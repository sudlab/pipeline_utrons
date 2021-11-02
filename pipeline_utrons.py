##############################################################################
#
#   MRC FGU CGAT
#
#   $Id$
#
#   Copyright (C) 2009 Andreas Heger
#
#   This program is free software; you can redistribute it and/or
#   modify it under the terms of the GNU General Public License
#   as published by the Free Software Foundation; either version 2
#   of the License, or (at your option) any later version.
#
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with this program; if not, write to the Free Software
#   Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
###############################################################################
"""
===========================
Pipeline template
===========================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python
:Updated for usage with python3 by: Cristina Alexandru

Overview
========

The Utrons pipeline is a tool for identifying introns in 3' UTR regions.

The pipeline performs the following:

    * Uses sequence alignment data from .bam files to assemble transcripts with StringTie.
    * Analyses and quantifies splice junctions in .bam files, filtering those which are unlinkely to be genuine.
    * Merges all assembled transcripts with StringTie.
    * Classifies and load transcripts into a databases.
    * Uses a separately custom-built pipeline to find utrons and their IDs.
    * Makes a Salmon index and quantifies transcripts with Salmon.
    * Merges all quantification files and uploads the outputs onto database.

Background
============

StringTie
---------

StringTie assembles transcripts from RNA seq reads aligned to the reference and performs quantification. 
It follows a netflow algorithm where it assembles and simultaneously quantifies the highly expressed 
transcripts, removing reads associated with them and repeating the process until all the reads are used. 
If provided with a reference annotation file, Stringtie uses it to construct assembly for low abundance genes.
Alternatively, the assembly of novel genes and transcripts may be skipped, using StringTie simply to quantify 
all the transcripts provided in an annotation file

Portcullis
---------

Portcullis stands for PORTable CULLing of Invalid Splice junctions from pre-aligned RNA-seq data. 
It is known that RNAseq mapping tools generate many invalid junction predictions, particularly in deep datasets 
with high coverage over splice sites. In order to address this, instead for creating a new RNAseq mapper, 
with a focus on SJ accuracy Portcullis takes in a BAM file generated by an RNAseq mapper of the user's own choice 
(e.g. Tophat2, Gsnap, STAR2 or HISAT2) as input (i.e. it's portable). It then, analyses and quantifies all splice 
junctions in the file before, filtering (culling) those which are unlikely to be genuine. Portcullis outputs 
junctions in a variety of formats, making it suitable for downstream analysis (such as differential splicing analysis
and gene modelling), without additional work. 
Portcullis can also filter the original BAM file, removing alignments associated with bad junctions. Both the filtered
junctions and BAM files are cleaner and more usable resources, which can more effectively be used to assist in 
downstream analyses such as gene prediction and genome annotation.

Usage
=====

See :ref:`PipelineSettingUp` and :ref:`PipelineRunning` on general
information how to use CGAT pipelines.

Configuration
-------------

The pipeline requires a configured :file:`pipeline.ini` file.
CGATReport report requires a :file:`conf.py` and optionally a
:file:`cgatreport.ini` file (see :ref:`PipelineReporting`).

Default configuration files can be generated by executing:

   python <srcdir>/pipeline_utrons.py config

Input files
-----------

1. RNA-se read alignments in a .bam format, in order to generate assembled transcripts in .gtf format.
2. A transcriptome in .fa format, in order to build a salmon index.
3. RNA-seq reads in fastq.1.gz and fastq.2.gz formats, for generating salmon quantification files in .sf format.

Requirements
------------

The pipeline requires the results from
:doc:`pipeline_annotations`. Set the configuration variable
:py:data:`annotations_database` and :py:data:`annotations_dir`.

On top of the default CGAT setup, the pipeline requires the following
software to be in the path:

.. Add any additional external requirements such as 3rd party software
   or R modules below:

Requirements:

* samtools >= 1.1

On top of the default CGAT setup, the pipeline requires the following
software to be in the path:

+--------------+----------+------------------------------------+
|*Program*     |*Version* |*Purpose*                           |
+--------------+----------+------------------------------------+
|gff2fasta     |          |bulding transcriptome in .fa format |
+--------------+----------+------------------------------------+
|salmon_       |>=0.7.2   |building an index                   |
+--------------+----------+------------------------------------+
|salmon_       |>=0.7.2   |alignment-free quantification       |
+--------------+----------+------------------------------------+

Pipeline output
===============

1. Assembled transcripts are generated in .gtf format.
2. A salmon index is built.
3. Salmon quantification files in .sf format are generted in quantification.dir.

Glossary
========

.. glossary::

salmon
      salmon_ - alignment-free quantification

.. _salmon: https://combine-lab.github.io/salmon/

###########################################################################


Code
====

"""
from ruffus import *
from ruffus.combinatorics import product
import sys
import os
import shutil
from gffutils import DataIterator as DataIterator
import sqlite3
import subprocess
import glob
from cgatcore import experiment as E
import cgat.Sra as Sra
from cgatcore import pipeline as P
import cgatpipelines.tasks.rnaseq as RnaSeq
import tempfile
from cgatpipelines.report import run_report

# load options from the config file
PARAMS = P.get_parameters(
    ["%s/pipeline.yml" % os.path.splitext(__file__)[0],
     "../pipeline.yml",
     "pipeline.yml"])

# add configuration values from associated pipelines
#
# 1. pipeline_annotations: any parameters will be added with the
#    prefix "annotations_". The interface will be updated with
#    "annotations_dir" to point to the absolute path names.
PARAMS.update(P.peek_parameters(
    PARAMS["annotations_dir"],
    'genesets',
    prefix="annotations_",
    update_interface=True,
    restrict_interface=True))

PARAMS["project_src"]=os.path.dirname(__file__)

# if necessary, update the PARAMS dictionary in any modules file.
# e.g.:
#
# import CGATPipelines.PipelineGeneset as PipelineGeneset
# PipelineGeneset.PARAMS = PARAMS
#
# Note that this is a hack and deprecated, better pass all
# parameters that are needed by a function explicitely.
RnaSeq.PARAMS = PARAMS

# -----------------------------------------------
# Utility functions
def connect():
    '''utility function to connect to database.

    Use this method to connect to the pipeline database.
    Additional databases can be attached here as well.
    Returns an sqlite3 database handle.
    '''

    dbh = sqlite3.connect(PARAMS["database_name"])
    statement = '''ATTACH DATABASE '%s' as annotations''' % (
        PARAMS["annotations_database"])
    cc = dbh.cursor()
    cc.execute(statement)
    cc.close()

    return dbh

@P.cluster_runnable
def create_index(table, columns):
    ''' utility function to create table indexes. Supply a list of columns
    to create a join index on those columns'''
    
    dbh = connect()
    cc = dbh.cursor()
    
    index_name = table + "_index_" + "_".join(columns)
    columns = ",".join(columns)
    cc.execute("CREATE INDEX %(index_name)s ON %(table)s(%(columns)s)" 
               % locals())
    
STRINGTIE_QUANT_FILES=["i_data.ctab", "e_data.ctab", "t_data.ctab",
                       "i2t.ctab", "e2t.ctab"]


# ---------------------------------------------------
@follows(mkdir("assembled_transcripts.dir"), mkdir("portcullis"))
@transform(["input_assemble.dir/*.bam",
            "input_assemble.dir/*.remote"],
           formatter(),
           add_inputs(os.path.join(
               PARAMS["annotations_dir"],
               PARAMS["annotations_interface_geneset_all_gtf"])),
           "assembled_transcripts.dir/{basename[0]}.gtf.gz")
def assembleWithStringTie(infiles, outfile):

    infile, reference = infiles
    basefile = os.path.basename(infile)
    job_threads = PARAMS["stringtie_threads"]
    job_memory = PARAMS["stringtie_memory"]
    tmpfile = P.get_temp_filename()
    if os.path.exists(tmpfile): 
    	os.unlink(tmpfile)
    
    statement =  '''
                    portcullis full 
			    -t 1
                            -o portcullis/%(basefile)s/
                            -r %(portcullis_bedref)s
                            -b
                            %(portcullis_fastaref)s
                            %(infile)s &&
                    mv portcullis/%(basefile)s/portcullis.filtered.bam %(tmpfile)s &&
                    rm -r portcullis/%(basefile)s/ &&
                    stringtie %(tmpfile)s
                           -p %(stringtie_threads)s
                           -G <(zcat %(reference)s)
                           %(stringtie_options)s
                           2> %(outfile)s.log
                   | gzip > %(outfile)s &&
                   rm %(tmpfile)s'''

    if infile.endswith(".remote"):
        token = glob.glob("gdc-user-token*")
        tmpfilename = P.get_temp_filename()
        if os.path.exists(tmpfilename):
            os.unlink(tmpfilename)    	
        if len(token) > 0:
            token = token[0]
        else:
            token = None

        s, infile = Sra.process_remote_BAM(
            infile, token, tmpfilename,
            filter_bed=os.path.join(
                PARAMS["annotations_dir"],
                PARAMS["annotations_interface_contigs_bed"]))

        infile = " ".join(infile)
        statement = "; ".join(
            ["mkdir -p %(tmpfilename)s",
             s,
             statement,
             "rm -r %(tmpfilename)s"])

    P.run(statement, job_condaenv="portcullis")


# ---------------------------------------------------
@follows(mkdir("final_genesets.dir"), assembleWithStringTie)
@merge([assembleWithStringTie,
       os.path.join(
           PARAMS["annotations_dir"],
           PARAMS["annotations_interface_geneset_all_gtf"])],
       "final_genesets.dir/agg-agg-agg.gtf.gz")
def mergeAllAssemblies(infiles, outfile):

    infiles = ["<(zcat %s)" % infile for infile in infiles]
    infiles, reference = infiles[:-1], infiles[-1]

    job_threads = PARAMS["stringtie_merge_threads"]
    job_memory = PARAMS["stringtie_merge_memory"]

    infiles = " ".join(infiles)

    statement = '''stringtie --merge
                             -G %(reference)s
                             -p %(stringtie_merge_threads)s
                             %(stringtie_merge_options)s
                             %(infiles)s
                            2> %(outfile)s.log
                   | cgat gtf2gtf --method=sort
                           --sort-order=gene+transcript
                            -S %(outfile)s -L %(outfile)s.log'''

    P.run(statement) 

@follows(assembleWithStringTie)
@collate([glob.glob("assembled_transcripts.dir/*%s*.gtf.gz" % PARAMS["stringtie_groups"][x]) for x in range(0, len(PARAMS["stringtie_groups"]))],
         regex("(.+)/(.+)-(.+).gtf.gz"),
         add_inputs(os.path.join(
            PARAMS["annotations_dir"],
            PARAMS["annotations_interface_geneset_all_gtf"])),
          r"final_genesets.dir/\2-agg-agg-agg.gtf.gz")
def merge_by_tissue(infiles, outfile):
    job_threads = PARAMS["stringtie_merge_threads"]
    job_memory= PARAMS["stringtie_merge_memory"]

    reference = "<(zcat %s)" % infiles[0][1]
    infiles = ["<(zcat %s)" % infile for infile in infiles[0][0]]
    infiles = " ".join(infiles)
    
    statement = '''stringtie --merge
                    -G %(reference)s
                    -p %(stringtie_merge_threads)s
                    %(stringtie_merge_options)s
                    %(infiles)s
                    2> %(outfile)s.log
               | cgat gtf2gtf --method=sort
                   --sort-order=gene+transcript
                   -S %(outfile)s -L %(outfile)s.log'''
    P.run(statement)


@follows(mergeAllAssemblies, merge_by_tissue)
def Assembly():
    pass


# ---------------------------------------------------
@transform([assembleWithStringTie, mergeAllAssemblies],
           suffix(".gtf.gz"),
           add_inputs(os.path.join(
               PARAMS["annotations_dir"],
               PARAMS["annotations_interface_geneset_all_gtf"])),
           ".class.gz")
def classifyTranscripts(infiles, outfile):
    '''classify transcripts.
    '''
    to_cluster = True

    infile, reference = infiles

    counter = PARAMS['gtf2table_classifier']

    job_memory = "16G"

    statement = '''
    zcat %(infile)s
    | cgat gtf2table
           --counter=%(counter)s
           --reporter=transcripts
           --gff-file=%(reference)s
           --log=%(outfile)s.log
    | gzip
    > %(outfile)s
    '''
    P.run(statement)


# ---------------------------------------------------
@follows(mkdir("database_load"), classifyTranscripts)
@merge(classifyTranscripts,
       "database_load/transcript_class.load")
def loadTranscriptClassification(infiles, outfile):

    P.concatenate_and_load(infiles, outfile,
                         regex_filename=".+/(.+).class.gz",
                         options="-i transcript_id -i gene_id"
                         " -i match_gene_id -i match_transcript_id"
                         " -i source",
                         job_memory="64G")


# ---------------------------------------------------
@follows(mkdir("utron_beds.dir"), classifyTranscripts)
@subdivide([assembleWithStringTie, mergeAllAssemblies],
           regex("(.+)/(.+).gtf.gz"),
           add_inputs(os.path.join(
               PARAMS["annotations_dir"],
               PARAMS["annotations_interface_geneset_all_gtf"]),
                      r"\1/\2.class.gz"),
           [r"utron_beds.dir/\2.all_utrons.bed.gz",
            r"utron_beds.dir/\2.partnered_utrons.bed.gz",
            r"utron_beds.dir/\2.novel_utrons.bed.gz"])
def find_utrons(infiles, outfiles):

    infile, reference, classfile = infiles
    job_threads=2
    job_memory="16G"

    all_out, part_out, novel_out = outfiles

    track = P.snip(all_out, ".all_utrons.bed.gz")
    current_file = __file__ 
    pipeline_path = os.path.abspath(current_file)
    pipeline_directory = os.path.dirname(pipeline_path)
    script_path = "pipeline_utrons/find_utrons.py"
    full_utron_path = os.path.join(pipeline_directory, script_path)  
    statement = '''cgat gtf2gtf -I %(infile)s
                             --method=sort
                             --sort-order=gene+transcript
                              -L %(track)s.log
                 | python full_utron_path 
                             --reffile=%(reference)s
                             --class-file=%(classfile)s
                             --outfile %(all_out)s
                             --partfile=%(part_out)s
                             --novel-file=%(novel_out)s
                              -L %(track)s.log'''

    P.run(statement)


# ---------------------------------------------------
@transform(find_utrons,
           suffix(".bed.gz"),
           ".ids.gz")
def getUtronIds(infile, outfile):

    statement = '''zcat %(infile)s 
                 | cut -f 4
                 | sed 's/:/\\t/g'
                 | sort -u
                 | gzip > %(outfile)s'''

    P.run(statement)


# ---------------------------------------------------
@collate(getUtronIds,
         regex(".+/(.+)\.(.+).ids.gz"),
         r"database_load/\2_ids.load")
def loadUtronIDs(infiles, outfile):
    job_threads=3

    header = "track,transcript_id"
    options = "-i track -i transcript_id"

    if not outfile == "all_utrons_ids.load":
        header += ",match_transcript_id"
        options += " -i match_transcript_id"

    P.concatenate_and_load(infiles, outfile,
                         regex_filename=".+/(.+)\..+\.ids.gz",
                         has_titles=False,
                         cat="track",
                         header=header,
                         options=options,
                         job_memory="64G")


@follows(loadUtronIDs, loadTranscriptClassification)
def AnnotateAssemblies():
    pass

# ---------------------------------------------------
@follows(mkdir("export/indexed_gtfs.dir"))
@transform([assembleWithStringTie,
            mergeAllAssemblies,
            merge_by_tissue,
            "final_genesets.dir/*.gtf.gz"],
           regex("(.+)/(.+).gtf.gz"),
           r"export/indexed_gtfs.dir/\2.gtf.gz")
def exportIndexedGTFs(infile, outfile):

    statement = '''zcat %(infile)s
                 | sort -k1,1 -k4,4n
                 | bgzip > %(outfile)s;
                 tabix -p gff %(outfile)s;
                 if [[ %(outfile)s == *"agg-agg-agg"* ]]; then cp %(outfile)s %(outfile)s.tbi export/; fi'''

    P.run(statement)


@follows(exportIndexedGTFs)
def export():
    pass


# ---------------------------------------------------
@follows(mkdir("salmon_index"),exportIndexedGTFs)
@transform("final_genesets.dir/agg-agg-agg.gtf.gz",
           formatter(),
           "salmon_index/{basename[0]}.salmon.index")
def makeSalmonIndex(infile,outfile):
    # Long transcripts cause indexing to use lots of memory?
    job_memory="64G"
    job_threads=1

    gtf_basename = P.snip(os.path.basename(infile), ".gtf.gz")
    transcript_fasta = "salmon_index/" + gtf_basename + "transcripts.fa"
    fastaref =PARAMS["portcullis_fastaref"]
    index_options=PARAMS["salmon_indexoptions"]
    tmpfile = P.get_temp_filename()
    


    statement = '''
    gunzip -c %(infile)s > %(tmpfile)s;
    gffread %(tmpfile)s -g %(fastaref)s -w %(transcript_fasta)s;
    salmon index
      -p %(job_threads)s
      %(index_options)s
      -t %(transcript_fasta)s
      -i %(outfile)s
      --perfectHash;
    rm %(tmpfile)s
    '''
    P.run(statement)


#---------------------------------------------------------

if not os.path.exists("temp_bams/"):
    os.makedirs("temp_bams/")

@follows(mkdir("quantification.dir"), mkdir("sorted_bams"),
         mergeAllAssemblies, makeSalmonIndex)
@product(["input_assemble.dir/*.bam",
          "input_assemble.dir/*.remote"],
         formatter(".+/(?P<TRACK>.+).([bam|remote])"),
         mergeAllAssemblies,
         formatter(".+/(?P<GENESET>.+).gtf.gz"),
         "quantification.dir/{TRACK[0][0]}_{GENESET[1][0]}")
def quantifyWithSalmon(infiles, outfile):
    '''Quantify existing samples against genesets'''
    job_threads=2
    job_memory="24G"

    infile, gtffile = infiles
    outdir = P.snip(outfile, ".sf")
    basefile = os.path.basename(infile)
    sample_name = basefile.split(os.extsep, 1)
    gtfbase = P.snip(os.path.basename(gtffile), ".gz")
    salmonIndex = "salmon_index/" + gtfbase + ".salmon.index"
    salmon_options=PARAMS["salmon_quantoptions"]

    sorted_bam="sorted_bams/" + sample_name[0] + "_sorted.bam"
    fastq1 = P.snip(outfile, "_agg-agg-agg.sf")+".1.fastq"
    fastq2 = P.snip(outfile, "_agg-agg-agg.sf")+".2.fastq"
    fastq0 = P.snip(outfile, "_agg-agg-agg.sf")+".0.fastq"

    statement = '''
    samtools sort -n %(infile)s -o %(sorted_bam)s &&
    samtools fastq
         -1 %(fastq1)s
         -2 %(fastq2)s
         -0 %(fastq0)s -s /dev/null -n -F 0x900
         %(sorted_bam)s &&
    paired_end_reads=$(samtools view -c -f 1 %(sorted_bam)s) &&
    if [ $paired_end_reads = 0 ]; then
        salmon quant -i %(salmonIndex)s
            --libType IU
            -r %(fastq0)s
            -o %(outdir)s
            %(salmon_options)s &&
    else
        salmon quant -i %(salmonIndex)s
            --libType IU
            -1 %(fastq1)s
            -2 %(fastq2)s
            -o %(outdir)s
            %(salmon_options)s 
    fi &&
    mv %(outdir)s/quant.sf %(outfile)s && 
    rm %(fastq1)s && rm %(fastq2)s && rm %(fastq0)s && rm %(sorted_bam)s 
    '''

    if infile.endswith(".remote"):
        token = glob.glob("gdc-user-token*")
        filename = "temp_bams/%s" % basefile
        tmpfilename = P.get_temp_filename()
        if os.path.exists(tmpfilename):
            os.unlink(tmpfilename)

        if len(token) > 0:
            token = token[0]
        else:
            token = None

        s, infile = Sra.process_remote_BAM(
            infile, token, filename,
            filter_bed=os.path.join(
                PARAMS["annotations_dir"],
                PARAMS["annotations_interface_contigs_bed"]))

        infile = " ".join(infile)
        statement = "&& ".join(
            ["mkdir %(filename)s",
             s,
             statement,
             "rm -r %(filename)s"])



    P.run(statement)

#----------------------------------------------------
@follows(quantifyWithSalmon)
@merge("quantification.dir/*_agg_agg_agg/logs/salmon_quant.log",
       "mapping_rates.txt")
def merge_mapping_rates(infiles, outfile):
    '''Find the mapping rates in the salmon logs'''

    with open(outfile, "w") as outf:
        for inf in infile:
            for line in open(infile):
                if re.match("Mapping rate", line):
                    outf.write(line)

#----------------------------------------------------
# ***********
# NOTE - last time I used this, the output database wasn't indexed correctly
# ***********
@merge(quantifyWithSalmon, "database_load/salmon_quant.load")
def mergeAllQuants(infiles, outfile):
    job_threads=3

    P.concatenate_and_load(infiles, outfile,
                         regex_filename="quantification.dir/(.*)_agg-agg-agg.sf",
                         options="-i Name -i Length -i EffectiveLength"
                         " -i TPM -i NumReads -i track"
                         " -i source",
                         job_memory="64G")


#---------------------------------------------------------------------------------------------------------------
###### Export all_utrons, novel_utrons ids and tx2gene text files from utrons database

@follows(mergeAllQuants, mkdir("expression.dir", "expression.dir/csvdb_files"))
def CSVDBfiles():
    '''utility function to connect to database.

    Use this method to connect to the pipeline database.
    Export all_utrons_ids, novel_utrons_ids and tx2gene data in :term:`txt` format
    to be used further on in the Rscript
    '''

    subprocess.call(["sqlite3", PARAMS["database_name"],
                     ".headers on", ".mode tab", ".output expression.dir/csvdb_files/tx2gene.txt",
                     "select transcript_id, match_gene_id from transcript_class where track = 'agg-agg-agg'"])

    subprocess.call(["sqlite3", PARAMS["database_name"],
                     ".headers on", ".mode tab", ".output expression.dir/csvdb_files/all_utrons_ids.txt", 
                     "select * from all_utrons_ids where track = 'agg-agg-agg'"])

    subprocess.call(["sqlite3", PARAMS["database_name"],
                     ".headers on", ".mode tab", ".output expression.dir/csvdb_files/partnered_utrons_ids.txt",
                     "select * from partnered_utrons_ids where track = 'agg-agg-agg'"])

    subprocess.call(["sqlite3", PARAMS["database_name"],
                     ".headers on", ".mode tab", ".output expression.dir/csvdb_files/novel_utrons_ids.txt",
                     "select * from novel_utrons_ids where track = 'agg-agg-agg'"])
    
    shutil.copy(PARAMS["database_false_positives"], 
                "expression.dir/csvdb_files/60db_novel_utrons_ids.txt")

###### Identify splice sites

@follows(find_utrons, mergeAllAssemblies, mergeAllQuants)
@transform("utron_beds.dir/*.bed.gz", 
           regex("(.+)/agg-agg-agg.(.+)_utrons.bed.gz"), 
           r"expression.dir/\2_splice_sites.txt")
def identify_splice_sites(infiles, outfiles):
    infile=infiles
    outfile, outfile_load = outfiles
    current_file = __file__
    pipeline_path = os.path.abspath(current_file)
    pipeline_directory = os.path.dirname(pipeline_path)
    script_path = os.path.join(pipeline_directory, 
                               "/splicesites_start_end_sizes.py")
    fastaref =PARAMS["portcullis_fastaref"]

    statement = ''' python %(script_path)s -g %(fastaref)s 
                           -I %(infile)s %(outfile)s'''
    P.run(statement)   

@transform(identify_splice_sites, suffix(".tsv"), ".load")
def load_splice_sites(infile, outfile):
    P.load(infile, 
           outfile, 
           options = "-i transcript_id "
                     "-i ss5 "
                     "-i ss3 "
                     "-i splice_site_start "
                     "-i splice_site_end "
                     "-i utron_size", 
            job_memory="16G")


@follows(find_utrons, mergeAllAssemblies, mergeAllQuants)
@transform(PARAMS["annotations_interface_geneset_all_gtf"],
           regex("(.+).gtf.gz"),
           "expression.dir/gtf_stop_codons.txt")
def gtf_stop_codons(infile, gtf):
    '''Extract stop codon positions from reference gtf'''

    outfile = open(gtf, "w")
    for line in DataIterator(infile):
        if line.featuretype == "stop_codon":
            if line.strand == "+":
                outfile.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % 
                     (line.attributes['gene_id'],
                      line.attributes['transcript_id'],
                      line.attributes['gene_name'],
                      line.strand, 
                      line.featuretype, 
                      line.seqid, 
                      line.start, 
                      line.end))
            elif line.strand =="-":
                 outfile.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % 
                     (line.attributes['gene_id'],
                      line.attributes['transcript_id'],
                      line.attributes['gene_name'],
                      line.strand, 
                      line.featuretype, 
                      line.seqid, 
                      line.end, 
                      line.start))
        else:
            pass
    outfile.close()

######### Rscript for generating tpm expression values for transcript and 
# gene level, as well as fraction expression ###########

@follows(CSVDBfiles, quantifyWithSalmon)
@merge("quantification.dir/*.sf",
       ["database_load/utrons_expression.load", 
        "database_load/partnered_utrons_expression.load", 
        "database_load/novel_utrons_expression.load"])
def utrons_expression(infiles, outfiles):
    outfile_load, partnered_load, novel_load = outfiles
    job_threads = 4
    job_memory = "48G"
    script_file = os.path.join(P.snip(__file__, ".py"), "utrons_Rscript.R")
    outfile = "expression.dir/utrons_expression.txt"
    statement = '''Rscript %(script_file)s '''
    P.run(statement)
    
    P.load(outfile, 
           outfile_load, 
           options = "-i Sample -i transcript_id -i gene_id -i tr_expr "
                     "-i gene_expr -i fract_expr",
           job_memory="16G")

    partnered = "expression.dir/partnered_utrons_expression.txt"
    
    connect().cursor().executescript(
        '''DROP TABLE IF EXISTS partenered_utron_expression;
           CREATE TABLE partnered_utron_expression AS
           SELECT * 
             FROM utrons_expression
             WHERE transcript_id in 
             (SELECT transcript_id 
               FROM partnered_utrons_ids
               WHERE track='agg-agg-agg')''')
    
    for column in ("Sample", "transcript_id", "gene_id", "tr_expr", 
                   "gene_expr", "fract_expr"):
        create_index("partnered_utrons_expression", column, submit=True,
                     job_memory="16G")
        
    connect().cursor().executescript(
        '''DROP TABLE IF EXISTS partenered_utron_expression;
           CREATE TABLE partnered_utron_expression AS
           SELECT * 
             FROM 
                  utrons_expression
             WHERE 
                  transcript_id in 
                  (SELECT transcript_id 
                   FROM partnered_utrons_ids
                   WHERE track='agg-agg-agg')''')
    
    for column in ("Sample", "transcript_id", "gene_id", "tr_expr", 
                   "gene_expr", "fract_expr"):
        create_index("novel_utrons_expression", column, submit=True,
                     job_memory="16G")

# --------------------------------------------------------------------------------------------------------------------------------------------
# Generic pipeline tasks
@follows(Assembly, AnnotateAssemblies, export, mergeAllQuants, 
         CSVDBfiles, utrons_expression, identify_splice_sites, gtf_stop_codons)
def full():
    pass

####################################################

@follows(mkdir("MultiQC_report.dir"))
@originate("MultiQC_report.dir/multiqc_report.html")
def renderMultiqc(infile):
    '''build mulitqc report'''

    statement = (
        "export LANG=en_GB.UTF-8 && "
        "export LC_ALL=en_GB.UTF-8 && "
        "multiqc . -f && "
        "mv multiqc_report.html MultiQC_report.dir/")

    P.run(statement)


@follows(renderMultiqc)

################################################### added #

@follows(mkdir("report"))
def build_report():
    '''build report from scratch.

    Any existing report will be overwritten.
    '''

    E.info("starting report build process from scratch")
    run_report(clean=True)


@follows(mkdir("report"))
def update_report():
    '''update report.

    This will update a report with any changes inside the report
    document or code. Note that updates to the data will not cause
    relevant sections to be updated. Use the cgatreport-clean utility
    first.
    '''

    E.info("updating report")
    run_report(clean=False)


@follows(update_report)
def publish_report():
    '''publish report in the CGAT downloads directory.'''

    E.info("publishing report")
    P.publish_report()


def main(argv=None):
    if argv is None:
        argv = sys.argv
    P.main(argv)


if __name__ == "__main__":
    sys.exit(P.main(sys.argv))
