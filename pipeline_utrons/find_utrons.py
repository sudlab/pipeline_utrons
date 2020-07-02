'''
find_utrons.py
====================================================

:Author: Isabel
:Release: $1.0$
:Date: |today|
:Tags: Python

Purpose
-------

.. To find utrons (introns in 3' UTRs).

Usage
-----

.. Example use case

Example::

   python cgat_script_template.py

Type::

   python find_utrons.py --help

for command line help.

Command line options
--------------------

'''

import sys
from cgatcore import experiment as E
import cgat.GTF as GTF
import cgat.Bed as Bed
import cgatcore.iotools as IOTools
import itertools
from cgatcore import database as Database
from collections import defaultdict
import pandas


def getGeneTable(reffile):
    E.info("Loading reference")
    table = defaultdict(dict)
    for ens_gene in GTF.gene_iterator(GTF.iterator(IOTools.open_file(reffile))):
        geneid = ens_gene[0][0].gene_id
        table[geneid]["models"] = dict()
        table[geneid]["start_codons"] = defaultdict(list)
        
        for transcript in ens_gene:

            transcript_id = transcript[0].transcript_id
            table[geneid]["models"][transcript_id] = transcript
            
            CDS = GTF.asRanges(transcript, "start_codon")
            if len(CDS) == 0:
                continue

            if transcript[0].strand == "-":
                start_codon = max(e[1] for e in CDS)
            else:
                start_codon = min(e[0] for e in CDS)

            table[geneid]["start_codons"][start_codon].append(transcript_id)

    E.info("Reference Loaded")
    return table


def main(argv=None):
    """script main.
    parses command line options in sys.argv, unless *argv* is given.
    """

    if argv is None:
        argv = sys.argv

    # setup command line parser
    parser = E.OptionParser(version="%prog version: $1.0$",
                            usage=globals()["__doc__"])

    parser.add_option("-r", "--reffile", dest="reffile", type="string",
                      help="Supply reference gtf file name")

    parser.add_option("-d", "--class-file", dest="classfile", type="string",
                      help="Supply database name")

    parser.add_option("-o", "--outfile", dest="outfile", type="string",
                      help="Supply output bed file name")

    parser.add_option("-u", "--indivfile", dest="indivfile", type="string",
                      help="Supply output bed file name for individual utrons")

    parser.add_option("-p", "--partfile", dest="partfile", type="string",
                      help="Supply output bed file name for partnered utrons")
    parser.add_option("-q", "--indivpartfile", dest="indivpartfile", type="string",
                      help="Supply output bed file name for individual partnered utrons")
    parser.add_option("-n", "--novel-file", dest="novelfile", type="string",
                      help="Supply output bed file name for novel introns")
    parser.add_option("--novel-transcript", dest="novel_id", type="string",
                      help="DEBUG: Output info for this transcript from the STDIN")
    parser.add_option("--target-transcript", dest="target_id", type="string",
                      help="DEBUG: Output info for this transcript from ref-file")
    
    # add common options (-h/--help, ...) and parse command line
    (options, args) = E.start(parser, argv=argv)

    outlines = []
    individuals = []
    partnered = []
    individualpartnered = []
    novel = []

    db = pandas.read_csv(options.classfile, sep="\t")

    # This keeps just one entry per-transcript - why? 
    #db = db.groupby("transcript_id").first()
    db = db.set_index("transcript_id")
    enshashtable = getGeneTable(options.reffile)
    
    for novel_transcript in GTF.transcript_iterator(GTF.iterator(options.stdin)):

        # Why do it on a gene by gene basis rather than transcript by transcript basis?
        transcript_id = novel_transcript[0].transcript_id

        if transcript_id == options.novel_id:
            output_novel = True
        else:
            output_novel = False
        
        try:
            geneid = db.loc[transcript_id].match_gene_id
        except KeyError:
            if output_novel:
                E.debug("Transcript %s not in class table" % transcript_id)
            continue

        if pandas.isnull(geneid):
            if output_novel:
                E.debug("Transcript %s matches no gene in class table" % transcript_id)
            continue

        ens_gene = enshashtable[geneid]
            
        all_ref_introns = set()
        novel_transcript_exons = GTF.asRanges(novel_transcript, "exon")
        novel_transcript_introns = GTF.toIntronIntervals(novel_transcript)      
        for ref_transcript in ens_gene["models"].values():
            ref_introns = GTF.toIntronIntervals(ref_transcript)
            all_ref_introns.update(ref_introns)


        #Identify comparison set
        def _in_exon(position, exons):
            return any(e[0] <= position <= e[1] for e in exons)

        # check if this ever gets the wrong start_codon. 
        filtered_starts = [s for s in ens_gene["start_codons"] if
                           _in_exon(s, novel_transcript_exons)]

        if len(filtered_starts) == 0:
            if output_novel:
                E.debug("No starts found for %s" % transcript_id)
            continue
        
        #if novel_transcript[0].strand == "-":
        #    selected_start = max(filtered_starts)
        #else:
        #    selected_start = min(filtered_starts)

        selected_models = list()
        for startc in filtered_starts:
            selected_models.extend(ens_gene["start_codons"][startc])

        if output_novel:
            E.debug("Transcripts with compatible starts are %s" % selected_models)
            
        for ref_transcript_id in selected_models:

            if output_novel and ref_transcript_id == options.target_id:
                output_ref=True
            else:
                output_ref=False
                
            second = ens_gene["models"][ref_transcript_id]
            ens_CDS = GTF.asRanges(second, "CDS")
            
            if len(ens_CDS) == 0:
                if output_ref:
                    E.debug("%s is not coding") # ensure only protein-coding transcripts
                continue

            ens_exons = GTF.asRanges(second, "exon")
            
            first_introns = set(novel_transcript_introns)
            second_introns = set(GTF.toIntronIntervals(second))

            first_CDSintrons = [intron for intron in first_introns if
                                (intron[0] > ens_CDS[0][0] and
                                 intron[1] < ens_CDS[-1][1])]

            second_CDSintrons = [intron for intron in second_introns if
                                 (intron[0] > ens_CDS[0][0] and
                                  intron[1] < ens_CDS[-1][1])]

            first_CDSintrons = set(first_CDSintrons)
            second_CDSintrons = set(second_CDSintrons)

            if not first_CDSintrons == second_CDSintrons:
                if output_ref:
                    E.debug("CDS chains do not match. Chains are:")
                    first_CDSintrons = sorted(list(first_CDSintrons))
                    second_CDSintrons = sorted(list(second_CDSintrons))
                    output = "\n".join(map(str, zip(first_CDSintrons, second_CDSintrons)))
                    E.debug(output)
                continue                           # match CDS intron chain

                      
            firstUTRintrons = first_introns - first_CDSintrons

            if len(firstUTRintrons) == 0:
                if output_ref:
                    E.debug("No UTR introns")
                continue

            secondUTRintrons = second_introns - second_CDSintrons

            found = False
            for intron in first_introns:
                if (intron[0] < ens_CDS[-1][1] and
                    intron[1] > ens_CDS[-1][1]) or \
                    (intron[0] < ens_CDS[0][0] and
                     intron[1] > ens_CDS[0][0]):

                    found=True
                    break      # ensure pruned transcript doesn't have
                        # introns overlapping start or stop codons in ensembl
                        # transcript
            if found:
                if output_ref:
                    E.debug("Start or stop in intron")
                continue
            
            if second[0].strand == "+":
                ens_stop = ens_CDS[-1][1]
                UTR3introns = [intron for intron in firstUTRintrons if
                               intron[0] >= ens_CDS[-1][1] and
                               intron[1] < ens_exons[-1][1]]
                secondUTR3introns = [intron for intron in secondUTRintrons if
                                     intron[0] >= ens_CDS[-1][1] and
                                     intron[1] < ens_exons[-1][1]]
            else:
                ens_stop = ens_CDS[0][0]
                UTR3introns = [intron for intron in firstUTRintrons if
                               intron[1] <= ens_CDS[0][0] and
                               intron[0] > ens_exons[0][0]]
                secondUTR3introns =  [intron for intron in secondUTRintrons if
                                      intron[1] <= ens_CDS[0][0] and
                                      intron[0] > ens_exons[0][0]]

            if len(UTR3introns) == 0:
                if output_ref:
                    E.debug("No UTR introns")
                continue

            outbed = Bed.Bed()
            outbed.fields = ['.', '.', '.', '.', '.', '.', '.', '.', '.']
            outbed.fromIntervals(UTR3introns)
            outbed.contig = novel_transcript[0].contig
            outbed["name"] = novel_transcript[0].transcript_id
            outbed["strand"] = novel_transcript[0].strand
            outlines.append(outbed)        # get output for each transcript
                    
            for item in UTR3introns:
                outbed2 = Bed.Bed()
                outbed2.fields = ['.', '.', '.', '.']
                outbed2.fromIntervals([item])
                outbed2.contig = novel_transcript[0].contig
                outbed2['name'] = novel_transcript[0].transcript_id
                outbed2["strand"] = novel_transcript[0].strand
                outbed2["thickStart"] = ens_stop
                individuals.append(outbed2)  # get output for each intron
          
            UTR3introns = set(UTR3introns)
            secondUTR3introns = set(secondUTR3introns)
            extraUTR3introns = list(UTR3introns - secondUTR3introns)

            
            if output_ref and len(secondUTR3introns -  UTR3introns) > 0:
                E.debug("Following introns in UTR of %s but not %s" % (options.target_id, options.novel_id))
                E.debug(secondUTRintrons - UTR3introns)
                
            # get only introns that are not in matched transcript
            if len(extraUTR3introns) != 0 and len(secondUTR3introns - UTR3introns) == 0:
                outbed3 = Bed.Bed()
                outbed3.fields = ['.'] * 9
                outbed3.fromIntervals(extraUTR3introns)
                outbed3.contig = novel_transcript[0].contig
                outbed3["name"] = novel_transcript[0].transcript_id + ":" + second[0].transcript_id
                outbed3["strand"] = novel_transcript[0].strand
                partnered.append(outbed3)
                        
                for item in extraUTR3introns:
                    outbed4 = Bed.Bed()
                    outbed4.fields = ['.', '.', '.', '.']
                    outbed4.fromIntervals([item])
                    outbed4.contig = novel_transcript[0].contig
                    outbed4["name"] = novel_transcript[0].transcript_id + ":" + second[0].transcript_id
                    outbed4["strand"] = novel_transcript[0].strand
                    outbed4["thickStart"] = ens_stop
                    individualpartnered.append(outbed4)
                
            if len(all_ref_introns) == 0:
                ens_starts, ens_ends = [], []
            else:
                ens_starts, ens_ends = zip(*all_ref_introns)

            novelEvents = [i for i in UTR3introns if
                           i[0] not in ens_starts and
                           i[1] not in ens_ends]
                
            for item in novelEvents:
                outbed5 = Bed.Bed()
                outbed5.fields = ['.']*4
                outbed5.fromIntervals([item])
                outbed5.contig = novel_transcript[0].contig
                outbed5["name"] = novel_transcript[0].transcript_id + ":" + second[0].transcript_id
                outbed5["strand"] = novel_transcript[0].strand
                outbed5["thickStart"] = ens_stop
                novel.append(outbed5)

    with IOTools.open_file(options.outfile, "w") as outf:
        for line in outlines:
            outf.write(str(line)+"\n")
  
    if options.indivfile is not None:
        with IOTools.open_file(options.indivfile, "w") as outf2:
            for line in individuals:
                outf2.write(str(line)+"\n")

    if options.partfile is not None:
        with IOTools.open_file(options.partfile, "w") as outf3:
            for line in partnered:
                outf3.write(str(line)+"\n")
        
    if options.indivpartfile is not None:
        with IOTools.open_file(options.indivpartfile, "w") as outf4:
            for line in individualpartnered:
                outf4.write(str(line)+"\n")

    if options.novelfile is not None:
        with IOTools.open_file(options.novelfile, "w") as outf5:
            for line in novel:
                outf5.write(str(line)+"\n")
    # write footer and output benchmark information.
    E.stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))

