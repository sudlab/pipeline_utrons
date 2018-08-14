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
import CGAT.Experiment as E
import CGAT.GTF as GTF
import CGAT.Bed as Bed
import CGAT.IOTools as IOTools
import itertools
import CGAT.Database as PUtils

import pandas


def getGeneTable(reffile):
    table = dict()
    for ens_gene in GTF.gene_iterator(GTF.iterator(IOTools.openFile(reffile))):
        geneid = ens_gene[0][0].gene_id
        table[geneid] = ens_gene
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

    # add common options (-h/--help, ...) and parse command line
    (options, args) = E.Start(parser, argv=argv)

    outlines = []
    individuals = []
    partnered = []
    individualpartnered = []
    novel = []

    db = pandas.read_csv(options.classfile, sep="\t")

    db = db.groupby("transcript_id").first()

    enshashtable = getGeneTable(options.reffile)
    
    for gene in GTF.gene_iterator(GTF.iterator(options.stdin)):
   
        transcript_ids = [transcript[0].transcript_id for transcript in gene]

        try:
            geneids = set(db.loc[transcript_ids].match_gene_id.dropna())
        except KeyError:
            continue
        
        geneids.discard(None)

        for geneid in geneids:  # get ensembl gene object from hashtable

            ens_gene = enshashtable[geneid]

            ens_introns = set()
            for transcript in ens_gene:
                ens_introns.update(GTF.toIntronIntervals(transcript))

            for first, second in itertools.product(gene, ens_gene):

                ens_CDS = GTF.asRanges(second, "CDS")
                if len(ens_CDS) == 0:  # ensure only protein-coding transcripts
                    continue
                ens_exons = GTF.asRanges(second, "exon")

                if second[0].strand == "+":
                    if first[0].start > ens_CDS[0][0]:
                        # ensure pruned transcript starts before ensembl CDS
                        continue
                else:
                    if first[-1].end < ens_CDS[-1][1]:
                        continue

                first_introns = set(GTF.toIntronIntervals(first))
                second_introns = set(GTF.toIntronIntervals(second))

                first_CDSintrons = [intron for intron in first_introns if
                                    (intron[0] > ens_CDS[0][0] and
                                     intron[1] < ens_CDS[-1][1])]

                second_CDSintrons = [intron for intron in second_introns if
                                     (intron[0] > ens_CDS[0][0] and
                                      intron[1] < ens_CDS[-1][1])]

                first_CDSintrons = set(first_CDSintrons)
                second_CDSintrons = set(second_CDSintrons)

                if not (len(first_CDSintrons - second_CDSintrons) == 0 and
                        len(second_CDSintrons - first_CDSintrons) == 0):
                    continue                           # match CDS intron chain

                firstUTRintrons = first_introns - first_CDSintrons

                if len(firstUTRintrons) == 0:
                    continue
                
                secondUTRintrons = second_introns - second_CDSintrons

                if second[0].strand == "+":
                    firstUTR5introns = set([intron for intron in firstUTRintrons
                                            if intron[1] < ens_CDS[0][0]])
                    secondUTR5introns = set([intron for intron in secondUTRintrons 
                                             if intron[1] < ens_CDS[0][0]])
                    if not ((first[0].start >= ens_exons[0][0]) or 
                            (len(firstUTR5introns-secondUTR5introns) == 0 and
                             len(secondUTR5introns-firstUTR5introns) == 0)):
                        continue                # ensure pruned transcript either
        # starts at same place or after ensembl transcript or has matching 5'UTR introns
                else: 
                    firstUTR5introns = set([intron for intron in firstUTRintrons 
                                            if intron[0] > ens_CDS[-1][1]])
                    secondUTR5introns = set([intron for intron in secondUTRintrons 
                                             if intron[0] > ens_CDS[-1][1]])
                    if not ((first[-1].end <= ens_exons[-1][1]) or
                            (len(firstUTR5introns-secondUTR5introns) == 0 and
                             len(secondUTR5introns-firstUTR5introns) == 0)):
                        continue

                for intron in firstUTRintrons:
                    if (intron[0] < ens_CDS[-1][1] and
                       intron[1] > ens_CDS[-1][1]) or \
                       (intron[0] < ens_CDS[0][0] and
                       intron[1] > ens_CDS[0][0]):
                        
                        continue        # ensure pruned transcript doesn't have
                        # introns overlapping start or stop codons in ensembl
                        # transcript
                
                if second[0].strand == "+":
                    ens_stop = ens_CDS[-1][1]
                    UTR3introns = [intron for intron in firstUTRintrons if
                                   intron[0] > ens_CDS[-1][1] and
                                   intron[1] < ens_exons[-1][1]]
                else:
                    ens_stop = ens_CDS[0][0]
                    UTR3introns = [intron for intron in firstUTRintrons if
                                   intron[1] < ens_CDS[0][0] and
                                   intron[0] > ens_exons[0][0]]

                if UTR3introns == []:
                    continue

                outbed = Bed.Bed()
                outbed.fields = ['.', '.', '.', '.', '.', '.', '.', '.', '.']
                outbed.fromIntervals(UTR3introns)
                outbed.contig = gene[0][0].contig
                outbed["name"] = first[0].transcript_id
                outbed["strand"] = gene[0][0].strand
                outlines.append(outbed)        # get output for each transcript
                    
                for item in UTR3introns:
                    outbed2 = Bed.Bed()
                    outbed2.fields = ['.', '.', '.', '.']
                    outbed2.fromIntervals([item])
                    outbed2.contig = gene[0][0].contig
                    outbed2['name'] = first[0].transcript_id
                    outbed2["strand"] = gene[0][0].strand
                    outbed2["thickStart"] = ens_stop
                    individuals.append(outbed2)  # get output for each intron
          
                UTR3introns = set(UTR3introns)
                extraUTR3introns = list(UTR3introns - secondUTRintrons)
                # get only introns that are not in matched transcript
                if len(extraUTR3introns) != 0:
                    outbed3 = Bed.Bed()
                    outbed3.fields = ['.'] * 9
                    outbed3.fromIntervals(extraUTR3introns)
                    outbed3.contig = gene[0][0].contig
                    outbed3["name"] = first[0].transcript_id + ":" + second[0].transcript_id
                    outbed3["strand"] = gene[0][0].strand
                    partnered.append(outbed3)
                        
                    for item in extraUTR3introns:
                        outbed4 = Bed.Bed()
                        outbed4.fields = ['.', '.', '.', '.']
                        outbed4.fromIntervals([item])
                        outbed4.contig = gene[0][0].contig
                        outbed4["name"] = first[0].transcript_id + ":" + second[0].transcript_id
                        outbed4["strand"] = gene[0][0].strand
                        outbed4["thickStart"] = ens_stop
                        individualpartnered.append(outbed4)
                
                if len(ens_introns) == 0:
                    ens_starts, ens_ends = [], []
                else:
                    ens_starts, ens_ends = zip(*ens_introns)

                novelEvents = [i for i in UTR3introns if
                               i[0] not in ens_starts and
                               i[1] not in ens_ends]
                
                for item in novelEvents:
                    outbed5 = Bed.Bed()
                    outbed5.fields = ['.']*4
                    outbed5.fromIntervals([item])
                    outbed5.contig = gene[0][0].contig
                    outbed5["name"] = first[0].transcript_id + ":" + second[0].transcript_id
                    outbed5["strand"] = gene[0][0].strand
                    outbed5["thickStart"] = ens_stop
                    novel.append(outbed5)

    with IOTools.openFile(options.outfile, "w") as outf:
        for line in outlines:
            outf.write(str(line)+"\n")
  
    if options.indivfile is not None:
        with IOTools.openFile(options.indivfile, "w") as outf2:
            for line in individuals:
                outf2.write(str(line)+"\n")

    if options.partfile is not None:
        with IOTools.openFile(options.partfile, "w") as outf3:
            for line in partnered:
                outf3.write(str(line)+"\n")
        
    if options.indivpartfile is not None:
        with IOTools.openFile(options.indivpartfile, "w") as outf4:
            for line in individualpartnered:
                outf4.write(str(line)+"\n")

    if options.novelfile is not None:
        with IOTools.openFile(options.novelfile, "w") as outf5:
            for line in novel:
                outf5.write(str(line)+"\n")
    # write footer and output benchmark information.
    E.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))

