'''
Finds the splice sites of the each utron sequence, as well as the size of each,
and its start and end (strand corrected)

usage:

splicesites_start_end_sizes.py [OPTIONS] -I INFILE -S OUTFILE -O PER_UTRON_OUT
'''
import sys

from cgat import Bed
from cgat import IndexedFasta
from cgatcore import iotools
from cgat import Genomics

from cgatcore import expriment as E

parser = E.OptionParser(version="%prog version: $1.0$",
                           usage=globals()["__doc__"])
parser.add_option("-g", "--genome", dest="genome",
                  help="index fasta genome sequence")
parser.add_option("-O", "--per-utron-out", dest="outfile",
                  help="File name for output file that will contain one row"
                       "per entry in the input")
options, args = E.start(parser, sys.argv)
                  
genome = IndexedFasta.IndexedFasta(options["genome"])


bedfile = Bed.iterator(options.stdin)
splice_site_dict = dict()
outfile = iotools.open_file(options["outfile"], "w")
outfile.write("\t".join("transcript_id",
                        "strand",
                        "ss5",
                        "ss3",
                        "contig",
                        "splice_site_start",
                        "splice_site_end",
                        "utron_size"))
for utron in bedfile:
    
    ss5_sequence = genome.getSequence(utron.contig, "+", utron.start, utron.start+2)
    ss3_sequence = genome.getSequence(utron.contig, "+", utron.end-2, utron.end)
    if utron.strand == "+":
        splice_site_dict[utron.name] = (ss5_sequence, ss3_sequence)
        if ":" in utron.name:
            transcript_id = utron.name.split(":")[0]
            match_transcript_id = utron.name.split(":")[1]
            outfile.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (transcript_id, utron.strand, ss5_sequence, ss3_sequence, utron.contig, utron.start, utron.end, utron.end-utron.start))
        else:
            outfile.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (utron.name, utron.strand, ss5_sequence, ss3_sequence, utron.contig, utron.start, utron.end, utron.end-utron.start))

    elif utron.strand == "-":
        ss5_sequence = Genomics.reverse_complement(ss5_sequence)
        ss3_sequence = Genomics.reverse_complement(ss3_sequence)
        splice_site_dict[utron.name] = (ss3_sequence, ss5_sequence)
        if ":" in utron.name:
            transcript_id = utron.name.split(":")[0]
            match_transcript_id = utron.name.split(":")[1]
            outfile.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (transcript_id, utron.strand, ss3_sequence, ss5_sequence, utron.contig, utron.end, utron.start, utron.end-utron.start))
        else:
            outfile.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (utron.name, utron.strand, ss3_sequence, ss5_sequence, utron.contig, utron.end, utron.start, utron.end-utron.start))

outfile.close()

from collections import defaultdict

counts = defaultdict(int)
for name, ss in splice_site_dict.items():
    counts[(ss[0].upper(), ss[1].upper())] += 1

sorted_keys = sorted(counts.keys(), key = lambda x: counts[(x)])    
for key in sorted_keys:
    print(key, counts[key])

