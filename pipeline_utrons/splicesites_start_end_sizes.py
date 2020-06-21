
# coding: utf-8

# In[2]:

import sys

from cgat import Bed
from cgat import IndexedFasta
from cgatcore import iotools
from cgat import Genomics


# In[3]:


genome = IndexedFasta.IndexedFasta("/shared/sudlab1/General/mirror/genomes/plain/hg38.fasta")


# In[7]:


bedfile = Bed.iterator(iotools.open_file(sys.argv[1]))
splice_site_dict = dict()
outfile = iotools.open_file(sys.argv[2], "w")
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

# In[17]:


from collections import defaultdict

counts = defaultdict(int)
for name, ss in splice_site_dict.items():
    counts[(ss[0].upper(), ss[1].upper())] += 1

sorted_keys = sorted(counts.keys(), key = lambda x: counts[(x)])    
for key in sorted_keys:
    print(key, counts[key])

