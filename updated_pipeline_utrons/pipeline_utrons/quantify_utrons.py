'''
quantify_utrons.py
====================================================

:Author: Isabel
:Release: $1.0$
:Date: |today|
:Tags: Python

Purpose
-------

.. To use files produced by find_utrons.py and a expression db to produce a table 
with utron and quantification data for all transcripts.

Usage
-----

.. Example use case

Example::

   python cgat_script_template.py

Type::

   python quantify_utrons.py --help

for command line help.

Command line options
--------------------

'''

import sys
from CGATCore import Experiment as E
import os
import CGAT.GTF as GTF
import CGAT.Bed as Bed
import CGATCore.IOTools as IOTools
import CGAT.Intervals as Intervals
import itertools
import CGATCore.PipelineUtilities as PUtils
import numpy as np
import pandas as pd



def getStopDistdf(row):
    if row['strand'] == '+':
        return row['start'] - row['stop']
    else:
        return row['stop'] - row['end']
    
def getOverUnder50(row):
    if row > 50:
        return "Over"
    else:
        return "Under"

def isUtron(row):
    if row['over_under_50'] == 'Over' or row['over_under_50'] == 'Under':
        return 'Yes'
    else:
        return 'No'

def insertYesCol(row):
    return 'Yes'

def get_tcons(row):
    split = row['name'].split(':')
    return split[0]

def get_enst(row):
    split = row['name'].split(':')
    return split[1]

def get_tcons_from_ens(row):
    ens = row['partner_id']
    try:
        db1 = PUtils.fetch_DataFrame("SELECT transcript_id, transfrag_id FROM agg_agg_agg_cuffcompare_transcripts WHERE transcript_id='%s'" %ens, db)
        db1 = db1.set_index('transcript_id').drop_duplicates()
        tcons = db1.loc[ens]
        return tcons
    except KeyError:
        return "No_id"

def label_treatment(row):
    if "Pre" in row['track']:
        return 'Pre'
    else:
        return 'Post'
    
def label_patient(row):
    split = row['track'].split('-')
    return split[0] + split[2]



##########
def main(argv=None):
    """script main.
    parses command line options in sys.argv, unless *argv* is given.
    """

    if argv is None:
        argv = sys.argv

    # setup command line parser
    parser = E.OptionParser(version="%prog version: $1.0$",
                            usage=globals()["__doc__"])

    parser.add_option("-d", "--database", dest="database", type="string",
                      help="Supply database name")

    parser.add_option("-u", "--indivfile", dest="indivfile", type="string",
                      help="Supply input bed file name for individual utrons")

    parser.add_option("-p", "--partfile", dest="partfile", type="string",
                      help="Supply input bed file name for partnered utrons")

    parser.add_option("-n", "--novelfile", dest="novelfile", type="string",
                      help="Supply input bed file name for novel utrons")

    parser.add_option("-t", "--targetfile", dest="targetfile", type="string",
                      help="Supply input bed file name for miRNA TSs")

    parser.add_option("-o", "--outfile", dest="outfile", type="string",
                      help="Supply output csv file name")


    # add common options (-h/--help, ...) and parse command line
    (options, args) = E.Start(parser, argv=argv)

    global db
    db = options.database

#get expressions files
    expressions = PUtils.fetch_DataFrame("SELECT track, match_gene_id, transfrag_id, fpkm FROM agg_agg_agg_cuffcompare_transcripts CROSS JOIN agg_agg_agg_class WHERE transfrag_id = agg_agg_agg_class.transcript_id AND fpkm > 0", options.database)
    expressions = expressions.set_index(["track", "match_gene_id", "transfrag_id"])

    grouped_expression = expressions["fpkm"].groupby(level=["track","match_gene_id"])
    ex_fracts = grouped_expression.apply(lambda x: x/x.sum())
    ex_fracts.to_csv("pruned_expressionfractions.csv") 
    ex_sums = grouped_expression.apply(lambda x: x.sum())
    ex_sums.to_csv("pruned_expressionsums.csv") 

    ex_sums = pd.read_csv("pruned_expressionsums.csv", names=['track','match_gene_id','exp_sum'])
    ex_sums = ex_sums.set_index(['match_gene_id','track'])
    ex_fracts = pd.read_csv("pruned_expressionfractions.csv", names=['track','match_gene_id','transfrag_id','exp_fract'])
    ex_fracts = ex_fracts.set_index(['track','match_gene_id','transfrag_id'])
    fpkm_ex_fracts = ex_fracts.join(expressions, how='inner')
    fpkm_ex_fracts = fpkm_ex_fracts.reset_index()
    fpkm_ex_fracts = fpkm_ex_fracts.set_index(['match_gene_id','track'])
    ex_all = fpkm_ex_fracts.join(ex_sums, how='inner')
    ex_all = ex_all.reset_index()
    ex_all.to_csv("pruned_expression_all.csv")  
    ex_all = pd.read_csv("pruned_expression_all.csv")
    ex_all = ex_all.set_index('transfrag_id')

#stop distances
    ind_utrons = pd.read_table(options.indivfile, header=0, sep='\t', names=["chrom","start","end","name","score","strand","stop"], usecols=["start","end","name","strand","stop"], compression='gzip')

    ind_utrons['dist'] = ind_utrons.apply(lambda row: getStopDistdf(row), axis=1)
    ind_utrons = ind_utrons.set_index('name')
    grouped_stopdist = ind_utrons.groupby(level='name')
    transcript_dist = grouped_stopdist.apply(lambda group: group['dist'].max())
    transcript_dist.name = 'dist'
    transcript_over_under_50 = transcript_dist.apply(lambda row: getOverUnder50(row))
    transcript_over_under_50.name = 'over_under_50'

    ex_all_dist = ex_all.join(transcript_over_under_50, how='left')
    ex_all_dist = ex_all_dist.join(transcript_dist, how='left')
    ex_all_dist['utron'] = ex_all_dist.apply(lambda row: isUtron(row), axis = 1)

#novel utrons

    novel_utrons = pd.read_table(options.novelfile, header=0, sep='\t', names=["chrom","start","end","name","score","strand","a","b","c","d","e","f"], usecols=["start","end","name"], compression='gzip')
    novel_utrons = novel_utrons.set_index(novel_utrons["name"])
    novel_utrons = novel_utrons.drop_duplicates(subset="name")  #excludes entries with different start/end utron coordinates in same transcript

    novel_utrons['novel_utron'] = novel_utrons.apply(lambda row: insertYesCol(row), axis = 1)
    novel_utrons = novel_utrons.drop(['start','end','name'], axis=1)
    ex_all_dist_nov = ex_all_dist.join(novel_utrons, how='left')

#TSs

    utron_TSs = pd.read_table(options.targetfile, header=0, sep='\t', names=["chrom","start","end","name","score","strand","stop"], usecols=["start","end","name","strand","stop"], compression='gzip')
    utron_TSs['miRNA_TS'] = utron_TSs.apply(lambda row: insertYesCol(row), axis=1)
    utron_TSs = utron_TSs.drop(["start","end","strand","stop"],axis=1).drop_duplicates()
    utron_TSs = utron_TSs.set_index(["name"])
    ex_all_dist_nov_TS = ex_all_dist_nov.join(utron_TSs, how='left')

#extra utrons

    tcons_ens = pd.read_table(options.partfile, header=0, sep='\t', names=["chrom","start","end","name","score","strand","a",'b','c','d','e','f'], usecols=["start","end","name","strand"], compression='gzip')
    tcons_ens['TCONS_id'] = tcons_ens.apply(lambda row: get_tcons(row), axis=1)
    tcons_ens['partner_id'] = tcons_ens.apply(lambda row: get_enst(row), axis=1)
    tcons_ens = tcons_ens.set_index('TCONS_id')

    tcons_ens['partner_id_TCONS'] = tcons_ens.apply(lambda row: get_tcons_from_ens(row), axis=1)
    tcons_ens = tcons_ens.drop_duplicates()
    tcons_ens['extra_utron'] = tcons_ens.apply(lambda row: insertYesCol(row), axis=1)

    partners = tcons_ens[['name','partner_id_TCONS']]
    partners = partners[partners['partner_id_TCONS']!='No_id']
    partners = partners.set_index('partner_id_TCONS')
    utrons_and_partners = tcons_ens.append(partners)
    utrons_and_partners = utrons_and_partners.join(ex_all_dist_nov, how='inner')
    utrons_and_partners = utrons_and_partners.reset_index().drop_duplicates(subset=['match_gene_id','track','index'])
    utrons_and_partners = utrons_and_partners.set_index(['match_gene_id','track'])
    groups = utrons_and_partners.groupby(level=['match_gene_id','track'])
    sums=groups.apply(lambda group: sum(group['fpkm']))
    utrons_and_partners['partner_exp_sum'] = sums
    utrons_and_partners['partner_exp_fract']=utrons_and_partners.apply(lambda row: row['fpkm']/row['partner_exp_sum'], axis=1)
    only_utrons = utrons_and_partners[utrons_and_partners['extra_utron']=='Yes']
    only_utrons = only_utrons[['index','extra_utron', 'partner_exp_sum', 'partner_exp_fract', 'partner_id_TCONS', 'partner_id']]
    only_utrons = only_utrons.reset_index()
    only_utrons = only_utrons.dropna(subset=['match_gene_id','track','index'])
    only_utrons = only_utrons.set_index(['match_gene_id','track','index'])

    ex_all_dist_nov_TS = ex_all_dist_nov_TS.reset_index()
    ex_all_dist_nov_TS = ex_all_dist_nov_TS.set_index(['match_gene_id','track','index'])
    ex_all_dist_nov_TS_ext = ex_all_dist_nov_TS.join(only_utrons, how='left')

#patients and treatment
    final = ex_all_dist_nov_TS_ext.reset_index()
    final['treatment'] = final.apply(lambda row: label_treatment(row),axis=1)
    final['patient'] = final.apply(lambda row: label_patient(row), axis=1)


    final.to_csv(options.outfile)

    
    # write footer and output benchmark information.
    E.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))

