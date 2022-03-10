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
Pipeline rMATs Simple
===========================

:Author: Andreas Heger
:Release: $Id$
:Date: |today|
:Tags: Python




Code
====

"""
from ruffus import *
import sys
import os
import shutil
import sqlite3
from glob import glob
from cgatcore import experiment as E
import cgat.Sra as Sra
from cgatcore import pipeline as P
import re

# load options from the config file
PARAMS = P.get_parameters(
    ["%s/pipeline.yml" % os.path.splitext(__file__)[0],
     "../pipeline.yml",
     "pipeline.yml"])

# -----------------------------------------------
# Utility functions
def connect():
    '''utility function to connect to database.

    Use this method to connect to the pipeline database.
    Additional databases can be attached here as well.
    Returns an sqlite3 database handle.
    '''

    dbh = sqlite3.connect(PARAMS["database_name"])

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
    

# ---------------------------------------------------
@transform(["*.bam",
           "*.remote"],
           formatter(),
           add_inputs(PARAMS["geneset"]),
           r"rmats.dir/prep_files.dir/{basename[0]}.rmats",
           r"{basename[0]}")
def run_rmats_pre(infiles, outfile, track):

    infile, gtffile = infiles
    
    od = os.path.abspath("rmats.dir")
    
    statement = '''rmats.py
                   --task prep
                   --tmp rmats.dir/%(track)s.dir
                   --gtf <(zcat %(gtffile)s)
                   --readLength %(rmats_readLength)s
                   -t %(rmats_paired)s
                   --od rmats.dir
                   --b1 <(echo %(infile)s)
                   &> rmats.dir/%(track)s.prep.log
                   '''                          
                   
    if infile.endswith(".remote"):
        token = glob("gdc-user-token*")
        tmpfilename = P.get_temp_filename()
        if os.path.exists(tmpfilename):
            os.unlink(tmpfilename)    	
        if len(token) > 0:
            token = token[0]
        else:
            token = None

        s, infile = Sra.process_remote_BAM(
            infile, token, tmpfilename,
            filter_bed=PARAMS["contigs_bed"])
        s = re.sub(";\n", " &&\n", s)

        infile = ",".join(infile)
        statement = " && ".join(
            ["mkdir -p %(tmpfilename)s",
             s,
             statement,
             "rm -r %(tmpfilename)s"])

    P.run(statement, 
          job_condaenv=PARAMS["rmats_env"],
          job_memory = PARAMS["rmats_prep_memory"])

    rmats_counter = ""    
    for f_path in glob("rmats.dir/%(track)s.dir/*.rmats" % locals()):
        shutil.copy(f_path, 
                    P.snip(outfile, ".rmats") + rmats_counter + ".rmats")
        if rmats_counter == "":
            rmats_counter = 1
        else:
            rmats_counter += 1
            
P.main(sys.argv)