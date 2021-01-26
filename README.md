# pipeline_utrons

## The pipeline performs the following:
   * Uses sequence alignment data from `.bam` files to assemble transcripts with StringTie.
   * Analyses and quantifies splice junctions in .bam files, filtering those which are 
   unlinkely to be genuine.
   * Merges all assembled transcripts with StringTie.
   * Classifies and load transcripts into a databases.
   * Uses a separately custom-built pipeline to find utrons and their IDs.
   * Makes a Salmon index and quantifies transcripts with Salmon.
   * Merges all quantification files and uploads the outputs onto database.
   
## Inputs needed          
1. RNA-se read alignments in a .bam format, in order to generate assembled transcripts in .gtf format.
2. A transcriptome in .fa format, in order to build a salmon index.
3. RNA-seq reads in fastq.1.gz and fastq.2.gz formats, for generating salmon quantification files in .sf format.

The updated `GenomeAnalysis.pyx` and `IndexedGenome.py` should replace the old versions in `cgat/cgat-apps/CGAT`.

## Configuration
The pipeline requires a configured :file: `pipeline.yml` file.

Make a directory with your project name.
Configure the pipeline with `python [path_to_repo]/pipeline_utrons.py config`.
A pipeline.log and pipeline.yml file(s) will be added to your new directory.
Modify the pipeline.yml according to your project (specify annotation database and directory, database for uploading the outputs; specify options for Salmon quantification).

## Pipeline use
Run the pipeline with `python [path_to_repo]/pipeline_utrons.py config make full`.

For running the pipeline on a large set of samples, submit the pipeline onto the cluster (sharc), using a submit_pipeline_cgtaflow custom script.



