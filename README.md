# trichoCODE version 1.0e

1.  Installation
2.  Softwares
3.  Configuration file
4.  Input files
5.  Running trichoCODE
6.  Output files

================
1.  Installation
================
trichoCODE is a collection of Python libraries that do the job.
Unpack somewhere in your home directory.

trichoCODE was developed and tested on RPM-based Linux, but should be worked on any
*unix system.

=========
2.  Softwares
=========
trichoCODE leverage various third party softwares. Before run it, please install all of the 
necessary softwares. Please make sure the program packages with the indicated versions. All these programs should be accessible through your system PATH variable.

**Required softwares**
- Perl 5
- Augustus 2.7          (http://bioinf.uni-greifswald.de/augustus/binaries/)
- GeneMark-ET 4.21      (http://opal.biology.gatech.edu/GeneMark/)
- SNAP 2006-07-28       (http://korflab.ucdavis.edu/software.html)
- Exonerate 2.2.0       (http://www.ebi.ac.uk/~guy/exonerate/)
- Blat v.35             (http://hgdownload.cse.ucsc.edu/admin/exe/)
- TransposonPSI 2.2.26  (http://transposonpsi.sourceforge.net/)
- EVM r2012-06-25       (http://evidencemodeler.github.io/)
- PASA r20130605        (http://pasapipeline.github.io/)

**Required Python modules**
- Biopython 1.6         (http://biopython.org/)

======================
3.  Configuration file
======================
The values set in the file "config.cfg" are used as defaults. Any of these can be 
left empty and set on the command line using the same name. If a required value (indicated by
"required" on the comment line in ""config.cfg") is set neither in CONFIG nor on the 
command line, trichoCODE will complain and exit.

trichoCODE saves a CONFIG file with the values for each project in the project root 
directory. Besides providing a record, this file allows you to restart a project 
run without re-entering any parameters.

================
4.  Input files
================
The minimum required input files are: genome sequence, masked genome sequence and homology protein sequence. All of the sequence file should be fasta format.

If RNA-Seq data is available, user can supply trichoCODE with assembled transcript fasta file, intron hint in gff format.

======================
5.  Running trichoCODE
======================

config.cfg configuration:

To run trichoCODE, user need to copy the config.cfg to some folder and modified the values.
The config.cfg contains three blocks. The necessary values should be indicated, such as project_names,
genome, masked_genome and homology_protein, etc in [Project] block. MySQL database is necessary for PASA,
if user want to use the assembled transcript data. mysql_user_name and mysql_user_passwd are required.

For the results filtering, the default value for the maximum intron is 1000, and minimum protein sequence is 50.
User may change the threshold based their strong evidence, if not we suggest to keep the default for fungal genomes.

For the RNA-Seq, if user didn't have RNA-Seq data for their own organism, trichoCODE is capable dealing with related cross-species transcript alignments. Generally, as little as 70 to 80 percent nucleotide identity have proved useful by BROAD institute(Haas, 2011). In such case, filter_intron_max may be set around 70 and min_avg_per_id is around 75 in config.cfg.

parameters:
trichoCODE is designed for user-control annotation. In another word, user cna choose to run the whole pipeline in a time or run the process stage by stage. This is very useful when user first run the whole process, but not satisfied with the results. Then re-run the pipeline at the beginning. There are four stages user can control through the -s parameter and -m. -s means which stage to run. If without -m, it will start at some stage to the end. For example, 

$ trichoCODE.py -c config.cfg -s prediction

It means from the prediction stage to the end. But if user just want to run the prediction stage only. -m is wanted. Like below:

$ trichoCODE.py -c config.cfg -s prediction -m

In such case, the program will only run the prediction stage and stop.

Taking full advantage of this method, user may have the flexibility to obtain the best annotation.

An example data set containing all the input files needed to predict the genes 
on scaffold_6, 7 and 8 of Trichoderma reesei. They are available for download at 
http://somewhere.


================
6.  Output files
================
The annotation results are in project/products folder. the project is the name in config.cfg. In this folder, evm.gff3 is the output file of EVM. evm_filter.gff3 contains the filtering results based on the configuration in config.cfg. If assembled transcript were supplied, pasa_update.gff will be produced. It contained the updated gene models with UTRs, alternative spliced, exon corrections. Then all the fasta files were produced for genes, proteins, cDNA, CDS.

In project/reports folder, user will find all the statistics of each important file, which supply user the hints for good or bad annotation.

If something wrong with the running process. User can check the logs in project/log folder.



 
