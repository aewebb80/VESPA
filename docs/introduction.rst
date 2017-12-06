*********************
Introduction to VESPA
*********************

Introduction
============

The VESPA (Very large-scale Evolution and Selective Pressure analyses) toolkit is a collection of commands designed to simplify molecular evolutionary analyses. The major motivation behind the development of VESPA was minimizing potential sources of error in large-scale selective pressure analyses using codeML from the PAML package [Yang 2007].

Assessing selective pressure variation on a large scale using protein coding DNA sequences requires a complex pipeline composed of numerous independent analyses, including: ortholog identification, multiple sequence alignment, phylogenetic reconstruction, and assessment of codon-based models of evolution. The pipeline requires multiple data manipulation steps to combine the output of each of these different stages (such as parsing BLAST result files to identify homologs and assessment of the suitability of the phylogenetic tree for selective pressure analysis). For researchers new to bioinformatics, manual data manipulation is prone to error, is potentially unstandardized, and is difficult to reproduce. VESPA eliminates the need for manual data manipulation by creating functions that automatically complete the majority of data manipulation steps using a standardized approach. In addition, the use of VESPA should minimize the requirements for users to create their own programs.

While the procedures within each stage of a selective pressure analysis are independent, there are requirements on the order in which the phases are carried out. VESPA creates a standardize pipeline of analyses with a specific ordering of phases in the process. In addition, the package encompasses multiple specialized pipelines to accurately assess selective pressure and reduce potential false positives (such as those caused by alignment error [Fletcher and Yang, 2010]).

Lastly, VESPA was designed to increase user productivity by automating labor intensive or highly repetitive tasks: i) automation by recursion – used to repeat an analysis on a number of files  (e.g. cleaning and translating a directory of genomes), and ii) automation of analysis methods – used to complete tasks that are normally demanding but invariable in execution (e.g. identifying homologs within BLAST output data). Automating these procedures within VESPA has created an analysis package that is highly scalable (i.e. from single gene to whole genome analyses) and that flexible enough to be useful for many levels of expertise and many alternative purposes/endpoints.


Phase and Pipeline Structure
=================================

The VESPA toolkit is separated into five separate analysis phases (Figure 1). The rationale behind the “phase” system was primarily to aid users in understanding the distinct procedures involved in selective pressure analysis and to provide more advanced users with a flexible and adaptable pipeline. Functions within a phase also analyze the same input type (e.g. sequences, BLAST output, etc.). 
The output of each phase in the VESPA toolit requires an analysis step that must be completed by the user with third-party software (Figure 1). These analyses are not automated by VESPA for three reasons: i) these analyses are far too computationally intensive, and ii) the submission process for these programs may differ from user to user, and iii) software updates may create bugs within the pipeline.  
The software package also incorporates two analysis pipelines, a basic pipeline for single gene orthologs (SGOs) and an advanced pipeline for both SGOs and multi-gene families (MGFs). The basic pipeline was designed to bypass the phylogenetic reconstruction techniques (phase 3) by inferring a gene phylogeny from a user-defined species phylogeny. Usage of the basic pipeline is only recommended if the genes are confirmed SGOs. See Figure 1 for additional details on the pipelines.

