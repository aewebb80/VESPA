************
Installation
************

Dependencies
============

The VESPA software package is designed to minimize potential software dependencies, as additional software requirements may be difficult for users to install on their systems. Currently, the non-standard python library dendropy [Sukumaran et al., 2010] is the only dependency that remains in VESPA. Dendropy incorporates numerous functions for storing phylogenetic information and simplifying tree-based analyses. Dendropy can be installed from here (https://pythonhosted.org/DendroPy/). Removal of dendropy would require substantial development time and the design of numerous core functions. However, installation of dendropy is simple and only requires a single command to be invoked by the user. If the user invokes a dendropy-dependent function, VESPA is designed to print a warning message detailing the installation process of dendropy if the software is not installed.

1.12: Third-Party Programs
The VESPA software package is designed to interface with a number of third-party programs. It should be noted that some functions were designed to interface with specific versions of these third-party programs and future updates may require updates to VESPA as well. Details on these third-party programs can be found below.

Program
Version
Program Link
BLAST
2.2.30+
ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/
DendroPy
4.0
https://pythonhosted.org/DendroPy/#installing
MetAL
1.1
http://kumiho.smith.man.ac.uk/blog/whelanlab/?page_id=396
MrBayes
3.2.3
http://mrbayes.sourceforge.net/
MUSCLE
3.8.21
http://www.drive5.com/muscle/downloads.htm
NoRMD
1.3
ftp://ftp-igbmc.u-strasbg.fr/pub/NORMD/
PAML
4.4e
http://abacus.gene.ucl.ac.uk/software/paml.html
ProtTest3
3.4
https://code.google.com/p/prottest3/

