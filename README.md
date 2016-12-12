# VESPA
VESPA: Very large-scale Evolutionary and Selective Pressure Analyses

Thanks for taking an interest in our pipeline. Everything you need to use this software is on this page. 
We hope you find the resources we provide here useful in getting you set up to analyse and interpret your data.

To reference VESPA:https://peerj.com/preprints/1895/

The VESPA software is released under the GNU general public license.


Installing VESPA:

VESPA1.0b is written in python and therefore should be compatible with most systems.

The VESPA software is available to download on GitHub from the following link:

Once downloaded, it can be installed as follows:

$ tar -xvzf VESPA.tar.gz

$ cd VESPA

$ chmod +x vespa.py

$ sudo mv vespa.py /usr/local/bin


Dependencies:

Perl dependencies:

VESPA requires users to install multiple perl scripts and modules to be fully operational. 
The perl dependencies may be found alongside vespa.py within the program tarball. 

Once downloaded, it can be installed as follows:

$ tar -xvzf VESPA.tar.gz

$ cd VESPA

$ chmod +x *Codeml*.pl

$ sudo mv *Codeml*.pl /usr/local/bin


$ sudo mv CodemlWrapper/ /Library/Perl/5.XX/
Note: Replace “5.XX” in the following command to the version of perl used by your system. (perl –v)

DendroPy: 

VESPA requires users to install the DendroPy python library (version 4.0). 
Instructions to installing DendroPy can be found at the following link: https://pythonhosted.org/DendroPy/#installing

 

Using VESPA:

We have written a manual that will take you through the process step by step. 
Please download the manual from here “Your_VESPA_Manual“. And follow the steps. 
If you have any difficulty please let us know. To get started we highly recommend 
you take our tutorial located here “VESPA_Tutorials.tar.gz“, which is specially
designed for the new user.
