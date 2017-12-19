# VESPA: Very large-scale Evolutionary and Selective Pressure Analyses
Thanks for taking an interest in our pipeline. We hope you find the resources we provide here useful in getting you set up to analyse and interpret your data.



To reference VESPA: https://peerj.com/preprints/1895/

Documentation is now hosted on [ReadTheDocs](http://vespa-evol.readthedocs.io/en/latest/)



## Installation

The main VESPA wrapper is written in Python, while the CodeML wrapper is written in Perl.

To begin, we recommend downloading a versioned tarball from the releases page.

Once downloaded, install as follows:

```
$ tar -xvzf VESPA.tar.gz
$ cd VESPA
$ chmod +x vespa.py
$ sudo mv vespa.py /usr/local/bin
```




### Dependencies

**Perl dependencies**

VESPA requires several Perl scripts and modules to be fully operational. 
The Perl dependencies may be found alongside `vespa.py` within the tarball. 

The following will install these modules on most systems.

```
$ tar -xvzf VESPA.tar.gz
$ cd VESPA
$ chmod +x \*Codeml\*.pl
$ sudo mv \*Codeml\*.pl /usr/local/bin
$ sudo mv CodemlWrapper/ /Library/Perl/5.XX/
Note: Replace “5.XX” in the following command to the version of perl used by your system. (perl –v)
```



**DendroPy** 

VESPA requires users to install the DendroPy Python library (version 4.0). 
[Instructions for installing DendroPy](https://pythonhosted.org/DendroPy/#installing)

 

## Using VESPA

We have written [a manual](http://vespa-evol.readthedocs.io/en/latest/) that will take you through the process step by step. 
This is constantly evolving with feedback with feedback from our users, and the latest version is hosted as [vespa-evol on ReadTheDocs](http://vespa-evol.readthedocs.io/en/latest/).



## Issues

If you encounter an issue with VESPA, please open a GitHub issue and we'll do our best to help.