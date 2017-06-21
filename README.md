# VESPA
## VESPA: Very large-scale Evolutionary and Selective Pressure Analyses

Thanks for taking an interest in our pipeline. Everything you need to use this software is on this page.
We hope you find the resources we provide here useful in getting you set up to analyse and interpret your data.

To reference VESPA: https://peerj.com/preprints/1895/

The VESPA software is released under the GNU General Public License.

### Installing VESPA:

VESPA 1.0b is written in Python and therefore should be compatible with most systems.

The VESPA software is available to download on GitHub from the following link:

Once downloaded, it can be installed as follows:

```
$ tar -xvzf VESPA.tar.gz

$ cd VESPA

$ chmod +x vespa.py

$ sudo mv vespa.py /usr/local/bin
```

### Dependencies:

#### Perl dependencies:

VESPA requires users to install multiple Perl scripts and modules to be fully operational.
The Perl dependencies may be found alongside vespa.py within the program tarball.

Once downloaded, it can be installed as follows:

```
$ tar -xvzf VESPA.tar.gz

$ cd VESPA

$ chmod +x \*Codeml\*.pl

$ sudo mv \*Codeml\*.pl /usr/local/bin

$ sudo mv CodemlWrapper/ /Library/Perl/5.XX/
```
*Note: Replace “5.XX” in the following command to the version of perl used by your system. (perl –v)*

#### Python dependencies:

VESPA requires users to install the DendroPy Python library (version 4.0).
Instructions to installing DendroPy can be found at the following link: https://pythonhosted.org/DendroPy/#installing

### Using VESPA:

We have written a manual that will take you through the process step by step.
Please download the manual from here "[Webb_VESPA_Manual][1]". And follow the steps.
If you have any difficulty please let us know. To get started we highly recommend
you take our tutorial located here "[VESPA_Tutorials.tar.gz][2]", which is specially
designed for the new user.

#### A little extra help:

If you would like a little extra help we have generated the following videos that take you through some of the more involved stages of the pipeline.
The videos are listed below (note our name for VESPA was Lamp when we made these videos):

1. Phase 1 tutorial video: [phase1][3]
2. Phase 2 tutorial video: [phase2][4]
3. Phase 3 tutorial video: [phase3][5]
4. Phase 4 tutorial video: [phase4][6]
5. Phase 5 tutorial video: [phase5][7]

[1]: http://mol-evol.org/wp-content/uploads/2015/09/Webb_VESPA_Manual.docx
[2]: http://mol-evol.org/wp-content/uploads/2015/12/VESPA_Tutorials.tar.gz
[3]: http://mol-evol.org/vespa/attachment/phase1
[4]: http://mol-evol.org/wp-content/uploads/2015/09/phase2.mp4
[5]: http://mol-evol.org/wp-content/uploads/2015/09/phase3.mp4
[6]: http://mol-evol.org/wp-content/uploads/2015/09/phase4.mov
[7]: http://mol-evol.org/wp-content/uploads/2015/09/phase5.mp4
