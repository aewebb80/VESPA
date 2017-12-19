************
Installation
************

::

    $ curl -o VESPA.tar.gz https://github.com/aewebb80/VESPA/archive/1.0.0.tar.gz
    $ tar -xzf VESPA.tar.gz
    $ cd VESPA
    $ chmod +x vespa.py
    $ sudo mv vespa.py /usr/local/bin

**Perl Dependencies**

VESPA requires users to install multiple Perl scripts and modules to be fully operational. These may be found alongside vespa.py within the program tarball.

Once downloaded, it can be installed as follows:
::

    $ tar -xvzf VESPA.tar.gz
    $ cd VESPA
    $ chmod +x *Codeml*.pl
    $ sudo mv *Codeml*.pl /usr/local/bin
    $ sudo mv CodemlWrapper/ /Library/Perl/5.XX/`

.. note:: Replace :code:`5.XX` with the version of Perl used by your system. (determined by executing :code:`$ perl –v`)

**DendroPy**

VESPA requires users to install the DendroPy python library (version 4.0). Instructions to installing DendroPy can be found at the following link: https://pythonhosted.org/DendroPy/#installing




Third party software
====================

The VESPA software package is designed to interface with a number of third-party programs. It should be noted that some functions were designed to interface with specific versions of these third-party programs and future updates may require updates to VESPA as well. Details on these third-party programs can be found below.

+-----------+---------+-------------------------------------------------------------------+
| Program   | Version | URL                                                               |
+===========+=========+===================================================================+
| BLAST     | 2.2.30\+| `<ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST>`_   |
+-----------+---------+-------------------------------------------------------------------+
| DendroPy  | 4.0     | `<https://pythonhosted.org/DendroPy/#installing>`_                |
+-----------+---------+-------------------------------------------------------------------+
| MetAL     | 1.1     | `<http://kumiho.smith.man.ac.uk/blog/whelanlab/?page_id=396>`_    |
+-----------+---------+-------------------------------------------------------------------+
| MrBayes   | 3.2.3   | `<http://mrbayes.sourceforge.net/>`_                              |
+-----------+---------+-------------------------------------------------------------------+
| MUSCLE    | 3.8.21  | `<http://www.drive5.com/muscle/downloads.htm>`_                   |
+-----------+---------+-------------------------------------------------------------------+
| NoRMD     | 1.3     | `<ftp://ftp-igbmc.u-strasbg.fr/pub/NORMD/>`_                      |
+-----------+---------+-------------------------------------------------------------------+
| PAML      | 4.4e    | `<http://abacus.gene.ucl.ac.uk/software/paml.html>`_              |
+-----------+---------+-------------------------------------------------------------------+
| ProtTest3 | 3.4     | `<https://github.com/ddarriba/prottest3>`_                        |
+-----------+---------+-------------------------------------------------------------------+