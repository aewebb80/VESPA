*************************
Phase 2: homology search
*************************

Command structure
=================

The VESPA software package was written in python (v2.7) and requires a UNIX environment to operate. VESPA may be invoked as follows: 
usr$ python vespa.py
The VESPA help screen will then be displayed by default. If desired, the help screen may also be displayed using the following commands.
usr$ python vespa.py help
usr$ python vespa.py h
In addition to the basic help screen, VESPA has the option to display basic help information for each VESPA command. If desired, the help information may be displayed by specifying the command of interest subsequent to the help screen call (please note the space):
usr$ python vespa.py help translate
Commands in VESPA are specified after the program call (i.e. python vespa.py) on the UNIX command-line. Please note a space is required between the program call and the desired command. For example, the translate command would be invoked as shown below:
usr$ python vespa.py translate
Commands also require specific options to be invoked to function correctly. Options are specified after the command and begin with a dash symbol (-) and end with an equal sign (=) followed by either a user-specified file or Boolean value (i.e. True/False). For example, the translate command requires the user to specify the input (here “user_data.txt”) as follows:
usr$ python vespa.py translate -input=user_data.txt
Please note the space between the command and option, it should also be noted that there is no space separating the option (i.e. “-input=”) and the user-specification (i.e. “user_data.txt”). Multiple options may be invoked on the same command-line as shown below and are separated by a space:
usr$ python vespa.py translate -input=user_data.txt -cleave_terminal=False
A comprehensive list of the commands supported in VESPA may be found on Pg. 10 of this manual. 
