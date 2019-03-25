#!/usr/bin/env python
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### __      ________  _____ _____        
### \ \    / /  ____|/ ____|  __ \ /\    
###  \ \  / /| |__  | (___ | |__) /  \   
###   \ \/ / |  __|  \___ \|  ___/ /\ \  
###    \  /  | |____ ____) | |  / ____ \ 
###     \/   |______|_____/|_| /_/    \_\
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### Authors: Andrew E. Webb, Thomas A. Walsh  & Mary J. O'Connell
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
def help_message(help_requested):
    global bme_command_table
    if help_requested in bme_command_table:
        third_party = False
        help_seperator = '----------------------------------------------------------------------------------'
        command_str = 'Command: {0}'.format(help_requested)
        command_text = '|||||{0}|||||'.format(command_str.center(72, ' '))
        print '\n{0}\n{1}\n{2}'.format(help_seperator, command_text, help_seperator)
        
        if help_requested == 'clean':
            print '''Details: QC filter for downloaded nucleotide sequences and/or genomes.
Basic usage: vespa.py clean -input=USR_INPUT
Supported file format(s): -input option: fasta formatted files'''
        
        elif help_requested == 'ensembl_clean':
            print'''Details: QC filter for identifying the longest nucleotide (canonical) transcript
         within an Ensembl nucleotide genome.
Basic usage: vespa.py ensembl_clean -input=USR_INPUT
Supported file format(s): -input option: fasta formatted files'''

        elif help_requested == 'translate':
            print'''Details: Translates nucleotide sequences that passed the QC filter of either clean
         function into amino acid sequences
Basic usage: vespa.py translate -input=USR_INPUT
Supported file format(s): -input option: fasta formatted files'''
        
        elif help_requested == 'create_database':
            print'''Details: Concatenates multiple genomes into the single database file.
Basic usage: vespa.py create_database -input=USR_INPUT
Supported file format(s): -input option: fasta formatted files'''
        
        elif help_requested == 'gene_selection':
            print '''Details: Searches a sequence database for gene identifiers specified within a
         separate csv file.
Basic usage: vespa.py gene_selection -input=USR_INPUT -selection_csv=USR_INPUT
Supported file format(s): -input option: fasta formatted files,
                          -selection_csv option: csv formatted files'''

        elif help_requested == 'similarity_groups':
            print '''Details: Construct sequence similarity groups with either non-reciprocal and
            reciprocal connections.
Basic usage: vespa.py similarity_groups -input=USR_INPUT -format=blast -database=USR_DB
Supported file format(s): -input option: BLAST tabular and HMMER output files
                          -database option: fasta formatted files
                          -format option: blast or hmmer'''

        elif help_requested == 'reciprocal_groups':
            print '''Details: Construct sequence similarity groups with only reciprocal connections.
Basic usage: vespa.py reciprocal_groups -input=USR_INPUT -format=blast -database=USR_DB
Supported file format(s): -input option: BLAST tabular and HMMER output files
                          -database option: fasta formatted files
                          -format option: blast or hmmer'''

        elif help_requested == 'best_reciprocal_groups':
            print '''Details: Construct sequence similarity groups with only reciprocal connections
         that share the best E-value for each species.
Basic usage: vespa.py best_reciprocal_groups -input=USR_INPUT -format=blast -database=USR_DB
Supported file format(s): -input option: BLAST tabular and HMMER output files
                          -database option: fasta formatted files
                          -format option: blast or hmmer'''

        elif help_requested == 'metal_compare':
            third_party = True
            print '''Details: Automates MSA comparison, scoring, and selection using the third-party
         programs MetAl and noRMD.
Basic usage: vespa.py metal_compare -input=USR_INPUT -compare=USR_INPUT
Supported file format(s): -input and -compare options: fasta formatted files'''

        elif help_requested == 'prottest_setup':
            third_party = True
            print '''Details: Automates the identification the best-fit model of amino acid replacement
         for a protein MSAs using the third-party program ProtTest3.
Basic usage: vespa.py prottest_setup -input=USR_INPUT
Supported file format(s): -input option: fasta formatted files'''

        elif help_requested == 'prottest_reader':
            third_party = True
            print '''Details: Automates the process of reading the output of the third-party program
         ProtTest3.
Basic usage: vespa.py prottest_reader -input=USR_INPUT
Supported file format(s): -input option: ProtTest3 output format'''

        elif help_requested == 'mrbayes_setup':
            third_party = True
            print '''Details: Simplifies phylogenetic reconstruction using the third-party program MrBayes
         by creating NEXUS formatted files with MrBayes command blocks.
Basic usage: vespa.py mrbayes_setup -input=USR_INPUT -model_list=MODEL_DATA
Supported file format(s): -input option: fasta formatted files
                          -model_list prottest_reader supported_output files'''

        elif help_requested == 'map_alignments':
            print '''Details: Automates the conversion of protein MSAs to nucleotide (codon) MSAs required
         for codeML.
Basic usage: vespa.py map_alignments -input=USR_INPUT -database=USR_DB
Supported file format(s): -input and -database options: fasta formatted files'''

        elif help_requested == 'infer_genetree':
            print '''Details: Automates the creation of the corresponding gene tree for a MSA using a
         user-specified species tree.
Basic usage: vespa.py infer_genetree -input=USR_INPUT -species_tree=USR_INPUT
Supported file format(s): -input option: fasta formatted files
                          -species_tree option: newick formatted files'''

        elif help_requested == 'codeml_setup':
            print '''Details: Automates the creation of the complex codeML directory structure.
Basic usage: vespa.py setup_codeml -input=USR_INPUT
Supported file format(s): -input option: fasta formatted files with corresponding
                                         newick formatted files'''

        elif help_requested == 'mrbayes_reader':
            print '''Details: Automates the conversion of nexus-formatted phylogenies into the
         newick format.
Basic usage: vespa.py mrbayes_reader -input=USR_INPUT
Supported file format(s): -input option: MrBayes converged NEXUS output files'''

        elif help_requested == 'create_subtrees':
            print '''Details: Simplifies pruning large multigene phylogenies into smaller
         sub-phylogenies.
Basic usage: vespa.py create_subtrees -input=USR_INPUT
Supported file format(s): -input option: newick formatted files'''

        elif help_requested == 'create_branch':
            print '''Details: Simplify the creation of the branch-label table required for
         the branch-site models of codeML.
Basic usage: vespa.py create_branch -input=USR_INPUT
Supported file format(s): -input option: newick formatted files'''

        elif help_requested == 'codeml_reader':
            print '''Details: Parses the complex codeML directory structure and create
         simplified results.
Basic usage: vespa.py codeml_reader -input=USR_INPUT
Supported file format(s): -input option: VESPA formmated codeML output files'''    
        
        if third_party:
            print '\nSee program manual for third-party program citations and additional options\n'
        else:
            print '\nSee program manual for additional options\n'
    elif help_requested:
        print 'Command not found, please confirm commands in help file.\n'
    else: 
        print '''VESPA v1.0b - [La]rge-scale [M]olecular evolution and selective pressure [P]ipeline
Authors: Andrew E. Webb, Thomas A. Walsh & Mary J. O'Connell
__      ________  _____ _____        
\ \    / /  ____|/ ____|  __ \ /\    
 \ \  / /| |__  | (___ | |__) /  \   
  \ \/ / |  __|  \___ \|  ___/ /\ \  
   \  /  | |____ ____) | |  / ____ \ 
    \/   |______|_____/|_| /_/    \_\

----------------------------------------------------------------------------------
|||||                              Command Help                             ||||||
----------------------------------------------------------------------------------
Specify the command of interest after invoking help/h.
For example: vespa.py help clean
----------------------------------------------------------------------------------

----------------------------------------------------------------------------------
|||||                        Phase 1: Data Preparation                       |||||
----------------------------------------------------------------------------------
This phase included for users new to bioinformatics. The phase prepares downloaded
genomes for homology searching using the two VESPA supported homology search tools:
BLAST and HMMER.

Commands: clean, ensemble_clean, translate, create_database, gene_selection
----------------------------------------------------------------------------------

----------------------------------------------------------------------------------
|||||                       Phase 2: Homology Searching                      |||||
----------------------------------------------------------------------------------
Details: This phase is concerned with identifying groups of similar sequences from
either BLAST or HMMER homology searches.

Commands: similarity_groups, reciprocal_groups, best_reciprocal_groups
----------------------------------------------------------------------------------

----------------------------------------------------------------------------------
|||||        Phase 3: Alignment Assessment & Phylogeny Reconstruction        |||||
----------------------------------------------------------------------------------
Details: This phase combines multiple third-party programs to automate the
assessment, selection, and phylogenetic reconstruction of protein MSAs.

Commands: metal_compare, prottest_setup, prottest_reader, mrbayes_setup
----------------------------------------------------------------------------------

----------------------------------------------------------------------------------
|||||                 Phase 4: Selection Analysis Preparation                |||||
----------------------------------------------------------------------------------
Details: This phase automates large-scale selective pressure analysis using codeML
from the PAML package.

Commands: map_alignments, infer_genetree, mrbayes_reader, link_input, codeml_setup,
          create_subtrees, create_branch
----------------------------------------------------------------------------------

----------------------------------------------------------------------------------
|||||                 Phase 5: Selection Analysis Assessment                 |||||
----------------------------------------------------------------------------------
Details: This phase automatically parses the codeML directory structure
and create simplified summary files 

Command: codeml_reader
----------------------------------------------------------------------------------

See program manual for additional information.
'''


### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
###   _____ _                         
###  / ____| |                        
### | |    | | __ _ ___ ___  ___  ___ 
### | |    | |/ _` / __/ __|/ _ \/ __|
### | |____| | (_| \__ \__ \  __/\__ \
###  \_____|_|\__,_|___/___/\___||___/                                
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### class: sequence - defines base sequence information and has basic functions                                                 ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
class sequence_data(object):
    def __init__(self, header, sequence):
        self.header = header
        self.sequence = sequence.upper()
        self.filename = header
        self.length = len(''.join(sequence))
        self.type = ''
    
    def __len__(self):
        return len(''.join(self.sequence))
    
    def __str__(self):
        prtSeq = ''.join(self.header) + '\n'.join([self.sequence[i:i+60] for i in range(0, len(self.sequence), 60)])
        return prtSeq
    
    def seq_filename (self, definement):
        self.filename = str(definement)
    
    def seq_translate (self):
        codonTbl = {'ATT':'I','ATC':'I','ATA':'I','CTT':'L',
                    'CTC':'L','CTA':'L','CTG':'L','TTA':'L',
                    'TTG':'L','GTT':'V','GTC':'V','GTA':'V',
                    'GTG':'V','TTT':'F','TTC':'F','ATG':'M',
                    'TGT':'C','TGC':'C','GCT':'A','GCC':'A',
                    'GCA':'A','GCG':'A','GGT':'G','GGC':'G',
                    'GGA':'G','GGG':'G','CCT':'P','CCC':'P',
                    'CCA':'P','CCG':'P','ACT':'T','ACC':'T',
                    'ACA':'T','ACG':'T','TCT':'S','TCC':'S',
                    'TCA':'S','TCG':'S','AGT':'S','AGC':'S',
                    'TAT':'Y','TAC':'Y','TGG':'W','CAA':'Q',
                    'CAG':'Q','AAT':'N','AAC':'N','CAT':'H',
                    'CAC':'H','GAA':'E','GAG':'E','GAT':'D',
                    'GAC':'D','AAA':'K','AAG':'K','CGT':'R',
                    'CGC':'R','CGA':'R','CGG':'R','AGA':'R',
                    'AGG':'R','TAA':'*','TAG':'*','TGA':'*'}
        strSeq = ''.join(self.sequence).strip()
        self.sequence =  ''.join([codonTbl.get(strSeq[3*n:3*n+3], 'X') for n in range(len(strSeq)//3)])
        self.length = len(self.sequence)
        return self
    
    
    def seq_revcomp (self):
        from string import maketrans
        strRev = ''.join(self.sequence).strip()[::-1]
        self.sequence = strRev.translate(maketrans('ACTGUBVDHKMRYNSW','TGACAVBHDMKYRNSW'))
        return self
    
    def seq_type (self):
        type_status = False 
        for unique_amino_characters in ['q', 'e', 'i', 'l', 'f', 'p']:
            if unique_amino_characters in ''.join(self.sequence).strip().lower():
                self.type = 'protein'
                type_status = True
        if not type_status:
            self.type = 'DNA'
        
    def internal_stop(self):
        def internal_translate(strSeq):
            codonTbl = {'ATT':'I','ATC':'I','ATA':'I','CTT':'L',
                        'CTC':'L','CTA':'L','CTG':'L','TTA':'L',
                        'TTG':'L','GTT':'V','GTC':'V','GTA':'V',
                        'GTG':'V','TTT':'F','TTC':'F','ATG':'M',
                        'TGT':'C','TGC':'C','GCT':'A','GCC':'A',
                        'GCA':'A','GCG':'A','GGT':'G','GGC':'G',
                        'GGA':'G','GGG':'G','CCT':'P','CCC':'P',
                        'CCA':'P','CCG':'P','ACT':'T','ACC':'T',
                        'ACA':'T','ACG':'T','TCT':'S','TCC':'S',
                        'TCA':'S','TCG':'S','AGT':'S','AGC':'S',
                        'TAT':'Y','TAC':'Y','TGG':'W','CAA':'Q',
                        'CAG':'Q','AAT':'N','AAC':'N','CAT':'H',
                        'CAC':'H','GAA':'E','GAG':'E','GAT':'D',
                        'GAC':'D','AAA':'K','AAG':'K','CGT':'R',
                        'CGC':'R','CGA':'R','CGG':'R','AGA':'R',
                        'AGG':'R','TAA':'*','TAG':'*','TGA':'*'} 
            return ''.join([codonTbl.get(strSeq[3*n:3*n+3], 'X') for n in range(len(strSeq)//3)])
        test_sequence = ''
        if self.type == 'unknown':
            self.seq_type()
        if self.type == 'DNA':
            test_sequence = internal_translate(''.join(self.sequence).strip())
        else:
            test_sequence = ''.join(self.sequence).strip()
        if not test_sequence.count('*'):
            return False
        else:
            if test_sequence.count('*') > 1:
                return True
            elif not test_sequence.endswith('*'):
                return True
            else:
                return False

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### class: sequence_reader - sequence reading capabilties                                                                       ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
class sequence_reader(object):
    def __init__(self, filename):
        self.filename = filename
        self.type = ''
        
    def sequence_format (self):
        with open(self.filename, 'rU') as unknown_sequence_data:
            format_data = unknown_sequence_data.readline().strip()
            
        list_format_data = format_data.split()
        if len(list_format_data) > 1:
            chk_num_format_data = [current_format_data.isdigit() for current_format_data in list_format_data]

        if format_data.startswith('>'):
            self.type = 'fasta'
        elif '#NEXUS' in format_data:
            self.type = 'nexus'
        elif len(list_format_data) == 2 and False not in chk_num_format_data:
            self.type = 'phylip'
        else:
            self.type = ''
    
    
    def read (self):
        from collections import defaultdict
        import sys
        sequence_reader.sequence_format(self)
        
        if self.type == 'fasta':
            def sequnce_reader (sequence_file):
                from itertools import groupby
                sequenceGroups = (entry[1] for entry in groupby(open(sequence_file, 'rU'), lambda line: line.startswith('>')))
                for header in sequenceGroups:
                    header = header.next()
                    seq = "".join(sequence.strip() for sequence in sequenceGroups.next())
                    yield header, seq
            for read_header, read_sequence in sequnce_reader(self.filename):
                yield sequence_data(read_header, read_sequence)
        
        elif self.type == 'phylip':
            with open(self.filename, 'rU') as phylip_sequence_data:
                phylip_stats = []
                phylip_dict = defaultdict(list)
                sequence_count = 0
                for phylip_sequence_lines in phylip_sequence_data:
                    if not phylip_stats:
                        sequence_stats = phylip_sequence_lines.split()
                        phylip_stats = [int(current_sequence_stat) for current_sequence_stat in sequence_stats]
                    elif not phylip_sequence_lines.strip():
                        sequence_count = 0
                    else:
                        if len(phylip_sequence_lines.split()) == 2:
                            phylip_split = phylip_sequence_lines.strip().split()
                            phylip_dict[sequence_count] = [phylip_split[0], phylip_split[1]]
                        else:
                            phylip_dict[sequence_count][1] += phylip_sequence_lines.strip()
                        sequence_count += 1   
            for phylip_key, (phylip_header, phylip_sequence) in phylip_dict.items():
                sequence_temp = sequence_data('>{0}\n'.format(phylip_header), phylip_sequence)
                del phylip_dict[phylip_key]
                yield sequence_temp
        
        elif self.type == 'nexus':
            with open(self.filename, 'rU') as nexus_sequence_data:
                nexus_stats = []
                nexus_dict = defaultdict(str)
                nexus_data_chk, nexus_matrix_chk = False, False
                for nexus_sequence_lines in nexus_sequence_data:
                    if nexus_data_chk and nexus_matrix_chk:
                        if nexus_sequence_lines.strip() == ';':
                            break
                        elif nexus_sequence_lines.strip() != '':
                            nexus_split = nexus_sequence_lines.strip().split()
                            print nexus_split
                            if not nexus_dict.has_key(nexus_split[0]):
                                nexus_dict[nexus_split[0]] = nexus_split[1]
                            else:
                                nexus_dict[nexus_split[0]] += nexus_split[1]
                    if 'Begin data;' in nexus_sequence_lines:
                        nexus_data_chk = True
                    if 'Matrix' in nexus_sequence_lines:
                        nexus_matrix_chk = True
                    if nexus_data_chk and not nexus_matrix_chk:
                        if 'ntax' in nexus_sequence_lines:
                            nexus_stats = [int(nexus_stat.split('=')[1]) for nexus_stat in nexus_sequence_lines.strip().replace(';','').split() if '=' in nexus_stat]
            for nexus_header, nexus_sequence in nexus_dict.items():
                sequence_temp = sequence_data('>{0}\n'.format(nexus_header), nexus_sequence)
                del nexus_dict[nexus_header]
                yield sequence_temp
        else:
            print 'Sequence format not recognized. Please confirm that input is correctly formmated.'
            sys.exit()

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### class: command_line_data - assigns multiple varibles defaults and allows for user defined options                           ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
class command_line_data(object):
    def __init__(self, define_location_string, define_if_input_in_dir):
        #assignment variables
        self.input_in_dir = define_if_input_in_dir
        self.dir_location = ''
        if self.input_in_dir:
            self.dir_location = define_location_string.rsplit('/',1)[0]
        self.input_location = define_location_string
        
        #current variables
        self.current_input = define_location_string
        self.current_input_filename = define_location_string
        self.current_input_dir = ''
        if self.input_in_dir:
            self.current_input_filename = define_location_string.rsplit('/',1)[-1]
            self.current_input_dir = define_location_string.rsplit('/',1)[0]
        self.current_output = ''
        self.current_output_filename = ''
        self.current_output_dir = ''
        self.current_output_singlefile = ''
        
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
###   _____                           _   ______                _   _                 
###  / ____|                         | | |  ____|              | | (_)                
### | |  __  ___ _ __   ___ _ __ __ _| | | |__ _   _ _ __   ___| |_ _  ___  _ __  ___ 
### | | |_ |/ _ \ '_ \ / _ \ '__/ _` | | |  __| | | | '_ \ / __| __| |/ _ \| '_ \/ __|
### | |__| |  __/ | | |  __/ | | (_| | | | |  | |_| | | | | (__| |_| | (_) | | | \__ \
###  \_____|\___|_| |_|\___|_|  \__,_|_| |_|   \__,_|_| |_|\___|\__|_|\___/|_| |_|___/
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### Function: check_output:
### Details: handles general output assignment for VESPA functions
def check_output(file_to_verify, default_output):
    if file_to_verify.input_in_dir:
        if not file_to_verify.current_output_dir:
            file_to_verify.current_output_dir = default_output + '_' + file_to_verify.current_input_dir
        if not os.path.exists(file_to_verify.current_output_dir):
            os.makedirs(file_to_verify.current_output_dir)
        if not file_to_verify.current_output_filename:
            file_to_verify.current_output_filename = file_to_verify.current_input_filename
        else:
            file_to_verify.current_output_filename = default_output + '_' + file_to_verify.current_output_filename
        file_to_verify.current_output = '{0}/{1}'.format(file_to_verify.current_output_dir, file_to_verify.current_output_filename)
    else:
        if not file_to_verify.current_output:
            file_to_verify.current_output = '{0}_{1}'.format(default_output,file_to_verify.current_input)
            file_to_verify.current_output_dir = '{0}_{1}'.format(default_output,remove_extension(file_to_verify.current_input))
    return file_to_verify.current_output_dir, file_to_verify.current_output_filename, file_to_verify.current_output


### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### Function: check_output_dir:
### Details: handles directory output assignment for VESPA functions
def check_output_dir(default_output):
    global bme_output_directory
    if not bme_output_directory:
        bme_output_directory = default_output
    if not os.path.exists(bme_output_directory):
        os.makedirs(bme_output_directory)

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### Function: return_filename:
### Details: Returns the filename without any file extention of the user specified filepath
def return_filename_wo_ext(file_path):
    return_file = ''
    if '/' in file_path:
        return_file = file_path.split('/')[-1]
        if '.' in return_file:
            return_file = return_file.split('.',1)[0]
    else:
        return_file = file_path
        if '.' in return_file:
            return_file = return_file.split('.',1)[0]
    return return_file
    
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### Function: return_filename:
### Details: Returns the filename of the user specified filepath
def return_filename(file_path):
    return_file = ''
    if '/' in file_path:
        return_file = file_path.split('/')[-1]
    else:
        return_file = file_path
    return return_file

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### Function: remove_extension:
### Details: Removes file extension from the user specified filepath
def remove_extension(file_path):
    return_file = ''
    if '.' in file_path:
        return_file = file_path.split('.',1)[0]
    else:
        return_file = file_path
    return return_file

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### Function: return_extension:
### Details: Returns the file extension from the user specified filepath
def return_extension(file_path):
    return_file = ''
    if '.' in file_path:
        return_file = file_path.rsplit('.')[-1]
    else:
        return_file = ''
    return return_file

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### Function: return_directory:
### Details: Returns the current directory from the user specified filepath
def return_directory(file_path):
    if '/' in file_path: 
        return file_path.rsplit('/', 1)[0]
    else:
        return False

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### Function: verify_sequence_file:
### Details: Verfies that the file is a sequence file
def verify_sequence_file(file_path):
    verify_sequence = sequence_reader(file_path)
    verify_sequence.sequence_format()
    return verify_sequence.type

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### Function: verify_alignment:
### Details: Verfies that sequence file is an alignment
def verify_alignment(file_path):
    len_list = []
    for check_sequence in sequence_reader(file_path).read():
        len_list.append(len(check_sequence))
    if len(set(len_list)) == 1:
        return True
    else:
        return False

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### Function: ensembl_infer:
### Details: Returns species for ENSEMBL sequences
def ensembl_infer(query_header):
    return_species = ''
    abnormal_ensembl_table = {'FB':'Fruitfly','ENSG':'Human', 'ENSCSAVG':'C_savignyi',}
    ensembl_table = {'ENSAMEG':'Panda', 'ENSAPLG':'Duck', 'ENSACAG':'Anole_lizard',
                     'ENSAMXG':'Cave_fish', 'ENSBTAG':'Cow', 'ENSCELG':'C_elegans',
                     'ENSCJAG':'Marmoset', 'ENSCAFG':'Dog', 'ENSCPOG':'Guinea_Pig',
                     'ENSCHOG':'Sloth', 'ENSCING':'C_intestinalis','ENSXMAG':'Platyfish', 
                     'ENSDARG':'Zebrafish', 'ENSDNOG':'Armadillo', 'ENSDORG':'Kangaroo_rat',
                     'ENSETEG':'Lesser_hedgehog_tenrec', 'ENSECAG':'Horse', 'ENSEEUG':'Hedgehog',
                     'ENSFCAG':'Cat', 'ENSFALG':'Flycatcher', 'ENSGMOG':'Cod',
                     'ENSGALG':'Chicken', 'ENSGACG':'Stickleback', 'ENSGGOG':'Gorilla',
                     'ENSSTOG':'Squirrel', 'ENSLACG':'Coelacanth', 'ENSLOCG':'Spotted_gar',
                     'ENSLAFG':'Elephant', 'ENSMMUG':'Macaque', 'ENSMEUG':'Wallaby',
                     'ENSMGAG':'Turkey', 'ENSMICG':'Mouse_Lemur', 'ENSMODG':'Opossum',
                     'ENSMUSG':'Mouse', 'ENSMPUG':'Ferret', 'ENSMLUG':'Microbat',
                     'ENSNLEG':'Gibbon', 'ENSOPRG':'Pika', 'ENSONIG':'Tilapia',
                     'ENSOANG':'Platypus', 'ENSOCUG':'Rabbit', 'ENSORLG':'Medaka',
                     'ENSOGAG':'Bushbaby', 'ENSOARG':'Sheep', 'ENSPTRG':'Chimpanzee',
                     'ENSPSIG':'Chinese_softshell_turtle', 'ENSPMAG':'Lamprey', 'ENSPPYG':'Orangutan',
                     'ENSPCAG':'Hyrax', 'ENSPVAG':'Megabat', 'ENSRNOG':'Rat',
                     'ENSSCEG':'S_cerevisiae', 'ENSSHAG':'Tasmanian_devil', 'ENSSARG':'Shrew',
                     'ENSSSCG':'Pig', 'ENSTGUG':'Zebra Finch', 'ENSTRUG':'Fugu',
                     'ENSTSYG':'Tarsier', 'ENSTNIG':'Tetraodon', 'ENSTBEG':'Tree Shrew',
                     'ENSTTRG':'Dolphin', 'ENSVPAG':'Alpaca', 'ENSXETG':'Xenopus'}
    
    if ensembl_table.has_key(query_header[1:8]):
        return_species = ensembl_table[query_header[1:8]]
    else:
        for check_species in abnormal_ensembl_table.keys():
            if check_species in query_header:
                return_species = abnormal_ensembl_table[check_species]
    return return_species
    
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### Function: create_unique_file:
### Details: Creates a file, if file already exits a number is added until the filename is unique
def create_unique_file(requested_filename):
    import os
    if '.' in requested_filename:
        filename = requested_filename.strip().rsplit('.',1)[0]
        filename_extension = requested_filename.strip().rsplit('.',1)[1]
        if os.path.isfile('{0}.{1}'.format(filename, filename_extension)):
            counter = 1
            while os.path.isfile('{0}_{1}.{2}'.format(filename, counter, filename_extension)):
                counter += 1
            return open('{0}_{1}.{2}'.format(filename, counter, filename_extension), 'w')
        else:
            return open('{0}.{1}'.format(filename, filename_extension), 'w')
    else:
        filename = requested_filename.strip()
        if os.path.isfile('{0}'.format(filename)):
            counter = 1
            while os.path.isfile('{0}_{1}'.format(filename, counter)):
                counter += 1
            return open('{0}_{1}'.format(filename, counter), 'w')
        else:
            return open('{0}'.format(filename), 'w')

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### Function: check_if_input_directory:
### Details: Checks if specified filepath exists and if so, check if the specified filepath is a file or directory.
def check_if_input_directory (input_varible):
    if os.path.isdir(input_varible):
        return True
    elif os.path.isfile(input_varible):
        return False
    else:
        print 'Input file or directory ({0}) does not exist - exiting program'.format(input_varible)
        sys.exit()
        
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
###  _____ _   _       _____  _                    
### | ____| | | |     |  __ \| |                   
### | |__ | |_| |__   | |__) | |__   __ _ ___  ___ 
### |___ \| __| '_ \  |  ___/| '_ \ / _` / __|/ _ \
###  ___) | |_| | | | | |    | | | | (_| \__ \  __/
### |____/ \__|_| |_| |_|    |_| |_|\__,_|___/\___|
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### Function: vespa_codeml_reader:
### Details: VESPA codeML output reader
def vespa_codeml_reader (input_files):
    print 'VESPA: CodeML Reader'
    import os, sys, glob
    from collections import defaultdict
    
    def codeml_raw_reader(input_files):
        import subprocess, os
        codeml_reader_output = create_unique_file('codeml_reader.log')
        report_list, return_list = [], []
        
        for raw_codeml_input in input_files:
            omega_position = 0
            split_codeml = raw_codeml_input.current_input.strip().split('/')
            for pos, string in enumerate(split_codeml):
                if 'Omega' in string:
                    omega_position = len(split_codeml) - (pos - 2)
            if omega_position != 0:
                report_list.append(raw_codeml_input.current_input.strip().rsplit('/',omega_position)[0])
        report_list = list(set(report_list))
        
        for orginal_path in report_list:
            report_split = orginal_path.split('/',1)
            codeml_wrapper_check, codeml_wrapper_bin = False, False  
            if not bme_output_directory:
                report_output_dir = 'Report_{0}'.format(report_split[0])
            check_output_dir(report_output_dir)
            report_path = '{0}/{1}'.format(bme_output_directory,report_split[1])
            if not os.path.exists(report_path):
                os.makedirs(report_path)
            if not codeml_wrapper_check:
                try:
                    codeml_wrapper_call = subprocess.Popen(['CreateCodemlReports.pl', orginal_path, report_path, '-overwrite=Yes'], stdout=subprocess.PIPE,stderr=subprocess.PIPE, stdin=subprocess.PIPE)
                    codeml_wrapper_check = True
                except:
                    codeml_wrapper_bin = True
            if codeml_wrapper_bin:
                try:
                    codeml_wrapper_call = subprocess.Popen(['perl', 'CreateCodemlReports.pl', orginal_path, report_path, '-overwrite=Yes'], stdout=subprocess.PIPE,stderr=subprocess.PIPE, stdin=subprocess.PIPE)
                    codeml_wrapper_check = True
                except:
                    pass
            if codeml_wrapper_check:
                wrapper_out, wrapper_error = codeml_wrapper_call.communicate()
                if not wrapper_error:
                    codeml_reader_output.write('Currently Parsing: {0}\n'.format(orginal_path))
                    codeml_reader_output.write(wrapper_out)
                    for path, sub_dirs, file_list in os.walk(report_path):
                        for files in file_list:
                            if not files.startswith('.'):
                                return_list.append(os.path.join(path, files))
                else:
                    print wrapper_error
                    codeml_reader_output.write('Error running CreateCodemlReports.pl, please confirm that the script and all modules are installed.\n')
            else:
                codeml_reader_output.write('Error running CreateCodemlReports.pl, please confirm that the script and all modules are installed.\n')

        codeml_reader_output.close()
        return return_list

    def append_character_data(codeml_matched_site_apd, extand_lineage_seq_apd, ps_char_seq_apd, seq_type_apd):
        if seq_type_apd == 'DNA':
            ps_char_seq_apd[int(codeml_matched_site_apd) * 3] = extand_lineage_seq_apd[int(codeml_matched_site_apd) * 3]
            ps_char_seq_apd[(int(codeml_matched_site_apd) * 3) + 1] = extand_lineage_seq_apd[(int(codeml_matched_site_apd) * 3) + 1]
            ps_char_seq_apd[(int(codeml_matched_site_apd) * 3) + 2] = extand_lineage_seq_apd[(int(codeml_matched_site_apd) * 3) + 2]
        elif seq_type_apd == 'protein':
            ps_char_seq_apd[int(codeml_matched_site_apd)] = extand_lineage_seq_apd[int(codeml_matched_site)]
        return ps_char_seq_apd
    
    def append_site_data(codeml_matched_site_apd, ps_site_seq_apd, seq_type_apd):
        if seq_type_apd == 'DNA':
            ps_site_seq_apd[int(codeml_matched_site_apd) * 3] = 'N'
            ps_site_seq_apd[(int(codeml_matched_site_apd) * 3) + 1] = 'N'
            ps_site_seq_apd[(int(codeml_matched_site_apd) * 3) + 2] = 'N'
        elif seq_type_apd == 'protein':
            ps_site_seq_apd[int(codeml_matched_site_apd)] = 'X'
        return ps_site_seq_apd
    
    global bme_branch_label_table, bme_alignment_path, bme_output_directory
    branch_models = False
    
    sequence_type = ''
    sequence_input_files = []
    
    if bme_branch_label_table:
        if os.path.isfile(bme_branch_label_table):
            branch_models = True
            ancestral_lineage_matcher = defaultdict(list)
            with open(bme_branch_label_table, 'rU') as label_data:
                for label_lines in label_data:
                    if ':' in label_lines:
                        ancestral_lineage_matcher[label_lines.strip().split(': ')[0]].extend([temp_labels.strip() for temp_labels in label_lines.strip().split(': ')[-1].split(',')])
                    else:
                        ancestral_lineage_matcher[label_lines.strip()].append(label_lines.strip())
        else:
            print 'Could not locate specified label table file, please check'
            sys.exit()
    else:
        print 'Branch labels not specified, please include the option: -label_table=USR_TBL. Only site models alignments will be created'
    
    if bme_alignment_path:
        if check_if_input_directory(bme_alignment_path):
            for path, sub_dirs, file_list in os.walk(bme_alignment_path):
                for files in file_list:
                    if not files.startswith('.'):
                        sequence_input_files.append(os.path.join(path, files))
        else:
            sequence_input_files.append(bme_alignment_path)
    else:
        print 'Alignments path not specified, please use -alignment_path='
        sys.exit()
    
    report_files = codeml_raw_reader(input_files)
    for sequence_input in sequence_input_files:
        sequence_type, current_sequence_file = '', return_filename_wo_ext(sequence_input)
        current_sequences = {}
        for working_sequence in sequence_reader(sequence_input).read():
            current_sequences[working_sequence.header.split('|')[0][1:]] = working_sequence
            if not sequence_type:
                working_sequence.seq_type()
                sequence_type = working_sequence.type
        for codeml_input in report_files:
            if current_sequence_file in codeml_input.split(os.sep):
                if 'PosSites' in codeml_input:
                    if branch_models:
                        if 'modelA.txt' in codeml_input:
                            current_lineage = []
                            for identifier_labels in ancestral_lineage_matcher.keys():
                                if identifier_labels in return_filename(codeml_input):
                                    ps_character_seqeunce_dict = {}
                                    get_length = []
                                    for extant_species in ancestral_lineage_matcher[identifier_labels]:
                                        if current_sequences.has_key(extant_species):
                                            get_length.append(len(current_sequences[extant_species]))
                                            ps_character_seqeunce_dict[extant_species] = list('-' * len(current_sequences[extant_species]))
                                    if len(set(get_length)) == 1:
                                        ps_site_seqeunce = list('-' * get_length[0])
                                        with open(codeml_input, 'rU') as codeml_data:
                                            for codeml_lines in codeml_data:
                                                if 'P(w>1) > 0.5' not in codeml_lines:
                                                    codeml_site_data = codeml_lines.strip().split()
                                                    codeml_matched_site = int(codeml_site_data[0]) - 1
                                                    ps_site_seqeunce = append_site_data(codeml_matched_site, ps_site_seqeunce, sequence_type)
                                                    for extant_species in ancestral_lineage_matcher[identifier_labels]:
                                                        if current_sequences.has_key(extant_species):
                                                            ps_character_seqeunce_dict[extant_species] = append_character_data(codeml_matched_site, current_sequences[extant_species].sequence, ps_character_seqeunce_dict[extant_species], sequence_type)
                                        codeml_alignment = open('{0}.fasta'.format(remove_extension(codeml_input)), 'w')
                                        codeml_alignment.write(str(sequence_data('>PS_Sites|{0}\n'.format(identifier_labels),''.join(ps_site_seqeunce))) + '\n')
                                        for keys in current_sequences.keys():
                                            if keys in ancestral_lineage_matcher[identifier_labels]:
                                                codeml_alignment.write(str(sequence_data('>PS_Characters|' + current_sequences[keys].header[1:],''.join(ps_character_seqeunce_dict[keys]))) + '\n')
                                                codeml_alignment.write(str(current_sequences[keys]) + '\n')
                                        for keys in current_sequences.keys():
                                            if keys not in ancestral_lineage_matcher[identifier_labels]:
                                                codeml_alignment.write(str(current_sequences[keys]) + '\n')
                                        codeml_alignment.close()
                                    else:
                                        print 'Error: Not specified'
                        elif '.txt' in codeml_input:
                            get_length = [len(current_sequences[keys]) for keys in current_sequences.keys()]
                            if len(set(get_length)) == 1:
                                ps_site_seqeunce = list('-' * get_length[0])
                                with open(codeml_input, 'rU') as codeml_data:
                                    for codeml_lines in codeml_data:
                                        if 'P(w>1) > 0.5' not in codeml_lines:
                                            codeml_site_data = codeml_lines.strip().split()
                                            codeml_matched_site = int(codeml_site_data[0]) - 1
                                            ps_site_seqeunce = append_site_data(codeml_matched_site, ps_site_seqeunce, sequence_type)
                                codeml_alignment = open('{0}.fasta'.format(remove_extension(codeml_input)), 'w')
                                codeml_alignment.write(str(sequence_data('>PS_Sites\n',''.join(ps_site_seqeunce))) + '\n')
                                for keys in current_sequences.keys():
                                    codeml_alignment.write(str(current_sequences[keys]) + '\n')
                                codeml_alignment.close()

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
###  _  _   _   _       _____  _                    
### | || | | | | |     |  __ \| |                   
### | || |_| |_| |__   | |__) | |__   __ _ ___  ___ 
### |__   _| __| '_ \  |  ___/| '_ \ / _` / __|/ _ \
###    | | | |_| | | | | |    | | | | (_| \__ \  __/
###    |_|  \__|_| |_| |_|    |_| |_|\__,_|___/\___|
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### Function: vespa_branch_table:
### Details: VESPA function that creates a branch table from species tree
def vespa_branch_table(input_files):
    print 'VESPA: Branch Table Creator'
    import os, sys, dendropy, re, copy
    from collections import defaultdict
    
    def select_subtree(sent_tree):
        pre_edited_tree = dendropy.Tree.get_from_string(sent_tree,"newick")
        finished_node_selection = False
        return_nodes = defaultdict(list)
        while not finished_node_selection:
            user_selected_node = raw_input('Please select a node for selection (numerical values): ') 
            if user_selected_node in [internal_nodes.label for internal_nodes in pre_edited_tree.internal_nodes()]:
                user_confirm = raw_input('Node {0} found. Please confirm (y / n): '.format(user_selected_node))
                if user_confirm.lower().startswith('y'):
                    edited_tree = copy.deepcopy(pre_edited_tree)
                    edited_node = edited_tree.find_node_with_label(user_selected_node)
                    edited_leafs =  [leaf_nodes.taxon for leaf_nodes in edited_node.leaf_iter()]
                    edited_leaf_taxa = [leaf_nodes.label for leaf_nodes in edited_leafs]
                    edited_tree.retain_taxa(edited_leafs)
                    confirm_ID = False
                    while not confirm_ID:
                        node_ID = raw_input('Please specify a label for this branch (used for codeml analysis): ')
                        if node_ID: 
                            id_confirm = raw_input('Label: {0}? Please confirm (y / n): '.format(node_ID))
                            if id_confirm.lower().startswith('y'):
                                confirm_ID = True
                                return_nodes[node_ID] = copy.deepcopy(edited_leaf_taxa)
            another_node_check = raw_input('Select an additional branch? Please confirm (y / n): ')
            if another_node_check.lower().startswith('n'):
                finished_node_selection = True
        return return_nodes

    def select_leaf(sent_taxa):
        finished_leaf_selection = False
        return_list = []
        while not finished_leaf_selection:
            user_selected_leafs = raw_input('Please select a leaf (taxa) for selection (if mutiple, seperate with comma): ')
            leaf_list = [current_leaf.strip() for current_leaf in user_selected_leafs.split(',')]         
            if len(leaf_list) == len(list(set(leaf_list) & set(sent_taxa))):
                user_confirm = raw_input('Leaf(s) {0} found. Please confirm (y / n): '.format(', '.join(leaf_list)))
                if user_confirm.lower().startswith('y'):
                    finished_leaf_selection = True
                    return_list = copy.deepcopy(leaf_list)
            else:
                unknown_leafs = list(set(leaf_list) - set(set(leaf_list) & set(sent_taxa)))
                print 'Leaf(s) {0} not found. Please confirm taxa spelling and formatting.'.format(', '.join(unknown_leafs))
        return return_list

    for species_tree in input_files:
        def tree_labeler (sent_tree):
            labeler = list(sent_tree)
            for pos, matched_nodes in enumerate([match.end() for match in re.finditer('\)', sent_tree)][::-1]):
                labeler.insert(matched_nodes,str(pos))
            return ''.join(labeler) + ';'
        
        branch_file = create_unique_file('branch_table.txt')
        
        saved_leafs = []
        saved_braches = []
        data_tree = dendropy.Tree.get_from_path(species_tree.current_input,"newick")
        original_tree = data_tree.as_string(schema="newick")
        original_taxa = [str(taxa_ids).replace("'",'') for taxa_ids in data_tree.taxon_namespace]
        
        screen_tree = dendropy.Tree.get_from_string(tree_labeler(original_tree),"newick")
        screen_tree.print_plot(show_internal_node_labels=True)
        screen_string = screen_tree.as_string(schema="newick")
        
        command_dict = {'0':'Finished', '1':'Species Selection', '2':'Ancestral Lineage Selection'}
        finished_selection = False
        while not finished_selection:
            print 'Possible actions\n____________________\n1. Species Selection\n2. Ancestral Lineage Selection\n0. Finished\n'
            user_request = raw_input('Please select an action ( 0 / 1 / 2 ): ')
            if user_request in ['0', '1', '2']:
                user_confirm = raw_input(command_dict[user_request] + '. Please confirm (y / n): ')
                if user_confirm.lower().startswith('y'):
                    if user_request == '1':
                        saved_leafs.extend(select_leaf(original_taxa))
                    elif user_request == '2':
                        saved_braches.append(select_subtree(screen_string))
                    elif user_request == '0':
                        finished_selection = True
        
        for current_leaf in saved_leafs:
            branch_file.write(current_leaf + '\n')
        for current_dict in saved_braches:
            for label, species in current_dict.items():
                branch_file.write('{0}: {1}\n'.format(label,', '.join(species)))
        branch_file.close()

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### Function: vespa_map_protein_gaps:
### Details: Creates nucleotide alignments using protein alignments and a nucleotide sequence database
def vespa_map_protein_gaps(input_files):
    print 'VESPA: Mapping Nucleotide Alignments'
    from collections import defaultdict
    global bme_sequence_database_location
    alignment_counter, alignment_map, protein_map = (defaultdict(list), defaultdict(list), defaultdict(list))

    def cleave_stop_codon(sequence):
        stop_codon_list = ['TAA','TAG','TGA']
        if sequence[-3:] in stop_codon_list:
            return sequence[:-3]
        else:
            return sequence
    
    if not bme_sequence_database_location:
        print 'Sequence database not found'
        sys.exit()
        
    for sequence_input in input_files:
        (mapped_output_directory, mapped_output_filename, mapped_output_path) = check_output(sequence_input, 'Map_Gaps')
        for working_sequence in sequence_reader(sequence_input.current_input).read():
            mapped_output_path = '{0}/{1}'.format(mapped_output_directory, return_filename(sequence_input.current_input))
            amino_acid_positions = []
            if working_sequence.sequence[-1] == '*':
                working_sequence.sequence = working_sequence.sequence[:-1]
            for position, amino_acid in enumerate(working_sequence.sequence):
                if amino_acid != '-':
                    amino_acid_positions.extend([position * 3, (position * 3) + 1, (position * 3) + 2])
            protein_map[working_sequence.header] = amino_acid_positions
            alignment_map[working_sequence.header] = [mapped_output_path, len(working_sequence) * 3, len(working_sequence.sequence.replace('-','')) * 3]
            alignment_counter[mapped_output_path].append(working_sequence.header)
       
    for mapped_alignments in alignment_counter.keys():
        try:
            with open(mapped_alignments): os.remove(mapped_alignments)
        except IOError:
            pass
        
    verify_database = False
    for genome_sequence in sequence_reader(bme_sequence_database_location).read():
        if not verify_database:
            verify_database = True
            genome_sequence.seq_type()
            if genome_sequence.type == 'protein':
                print 'Protein database detected, please use a nucleotide database'
                sys.exit()
        if alignment_map.has_key(genome_sequence.header):
            alignment_output_filename = alignment_map[genome_sequence.header][0]
            protein_adjested_sequence = cleave_stop_codon(genome_sequence.sequence)
            if len(protein_adjested_sequence) == alignment_map[genome_sequence.header][2]:
                nucleotide_alignment = list(alignment_map[genome_sequence.header][1] * '-')
                for nucleotide, position in zip(protein_adjested_sequence, protein_map[genome_sequence.header]):
                    nucleotide_alignment[position] = nucleotide
                sequence_output = open(alignment_output_filename, 'a')
                sequence_output.write(str(sequence_data(genome_sequence.header, ''.join(nucleotide_alignment))) + '\n')
                sequence_output.close()
            else:
                print 'Sequence length differences identified, please verify that the protein sequences have not been altered since translation'
                sys.exit()
           
    for alignment_files in alignment_counter.keys():
        written_header_list = []
        with open(alignment_files) as alignment_data:
            for alignment_lines in alignment_data:
                if alignment_lines.startswith('>'):
                    written_header_list.append(alignment_lines)
        if len(written_header_list) != len(alignment_counter[alignment_files]):
            convert_list = []
            for written_header in written_header_list:
                convert_list.append(written_header.replace(' ', '_').replace(':', '_').replace('(', '_').replace(')', '_').replace(',', '_'))
            print 'Error in:', alignment_files + ',' + ','.join([print_id.strip() for print_id in list(set(alignment_counter[alignment_files]) - set(convert_list))])

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### Function: vespa_mrbayes_reader:
### Details: Reads mrbayes output and converts to newick
def vespa_mrbayes_reader (input_files):
    print 'VESPA: Reading MrBayes Output'
    import re
    for alignment_input in input_files:
        (mrbayes_output_dir, mrbayes_output_filename, mrbayes_output) = check_output(alignment_input, 'MrBayes_Reader')
        newick_output = open(remove_extension(mrbayes_output) + '.tre', 'w')
        with open(alignment_input.current_input) as nexus_data:
            check_nexus = nexus_data.readline()
            tree_list = []
            if check_nexus.strip() == '#NEXUS':
                parsing_trees = False
                for nexus_lines in nexus_data:
                    if 'end;' in nexus_lines:
                        parsing_trees = False
                    elif parsing_trees:
                        if nexus_lines.strip().startswith('tree'):
                            current_tree = nexus_lines.strip().split('=')[1].strip()
                            tree_list.append(re.sub('\)\d.\d+', ')',re.sub(':\d.\d+', '', current_tree)))
                    elif 'begin trees;' in nexus_lines:
                        parsing_trees = True
            if len(set(tree_list)):
                newick_output.write('{0}\n'.format(list(set(tree_list))[0]))
            else:
                print 'Error parsing NEXUS file. Please comfirm file is from MrBayes.'
        newick_output.close()

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### Function: vespa_infer_genetree:
### Details: Function that generates gene trees using a specififed spcies tree
def vespa_infer_genetree (input_files):
    print 'VESPA: Inferring GeneTree'
    import dendropy, re, os, shutil
    global bme_species_tree, bme_in_paralogs
    if bme_species_tree:
        if os.path.isfile(bme_species_tree):
            try:
                species_tree = dendropy.Tree.get(path=bme_species_tree, schema="newick", preserve_underscores=True)
            except:
                print 'Error reading: {0}. Please verify that species tree is in newick format'.format(bme_species_tree)
                sys.exit()
            taxa_list = [str(taxa_ids).replace("'",'') for taxa_ids in species_tree.taxon_namespace]
            error_log = create_unique_file('infer_genetree.log')
            for alignment_input in input_files:
                sequence_counter, in_paralog_counter = 0, 0
                if verify_alignment(alignment_input.current_input):
                    convert_tree = dendropy.Tree(species_tree)
                    taxa_check = dict([(taxa,'') for taxa in taxa_list])
                    for working_sequence in sequence_reader(alignment_input.current_input).read():
                        sequence_counter += 1
                        if working_sequence.header.split('|')[0][1:] in taxa_list:
                            if taxa_check[working_sequence.header.split('|')[0][1:]]:
                                if bme_in_paralogs:
                                    in_paralog_counter += 1
                                    if '(' in taxa_check[working_sequence.header.split('|')[0][1:]]:
                                        current_leaf = taxa_check[working_sequence.header.split('|')[0][1:]][1:-1]
                                    else:
                                        current_leaf = taxa_check[working_sequence.header.split('|')[0][1:]]
                                    taxa_check[working_sequence.header.split('|')[0][1:]] = '({0},{1})'.format(current_leaf, working_sequence.header.strip()[1:])
                                else:
                                    print 'Duplication detected: {0}. Please check files'.format(working_sequence.header.strip()[1:])
                                    error_log.write('Duplication detected: {0}. Please check files\n'.format(working_sequence.header.strip()[1:]))
                            else:  
                                taxa_check[working_sequence.header.split('|')[0][1:]] = working_sequence.header.strip()[1:]
                        else:
                            print 'Non-tree species detected: {0}. Please check files'.format(working_sequence.header.strip()[1:])
                            error_log.write('Non-tree species detected: {0}. Please check files\n'.format(working_sequence.header.strip()[1:]))
                    if '' in taxa_check.values():
                        nodes_to_prune = [convert_tree.find_node_with_taxon_label(taxa).taxon for taxa in taxa_check.keys() if taxa_check[taxa] == '']
                        convert_tree.prune_taxa(nodes_to_prune)
                    pruned_tree = convert_tree.as_string(schema="newick").replace("'",'')
                    pruned_taxa = dendropy.Tree.get(data=pruned_tree, schema="newick").taxon_namespace.labels()
                    if len(pruned_taxa) == (sequence_counter - in_paralog_counter):
                        for raw_taxa, gene in taxa_check.items():
                            if gene:
                                replace_finder = re.search(raw_taxa +"(\)|\,|;)", pruned_tree)
                                pruned_tree = pruned_tree.replace(replace_finder.group(), replace_finder.group().replace(raw_taxa, gene))
                        (inferred_output_dir, inferred_output_filename, inferred_output) = check_output(alignment_input, 'Inferred_Genetree')
                        current_filename = return_filename_wo_ext(alignment_input.current_input)
                        sub_output_path = '{0}/{1}'.format(inferred_output_dir,current_filename)
                        if not os.path.exists(sub_output_path):
                            os.makedirs(sub_output_path)
                        shutil.copy(alignment_input.current_input,sub_output_path)
                        tree_file = open('{0}/{1}.tre'.format(sub_output_path,current_filename), 'w')
                        tree_file.write(pruned_tree)    
                        tree_file.close()
                    else:
                        error_log.write('Error pruning: {0}.tre\n'.format(remove_extension(alignment_input.current_input)))
                else:
                    print 'Unaligned seqeunces detected in: {0}'.format(return_filename(alignment_input.current_input))
                    error_log.write('Unaligned seqeunces detected in: {0}\n'.format(return_filename(alignment_input.current_input)))
            error_log.close()
        else:
            print 'Could not locate specified tree file, please check'
            sys.exit()
    else:
        print 'Species tree not specified, please include the option: -species_tree=USR_TRE'
        sys.exit()

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### Function: vespa_link_input:
### Details: Combines alignment and tree data for Codeml
def vespa_link_input(input_files):
    print 'VESPA: Linking Input'
    import dendropy, os, sys, shutil
    global bme_alignment_path
    if bme_alignment_path:
        error_log = create_unique_file('link_input.log')
        for tree_files in input_files:
            taxa_list = []
            try:
                gene_tree = dendropy.Tree.get_from_path(tree_files.current_input, schema="newick",preserve_underscores=True)
                taxa_list = [str(taxa_ids).replace("'",'') for taxa_ids in gene_tree.taxon_namespace]
            except IOError:
                print 'Non-Tree file detected: {0}'.format(return_filename(tree_files.current_input))
                error_log.write('Non-Tree file detected: {0}\n'.format(return_filename(tree_files.current_input)))
            if taxa_list:
                alignment_files = []
                if check_if_input_directory(bme_alignment_path):
                    for path, sub_dirs, file_list in os.walk(bme_alignment_path):
                        for files in file_list:
                            if not files.startswith('.'):
                                alignment_files.append(os.path.join(path, files))
                else:
                    alignment_files.append(bme_alignment_path)
                
                for alignment_input in alignment_files:
                    alignment_taxa = []
                    if verify_alignment(alignment_input):
                        if return_filename_wo_ext(tree_files.current_input) == return_filename_wo_ext(alignment_input):
                            for working_sequence in sequence_reader(alignment_input).read():
                                alignment_taxa.append(working_sequence.header.strip()[1:])
                            if set(taxa_list) == set(alignment_taxa):
                                (linked_dir, linked_file, linked_output) = check_output(tree_files, 'Linked')
                                if not os.path.exists(linked_dir):
                                    os.makedirs(linked_dir)
                                if tree_files.input_in_dir:
                                    sub_output_dir = '{0}/{1}'.format(linked_dir,remove_extension(return_filename(tree_files.current_input)))
                                    if not os.path.exists(sub_output_dir):
                                        os.makedirs(sub_output_dir)
                                    shutil.copy(tree_files.current_input,sub_output_dir)
                                    shutil.copy(alignment_input,sub_output_dir)
                                else:
                                    shutil.copy(tree_files.current_input,linked_dir)
                                    shutil.copy(alignment_input,linked_dir)
                    else:
                        print 'Unaligned seqeunces detected in: {0}'.format(return_filename(alignment_input))
                        error_log.write('Unaligned seqeunces detected in: {0}\n'.format(return_filename(alignment_input)))
        error_log.close()   
    else:
        print 'Alignment directory not specified, please include the option: -alignment_path=USR_DIR'
        sys.exit()
                
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### Function: vespa_codeml_setup:
### Details: Takes combined alignment and tree data and prepares the data for Codeml
def vespa_codeml_setup(input_files):
    print 'VESPA: CodeML Setup'
    import dendropy, os, sys, subprocess, copy
    from collections import defaultdict
    global bme_in_paralogs, bme_branch_label_table
    error_log = create_unique_file('codeml_setup.log')
    taskfarm_log = create_unique_file('codeml_taskfarm.txt')
    branch_models = False
    alignment_found_check, tree_found_check = False, False
    
    if bme_branch_label_table:
        if os.path.isfile(bme_branch_label_table):
            branch_models = True
            label_dict = defaultdict(list)
            with open(bme_branch_label_table, 'rU') as label_data:
                for label_lines in label_data:
                    if ':' in label_lines:
                        label_dict['node'].append([label_lines.strip().split(': ')[0], [temp_labels.strip() for temp_labels in label_lines.strip().split(': ')[-1].split(',')]])
                    else:
                        label_dict['leaf'].append(label_lines.strip())
        else:
            print 'Could not locate specified label table file, please check'
            sys.exit()
    else:
        print 'Branch labels not specified, please include the option: -label_table=USR_TBL. Setup will be site models only'
        error_log.write('Branch labels not specified, please include the option: -label_table=USR_TBL. Setup will be site models only\n')

    for alignment_input in input_files:
        if verify_sequence_file(alignment_input.current_input) != '':
            species_convert = defaultdict(list)                       
            if verify_alignment(alignment_input.current_input):
                alignment_found_check = True
                taxa_list = []
                for tree_files in input_files:
                    if remove_extension(tree_files.current_input) == remove_extension(alignment_input.current_input):
                        if tree_files.current_input != alignment_input.current_input:
                            try:
                                tree_input = tree_files.current_input
                                original_tree = dendropy.Tree.get_from_path(tree_input, schema="newick",preserve_underscores=True)
                                taxa_list = [str(taxa_ids).replace("'",'') for taxa_ids in original_tree.taxon_namespace]
                            except IOError:
                                pass
                if taxa_list:
                    tree_found_check = True
                    if bme_in_paralogs:
                        for taxa in taxa_list:
                            species_convert[taxa.split('|')[0]].append(taxa)
                    else:
                        if len(set([taxa.split('|')[0] for taxa in taxa_list])) == len(taxa_list):
                            for taxa in taxa_list:
                                species_convert[taxa.split('|')[0]].append(taxa)
                        else:
                            print 'Duplication detected: {0}. Please check file'.format(alignment_input.current_input)
                            error_log.write('Duplication detected: {0}. Please check file\n'.format(alignment_input.current_input))
                            continue
                    codeml_trees = [tree_input]
                    if branch_models:
                        if label_dict.has_key('leaf'):
                            for species in label_dict['leaf']:
                                if species_convert.has_key(species):
                                    tree_file_name = '{0}_{1}.tre'.format(remove_extension(tree_input), species)
                                    if len(species_convert[species]) > 1:
                                        working_tree = copy.deepcopy(original_tree)
                                        mrca_node = working_tree.mrca(taxon_labels=species_convert[species])
                                        mrca_leafs =  [leaf_nodes.taxon for leaf_nodes in mrca_node.leaf_iter()]
                                        mrca_leaf_taxa = [leaf_nodes.label for leaf_nodes in mrca_leafs]
                                        mrca_tree = dendropy.Tree(working_tree)
                                        mrca_tree.retain_taxa(mrca_leafs)
                                        mrca_newick = mrca_tree.as_string(schema="newick", suppress_rooting=True).strip().replace(';','')
                                        if set(mrca_leaf_taxa) == set(species_convert[species]):                                      
                                            label_tree = original_tree.as_string(schema="newick", suppress_rooting=True)
                                            label_tree = label_tree.replace("'",'').replace(mrca_newick.replace("'",''), mrca_newick.replace("'",'') + "'#1'")
                                            tree_file = open(tree_file_name, 'w')
                                            tree_file.write(label_tree)
                                            tree_file.close()
                                            codeml_trees.append(tree_file_name)
                                        else:
                                            print 'Unable to label: {0}. Additional genes present'.format(tree_file_name)
                                            error_log.write('Unable to label: {0}. Additional genes present\n'.format(tree_file_name))
                                        
                                        for gene_conversions in species_convert[species]:
                                            paralog_file_name = '{0}_{1}.tre'.format(remove_extension(tree_input), gene_conversions.split('|')[1].strip())
                                            label_tree = original_tree.as_string(schema="newick", suppress_rooting=True)
                                            current_tree = label_tree.replace(gene_conversions, gene_conversions + '#1').replace("'",'')
                                            tree_file = open(paralog_file_name, 'w')
                                            tree_file.write(current_tree)
                                            tree_file.close()
                                            codeml_trees.append(paralog_file_name)
                                    else:
                                        for gene_conversions in species_convert[species]:
                                            label_tree = original_tree.as_string(schema="newick", suppress_rooting=True)
                                            current_tree = label_tree.replace(gene_conversions, gene_conversions + '#1').replace("'",'')
                                            tree_file = open(tree_file_name, 'w')
                                            tree_file.write(current_tree)
                                            tree_file.close()
                                            codeml_trees.append(tree_file_name)
                                else:
                                    print 'Unable to label {0} for {1}. Cannot find {0} in genetree. {0} will not be included in codeML analysis for {1}.'.format(species, return_filename(tree_input))
                                    error_log.write('Unable to label {0} for {1}. Cannot find {0} in genetree. {0} will not be included in codeML analysis for {1}.\n'.format(species, return_filename(tree_input)))
                        if label_dict.has_key('node'):
                            for node_data in label_dict['node']:
                                tree_file_name = remove_extension(tree_input) + '_' + node_data[0] + '.tre'
                                species_list = node_data[1]
                                node_genes, nodes_found = ([], [])
                                paralog_counter = 0
                                for node_species in species_list:
                                    if species_convert.has_key(node_species):
                                        if bme_in_paralogs:
                                            if len(species_convert[node_species]) > 1:
                                                paralog_counter += (len(species_convert[node_species]) - 1)
                                        for append_genes in species_convert[node_species]: 
                                            node_genes.append(append_genes)
                                        nodes_found.append(node_species)
                                    else:
                                        print 'Unable to label {0} for {1}. Cannot find {2} in genetree. {0} will not be included in codeML analysis for {1}.'.format(node_data[0], return_filename(tree_input), node_species)
                                        error_log.write('Unable to label {0} for {1}. Cannot find {2} in genetree. {0} will not be included in codeML analysis for {1}.\n'.format(node_data[0], return_filename(tree_input), node_species))
                                if (len(node_genes) - paralog_counter) == len(species_list):
                                    working_tree = copy.deepcopy(original_tree)
                                    mrca_node = working_tree.mrca(taxon_labels=node_genes)
                                    mrca_leafs =  [leaf_nodes.taxon for leaf_nodes in mrca_node.leaf_iter()]
                                    mrca_leaf_taxa = [leaf_nodes.label for leaf_nodes in mrca_leafs]
                                    mrca_tree = dendropy.Tree(working_tree)
                                    mrca_tree.retain_taxa(mrca_leafs)
                                    mrca_newick = mrca_tree.as_string(schema="newick", suppress_rooting=True).strip().replace(';','')
                                    if set(mrca_leaf_taxa) == set(node_genes):
                                        label_tree = original_tree.as_string(schema="newick", suppress_rooting=True)
                                        label_tree = label_tree.replace("'",'').replace(mrca_newick.replace("'",''), mrca_newick.replace("'",'') + "'#1'")
                                        tree_file = open(tree_file_name, 'w')
                                        tree_file.write(label_tree)
                                        tree_file.close()
                                        codeml_trees.append(tree_file_name)
                                    else:
                                        print 'Unable to label: {0}. Additional genes present'.format(node_data)
                                        error_log.write('Unable to label: {0}. Additional genes present\n'.format(node_data))
                                else:
                                    reduced_list_not_already_labeled = False
                                    for check_nodes in label_dict['node']:
                                        if node_data[0] != check_nodes[0]:
                                            if set(nodes_found) == set(check_nodes[1]):
                                               reduced_list_not_already_labeled = True
                                    if len(nodes_found) == 1:
                                        for check_leaves in label_dict['leaf']:
                                            if str(nodes_found[0]) == str(check_leaves):
                                                reduced_list_not_already_labeled = True
                                    if not reduced_list_not_already_labeled and len(nodes_found) != 0:
                                        working_tree = copy.deepcopy(original_tree)
                                        mrca_node = working_tree.mrca(taxon_labels=node_genes)
                                        mrca_leafs =  [leaf_nodes.taxon for leaf_nodes in mrca_node.leaf_iter()]
                                        mrca_leaf_taxa = [leaf_nodes.label for leaf_nodes in mrca_leafs]
                                        mrca_tree = dendropy.Tree(working_tree)
                                        mrca_tree.retain_taxa(mrca_leafs)
                                        mrca_newick = mrca_tree.as_string(schema="newick", suppress_rooting=True).strip().replace(';','')
                                        if set(node_genes) == set(mrca_leaf_taxa):
                                            label_tree = original_tree.as_string(schema="newick", suppress_rooting=True)
                                            label_tree = label_tree.replace("'",'').replace(mrca_newick.replace("'",''), mrca_newick.replace("'",'') + "'#1'")
                                            tree_file = open(tree_file_name, 'w')
                                            tree_file.write(label_tree)
                                            tree_file.close()
                                            codeml_trees.append(tree_file_name)
                                        else:
                                            print 'Unable to label: {0}. Additional genes present'.format(node_data)
                                            error_log.write('Unable to label: {0}. Additional genes present\n'.format(node_data))
                    (codeml_dir, codeml_file, codeml_output) = check_output(alignment_input, 'Codeml_Setup')
                    main_output_dir = codeml_dir
                    codeml_wrapper_check, codeml_wrapper_bin = False, False
                    if alignment_input.input_in_dir: 
                        setup_list = ['GenerateCodemlWorkspace.pl',alignment_input.current_input] + codeml_trees + [codeml_dir]
                        try:
                            codeml_wrapper_call = subprocess.Popen(setup_list, stdout=subprocess.PIPE,stderr=subprocess.PIPE, stdin=subprocess.PIPE)
                            codeml_wrapper_check = True
                        except:
                            codeml_wrapper_bin = True
                        if codeml_wrapper_bin:
                            setup_list = ['perl'] + setup_list
                            try:
                                codeml_wrapper_call = subprocess.Popen(setup_list, stdout=subprocess.PIPE,stderr=subprocess.PIPE, stdin=subprocess.PIPE)
                                codeml_wrapper_check = True
                            except:
                                pass
                        if codeml_wrapper_check:
                            wrapper_out, wrapper_error = codeml_wrapper_call.communicate()
                            if not wrapper_error:
                                for path, sub_dirs, file_list in os.walk(codeml_dir):
                                    if 'Omega' in path:
                                        taskfarm_log.write('cd {0}; codeml\n'.format(path))
                            else:
                                if 'number of sequences' in wrapper_error:
                                    error_string = wrapper_error.strip().split(': ')[1].split(',')[0]
                                    error_string = error_string[0].upper() + error_string[1:]
                                    print error_string
                                    error_log.write('Error running GenerateCodemlWorkspace.pl.\n')
                                    error_log.write('Error Reported: {0}.\n'.format(error_string))   
                                else:
                                    print 'Error running GenerateCodemlWorkspace.pl. Please check log file for details.'
                                    error_log.write('Error running GenerateCodemlWorkspace.pl.\n')
                                    error_log.write('{0}.\n'.format(wrapper_error))    
                        else:
                            print 'Error running GenerateCodemlWorkspace.pl, please confirm that the script and all modules are installed.'
                            error_log.write('Error running GenerateCodemlWorkspace.pl, please confirm that the script and all modules are installed.\n')
                    else:
                        setup_list = ['GenerateCodemlWorkspace.pl',alignment_input.current_input] + codeml_trees + [codeml_dir]
                        try:
                            codeml_wrapper_call = subprocess.Popen(setup_list, stdout=subprocess.PIPE,stderr=subprocess.PIPE, stdin=subprocess.PIPE)
                            codeml_wrapper_check = True
                        except:
                            codeml_wrapper_bin = True
                        if codeml_wrapper_bin:
                            setup_list = ['perl'] + setup_list
                            try:
                                codeml_wrapper_call = subprocess.Popen(setup_list, stdout=subprocess.PIPE,stderr=subprocess.PIPE, stdin=subprocess.PIPE)
                                codeml_wrapper_check = True
                            except:
                                pass
                        if codeml_wrapper_check:
                            wrapper_out, wrapper_error = codeml_wrapper_call.communicate()
                            if not wrapper_error:
                                for path, sub_dirs, file_list in os.walk(codeml_dir):
                                    if 'Omega' in path:
                                        taskfarm_log.write('cd {0}; codeml\n'.format(path))
                            else:
                                if 'number of sequences' in wrapper_error:
                                    error_string = wrapper_error.strip().split(': ')[1].split(',')[0]
                                    error_string = error_string[0].upper() + error_string[1:]
                                    print error_string
                                    error_log.write('Error running GenerateCodemlWorkspace.pl.\n')
                                    error_log.write('Error Reported: {0}.\n'.format(error_string))   
                                else:
                                    print 'Error running GenerateCodemlWorkspace.pl. Please check log file for details.'
                                    error_log.write('Error running GenerateCodemlWorkspace.pl.\n')
                                    error_log.write('{0}.\n'.format(wrapper_error))   
                        else:
                            print 'Error running GenerateCodemlWorkspace.pl, please confirm that the script and all modules are installed.'
                            error_log.write('Error running GenerateCodemlWorkspace.pl, please confirm that the script and all modules are installed.\n')
            else:
                with open(alignment_input.current_input) as test_data:
                    if test_data.readline().startswith('>'):
                        print 'Unaligned seqeunces detected in:{0}'.format(return_filename(alignment_input.current_input))
                        error_log.write('Unaligned seqeunces detected in:{0}\n'.format(return_filename(alignment_input.current_input)))
    if not alignment_found_check:
        print 'No alignment(s) specified by input option, please verify that each tree file is accompanied by an alignment'
    if not tree_found_check:
        print 'No tree file(s) specified by input option, please verify that each tree file is accompanied by an alignment'
    taskfarm_log.close()
    error_log.close()


### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### Function: vespa_subtrees:
### Details: Creates subtrees from newick tree data
def vespa_subtrees(input_files):
    print 'VESPA: Create SubTrees'
    import os, sys, dendropy, re
    from collections import defaultdict
    
    class nodes_picker_data(object):
        def __init__(self, current_tree, log_data):
            self.job_name = current_tree.strip()
            self.requested_command = log_data[0]
            self.original_tree = log_data[1]
            self.labeled_tree = ''
            self.edited_tree = log_data[2]
            self.outgroup_request = True
            self.outgroups = []
            
        def __str__(self):
            return '#' + self.job_name + '\nEdit_Method:' + self.requested_command + '\nOrignal_Tree:' + self.original_tree + '\nEdited_Tree:' + self.edited_tree + '\nOutgroup(s):' + ', '.join(self.outgroups) + '\n'

    def np_log_reader(np_log_file):
        from itertools import groupby
        return_trees = []
        tree_groups = (entry[1] for entry in groupby(open(np_log_file, 'rU'), lambda line: line.startswith('#')))
        for current_tree in tree_groups:
            return_trees.append(nodes_picker_data(current_tree.next().strip()[1:], [tree_lines.strip().split(':')[1] for tree_lines in tree_groups.next()]))
        return return_trees

    def select_subtree(tree_to_modify):
        unedited_tree = dendropy.Tree.get_from_string(tree_to_modify.labeled_tree,"newick")
        while not tree_to_modify.edited_tree:
            user_selected_node = raw_input('Please select a node for subtree creation: ')
            if user_selected_node in [internal_nodes.label for internal_nodes in unedited_tree.internal_nodes()]:
                user_confirm = raw_input('Node ' + user_selected_node + ' found. Please confirm (y / n): ')
                if user_confirm.lower().startswith('y'):
                    edited_node = unedited_tree.find_node_with_label(user_selected_node)
                    edited_leafs = [leaf_nodes.taxon for leaf_nodes in edited_node.leaf_iter()]
                    edited_leaf_taxa = [leaf_nodes.label for leaf_nodes in edited_leafs]
                    edited_tree = dendropy.Tree(unedited_tree)
                    edited_tree.retain_taxa(edited_leafs)
                    edited_tree_string = edited_tree.as_string(schema="newick", suppress_rooting=True).strip()
                    for label_match in [match.group() for match in re.finditer('\)\d+',edited_tree_string)]:
                        edited_tree_string = edited_tree_string.replace(label_match, ')')
                    tree_to_modify.edited_tree = edited_tree_string.strip()

    def remove_node(tree_to_modify):
        edited_tree = dendropy.Tree.get_from_string(tree_to_modify.labeled_tree,"newick")
        while not tree_to_modify.edited_tree:
            user_selected_node = raw_input('Please select a node for removal: ')
            if user_selected_node in [internal_nodes.label for internal_nodes in edited_tree.internal_nodes()]:
                user_confirm = raw_input('Node ' + user_selected_node + ' found. Please confirm (y / n): ')
                if user_confirm.lower().startswith('y'):
                    edited_node = edited_tree.find_node_with_label(user_selected_node)
                    edited_tree.prune_subtree(edited_node)
                    edited_tree_string = edited_tree.as_string(schema="newick", suppress_rooting=True).strip()
                    for label_match in [match.group() for match in re.finditer('\)\d+',edited_tree_string)]:
                        edited_tree_string = edited_tree_string.replace(label_match, ')')
                    tree_to_modify.edited_tree = edited_tree_string.strip()

    def remove_leaf(tree_to_modify):
        edited_tree = dendropy.Tree.get_from_string(tree_to_modify.labeled_tree,"newick")
        leaf_list = [str(taxa_ids).replace("'",'') for taxa_ids in edited_tree.taxon_namespace]
        while not tree_to_modify.edited_tree:
            user_selected_leafs = raw_input('Please select a leaf (taxa) for removal (if mutiple, seperate with comma): ')
            user_leafs = [current_leaf.strip() for current_leaf in user_selected_leafs.split(',')]
            if len(user_leafs) == len(list(set(user_leafs) & set(leaf_list))):
                user_confirm = raw_input('Leaf(s) ' + ', '.join(user_leafs) + ' found. Please confirm (y / n): ')
                if user_confirm.lower().startswith('y'):
                    remove_list = [edited_tree.find_node_with_taxon_label(current_taxa).taxon for current_taxa in leaf_list if current_taxa in user_leafs]
                    edited_tree.prune_taxa(remove_list)
                    edited_tree_string = edited_tree.as_string(schema="newick").replace("'",'')
                    for label_match in [match.group() for match in re.finditer('\)\d+',edited_tree_string)]:
                        edited_tree_string = edited_tree_string.replace(label_match, ')')
                    tree_to_modify.edited_tree = edited_tree_string.strip()

    def current_tree_user_request(requesting_job):
        
        def tree_labeler(unlabeled_tree):
            labeled_tree = list(unlabeled_tree.original_tree)
            for pos, matched_nodes in enumerate([match.end() for match in re.finditer('\)', unlabeled_tree.original_tree)][::-1]):
                labeled_tree.insert(matched_nodes,str(pos))
            unlabeled_tree.labeled_tree = ''.join(labeled_tree)
            
        def command_reqeust(command_job):
            command_dict = {'1':'Subtree Selection', '2':'Node Removal', '3':'Leaf (Taxa) Removal', '4':'Keep Original'}
            while not command_job.requested_command:
                user_request = raw_input('Please select an action ( 1 / 2 / 3 / 4 ): ')
                if user_request in ['1', '2', '3', '4']:
                    user_confirm = raw_input(command_dict[user_request] + '. Please confirm (y / n): ')
                    if user_confirm.lower().startswith('y'):
                        command_job.requested_command = command_dict[user_request]
                        if command_job.requested_command == 'Subtree Selection':
                            select_subtree(command_job)
                            command_job.job_status = 'Finished'
                        elif command_job.requested_command == 'Node Removal':
                            remove_node(command_job)
                            command_job.job_status = 'Finished'
                        elif command_job.requested_command == 'Leaf (Taxa) Removal':
                            remove_leaf(command_job)
                            command_job.job_status = 'Finished'
                        elif command_job.requested_command == 'Keep Original':
                            command_job.job_status = 'Finished'
                            command_job.edited_tree = command_job.original_tree
                        
        def outgroup_reqeust(command_job):
            taxa_tree = dendropy.Tree.get_from_string(command_job.original_tree,"newick")
            taxa_list = [str(taxa_ids).replace("'",'') for taxa_ids in taxa_tree.taxon_namespace]
            while command_job.outgroup_request and not command_job.outgroups:
                verify_outgroups = raw_input('Additional outgroup(s) required? (y / n): ')
                if 'y' in verify_outgroups.lower():
                    while not command_job.outgroups:
                        selected_outgroups = raw_input('Please indicate outgroup(s) to add (if mutiple, seperate with comma): ')
                        outgroup_list = [current_outgroup.strip() for current_outgroup in selected_outgroups.split(',')]
                        user_confirm = raw_input('Following outgroup(s) selected: ' + ', '.join(outgroup_list) + '. Please confirm (y / n): ')
                        if user_confirm.lower().startswith('y'):
                            if len(outgroup_list) == len(list(set(outgroup_list) & set(taxa_list))):
                                command_job.outgroups = outgroup_list
                            else:
                                print 'Unknown outgroup detected'
                else:
                    command_job.outgroup_request = False
                    
        tree_labeler(requesting_job)
        screen_tree = dendropy.Tree.get_from_string(requesting_job.labeled_tree,"newick")
        screen_tree.print_plot(show_internal_node_labels=True) 
        print 'Current Tree: ' + requesting_job.job_name + '\n'
        print 'Possible actions\n____________________\n1. Subtree Selection\n2. Node Removal\n3. Leaf (Taxa) Removal\n4. Keep Original\n'
        command_reqeust(requesting_job)
        if requesting_job.requested_command != 'Keep Original': 
            outgroup_reqeust(requesting_job)
    
    global bme_sequence_database_location, bme_output_directory
    if not bme_sequence_database_location:
        print 'Sequence database not found'
        sys.exit()
    
    check_output_dir('Subtrees')
        
    current_log_data = []
    if os.path.isfile(bme_subtree_log_file):
        current_log_data = np_log_reader(bme_subtree_log_file)
        completed_jobs = [completed_jobs.job_name for completed_jobs in current_log_data]
        updated_log = open(bme_subtree_log_file,'a')
        for tree_files in input_files:
            if tree_files.current_input not in completed_jobs:
                newick_tree = dendropy.Tree.get_from_path(tree_files.current_input,"newick").as_string(schema="newick", suppress_rooting=True).strip()
                current_tree = nodes_picker_data(tree_files.current_input,['', newick_tree, ''])
                current_tree_user_request(current_tree)
                current_log_data.append(current_tree)
                updated_log.write(str(current_tree))
        updated_log.close()
    else:
        new_log = open(bme_subtree_log_file,'a')
        for tree_files in input_files:
            newick_tree = dendropy.Tree.get_from_path(tree_files.current_input,"newick").as_string(schema="newick", suppress_rooting=True).strip()
            current_tree = nodes_picker_data(tree_files.current_input, ['', newick_tree, ''])
            current_tree_user_request(current_tree)
            current_log_data.append(current_tree)
            new_log.write(str(current_tree))
        new_log.close()
    
    sequence_table = defaultdict(list)
    for tree_data in current_log_data:
        sequence_filename = '{0}/{1}.{2}'.format(bme_output_directory, return_filename_wo_ext(tree_data.job_name), return_extension(bme_sequence_database_location))
        log_tree = dendropy.Tree.get_from_string(tree_data.edited_tree,"newick")
        for taxa in [str(taxa_ids).replace("'",'') for taxa_ids in log_tree.taxon_namespace]:
            sequence_table[taxa.strip().replace('#1','')].append(sequence_filename)
        if tree_data.outgroups:
            for current_outgroup in tree_data.outgroups:
                sequence_table[current_outgroup.strip().replace('#1','')].append(sequence_filename)
        try:
            with open(sequence_filename): os.remove(sequence_filename)
        except IOError:
            pass

    for working_sequence in sequence_reader(bme_sequence_database_location).read():
        for seqeunce_keys in sequence_table.keys():
            if seqeunce_keys in working_sequence.header:
                for current_sequence_file in sequence_table[seqeunce_keys]: 
                    sequence_output = open(current_sequence_file, 'a')
                    sequence_output.write(str(working_sequence) + '\n')
                    sequence_output.close()

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
###  ____          _   _____  _                    
### |___ \        | | |  __ \| |                   
###   __) |_ __ __| | | |__) | |__   __ _ ___  ___ 
###  |__ <| '__/ _` | |  ___/| '_ \ / _` / __|/ _ \
###  ___) | | | (_| | | |    | | | | (_| \__ \  __/
### |____/|_|  \__,_| |_|    |_| |_|\__,_|___/\___|
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### Function: vespa_metAl_compare:
### Details: Compares two sets of files
def vespa_metAl_compare (input_files):
    print 'VESPA: metAl compare'
    from collections import defaultdict
    from subprocess import Popen, PIPE
    import sys, shutil
    global bme_metAl_compare_files, bme_metAl_compare_dir
    
    def check_noRMD(location_list):
        noRMD_values = {}
        for pos, alingment_files in enumerate(location_list):
            noRMD_program_check, noRMD_program_bin = False, False
            try:
                noRMD_program_call = Popen(['normd', alingment_files], stdout=PIPE, stderr=PIPE)
                noRMD_program_check = True
            except:
                noRMD_program_bin = True
            if noRMD_program_bin:
                try:
                    noRMD_program_call = Popen(['./normd', alingment_files], stdout=PIPE, stderr=PIPE)
                    noRMD_program_check = True
                except:
                    pass
            if noRMD_program_check:
                noRMD_output, noRMD_error = noRMD_program_call.communicate()
                if not noRMD_error:
                    noRMD_values[location_list[pos]] = float(noRMD_output.strip())
                else:
                    print 'Error detected with noRMD. Please confirm the program is correctly compiled'
                    sys.exit(0)
            else:
                print 'Error running noRMD. Please confirm the program is installed'
                sys.exit(0)
        return noRMD_values
    
    def scoreMetAl (sent_compare_data, sent_output_dir):
        global bme_metAl_cutoff
        for alignment_ID, alignment_locations in sent_compare_data.items():
            if len(alignment_locations) == 2:
                metal_program_check, metal_program_bin = False, False
                try:
                    metAl_program_call = Popen(['metal', alignment_locations[0], alignment_locations[1]], stdout=PIPE, stderr=PIPE)
                    metal_program_check = True
                except:
                    metal_program_bin = True
                if metal_program_bin:
                    try:
                        metAl_program_call = Popen(['./metal', alignment_locations[0], alignment_locations[1]], stdout=PIPE, stderr=PIPE)
                        metal_program_check = True
                    except:
                        pass
                if metal_program_check:
                #metAl_program_call = Popen(['metal', alignment_locations[0], alignment_locations[1]], stdout=PIPE, stderr=PIPE)
                    metAl_output, metAl_error = metAl_program_call.communicate()
                    if not metAl_error:
                        split_metAl = metAl_output.strip().split('= ')
                        if float(split_metAl[1]) < bme_metAl_cutoff:
                            metal_compare_results.write('{0},{1},{2},{3},{4},{5}\n'.format(alignment_ID, split_metAl[0].strip(), split_metAl[1], 'null', 'null', alignment_locations[0]))
                            shutil.copy(alignment_locations[0], sent_output_dir)
                        else:
                            returned_values_dict = check_noRMD(alignment_locations)
                            alignment_compare = returned_values_dict.keys()
                            if returned_values_dict[alignment_compare[0]] == returned_values_dict[alignment_compare[1]]:
                                metal_compare_results.write('{0},{1},{2},{3},{4},{5}\n'.format(alignment_ID, split_metAl[0].strip(), split_metAl[1], returned_values_dict[alignment_compare[0]], returned_values_dict[alignment_compare[1]], alignment_locations[0]))
                                shutil.copy(alignment_locations[0], sent_output_dir)
                           
                            elif returned_values_dict[alignment_compare[0]] > returned_values_dict[alignment_compare[1]]:
                                if alignment_compare[0] == alignment_locations[0]:
                                    metal_compare_results.write('{0},{1},{2},{3},{4},{5}\n'.format(alignment_ID, split_metAl[0].strip(), split_metAl[1], returned_values_dict[alignment_compare[0]], returned_values_dict[alignment_compare[1]], alignment_compare[0]))
                                if alignment_compare[1] == alignment_locations[0]:
                                    metal_compare_results.write('{0},{1},{2},{3},{4},{5}\n'.format(alignment_ID, split_metAl[0].strip(), split_metAl[1], returned_values_dict[alignment_compare[1]], returned_values_dict[alignment_compare[0]], alignment_compare[0]))
                                shutil.copy(alignment_compare[0], sent_output_dir)
                        
                            elif returned_values_dict[alignment_compare[0]] < returned_values_dict[alignment_compare[1]]:
                                if alignment_compare[0] == alignment_locations[0]:
                                    metal_compare_results.write('{0},{1},{2},{3},{4},{5}\n'.format(alignment_ID, split_metAl[0].strip(), split_metAl[1], returned_values_dict[alignment_compare[0]], returned_values_dict[alignment_compare[1]], alignment_compare[1]))
                                if alignment_compare[1] == alignment_locations[0]:
                                    metal_compare_results.write('{0},{1},{2},{3},{4},{5}\n'.format(alignment_ID, split_metAl[0].strip(), split_metAl[1], returned_values_dict[alignment_compare[1]], returned_values_dict[alignment_compare[0]], alignment_compare[1]))
                                shutil.copy(alignment_compare[1], sent_output_dir)
                    else:
                        print 'Error detected with metAl. Please confirm the program is correctly compiled'
                        sys.exit(0)
                else:
                        print 'Error running metAl. Please confirm the program is installed'
                        sys.exit(0)
            else:
                print 'Cannot find comparison alignment for: {0}'.format(alignment_ID)

    metal_compare_results = create_unique_file('metAl_compare.csv')
    metal_compare_results.write('{0},{1},{2},{3},{4},{5}\n'.format('Alignment_ID', 'metAL_d_pos', 'metAL_score', 'noRMD_input', 'noRMD_compare', 'Selected_Alignment'))
    metAl_output_dir, metAl_output_filename, metAl_output = '', '', ''
    compare_dict = defaultdict(list)

    for sequence_input in input_files:
        if verify_alignment(sequence_input.current_input):
            (metAl_output_dir, metAl_output_filename, metAl_output) = check_output(sequence_input, 'metAl_compare')
            compare_dict[return_filename_wo_ext(sequence_input.current_input)].append(sequence_input.current_input)
        else:
            print return_filename(sequence_input.current_input) + ': Not an alignment file'
            
    for sequence_compare_input in bme_metAl_compare_files:
        if verify_alignment(sequence_compare_input):
            compare_dict[return_filename_wo_ext(sequence_compare_input)].append(sequence_compare_input)
        else:
            print return_filename(sequence_compare_input) + ': Not an alignment file'
    
    scoreMetAl(compare_dict, metAl_output_dir)
    metal_compare_results.close()

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### Function: vespa_setup_prottest:
### Details: Creates the needed inputfile files for a prottest run
def vespa_setup_prottest(input_files):
    print 'VESPA: ProtTest Setup'
    import shutil
    prottest_file = create_unique_file('setup_prottest_taskfarm')
    for sequence_input in input_files:
        if verify_alignment(sequence_input.current_input):
            (prottest_output_dir, prottest_output_filename, prottest_output) = check_output(sequence_input, 'ProtTest_Setup')
            shutil.copy(sequence_input.current_input, prottest_output_dir)
            prottest_file.write('java -jar prottest.jar -i ' + prottest_output + ' -o ' + prottest_output_dir + '/' + remove_extension(prottest_output_filename) + '.models -all-distributions\n')
        else:
            print return_filename(sequence_input) + ': Not an alignment file' 
    prottest_file.close()

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### Function: vespa_prottest_reader:
### Details: Reads prottest output and creates generic model output (best/supported by MrBayes)
def vespa_prottest_reader (input_files):
    print 'VESPA: ProtTest Results Reader'
        
    def protest_verify(protest_output):
        verify_model = False
        alignment_file = ''
        with open(protest_output) as model_file:
            supported_by_mrbayes = ['Dayhoff', 'JTT', 'Blosum62', 'VT', 'WAG']
            for check_model_file in [model_file.next() for x in xrange(3)]:
                if 'ProtTest' in check_model_file:
                    verify_model = True
                    break
            if verify_model:
                for model_lines in model_file:
                    if 'Alignment file' in model_lines:
                        alignment_file = return_filename(model_lines.split(':')[-1].strip())
        return alignment_file
    
    def protest_output_reader(protest_output):
        return_best_model, return_best_supported_model = ('', '')
        with open(protest_output) as model_file:
            supported_by_mrbayes = ['Dayhoff', 'JTT', 'Blosum62', 'VT', 'WAG']
            data_block_test, end_of_block = (False, False)
            for model_lines in model_file:
                if data_block_test and ('-' * 75) == model_lines.strip():
                    end_of_block = True
                elif data_block_test and not end_of_block:
                    split_model = model_lines.strip().split()
                    if not return_best_model:
                        return_best_model = split_model[0].split('+')
                    if split_model[0].split('+')[0] in supported_by_mrbayes and not return_best_supported_model:
                        return_best_supported_model = split_model[0].split('+')
                elif ('-' * 75) == model_lines.strip():
                    data_block_test = True
        return return_best_model, return_best_supported_model

    report_best_models = create_unique_file('prottest_reader.best_models')
    report_best_supported_models = create_unique_file('prottest_reader.best_supported_models')

    for model_input in input_files:
        best_model, best_supported_model = ('', '')
        alignment_input = protest_verify(model_input.current_input)
        if alignment_input:
            best_model, best_supported_model = protest_output_reader(model_input.current_input)
            report_best_models.write(alignment_input + ',' + '+'.join(best_model) + '\n')
            report_best_supported_models.write(alignment_input + ',' + '+'.join(best_supported_model) + '\n')
    report_best_models.close()
    report_best_supported_models.close()
    
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### Function: vespa_setup_mrbayes:
### Details: Reads create Nexus files for MrBayes given fasta files and prottest supported file
def vespa_setup_mrbayes (input_files):
    print 'VESPA: MrBayes Setup'
    global bme_supported_model_list, bme_mrbayes_mcmc_gen, bme_mrbayes_mcmc_chains, bme_mrbayes_mcmc_temp, bme_mrbayes_mcmc_burnin
    import os, sys
    header_warning = False

    def format_nexus (alignment_list, nexus_filename):
        sequence_for_info = alignment_list[0]
        sequence_for_info.seq_type()
        sequence_length = len(sequence_for_info)
        nexus_file = open(nexus_filename, 'w')
        nexus_file.write('\n'.join(['#NEXUS', '', 'BEGIN DATA;', 'DIMENSIONS NTAX=' + str(len(alignment_list)) + ' NCHAR=' +
              str(sequence_length) + ';', 'FORMAT DATATYPE=' + sequence_for_info.type + ' MISSING=- INTERLEAVE;', '', 'MATRIX']) + '\n')
        current_sequence_position = 0
        while current_sequence_position < sequence_length:
            for alignment_seqeunce in alignment_list:
                alignment_header = alignment_seqeunce.header.strip()[1:] 
                nexus_file.write(alignment_header + (' ' * (22 - len(alignment_header))) + ' '.join([alignment_seqeunce.sequence[seqeunce_block:seqeunce_block + 20] for seqeunce_block in range(current_sequence_position, current_sequence_position + 100, 20)]) + '\n')
            current_sequence_position += 100
            nexus_file.write('\n')
        nexus_file.write('\n'.join([';', 'END;\n']))
        nexus_file.close()
    
    
    def convert_for_mrbayes (orginal_model):
        if '+' in orginal_model:
            command_list = orginal_model.split('+')
        else:
            command_list = [orginal_model,'']
        convert_model = {'Dayhoff':'dayhoff', 'JTT':'jones', 'Blosum62':'blosum', 'VT':'vt','WAG':'wag'}
        convert_options = {'I':'propinv', 'G':'gamma', 'IG':'invgamma', '':'equal'}
        return convert_model[command_list[0]], convert_options[''.join(command_list[1:])]
    
    
    if bme_supported_model_list:
        supported_model_dict = {}
        with open(bme_supported_model_list) as model_data:
            for model_lines in model_data:
                model_split = model_lines.strip().split(',')
                supported_model_dict[model_split[0]] = model_split[1]
        
        for sequence_input in input_files:
            if verify_alignment(sequence_input.current_input):
                if supported_model_dict.has_key(return_filename(sequence_input.current_input)):
                    mrbayes_lines = ['\nbegin mrbayes;', 'log start filename=Logs/' + return_filename_wo_ext(sequence_input.current_input) + '.log replace;', 'set autoclose=yes;']
                    (mrbayes_output_dir, mrbayes_output_filename, mrbayes_output) = check_output(sequence_input, 'MrBayes_Setup')
                    
                    model_input, rate_input = convert_for_mrbayes(supported_model_dict[return_filename(sequence_input.current_input)])
                    mrbayes_lines.extend(['lset applyto=(all) rates=' + rate_input + ';', 'prset aamodelpr=fixed(' + model_input + ');',
                                        'mcmcp ngen={0} printfreq=2000 samplefreq=200 nchains={1} temp={2} savebrlens=yes relburnin=yes burninfrac={3};'.format(bme_mrbayes_mcmc_gen, bme_mrbayes_mcmc_chains, bme_mrbayes_mcmc_temp, bme_mrbayes_mcmc_burnin),
                                        'mcmc;', 'sumt;', 'sump;', 'log stop;', 'end;'])

                    nexus_sequence_input = []
                    for working_sequence in sequence_reader(sequence_input.current_input).read():
                        if len(working_sequence.header) > 22:
                            if not header_warning:
                                header_warning = True
                                print 'Warning: Sequence headers too long for NEXUS format - Editing headers for length (Manual editing beforehand is recommended)'
                            working_sequence.header = '{0}\n'.format(working_sequence.header[:22])
                        nexus_sequence_input.append(working_sequence)
                    
                    coverted_filename = '{0}.nex'.format(remove_extension(mrbayes_output))
                    format_nexus(nexus_sequence_input, coverted_filename)
                    append_mrbayes_block = open(coverted_filename, 'a')
                    append_mrbayes_block.write('\n'.join(mrbayes_lines))
                    append_mrbayes_block.close()
                        
    else:
        print 'No ProtTest model table specified. Please specify using -model_table='

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
###  ___            _   _____  _                    
### |__ \          | | |  __ \| |                   
###    ) |_ __   __| | | |__) | |__   __ _ ___  ___ 
###   / /| '_ \ / _` | |  ___/| '_ \ / _` / __|/ _ \
###  / /_| | | | (_| | | |    | | | | (_| \__ \  __/
### |____|_| |_|\__,_| |_|    |_| |_|\__,_|___/\___|                                               
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### Functions: vespa_setup_reciprocal_input
### Details: creates an input database for a reciprocal run
def vespa_setup_reciprocal_input(input_files):
    global bme_sequence_database_location
    csv_list = []
    for similarity_file in input_files:
        with open(similarity_file.current_input) as similarity_data:
            for similarity_lines in similarity_data:
                if similarity_lines.strip().split()[1] not in csv_list:
                    csv_list.append(similarity_lines.strip().split()[1])
    reciprocal_output = open(create_unique_file('Reciprocal_Input.' + bme_sequence_database_location.split('.')[-1]), 'w')
    for working_sequence in sequence_reader(bme_sequence_database_location).read():
        for csv_entries in csv_list:
            if csv_entries in working_sequence.header:
                reciprocal_output.write(str(working_sequence) + '\n')
    reciprocal_output.close()
    
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### Functions: create_similarity_groups
### Details: create_similarity_groups: Creates connected components and creates sequence files
def create_similarity_groups(graph_data):
    def merge_groups(merge_data):
        merge_check = True
        while merge_check:
            merge_check = False
            current_mergers = []
            while merge_data:
                query, subject_list = merge_data[0], merge_data[1:]
                merge_data = []
                for subject in subject_list:
                    if subject.isdisjoint(query):
                        merge_data.append(subject)
                    else:
                        merge_check = True
                        query.update(subject)
                current_mergers.append(query)
            merge_data = current_mergers
        return merge_data
    
    import os
    from collections import defaultdict    
    
    global bme_output_directory, bme_sequence_database_location
    sequence_table = defaultdict(str)
    check_output_dir('Similarity_Groups')
    
    merged_graph_data = merge_groups(graph_data)
    total_number_of_files = len(str(len(merged_graph_data)))
    for file_counter, sequence_list in enumerate(merged_graph_data):
        similarity_filename = bme_output_directory + '/similarity_group_' + ('0' * (total_number_of_files - len(str(file_counter)))) + str(file_counter) + '.fasta'
        for sequences in sequence_list:
            sequence_table[sequences] = similarity_filename
        try:
            with open(similarity_filename): os.remove(similarity_filename)
        except IOError:
            pass

    for working_sequence in sequence_reader(bme_sequence_database_location).read():
        if sequence_table.has_key(working_sequence.header[1:].strip()):
            sequence_output = open(sequence_table.pop(working_sequence.header[1:].strip()), 'a')
            sequence_output.write(str(working_sequence) + '\n')
            sequence_output.close()
    
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### Functions: assign_edges
### Details: assign_edges: Checks for thresholds and assigns edges using the correct format template
def assign_connections (assignment_graph, assign_data):
    global bme_similarity_e_value_cutoff, bme_similarity_alignment_length_cutoff, bme_similarity_percent_identity_cutoff, bme_similarity_data_format
    global blast_alignment_length_warn, hmmer_percent_identity_warn
    pass_thresholds = True
    assign_e_value, assign_percent_identity, assign_alignment_length = (False, False, False)
    
    if bme_similarity_data_format == 'blast':
        assign_e_value, assign_percent_identity = (float(assign_data[10]), float(assign_data[2]))
    if bme_similarity_data_format == 'hmmer':
        assign_e_value = float(assign_data[6])
        assign_alignment_length = float(assign_data[2]) / float(assign_data[5])
    
    if bme_similarity_e_value_cutoff:
        if float(bme_similarity_e_value_cutoff) < assign_e_value:
            pass_thresholds = False
    if bme_similarity_percent_identity_cutoff:
        if bme_similarity_data_format == 'hmmer':
            if not hmmer_percent_identity_warn:
                hmmer_percent_identity_warn = True
                print 'Percent Identity: HMMER does not compute use percent identity - command ignored'
        else:
            if float(bme_similarity_percent_identity_cutoff) > assign_percent_identity:
                pass_thresholds = False
    if bme_similarity_alignment_length_cutoff:
        if assign_alignment_length:
            if float(bme_similarity_alignment_length_cutoff) > assign_alignment_length:
                pass_thresholds = False
        else:
            if not blast_alignment_length_warn:
                blast_alignment_length_warn = True
                print 'Alignment length: cannot compute from BLAST output alone'
    if pass_thresholds:
        if assignment_graph.has_key((assign_data[1],assign_data[0])):
            assignment_graph[(assign_data[1],assign_data[0])] = True
        else:
            assignment_graph[(assign_data[0],assign_data[1])] = False
    
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### Functions: vespa_best_reciprocal_similarity_groups
### Details: Identifies best reciprocals between species
def vespa_best_reciprocal_similarity_groups(input_files):
    
    def format_splitter(unsplit_similarity_line, similarity_format):
        return_line = []
        if similarity_format == 'blast':
            return_line = unsplit_similarity_line.strip().split('\t')
        if similarity_format == 'hmmer':
            return_line = unsplit_similarity_line.strip().split()
        return return_line
    
    def return_compare_data(umcompared_line, similarity_format):
        return_query_sequence, return_query_species  = '', ''
        return_subject_sequence, return_subject_species  = '', ''
        return_e_value = 1.0
        if similarity_format == 'blast':
            return_query_sequence, return_subject_sequence = umcompared_line[0], umcompared_line[1]
            return_query_species, return_subject_species = return_query_sequence.split('|')[0], return_subject_sequence.split('|')[0]
            return_e_value = float(umcompared_line[10])
        if similarity_format == 'hmmer':
            pass
        return  return_query_sequence, return_query_species, return_subject_sequence, return_subject_species, return_e_value
    
    from collections import defaultdict
    global bme_similarity_data_format
    
    print 'VESPA: Best-Reciprocal Groups'
    reciprocality_table = defaultdict(dict)
    reciprocality_graph = {}
    for similarity_file in input_files:
        with open(similarity_file.current_input) as similarity_data:
            for similarity_lines in similarity_data:
                split_lines = format_splitter(similarity_lines, bme_similarity_data_format)
                query_sequence, query_species, subject_sequence, subject_species, e_value = return_compare_data(split_lines, bme_similarity_data_format)
                if query_species != subject_species:
                    if reciprocality_table[query_sequence].has_key(subject_species):
                        if e_value < reciprocality_table[query_sequence][subject_species][1]:
                            reciprocality_table[query_sequence][subject_species] = [subject_sequence, e_value]
                    else:
                        reciprocality_table[query_sequence][subject_species] = [subject_sequence, e_value]
    for similarity_file in input_files:
        with open(similarity_file.current_input) as similarity_data:
            for similarity_lines in similarity_data:
                split_lines = format_splitter(similarity_lines, bme_similarity_data_format)
                query_sequence, query_species, subject_sequence, subject_species, e_value = return_compare_data(split_lines, bme_similarity_data_format)
                if reciprocality_table[query_sequence].has_key(subject_species):
                    if reciprocality_table[query_sequence][subject_species][0] == subject_sequence:
                        assign_connections(reciprocality_graph, split_lines)
    
    sub_graph = []
    for connection, reciprocality_confirmation in reciprocality_graph.items():
        if reciprocality_confirmation:
            sub_graph.append(set(connection))
    reciprocality_graph.clear()
    create_similarity_groups(sub_graph)
    
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### Functions: vespa_similarity_groups
### Details: Identifies either simple or reciprocal connections within file
def vespa_similarity_groups(input_files, reciprocality_check):
    
    def format_splitter(unsplit_similarity_line, similarity_format):
        return_line = []
        if similarity_format == 'blast':
            return_line = unsplit_similarity_line.strip().split('\t')
        if similarity_format == 'hmmer':
            return_line = unsplit_similarity_line.strip().split()
        return return_line
    
    def test_selfhit(check_line, similarity_format):
        return_check = True
        if similarity_format == 'blast':
            if check_line[0] == check_line[1]:
                return_check = False
        if similarity_format == 'hmmer':
            pass
        return return_check
    
    if reciprocality_check:
        print 'VESPA: Reciprocal Groups'
    else:
        print 'VESPA: Similarity Groups'
        
    global bme_similarity_data_format
    reciprocality_graph = {}
    for similarity_file in input_files:
        with open(similarity_file.current_input) as similarity_data:
            for similarity_lines in similarity_data:
                split_lines = format_splitter(similarity_lines, bme_similarity_data_format)
                if test_selfhit(split_lines, bme_similarity_data_format):     
                    assign_connections(reciprocality_graph, split_lines)
    
    if not reciprocality_check:
        graph = []
        for connection in reciprocality_graph.keys():
            graph.append(set(connection))
        reciprocality_graph.clear()
        create_similarity_groups(graph)
    else:
        sub_graph = []
        for connection, reciprocality_confirmation in reciprocality_graph.items():
            if reciprocality_confirmation:
                sub_graph.append(set(connection))
        reciprocality_graph.clear()
        create_similarity_groups(sub_graph)



### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
###  __     _     _____  _                    
### /_ |   | |   |  __ \| |                   
###  | |___| |_  | |__) | |__   __ _ ___  ___ 
###  | / __| __| |  ___/| '_ \ / _` / __|/ _ \
###  | \__ \ |_  | |    | | | | (_| \__ \  __/
###  |_|___/\__| |_|    |_| |_|\__,_|___/\___|
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### Function: vespa_clean:
### Details: Removes sequences from the specified filepath that are not divisible by three
def vespa_clean (input_files):
    print 'VESPA: Cleaning sequences'
    removed_in_cleanfile = create_unique_file('cleaned_genes_removed.log')
    global bme_remove_internal_stop, bme_label_with_filename, bme_infer_labels_ensembl
    for sequence_input in input_files:
        removed_in_cleanfile.write('Cleaning File: {0}\n'.format(return_filename(sequence_input.current_input)))
        (cleaned_output_dir, cleaned_output_filename, cleaned_ouput) = check_output(sequence_input, 'Cleaned')
        bme_clean_file = open(cleaned_ouput, 'w')
        for working_sequence in sequence_reader(sequence_input.current_input).read():
            if len(working_sequence) % 3 == 0:
                remove_check = False
                if bme_remove_internal_stop:
                    working_sequence.type = 'DNA'
                    if working_sequence.internal_stop():
                        remove_check = True
                        removed_in_cleanfile.write('Gene removed - Internal stop-codon found: {0}'.format(working_sequence.header))
                if not remove_check:
                    if bme_label_with_filename:
                        working_sequence.header = '>{0}|{1}'.format(return_filename_wo_ext(sequence_input.current_input), working_sequence.header[1:])
                    if bme_infer_labels_ensembl:
                        species_header = ensembl_infer(working_sequence.header)
                        if species_header:
                            working_sequence.header = '>{0}|{1}'.format(species_header, working_sequence.header[1:]) 
                    bme_clean_file.write(str(working_sequence) + '\n')
            else:
                removed_in_cleanfile.write('Gene removed - Abnormal sequence length: {0}'.format(working_sequence.header))
        bme_clean_file.close()
    removed_in_cleanfile.close()
    
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### Function: vespa_clean_ensembl:
### Details: Cleans ensembl genome, returns longest transcripts divisible by three
def vespa_clean_ensembl (input_files):
    print 'VESPA: Cleaning ENSEBML sequences'
    global bme_remove_internal_stop, bme_label_with_filename, bme_infer_labels_ensembl
    removed_in_cleanfile = create_unique_file('cleaned_ensembl_removed.log')
    from collections import defaultdict
    for sequence_input in input_files:
        removed_in_cleanfile.write('Cleaning File: {0}\n'.format(return_filename(sequence_input.current_input)))
        (cleaned_output_dir, cleaned_output_filename, cleaned_output) = check_output(sequence_input, 'Cleaned')
        bme_clean_file = open(cleaned_output, 'w')
        geneDict = defaultdict(list)
        geneKeys = []
        for working_sequence in sequence_reader(sequence_input.current_input).read():
            geneData = working_sequence.header.strip().split('|')
            if len(working_sequence) % 3 == 0:
                remove_check = False
                if bme_remove_internal_stop:
                    working_sequence.type = 'DNA'
                    if working_sequence.internal_stop():
                        remove_check = True
                        removed_in_cleanfile.write('Gene removed - Internal stop codon found: {0}'.format(working_sequence.header))
                if not remove_check:
                    if geneDict.has_key(geneData[0]):
                        if geneDict[geneData[0]][1] < len(working_sequence):
                            removed_in_cleanfile.write('Gene removed - Longer transcript found: {0}'.format(geneDict[geneData[0]][0]))
                            geneDict[geneData[0]] = [working_sequence.header.strip(), len(working_sequence)]
                        else:
                            removed_in_cleanfile.write('Gene removed - Longer transcript found: {0}'.format(working_sequence.header))
                    else:
                        geneDict[geneData[0]] = [working_sequence.header.strip(), len(working_sequence)]
            else:
                removed_in_cleanfile.write('Gene removed - Abnormal sequence length: {0}'.format(working_sequence.header))
        geneKeys = [longest_gene[0] for longest_gene in geneDict.values()]  
        for working_sequence in sequence_reader(sequence_input.current_input).read():
            if working_sequence.header.strip() in geneKeys:
                if bme_label_with_filename:
                    working_sequence.header = '>{0}|{1}'.format(return_filename_wo_ext(sequence_input.current_input), working_sequence.header[1:])
                if bme_infer_labels_ensembl:
                    species_header = ensembl_infer(working_sequence.header)
                    if species_header:
                        working_sequence.header = '>{0}|{1}'.format(species_header, working_sequence.header[1:])
                bme_clean_file.write(str(working_sequence) + '\n')
        bme_clean_file.close()
    removed_in_cleanfile.close()

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### Function: vespa_reverse_complement:
### Details: Returns the reverse complement of the sequence in the specified filepath
def vespa_reverse_complement (input_files):
    print 'VESPA: Reverse Complementing Sequences'
    global bme_label_with_filename, bme_infer_labels_ensembl
    for sequence_input in input_files:
        (reversed_output_dir, reversed_output_filename, reversed_output) = check_output(sequence_input, 'RevComp')  
        bme_revcomp_file = open(reversed_output, 'w')
        for working_sequence in sequence_reader(sequence_input.current_input).read():
            if bme_label_with_filename:
                    working_sequence.header = '>{0}|{1}'.format(return_filename_wo_ext(sequence_input.current_input), working_sequence.header[1:])
            if bme_infer_labels_ensembl:
                species_header = ensembl_infer(working_sequence.header)
                if species_header:
                    working_sequence.header = '>{0}|{1}'.format(species_header, working_sequence.header[1:])
            bme_revcomp_file.write(str(working_sequence.seq_revcomp()) + '\n')
        bme_revcomp_file.close()

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### Function: vespa_translate:
### Details: Translate the specified filepath from DNA to Protien
def vespa_translate (input_files):
    print 'VESPA: Translating Sequences'
    global bme_remove_internal_stop, bme_remove_terminal_stop, bme_label_with_filename, bme_infer_labels_ensembl
    clean_warn = False
    removed_in_transfile = create_unique_file('translated_genes_removed.log')
    for sequence_input in input_files:
        removed_in_transfile.write('Translating File: {0}\n'.format(sequence_input.current_input))
        (translated_output_dir, translated_output_filename, translated_output) = check_output(sequence_input, 'Translated')
        bme_translate_file = open(translated_output, 'w')
        for working_sequence in sequence_reader(sequence_input.current_input).read():
            if len(working_sequence) % 3 == 0:
                working_sequence.seq_translate()
                remove_check = False
                if bme_remove_internal_stop:
                    working_sequence.type = 'protein'
                    if working_sequence.internal_stop():
                        removed_in_transfile.write('Gene removed - Internal stop codon found: {0}'.format(working_sequence.header))
                        remove_check = True
                if not remove_check:
                    if bme_remove_terminal_stop:
                        if working_sequence.sequence[-1] == '*':
                            working_sequence.sequence = working_sequence.sequence[:-1]
                    if bme_label_with_filename:
                        working_sequence.header = '>{0}|{1}'.format(return_filename_wo_ext(sequence_input.current_input), working_sequence.header[1:])
                    if bme_infer_labels_ensembl:
                        species_header = ensembl_infer(working_sequence.header)
                        if species_header:
                            working_sequence.header = '>{0}|{1}'.format(species_header, working_sequence.header[1:])
                    bme_translate_file.write(str(working_sequence) + '\n')
            else:
                if not clean_warn:
                    clean_warn= True
                    print 'Warning: Abnormal sequence length detected. Please confirm seqeunces are DNA and have been cleaned'
                removed_in_transfile.write('Gene removed - Abnormal sequence length: {0}'.format(working_sequence.header))
        bme_translate_file.close()
    removed_in_transfile.close()

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### Function: vespa_create_database:
### Details: Returns a database of all sequences in the specified filepath
def vespa_create_database (input_files):
    def assign_database_filename(assign_filename):
        if assign_filename:
            return create_unique_file(assign_filename)
        else:
            return create_unique_file('database.fas')  
    
    print 'VESPA: Creating Database'
    global bme_format_blast_database, bme_output_filename
    import subprocess
    sequence_type = ''
    database_created = False
    
    for sequence_input in input_files:
        if not database_created:
            bme_database_file = assign_database_filename(bme_output_filename)
            database_created = True
        for working_sequence in sequence_reader(sequence_input.current_input).read():
            bme_database_file.write(str(working_sequence) + '\n')
            if not sequence_type:
               working_sequence.seq_type()
               sequence_type = working_sequence.type
    bme_database_file.close()

    if bme_format_blast_database:
        type_convert = {'protein':'prot','DNA':'nucl'}
        try:
            blast_test = subprocess.Popen(['makeblastdb', '-dbtype', type_convert[sequence_type], '-in', bme_database_file.name], stdout=subprocess.PIPE,stderr=subprocess.PIPE, stdin=subprocess.PIPE)
            blast_out, blast_error = blast_test.communicate()
            if not blast_error:
                print 'VESPA: Formatting BLAST Database'
            else:
                print 'Error with makeblastdb function. Aborting format'
        except:
            print 'Cannot locate makeblastdb function. Aborting format'
    
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### Function: vespa_individual_sequences:
### Details: Returns single files of sequences in the specified filepath
def vespa_individual_sequences (input_files):
    def return_sequence_filename(sequence_header):
        return_seq_filename = ''
        for header_characters in sequence_header:
            if return_seq_filename:
                if header_characters == '|':
                    return_seq_filename += '_'
                elif header_characters == '_':
                    return_seq_filename += header_characters
                elif header_characters == ' ':
                    pass
                elif not header_characters.isalnum():
                    break
            if header_characters.isalnum():
                return_seq_filename += header_characters
        return return_seq_filename
        
    import os
    print 'VESPA: Creating Individual Sequences'
    global bme_output_directory
    directory_created = False
    
    for sequence_input in input_files:
        if not directory_created:
            (individual_output_dir, individual_output_filename, individual_output) = check_output(sequence_input, 'Individual')
            check_output_dir(individual_output_dir)
            directory_created = True
        for working_sequence in sequence_reader(sequence_input.current_input).read():
            individual_filename = '{0}/{1}'.format(bme_output_directory,return_sequence_filename(working_sequence.header.strip()))
            if '.' in sequence_input.current_input:
                individual_filename += '.{0}'.format(sequence_input.current_input.split('.')[-1])
            bme_individual_file = open(individual_filename, 'w')
            bme_individual_file.write(str(working_sequence) + '\n')
            bme_individual_file.close()

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### Function: vespa_split_in_groups:
### Details: Returns files that each contain multiple sequences
def vespa_split_in_groups (input_files):
    import os
    print 'VESPA: Creating sequence groups'
    global bme_split_number_in_groups, bme_output_directory
    total_sequences, total_files, sequence_counter, group_counter = 0, 0, 0, 0
    
    for sequence_input in input_files:
        with open(sequence_input.current_input) as data_for_totals:
            for lines_for_totals in data_for_totals:
                if '>' in lines_for_totals:
                    total_sequences += 1
    
    total_files = len(str(total_sequences/bme_split_number_in_groups))
    initial_loop = True
    
    print bme_split_number_in_groups
    
    for sequence_input in input_files:
        for working_sequence in sequence_reader(sequence_input.current_input).read():
            if initial_loop:
                (split_output_dir, split_output_filename, split_output) = check_output(sequence_input, 'Split')
                check_output_dir(split_output_dir)
                split_filename = '{0}/sequence_group_{1}'.format(bme_output_directory, ('0' * (total_files - len(str(group_counter)))) + str(group_counter))
                if '.' in sequence_input.current_input:
                    split_filename += '.{0}'.format(return_extension(sequence_input.current_input))
                bme_split_file = open(split_filename, 'w')
                initial_loop = False
            
            if sequence_counter == bme_split_number_in_groups:
                group_counter += 1
                bme_split_file.close()
                split_filename = '{0}/sequence_group_{1}'.format(bme_output_directory, ('0' * (total_files - len(str(group_counter)))) + str(group_counter))
                if '.' in sequence_input.current_input:
                    split_filename += '.{0}'.format(return_extension(sequence_input.current_input))
                bme_split_file = open(split_filename, 'w')
                bme_split_file.write(str(working_sequence) + '\n')
                sequence_counter = 0
            else:
                bme_split_file.write(str(working_sequence) + '\n')
            sequence_counter += 1
    bme_split_file.close()

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### Function: vespa_gene_selection:
### Details: Returns single files of sequences in the specified filepath if present within a csv file
def vespa_gene_selection (input_files):
    def return_sequence_filename(sequence_header):
        return_seq_filename = ''
        for header_characters in sequence_header:
            if return_seq_filename:
                if header_characters == '|':
                    return_seq_filename += '_'
                elif header_characters == '_':
                    return_seq_filename += header_characters
                elif header_characters == ' ':
                    pass
                elif not header_characters.isalnum():
                    break
            if header_characters.isalnum():
                return_seq_filename += header_characters
        return return_seq_filename
    
    import csv
    print 'VESPA: Gene selection'
    global bme_selection_csv, bme_output_directory
    csv_found, csv_missing, csv_list = ([], [], [row[0].strip() for row in csv.reader(open(bme_selection_csv, 'rU'))])
    directory_created = False
    
    for sequence_input in input_files:
        if not directory_created:
            (selected_output_dir, selected_output_filename, selected_output) = check_output(sequence_input, 'Selected')
            check_output_dir(selected_output_dir)
            directory_created = True
        for working_sequence in sequence_reader(sequence_input.current_input).read():
            for csv_entries in csv_list:
                if csv_entries in working_sequence.header:
                    csv_found.append(csv_entries)
                    selection_filename = '{0}/{1}'.format(bme_output_directory, csv_entries)
                    if '.' in sequence_input.current_input:
                        selection_filename += '.{0}'.format(return_extension(sequence_input.current_input))
                    bme_selection_file = open(selection_filename, 'w')
                    bme_selection_file.write(str(working_sequence) + '\n')
                    bme_selection_file.close()
    
    csv_missing = list(set(csv_list) - set(csv_found))
    
    if csv_missing:
        print '{0} genes not found, creating file: missing_genes.log'.format(len(csv_missing))
        bme_missing_file = create_unique_file('missing_genes.log')
        for missing_entries in csv_missing:
            bme_missing_file.write(missing_entries.strip() + '\n')
        bme_missing_file.close()

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### Function: vespa_check_SGO:
### Details: Checks for SGOs using seqeunce headers.
def vespa_check_SGO (input_files):        
    import os
    from collections import defaultdict
    print 'VESPA: Checking SGO status'
    check_SGO_log = create_unique_file('SGO_Check.log')
    for sequence_input in input_files:
        sgo_status = True
        species_counter = defaultdict(int)
        for working_sequence in sequence_reader(sequence_input.current_input).read():
            current_species = working_sequence.header[1:].strip().split('|')[0]
            if species_counter.has_key(current_species):
                species_counter[current_species] += 1
            else:
                species_counter[current_species] = 1
        for species_counts in species_counter.values():
            if species_counts != 1:
                sgo_status = False
        if sgo_status:
            check_SGO_log.write('{0},PASS\n'.format(sequence_input.current_input))
        else:
            check_SGO_log.write('{0},FAIL\n'.format(sequence_input.current_input))
    check_SGO_log.close()
    
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### Function: vespa_reduce_ensembl:
### Details: Reduces the length of Ensembl ID headers 
def vespa_reduce_ensembl (input_files):
    print 'VESPA: Reduce Ensembl Created Headers'
    reduced_conversion = create_unique_file('reduced_conversion.log')
    for sequence_input in input_files:
        (reduced_output_dir, reduced_output_filename, reduced_output) = check_output(sequence_input, 'Reduced')  
        bme_reduced_file = open(reduced_output, 'w')
        for working_sequence in sequence_reader(sequence_input.current_input).read():
            original_list, reduced_list = working_sequence.header.strip()[1:].split('|'), []
            for header_entries in original_list:
                check_header = header_entries.upper()
                if check_header.startswith('ENS'):
                    if check_header.split('0',1)[0].endswith('G'): 
                        reduced_list.append(header_entries)
                else:
                    reduced_list.append(header_entries)
            reduced_header = '>{0}\n'.format('|'.join(reduced_list))
            reduced_conversion.write('{0},{1}\n'.format(working_sequence.header.strip(),reduced_header.strip()))
            working_sequence.header = reduced_header
            bme_reduced_file.write(str(working_sequence) + '\n')
        bme_reduced_file.close()

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
###   _____ _       _           _  __      __        _       _     _           
###  / ____| |     | |         | | \ \    / /       (_)     | |   | |          
### | |  __| | ___ | |__   __ _| |  \ \  / /_ _ _ __ _  __ _| |__ | | ___  ___ 
### | | |_ | |/ _ \| '_ \ / _` | |   \ \/ / _` | '__| |/ _` | '_ \| |/ _ \/ __|
### | |__| | | (_) | |_) | (_| | |    \  / (_| | |  | | (_| | |_) | |  __/\__ \
###  \_____|_|\___/|_.__/ \__,_|_|     \/ \__,_|_|  |_|\__,_|_.__/|_|\___||___/
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### Global variables
bme_command_table = ['ensembl_clean', 'clean', 'translate', 'rev_complement', 'create_database', 'individual_sequences', 'gene_selection',
                    'sgo_check', 'best_reciprocal_groups', 'reciprocal_groups', 'similarity_groups', 'create_subtrees', 'map_alignments',
                    'infer_genetree', 'codeml_setup', 'split_sequences', 'codeml_reader', 'setup_reciprocal_input', 'prottest_setup',
                    'prottest_reader', 'metal_compare', 'mrbayes_setup', 'create_branch', 'mrbayes_reader', 'link_input', 'reduce_ensembl',
                    'h', 'help']
bme_sequence_database_location = ''
bme_assign_database_filename = ''
bme_output_directory = ''
bme_output_filename = ''

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### General sequence tool variables
bme_remove_internal_stop = True
bme_remove_terminal_stop = True
bme_label_with_filename = False
bme_infer_labels_ensembl = False
bme_split_number_in_groups = 100
bme_selection_csv = ''
bme_format_blast_database = False

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### General similarity variables
bme_similarity_bme_format_location = ''
bme_similarity_data_format = 'blast'
bme_similarity_e_value_cutoff = ''
bme_similarity_percent_identity_cutoff = ''
bme_similarity_alignment_length_cutoff = ''

blast_alignment_length_warn = False
hmmer_percent_identity_warn = False

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### General 3rd phase variables
bme_metAl_compare_files = []
bme_metAl_compare_dir = False
bme_metAl_cutoff = 0.05
bme_supported_model_list = ''
bme_mrbayes_mcmc_gen = 200000
bme_mrbayes_mcmc_chains = 4
bme_mrbayes_mcmc_temp = 0.2
bme_mrbayes_mcmc_burnin = 0.25


### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### General tree variables
bme_subtree_log_file = 'vespa_subtrees.log'

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### General codeML variables
bme_species_tree = ''
bme_branch_label_table = ''       
bme_in_paralogs = False
bme_alignment_path = ''
bme_main_output = ''



### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
###   _____                                          _ _      _            
###  / ____|                                        | | |    (_)           
### | |     ___  _ __ ___  _ __ ___   __ _ _ __   __| | |     _ _ __   ___ 
### | |    / _ \| '_ ` _ \| '_ ` _ \ / _` | '_ \ / _` | |    | | '_ \ / _ \
### | |___| (_) | | | | | | | | | | | (_| | | | | (_| | |____| | | | |  __/
###  \_____\___/|_| |_| |_|_| |_| |_|\__,_|_| |_|\__,_|______|_|_| |_|\___|   
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### Function: command_line:
### Details: Handles user commands and assigning variables
def command_line():
    global bme_sequence_database_location, bme_assign_database_filename, bme_output_directory, bme_output_filename, bme_command_table
    global bme_remove_internal_stop, bme_remove_terminal_stop, bme_label_with_filename, bme_infer_labels_ensembl, bme_split_number_in_groups, bme_selection_csv, bme_format_blast_database
    global bme_similarity_e_value_cutoff, bme_similarity_percent_identity_cutoff, bme_similarity_alignment_length_cutoff
    global bme_metAl_compare_files, bme_metAl_compare_dir, bme_metAl_cutoff, bme_supported_model_list, bme_mrbayes_mcmc_gen, bme_mrbayes_mcmc_chains, bme_mrbayes_mcmc_temp, bme_mrbayes_mcmc_burnin
    global bme_species_tree, bme_branch_label_table, bme_in_paralogs, bme_alignment_path
    

    def command_splitter (command_array):
        command_return, option_return = '', []
        for command_groups in command_array:
            if command_return and not command_groups.startswith('-'):
                yield command_return, option_return
                command_return = command_groups
                option_return = []
            if not command_return and not command_groups.startswith('-'):
                command_return = command_groups
            if command_groups.startswith('-'):
                option_return.append(command_groups)
        yield command_return, option_return
    
    def assign_input_files (input_directory, input_varible):
        import os, sys
        input_to_return = []
        if input_directory:
            for path, sub_dirs, file_list in os.walk(input_varible):
                for files in file_list:
                    if not files.startswith('.'):
                        input_to_return.append(command_line_data(os.path.join(path, files), True))
        else:
            input_to_return.append(command_line_data(input_varible, False))
        return input_to_return
    
    def assign_compare_files(input_varible):
        import os, sys
        if check_if_input_directory(input_varible):
            for path, sub_dirs, file_list in os.walk(input_varible):
                for files in file_list:
                    if not files.startswith('.'):
                        bme_metAl_compare_dir = True
                        bme_metAl_compare_files.append(os.path.join(path, files))
        else:
            bme_metAl_compare_files.append(input_varible)

        
    import sys, os
    command_input = []
    input_directory_check = False
    
    if 'h' in sys.argv[1:]:
        if len(sys.argv[1:]) > 1:
            for help_request in sys.argv[2:]:
                help_message(help_request) 
        else:
            help_message('')
        sys.exit()
    elif 'help' in sys.argv[1:]:
        if len(sys.argv[1:]) > 1:
            for help_request in sys.argv[2:]:
                help_message(help_request) 
        else:
            help_message('')
        sys.exit()
    
    for current_command, options_list in command_splitter(sys.argv[1:]):
        for options in options_list:
            if options.startswith('-input='):
                input_directory_check = check_if_input_directory(options.split('=')[1])
                command_input = assign_input_files(input_directory_check, options.split('=')[1])
        for options in options_list:
            if options.startswith('-output='):
                bme_output_directory = options.split('=')[1]
                bme_output_filename = options.split('=')[1]
                if input_directory_check:
                    for input_in_dir in command_input:
                        input_in_dir.current_output_filename = input_in_dir.current_input.split('/')[-1]
                        input_in_dir.current_output_dir = options.split('=')[1]
                        input_in_dir.current_output = '{0}/{1}'.format(options.split('=')[1], input_in_dir.current_input.split('/')[-1])
                else:
                    for input_in_dir in command_input:
                        input_in_dir.current_output = options.split('=')[1]
            elif options.startswith('-label_filename='):
                if 'true' in options.split('=')[1].lower():
                    bme_label_with_filename = True
                else:
                    bme_label_with_filename = False
            elif options.startswith('-infer_ensembl_species='):
                if 'true' in options.split('=')[1].lower():
                    bme_infer_labels_ensembl = True
                else:
                    bme_infer_labels_ensembl = False
            elif options.startswith('-rm_internal_stop='):
                if 'true' in options.split('=')[1].lower():
                    bme_remove_internal_stop = True
                else:
                    bme_remove_internal_stopp = False
            elif options.startswith('-cleave_terminal='):
                if 'true' in options.split('=')[1].lower():
                    bme_remove_terminal_stop = True
                else:
                    bme_remove_terminal_stop = False
            elif options.startswith('-format_blast='):
                if 'true' in options.split('=')[1].lower():
                    bme_format_blast_database = True
                else:
                    bme_format_blast_database = False        
            elif options.startswith('-selection_csv='):
                bme_selection_csv = options.split('=')[1]
            elif options.startswith('-output_database='):
                bme_assign_database_filename = options.split('=')[1]
            elif options.startswith('-split_number='):
                bme_split_number_in_groups = int(options.split('=')[1])
            elif options.startswith('-subtree_log='):
                bme_subtree_log_file = options.split('=')[1]
            elif options.startswith('-species_tree='):
                bme_species_tree = options.split('=')[1]
            elif options.startswith('-branch_file='):
                bme_branch_label_table = options.split('=')[1]
            elif options.startswith('-allow_inparalogs='):
                if 'true' in options.split('=')[1].lower():
                    bme_in_paralogs = True
                else:
                    bme_in_paralogs = False
            elif options.startswith('-format='):
                bme_similarity_data_format = options.split('=')[1]
            elif options.startswith('-e_value='):
                bme_similarity_e_value_cutoff = float(options.split('=')[1])
            elif options.startswith('-percent_identity='):
                bme_similarity_percent_identity_cutoff = float(options.split('=')[1])
            elif options.startswith('-alignment_length='):
                bme_similarity_alignment_length_cutoff = float(options.split('=')[1])
            elif options.startswith('-database='):
                bme_sequence_database_location = options.split('=')[1]
            elif options.startswith('-alignment_path='):
                bme_alignment_path = options.split('=')[1]
            elif options.startswith('-compare='):
                assign_compare_files(options.split('=')[1])
            elif options.startswith('-metal_cutoff='):
                bme_metAl_cutoff = float(options.split('=')[1])
            elif options.startswith('-model_list='):
                bme_supported_model_list = options.split('=')[1]
            elif options.startswith('-mcmc_gen='):
                bme_mrbayes_mcmc_gen = int(options.split('=')[1])
            elif options.startswith('-mcmc_chains='):
                bme_mrbayes_mcmc_chains = int(options.split('=')[1])
            elif options.startswith('-mcmc_temp='):
                bme_mrbayes_mcmc_temp = float(options.split('=')[1])
            elif options.startswith('-mcmc_burnin='):
                bme_mrbayes_mcmc_burnin = float(options.split('=')[1])


        if current_command.lower() in bme_command_table:
            if command_input:
                #1st Phase
                if 'ensembl_clean' in current_command.lower():
                    vespa_clean_ensembl(command_input)
                elif 'clean' in current_command.lower():
                    vespa_clean(command_input)
                elif 'translate' in current_command.lower():
                    vespa_translate(command_input)
                elif 'rev_complement' in current_command.lower():
                    vespa_reverse_complement(command_input)
                elif 'create_database' in current_command.lower():
                    vespa_create_database(command_input,)
                elif 'individual_sequences' in current_command.lower():
                    vespa_individual_sequences(command_input)   
                elif 'split_sequences' in current_command.lower():
                    vespa_split_in_groups(command_input)
                elif 'gene_selection' in current_command.lower():
                    vespa_gene_selection(command_input)
                elif 'sgo_check' in current_command.lower():
                    vespa_check_SGO(command_input)   
                elif 'reduce_ensembl' in current_command.lower():
                    vespa_reduce_ensembl(command_input)    
                    
                #2nd Phase
                elif 'setup_reciprocal_input' in current_command.lower():
                    vespa_setup_reciprocal_input(command_input)
                elif 'best_reciprocal_groups' in current_command.lower():
                    vespa_best_reciprocal_similarity_groups(command_input)
                elif 'reciprocal_groups' in current_command.lower():
                    vespa_similarity_groups(command_input,True)
                elif 'similarity_groups' in current_command.lower():
                    vespa_similarity_groups(command_input,False)
                
                #3rd Phase
                elif 'metal_compare' in current_command.lower():
                    vespa_metAl_compare(command_input)
                elif 'prottest_setup' in current_command.lower():
                    vespa_setup_prottest(command_input)
                elif 'prottest_reader' in current_command.lower():
                    vespa_prottest_reader(command_input)
                elif 'mrbayes_setup' in current_command.lower():
                    vespa_setup_mrbayes(command_input)
                
                
                #4th Phase
                elif 'mrbayes_reader' in current_command.lower():
                    vespa_mrbayes_reader(command_input)
                elif 'create_branch' in current_command.lower():
                    vespa_branch_table(command_input)
                elif 'create_subtrees' in current_command.lower():
                    vespa_subtrees(command_input)
                elif 'map_alignments' in current_command.lower():
                    vespa_map_protein_gaps(command_input)
                elif 'infer_genetree' in current_command.lower():
                    vespa_infer_genetree(command_input)
                elif 'link_input' in current_command.lower():
                    vespa_link_input(command_input)    
                elif 'codeml_setup' in current_command.lower():
                    vespa_codeml_setup(command_input)
                    
                #5th Phase
                elif 'codeml_reader' in current_command.lower():
                    vespa_codeml_reader(command_input)
            else:
                print 'No input specified for command: {0}. Please check command-line input'.format(current_command)
        else:
            print 'Specified command ({0}) not found. Please check command-line input'.format(current_command) 
            

import sys, os
if len(sys.argv) > 1:
    command_line()
else:
    help_message('')