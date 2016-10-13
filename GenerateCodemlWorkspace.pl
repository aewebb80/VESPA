#!/usr/bin/perl -w
use strict;
use CodemlWrapper::Alignment;
use CodemlWrapper::Tree;
use CodemlWrapper::Job;
#GenerateCodemlWorkspace v1.01 by Tom Walsh (2012/07/10)

my @default_models = ('m0', 'm1Neutral', 'm2Selection', 'm3Discrtk2', 'm3Discrtk3', 'm7', 'm8', 'm8a', 'modelA', 'modelAnull');
my @default_omegas = (0, 1, 2, 10);

=head1 NAME

GenerateCodemlWorkspace.pl - Sets up a Codeml workspace using an alignment and one or more tree files. 

=head1 USAGE

perl GenerateCodemlWorkspace.pl alignment_file tree_files... codeml_workspace [-models=model_set] [-omegas=omega_set]

=head1 DESCRIPTION

GenerateCodemlWorkspace.pl takes as input an alignment file and one or more tree 
files. Using this input, it generates a Codeml workspace, an ordered directory 
structure in which each Codeml task is placed in its own directory, with its own 
Codeml control, alignment and tree files. 

=head1 ARGUMENTS

=over 4

=item • alignment_file

Nucleotide alignment file in PHYLIP sequential, PHYLIP interleaved or FASTA format.

=item • tree_files...

One or more tree files in 'PAML', 'PHYLIP' or Newick tree format. These Codeml 
trees must have identical topologies but be uniquely labelled (or unlabelled). 

=item • codeml_workspace

Name of the Codeml 'workspace', a directory structure containing a set of 
Codeml tasks associated with a specific Codeml job. 

=back

=head1 OPTIONS

=over 4

=item • -models=model_set

A model to include in the Codeml workspace, specified individually on the command 
line or as a list in a text file with extension '.txt'. If no models are specified,
the default models are used.

=item • -omegas=omega_set

An initial omega to include in the Codeml workspace, specified individually on 
the command line or as a list in a text file with extension '.txt'. If no omegas
are specified, the default omegas are used.

=back

=head1 DEPENDENCIES

=over 4

=item • L<CodemlWrapper::Job>

=item • L<CodemlWrapper::Alignment>

=item • L<CodemlWrapper::Tree>

=back

=cut 

#Accept arguments.
my (@arguments) = @ARGV;

#Init variables.
undef my @standard_arguments;
undef my @models;
undef my @omegas;

#Process each argument.
foreach my $arg (@arguments)
{
	#If marked as a model or model list file, add to list of models..
	if($arg=~/^-models=(.+)/i)
	{
		my $m = $1;
		
		my @list = ($m);
		if($m=~/\.txt$/)
		{
			open(LIST_FILE, '<', $m) or die($0." failed: can't open \'".$m."\' for reading, exiting;");
			@list = <LIST_FILE>;
			chomp(@list);
			close(LIST_FILE) or die($0." failed: can't close \'".$m."\' after reading, exiting;");
		}
		
		foreach my $item (@list)
		{
			#Trim leading and trailing whitespace.
			$item=~s/(^\s*|\s*$)//g;
		
			if( !(grep {$item eq $_} @models) )
			{push(@models, $item);}
		}
	}
	#..otherwise, if marked as an omega or omega list file, add to list of initial omegas..
	elsif($arg=~/^-omegas=(.+)/i)
	{
		my $o = $1;
		
		my @list = ($o);
		if($o=~/\.txt$/)
		{
			open(LIST_FILE, '<', $o) or die($0." failed: can't open \'".$o."\' for reading, exiting;");
			@list = <LIST_FILE>;
			chomp(@list);
			close(LIST_FILE) or die($0." failed: can't close \'".$o."\' after reading, exiting;");
		}
		
		foreach my $item (@list)
		{
			#Trim leading and trailing whitespace.
			$item=~s/(^\s*|\s*$)//g;		
		
			if( !(grep {$item eq $_} @omegas) )
			{push(@omegas, $item);}
		}
	}
	#..otherwise, it's a standard argument.
	else
	{
		push(@standard_arguments, $arg);
	}
}

#If no models specified, set default.
@models = @default_models if(@models==0);

#If no initial omegas specified, set default.
@omegas = @default_omegas if(@omegas==0);

#Verify that there are at least three arguments.
if(@standard_arguments < 3)
{die($0." failed: must have at least 3 arguments, exiting;");}

#Set variables from arguments.
my $codeml_workspace = pop(@standard_arguments);
my ($alignment_file, @tree_files) = @standard_arguments;

#If existing Codeml workspace directory can't be removed, die.
if(-d $codeml_workspace)
{system("rm -rf ".$codeml_workspace)==0 or die($0." failed: can't remove existing directory \'".$codeml_workspace."\', exiting;");}

#Load given alignment from file.	
my $alignment = Alignment->new($alignment_file);

#Load each specified tree from file, and take the name of the file (without its
#extension) as the name of the tree to be added to the Codeml job.
undef my %trees;
foreach my $tree_file (@tree_files)
{
	my $tree_name = substr($tree_file, rindex($tree_file, "/") + 1 );	
	$tree_name=~s/\.[^\.]+$//;
		
	$trees{$tree_name} = Tree->new($tree_file);
}

#Create new CodemlWrapper::Job. 	
my $CodemlJob = Job->new();

#Set models for the Codeml job.
$CodemlJob->set_models(@models);

#Set initial ω-values for the Codeml job.
$CodemlJob->set_initial_omegas(@omegas);

#Add the given alignment to the Codeml job.
$CodemlJob->add_alignment($alignment);

#Add the given trees to the Codeml job.
foreach my $tree_name (keys(%trees))
{
	$CodemlJob->add_tree($tree_name, $trees{$tree_name});	
}

#Save Codeml job to Codeml workspace.
$CodemlJob->save($codeml_workspace);
