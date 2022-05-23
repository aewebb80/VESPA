#!/usr/bin/env perl
use strict;
use warnings;
use Carp qw(carp cluck croak confess);
use CodemlWrapper::Alignment;
use CodemlWrapper::Control;
use CodemlWrapper::Output;
use CodemlWrapper::Tree;
#CodemlWrapper::Job v1.02 by Tom Walsh (2012/07/12)

=head1 NAME

CodemlWrapper::Job - Wrapper for a Codeml job.

=head1 SYNOPSIS: PREPARATION

	use CodemlWrapper::Tree;
	use CodemlWrapper::Job;
		
	#SET MAIN VARIABLES.
	my $codeml_workspace = 'Workspace1';
	my $alignment_file = 'align.phy';
	my $tree_file = 'tree';
	
	my @models = ('m0', 'm1Neutral', 'm2Selection', 'm3Discrtk2', 'm3Discrtk3', 'm7', 'm8', 'm8a', 'modelA', 'modelAnull');	
	my @initial_omegas = (0, 1, 2, 10);				
	my @foreground_sequences = ('Human', 'Chimp', 'Gorilla');				
	
	#SET UP PRIMATE TREE.
	#Read in tree from file.
	my $tree = Tree->new($tree_file);					
	
	#Label foreground branches using the foreground sequence names.
	$tree->set_foreground_genes(@foreground_sequences);
		
	#SET UP CODEML WORKSPACE.
	#Create new CodemlWrapper::Job.
	my $CodemlJob = Job->new();
	
	#Set models.
	$CodemlJob->set_models(@models);
	
	#Set initial ω-values.
	$CodemlJob->set_initial_omegas(@initial_omegas);
	
	#Add alignment from file.
	$CodemlJob->add_alignment($alignment_file);
	
	#Add ape tree.
	$CodemlJob->add_tree('Apes', $tree);
		
	#Save workspace.
	$CodemlJob->save($codeml_workspace);

=head1 SYNOPSIS: POST-PROCESSING

	use CodemlWrapper::Job;
	
	#SET MAIN VARIABLES.
	my $codeml_workspace = 'Workspace1';			
		
	#CREATE CODEML JOB SUMMARY REPORTS.
	#Create new CodemlWrapper::Job.
	my $CodemlJob = Job->new();
	
	#Load Codeml job from Codeml workspace.
	$CodemlJob->load($codeml_workspace);
	
	#Create a summary report for the Codeml job.
	$CodemlJob->create_summary_report($codeml_workspace."/Codeml_Summary_Report.txt");

=head1 DESCRIPTION

Codeml is a component of Ziheng Yang's PAML package. Codeml performs selection 
analysis on gene sequences, allowing episodes of positive selection of a gene on 
evolutionary timescales to be inferred. For information on PAML see the PAML manual.

CodemlWrapper::Job is a Perl wrapper module for Codeml (PAML4.4e). It helps with the 
preparation and post-processing of a Codeml job, through manipulation of the Codeml 
'workspace': a directory structure containing the set of Codeml tasks in the Codeml 
job. CodemlWrapper::Job usually runs silently. To see full output, include 'Verbose' 
as an argument when creating a L</new> object.

Each Codeml task is located in its own directory within the workspace, which 
contains at least a Codeml control file ('codeml.ctl'), a nucleotide alignment 
and a Codeml tree. All tasks in a Codeml job take the same alignment and tree 
topology as input, although a Codeml tree may labelled differently (or unlabelled) 
depending on the lineage (or lineages) of interest. For a given labelled (or unlabelled)
tree, Codeml tasks can also differ in the evolutionary model and initial ω-value used. 
Note that although PAML I<does> allow for multiple trees and models in a Codeml 
task, for simplicity's sake CodemlWrapper::Job I<does not>.

Each Codeml job has four main attributes: an alignment, one or more Codeml trees, 
a set of evolutionary models and a set of initial ω-values. The alignment and Codeml 
tree (or trees) are input data, while the models and initial ω-values determine 
the set of tasks for each tree in the Codeml job. All the tasks in the Codeml 
job use the same alignment. Any number of trees may be used as input, up to an 
arbitrary limit. All but one of these must be labelled with a unique foreground, 
such that no two trees can have the same set of foreground genes. One tree is always 
kept as the topology tree, an unlabelled Codeml tree (named 'Sites' by default). 
The topology tree is used as input for homogeneous and site-specific models and 
is added automatically to every Codeml job. Each labelled Codeml tree is used 
as input for branch-specific and branch-site  models with respect to its foreground 
lineage(s). 

The 'active sets' of models and ω-values are used to determine which tasks the 
Codeml job will handle, whether at the preparation or post-processing stage. Whether 
a model or ω-value is in the active set determines whether it's 'visible' to the 
Codeml job. For example, adding the 'm1Neutral' model ensures that a saved Codeml 
workspace will contain 'm1Neutral' Codeml tasks, while removing '0' from the set 
of ω-values will mean tasks with an initial ω-value of zero are not loaded with 
a Codeml workspace, even if they are present. Not all models in the active set 
will be used for all trees: homogeneous and site-specific models can only be used 
with the topology tree, while branch-specific and branch-site models are used with 
labelled Codeml trees. Similarly, not all models will use the full set of ω-values: 
if a model requires that an ω-value be fixed (e.g. 'modelAnull'), only the fixed 
ω-value will be used for that model, and then only if that ω-value is in the active 
set. In general, the default ω-values should be used. 

Preparation of a Codeml job involves setup of data, models and ω-values and 
culminates in a call to L</save>, which saves the workspace directory to disk. After 
the Codeml job has been run, post-processing begins with a call to L</load>, which 
loads the contents of the Codeml workspace, including any output. To prevent the 
output of a Codeml job being associated with incorrect input, the data of a loaded 
workspace can't be directly modified: the alignment and trees are taken as read-only 
and can't be changed. Models and parameters may be changed in this context, but 
changing these doesn't change the actual contents of the workspace, just whether 
they are 'visible' for post-processing. Once a Codeml workspace is loaded, the 
results can be processed and reports can be created, including for example a summary 
report for the job as a whole, or a positive site report for a particular model. 

CodemlWrapper::Job creates workspace directories in the following order: tree, model, 
initial ω-value. For example, in a workspace called 'Codeml_Workspace', with a 
tree called 'Eutheria' (in which Eutherian taxa are labelled as foreground), 
a selection analysis with 'modelA' and initial ω-value of '0' would be located in 
'Codeml_Workspace/Eutheria/modelA/Omega0'. CodemlWrapper::Job can load a workspace even 
if it doesn't follow this convention; the only requirement is that all tasks for 
each tree share a unique common directory within the Codeml workspace. 

=head2 Input

=head3 Alignment

The nucleotide alignment is the same for all tasks in a Codeml job. Gaps, if they 
are present, should be within the reading frame (so that a given codon will have 
either three gaps or none). In the CodemlWrapper::Job module, the alignment can be added 
using L</add_alignment>, and removed using L</remove_alignment>. Most of the work 
done with the alignment is done using L<'CodemlWrapper::Alignment'|/CodemlWrapper::Alignment> 
objects. These can take as input in PHYLIP sequential, PHYLIP interleaved or FASTA 
format. The sequence names in the alignment should match the gene names in all 
trees. For more information see the L<'CodemlWrapper::Alignment'|/CodemlWrapper::Alignment> 
documentation.

=head3 Trees

The tree topology is identical for all tasks in a Codeml job, although a given 
job may contain many uniquely labelled Codeml trees. In the CodemlWrapper::Job module, 
a tree can be added using L</add_tree>, and removed using L</remove_tree>. Most 
of the work done with trees is done using L<'CodemlWrapper::Tree'|/CodemlWrapper::Tree> objects.
These can take as input files in 'PAML', 'PHYLIP' or standard Newick tree format. 
Trees should be unrooted, but need not be as CodemlWrapper::Job automatically unroots 
trees when they're loaded. The gene names in each tree should match all other trees 
and also the alignment. For more info see the L<'CodemlWrapper::Tree'|/CodemlWrapper::Tree> docs.

=head2 Settings

Two main settings can be modified through the CodemlWrapper::Job wrapper: Codeml's 
evolutionary models and the initial ω-values used. 

=head3 Models

A Codeml job can have between one and twelve active models. The active set of models 
can be set by calling L</set_models>, and can be obtained by calling L</get_models>. 
All twelve models are listed below. All models except '2ratios' and 'modelB' are 
active by default.

=over 4

=item • m0

=item • 2ratios

=item • m1Neutral

=item • m2Selection

=item • m3Discrtk2

=item • m3Discrtk3

=item • m7

=item • m8

=item • m8a

=item • modelA

=item • modelAnull

=item • modelB

=back

These twelve models are grouped into four model categories: Homogeneous ('m0'), 
Branch-specific ('2ratios'), Site-specific ('m1Neutral', 'm2Selection', 'm3Discrtk2', 
'm3Discrtk3', 'm7', 'm8', 'm8a') and Branch-site models ('modelA', 'modelAnull', 'modelB').

=head3 Initial ω-Values

A Codeml job can have 1 or more initial ω-values, such that each ω-value is a 
positive real number. The active set of ω-values can be modified by calling 
L</set_initial_omegas>, and can be obtained by calling L</get_initial_omegas>. 
The default ω-values are 0, 1, 2 and 10. For most purposes, the default ω-values 
are sufficient.

=head2 Output

CodemlWrapper::Job output comes in the form of tab-delimited reports. There are currently 
five types of 'report': a summary report, a positive site report, a performance report, 
an LRT report and a positive site MSA. 

=head3 Summary Report

The summary report is obtained by calling L</create_summary_report>, and shows the 
results of a Codeml job, optionally for a specific tree. By default, if no tree 
is specified, a summary report is produced for all trees in the Codeml job.

Every loaded task for the trees of interest is included in the report, except when 
two or more tasks were run for the same model (e.g. with different initial ω-values). 
In this case, only the task with the highest log-likelihood is retained. A likelihood 
ratio test (LRT) is performed on certain models (the alternative models in the list 
of LRTs below). For every LRT performed, the null model I<must> be present. The 
following LRTs are supported:

=over 4

=item • m0 vs m3Discrtk2

=item • m1Neutral vs 2ratios

=item • m3Discrtk2 vs m3Discrtk3

=item • m1Neutral vs m2Selection

=item • m7 vs m8

=item • m8a vs m8

=item • m1Neutral vs modelA

=item • modelAnull vs modelA

=item • m3Discrtk2 vs modelB

=back

If the null model is rejected for every LRT of a model and the model's parameter 
estimates indicate that at least one class of sites has ω > 1, then positive selection 
is inferred. 

The headings of a Codeml summary report are: 

=over 4

=item • Model

The name of the model.

=item • Tree

The tree on which the model was used in Codeml (e.g. 'Sites').

=item • Model Type

The category of model (e.g. 'Site-specific', 'Branch-site').

=item • p

The number of parameters in the model.

=item • w (t=0)

The initial ω-value used. This indicates which output file was selected for the 
report. For example, if 'modelA' had several tasks with different initial ω-values, 
and the task with an initial ω-value of '10' had the highest log-likelihood, that 
task's output would be used, and the value under this heading would be '10'. 

=item • lnL

The log-likelihood of the model. 

=item • LRT Result

The result of the likelihood ratio test (LRT). If this model is an alternative model 
for an LRT (or LRTs) and the null model is rejected, this model's name is shown. 
If an LRT fails to reject a null model, the null model is shown. If no LRT is performed 
on the model, the value under this heading is 'N/A'. 

=item • Parameter Estimates

Estimates of the model parameters. 

=item • Positive Selection

Indicates whether positive selection has occurred. For the models 'm1Neutral', 
'm7', 'm8a' and 'modelAnull', the value is 'Not Allowed'. For 'm0', '2ratios', 
'm3Discrtk2', 'm3Discrtk3', 'm2Selection', 'm8', 'modelA' and 'modelB', the value 
is 'Yes' or 'No', depending on whether or not there is evidence of positive selection. 
If the null model is rejected for every LRT of a model and the model's parameter 
estimates indicate that at least one class of sites has ω > 1, then positive selection 
is inferred and the value under this heading is 'Yes'.

=item • Positively Selected Site Summary

The column under this heading contains positive site summaries, where appropriate.
No positive site summary is produced if the model doesn't allow positive sites, 
or if positive selection hasn't been inferred.

In the report, the specific heading of this column depends on the dimensions of the 
alignment and the value of the 'positive site probability threshold'. In order to be 
reported, a positive site's probability must be greater than the positive site probability
threshold. This value can be set by calling L</set_positive_site_probability_threshold>, 
and can be obtained by calling L</get_positive_site_probability_threshold>, (default=0.5).

If positive selection occurs for a model and positive sites are found, a positive 
site summary is produced. Two positive site inference types are used by Codeml: 
Bayes Empirical Bayes (BEB) and Naive Empirical Bayes (NEB). CodemlWrapper::Job gives 
preference to BEB positive sites over NEB sites. When both are present, CodemlWrapper::Job 
ignores the NEB sites, so NEB sites are only reported in the absence of BEB sites.
All positive site summaries indicate the number and inference type (i.e. NEB or BEB) 
of positive sites for the given model. The site position information given in the 
positive site summary differs depending on the sequence(s) specified in the call 
to L</create_summary_report>. If no sequences are specified, positive site info 
is given with respect to the Codeml alignment. If one or more sequences are 
specified, positive site info is given with respect to the specified sequences. 
(Pos site info can also be output in terms of an original alignment - see L</create_summary_report>).

Note the distinction made between a positive site summary of an I<alignment> and 
that of a gene I<sequence>. In the case of the alignment, only site positions are 
reported. No amino acids are included, since these could differ between sequences 
at that site in the alignment. Where sequences are specified, positive site positions 
are reported with the amino acid at that site in that sequence. Sites for which a 
sequence has a gap character are excluded from the positive site summary, since 
they're not part of that sequence in any meaningful way. Instead, a tally of gap 
sites is included where appropriate. 

For clarity, not all positive site information is included in the positive site 
summary. For some models (e.g. m8), positive sites are accompanied by other information. 
This isn't included in the summary report, but can be obtained from a L</Positive Site Report>.

=back

=head3 Positive Site Report

The positive site report is obtained by calling L</create_pos_site_report>, 
and shows detailed positive site information for a given model in a given tree, 
with the option to select a sequence of interest. If no sequence of interest is 
given, the positive site report will relate to the alignment as a whole. Pos site 
info can also be output in terms of an original alignment - see L</create_pos_site_report>).
The positive site report is useful for obtaining full positive site information in cases where 
positive selection has been inferred. If positive sites aren't possible within the 
given model, the process fails. If positive sites are possible but haven't been
inferred, the process continues but no output file is produced.

Note that gap sites for the sequence of interest are completely excluded from the 
positive site report; a tally of gaps is given in the summary report. As with the 
Codeml summary report, positive sites of type BEB are given preference over NEB 
sites. When both are present, the NEB sites are ignored, and NEB sites are only 
reported in the absence of BEB sites.

The headings of a positive site report are: 

=over 4

=item • Site Info

The position of the site in the alignment/sequence (ranging from 1 to N, where N 
is the number of codon positions in the alignment/sequence). In the report, the 
specific heading of this column depends on the dimensions of the alignment and 
the positive site inference type.

=item • Amino Acid

If the positive site report is for a sequence of interest (as opposed to the alignment), 
the amino acid is given for each site. 

=item • P(w>1)

Probability that ω-value of the site is greater than 1. In the report, the heading 
of this column will include the 'positive site probability threshold'.

=item • Posterior Mean (w)

Posterior mean value of ω. Under some models the standard error (SE) is also given. 

=back

The positive site report is intended to complement the summary report, so that 
positive site reports can be created for each model with positively selected sites 
in the summary report. 

=head3 Positive Site MSA

The positive site MSA is obtained by calling L</create_pos_site_msa>, and is an
alignment showing positive sites highlighted. Site-specific pos sites are shown
in red, while branch-site specific pos sites are shown in red for foreground
sequences and blue for background sequences. Pos site info can also be shown on 
an original alignment - see L</create_pos_site_msa>).

If positive sites aren't possible within the given model, the process fails. 
If positive sites are possible but haven't been inferred, the process continues 
but no output file is produced.

Note that gap sites for the sequence of interest are completely excluded from the 
positive site MSA; a tally of gaps is given in the summary report. As with the 
Codeml summary report, positive sites of type BEB are given preference over NEB 
sites. When both are present, the NEB sites are ignored, and NEB sites are only 
reported in the absence of BEB sites.

The positive site MSA is intended to complement the pos site report and provide 
a quick way to visualise positive sites along the alignment.

=head3 Performance Report

The performance report is obtained by calling L</create_performance_report>, and 
shows performance information for all loaded models. This information is taken 
directly from the beginning and end of the Codeml output file.

The headings of a performance report are: 

=over 4

=item • Program Info

PAML program info, typically including version number and date.

=item • Model

PAML model used. 

=item • Omega

Initial ω-value used. 

=item • Sequence Count

Number of sequences in the input nucleotide alignment.

=item • Sequence Length

Number of codons in each sequence of the input nucleotide alignment.

=item • Time Used

Time (in seconds) used by Codeml to complete the task.

=back

=head1 OBJECT METHODS

=cut

package Job;
{

use constant TRUE  => 1;
use constant FALSE => 0;
use constant DEF   => 1;
use constant MODULE_NAME_LENGTH => 3;
use constant ARBITRARY_LIMIT => 1000;
use constant DEFAULT_PSP_THRESHOLD => 0.5;

my %model_categories = 
( 
	'Homogeneous' 		=> [ 'm0' ],
	'Branch-specific'	=> [ '2ratios' ],
	'Site-specific'   	=> [ 'm1Neutral', 'm2Selection', 'm3Discrtk2', 'm3Discrtk3', 'm7', 'm8', 'm8a' ],
	'Branch-site' 		=> [ 'modelA', 'modelAnull', 'modelB' ]		
);

my @supported_model_categories = ('Homogeneous', 'Branch-specific', 'Site-specific', 'Branch-site');

my @supported_models = 
(
	'm0',
	'2ratios',
	'm1Neutral', 
	'm2Selection', 
	'm3Discrtk2', 
	'm3Discrtk3', 
	'm7', 
	'm8', 
	'm8a',
	'modelA', 
	'modelAnull', 
	'modelB'
);

my @default_models = ('m0', 'm1Neutral', 'm2Selection', 'm3Discrtk2', 'm3Discrtk3', 'm7', 'm8', 'm8a', 'modelA', 'modelAnull');
my @default_omegas = (0, 1, 2, 10);

my %report_headings = 
(
	'performance-report'	=> ['Program Info', 'Model', 'omega', 'Sequence Count', 'Sequence Length', 'Time Used (secs)'],	
	'summary-report'		=> ['Model', 'Tree', 'Model Type', 'p', 'w (t=0)', 'lnL', 'LRT Result', 'Parameter Estimates', 'Positive Selection'], 
	'lrt-report'			=> ['Tree', 'LRT', 'Null Model lnL', 'Alternative Model lnL', 'D', 'Critical Value', 'Null Rejected?']
);

my @supported_model_features = ('category', 'positive-site-capability', 'free-parameters', 'null-models', 'fixed-omega');

my %model_feature = 
(
	'm0' 			=> 	{ 
							'category'					=>	'Homogeneous', 
							'positive-site-capability'	=>	'N/A', 
							'free-parameters'			=>	1							
						},	
						
	'2ratios'		=> 	{ 
							'category'					=>	'Branch-specific',
							'positive-site-capability'	=>	'N/A', 
							'free-parameters'			=>	2, 
							'null-models'				=> 	['m1Neutral']												
						},
						
	'm1Neutral' 	=> 	{ 
							'category'					=>	'Site-specific', 
							'positive-site-capability'	=>	'Not Allowed', 
							'free-parameters'			=>	1							
						},
	
	'm2Selection' 	=> 	{ 
							'category'					=>	'Site-specific', 
							'positive-site-capability'	=>	'Allowed', 
							'free-parameters'			=>	2, 
							'null-models'				=> 	['m1Neutral']
						},
	
	'm3Discrtk2' 	=> 	{	
							'category'					=>	'Site-specific',
							'positive-site-capability'	=>	'Allowed', 
							'free-parameters'			=>	3, 
							'null-models'				=>	['m0']
						},
	
	'm3Discrtk3' 	=> 	{ 	
							'category'					=>	'Site-specific', 
							'positive-site-capability'	=>	'Allowed', 
							'free-parameters'			=>	5, 
							'null-models'				=> 	['m3Discrtk2']	
						},
	
	'm7' 			=> 	{ 
							'category'					=>	'Site-specific', 
							'positive-site-capability'	=>	'Not Allowed', 
							'free-parameters'			=>	2
						},
	
	'm8'			=> 	{
							'category'					=>	'Site-specific', 
							'positive-site-capability'	=>	'Allowed', 
							'free-parameters'			=>	4, 
							'null-models'				=> 	['m7', 'm8a']
						},
	
	'm8a' 			=> 	{
							'category'					=>	'Site-specific', 
							'positive-site-capability'	=>	'Not Allowed', 
							'free-parameters'			=>	4, 
							'fixed-omega'			=>	1	
						},
	
	'modelA'		=> 	{ 
							'category'					=>	'Branch-site', 
							'positive-site-capability'	=>	'Allowed', 
							'free-parameters'			=>	3, 
							'null-models'				=>	['m1Neutral', 'modelAnull']
						},
	
	'modelAnull'	=> 	{ 
							'category'					=>	'Branch-site', 
							'positive-site-capability'	=>	'Not Allowed', 
							'free-parameters'			=>	3, 
							'fixed-omega'			=>	1	
						},
	
	'modelB'		=> 	{ 	
							'category'					=>	'Branch-site', 
							'positive-site-capability'	=>	'Allowed', 
							'free-parameters'			=>	5, 
							'null-models'				=>	['m3Discrtk2']	
						}
);
	
my %LRT_params = 
(
		'm0 vs m3Discrtk2' 		   => {'AdjFactor'=>2, 'CriticalValue'=>5.99}, 
		'm1Neutral vs 2ratios'     => {'AdjFactor'=>2, 'CriticalValue'=>5.99},
		'm3Discrtk2 vs m3Discrtk3' => {'AdjFactor'=>1, 'CriticalValue'=>1.00},
		'm1Neutral vs m2Selection' => {'AdjFactor'=>2, 'CriticalValue'=>5.99},
		'm7 vs m8' 				   => {'AdjFactor'=>2, 'CriticalValue'=>5.99},
		'm8a vs m8' 			   => {'AdjFactor'=>2, 'CriticalValue'=>2.71},
		'm1Neutral vs modelA'      => {'AdjFactor'=>2, 'CriticalValue'=>5.99},
		'modelAnull vs modelA'     => {'AdjFactor'=>2, 'CriticalValue'=>3.84},
		'm3Discrtk2 vs modelB'     => {'AdjFactor'=>2, 'CriticalValue'=>5.99}			
);	

=head2 new

 Usage   : $obj = CodemlWrapper::Job->new();
           $obj = CodemlWrapper::Job->new($Codeml_workspace);
 Function: Creates a new CodemlWrapper::Job object. If a Codeml workspace directory
           is specified, loads CodemlWrapper::Job object from that workspace.
 Returns : $obj (CodemlWrapper::Job object)
 Args    : $Codeml_workspace (optional directory path)

=cut

sub new
{
	my($self, @arguments) = @_;
	undef my %Job;
	$self = \%Job;
	
	bless($self, 'Job');
			
	my $sub = $self->get_nickname().substr((caller(0))[3], MODULE_NAME_LENGTH);		
			
	$self->reset_job();
	
	if(@arguments > 2)
	{Carp::croak($sub." failed: too many arguments (max=2), exiting;");}
	
	my %valid_settings = ('Verbose'=>1);
	my @setting_arguments = grep { defined($valid_settings{$_}) } @arguments;
	my @standard_arguments = grep { !defined($valid_settings{$_}) } @arguments;
	
	foreach my $set_arg (@setting_arguments)
	{	
		if($set_arg=~/^Verbose$/i)
		{$self->{'VERBOSE'} = TRUE;}
	}
	
	if(@standard_arguments == 1)
	{
		my ($Codeml_workspace) = @standard_arguments;
		
		$self->load($Codeml_workspace);
	}
	elsif(@standard_arguments > 1)
	{
		Carp::croak($sub." failed: invalid arguments: (\'".join("\', \'", @standard_arguments)."\'), exiting;");
	}
										
	return $self;
}

=head2 reset_job

 Usage   : $obj->reset_job();
 Function: Clears all data and resets default settings.
 Returns : none
 Args    : none

=cut

sub reset_job
{
	my($self) = @_;
	my $sub = $self->get_nickname().substr((caller(0))[3], MODULE_NAME_LENGTH);	
	
	delete $self->{'ALIGNMENT'};
	delete $self->{'TREES'};	
	delete $self->{'TASK-SET'};
	delete $self->{'RESULTS'};
	
	delete $self->{'LOADED-WORKSPACE-DIRECTORY'};
	delete $self->{'POSITIVE-SELECTION-FOUND'};
			
	%{$self->{'OMEGAS'}} = map { $_ => DEF } @default_omegas;	
	%{$self->{'MODELS'}} = map { $_ => DEF } @default_models;
		
	$self->{'TOPOLOGY-NAME'} = 'Sites';
	$self->{'POS-SITE-PROBABILITY-THRESHOLD'} = DEFAULT_PSP_THRESHOLD;
	
	print($sub.":...\n" ) if($self->{'VERBOSE'});
}

=head2 load

 Usage   : $obj->load($Codeml_workspace);
 Function: Loads Codeml workspace from the specified Codeml workspace directory, 
           processes any output.
 Returns : none
 Args    : $Codeml_workspace (directory)

=cut

sub load
{
	my($self, $dir) = @_;
	my $sub = $self->get_nickname().substr((caller(0))[3], MODULE_NAME_LENGTH);	
	
	if( !(-d $dir) )
	{Carp::croak($sub." failed: can't load directory \'".$dir."\', exiting;");}
	
	if(keys(%{$self->{'MODELS'}})==0)
	{Carp::croak($sub." failed: no model(s) set, exiting;");}
	
	if(keys(%{$self->{'OMEGAS'}})==0)
	{Carp::croak($sub." failed: no initial w-value(s) set, exiting;");}
	
	delete $self->{'ALIGNMENT'};
	delete $self->{'TREES'};	
	delete $self->{'TASK-SET'};	
	delete $self->{'POSITIVE-SELECTION-FOUND'};
	
	$self->{'LOADED-WORKSPACE-DIRECTORY'} = $dir;
	
	my $any_output = FALSE;			
	my $ctl_temp_outfile = $dir."/ctl_file_search.out";
	my $ctl_temp_errfile = $dir."/ctl_file_search.err";
	while( (-f $ctl_temp_outfile) || (-f $ctl_temp_errfile) )
	{
		my $r = int(rand(1000000));
		$ctl_temp_outfile = $dir."/ctl_file_search_".$r.".out";
		$ctl_temp_errfile = $dir."/ctl_file_search_".$r.".err";
	}
			
	my $ctl_file_search_status = system("find \"".$dir."\" -name \"codeml.ctl\" 1> ".$ctl_temp_outfile." 2> ".$ctl_temp_errfile);
	my @ctl_files = util_read_file($ctl_temp_outfile);
	chomp(@ctl_files);
	system("rm ".$ctl_temp_outfile)==0 or Carp::croak($sub." failed: can't remove temp file \'".$ctl_temp_outfile."\', exiting;");

	my $ctl_file_search_errors = util_read_file($ctl_temp_errfile);
	system("rm ".$ctl_temp_errfile)==0 or Carp::croak($sub." failed: can't remove temp file \'".$ctl_temp_errfile."\', exiting;");
			
	if($ctl_file_search_status!=0)
	{Carp::croak($sub." failed: \n".$ctl_file_search_errors);}

	if(@ctl_files==0)
	{Carp::croak($sub." failed: no Codeml control files in \'".$dir."\', exiting;");}
	
	my $next_temp_key = 1;
	undef my %temp;
			
	foreach my $ctl_file (@ctl_files)
	{	
		my $ctl_dir = substr($ctl_file, 0, rindex($ctl_file, "/")+1);
		
		my $ctl = Control->new($ctl_file);
		
		my $seq_file  = $ctl_dir.$ctl->get_parameter('seqfile');			
		my $tree_file = $ctl_dir.$ctl->get_parameter('treefile');
		my $out_file  = $ctl_dir.$ctl->get_parameter('outfile');	
		my $model 	  = $ctl->get_model();
		my $omega 	  = $ctl->get_parameter('omega');
		
		if( defined($self->{'MODELS'}{$model}) && defined($self->{'OMEGAS'}{$omega}) )
		{
			my $tree = Tree->new($tree_file);
			
			if($tree->is_rooted())
			{Carp::croak($sub." failed: rooted tree in \'".$tree_file."\', exiting;");}
			
			if( ( ($model_feature{$model}{'category'} eq 'Branch-specific') or ($model_feature{$model}{'category'} eq 'Branch-site') ) and ($tree->get_foreground_genes() == 0) )
			{Carp::croak($sub." failed: \'".$ctl_file."\' specifies \'".$model."\' but refers to unlabelled tree file \'".$tree_file."\', exiting;");}

			if( ( ($model_feature{$model}{'category'} eq 'Homogeneous') or ($model_feature{$model}{'category'} eq 'Site-specific') ) and ($tree->get_foreground_genes() > 0) )
			{Carp::croak($sub." failed: \'".$ctl_file."\' specifies \'".$model."\' but refers to labelled tree file \'".$tree_file."\', exiting;");}
			
			undef my $output;			
			if(-f $out_file)
			{
				$output = Output->new($out_file);
				$any_output = TRUE;
			}
			
			if(!$self->model_supported($model))
			{Carp::croak($sub." failed: model \'".$model."\' not supported, exiting;");}
						
			if( defined($output) && ($output->get_model() ne $model) )
			{Carp::croak($sub." failed: model given differs between Codeml control file (\'".$ctl_file."\') and output file (\'".$out_file."\'), exiting;");}
									
			if(defined($temp{'ALIGNMENT'}))
			{	
				if(!util_alignments_equal($seq_file, $temp{'ALIGNMENT'}))
				{Carp::croak($sub." failed: more than one Codeml dataset in \'".$dir."\', exiting;");}
			}
			else
			{
				$temp{'ALIGNMENT'} = Alignment->new($seq_file);
				
				if(!defined($temp{'GENE-NAMES'}))
				{
					$temp{'GENE-NAMES'} =  join("\t", sort {$a cmp $b} $temp{'ALIGNMENT'}->get_sequence_names());
				}
			}
						
			if(defined($temp{'TREES'}{ $self->{'TOPOLOGY-NAME'} }))
			{				
				if( !util_topologies_equal($tree, $temp{'TREES'}{ $self->{'TOPOLOGY-NAME'} }) )
				{Carp::croak($sub." failed: more than one Codeml dataset in \'".$dir."\', exiting;");}
			}
			else
			{
				my $topology = $tree->copy();
				$topology->remove_attributes('BRANCH-LENGTHS', 'BRANCH-LABELS');
				
				my $topology_leaves = join("\t", $topology->get_genes());
				
				if($topology_leaves ne $temp{'GENE-NAMES'})
				{Carp::croak($sub." failed: more than one Codeml dataset in \'".$dir."\', exiting;");}
				
				my $ladderised_topology_string = $topology->get_ladderised_tree_string();
				$temp{'TREE-KEYS'}{$ladderised_topology_string} = $self->{'TOPOLOGY-NAME'};
				
				$temp{'TREES'}{ $self->{'TOPOLOGY-NAME'} } = $topology;
			}
			
			undef my $tmp_tree_key;
			my $ladderised_tree_string = $tree->get_ladderised_tree_string('BRANCH-LABELS');						
			if( ($model_feature{$model}{'category'} eq 'Homogeneous') || ($model_feature{$model}{'category'} eq 'Site-specific') )
			{
				$tmp_tree_key = $self->{'TOPOLOGY-NAME'};
			}
			elsif(defined($temp{'TREE-KEYS'}{$ladderised_tree_string}))
			{
				$tmp_tree_key = $temp{'TREE-KEYS'}{$ladderised_tree_string};
			}
			else
			{
				$tmp_tree_key = $next_temp_key++;
				$temp{'TREE-KEYS'}{$ladderised_tree_string} = $tmp_tree_key;
				$temp{'TREES'}{$tmp_tree_key} = $tree;
			}
			
			if(defined($temp{'TASKS'}{$tmp_tree_key}{$model}{$omega}))
			{Carp::croak($sub." failed: results in \'".$ctl_dir."\' duplicate those in \'".$temp{'TASKS'}{$tmp_tree_key}{$model}{$omega}{"DIR"}."\', exiting;");}
			
			$temp{'TASKS'}{$tmp_tree_key}{$model}{$omega}{'CTL-FILE'} = $ctl_file;
			$temp{'TASKS'}{$tmp_tree_key}{$model}{$omega}{'CONTROL'} = $ctl;
			
			if(defined($output))
			{$temp{'TASKS'}{$tmp_tree_key}{$model}{$omega}{'OUTPUT'} = $output;}
			
			if(grep {$_ eq $ctl_dir} @{$temp{'DIRS'}{$tmp_tree_key}})
			{Carp::croak($sub." failed: directory \'".$ctl_dir."\' has multiple Codeml control files, exiting;");}
					
			push(@{$temp{'DIRS'}{$tmp_tree_key}}, $ctl_dir);
		}
	}
	
	if(keys(%{$temp{'TASKS'}})==0)
	{Carp::croak($sub." failed: no Codeml tasks could be loaded from \'".$dir."\' with current model and initial w-value settings, exiting;");}
	
	undef my %tmp2name;
	undef my @apparent_names;
	foreach my $tmp_key (keys(%{$temp{'TASKS'}}))
	{	
		if(@{$temp{'DIRS'}{$tmp_key}}==1)
		{Carp::croak($sub." failed: can't resolve task-set name for Codeml workspace \'".$dir."\' (just one task in task-set), exiting;");}
				
		my @dir_arrays = map { [split(/\//, $_)] } @{$temp{'DIRS'}{$tmp_key}};
		my $comparator = pop(@dir_arrays);
		undef my $apparent_name;
		
		#If only one model in task-set, check if it follows the convention used in this module.
		#If it does, use that. If it doesn't, croak.
		if(keys(%{$temp{'TASKS'}{$tmp_key}})==1)
		{
			my $model = (keys(%{$temp{'TASKS'}{$tmp_key}}))[0];
			my $CodemlJob_directory_structure = TRUE;		
			
			if( ($comparator->[0] ne $dir) || ($comparator->[2] ne $model) || ($comparator->[3]!~/Omega\d+/) || (@{$comparator}>4) )
			{$CodemlJob_directory_structure = FALSE;}
	
			if($CodemlJob_directory_structure != FALSE)
			{
				foreach my $dir_array (@dir_arrays)
				{
					if( ($dir_array->[0] ne $comparator->[0]) || ($dir_array->[1] ne $comparator->[1]) || ($dir_array->[2] ne $comparator->[2]) || ($comparator->[3]!~/Omega\d+/) || (@{$dir_array}>4) )
					{$CodemlJob_directory_structure = FALSE;}
				}
			}
			
			if($CodemlJob_directory_structure==FALSE)
			{Carp::croak($sub." failed: can't resolve task-set name for Codeml workspace \'".$dir."\' (just one model in task-set), exiting;");}
			
			$apparent_name = $comparator->[1];			
		}
		else
		{
			my $i=0;
			while( !(grep {$_->[$i] ne $comparator->[$i]} @dir_arrays) && ($i<@{$comparator}) )
			{
				$apparent_name = $comparator->[$i++];
			}
						
			if( $i==(@{$comparator}-1) )
			{Carp::croak($sub." failed: can't resolve task-set name for Codeml workspace \'".$dir."\' (no common directory for task-set), exiting;");}	
		}
		
		if(!defined($apparent_name))
		{Carp::croak($sub." failed: can't resolve task-set name for Codeml workspace \'".$dir."\' (name not inferred), exiting;");}
				
		if(grep {$_ eq $apparent_name} @apparent_names)
		{Carp::croak($sub." failed: duplicate task-set name (\'".$apparent_name."\') for Codeml workspace \'".$dir."\' , exiting;");}
			
		if(grep {$_ eq $apparent_name} keys(%model_categories))
		{Carp::croak($sub." failed: can't load Codeml workspace as task-set name \'".$apparent_name."\' clashes with a Codeml model category name, exiting;");}
	
		
		$tmp2name{$tmp_key} = $apparent_name;			
		push(@apparent_names, $apparent_name);
	}
	
	foreach my $k (keys(%tmp2name))
	{
		my $v = $tmp2name{$k};
		
		if($k eq $self->{'TOPOLOGY-NAME'})
		{$self->{'TOPOLOGY-NAME'} = $v;}
		
		$self->{'TASK-SET'}{$v} = $temp{'TASKS'}{$k};
		$self->{'TREES'}{$v} = $temp{'TREES'}{$k};
	}		
	
	$self->{'ALIGNMENT'} = $temp{'ALIGNMENT'};
	
	if($any_output)
	{
		#get all the lnLs..
		foreach my $tree (keys(%{$self->{'TASK-SET'}}))
		{
			foreach my $model (keys(%{$self->{'TASK-SET'}{$tree}}))
			{
				my $task = \%{$self->{'TASK-SET'}{$tree}{$model}};
				my @outputs = grep { defined($task->{$_}{'OUTPUT'}) } keys(%{$task});
				
				if(@outputs>0)
				{
					my @lnL_sorted_initial_omegas = sort { $task->{$a}{'OUTPUT'}->get_lnL() <=> $task->{$b}{'OUTPUT'}->get_lnL() } @outputs;
					my $w0 = pop(@lnL_sorted_initial_omegas);
					
					$self->{'RESULTS'}{$tree}{$model}{'Model'} = $model;
					$self->{'RESULTS'}{$tree}{$model}{'Tree'} = $tree;
					$self->{'RESULTS'}{$tree}{$model}{'Model Type'} = $model_feature{$model}{'category'};
					$self->{'RESULTS'}{$tree}{$model}{'p'} = $model_feature{$model}{'free-parameters'};
					$self->{'RESULTS'}{$tree}{$model}{'w (t=0)'} = $w0;	
					$self->{'RESULTS'}{$tree}{$model}{'lnL'} = $task->{$w0}{'OUTPUT'}->get_lnL();
										
					my @parameters = map {  $_."=".$task->{$w0}{'OUTPUT'}->get_parameter_value($_) } $task->{$w0}{'OUTPUT'}->get_model_parameters();	
					$self->{'RESULTS'}{$tree}{$model}{'Parameter Estimates'} = join(" ", @parameters);
				}
			}
		}
		
		#perform LRTs, get pos site info where relevant..
		foreach my $tree (keys(%{$self->{'RESULTS'}}))
		{
			$self->{'POSITIVE-SELECTION-FOUND'}{$tree} = FALSE;
			
			foreach my $model (keys(%{$self->{'RESULTS'}{$tree}}))
			{
				my $result = \%{$self->{'RESULTS'}{$tree}{$model}};
				my $w0 = $result->{'w (t=0)'};
				my $lnL = $result->{'lnL'};	
				my $output = \%{$self->{'TASK-SET'}{$tree}{$model}{$w0}{'OUTPUT'}};	
							
				$result->{'LRT Result'} = 'N/A';
				if(defined($model_feature{$model}{'null-models'}))
				{
					$result->{'LRT Result'} = $model;
					foreach my $null_model (@{$model_feature{$model}{'null-models'}})
					{
						my $null_tree = ($model_feature{$null_model}{'category'} eq 'Homogeneous') || ($model_feature{$null_model}{'category'} eq 'Site-specific') ? $self->{'TOPOLOGY-NAME'} : $tree ;
							
						if(!defined($self->{'RESULTS'}{$null_tree}{$null_model}))
						{Carp::croak($sub." failed: can't report on model \'".$model."\' without its null model \'".$null_model."\', exiting;");}
					
						my $LRT = $null_model.' vs '.$model;						
						my $null_lnL = $self->{'RESULTS'}{$null_tree}{$null_model}{'lnL'};	
						my $D = $LRT_params{$LRT}{'AdjFactor'}*($lnL-$null_lnL);
						my $c = $LRT_params{$LRT}{'CriticalValue'};
						my $resultLRT = ($D > $c) ? 'Yes' : 'No' ;
						
						$self->{'LRTs'}{$tree}{$LRT}{'Null Model lnL'} 			= $null_lnL;
						$self->{'LRTs'}{$tree}{$LRT}{'Alternative Model lnL'} 	= $lnL;
						$self->{'LRTs'}{$tree}{$LRT}{'D'} 						= $D;
						$self->{'LRTs'}{$tree}{$LRT}{'Critical Value'} 			= $c;
						$self->{'LRTs'}{$tree}{$LRT}{'Null Rejected?'} 			= $resultLRT;						
						
						if($D <= $c)
						{
							$result->{'LRT Result'} = ($result->{'LRT Result'} eq $model) ? $null_model : $result->{'LRT Result'}.", ".$null_model ;
						}
					}
				}
								
				undef my $positive_selection;			
				if($model_feature{$model}{'positive-site-capability'} ne 'Not Allowed')
				{
					if( ($result->{'LRT Result'} eq $model) && ($output->get_positive_selection_status()) )
					{		
						$self->{'POSITIVE-SELECTION-FOUND'}{$tree} = TRUE;
						$result->{'Positive Selection'} = 'Yes';

						foreach my $site ( $output->get_positive_sites() )
						{
							$result->{'Pos Sites'}{$site}{'P(w>1)'} = $output->get_pos_site_attribute($site, 'P(w>1)');
							$result->{'Pos Sites'}{$site}{'Posterior Mean (w)'} = $output->get_pos_site_attribute($site, 'Posterior Mean (w)');
							$result->{'Pos Sites Inference Type'} = $output->get_pos_site_inference_type();
						}
					}
					else
					{
						$result->{'Positive Selection'} = 'No';
					}
				}
				else
				{
					$result->{'Positive Selection'} = $model_feature{$model}{'positive-site-capability'};
				}
			}
		}
	}
	
	print($sub.": Codeml workspace loaded from \'".$dir."\'...\n" ) if($self->{'VERBOSE'});
}

=head2 save

 Usage   : $obj->save($Codeml_workspace);
 Function: Saves Codeml workspace to the specified Codeml workspace directory.
 Returns : none
 Args    : $Codeml_workspace (directory)

=cut

sub save
{
	my($self, $dir) = @_;
	my $sub = $self->get_nickname().substr((caller(0))[3], MODULE_NAME_LENGTH);	
			
	if(!defined($self->{'ALIGNMENT'}))
	{Carp::croak($sub." failed: no alignment set, exiting;");}
	
	if(!defined($self->{'TREES'}))
	{Carp::croak($sub." failed: no tree(s) set, exiting;");}
	
	if(keys(%{$self->{'MODELS'}})==0)
	{Carp::croak($sub." failed: no model(s) set, exiting;");}
	
	if(keys(%{$self->{'OMEGAS'}})==0)
	{Carp::croak($sub." failed: no initial w-value(s) set, exiting;");}
	
	if(!defined($self->{'TASK-SET'}))
	{Carp::croak($sub." failed: task-set has not been defined, exiting;");}
	
	if(-d $dir)
	{Carp::croak($sub." failed: can't overwrite existing directory \'".$dir."\', exiting;");}
	
	my %tree_foreground_info = map { $_ => join(' ', $self->{'TREES'}{$_}->get_foreground_genes()) } keys(%{$self->{'TREES'}});	
	my @unlabelled_trees = grep { ($tree_foreground_info{$_} eq '') && ($_ ne $self->{'TOPOLOGY-NAME'}) } keys(%tree_foreground_info);
	my @labelled_trees = grep { ($tree_foreground_info{$_} ne '') } keys(%tree_foreground_info);
	
	if(@unlabelled_trees>1)
	{Carp::croak($sub." failed: can't save Codeml workspace with more than one unlabelled tree, exiting;");}
			
	if(@unlabelled_trees==1)
	{
		my $new_topology = $unlabelled_trees[0];
		
		delete $self->{'TREES'}{ $self->{'TOPOLOGY-NAME'} };
		delete $self->{'TASK-SET'}{ $self->{'TOPOLOGY-NAME'} };
		delete $self->{'TASK-SET'}{ $new_topology };
		
		$self->{'TOPOLOGY-NAME'} = $new_topology;
		
		foreach my $model ($self->get_models($new_topology))
		{
			foreach my $omega ($self->get_initial_omegas($model))
			{					
				$self->{'TASK-SET'}{$new_topology}{$model}{$omega}{'CONTROL'} = Control->new();
				$self->{'TASK-SET'}{$new_topology}{$model}{$omega}{'CONTROL'}->set_model($model);
				$self->{'TASK-SET'}{$new_topology}{$model}{$omega}{'CONTROL'}->set_parameter('omega', $omega);
			}
		}
	}
	
	for(my $i=0 ; $i<(@labelled_trees-1) ; $i++)
	{			
		for(my $j=$i+1 ; $j<@labelled_trees ; $j++)
		{
			if($tree_foreground_info{ $labelled_trees[$i] } eq $tree_foreground_info{ $labelled_trees[$j] })
			{Carp::croak($sub." failed: trees \'".$i."\' and \'".$j."\' have duplicate foreground taxa, exiting;");}
		}	
	}
	
	mkdir $dir;
	chdir $dir;	
	
	foreach my $tree (keys(%{$self->{'TASK-SET'}}))
	{
		mkdir $tree;
		chdir $tree;
				
		foreach my $model (keys(%{$self->{'TASK-SET'}{$tree}}))
		{
			mkdir $model;
			chdir $model;
	
			foreach my $omega (keys(%{$self->{'TASK-SET'}{$tree}{$model}}))
			{
				mkdir 'Omega'.$omega;
				chdir 'Omega'.$omega;
											
				$self->{'TASK-SET'}{$tree}{$model}{$omega}{'CONTROL'}->save();						
				$self->{'TASK-SET'}{$tree}{$model}{$omega}{'CTL-FILE'} = $dir."/".$tree."/".$model."/Omega".$omega."/codeml.ctl";
										
				my $seqfile = $self->{'TASK-SET'}{$tree}{$model}{$omega}{'CONTROL'}->get_parameter('seqfile');
				$self->{'ALIGNMENT'}->save($seqfile);
				
				my $treefile = $self->{'TASK-SET'}{$tree}{$model}{$omega}{'CONTROL'}->get_parameter('treefile');
				$self->{'TREES'}{$tree}->save($treefile);
				
				my $outfile = $self->{'TASK-SET'}{$tree}{$model}{$omega}{'CONTROL'}->get_parameter('outfile');
				if(defined($self->{'TASK-SET'}{$tree}{$model}{$omega}{'OUTPUT'}))
				{$self->{'TASK-SET'}{$tree}{$model}{$omega}{'OUTPUT'}->save($outfile);}
				
				chdir "..";
			}
			chdir "..";
		}
		chdir "..";
	}
	
	print($sub.": Codeml workspace saved to \'".$dir."\'...\n" ) if($self->{'VERBOSE'});
}

=head2 add_alignment

 Usage   : $obj->add_alignment($alignment);
 Function: Adds the given alignment to the Codeml job. If the specified alignment
           is a CodemlWrapper::Alignment object, it's added directly to the Codeml job. 
           If it's a file path, the alignment is loaded from the specified file 
           and added to the Codeml job.
 Returns : none
 Args    : $alignment (CodemlWrapper::Alignment OR file path)

=cut

sub add_alignment
{
	my($self, $alignment) = @_;
	my $sub = $self->get_nickname().substr((caller(0))[3], MODULE_NAME_LENGTH);	
	
	if(defined($self->{'LOADED-WORKSPACE-DIRECTORY'}))
	{Carp::croak($sub." failed: can't add \'CodemlWrapper::Alignment\' after a Codeml workspace has been loaded, exiting;");}
			
	delete $self->{'ALIGNMENT'};
	
	if(ref($alignment) eq 'Alignment')
	{
		$self->{'ALIGNMENT'} = $alignment->copy();
	}
	elsif(-f $alignment)
	{
		$self->{'ALIGNMENT'} = Alignment->new($alignment);
	}
	else
	{
		Carp::croak($sub." failed: can only add \'CodemlWrapper::Alignment\' object or file, exiting;");
	}

	if( defined($self->{'TREES'}) && defined($self->{'MODELS'}) && defined($self->{'OMEGAS'}) && !defined($self->{'TASK-SET'}) )
	{
		$self->_setup_task_set();	
	}
	
	if(defined($self->{'TREES'}{ $self->{'TOPOLOGY-NAME'} }))
	{
		my $alignment_headers = join("\t", sort {$a cmp $b} $self->{'ALIGNMENT'}->get_sequence_names());
		my $topology_leaves = join("\t", $self->{'TREES'}{ $self->{'TOPOLOGY-NAME'} }->get_genes());
		
		if($alignment_headers ne $topology_leaves)
		{Carp::croak($sub." failed: gene names in \'CodemlWrapper::Alignment\' differ from those in tree(s), exiting;");}
	}
	
	print($sub.": alignment \'".$alignment."\' added...\n" ) if($self->{'VERBOSE'});
}

=head2 remove_alignment

 Usage   : $obj->remove_alignment();
 Function: Removes the current alignment from the Codeml job.
 Returns : none
 Args    : none

=cut

sub remove_alignment
{	
	my($self) = @_;	
	my $sub = $self->get_nickname().substr((caller(0))[3], MODULE_NAME_LENGTH);	
	
	if(defined($self->{'LOADED-WORKSPACE-DIRECTORY'}))
	{Carp::croak($sub." failed: can't remove \'CodemlWrapper::Alignment\' after a Codeml workspace has been loaded, exiting;");}
			
	delete $self->{'ALIGNMENT'};
	delete $self->{'TASK-SET'};
	
	print($sub.": alignment removed...\n" ) if($self->{'VERBOSE'});
}

=head2 add_tree

 Usage   : $obj->add_tree($tree_name, $tree);
 Function: Adds tree with specified tree name to the Codeml job. If the given 
           tree is a CodemlWrapper::Tree object, it's added directly to the Codeml job. 
           If it's a file path, the tree is loaded from the specified file and 
           added to the Codeml job. The only restriction on tree names is that 
           they can't clash with an existing tree name (including the topology 
           name) or Codeml model category names.
 Returns : none
 Args    : $tree_name (string), $tree (CodemlWrapper::Tree or file path)

=cut

sub add_tree
{
	my($self, $tree_name, $tree) = @_;
	my $sub = $self->get_nickname().substr((caller(0))[3], MODULE_NAME_LENGTH);	
	
	if(defined($self->{'LOADED-WORKSPACE-DIRECTORY'}))
	{Carp::croak($sub." failed: can't add \'CodemlWrapper::Tree\' after a Codeml workspace has been loaded, exiting;");}
			
	if(!defined($tree_name))
	{Carp::croak($sub." failed: no tree name given, exiting;");}
	
	if($tree_name=~/\W/)
	{Carp::croak($sub." failed: tree name \'".$tree_name."\' contains non-alphanumeric character(s), exiting;");}
		
	if($tree_name eq $self->{'TOPOLOGY-NAME'})
	{Carp::croak($sub." failed: tree name \'".$tree_name."\' is reserved for the tree topology, exiting;");}
	
	if(grep {$_ eq $tree_name} keys(%model_categories))
	{Carp::croak($sub." failed: tree name \'".$tree_name."\' can't be used as it clashes with a Codeml model category name, exiting;");}
	
	if(defined($tree_name) && defined($self->{'TREES'}{$tree_name}))
	{Carp::croak($sub." failed: tree name \'".$tree_name."\' duplicates that of a previously loaded tree, exiting;");}
		
	if(!defined($tree))
	{Carp::croak($sub." failed: no tree given, exiting;");}
	
	if( (keys(%{$self->{'TREES'}}) + 1) >= ARBITRARY_LIMIT )
	{Carp::croak($sub." failed: do you really want to add ".ARBITRARY_LIMIT." trees? Really?!");}
	
	if(ref($tree) eq 'Tree')
	{
		$self->{'TREES'}{$tree_name} = $tree->copy();
	}
	elsif(-f $tree)
	{
		$self->{'TREES'}{$tree_name} = Tree->new($tree);
	}
	else
	{
		Carp::croak($sub." failed: can only add \'CodemlWrapper::Tree\' object or file, exiting;");
	}
		
	$self->{'TREES'}{$tree_name}->remove_attributes('BRANCH-LENGTHS');
	$self->{'TREES'}{$tree_name}->unroot_tree();
				
	if(defined($self->{'TREES'}{ $self->{'TOPOLOGY-NAME'} }))
	{
		if(!util_topologies_equal($self->{'TREES'}{ $self->{'TOPOLOGY-NAME'} }, $self->{'TREES'}{$tree_name}))
		{Carp::croak($sub." failed: tree \'".$tree_name."\' has different topology from previously loaded tree(s), exiting;");}
	}
	else
	{	
		$self->{'TREES'}{ $self->{'TOPOLOGY-NAME'} } = $self->{'TREES'}{$tree_name}->copy();			
		$self->{'TREES'}{ $self->{'TOPOLOGY-NAME'} }->remove_attributes('BRANCH-LABELS');
		
		if(defined($self->{'ALIGNMENT'}))
		{
			my $alignment_headers = join("\t", sort {$a cmp $b} $self->{'ALIGNMENT'}->get_sequence_names());
			my $tree_leaves = join("\t", $self->{'TREES'}{$tree_name}->get_genes());
			
			if($tree_leaves ne $alignment_headers)
			{Carp::croak($sub." failed: gene names in \'".$tree_name."\' differ from those in alignment, exiting;");}
		}
	}
	
	if( defined($self->{'ALIGNMENT'}) && defined($self->{'MODELS'}) && defined($self->{'OMEGAS'}) )
	{
		if(defined($self->{'TASK-SET'}))
		{
			foreach my $model ($self->get_models($tree_name))
			{
				foreach my $omega ($self->get_initial_omegas($model))
				{					
					$self->{'TASK-SET'}{$tree_name}{$model}{$omega}{'CONTROL'} = Control->new();
					$self->{'TASK-SET'}{$tree_name}{$model}{$omega}{'CONTROL'}->set_model($model);
					$self->{'TASK-SET'}{$tree_name}{$model}{$omega}{'CONTROL'}->set_parameter('omega', $omega);
				}
			}			
		}
		else
		{
			$self->_setup_task_set();
		}
	}

	print($sub.": tree \'".$tree_name."\' added from \'".$tree."\'...\n" ) if($self->{'VERBOSE'});
}

=head2 remove_tree

 Usage   : $obj->remove_tree($tree_name);
 Function: Removes the tree with the given name from the Codeml job.
 Returns : none
 Args    : $tree_name (string)

=cut

sub remove_tree
{	
	my($self, $tree_name) = @_;	
	my $sub = $self->get_nickname().substr((caller(0))[3], MODULE_NAME_LENGTH);	
	
	if(defined($self->{'LOADED-WORKSPACE-DIRECTORY'}))
	{Carp::croak($sub." failed: can't remove \'CodemlWrapper::Tree\' after a Codeml workspace has been loaded, exiting;");}
	
	if(!defined($tree_name))
	{Carp::croak($sub." failed: no tree name given, exiting;");}
	
	if( defined($tree_name) && ($tree_name eq $self->{'TOPOLOGY-NAME'}) )
	{Carp::croak($sub." failed: can't remove tree topology (\'".$self->{'TOPOLOGY-NAME'}."\'), exiting;");}
	
	if(defined($tree_name) && !defined($self->{'TREES'}{$tree_name}))
	{Carp::croak($sub." failed: no tree named \'".$tree_name."\', exiting;");}
	
	delete $self->{'TREES'}{$tree_name};
	delete $self->{'TASK-SET'}{$tree_name};
	
	if( (keys(%{$self->{'TREES'}})==1) && (defined($self->{'TREES'}{ $self->{'TOPOLOGY-NAME'} })) ) 
	{
		delete $self->{'TREES'};
		delete $self->{'TASK-SET'};
	}
	
	print($sub.": tree \'".$tree_name."\' removed...\n" ) if($self->{'VERBOSE'});
}

=head2 set_models

 Usage   : $obj->set_models(@models);
 Function: Sets the given models as the active set.
 Returns : none
 Args    : @models (array of strings)

=cut

sub set_models
{
	my($self, @models) = @_;	
	my $sub = $self->get_nickname().substr((caller(0))[3], MODULE_NAME_LENGTH);	
			
	delete $self->{'MODELS'};
	
	if(@models==0)
	{Carp::croak($sub." failed: no models given, exiting;");}
	
	if(grep { !$self->model_supported($_) } @models)
	{Carp::croak($sub." failed: model not supported, exiting;");}
	
	%{$self->{'MODELS'}} = map { $_ => DEF } @models;
	
	foreach my $model (keys(%{$self->{'MODELS'}}))
	{
		if( defined($model_feature{$model}{'null-models'}) )
		{
			foreach my $null (@{$model_feature{$model}{'null-models'}})
			{
				if( !defined($self->{'MODELS'}{$null}) )
				{Carp::croak($sub." failed: can't set model \'".$model."\' without its null (\'".$null."\'), exiting;");}
			}
		}
	}
	
	if( defined($self->{'ALIGNMENT'}) &&  defined($self->{'TREES'}) && defined($self->{'OMEGAS'}) )
	{
		if(defined($self->{'TASK-SET'}))
		{
			foreach my $tree_name (keys(%{$self->{'TREES'}}))
			{					
				foreach my $model ($self->get_models($tree_name))
				{
					if( !(grep {$_ eq $model} @models) )
					{
						delete $self->{'TASK-SET'}{$tree_name}{$model};
					}
					
					if( (grep {$_ eq $model} @models) && !(defined($self->{'TASK-SET'}{$tree_name}{$model})) )
					{
						foreach my $omega ($self->get_initial_omegas($model))
						{					
							$self->{'TASK-SET'}{$tree_name}{$model}{$omega}{'CONTROL'} = Control->new();
							$self->{'TASK-SET'}{$tree_name}{$model}{$omega}{'CONTROL'}->set_model($model);
							$self->{'TASK-SET'}{$tree_name}{$model}{$omega}{'CONTROL'}->set_parameter('omega', $omega);
						}
					}					
				}		
			}	
		}
		else
		{
			$self->_setup_task_set();
		}
	}
	
	if(defined($self->{'LOADED-WORKSPACE-DIRECTORY'}))
	{$self->load($self->{'LOADED-WORKSPACE-DIRECTORY'});}	

	print($sub.": \'".join("\', \'", @models)."\'...\n" ) if($self->{'VERBOSE'});
}

=head2 set_initial_omegas

 Usage   : $obj->set_initial_omegas(@omegas);
 Function: Sets the given initial ω-values as the active set.
 Returns : none
 Args    : @omegas (array of numbers)

=cut

sub set_initial_omegas
{
	my($self, @omegas) = @_;
	my $sub = $self->get_nickname().substr((caller(0))[3], MODULE_NAME_LENGTH);	
			
	delete $self->{'OMEGAS'};
	
	if(@omegas==0)
	{Carp::croak($sub." failed: no initial w-values given, exiting;");}
	
	if(grep { $_ < 0.0 } @omegas)
	{Carp::croak($sub." failed: initial w-value out of bounds, exiting;");}
	
	%{$self->{'OMEGAS'}} = map { $_ => DEF } @omegas;
	
	if( defined($self->{'ALIGNMENT'}) &&  defined($self->{'TREES'}) && defined($self->{'MODELS'}) )
	{
		if(defined($self->{'TASK-SET'}))
		{
			foreach my $tree_name (keys(%{$self->{'TREES'}}))
			{
				foreach my $model ($self->get_models($tree_name))
				{
					foreach my $omega ($self->get_initial_omegas($model))
					{
						if(!grep {$_ eq $omega} @omegas)
						{
							delete $self->{'TASK-SET'}{$tree_name}{$model}{$omega};
						}
		
						if( (grep {$_ eq $omega} @omegas) && !(defined($self->{'TASK-SET'}{$tree_name}{$model}{$omega})) )
						{
							$self->{'TASK-SET'}{$tree_name}{$model}{$omega}{'CONTROL'} = Control->new();
							$self->{'TASK-SET'}{$tree_name}{$model}{$omega}{'CONTROL'}->set_model($model);
							$self->{'TASK-SET'}{$tree_name}{$model}{$omega}{'CONTROL'}->set_parameter('omega', $omega);
						}								
					}	
				}
			}	
		}
		else
		{
			$self->_setup_task_set();
		}
	}
	
	if(defined($self->{'LOADED-WORKSPACE-DIRECTORY'}))
	{$self->load($self->{'LOADED-WORKSPACE-DIRECTORY'});}	

	print($sub.": \'".join("\', \'", @omegas)."\'...\n" ) if($self->{'VERBOSE'});
}

=head2 set_topology_name

 Usage   : $obj->set_topology_name($topology_name);
 Function: Sets the topology tree name.
 Returns : none
 Args    : $topology_name (string)

=cut

sub set_topology_name
{
	my($self, $topology_name) = @_;
	my $sub = $self->get_nickname().substr((caller(0))[3], MODULE_NAME_LENGTH);	
	
	if(!defined($topology_name))
	{Carp::croak($sub." failed: topology name not specified, exiting;");}
		
	if(defined($self->{'TREES'}{$topology_name}))
	{Carp::croak($sub." failed: topology name \"".$topology_name."\" duplicates that of another tree, exiting;");}
	
	$self->{'TOPOLOGY-NAME'} = $topology_name;
	
	print($sub.": \'".$topology_name."\'...\n" ) if($self->{'VERBOSE'});
}

=head2 set_positive_site_probability_threshold

 Usage   : $obj->set_positive_site_probability_threshold($threshold);
 Function: Sets the probability threshold for reporting positive sites. Only 
           positively selected sites with posterior probabilities above this 
           threshold are reported.
 Returns : none
 Args    : $threshold (number)

=cut

sub set_positive_site_probability_threshold
{
	my($self, $threshold) = @_;
	my $sub = $self->get_nickname().substr((caller(0))[3], MODULE_NAME_LENGTH);	
	
	if(!defined($threshold))
	{Carp::croak($sub." failed: no positive site probability threshold given, exiting;");}
	
	if( ($threshold < 0.0) || ($threshold > 1.0) )
	{Carp::croak($sub." failed: positive site probability threshold out of bounds, exiting;");}
	
	$self->{'POS-SITE-PROBABILITY-THRESHOLD'} = $threshold;
	
	print($sub.": \'".$threshold."\'...\n" ) if($self->{'VERBOSE'});
}

=head2 create_performance_report

 Usage   : $obj->create_performance_report('OUTFILE' => $output_file);
 Function: Creates report of time used for each Codeml task.
 Returns : none
 Args    : 'OUTFILE' => $output_file (file name)

=cut

sub create_performance_report
{
	my($self, @paired_parameters) = @_;
	my $sub = $self->get_nickname().substr((caller(0))[3], MODULE_NAME_LENGTH);	
	
	if( !defined($self->{'LOADED-WORKSPACE-DIRECTORY'}) || !defined($self->{'TASK-SET'}) )
	{Carp::croak($sub." failed: no Codeml workspace loaded, exiting;");}
	
	if(@paired_parameters == 0)
	{Carp::croak($sub." failed: no parameters given, exiting;");}
	
	if(@paired_parameters % 2 != 0)
	{Carp::croak($sub." failed: invalid parameters, exiting;");}

	my %supported = ('OUTFILE' => {'min'=>1,'max'=>1});
	util_check_pairwise_parameter_constraints(\%supported);
	
	undef my %param;
	for(my $p=0 ; $p < @paired_parameters ; $p+=2)
	{
		my ($k, $v) = @paired_parameters[$p, $p+1];

		if(!defined($supported{$k}))
		{Carp::croak($sub." failed: unsupported parameter (\'".$k."\'), exiting;");}
		
		if($supported{$k}{'max'} == 1)
		{
			if(defined($param{$k}))
			{Carp::croak($sub." failed: parameter (\'".$k."\') can't be specified more than once, exiting;");}
			
			$param{$k} = $v;			
		}
		else
		{
			push(@{$param{$k}}, $v);
		}
	}
	
	foreach my $s (keys(%supported))
	{
		if(defined($param{$s}))
		{
			my @values = ($supported{$s}{'max'}==1) ? ($param{$s}) : @{$param{$s}};
			
			if(@values < $supported{$s}{'min'})
			{Carp::croak($sub." failed: too few values for parameter (\'".$s."\'), exiting;");}
			
			if(@values > $supported{$s}{'max'})
			{Carp::croak($sub." failed: too many values for parameter (\'".$s."\'), exiting;");}
		}
		elsif($supported{$s}{'min'} > 0)
		{
			Carp::croak($sub." failed: no value given for required parameter (\'".$s."\'), exiting;");
		}
	}

	my  @report_headings = @{$report_headings{'performance-report'}};
	undef my %report;
	my $index = 0;		
	
	foreach my $tree (keys(%{$self->{'TASK-SET'}}))
	{
		foreach my $model (keys(%{$self->{'TASK-SET'}{$tree}}))
		{
			foreach my $omega (keys(%{$self->{'TASK-SET'}{$tree}{$model}}))
			{
				if(defined($self->{'TASK-SET'}{$tree}{$model}{$omega}{'OUTPUT'}))
				{
					my $output = \%{$self->{'TASK-SET'}{$tree}{$model}{$omega}{'OUTPUT'}};	
				
					$report{$index}{'Program Info'} = $output->get_program_info();
					$report{$index}{'Model'} = $model;
					$report{$index}{'omega'} = $omega;
					$report{$index}{'Sequence Count'} = $output->get_sequence_count();
					$report{$index}{'Sequence Length'} = $output->get_sequence_length();
					$report{$index}{'Time Used (secs)'} = $output->get_execution_time();
					
					$index++;
				}				
			}				
		}
	}
	
	my $report_text = join("\t", @report_headings)."\n";
	foreach my $i (sort {$report{$b}{'Time Used (secs)'} <=> $report{$a}{'Time Used (secs)'}} keys(%report))
	{
		my @entry = map { $report{$i}{$_} } @report_headings;				
		$report_text .= join("\t", @entry)."\n";
	}
	
	util_write_file($param{'OUTFILE'}, $report_text);
	
	print($sub.": performance report for \'".$self->{'LOADED-WORKSPACE-DIRECTORY'}."\' saved to \'".$param{'OUTFILE'}."\'...\n" ) if($self->{'VERBOSE'});
}

=head2 create_summary_report

 Usage   : $obj->create_summary_report('OUTFILE' => $output_file);
           $obj->create_summary_report('OUTFILE' => $output_file, $parameter_name => $parameter_value, ...);
 Function: Creates a summary report for the Codeml job. Optional parameters can 
           be set as shown above, with each parameter name followed by its value.
           
           To set one or more trees of interest, include the following argument
           for each tree:
           
              'TREE' => $tree
           
           ...where $tree is the name of a tree of interest. If a tree other than
           the topology tree is specified, homogeneous and site-specific models
           are reported along with lineage-specific models for that tree. If the
           topology tree is specified, homogeneous and site-specific models are
           reported. If no tree is specified, a summary report is produced for all 
           trees in the Codeml job. 
           
           To set one or more sequences of interest, simply include the following 
           argument for each sequence:
           
              'SEQUENCE' => $sequence
           
           ...where $sequence is the name of a sequence of interest, or the string
           'FOREGROUND'. If 'FOREGROUND' is specified, positive sites are reported 
           for every sequence in the foreground. If particular sequences are given,
           positive sites are reported for those sequences. If no sequence is given,
           positive sites are reported as positions in the alignment, without amino
           acid symbols.
           
           In some cases the alignment input to Codeml differs from an original
           alignment and that original alignment has been stored in a HTML masked 
           alignment (such that sites to be deleted were masked out). In this case,
           positive site info can be output in terms of the original alignment by
           including the following argument:
           
              'ORIGINAL-ALIGNMENT' => $file
           
           ...where $file is the name of the masked alignment file. Note that this
           will only work with masked alignments created by the 'Codeml::Alignment'
           module.
 Returns : none
 Args    : 'OUTFILE' => $output_file (file name), $parameter_name => $parameter_value (optional paired parameters)

=cut

sub create_summary_report
{
	my($self, @paired_parameters) = @_;	
	my $sub = $self->get_nickname().substr((caller(0))[3], MODULE_NAME_LENGTH);	

	my $threshold = $self->{'POS-SITE-PROBABILITY-THRESHOLD'};

	if( !defined($self->{'LOADED-WORKSPACE-DIRECTORY'}) || !defined($self->{'TASK-SET'}) )
	{Carp::croak($sub." failed: no Codeml workspace loaded, exiting;");}
	
	if(@paired_parameters == 0)
	{Carp::croak($sub." failed: no parameters given, exiting;");}
	
	if(@paired_parameters % 2 != 0)
	{Carp::croak($sub." failed: invalid parameters, exiting;");}
	
	my $tree_num = $self->get_tree_names();
	my $seq_num = $self->{'ALIGNMENT'}->get_sequence_names();
	my %supported = ('OUTFILE' => {'min'=>1,'max'=>1}, 'TREE' => {'min'=>0,'max'=>$tree_num}, 'SEQUENCE' => {'min'=>0,'max'=>$seq_num}, 'ORIGINAL-ALIGNMENT' => {'min'=>0,'max'=>1});
	util_check_pairwise_parameter_constraints(\%supported);
	
	undef my %param;
	for(my $p=0 ; $p < @paired_parameters ; $p+=2)
	{
		my ($k, $v) = @paired_parameters[$p, $p+1];

		if(!defined($supported{$k}))
		{Carp::croak($sub." failed: unsupported parameter (\'".$k."\'), exiting;");}
		
		if($supported{$k}{'max'} == 1)
		{
			if(defined($param{$k}))
			{Carp::croak($sub." failed: parameter (\'".$k."\') can't be specified more than once, exiting;");}
			
			$param{$k} = $v;			
		}
		else
		{
			push(@{$param{$k}}, $v);
		}
	}
	
	foreach my $s (keys(%supported))
	{
		if(defined($param{$s}))
		{
			my @values = ($supported{$s}{'max'}==1) ? ($param{$s}) : @{$param{$s}};
			
			if(@values < $supported{$s}{'min'})
			{Carp::croak($sub." failed: too few values for parameter (\'".$s."\'), exiting;");}
			
			if(@values > $supported{$s}{'max'})
			{Carp::croak($sub." failed: too many values for parameter (\'".$s."\'), exiting;");}
		}
		elsif($supported{$s}{'min'} > 0)
		{
			Carp::croak($sub." failed: no value given for required parameter (\'".$s."\'), exiting;");
		}
	}
	
	undef my @trees;	
	if(defined($param{'TREE'}))
	{
		@trees = ($self->{'TOPOLOGY-NAME'});	
		my @non_tops = grep { $_ ne $self->{'TOPOLOGY-NAME'} } @{$param{'TREE'}};
		
		foreach my $tree (@non_tops)
		{
			if( !defined($self->{'TREES'}{$tree}) )
			{Carp::croak($sub." failed: no tree with name \'".$tree."\', exiting;");}	
		
			push(@trees, $tree);
		}
	}
	else
	{
		@trees = $self->get_tree_names();
	}	

	undef my $original_alignment;	
	my $alignment_handle = $self->{'ALIGNMENT'};
	if(defined($param{'ORIGINAL-ALIGNMENT'}))
	{
		$original_alignment = Alignment->new($param{'ORIGINAL-ALIGNMENT'});
		$alignment_handle = $original_alignment;
				
		my @mask_sites = $original_alignment->get_mask_sites();
		my $edited_alignment = $original_alignment->copy();
		$edited_alignment->remove_sites(@mask_sites);
		
		if( !util_alignments_equal($edited_alignment, $self->{'ALIGNMENT'}) )
		{Carp::croak($sub." failed: alignment in \'".$param{'ORIGINAL-ALIGNMENT'}."\' doesn't match this dataset, exiting;");}
	}

	my @seqnames = $alignment_handle->get_sequence_names();
	my $num_codons = $alignment_handle->get_codon_count();
	my $num_seqs = $alignment_handle->get_sequence_count();

	my @report_headings = @{$report_headings{'summary-report'}};
	my $pos_site_heading = "Positively Selected Sites in ".$num_seqs."x".$num_codons." Alignment (P(w>1) > ".$threshold."\)";
	push(@report_headings, $pos_site_heading) if(grep {$self->{'POSITIVE-SELECTION-FOUND'}{$_}} @trees);
				
	my $report_text = join("\t", @report_headings)."\n";
	foreach my $tree (@trees)
	{
		undef my @sequences;
		if(defined($param{'SEQUENCE'}))
		{
			undef my %nr_seqs;
			foreach my $seq (@{$param{'SEQUENCE'}})
			{
				if($seq=~/^FOREGROUND$/i)
				{
					foreach ($self->{'TREES'}{$tree}->get_foreground_genes())
					{$nr_seqs{$_} = DEF;}
				}
				else
				{
					if( !(grep {$_ eq $seq} @seqnames) )
					{Carp::croak($sub." failed: no sequence with name \'".$seq."\', exiting;");}					
				
					$nr_seqs{$seq} = DEF;
				}
			}		
			
			@sequences = sort {$a cmp $b} keys(%nr_seqs);
		}
		
		if(@sequences == 0)
		{
			@sequences = ('Alignment');
		}
		
		foreach my $model ($self->get_models($tree))
		{		
			if(defined($self->{'RESULTS'}{$tree}{$model}))
			{
				undef my @pos_site_info_list;
				my $result = \%{$self->{'RESULTS'}{$tree}{$model}};
				if( defined($result->{'Positive Selection'}) && ($result->{'Positive Selection'} eq 'Yes') && ($model_feature{$model}{'positive-site-capability'} eq 'Allowed') && defined($result->{'Pos Sites'}) )
				{
					my $inf_type = $result->{'Pos Sites Inference Type'};
					my @complete_sites = sort {$a <=> $b} keys(%{$result->{'Pos Sites'}});
					
					my @qualifying_sites = grep { $result->{'Pos Sites'}{$_}{'P(w>1)'} >= $threshold } @complete_sites;
					
					if(defined($original_alignment))
					{
						#$seqnames[0] is taken as an arbitrary sequence..it shouldn't matter, since mask sites should be the same in all sequences.
						@qualifying_sites = $original_alignment->map_sites('UNMASKED', 'ALIGNMENT', $seqnames[0], @qualifying_sites);
					}
					
					foreach my $sequence (@sequences)
					{		
						my $gap_count_info = '';
						undef my @positive_site_data;
						
						if(@qualifying_sites > 0)
						{
							if($sequence ne 'Alignment')
							{
								foreach my $site (@qualifying_sites)
								{
									my $pos = $alignment_handle->map_sites('ALIGNMENT', 'UNGAPPED', $sequence, $site);
									my $res = $alignment_handle->get_residue($sequence, $site);
									push(@positive_site_data, $pos.$res) if(defined($pos));
								}
								
								my $gap_count = @qualifying_sites - @positive_site_data;
								$gap_count_info = ($gap_count > 0) ? ", ".$gap_count." gap(s)" : "";
							}
							else
							{
								@positive_site_data = @qualifying_sites;
							}
						}
						
						my $pos_site_info = (@positive_site_data > 0) ? ": ".join(" ", @positive_site_data) : '';
						my $inf_type_info = (@positive_site_data > 0) ? " ".$inf_type : '';
						push(@pos_site_info_list, $sequence." (".@positive_site_data.$inf_type_info." sites".$gap_count_info.")".$pos_site_info);
					}
				}
				
				$result->{$pos_site_heading} = join("; ", @pos_site_info_list);

				my @entry = map { $result->{$_} } @report_headings;				
				$report_text .= join("\t", @entry)."\n";
			}
		}
	}	
	
	util_write_file($param{'OUTFILE'}, $report_text);

	print($sub.": summary report for \'".$self->{'LOADED-WORKSPACE-DIRECTORY'}."\' saved to \'".$param{'OUTFILE'}."\'...\n" ) if($self->{'VERBOSE'});
}

=head2 create_pos_site_report

 Usage   : $obj->create_pos_site_report('OUTFILE' => $output_file, 'MODEL' => $model);
           $obj->create_pos_site_report('OUTFILE' => $output_file, 'MODEL' => $model, $parameter_name => $parameter_value, ...);
 Function: Creates a report of positive site info under a specified model. To 
           specify the model of interest, include the following required
           argument:
           
              'MODEL' => $model
           
 		   ...where $model is the name of a supported Codeml model. If that model
 		   is present in more than one tree, the tree must also be specified as 
 		   follows:
 		   
 		      'TREE' => $tree
 
           ...where $tree is the name of the tree associated with the model of 
           interest. Pos sites are reported in terms of the alignment by default.
           To set pos sites to be reported in terms of a particular sequence, set
           the following argument:
           
              'SEQUENCE' => $sequence
        
           ...where $sequence is the name of the sequence of interest. 
           
           In some cases the alignment input to Codeml differs from an original
           alignment and that original alignment has been stored in a HTML masked 
           alignment (such that sites to be deleted were masked out). In this case,
           pos sites can be reported in terms of the original alignment by
           including the following argument:
           
              'ORIGINAL-ALIGNMENT' => $file
           
           ...where $file is the name of the masked alignment file. Note that this
           will only work with masked alignments created by the 'Codeml::Alignment'
           module.
 Returns : none
 Args    : 'OUTFILE' => $output_file (file name), 'MODEL' => $model (model name), $parameter_name => $parameter_value (optional paired parameters)

=cut

sub create_pos_site_report
{
	my($self, @paired_parameters) = @_;
	my $sub = $self->get_nickname().substr((caller(0))[3], MODULE_NAME_LENGTH);	

	my $threshold = $self->{'POS-SITE-PROBABILITY-THRESHOLD'};

	if( !defined($self->{'LOADED-WORKSPACE-DIRECTORY'}) || !defined($self->{'TASK-SET'}) )
	{Carp::croak($sub." failed: no Codeml workspace loaded, exiting;");}
	
	if(@paired_parameters == 0)
	{Carp::croak($sub." failed: no parameters given, exiting;");}
	
	if(@paired_parameters % 2 != 0)
	{Carp::croak($sub." failed: invalid parameters, exiting;");}
	
	my %supported = ('OUTFILE' => {'min'=>1,'max'=>1}, 'TREE' => {'min'=>0,'max'=>1}, 'MODEL' => {'min'=>1,'max'=>1}, 'SEQUENCE' => {'min'=>0,'max'=>1}, 'ORIGINAL-ALIGNMENT' => {'min'=>0,'max'=>1});
	util_check_pairwise_parameter_constraints(\%supported);
	
	undef my %param;
	for(my $p=0 ; $p < @paired_parameters ; $p+=2)
	{
		my ($k, $v) = @paired_parameters[$p, $p+1];

		if(!defined($supported{$k}))
		{Carp::croak($sub." failed: unsupported parameter (\'".$k."\'), exiting;");}
		
		if($supported{$k}{'max'} == 1)
		{
			if(defined($param{$k}))
			{Carp::croak($sub." failed: parameter (\'".$k."\') can't be specified more than once, exiting;");}
			
			$param{$k} = $v;			
		}
		else
		{
			push(@{$param{$k}}, $v);
		}
	}
	
	foreach my $s (keys(%supported))
	{
		if(defined($param{$s}))
		{
			my @values = ($supported{$s}{'max'}==1) ? ($param{$s}) : @{$param{$s}};
			
			if(@values < $supported{$s}{'min'})
			{Carp::croak($sub." failed: too few values for parameter (\'".$s."\'), exiting;");}
			
			if(@values > $supported{$s}{'max'})
			{Carp::croak($sub." failed: too many values for parameter (\'".$s."\'), exiting;");}
		}
		elsif($supported{$s}{'min'} > 0)
		{
			Carp::croak($sub." failed: no value given for required parameter (\'".$s."\'), exiting;");
		}
	}
	
	my $model = $param{'MODEL'};
	
	if($model_feature{$model}{'positive-site-capability'} ne 'Allowed')
	{Carp::croak($sub." failed: model \'".$model."\' doesn't support pos sites, exiting;");}

	undef my $tree;
	if( ($model_feature{$model}{'category'} eq 'Homogeneous') || ($model_feature{$model}{'category'} eq 'Site-specific') )
	{
		$tree = $self->{'TOPOLOGY-NAME'};
	}
	elsif(keys(%{$self->{'TREES'}}) == 2)
	{
		$tree = ( grep { $_ ne $self->{'TOPOLOGY-NAME'} } keys(%{$self->{'TREES'}}) ) [0];
	}
	elsif(defined($param{'TREE'}))
	{
		$tree = $param{'TREE'};
	}
	else
	{
		Carp::croak($sub." failed: no tree specified, exiting;");
	}	

	undef my $original_alignment;	
	my $alignment_handle = $self->{'ALIGNMENT'};
	if(defined($param{'ORIGINAL-ALIGNMENT'}))
	{
		$original_alignment = Alignment->new($param{'ORIGINAL-ALIGNMENT'});
		$alignment_handle = $original_alignment;
				
		my @mask_sites = $original_alignment->get_mask_sites();
		my $edited_alignment = $original_alignment->copy();
		$edited_alignment->remove_sites(@mask_sites);
		
		if( !util_alignments_equal($edited_alignment, $self->{'ALIGNMENT'}) )
		{Carp::croak($sub." failed: alignment in \'".$param{'ORIGINAL-ALIGNMENT'}."\' doesn't match this dataset, exiting;");}
	}	

	my @seqnames = $alignment_handle->get_sequence_names();
	my $num_codons = $alignment_handle->get_codon_count();
	my $num_seqs = $alignment_handle->get_sequence_count();	
	
	my $sequence = 'Alignment';
	if(defined($param{'SEQUENCE'}))
	{
		if( !(grep {$_ eq $param{'SEQUENCE'}} @seqnames) )
		{Carp::croak($sub." failed: no sequence with name \'".$param{'SEQUENCE'}."\', exiting;");}	
		
		$sequence = $param{'SEQUENCE'};
	}
	
	if(	!defined($self->{'RESULTS'}{$tree}{$model}) )
	{Carp::croak($sub." failed: no results found for model \'".$model."\' in tree \'".$tree."\', exiting;");}

	my $result = \%{$self->{'RESULTS'}{$tree}{$model}};
	my $inf_type = $result->{'Pos Sites Inference Type'};	
	
	my $site_heading = $inf_type." Site in ".$num_seqs."x".$num_codons." Alignment";
	substr($site_heading, 12, 0) =  "sequence \'".$sequence."\' of " if($sequence ne 'Alignment');	
	my $mean_heading = "P(w>1) > ".$threshold;
	
	my @report_headings = ($site_heading);
	push(@report_headings, 'Amino Acid') if($sequence ne 'Alignment');
	push(@report_headings, $mean_heading);
	
	undef my %report;

	if( defined($result->{'Positive Selection'}) && ($result->{'Positive Selection'} eq 'Yes') && ($model_feature{$model}{'positive-site-capability'} eq 'Allowed') && defined($result->{'Pos Sites'}) )
	{
		my @complete_sites = sort {$a <=> $b} keys(%{$result->{'Pos Sites'}});

		my @qualifying_sites = grep { $result->{'Pos Sites'}{$_}{'P(w>1)'}  >= $threshold } @complete_sites;	
		
		foreach my $site (@qualifying_sites)
		{
			my $pos = $site;
			
			if( defined($original_alignment) )
			{
				#$seqnames[0] is taken as an arbitrary sequence..it shouldn't matter, since mask sites should be the same in all sequences.
				$pos = $alignment_handle->map_sites('UNMASKED', 'ALIGNMENT', $seqnames[0], $pos);
			}
			
			if($sequence ne 'Alignment')
			{
				$pos = $alignment_handle->map_sites('ALIGNMENT', 'UNGAPPED', $sequence, $pos);
			}
			
			my $coord = $site;
			if( defined($original_alignment) )
			{
				#$seqnames[0] is taken as an arbitrary sequence..it shouldn't matter, since mask sites should be the same in all sequences.
				$coord = $original_alignment->map_sites('UNMASKED', 'ALIGNMENT', $seqnames[0], $coord);
			}			
			
			if( defined($pos) )
			{
				$report{$site}{ $site_heading } = $pos;			
				$report{$site}{ $mean_heading } = $result->{'Pos Sites'}{$site}{'P(w>1)'};
				$report{$site}{'Posterior Mean (w)'} = $result->{'Pos Sites'}{$site}{'Posterior Mean (w)'};
				$report{$site}{'Amino Acid'} = $alignment_handle->get_residue($sequence, $coord) if($sequence ne 'Alignment');			
			}
		}
	}

	if(keys(%report) == 0)
	{
		print("Pos site report failed: no positive sites in \'".$sequence."\' sequence with (P(w>1) >= ".$threshold.") for model \'".$model."\' in tree \'".$tree."\'...\n");
		return;	
	}
	
	if( grep { defined($report{$_}{'Posterior Mean (w)'}) && $report{$_}{'Posterior Mean (w)'} ne '' } keys(%report) )
	{push(@report_headings, 'Posterior Mean (w)');}
	
	my $report_text = join("\t", @report_headings)."\n";
	foreach my $site (sort {$a <=> $b} keys(%report))
	{
		my @entry = map { $report{$site}{$_} } @report_headings;				
		$report_text .= join("\t", @entry)."\n";
	}
			
	util_write_file($param{'OUTFILE'}, $report_text);
	
	print($sub.": pos site report for \'".$self->{'LOADED-WORKSPACE-DIRECTORY'}."\' saved to \'".$param{'OUTFILE'}."\'...\n" ) if($self->{'VERBOSE'});
}

=head2 create_pos_site_msa

 Usage   : $obj->create_pos_site_msa('OUTFILE' => $output_file, 'MODEL' => $model);
           $obj->create_pos_site_msa('OUTFILE' => $output_file, 'MODEL' => $model, $parameter_name => $parameter_value, ...);
 Function: Creates a HTML alignment showing positive site info under a specified 
           model. To specify the model of interest, include the following required
           argument:
           
              'MODEL' => $model
           
 		   ...where $model is the name of a supported Codeml model. If that model
 		   is present in more than one tree, the tree must also be specified as 
 		   follows:
 		   
 		      'TREE' => $tree
 
           ...where $tree is the name of the tree associated with the model of 
           interest. Pos sites are reported in terms of the alignment by default.
           To set pos sites to be reported in terms of a particular sequence, set
           the following argument:
           
              'SEQUENCE' => $sequence
        
           ...where $sequence is the name of the sequence of interest. 
           
           In some cases the alignment input to Codeml differs from an original
           alignment and that original alignment has been stored in a HTML masked 
           alignment (such that sites to be deleted were masked out). In this case,
           the original alignment can be saved with pos site info by including the
           following argument:
           
              'ORIGINAL-ALIGNMENT' => $file
           
           ...where $file is the name of the masked alignment file. Note that this
           will only work with masked alignments created by the 'Codeml::Alignment'
           module.
 Returns : none
 Args    : 'OUTFILE' => $output_file (file name), 'MODEL' => $model (model name), $parameter_name => $parameter_value (optional paired parameters)

=cut

sub create_pos_site_msa
{
	my($self, @paired_parameters) = @_;
	my $sub = $self->get_nickname().substr((caller(0))[3], MODULE_NAME_LENGTH);	
	
	my $threshold = $self->{'POS-SITE-PROBABILITY-THRESHOLD'};

	if( !defined($self->{'LOADED-WORKSPACE-DIRECTORY'}) || !defined($self->{'TASK-SET'}) )
	{Carp::croak($sub." failed: no Codeml workspace loaded, exiting;");}
	
	if(@paired_parameters == 0)
	{Carp::croak($sub." failed: no parameters given, exiting;");}
	
	if(@paired_parameters % 2 != 0)
	{Carp::croak($sub." failed: invalid parameters, exiting;");}

	my %supported = ('OUTFILE' => {'min'=>1,'max'=>1}, 'TREE' => {'min'=>0,'max'=>1}, 'MODEL' => {'min'=>1,'max'=>1}, 'ORIGINAL-ALIGNMENT' => {'min'=>0,'max'=>1});
	util_check_pairwise_parameter_constraints(\%supported);
	
	undef my %param;
	for(my $p=0 ; $p < @paired_parameters ; $p+=2)
	{
		my ($k, $v) = @paired_parameters[$p, $p+1];

		if(!defined($supported{$k}))
		{Carp::croak($sub." failed: unsupported parameter (\'".$k."\'), exiting;");}
		
		if($supported{$k}{'max'} == 1)
		{
			if(defined($param{$k}))
			{Carp::croak($sub." failed: parameter (\'".$k."\') can't be specified more than once, exiting;");}
			
			$param{$k} = $v;			
		}
		else
		{
			push(@{$param{$k}}, $v);
		}
	}
	
	foreach my $s (keys(%supported))
	{
		if(defined($param{$s}))
		{
			my @values = ($supported{$s}{'max'}==1) ? ($param{$s}) : @{$param{$s}};
			
			if(@values < $supported{$s}{'min'})
			{Carp::croak($sub." failed: too few values for parameter (\'".$s."\'), exiting;");}
			
			if(@values > $supported{$s}{'max'})
			{Carp::croak($sub." failed: too many values for parameter (\'".$s."\'), exiting;");}
		}
		elsif($supported{$s}{'min'} > 0)
		{
			Carp::croak($sub." failed: no value given for required parameter (\'".$s."\'), exiting;");
		}
	}
	
	my $model = $param{'MODEL'};

	if($model_feature{$model}{'positive-site-capability'} ne 'Allowed')
	{Carp::croak($sub." failed: model \'".$model."\' doesn't support pos sites, exiting;");}

	undef my $tree;
	if( ($model_feature{$model}{'category'} eq 'Homogeneous') || ($model_feature{$model}{'category'} eq 'Site-specific') )
	{
		$tree = $self->{'TOPOLOGY-NAME'};
	}
	elsif(keys(%{$self->{'TREES'}}) == 2)
	{
		$tree = ( grep { $_ ne $self->{'TOPOLOGY-NAME'} } keys(%{$self->{'TREES'}}) ) [0];
	}
	elsif(defined($param{'TREE'}))
	{
		$tree = $param{'TREE'};
	}
	else
	{
		Carp::croak($sub." failed: no tree specified, exiting;");
	}	

	undef my $original_alignment;	
	my $alignment_handle = $self->{'ALIGNMENT'};
	if(defined($param{'ORIGINAL-ALIGNMENT'}))
	{
		$original_alignment = Alignment->new($param{'ORIGINAL-ALIGNMENT'});
		$alignment_handle = $original_alignment;
				
		my @mask_sites = $original_alignment->get_mask_sites();
		my $edited_alignment = $original_alignment->copy();
		$edited_alignment->remove_sites(@mask_sites);
		
		if( !util_alignments_equal($edited_alignment, $self->{'ALIGNMENT'}) )
		{Carp::croak($sub." failed: alignment in \'".$param{'ORIGINAL-ALIGNMENT'}."\' doesn't match this dataset, exiting;");}
	}	
		
	my @seqnames = $alignment_handle->get_sequence_names();
	
	if(	!defined($self->{'RESULTS'}{$tree}{$model}) )
	{Carp::croak($sub." failed: no results found for model \'".$model."\' in tree \'".$tree."\', exiting;");}

	undef my @foreground_genes;
	my $site_count = 0;
	my $result = \%{$self->{'RESULTS'}{$tree}{$model}};
	if( defined($result->{'Positive Selection'}) && ($result->{'Positive Selection'} eq 'Yes') && ($model_feature{$model}{'positive-site-capability'} eq 'Allowed') && defined($result->{'Pos Sites'}) )
	{
		my @complete_sites = sort {$a <=> $b} keys(%{$result->{'Pos Sites'}});
		my @qualifying_sites = grep { $result->{'Pos Sites'}{$_}{'P(w>1)'} >= $threshold } @complete_sites;	

		if(defined($original_alignment))
		{
			#$seqnames[0] is taken as an arbitrary sequence..it shouldn't matter, since mask sites should be the same in all sequences.
			@qualifying_sites = $original_alignment->map_sites('UNMASKED', 'ALIGNMENT', $seqnames[0], @qualifying_sites);
		}
	
		@foreground_genes = $self->{'TREES'}{$tree}->get_foreground_genes();
		
		$alignment_handle->set_foreground_sequences(@foreground_genes) if(@foreground_genes > 0);;
		
		$alignment_handle->set_pos_sites(@qualifying_sites);
		
		$site_count = scalar(@qualifying_sites);
	}

	if($site_count == 0)
	{
		print("Pos site MSA failed: no positive sites at (P(w>1) >= ".$threshold.") for model \'".$model."\' in tree \'".$tree."\'...\n");
		return;
	}
	
	$alignment_handle->save($param{'OUTFILE'}, 'HTML');	
	
	$alignment_handle->set_foreground_sequences() if(@foreground_genes > 0);;
	
	print($sub.": pos site MSA for \'".$self->{'LOADED-WORKSPACE-DIRECTORY'}."\' saved to \'".$param{'OUTFILE'}."\'...\n" ) if($self->{'VERBOSE'});
}

=head2 create_lrt_report

 Usage   : $obj->create_lrt_report('OUTFILE' => $output_file);
 Function: Creates report of LRTs in the Codeml job.
 Returns : none
 Args    : 'OUTFILE' => $output_file (file name)

=cut

sub create_lrt_report
{
	my($self, @paired_parameters) = @_;
	my $sub = $self->get_nickname().substr((caller(0))[3], MODULE_NAME_LENGTH);	
	
	if( !defined($self->{'LOADED-WORKSPACE-DIRECTORY'}) || !defined($self->{'TASK-SET'}) )
	{Carp::croak($sub." failed: no Codeml workspace loaded, exiting;");}
	
	if(@paired_parameters == 0)
	{Carp::croak($sub." failed: no parameters given, exiting;");}
	
	if(@paired_parameters % 2 != 0)
	{Carp::croak($sub." failed: invalid parameters, exiting;");}
	
	my %supported = ('OUTFILE' => {'min'=>1,'max'=>1});
	util_check_pairwise_parameter_constraints(\%supported);
	
	undef my %param;
	for(my $p=0 ; $p < @paired_parameters ; $p+=2)
	{
		my ($k, $v) = @paired_parameters[$p, $p+1];

		if(!defined($supported{$k}))
		{Carp::croak($sub." failed: unsupported parameter (\'".$k."\'), exiting;");}
		
		if($supported{$k}{'max'} == 1)
		{
			if(defined($param{$k}))
			{Carp::croak($sub." failed: parameter (\'".$k."\') can't be specified more than once, exiting;");}
			
			$param{$k} = $v;			
		}
		else
		{
			push(@{$param{$k}}, $v);
		}
	}
	
	foreach my $s (keys(%supported))
	{
		if(defined($param{$s}))
		{
			my @values = ($supported{$s}{'max'}==1) ? ($param{$s}) : @{$param{$s}};
			
			if(@values < $supported{$s}{'min'})
			{Carp::croak($sub." failed: too few values for parameter (\'".$s."\'), exiting;");}
			
			if(@values > $supported{$s}{'max'})
			{Carp::croak($sub." failed: too many values for parameter (\'".$s."\'), exiting;");}
		}
		elsif($supported{$s}{'min'} > 0)
		{
			Carp::croak($sub." failed: no value given for required parameter (\'".$s."\'), exiting;");
		}
	}
	
	my @report_headings = @{$report_headings{'lrt-report'}};
	
	undef my %report;
	my $index = 0;		
	
	foreach my $tree ($self->get_tree_names())
	{
		if(defined($self->{'LRTs'}{$tree}))
		{
			foreach my $LRT (sort {$a cmp $b} keys(%{$self->{'LRTs'}{$tree}}))
			{				
				$report{$index}{'Tree'} = $tree;
				$report{$index}{"LRT"} = $LRT;
				
				$report{$index}{'Null Model lnL'} 		 = $self->{'LRTs'}{$tree}{$LRT}{'Null Model lnL'};
				$report{$index}{'Alternative Model lnL'} = $self->{'LRTs'}{$tree}{$LRT}{'Alternative Model lnL'};
				$report{$index}{'D'} 					 = $self->{'LRTs'}{$tree}{$LRT}{'D'};
				$report{$index}{'Critical Value'} 		 = $self->{'LRTs'}{$tree}{$LRT}{'Critical Value'};
				$report{$index}{'Null Rejected?'} 		 = $self->{'LRTs'}{$tree}{$LRT}{'Null Rejected?'};
							
				$index++;
			}
		}
	}
	
	my $report_text = join("\t", @report_headings)."\n";
	foreach my $i (sort {$a <=> $b} keys(%report))
	{
		my @entry = map { $report{$i}{$_} } @report_headings;				
		$report_text .= join("\t", @entry)."\n";
	}
	
	util_write_file($param{'OUTFILE'}, $report_text);

	print($sub.": LRT report for \'".$self->{'LOADED-WORKSPACE-DIRECTORY'}."\' saved to \'".$param{'OUTFILE'}."\'...\n" ) if($self->{'VERBOSE'});
}

=head2 get_tree_names

 Usage   : @tree_names = $obj->get_tree_names();
 Function: Gets a list of all trees in the Codeml job, with the topology tree 
           first and the remaining trees sorted in alphabetical order.
 Returns : @tree_names (array of strings)
 Args    : none

=cut

sub get_tree_names
{
	my($self) = @_;		
	my $sub = $self->get_nickname().substr((caller(0))[3], MODULE_NAME_LENGTH);	
	
	undef my @tree_list;
	
	if(keys(%{$self->{'TREES'}})>0)
	{
		my @non_topology_set = grep { $_ ne $self->{'TOPOLOGY-NAME'} } keys(%{$self->{'TREES'}});
		my @sorted_non_topology_set = sort {$a cmp $b} @non_topology_set;
		@tree_list = ( $self->{'TOPOLOGY-NAME'}, @sorted_non_topology_set );
	}
			
	print($sub.": \'".join("\', \'", @tree_list)."\'...\n" ) if($self->{'VERBOSE'} && ( !defined((caller(1))[3]) || (substr((caller(1))[3], 0, 5) ne 'Job::') ) );
			
	return @tree_list;
}

=head2 get_foreground

 Usage   : @foreground_taxa = $obj->get_foreground($tree);
 Function: Gets a list of foreground taxa for the specified tree in the Codeml job.
 Returns : @foreground_taxa (array of strings)
 Args    : $tree (string)

=cut

sub get_foreground
{
	my($self, $tree) = @_;		
	my $sub = $self->get_nickname().substr((caller(0))[3], MODULE_NAME_LENGTH);	
	
	if(!defined($tree))
	{Carp::croak($sub." failed: no tree name given, exiting;");}
	
	if(defined($tree) && !defined($self->{'TREES'}{$tree}))
	{Carp::croak($sub." failed: no tree named \'".$tree."\', exiting;");}
	
	my @foreground_taxa = $self->{'TREES'}{$tree}->get_foreground_genes();
			
	print($sub.": foreground of \'".$tree."\' - \'".join("\', \'", @foreground_taxa)."\'...\n" ) if($self->{'VERBOSE'} && ( !defined((caller(1))[3]) || (substr((caller(1))[3], 0, 5) ne 'Job::') ) );
			
	return @foreground_taxa;
}

=head2 get_models

 Usage   : @models = $obj->get_models();
           @models = $obj->get_models($grouping);
 Function: Gets a list of all active models, or if a grouping is given, models for 
           that grouping. A grouping can either be a tree name or a model category 
           name. If it's a tree name, active models associated with that tree are 
           returned. If it's a model category name, active models in that category
           are returned.
 Returns : @models (array of strings)
 Args    : $grouping (optional string)

=cut

sub get_models
{
	my($self, $grouping) = @_;
	my $sub = $self->get_nickname().substr((caller(0))[3], MODULE_NAME_LENGTH);	
	
	undef my @models;
	
	if(defined($grouping))
	{
		if(grep {$_ eq $grouping} keys(%model_categories))
		{
			@models = @{$model_categories{$grouping}};
		}
		elsif(defined($self->{'TREES'}{$grouping}))
		{
			if($grouping eq $self->{'TOPOLOGY-NAME'})
			{
				@models = grep { defined($self->{'MODELS'}{$_}) } (@{$model_categories{'Homogeneous'}}, @{$model_categories{'Site-specific'}});
			}
			else
			{
				@models = grep { defined($self->{'MODELS'}{$_}) } (@{$model_categories{"Branch-specific"}}, @{$model_categories{"Branch-site"}});
			}
		}
		else
		{
			Carp::croak($sub." failed: grouping \'".$grouping."\' is neither a tree nor a model category, exiting;");
		}
	}
	else
	{
		@models = keys(%{$self->{'MODELS'}});
	}
	
	print($sub.": \'".join("\', \'", @models)."\'...\n" ) if($self->{'VERBOSE'} && ( !defined((caller(1))[3]) || (substr((caller(1))[3], 0, 5) ne 'Job::') ) );
		
	return @models;
}

=head2 get_initial_omegas

 Usage   : @omegas = $obj->get_initial_omegas();
           @omegas = $obj->get_initial_omegas($model);
 Function: Gets a list of all active initial ω-values, or if a model is given, 
           only initial ω-values for that model.
 Returns : @omegas (array of numbers)
 Args    : $model (optional string)

=cut

sub get_initial_omegas
{
	my($self, $model) = @_;
	my $sub = $self->get_nickname().substr((caller(0))[3], MODULE_NAME_LENGTH);	
	
	undef my @initial_omegas;
	
	if( defined($model) && !$self->model_supported($model) )
	{Carp::croak($sub." failed: model \'".$model."\' not supported, exiting;");}
	
	if( defined($model) && defined($model_feature{$model}{"fixed-omega"}) )
	{
		my $fixed_omega = $model_feature{$model}{"fixed-omega"};
		
		if(defined($self->{'OMEGAS'}{$fixed_omega}))
		{
			@initial_omegas = ($fixed_omega);
		}
	}
	else
	{
		@initial_omegas = keys(%{$self->{'OMEGAS'}});
	}
	
	print($sub.": \'".join("\', \'", @initial_omegas)."\'...\n" ) if($self->{'VERBOSE'} && ( !defined((caller(1))[3]) || (substr((caller(1))[3], 0, 5) ne 'Job::') ) );
	
	return @initial_omegas;
}

=head2 get_topology_name

 Usage   : $topology_name = $obj->get_topology_name();
 Function: Gets the topology tree name.
 Returns : $topology_name (string)
 Args    : none

=cut

sub get_topology_name
{
	my($self) = @_;
	my $sub = $self->get_nickname().substr((caller(0))[3], MODULE_NAME_LENGTH);	

	print($sub.": \'".$self->{'TOPOLOGY-NAME'}."\'...\n" ) if($self->{'VERBOSE'});	
	
	return $self->{'TOPOLOGY-NAME'};
}

=head2 get_positive_site_probability_threshold

 Usage   : $threshold = $obj->get_positive_site_probability_threshold();
 Function: Gets the probability threshold for reporting positive sites.
 Returns : $threshold (number)
 Args    : none

=cut

sub get_positive_site_probability_threshold
{
	my($self) = @_;
	
	if($self->{'VERBOSE'})
	{ 
		my $sub = $self->get_nickname().substr((caller(0))[3], MODULE_NAME_LENGTH);		
		print($sub.": \'".$self->{'POS-SITE-PROBABILITY-THRESHOLD'}."\'...\n" );	
	}
			
	return $self->{'POS-SITE-PROBABILITY-THRESHOLD'};
}

=head2 get_supported_model_categories

 Usage   : @model_categories = $obj->get_supported_model_categories();
 Function: Gets a list of supported model categories.
 Returns : @model_categories (array of strings)
 Args    : none

=cut

sub get_supported_model_categories
{
	my($self) = @_;
	my $sub = $self->get_nickname().substr((caller(0))[3], MODULE_NAME_LENGTH);	

	if($self->{'VERBOSE'})
	{ 
		my $sub = $self->get_nickname().substr((caller(0))[3], MODULE_NAME_LENGTH);	
		print($sub.": \'".join("\', \'", @supported_model_categories)."\'...\n" );
	}
	
	return @supported_model_categories;
}

=head2 get_supported_model_features

 Usage   : @model_features = $obj->get_supported_model_features();
 Function: Gets a list of supported model features. 
 Returns : @model_features (array of strings)
 Args    : none

=cut

sub get_supported_model_features
{
	my($self) = @_;

	if($self->{'VERBOSE'})
	{ 
		my $sub = $self->get_nickname().substr((caller(0))[3], MODULE_NAME_LENGTH);	
		print($sub.": \'".join("\', \'", @supported_model_features)."\'...\n" );
	}
	
	return @supported_model_features;
}

=head2 get_supported_models

 Usage   : @models = $obj->get_supported_models();
 Function: Gets a list of supported models.
 Returns : @models (array of strings)
 Args    : none

=cut

sub get_supported_models
{
	my($self) = @_;
	
	if($self->{'VERBOSE'})
	{ 
		my $sub = $self->get_nickname().substr((caller(0))[3], MODULE_NAME_LENGTH);	
		print($sub.": \'".join("\', \'", @supported_models)."\'...\n" );
	}
	
	return @supported_models;
}

=head2 get_model_feature

 Usage   : $model_feature_value = $obj->get_model_feature($model, $feature);
 Function: Gets the value of a given feature for a given model.
 Returns : $model_feature_value (string OR number) 
 Args    : $model (string), $feature (string)

=cut

sub get_model_feature
{
	my($self, $model, $feature) = @_;
	my $sub = $self->get_nickname().substr((caller(0))[3], MODULE_NAME_LENGTH);	
	
	if(!defined($model))
	{Carp::croak($sub." failed: model not defined, exiting;");}
	
	if(!defined($feature))
	{Carp::croak($sub." failed: model feature not defined, exiting;");}
		
	if( !(grep { $model eq $_ } @supported_models) )
	{Carp::croak($sub." failed: model \'".$model."\' not supported, exiting;");}
	
	if( !(grep { $feature eq $_ } @supported_model_features) )
	{Carp::croak($sub." failed: model feature \'".$feature."\' not supported, exiting;");}
		
	print($sub.": feature \'".$feature."\' of model \'".$model."\' is \'".$model_feature{$model}{$feature}."\'...\n" ) if($self->{'VERBOSE'});
		
	return $model_feature{$model}{$feature};
}

=head2 get_result

 Usage   : $result = $obj->get_result($tree_name, $model, $result_name);
 Function: Gets the value of a given result for a given model in a given tree.
 Returns : $result (string)
 Args    : $tree_name (string), $model (string), $result_name (string)

=cut

sub get_result
{
	my($self, $tree_name, $model, $result_name) = @_;
	my $sub = $self->get_nickname().substr((caller(0))[3], MODULE_NAME_LENGTH);	
	
	if( !defined($self->{'LOADED-WORKSPACE-DIRECTORY'}) || !defined($self->{'TASK-SET'}) )
	{Carp::croak($sub." failed: no Codeml workspace loaded, exiting;");}
	
	if(!defined($tree_name))
	{Carp::croak($sub." failed: no tree name given, exiting;");}
	
	if(!defined($model))
	{Carp::croak($sub." failed: no model given, exiting;");}
	
	if(!defined($result_name))
	{Carp::croak($sub." failed: no result heading given, exiting;");}
	
	if( defined($tree_name) && !defined($self->{'RESULTS'}{$tree_name}) )
	{Carp::croak($sub." failed: no results for tree \'".$tree_name."\', exiting;");}
	
	if( defined($model) && !defined($self->{'RESULTS'}{$tree_name}{$model}) )
	{Carp::croak($sub." failed: no results for model \'".$model."\' in tree \'".$tree_name."\', exiting;");}
	
	if( !(grep { $result_name eq $_ } @{$report_headings{'summary-report'}}) ) 
	{Carp::croak($sub." failed: result name \'".$result_name."\' not supported, exiting;");}
	
	if( defined($model) && !defined($self->{'RESULTS'}{$tree_name}{$model}{$result_name}) )
	{Carp::croak($sub." failed: no result with name \'".$result_name."\' for model \'".$model."\' in tree \'".$tree_name."\', exiting;");}

	print($sub.": result \'".$result_name."\' for model \'".$model."\' in tree \'".$tree_name."\' is \'".$self->{'RESULTS'}{$tree_name}{$model}{$result_name}."\'...\n" ) if($self->{'VERBOSE'});
		
	return $self->{'RESULTS'}{$tree_name}{$model}{$result_name};
}
	
=head2 model_feature_supported

 Usage   : if($obj->model_feature_supported($model, $feature)) {...}
 Function: Tests if a model feature is supported.
 Returns : Boolean indicating if model feature is supported
 Args    : $model (string), $feature (string)

=cut

sub model_feature_supported
{
	my($self, $model, $feature) = @_;
	my $sub = $self->get_nickname().substr((caller(0))[3], MODULE_NAME_LENGTH);	
	
	if(!defined($model))
	{Carp::croak($sub." failed: no model given, exiting;");}
	
	if(!$self->model_supported($model))
	{Carp::croak($sub." failed: model \'".$model."\' not supported, exiting;");}
	
	if(!defined($feature))
	{Carp::croak($sub." failed: no model feature given, exiting;");}
	
	undef my $feature_supported;
	if(defined($model_feature{$model}{$feature}))
	{$feature_supported = TRUE;}
	else
	{$feature_supported = FALSE;}
	
	print($sub.": feature \'".$feature."\' supported in model \'".$model."\' - ".($feature_supported ? "TRUE" : "FALSE")."...\n" ) if($self->{'VERBOSE'} && ( !defined((caller(1))[3]) || (substr((caller(1))[3], 0, 5) ne 'Job::') ) );

	return $feature_supported;
}	

=head2 model_supported

 Usage   : if($obj->model_supported($model)) {...}
 Function: Tests if a model is supported.
 Returns : Boolean indicating if model is supported
 Args    : $model (string)

=cut

sub model_supported
{
	my($self, $model) = @_;
	my $sub = $self->get_nickname().substr((caller(0))[3], MODULE_NAME_LENGTH);	
	
	if(!defined($model))
	{Carp::croak($sub." failed: no model given, exiting;");}
	
	undef my $model_supported;
	if(defined($model_feature{$model}))
	{$model_supported = TRUE;}
	else
	{$model_supported = FALSE;}

	print($sub.": model \'".$model."\' supported - ".($model_supported ? "TRUE" : "FALSE")."...\n" ) if($self->{'VERBOSE'} && ( !defined((caller(1))[3]) || (substr((caller(1))[3], 0, 5) ne 'Job::') ) );

	return $model_supported;
}

sub DESTROY
{
	my($self) = @_;
			
	undef %$self;
}

sub get_nickname
{
	my ($self) = @_;

	#INTERNAL USE ONLY!
	if( !defined((caller(1))[3]) || (substr((caller(1))[3], 0, 5) ne 'Job::') )
	{Carp::croak((caller(0))[3]." failed: can only be called by \'CodemlWrapper::Job\' subroutines, exiting;");}
	
	my $nickname = 'Job';

	if( defined($self->{'LOADED-WORKSPACE-DIRECTORY'}) )
	{
		$nickname .= "(".$self->{'LOADED-WORKSPACE-DIRECTORY'}.")";
	}
	elsif( defined($self->{'TREES'}) )
	{
		my @non_tops = sort grep { $_ ne $self->{'TOPOLOGY-NAME'} } keys(%{$self->{'TREES'}});
		$nickname .= "(".join("/", @non_tops).")" if(@non_tops > 0);
	}

	return $nickname;
}

sub _setup_task_set
{
	my($self) = @_;	
	
	#INTERNAL USE ONLY!
	if( !defined((caller(1))[3]) || (caller(1))[3]!~/^Job::([^:]+)$/ )
	{Carp::croak((caller(0))[3]." failed: can only be called by \'CodemlWrapper::Job\' subroutines, exiting;");}
	
	if(defined($self->{'LOADED-WORKSPACE-DIRECTORY'}))
	{Carp::croak((caller(0))[3]." failed: can't change task-set of Codeml job after a Codeml workspace has been loaded, exiting;");}
	
	if(defined($self->{'TASK-SET'}))
	{Carp::croak((caller(0))[3]." failed: task-set already set up, exiting;");}
			
	if( !defined($self->{'ALIGNMENT'}) || !defined($self->{'TREES'}) || !defined($self->{'MODELS'}) || !defined($self->{'OMEGAS'}) )
	{Carp::croak((caller(0))[3]." failed: one or more of alignment, tree, models or initial w-values not set, exiting;");}
	
	foreach my $tree_name (keys(%{$self->{'TREES'}}))
	{
		foreach my $model ($self->get_models($tree_name))
		{
			foreach my $omega ($self->get_initial_omegas($model))
			{					
				$self->{'TASK-SET'}{$tree_name}{$model}{$omega}{'CONTROL'} = Control->new();
				$self->{'TASK-SET'}{$tree_name}{$model}{$omega}{'CONTROL'}->set_model($model);
				$self->{'TASK-SET'}{$tree_name}{$model}{$omega}{'CONTROL'}->set_parameter('omega', $omega);
			}
		}
	}
}

sub util_alignments_equal
{
	my(@alignments) = @_;
	my $alignment_equality = TRUE;
	my $n = @alignments;
	undef my @strings;
	undef my @refs;
	
	#INTERNAL USE ONLY!
	if( !defined((caller(1))[3]) || (caller(1))[3]!~/^Job::([^:]+)$/ )
	{Carp::croak((caller(0))[3]." failed: can only be called by \'CodemlWrapper::Job\' subroutines, exiting;");}
			
	if($n<2)
	{Carp::croak((caller(1))[3]." failed: can only compare two or more alignments, exiting;");}
		
	foreach my $alignment (@alignments)
	{			
		if(ref($alignment) eq 'Alignment')
		{
			push(@strings, $alignment->get_file_string());
		}
		elsif(-f $alignment)
		{
			push(@strings, util_read_file($alignment));
		}
		else
		{
			Carp::croak((caller(1))[3]." failed: can only compare \'CodemlWrapper::Alignment\' objects and/or alignment files, exiting;");
		}
	}
			
	my @string_mismatches = grep { $_ ne $strings[0] } @strings[1..($n-1)];
			
	if(@string_mismatches>0)
	{	
		foreach my $alignment (@alignments)
		{		
			if(ref($alignment) eq 'Alignment')
			{
				push(@refs, $alignment);
			}
			elsif(-f $alignment)
			{
				push(@refs, Alignment->new($alignment));
			}
		}
		
		my @alignment_mismatches = grep { !($_->alignment_match($refs[0])) } @refs[1..($n-1)];
		
		if(@alignment_mismatches>0)
		{
			$alignment_equality = FALSE;
		}
	}	
	
	return $alignment_equality;
}

sub util_topologies_equal
{
	my(@trees) = @_;
	my $tree_equality = TRUE;
	my $n = @trees;
	undef my @strings;
	undef my @refs;
	
	#INTERNAL USE ONLY!
	if( !defined((caller(1))[3]) || (caller(1))[3]!~/^Job::([^:]+)$/ )
	{Carp::croak((caller(0))[3]." failed: can only be called by \'CodemlWrapper::Job\' subroutines, exiting;");}
		
	if($n<2)
	{Carp::croak((caller(1))[3]." failed: can only compare two or more trees, exiting;");}
		
	foreach my $tree (@trees)
	{			
		if(ref($tree) eq 'Tree')
		{
			push(@strings, $tree->get_file_string());
		}
		elsif(-f $tree)
		{
			push(@strings, util_read_file($tree));
		}
		else
		{
			Carp::croak((caller(1))[3]." failed: can only compare \'CodemlWrapper::Tree\' objects and/or tree files, exiting;");
		}
	}
			
	my @string_mismatches = grep { $_ ne $strings[0] } @strings[1..($n-1)];
			
	if(@string_mismatches>0)
	{	
		foreach my $tree (@trees)
		{		
			if(ref($tree) eq 'Tree')
			{
				push(@refs, $tree);
			}
			elsif(-f $tree)
			{
				push(@refs, Tree->new($tree));
			}
		}
		
		my @tree_mismatches = grep { !($_->tree_match($refs[0])) } @refs[1..($n-1)];
		
		if(@tree_mismatches>0)
		{
			$tree_equality = FALSE;
		}
	}	
	
	return $tree_equality;
}
			
sub util_write_file
{
	my($file, @string_data) = @_;
	
	#INTERNAL USE ONLY!
	if( !defined((caller(1))[3]) || (caller(1))[3]!~/^Job::([^:]+)$/ )
	{Carp::croak((caller(0))[3]." failed: can only be called by \'CodemlWrapper::Job\' subroutines, exiting;");}
	
	my $string = "";
	
	if(@string_data==0)
	{Carp::croak((caller(1))[3]." failed: can't write empty array to \'".$file."\', exiting;");}
	
	if((@string_data==1) && ($string_data[0] eq ""))
	{Carp::croak((caller(1))[3]." failed: can't write empty string to \'".$file."\', exiting;");}
	
	if(@string_data==1)
	{
		$string = $string_data[0];
	}
	else
	{
		$string = join("\n", @string_data);
	}
	
	open(FILE_HANDLE, '>', $file) or Carp::croak((caller(1))[3]." failed: can't open \'".$file."\' for writing, exiting;");
	print FILE_HANDLE $string;
	close(FILE_HANDLE) or Carp::croak((caller(1))[3]." failed: can't close \'".$file."\' after writing, exiting;");
}

sub util_read_file
{
	my($file) = @_;
	
	#INTERNAL USE ONLY!
	if( !defined((caller(1))[3]) || (caller(1))[3]!~/^Job::([^:]+)$/ )
	{Carp::croak((caller(0))[3]." failed: can only be called by \'CodemlWrapper::Job\' subroutines, exiting;");}
		
	open(FILE_HANDLE, '<', $file) or Carp::croak((caller(1))[3]." failed: can't open \'".$file."\' for reading, exiting;");
	my @string_data = <FILE_HANDLE>;
	my $string = join("", @string_data);
	chomp(@string_data);
	close(FILE_HANDLE) or Carp::croak((caller(1))[3]." failed: can't close \'".$file."\' after reading, exiting;");
	
	if(wantarray)
	{
		return @string_data;
	}
	else
	{
		return $string;
	}
}

sub util_check_pairwise_parameter_constraints
{
	my($ref) = @_;
	
	foreach my $p (keys(%{$ref}))
	{
		my ($min, $max) = ($ref->{$p}{'min'}, $ref->{$p}{'max'});
		
		if($max < 1)
		{Carp::croak((caller(1))[3]." failed: \'max\' constraint for  \'".$p."\' shouldn't be less than one, exiting;");}

		if($min < 0)
		{Carp::croak((caller(1))[3]." failed: \'min\' constraint for  \'".$p."\' shouldn't be less than zero, exiting;");}		
	
		if($max < $min)
		{Carp::croak((caller(1))[3]." failed: \'max\' constraint for  \'".$p."\' shouldn't be less than \'min\' constraint, exiting;");}
	}
}

}

1
