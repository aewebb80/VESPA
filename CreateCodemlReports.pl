#!/usr/bin/env perl
use strict;
use warnings;
use CodemlWrapper::Job;
#CreateCodemlReports v1.0 by Tom Walsh (2012/03/28)

=head1 NAME

CreateCodemlReports.pl - Creates a set of reports for a Codeml job. 

=head1 USAGE

perl CreateCodemlReports.pl codeml_workspace report_directory [-alignment=original_alignment] [-tree=tree_of_interest] [-sequence=sequence_of_interest] [-models=model_set] [-omegas=omega_set] [-threshold=pos_site_threshold] [-overwrite=overwrite_option]

=head1 DESCRIPTION

CreateCodemlReports.pl searches a Codeml workspace and processes any Codeml output 
found. All reports are saved within the specified report directory.

In the report directory, Codeml summary reports are produced. If trees are specified, 
summary reports are created for those trees. Otherwise, reports are created for all 
trees. An LRT report is also created in each report directory, showing the intermediate
values of the LRT calculations. In addition, a pos site report and pos site MSA are 
produced for every model with inferred positive sites. If one or more sequences are 
specified, pos site reports are created for each sequence in turn. An original 
alignment can be specified if it's desired to map pos sites back to that original 
alignment. If an original alignment is specified, any pos site MSAs will be based  
on the original alignment.

If the specified report directory already exists and the 'overwrite' option isn't 
set, the user is prompted to confirm that (s)he wishes to write reports to that 
directory. Only set an existing directory as the report directory if you know what 
you're doing. 

=head1 ARGUMENTS

=over 4

=item • codeml_workspace

Codeml workspace: a directory structure containing a set of Codeml tasks 
- and in this case, output - associated with a specific Codeml job. 

=item • report_directory

The directory within which Codeml reports are to be saved. 

=back

=head1 OPTIONS

=over 4

=item • -alignment=original_alignment

The name of a masked alignment file in 'HTML' format, in which sites that were 
removed before Codeml analysis are instead masked. This allows sites to be mapped 
to their original positions. Note that this will only work with masked alignments 
created by the 'CodemlWrapper::Alignment' module.

=item • -tree=tree_of_interest

The name of a tree of interest. If a tree other than the topology tree is specified, 
homogeneous and site-specific models are reported along with lineage-specific models 
for that tree. If the topology tree is specified, homogeneous and site-specific models 
are reported. If no tree is specified, reports are produced for all trees in the 
Codeml job. 

=item • -sequence=sequence_of_interest

The name of a sequence of interest, or the string 'FOREGROUND'. If 'FOREGROUND' 
is specified, positive sites are reported for every sequence in the foreground. 
If particular sequences are given, positive sites are reported for those sequences. 
If no sequence is given, positive sites are reported as positions in the alignment, 
without amino acid symbols.

=item • -models=model_set

A model to include in the Codeml reports, specified individually on the command 
line or as a list in a text file with extension '.txt'. If no models are specified,
the default models are used.

=item • -omegas=omega_set

An initial omega to include in the Codeml reports, specified individually on 
the command line or as a list in a text file with extension '.txt'. If no omegas
are specified, the default omegas are used.

=item • -threshold=pos_site_threshold

Pos site probability threshold. Pos sites with a probability lower than the 
threshold will be ignored. If no threshold is specified, the default (0.5) will 
be used.

=item • -overwrite=overwrite_option

The overwrite option indicates whether an existing directory should be used as a
report directory. If 'YES', existing directories will be used as a report directory.
If 'NO', specifying an existing directory as a report directory will cause the 
script to fail.

=back

=head1 DEPENDENCIES

=over 4

=item • L<CodemlWrapper::Job>

=back

=cut 

#Accept arguments.
my (@arguments) = @ARGV;

#Init variables.
undef my @standard_arguments;
undef my $alignment;
undef my @sequences;
undef my @trees;
undef my @models;
undef my @omegas;
undef my $threshold;
undef my $overwrite;

#Process each argument.
foreach my $arg (@arguments)
{	
	#If marked as an original alignment, set original alignment..
	if($arg=~/^-alignment=(.+)/i)
	{
		my $a = $1;
		
		if( defined($alignment) )
		{die($0." failed: multiple original alignments specified, exiting;");}
		
		$alignment = $a;
	}	
	#..otherwise, if marked as a tree name or tree list file, add to list of trees of interest..
	elsif($arg=~/^-tree=(.+)/i)
	{
		my $t = $1;
		
		my @list = ($t);
		if($t=~/\.txt$/)
		{
			open(LIST_FILE, '<', $t) or die($0." failed: can't open \'".$t."\' for reading, exiting;");
			@list = <LIST_FILE>;
			chomp(@list);
			close(LIST_FILE) or die($0." failed: can't close \'".$t."\' after reading, exiting;");
		}
		
		foreach my $item (@list)
		{
			#Trim leading and trailing whitespace.
			$item=~s/(^\s*|\s*$)//g;		
		
			if( !(grep {$item eq $_} @trees) )
			{push(@trees, $item);}
		}
	}
	#..otherwise, if marked as a sequence name or sequence list file, add to list of sequences of interest..
	elsif($arg=~/^-sequence=(.+)/i)
	{
		my $s = $1;
		
		my @list = ($s);
		if($s=~/\.txt$/)
		{
			open(LIST_FILE, '<', $s) or die($0." failed: can't open \'".$s."\' for reading, exiting;");
			@list = <LIST_FILE>;
			chomp(@list);
			close(LIST_FILE) or die($0." failed: can't close \'".$s."\' after reading, exiting;");
		}
		
		foreach my $item (@list)
		{
			#Trim leading and trailing whitespace.
			$item=~s/(^\s*|\s*$)//g;		
		
			if( !(grep {$item eq $_} @sequences) )
			{push(@sequences, $item);}
		}
	}		
	#..otherwise, if marked as a model or model list file, add to list of models..
	elsif($arg=~/^-models=(.+)/i)
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
	#..otherwise, if marked as a pos site probability threshold, set pos site probability threshold..
	elsif($arg=~/^-threshold=(.+)/i)
	{
		my $t = $1;
		
		if( defined($threshold) )
		{die($0." failed: multiple pos site probability thresholds specified, exiting;");}
		
		$threshold = $t;
	}
	#..otherwise, if marked as an overwrite option, set overwrite option..
	elsif($arg=~/^-overwrite=(.+)/i)
	{
		my $o = $1;
		my $u = uc($o);
		
		if( defined($overwrite) )
		{die($0." failed: multiple overwrite options specified, exiting;");}
		
		if( ($u ne 'YES') && ($u ne 'NO') )
		{die($0." failed: invalid overwrite option (\'".$o."\'), exiting;");}		
		
		$overwrite = $u;
	}		
	#..otherwise, it's a standard argument.
	else
	{
		push(@standard_arguments, $arg);
	}
}

#Verify that there are two standard arguments.
if(@standard_arguments != 2)
{die($0." failed: must have 2 standard arguments, exiting;");}

#Set variables from arguments.
my ($codeml_workspace, $report_directory) = @standard_arguments;

#If Codeml workspace directory can't be found, die.
if( !(-d $codeml_workspace) )
{die($0." failed: couldn't find Codeml workspace \'".$codeml_workspace."\', exiting;");}

#If the report directory already exists, decide whether to write into that directory..
if(-d $report_directory)
{
	#If the overwrite option has been set, decide based on that..
	if(defined($overwrite))
	{
		if($overwrite eq 'NO')
		{
			print("ERROR: Can't write to existing directory \'".$report_directory."\', exiting;\n");
			exit(1);
		}
	}
	#..otherwise, prompt the user to decide.
	else
	{
		undef my $decision;
		
		#Prompt the user to confirm whether existing directory is to be used as report directory.
		print("Report directory \'".$report_directory."\' already exists.\n");
		print("Are you sure you want to create reports in \'".$report_directory."\'? (y|n): ");
		
		#Accept input.
		$decision = <STDIN>;
		
		#While user input is invalid, prompt the user again.
		while($decision!~/^(y|n)$/i)
		{	
			#Prompt user.
			print("Please type \'y\' if you wish to set \'".$report_directory."\' as report directory, or \'n\' to cancel and exit.\n");
			print("Create Codeml reports in \'".$report_directory."\'? (y|n): ");		
			
			#Accept input.
			$decision = <STDIN>;
		}
		
		#If user decision is 'n', exit.	
		if($decision=~/^n$/i)
		{
			print($0." terminated by user.\n");
			exit(0);
		}
	}
}
#..otherwise, make report directory.
else
{
	mkdir $report_directory or die($0." failed: couldn't make directory \'".$report_directory."\', exiting;");
}

print "Creating Codeml reports in \'".$report_directory."\' using Codeml output from \'".$codeml_workspace."\'...\n";

#Create new Codeml::Job. 
my $Codeml_Job = Job->new();
	
#Set models if specified.
$Codeml_Job->set_models(@models) if(@models > 0);

#Set initial omegas if specified.
$Codeml_Job->set_initial_omegas(@omegas) if(@omegas > 0);

#Set pos site probability threshold if specified.
$Codeml_Job->set_positive_site_probability_threshold($threshold) if(defined($threshold));

#Load from Codeml workspace.
$Codeml_Job->load($codeml_workspace);

#Get topology tree name.
my $topology = $Codeml_Job->get_topology_name();	

#Create summary report for Codeml workspace.
my $summary_report_file = $report_directory."/SummaryReport.txt";
my @summary_report_args = ('OUTFILE' => $summary_report_file);
push(@summary_report_args, ('ORIGINAL-ALIGNMENT' => $alignment)) if(defined($alignment));
push(@summary_report_args, ('SEQUENCE' => $_)) foreach (@sequences);
push(@summary_report_args, ('TREE' => $_)) foreach (@trees);
$Codeml_Job->create_summary_report(@summary_report_args);

print("Creating Codeml summary report \'".$summary_report_file."\'...\n");
print("Summary report pos sites mapped to original alignment: \'".$alignment."\'...\n") if(defined($alignment));	

#Set trees if not specified.
@trees = $Codeml_Job->get_tree_names() if(@trees==0);

#Process output for each tree in the Codeml job.
foreach my $tree (@trees)
{
	#Get the sequences of interest for this tree, including foreground if specified.
	my @SOIs = map { if($_ eq 'FOREGROUND'){ $Codeml_Job->get_foreground($tree) }else{ $_ } } @sequences;

	#Process output for each model in the tree.
	foreach my $model ($Codeml_Job->get_models($tree))
	{
		#Get positive selection status of the model.
		my $positive_selection_status = $Codeml_Job->get_result($tree, $model, 'Positive Selection');
		
		#Get positive site capability of the model.
		my $positive_site_capability = $Codeml_Job->get_model_feature($model, 'positive-site-capability');
		
		#If positive selection is inferred by the model and positive sites are 
		#allowed, create a positive site report and pos site MSA for this model.
		if( ($positive_selection_status eq 'Yes') && ($positive_site_capability eq 'Allowed') )
		{
			#If no sequences specified, create pos site report for given model in tree wrt alignment..
			if(@SOIs == 0)
			{
				#Create pos site report for given model in tree wrt alignment.
				my $pos_sites_report_file = $report_directory."/PosSites_".(($tree eq $topology) ? $model : $tree."_".$model).".txt";
				my @pos_sites_args = ('OUTFILE' => $pos_sites_report_file, 'TREE' => $tree, 'MODEL' => $model);
				push(@pos_sites_args, ('ORIGINAL-ALIGNMENT' => $alignment)) if(defined($alignment));
				$Codeml_Job->create_pos_site_report(@pos_sites_args);	
				
				print("Creating Codeml pos site report \'".$pos_sites_report_file."\'...\n");
				print("Pos sites mapped to original alignment: \'".$alignment."\'...\n") if(defined($alignment));
			}
			#..otherwise, create pos site reports for given model in tree wrt specified sequences.
			else
			{
				foreach my $SOI (@SOIs)
				{
					#Create pos site report for given model in tree wrt sequence of interest.
					my $pos_sites_report_file = $report_directory."/PosSites_".(($tree eq $topology) ? $model : $tree."_".$model)."_".$SOI.".txt";
					my @pos_sites_args = ('OUTFILE' => $pos_sites_report_file, 'TREE' => $tree, 'MODEL' => $model, 'SEQUENCE' => $SOI);
					push(@pos_sites_args, ('ORIGINAL-ALIGNMENT' => $alignment)) if(defined($alignment));
					$Codeml_Job->create_pos_site_report(@pos_sites_args);
					
					print("Creating Codeml pos site report \'".$pos_sites_report_file."\'...\n");
					print("Pos sites mapped to original alignment: \'".$alignment."\'...\n") if(defined($alignment));					
				}
			}
		
			#Create pos site MSA for given model in tree.
			my $pos_sites_msa_file = $report_directory."/PosSites_".(($tree eq $topology) ? $model : $tree."_".$model).".html";
			my @pos_msa_args = ('OUTFILE' => $pos_sites_msa_file, 'TREE' => $tree, 'MODEL' => $model);
			push(@pos_msa_args, ('ORIGINAL-ALIGNMENT' => $alignment)) if(defined($alignment));
			$Codeml_Job->create_pos_site_msa(@pos_msa_args);
			
			print("Creating Codeml pos site MSA \'".$pos_sites_msa_file."\'...\n");
			print("Pos site MSA based on original alignment: \'".$alignment."\'...\n") if(defined($alignment));	
		}
	}
}

#Create LRT report for Codeml Job.
my $lrt_report_file = $report_directory."/LRTs.txt";
$Codeml_Job->create_lrt_report('OUTFILE' => $lrt_report_file);
print "Creating LRT report \'".$lrt_report_file."\'...\n\n";


