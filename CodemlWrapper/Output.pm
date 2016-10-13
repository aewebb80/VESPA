#!/usr/bin/perl -w
use strict;
use Carp qw(carp cluck croak confess);
#CodemlWrapper::Output v1.0 by Tom Walsh (2012/03/28)

=head1 NAME

CodemlWrapper::Output - Wrapper for a Codeml output file.

=head1 SYNOPSIS

	use CodemlWrapper::Output;
	
	#Set up variables.
	my $null_rejected = FALSE;
	undef my @pos_sites_m8;
	
	#Create new CodemlWrapper::Output objects for 'm8' and its null model, 'm8a'.
	my $output_m8  = Output->new("Workspace1/Eutheria/m8/Omega0");
	my $output_m8a = Output->new("Workspace1/Eutheria/m8a/Omega1");
	
	#Get the log-likelihood for each model.
	my $lnL_m8  = $output_m8->get_lnL();
	my $lnL_m8a = $output_m8a->get_lnL();
	
	#Do a likelihood ratio test; if the likelihood of m8 is significantly 
	#greater than that of m8a, m8a is rejected.
	if( 2*($lnL_m8-$lnL_m8a) > 2.71)
	{
	   $null_rejected = TRUE;
	}
	
	#If the null model (m8a) is rejected, get the 
	#positively selected sites for m8. 
	if($null_rejected)
	{
	   @pos_sites_m8 = $output_m8->get_positive_sites();
	}

=head1 DESCRIPTION

Codeml is a component of Ziheng Yang's PAML package. Codeml performs selection 
analysis on gene sequences, allowing episodes of positive selection of a gene on 
evolutionary timescales to be inferred. For information on PAML see the PAML manual.

CodemlWrapper::Output is a Perl wrapper module for the standard output file of Codeml 
(PAML4.4e). It helps with processing of Codeml output files and is used mainly 
in the context of the CodemlWrapper::Job module to handle the Codeml output file for 
each completed task in a Codeml job. Codeml output files are written by PAML and 
are treated as read-only by the CodemlWrapper::Output module. They can be read from disk 
using L</load>. 

CodemlWrapper::Output currently doesn't process all parts of the Codeml output file, 
just those that are most relevant to interpreting the results of a Codeml task. 
These can be obtained using a number of subroutines (listed below), each related to a particular 
relevant feature of the output. 

CodemlWrapper::Output usually runs silently. To see full output, include 
'Verbose' as an argument when creating a L</new> object.

=head1 OBJECT METHODS

=cut

package Output;
{

use constant TRUE  => 1;
use constant FALSE => 0;
use constant MODULE_NAME_LENGTH => 6;

=head2 new

 Usage   : $obj = CodemlWrapper::Output->new();
           $obj = CodemlWrapper::Output->new($file);
           $obj = CodemlWrapper::Output->new('Verbose');
           $obj = CodemlWrapper::Output->new($file, 'Verbose');
 Function: Creates a new CodemlWrapper::Output object. If a file name is specified,  
           loads CodemlWrapper::Output object from that file. If 'Verbose' is 
           specified, the CodemlWrapper::Output object will have verbose output.
 Returns : $obj (CodemlWrapper::Output object)
 Args    : $file (optional file name)

=cut

sub new
{
	my($self, @arguments) = @_;
	undef my %Output;
	$self = \%Output;
	
	bless($self, 'Output');
	
	my $sub = $self->get_nickname().substr((caller(0))[3], MODULE_NAME_LENGTH);
	
	#Process arguments.
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
		my ($file) = @standard_arguments;
		
		$self->load($file);
	}
	elsif(@standard_arguments > 1)
	{
		Carp::croak($sub." failed: invalid arguments: (\'".join("\', \'", @standard_arguments)."\'), exiting;");
	}
	
	return $self;
}

=head2 load

 Usage   : $obj->load($file);
 Function: Loads CodemlWrapper::Output object from the specified file.
 Returns : none
 Args    : $file (file name)

=cut

sub load
{
	my($self, $file) = @_;
	my $sub = $self->get_nickname().substr((caller(0))[3], MODULE_NAME_LENGTH);
	
	if(!defined($file))
	{Carp::croak($sub." failed: no file name given, exiting;");}
	
	if(!(-f $file))
	{Carp::croak($sub." failed: file \'".$file."\' not found, exiting;");}
	
	delete $self->{'OUTPUT'};
	$self->{'OUTPUT'}{'FILE'} = $file;
	
	undef my $top;
	undef my $header_string;
	undef my $model_cue;
	undef my $omega_fixed;
	undef my $model;
	
	undef my $bottom;
	undef my $time_used;
	
	my $file_string = util_read_file($file);
	my $total = length($file_string);

	if($file_string=~/^((.|\n)*)(CODONML((.|\n)+))(ns =(\s+)(\d+)  ls =(\s+)(\d+)\n\n)/)
	{
		$top = length($&);
		$header_string = $3;			
		$self->{'OUTPUT'}{'SEQ-COUNT'}  = $8;
		$self->{'OUTPUT'}{'SEQ-LENGTH'} = $10;
	}
	else
	{
		Carp::croak($sub." failed: no Codeml output file header in \'".$file."\', exiting;");
	}
	
	if($header_string=~/^(\S+) \(in (.+)\)\s+(\S.+)\nModel: (.+) for branches/)
	{
		$self->{'OUTPUT'}{'PROGRAM'}  = $1;
		$self->{'OUTPUT'}{'VERSION'}  = $2;
		$self->{'OUTPUT'}{'SEQ-FILE'} = $3;			
		$model_cue = $4;
	}
	else
	{
		Carp::croak($sub." failed: invalid Codeml output file header in \'".$file."\', exiting;");
	}
	
	if($header_string=~/omega =\s+(\S+)\s+fixed/)
	{$omega_fixed = $1;}
	
	if($header_string=~/Site-class models:\s+([^\n]+)/)
	{
		my $siteclass_model = $1;
		
		if($model_cue eq "One dN/dS ratio")
		{
			if($siteclass_model 	eq "NearlyNeutral")
			{$model = "m1Neutral";}
			elsif($siteclass_model 	eq "PositiveSelection")
			{$model = "m2Selection";}
			elsif($siteclass_model 	eq "discrete \(2 categories\)")
			{$model = "m3Discrtk2";}
			elsif($siteclass_model 	eq "discrete \(3 categories\)")
			{$model = "m3Discrtk3";}	
			elsif($siteclass_model 	eq "beta \(10 categories\)")
			{$model = "m7";}
			elsif(($siteclass_model eq "beta&w>1 \(11 categories\)") && (!defined($omega_fixed)))
			{$model = "m8"}
			elsif(($siteclass_model eq "beta&w>1 \(11 categories\)") && (defined($omega_fixed)) && ($omega_fixed==1.000))
			{$model = "m8a"}
		}
		elsif($model_cue eq "several dN/dS ratios for branches")
		{
			if(($siteclass_model 	eq "PositiveSelection") && (!defined($omega_fixed)))
			{$model = "modelA";}
			elsif(($siteclass_model eq "PositiveSelection") && (defined($omega_fixed)) && ($omega_fixed==1.000))
			{$model = "modelAnull";}
			elsif($siteclass_model 	eq "discrete (4 categories)")
			{$model = "modelB";}
		}
	}
	else
	{
		if($model_cue 	 eq "One dN/dS ratio")
		{
			$model = "m0";
		}
		elsif($model_cue eq "several dN/dS ratios for branches")
		{
			$model = "2ratios";	
		}
	}
	
	if(defined($model))
	{$self->{'OUTPUT'}{'MODEL'} = $model;}
	else
	{Carp::croak($sub." failed: unsupported model in Codeml output file \'".$file."\', exiting;");}
	
	if($file_string=~/\nTime used:(\s+)(\S+)\n$/)
	{
		$bottom = length($&);
		
		my $raw_time = $2;
		my $time = 0;
		
		if($raw_time=~/((\d+):)?(\d+):(\d+)/)
		{
			my $hours = defined($2) ? $2 : 0 ;
			my $minutes = $3;
			my $seconds = $4;
			
			$time = ($hours*3600) + ($minutes*60) + $seconds;
		}
		
		$self->{'OUTPUT'}{'TIME'}  = $time;
	}
	else
	{
		Carp::croak($sub." failed: incomplete output file \'".$file."\', exiting;");
	}	
	
	my $body_text = substr($file_string, $top, $total-($top+$bottom));
	
	undef my $data_text;
	if($body_text=~/\n\n(TREE # (\s+\d+):  )/)
	{$data_text = $1.$';}
	else
	{Carp::croak($sub." failed: incomplete output file \'".$file."\', exiting;");}
			
	if($data_text=~/\nlnL\(ntime:\s*\d+  np:\s*\d+\):\s+(\S+)\s+\S+\n/)
	{$self->{'OUTPUT'}{'lnL'} = $1;}
	else
	{Carp::croak($sub." failed: can't find \'lnL\' value in Codeml output file \'".$file."\', exiting;");}
	
	if($data_text=~/kappa \(ts\/tv\) =(\s*)(\d+\.\d+)/)
	{$self->{'OUTPUT'}{'kappa'}=$';}
	else
	{Carp::croak($sub." failed: can't find \'kappa\' value in Codeml output file \'".$file."\', exiting;");}
	
	my $parameters_not_found = FALSE;
	if($model eq 'm0')
	{
		if($data_text=~/omega \(dN\/dS\) =  (\d+\.\d+)/)
		{
			$self->{'OUTPUT'}{'PARAMETERS'}{'w'}=$1;
		}
		else
		{$parameters_not_found = TRUE;}
	}
	elsif($model eq '2ratios')
	{	
		if($data_text=~/w \(dN\/dS\) for branches:\s*(\d+\.\d{5})\s*(\d+\.\d{5})/)
		{
			$self->{'OUTPUT'}{'PARAMETERS'}{'w0'}=$1;
			$self->{'OUTPUT'}{'PARAMETERS'}{'w1'}=$2;
		}
		else
		{$parameters_not_found = TRUE;}			
	}
	elsif( ($model eq 'm1Neutral') || ($model eq 'm2Selection') || ($model eq 'm3Discrtk2')|| ($model eq 'm3Discrtk3') )
	{
		if($data_text=~/dN\/dS( \(w\))? for site classes \(K=(\d+)\)\n\np:\s+\b(.+)\nw:\s+\b(.+)\n/)
		{
			my( $class_tally, $prop_string, $omega_string ) = ($2, $3, $4);
			my @prop_array = ($prop_string=~/\d{1,3}\.\d{5}/g);
			my @omega_array = ($omega_string=~/\d{1,3}\.\d{5}/g);
			
			if( (@prop_array!=$class_tally) || (@omega_array!=$class_tally) )
			{Carp::croak($sub." failed: can't process detailed parameters in Codeml output file \'".$file."\', exiting;");}
			
			for(my $i=0 ; $i<@prop_array ; $i++)
			{
				$self->{'OUTPUT'}{'PARAMETERS'}{'p'.$i} = $prop_array[$i];
				$self->{'OUTPUT'}{'PARAMETERS'}{'w'.$i} = $omega_array[$i];
			}
		}
		else
		{$parameters_not_found = TRUE;}
	}
	elsif( ($model eq 'm7') || ($model eq 'm8') || ($model eq 'm8a') )
	{
		if($data_text=~/Parameters in M(7|8).+:\n((.+\n)+)\n/)
		{
			my $param_string = $2;
			$param_string=~s/[\s\n\(\)]//g;
			my @param_array = ($param_string=~/\w+=\d+\.\d+/g);
			%{$self->{'OUTPUT'}{'PARAMETERS'}} = map { if($_=~/=/){ $` => $' }else{} } @param_array ;
		}
		else
		{$parameters_not_found = TRUE;}
	}
	elsif( ($model eq 'modelA') || ($model eq 'modelAnull') || ($model eq 'modelB') )
	{
		if($data_text=~/dN\/dS( \(w\))? for site classes \(K=(\d+)\)\n\nsite class\s+\b(.+)\nproportion\s+\b(.+)\nbackground w\s+\b(.+)\nforeground w\s+\b(.+)\n/)
		{
			my( $class_tally, $classes, $prop_string, $background_string, $foreground_string ) = ($2, $3, $4, $5, $6);
			my @prop_array = ($prop_string=~/\d{1,3}\.\d{5}/g);
			my @background_array = ($background_string=~/\d{1,3}\.\d{5}/g);
			my @foreground_array = ($foreground_string=~/\d{1,3}\.\d{5}/g);
			
			if( (@prop_array!=$class_tally) || (@background_array!=$class_tally) || (@foreground_array!=$class_tally) )
			{Carp::croak($sub." failed: can't process detailed parameters in Codeml output file \'".$file."\', exiting;");}
			
			my $expected_props = 4;
			for(my $i=0 ; $i<$expected_props ; $i++)
			{$self->{'OUTPUT'}{'PARAMETERS'}{'p'.$i} = $prop_array[$i];}
		
			my $expected_omegas = 3;
			for(my $i=0 ; $i<$expected_omegas ; $i++)
			{$self->{'OUTPUT'}{'PARAMETERS'}{'w'.$i} = $foreground_array[$i];}
		}
		else
		{$parameters_not_found = TRUE;}
	}
	
	if($parameters_not_found)
	{Carp::croak($sub." failed: can't find detailed parameters in Codeml output file \'".$file."\', exiting;");}
	
	undef my $pos_site_string;
	if( $data_text =~ /Bayes Empirical Bayes (.+\n)Positively selected sites(.*\n){5}((.+\n)+)/ )
	{ 
		$pos_site_string = $3;
		$self->{'OUTPUT'}{'SITES-INFERENCE-TYPE'} = 'BEB';
	}
	elsif( $data_text =~ /Bayes Empirical Bayes (.+\n)Positive sites for foreground lineages(.+\n)((.+\n)+)/ )
	{
		$pos_site_string = $3;
		$self->{'OUTPUT'}{'SITES-INFERENCE-TYPE'} = 'BEB';
	}
	elsif( $data_text =~ /Naive Empirical Bayes (.+\n)Positively selected sites(.*\n){5}((.+\n)+)/ )
	{
		$pos_site_string = $3;
		$self->{'OUTPUT'}{'SITES-INFERENCE-TYPE'} = 'NEB';
	}
	elsif( $data_text =~ /Naive Empirical Bayes (.+\n)Positive sites for foreground lineages(.+\n)\n((.+\n)+)/ )
	{
		$pos_site_string = $3;
		$self->{'OUTPUT'}{'SITES-INFERENCE-TYPE'} = 'NEB';
	}
	
	if(defined($pos_site_string))
	{
		while($pos_site_string=~/\s*(\d+)\s(\S)\s+(\S+)(\s+(\d+\.\d+( \+- \d+\.\d+)?))?/g)
		{
			my ( $position, $amino_acid, $probability, $posterior_mean ) = ( $1, $2, $3, $5 );
			
			$probability=~s/\*{1,2}//;	#remove asterisks from probability
			
			$self->{'OUTPUT'}{'SITES'}{$position}{'Amino Acid'}	 = $amino_acid;
			$self->{'OUTPUT'}{'SITES'}{$position}{'P(w>1)'} 	 = $probability;
			$self->{'OUTPUT'}{'SITES'}{$position}{'Posterior Mean (w)'} = $posterior_mean if(defined($posterior_mean));
		}
	}
	
	print($sub.": Codeml output loaded from with \'".$file."\'...\n" ) if($self->{'VERBOSE'});	
}

=head2 get_lnL

 Usage   : $lnL = $obj->get_lnL();
 Function: Gets the log-likelihood estimate of the model used.
 Returns : $lnL (number)
 Args    : none

=cut

sub get_lnL
{
	my($self) = @_;
	my $sub = $self->get_nickname().substr((caller(0))[3], MODULE_NAME_LENGTH);
	
	if(!defined($self->{'OUTPUT'}))
	{Carp::croak($sub." failed: no Codeml output loaded, exiting;");}
	
	print($sub.": \'".$self->{'OUTPUT'}{'lnL'}."\'...\n" ) if($self->{'VERBOSE'});	
	
	return $self->{'OUTPUT'}{'lnL'};	
}

=head2 get_model_parameters

 Usage   : @parameters = $obj->get_model_parameters();
 Function: Gets a list of the parameters of the model used.
 Returns : @parameters (array of strings)
 Args    : none

=cut

sub get_model_parameters
{
	my($self) = @_;
	my $sub = $self->get_nickname().substr((caller(0))[3], MODULE_NAME_LENGTH);
	
	if(!defined($self->{'OUTPUT'}))
	{Carp::croak($sub." failed: no Codeml output loaded, exiting;");}	
	
	my @parameter_list = sort( keys(%{$self->{'OUTPUT'}{'PARAMETERS'}}) );
	
	print($sub.": ".@parameter_list." parameters retrieved...\n" ) if($self->{'VERBOSE'});	
	
	return @parameter_list;	
}

=head2 get_parameter_value

 Usage   : $parameter_value = $obj->get_parameter_value($parameter);
 Function: Gets the value of the given parameter.
 Returns : $parameter_value (number)
 Args    : $parameter (string)

=cut

sub get_parameter_value
{
	my($self, $parameter) = @_;
	my $sub = $self->get_nickname().substr((caller(0))[3], MODULE_NAME_LENGTH);
	
	if(!defined($self->{'OUTPUT'}))
	{Carp::croak($sub." failed: no Codeml output loaded, exiting;");}
	
	if(!defined($self->{'OUTPUT'}{'PARAMETERS'}{$parameter}))
	{Carp::croak($sub." failed: parameter \'".$parameter."\' not defined, exiting;");}	
	
	print($sub.": parameter \'".$parameter."\' has value \'".$self->{'OUTPUT'}{'PARAMETERS'}{$parameter}."\'...\n" ) if($self->{'VERBOSE'} && ( !defined((caller(1))[3]) || (substr((caller(1))[3], 0, 8) ne 'Output::') ) );	
	
	return $self->{'OUTPUT'}{'PARAMETERS'}{$parameter};		
}

=head2 get_positive_sites

 Usage   : @positive_site_positions = $obj->get_positive_sites();
 Function: Gets a list of the positions of inferred positive sites. When output 
           from Codeml, positive site positions are given in terms of the alignment.
 Returns : @positive_site_positions (array of numbers)
 Args    : none

=cut

sub get_positive_sites
{
	my($self) = @_;
	my $sub = $self->get_nickname().substr((caller(0))[3], MODULE_NAME_LENGTH);
	
	if(!defined($self->{'OUTPUT'}))
	{Carp::croak($sub." failed: no Codeml output loaded, exiting;");}
	
	if(!defined($self->{'OUTPUT'}{'SITES'}))
	{Carp::croak($sub." failed: no positive site data, exiting;");}	
	
	my @pos_site_positions = sort {$a <=> $b} ( keys(%{$self->{'OUTPUT'}{'SITES'}}) );
	
	print($sub.": ".@pos_site_positions." positive sites retrieved...\n" ) if($self->{'VERBOSE'});
	
	return @pos_site_positions;
}	

=head2 get_pos_site_inference_type

 Usage   : $pos_site_inference_type = $obj->get_pos_site_inference_type();
 Function: Gets the positive site inference type. Possible values are 'NEB'
           (Naive Empirical Bayes) and 'BEB' (Bayes Empirical Bayes). 
 Returns : $pos_site_inference_type (string)
 Args    : none

=cut

sub get_pos_site_inference_type
{
	my($self) = @_;
	my $sub = $self->get_nickname().substr((caller(0))[3], MODULE_NAME_LENGTH);
	
	if(!defined($self->{'OUTPUT'}))
	{Carp::croak($sub." failed: no Codeml output loaded, exiting;");}

	if(!defined($self->{'OUTPUT'}{'SITES-INFERENCE-TYPE'}))
	{Carp::croak($sub." failed: positive site inference type not defined, exiting;");}		
	
	print($sub.": \'".$self->{'OUTPUT'}{'SITES-INFERENCE-TYPE'}."\'...\n" ) if($self->{'VERBOSE'});
	
	return $self->{'OUTPUT'}{'SITES-INFERENCE-TYPE'};
}

=head2 get_pos_site_attribute

 Usage   : $pos_site_attribute = $obj->get_pos_site_attribute($position, $attribute_name);
 Function: Gets the value of the specified attribute for the given positive site.
 Returns : $pos_site_attribute (string)
 Args    : $position (number), $attribute_name (string)

=cut

sub get_pos_site_attribute
{
	my($self, $position, $attribute_name) = @_;
	my $sub = $self->get_nickname().substr((caller(0))[3], MODULE_NAME_LENGTH);
	
	if(!defined($self->{'OUTPUT'}))
	{Carp::croak($sub." failed: no Codeml output loaded, exiting;");}
	
	if(!defined($self->{'OUTPUT'}{'SITES'}))
	{Carp::croak($sub." failed: no positive site data, exiting;");}	
	
	if(!defined($position))
	{Carp::croak($sub." failed: no positive site position given, exiting;");}
	
	if(!defined($self->{'OUTPUT'}{'SITES'}{$position}))
	{Carp::croak($sub." failed: no positive site at position ".$position.", exiting;");}	
	
	if(!defined($attribute_name))
	{Carp::croak($sub." failed: no positive site attribute name given, exiting;");}
	
	if(! grep { $_ eq $attribute_name } ('Amino Acid', 'P(w>1)', 'Posterior Mean (w)') )
	{Carp::croak($sub." failed: \'".$attribute_name."\' is not a valid positive site attribute, exiting;");}

	my $pos_site_attribute = $self->{'OUTPUT'}{'SITES'}{$position}{$attribute_name};
	
	print($sub.": ".$attribute_name." of pos site at codon ".$position." is \'".$pos_site_attribute."\'...\n" ) if($self->{'VERBOSE'});
	
	return $pos_site_attribute;		
}

=head2 get_positive_selection_status

 Usage   : if($obj->get_positive_selection_status()) {...}
 Function: Tests if Codeml has inferred positive selection.
 Returns : Boolean indicating if Codeml has inferred positive selection
 Args    : none

=cut

sub get_positive_selection_status
{
	my($self) = @_;
	my $sub = $self->get_nickname().substr((caller(0))[3], MODULE_NAME_LENGTH);
	
	if(!defined($self->{'OUTPUT'}))
	{Carp::croak($sub." failed: no Codeml output loaded, exiting;");}	
	
	if(!defined($self->{'OUTPUT'}{'PARAMETERS'}))
	{Carp::croak($sub." failed: parameters not defined, exiting;");}	
	
	my $positive_selection = FALSE;
	my $model = $self->get_model();		
	
	if($model eq "m0")
	{
		my $w = $self->get_parameter_value("w");
		
		if($w>1.0)
		{$positive_selection = TRUE;}
	}
	elsif($model eq "2ratios")
	{
		my $w0 = $self->get_parameter_value("w0");
		my $w1 = $self->get_parameter_value("w1");
		
		if( ($w0<=1.0) && ($w1>1.0) )
		{$positive_selection = TRUE;}		
	}
	elsif($model eq "m2Selection")
	{
		my $w2 = $self->get_parameter_value("w2");
		my $p2 = $self->get_parameter_value("p2");
		
		if( ($w2>1.0) && ($p2>0.0) )
		{$positive_selection = TRUE;}
	}
	elsif( ($model eq "m3Discrtk2") || ($model eq "m3Discrtk3") )
	{
		my $x = 0;
		while($self->parameter_defined("w".$x))
		{
			my $w = $self->get_parameter_value("w".$x);
			my $p = $self->get_parameter_value("p".$x);
			
			if( ($w>1.0) && ($p>0.0) )
			{$positive_selection = TRUE;}
			
			$x++;
		}
	}
	elsif($model eq "m8")
	{
		my $w = $self->get_parameter_value("w");
		my $p1 = $self->get_parameter_value("p1");
		
		if( ($w>1.0) && ($p1>0.0) )
		{$positive_selection = TRUE;}
	}
	elsif( ($model eq "modelA") || ($model eq "modelB") )
	{
		my $w1 = $self->get_parameter_value("w1");
		my $w2 = $self->get_parameter_value("w2");
		my $p2 = $self->get_parameter_value("p2");
		my $p3 = $self->get_parameter_value("p3");
		
		if( ($w1<=1.0) && ($w2>1.0) && (($p2+$p3)>0.0) )
		{$positive_selection = TRUE;}
	}
	
	print($sub.": ".($positive_selection ? "TRUE" : "FALSE")."...\n" ) if($self->{'VERBOSE'});
		
	return $positive_selection;
}

=head2 get_kappa

 Usage   : $kappa = $obj->get_kappa();
 Function: Gets the transition/transversion rate ratio (kappa).
 Returns : $kappa (number)
 Args    : none

=cut

sub get_kappa
{
	my($self) = @_;
	my $sub = $self->get_nickname().substr((caller(0))[3], MODULE_NAME_LENGTH);
	
	if(!defined($self->{'OUTPUT'}))
	{Carp::croak($sub." failed: no Codeml output loaded, exiting;");}
	
	print($sub.": \'".$self->{'OUTPUT'}{'kappa'}."\'...\n" ) if($self->{'VERBOSE'});
	
	return $self->{'OUTPUT'}{'kappa'};	
}

=head2 get_program_info

 Usage   : $program_info = $obj->get_program_info();
 Function: Gets info about the version of PAML used for the Codeml task.
 Returns : $program_info (string)
 Args    : none

=cut

sub get_program_info
{
	my($self) = @_;
	my $sub = $self->get_nickname().substr((caller(0))[3], MODULE_NAME_LENGTH);
	
	if(!defined($self->{'OUTPUT'}))
	{Carp::croak($sub." failed: no Codeml output loaded, exiting;");}
	
	print($sub.": ".$self->{'OUTPUT'}{'PROGRAM'}." (in ".$self->{'OUTPUT'}{'VERSION'}.")...\n" ) if($self->{'VERBOSE'});
	
	return $self->{'OUTPUT'}{'PROGRAM'}." (in ".$self->{'OUTPUT'}{'VERSION'}.")";
}

=head2 get_model

 Usage   : $model = $obj->get_model();
 Function: Gets the model used in the Codeml task.
 Returns : $model (string)
 Args    : none

=cut

sub get_model
{
	my($self) = @_;
	my $sub = $self->get_nickname().substr((caller(0))[3], MODULE_NAME_LENGTH);
	
	if(!defined($self->{'OUTPUT'}))
	{Carp::croak($sub." failed: no Codeml output loaded, exiting;");}
	
	print($sub.": \'".$self->{'OUTPUT'}{'MODEL'}."\'...\n" ) if($self->{'VERBOSE'} && ( !defined((caller(1))[3]) || (substr((caller(1))[3], 0, 8) ne 'Output::') ) );
	
	return $self->{'OUTPUT'}{'MODEL'};	
}

=head2 get_sequence_count

 Usage   : $sequence_count = $obj->get_sequence_count();
 Function: Gets the number of sequences in the input Codeml alignment.
 Returns : $sequence_count (number)
 Args    : none

=cut

sub get_sequence_count
{
	my($self) = @_;
	my $sub = $self->get_nickname().substr((caller(0))[3], MODULE_NAME_LENGTH);
	
	if(!defined($self->{'OUTPUT'}))
	{Carp::croak($sub." failed: no Codeml output loaded, exiting;");}

	print($sub.": \'".$self->{'OUTPUT'}{'SEQ-COUNT'}."\'...\n" ) if($self->{'VERBOSE'});

	return $self->{'OUTPUT'}{'SEQ-COUNT'};
}

=head2 get_sequence_length

 Usage   : $sequence_length = $obj->get_sequence_length();
 Function: Gets the sequence length (in codons) of the input Codeml alignment.
 Returns : $sequence_length (number)
 Args    : none

=cut

sub get_sequence_length
{
	my($self) = @_;
	my $sub = $self->get_nickname().substr((caller(0))[3], MODULE_NAME_LENGTH);
	
	if(!defined($self->{'OUTPUT'}))
	{Carp::croak($sub." failed: no Codeml output loaded, exiting;");}

	print($sub.": \'".$self->{'OUTPUT'}{'SEQ-LENGTH'}."\'...\n" ) if($self->{'VERBOSE'});	

	return $self->{'OUTPUT'}{'SEQ-LENGTH'};
}

=head2 get_execution_time

 Usage   : $execution_time = $obj->get_execution_time();
 Function: Gets the execution time (in seconds) of the Codeml task.
 Returns : $execution_time (number)
 Args    : none

=cut

sub get_execution_time
{
	my($self) = @_;
	my $sub = $self->get_nickname().substr((caller(0))[3], MODULE_NAME_LENGTH);
	
	if(!defined($self->{'OUTPUT'}))
	{Carp::croak($sub." failed: no Codeml output loaded, exiting;");}

	print($sub.": \'".$self->{'OUTPUT'}{'TIME'}."\'...\n" ) if($self->{'VERBOSE'});	
	
	return $self->{'OUTPUT'}{'TIME'};
}

=head2 parameter_defined

 Usage   : if($obj->parameter_defined($parameter)) {...}
 Function: Tests if a parameter is defined under the model used.
 Returns : Boolean indicating if parameter is defined
 Args    : $parameter (string)

=cut

sub parameter_defined
{
	my($self, $parameter) = @_;
	my $sub = $self->get_nickname().substr((caller(0))[3], MODULE_NAME_LENGTH);
	
	if(!defined($self->{'OUTPUT'}))
	{Carp::croak($sub." failed: no Codeml output loaded, exiting;");}
	
	my $parameter_status = defined($self->{'OUTPUT'}{'PARAMETERS'}{$parameter});
	
	print($sub.": parameter \'".$parameter."\' defined - ".($parameter_status ? "TRUE" : "FALSE")."...\n" ) if($self->{'VERBOSE'} && ( !defined((caller(1))[3]) || (substr((caller(1))[3], 0, 8) ne 'Output::') ) );	
	
	return $parameter_status;
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
	if( !defined((caller(1))[3]) || (substr((caller(1))[3], 0, 8) ne 'Output::') )
	{Carp::croak((caller(0))[3]." failed: can only be called by \'CodemlWrapper::Output\' subroutines, exiting;");}
	
	my $nickname = 'Output';

	if( defined($self->{'OUTPUT'}{'FILE'}) )
	{
		$nickname .= "(".$self->{'OUTPUT'}{'FILE'}.")";
	}

	return $nickname;
}	

sub util_read_file
{
	my($file) = @_;
	
	#INTERNAL USE ONLY!
	if( !defined((caller(1))[3]) || (caller(1))[3]!~/^Output::([^:]+)$/ )
	{Carp::croak((caller(0))[3]." failed: can only be called by \'CodemlWrapper::Output\' subroutines, exiting;");}
		
	open(FILE_HANDLE, '<', $file) or Carp::croak((caller(1))[3]."\' failed: can't open \'".$file."\' for reading, exiting;");
	my @string_data = <FILE_HANDLE>;
	my $string = join("", @string_data);
	chomp(@string_data);
	close(FILE_HANDLE) or Carp::croak((caller(1))[3]."\' failed: can't close \'".$file."\' after reading, exiting;");
	
	if(wantarray)
	{
		return @string_data;
	}
	else
	{
		return $string;
	}
}

}

1
