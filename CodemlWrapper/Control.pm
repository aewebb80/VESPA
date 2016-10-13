#!/usr/bin/perl -w
use strict;
use Carp qw(carp cluck croak confess);
#CodemlWrapper::Control v1.0 by Tom Walsh (2012/03/28)

=head1 NAME

CodemlWrapper::Control - Wrapper for a Codeml control file.

=head1 SYNOPSIS

	use CodemlWrapper::Control;
		
	#Create new CodemlWrapper::Control object.
	my $ctl = Control->new();	
	
	#Set the model to 'modelA'.
	$ctl->set_model("modelA") if $ctl->model_supported("modelA");
	
	#Set the initial ω-value to 2.
	$ctl->set_parameter("omega", 2) if $ctl->parameter_supported("omega");
	
	#Set the nucleotide sequence alignment file name to 'align.phy'.
	$ctl->set_parameter("seqfile", 	"align.phy");
	
	#Set the tree file name to 'tree'.
	$ctl->set_parameter("treefile", "tree");
	
	#Set the output file name to 'out'.
	$ctl->set_parameter("outfile", 	"out");
	
	#Save the Codeml control file to the current directory. 
	$ctl->save();	

=head1 DESCRIPTION

Codeml is a component of Ziheng Yang's PAML package. Codeml performs selection 
analysis on gene sequences, allowing episodes of positive selection of a gene on 
evolutionary timescales to be inferred. For information on PAML see the PAML manual.

CodemlWrapper::Control is a Perl wrapper module for the control file of Codeml (PAML4.4e). 
It helps with processing of Codeml control files and is used mainly in the context 
of the CodemlWrapper::Job module to handle the Codeml control file for each task in a 
Codeml job. A control file is written to disk by calling L</save>, and can be read 
from disk using L</load>. 

A Codeml control file is essentially just a list of parameter settings that determine 
the way in which Codeml is run, including the names of input and output files, 
model types and other parameters. CodemlWrapper::Control allows the user to set the values
of the control file in two ways: by setting parameters individually or by specifying 
a preset model. 

The individual parameters that are supported by the CodemlWrapper::Control module are 
listed below, with their default values in parentheses. For more info on these 
see the PAML manual. A parameter can be set by calling L</set_parameter>, and its 
value can be obtained by calling L</get_parameter>. Calling L</reset_parameters> 
resets each parameter to its default value. 

=over 4

=item • seqfile     ('align.phy')

=item • treefile    ('tree')

=item •	outfile     ('out')

=item • noisy       (3)

=item • verbose     (0)

=item • runmode     (0)

=item • seqtype     (1)

=item • CodonFreq   (2)

=item • aaDist      (0)

=item • aaRatefile  ('wag.dat')

=item • model       (0)

=item • NSsites     (0)

=item • icode       (0)

=item • fix_kappa   (0)

=item • kappa       (3)

=item • fix_omega   (0)

=item • omega       (1)

=item • fix_alpha   (1)

=item • alpha       (0)

=item • Malpha      (0)

=item • ncatG       (1)

=item • clock       (0)

=item • getSE       (0)

=item • RateAncestor (0)

=item • Small_Diff  (.5e-6)

=back

Codeml control files can also be set to specify a preset model (listed below). 
These can be set by calling L</set_model>, and the model specified by the control 
file can be obtained using L</get_model>.

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

CodemlWrapper::Control usually runs silently. To see full output, include 'Verbose' as 
an argument when creating a L</new> object.

=head1 OBJECT METHODS

=cut

package Control;
{  	
	
use constant TRUE  => 1;
use constant FALSE => 0;
use constant MODULE_NAME_LENGTH => 7;

my %model_settings = 
(
	'm0' 			=> { 'model'=>0, 'NSsites'=>0, 'fix_omega'=>0, 'ncatG'=>1 	},	
	'2ratios'		=> { 'model'=>2, 'NSsites'=>0, 'fix_omega'=>0, 'ncatG'=>2 	},
	'm1Neutral' 	=> { 'model'=>0, 'NSsites'=>1, 'fix_omega'=>0, 'ncatG'=>2 	},
	'm2Selection' 	=> { 'model'=>0, 'NSsites'=>2, 'fix_omega'=>0, 'ncatG'=>3 	},
	'm3Discrtk2' 	=> { 'model'=>0, 'NSsites'=>3, 'fix_omega'=>0, 'ncatG'=>2 	},
	'm3Discrtk3' 	=> { 'model'=>0, 'NSsites'=>3, 'fix_omega'=>0, 'ncatG'=>3 	},
	'm7' 			=> { 'model'=>0, 'NSsites'=>7, 'fix_omega'=>0, 'ncatG'=>10 	},
	'm8'			=> { 'model'=>0, 'NSsites'=>8, 'fix_omega'=>0, 'ncatG'=>10 	},
	'm8a' 			=> { 'model'=>0, 'NSsites'=>8, 'fix_omega'=>1, 'ncatG'=>10,	'omega'=>1	},
	'modelA'		=> { 'model'=>2, 'NSsites'=>2, 'fix_omega'=>0, 'ncatG'=>4	},
	'modelAnull'	=> { 'model'=>2, 'NSsites'=>2, 'fix_omega'=>1, 'ncatG'=>4,	'omega'=>1	},
	'modelB'		=> { 'model'=>2, 'NSsites'=>3, 'fix_omega'=>0, 'ncatG'=>4	}
);
	
my @required_parameters = 
( 
	'seqfile', 
	'treefile', 
	'outfile',
	'seqtype',
	'model',
	'NSsites'
);

my @parameter_categories = 
( 
	[ my @FILE 			= ( 'seqfile', 'treefile', 'outfile' ) ],
	[ my @MODES 		= ( 'noisy', 'verbose', 'runmode' ) ],
	[ my @SEQ 			= ( 'seqtype', 'CodonFreq', 'aaDist', 'aaRatefile' ) ],
	[ my @MODEL_TYPE	= ( 'model' ) ],
	[ my @MODEL 		= ( 'NSsites' ) ],
	[ my @CODE 			= ( 'icode' ) ],
	[ my @RATIOS 		= ( 'fix_kappa', 'kappa', 'fix_omega', 'omega' ) ],
	[ my @ALPHA 		= ( 'fix_alpha', 'alpha', 'Malpha', 'ncatG' ) ],
	[ my @CLOCK 		= ( 'clock', 'getSE', 'RateAncestor' ) ],
	[ my @MISC 			= ( 'Small_Diff' ) ]
);
	
my %parameter_defaults = 
( 
	'seqfile' 		=> 'align.phy', 
	'treefile' 		=> 'tree', 
	'outfile' 		=> 'out', 
	'noisy' 		=> '3',
	'verbose' 		=> '0',
	'runmode' 		=> '0',
	'seqtype' 		=> '1',
	'CodonFreq' 	=> '2',
	'aaDist' 		=> '0',
	'aaRatefile' 	=> 'wag.dat',
	'model' 		=> '0',
	'NSsites' 		=> '0',
	'icode' 		=> '0',
	'fix_kappa' 	=> '0',
	'kappa' 		=> '3',
	'fix_omega' 	=> '0',
	'omega' 		=> '1',
	'fix_alpha' 	=> '1',
	'alpha' 		=> '0',
	'Malpha' 		=> '0',
	'ncatG' 		=> '1',
	'clock' 		=> '0',
	'getSE' 		=> '0',
	'RateAncestor' 	=> '0',
	'Small_Diff' 	=> '.5e-6'		
);

my %parameter_comments = 
( 
	'seqfile' 		=> [('* sequence data file name')], 
	'treefile' 		=> [('* tree structure file name')], 
	'outfile' 		=> [('* main result file name')], 
	'noisy' 		=> [('* 0,1,2,3,9: how much rubbish on the screen')],
	'verbose' 		=> [('* 0: concise; 1: detailed, 2: too much')],
	'runmode' 		=> [('* 0: user tree;  1: semi-automatic;  2: automatic', '* 3: StepwiseAddition; (4,5):PerturbationNNI; -2: pairwise')],
	'seqtype' 		=> [('* 1:codons; 2:AAs; 3:codons-->AAs')],
	'CodonFreq' 	=> [('* 0:1/61 each, 1:F1X4, 2:F3X4, 3:codon table')],
	'aaDist' 		=> [('* 0:equal, +:geometric; -:linear, 1-6:G1974,Miyata,c,p,v,a')],
	'aaRatefile' 	=> [('* only used for aa seqs with model=empirical(_F)', '* dayhoff.dat, jones.dat, wag.dat, mtmam.dat, or your own')],
	'model' 		=> [('* models for codons:', '* 0:one, 1:b, 2:2 or more dN/dS ratios for branches', '* models for AAs or codon-translated AAs:', '* 0:poisson, 1:proportional, 2:Empirical, 3:Empirical+F', '* 6:FromCodon, 7:AAClasses, 8:REVaa_0, 9:REVaa(nr=189)')],
	'NSsites'	 	=> [('* 0:one w;1:neutral;2:selection; 3:discrete;4:freqs;', '* 5:gamma;6:2gamma;7:beta;8:beta&w;9:beta&gamma;', '* 10:beta&gamma+1; 11:beta&normal>1; 12:0&2normal>1;', '* 13:3normal>0')],
	'icode' 		=> [('* 0:universal code; 1:mammalian mt; 2-10:see PAML manual')],
	'fix_kappa' 	=> [('* 1: kappa fixed, 0: kappa to be estimated')],
	'kappa' 		=> [('* initial or fixed kappa')],
	'fix_omega' 	=> [('* 1: omega or omega_1 fixed, 0: estimate')],
	'omega' 		=> [('* initial or fixed omega, for codons or codon-based AAs')],
	'fix_alpha' 	=> [('* 0: estimate gamma shape parameter; 1: fix it at alpha')],
	'alpha' 		=> [('* initial or fixed alpha, 0:infinity (constant rate)')],
	'Malpha' 		=> [('* different alphas for genes')],
	'ncatG' 		=> [('* # of categories in dG of NSsites models')],
	'clock' 		=> [('* 0:no clock, 1:global clock; 2:local clock; 3:TipDate')],
	'getSE' 		=> [('* 0: don\'t want them, 1: want S.E.s of estimates')],
	'RateAncestor' 	=> [('* (0,1,2): rates (alpha>0) or ancestral states (1 or 2)')],		
	'Small_Diff' 	=> [('* used in the difference approximation of derivatives')]
);

=head2 new

 Usage   : $obj = CodemlWrapper::Control->new();
           $obj = CodemlWrapper::Control->new($file);
 Function: Creates a new CodemlWrapper::Control object. If a file name is specified,  
           loads CodemlWrapper::Control object from that file.
 Returns : $obj (CodemlWrapper::Control object)
 Args    : file (optional file name)

=cut

sub new
{
	my( $self, @arguments) = @_;
	undef my %Control;
	$self = \%Control;
	
	bless( $self, 'Control');		
	
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
	elsif(@standard_arguments == 0)
	{
		$self->reset_parameters();
	}
	elsif(@standard_arguments > 1)
	{
		Carp::croak($sub." failed: invalid arguments: (\'".join("\', \'", @standard_arguments)."\'), exiting;");
	}

	return $self;
}

=head2 load

 Usage   : $obj->load();
           $obj->load($filepath);
 Function: Loads CodemlWrapper::Control object from a specified file named 'codeml.ctl'
           or from a specified directory containing a file so named. If no file 
           path is specified, the CodemlWrapper::Control object is loaded from a file 
           named 'codeml.ctl' in the current working directory.
 Returns : none
 Args    : $filepath (optional file name OR directory)

=cut

sub load
{
	my($self, $filepath) = @_;		
	my $sub = $self->get_nickname().substr((caller(0))[3], MODULE_NAME_LENGTH);

	if(defined($filepath))
	{
		if($filepath!~/codeml\.ctl$/)
		{	
			$filepath=~s/(\/|)$/\/codeml\.ctl/;
		}			
	}
	else
	{
		$filepath = "codeml.ctl";
	}
	
	if(!(-f $filepath))
	{Carp::croak($sub." failed: file \'".$filepath."\' not found, exiting;");}
	
	undef my %ctl;
	
	my @records = util_read_file($filepath);

	my $l = 1;
	foreach my $record (@records)
	{
		$record=~s/(\*.*| |\t)//g;
		
		if($record=~/^([^=]+)=([^=]+)$/)
		{						
			my ($p, $v) = ($1, $2);
			
			if(defined $ctl{$p})
			{
				Carp::croak($sub." failed: parameter \'".$p."\' in line ".$l." of \'".$filepath."\' previously set in line ".$ctl{$p}{"L"}.", exiting;");
			}
			else
			{
				#Trim leading and trailing whitespace.
				$v=~s/(^\s+|\s+$)//g;
				$ctl{$p}{"V"} = $v;
				$ctl{$p}{"L"} = $l;					
			}
		}
		elsif($record ne "")
		{
			Carp::croak($sub." failed: no valid parameter setting in line ".$l." of \'".$filepath."\', exiting;");
		}
	
		$l++;
	}
					
	foreach my $req (@required_parameters)
	{
		if(!defined($ctl{$req}{"V"}))
		{
			Carp::croak($sub." failed: required parameter \'".$req."\' is absent from \'".$filepath."\', exiting;");
		}
	}
	
	foreach my $p (keys(%ctl))
	{
		if($self->parameter_supported($p)==0)
		{
			Carp::croak($sub." failed: unsupported parameter \'".$p."\' in line ".$ctl{$p}{"L"}." of \'".$filepath."\', exiting;");
		}
		
		if(util_parameter_sanity_check($p, $ctl{$p}{"V"})==0)
		{
			Carp::croak($sub." failed: parameter \'".$p."\' has invalid value \'".$ctl{$p}{"V"}."\' in line ".$ctl{$p}{"L"}." of \'".$filepath."\', exiting;");
		}
	}
					
	foreach my $k (keys(%ctl))
	{$self->{$k} = $ctl{$k}{"V"};}
	
	$self->{'CODEML-WRAPPER-MODEL'} = $self->get_model();
	$self->{'CODEML-WRAPPER-OMEGA'} = $self->get_parameter('omega');
	
	print($sub.": Codeml control file \'".$filepath."\' loaded (model=\'".$self->{'CODEML-WRAPPER-MODEL'}."\', ω=".$self->{'CODEML-WRAPPER-OMEGA'}.")...\n" ) if($self->{'VERBOSE'});	
}

=head2 save

 Usage   : $obj->save($filepath);
 Function: Saves CodemlWrapper::Control object to a specified file named 'codeml.ctl'
           or to a file named 'codeml.ctl' in a specified directory. If no file 
           path is specified, the CodemlWrapper::Control object is saved to a file 
           named 'codeml.ctl' in the current working directory.
 Returns : none
 Args    : $filepath (optional file name OR directory)

=cut

sub save
{
	my($self, $filepath) = @_;
	my $sub = $self->get_nickname().substr((caller(0))[3], MODULE_NAME_LENGTH);
		
	undef my $ctl_string;
	
	if(defined($filepath))
	{
		if($filepath!~/codeml\.ctl$/)
		{	
			$filepath=~s/(\/|)$/\/codeml\.ctl/;
		}			
	}
	else
	{
		$filepath = "codeml.ctl";
	}
	
	my $pn_width = 0;
	my $pv_width = 0;		
	foreach my $k (keys %{$self})
	{
		my $pn_length = length($k);
		my $pv_length = length($self->{$k});
		
		if($pn_length>$pn_width)
		{$pn_width = $pn_length;}
		
		if($pv_length>$pv_width)
		{$pv_width = $pv_length;}
	}
	
	my $pn_space = " " x ($pn_width + 2);
	my $pv_space = " " x ($pv_width + 2);
	my $pc_space = " " x ($pn_width + $pv_width + 3);
	
	foreach my $category (@parameter_categories)
	{
		foreach my $parameter (@{$category})
		{
			if(defined($self->{$parameter}))
			{
				my $pn = $parameter;
				my $pv = $self->{$parameter};
									
				my $pn_blank = substr($pn_space, length($pn));
				my $pv_blank = substr($pv_space, length($pv));
				my $pc_blank = $pc_space;					
				
				my $entry    = $pn_blank.$pn." = ".$pv.$pv_blank;
				my $comment  = join( "\n\t".$pc_blank, @{$parameter_comments{$pn}});
				
				$ctl_string .= $entry.$comment."\n";
			}
		}
		
		$ctl_string .= "\n";
	}
	
	util_write_file($filepath, $ctl_string);
	
	print($sub.": Codeml control file \'".$filepath."\' saved (model=\'".$self->{'CODEML-WRAPPER-MODEL'}."\', ω=".$self->{'CODEML-WRAPPER-OMEGA'}.")...\n" ) if($self->{'VERBOSE'});	
}

=head2 reset_parameters

 Usage   : $obj->reset_parameters();
 Function: Sets each parameter to its default value.
 Returns : none
 Args    : none

=cut

sub reset_parameters
{
	my ($self) = @_;
	my $sub = $self->get_nickname().substr((caller(0))[3], MODULE_NAME_LENGTH);
	
	foreach my $p (keys(%parameter_defaults))
	{$self->{$p} = $parameter_defaults{$p};}	
	
	$self->{'CODEML-WRAPPER-MODEL'} = 'm0';
	$self->{'CODEML-WRAPPER-OMEGA'} = 1;
	
	print($sub.": parameters reset (model=\'m0\', ω=1)...\n" ) if($self->{'VERBOSE'});
}

=head2 set_model

 Usage   : $obj->set_model($model);
 Function: Sets the CodemlWrapper control file to specify the given model.
 Returns : none
 Args    : $model (string)

=cut

sub set_model
{
	my($self, $model) = @_;
	my $sub = $self->get_nickname().substr((caller(0))[3], MODULE_NAME_LENGTH);
	
	if(!defined($model))
	{Carp::croak($sub." failed: no model given, exiting;");}
			
	if($self->model_supported($model)==FALSE)
	{Carp::croak($sub." failed: unsupported model \'".$model."\', exiting;");}
	
	$self->{"model"}	 = $model_settings{$model}->{"model"};
	$self->{"NSsites"}   = $model_settings{$model}->{"NSsites"};
	$self->{"ncatG"} 	 = $model_settings{$model}->{"ncatG"};
	$self->{"fix_omega"} = $model_settings{$model}->{"fix_omega"};
	$self->{"omega"} 	 = $model_settings{$model}->{"omega"} if(defined($model_settings{$model}->{"omega"}));

	$self->{'CODEML-WRAPPER-MODEL'} = $model;
	print($sub.": model set to \'".$model."\'...\n" ) if($self->{'VERBOSE'} && ( !defined((caller(1))[3]) || (substr((caller(1))[3], 0, 9) ne 'Control::') ) );
}

=head2 get_model

 Usage   : $model = $obj->get_model();
 Function: Gets the model specified by the CodemlWrapper control file, 
           if that model is supported.
 Returns : $model (string)
 Args    : none

=cut

sub get_model
{
	my($self) = @_;
	my $sub = $self->get_nickname().substr((caller(0))[3], MODULE_NAME_LENGTH);
	
	undef my $ctl_model;
	
	my $model     = $self->{"model"};
	my $NSsites   = $self->{"NSsites"};
	my $ncatG     = $self->{"ncatG"};
	my $fix_omega = $self->{"fix_omega"};
	my $omega     = $self->{"omega"};
	
	if( !defined($model) || !defined($NSsites) || !defined($ncatG) || !defined($fix_omega) || !defined($omega) )
	{Carp::croak($sub." failed: critical parameters undefined, exiting;");}
	
	if($model==0)
	{
		if($NSsites==0)
		{
			$ctl_model = 'm0';
		}
		elsif($NSsites==1)
		{
			$ctl_model = 'm1Neutral';
		}
		elsif($NSsites==2)
		{
			$ctl_model = 'm2Selection';
		}
		elsif($NSsites==3)
		{
			if($ncatG==2)
			{$ctl_model = 'm3Discrtk2';}
			elsif($ncatG==3)
			{$ctl_model = 'm3Discrtk3';}
		}
		elsif($NSsites==7)
		{
			$ctl_model = 'm7';
		}
		elsif($NSsites==8)
		{
			if( $fix_omega==1 && $omega==1 )
			{$ctl_model = 'm8a';}
			elsif($fix_omega==0)
			{$ctl_model = 'm8';}
		}
	}
	elsif($model==2)
	{
		if($NSsites==0)
		{
			$ctl_model = '2ratios';
		}
		elsif($NSsites==2)
		{
			if( $fix_omega==1 && $omega==1 )
			{$ctl_model = 'modelAnull';}
			elsif($fix_omega==0)
			{$ctl_model = 'modelA';}
		}
		elsif($NSsites==3)
		{
			$ctl_model = 'modelB';
		}
	}
	
	if(!defined($ctl_model))
	{Carp::croak($sub." failed: model unsupported, exiting;");}
	
	print($sub.": model is \'".$ctl_model."\'...\n" ) if($self->{'VERBOSE'} && ( !defined((caller(1))[3]) || (substr((caller(1))[3], 0, 9) ne 'Control::') ) );
	
	return $ctl_model;
}

=head2 set_parameter

 Usage   : $obj->set_parameter($parameter, $value);
 Function: Sets the value of the given parameter.
 Returns : none
 Args    : $parameter (string), $value (string OR number)

=cut

sub set_parameter
{
	my($self, $parameter, $value) = @_;
	my $sub = $self->get_nickname().substr((caller(0))[3], MODULE_NAME_LENGTH);
	
	if(!defined($parameter))
	{Carp::croak($sub." failed: no Codeml control parameter given, exiting;");}
	
	if($self->parameter_supported($parameter)==FALSE)
	{Carp::croak($sub." failed: unsupported parameter \'".$parameter."\', exiting;");}
	
	if(!defined($value))
	{Carp::croak($sub." failed: no value given for Codeml control parameter \'".$parameter."\', exiting;");}
		
	if(util_parameter_sanity_check($parameter, $value)==0)
	{Carp::croak($sub." failed: invalid value (\'".$value."\') for Codeml control parameter \'".$parameter."\', exiting;");}
	
	$self->{$parameter} = $value;	
	
	$self->{'CODEML-WRAPPER-OMEGA'} = $value if($parameter eq 'omega');
	
	print($sub.": parameter \'".$parameter."\' set to \'".$value."\'...\n" ) if($self->{'VERBOSE'} && ( !defined((caller(1))[3]) || (substr((caller(1))[3], 0, 9) ne 'Control::') ) );
}

=head2 get_parameter

 Usage   : $value = $obj->get_parameter($parameter);
 Function: Gets the value of the given parameter.
 Returns : $value (string OR number)
 Args    : $parameter (string)

=cut

sub get_parameter
{
	my($self, $parameter) = @_;
	my $sub = $self->get_nickname().substr((caller(0))[3], MODULE_NAME_LENGTH);

	if(!defined($parameter))
	{Carp::croak($sub." failed: no Codeml control parameter given, exiting;");}
	
	if($self->parameter_supported($parameter)==FALSE)
	{Carp::croak($sub." failed: unsupported parameter \'".$parameter."\', exiting;");}
	
	if(!defined($self->{$parameter}))
	{Carp::croak($sub." failed: Codeml control parameter \'".$parameter."\' not defined, exiting;");}
	
	print($sub.": parameter \'".$parameter."\' has value \'".$self->{$parameter}."\'...\n" ) if($self->{'VERBOSE'} && ( !defined((caller(1))[3]) || (substr((caller(1))[3], 0, 9) ne 'Control::') ) );
	
	return $self->{$parameter};
}

=head2 model_supported

 Usage   : if($obj->model_supported($model)) {...}
 Function: Tests if the given model is supported.
 Returns : Boolean indicating if the given model is supported
 Args    : $model (string)

=cut

sub model_supported
{
	my($self, $model) = @_;
	my $sub = $self->get_nickname().substr((caller(0))[3], MODULE_NAME_LENGTH);
	
	if(!defined($model))
	{Carp::croak($sub." failed: no model given, exiting;");}
			
	undef my $model_supported;
	if(defined($model_settings{$model}))
	{$model_supported = TRUE;}
	else
	{$model_supported = FALSE;}
	
	print($sub.": model \'".$model."\' supported - ".($model_supported ? "TRUE" : "FALSE")."...\n" ) if($self->{'VERBOSE'} && ( !defined((caller(1))[3]) || (substr((caller(1))[3], 0, 9) ne 'Control::') ) );
	
	return $model_supported;
}

=head2 parameter_supported

 Usage   : if($obj->parameter_supported($parameter)) {...}
 Function: Tests if the given Codeml control parameter is supported.
 Returns : Boolean indicating if the given parameter is supported
 Args    : $parameter (string)

=cut

sub parameter_supported
{
	my( $self, $parameter ) = @_;	
	my $sub = $self->get_nickname().substr((caller(0))[3], MODULE_NAME_LENGTH);
	
	if(!defined($parameter))
	{Carp::croak($sub." failed: no Codeml control parameter given, exiting;");}
	
	undef my $parameter_supported;
	if(defined($parameter_defaults{$parameter}))
	{$parameter_supported = TRUE;}
	else
	{$parameter_supported = FALSE;}

	print($sub.": parameter \'".$parameter."\' supported - ".($parameter_supported ? "TRUE" : "FALSE")."...\n" ) if($self->{'VERBOSE'} && ( !defined((caller(1))[3]) || (substr((caller(1))[3], 0, 9) ne 'Control::') ) );
	
	return $parameter_supported;
}

sub DESTROY
{
	my ($self) = @_;
			
	undef %$self;
}

sub get_nickname
{
	my ($self) = @_;

	#INTERNAL USE ONLY!
	if( !defined((caller(1))[3]) || (substr((caller(1))[3], 0, 9) ne 'Control::') )
	{Carp::croak((caller(0))[3]." failed: can only be called by \'CodemlWrapper::Control\' subroutines, exiting;");}
	
	my $nickname = 'Control';

	if( defined($self->{'CODEML-WRAPPER-MODEL'}) && defined($self->{'CODEML-WRAPPER-OMEGA'}) )
	{
		$nickname .= "(".$self->{'CODEML-WRAPPER-MODEL'}.",ω=".$self->{'CODEML-WRAPPER-OMEGA'}.")";
	}
	
	return $nickname;
}

sub util_parameter_sanity_check
{
	my($parameter, $value) = @_;
			
	#INTERNAL USE ONLY!
	if( !defined((caller(1))[3]) || (caller(1))[3]!~/^Control::([^:]+)$/ )
	{Carp::croak((caller(0))[3]." failed: can only be called by \'CodemlWrapper::Control\' subroutines, exiting;");}
	
	#NB*: the filenames "seqfile", "treefile" and "outfile" automatically pass the sanity check.
	
	if( ($parameter eq "noisy") && ($value!~/^[01239]$/) )										#noisy
	{return FALSE;}
	
	if( ($parameter eq "runmode") && ($value!~/^(0|1|2|3|4|5|-2)$/) )							#runmode
	{return FALSE;}
	
	if( ($parameter eq "seqtype") && ($value!=1) )												#seqtype
	{return FALSE;}
	
	if( ($parameter eq "aaDist") && ($value!~/^(0|\+|-|1|2|3|4|5|6)$/) )						#aaDist
	{return FALSE;}
	
	if( ($parameter eq "aaRatefile") && ($value!~/\.dat$/) )									#aaRatefile
	{return FALSE;}
	
	if( ($parameter eq "model") && ($value!~/^[01236789]$/) )									#model
	{return FALSE;}
	
	if( ($parameter eq "NSsites") && ($value!~/^(0|1|2|3|4|5|6|7|8|9|10|11|12|13)$/) )			#NSsites
	{return FALSE;}
	
	if( ($parameter eq "icode") && ($value!~/^(0|1|2|3|4|5|6|7|8|9|10)$/) )						#icode
	{return FALSE;}
	
	if( ($parameter eq "RateAncestor") && ($value!~/^[012]$/) )									#RateAncestor
	{return FALSE;}
	
	if( ($parameter eq "verbose") && ($value!~/^[012]$/) )										#verbose
	{return FALSE;}	
	
	if( ($parameter eq "CodonFreq") && ($value!~/^[0123]$/) )									#CodonFreq
	{return FALSE;}
	
	if( ($parameter eq "clock") && ($value!~/^[0123]$/) )										#clock
	{return FALSE;}
	
	if( ($parameter eq "getSE") && ($value!~/^[01]$/) )											#getSE
	{return FALSE;}
	
	if( ($parameter eq "fix_omega") && ($value!~/^[01]$/) )										#fix_omega
	{return FALSE;}	

	if( ($parameter eq "fix_kappa") && ($value!~/^[01]$/) )										#fix_kappa
	{return FALSE;}	
	
	if( ($parameter eq "fix_alpha") && ($value!~/^[01]$/) )										#fix_alpha
	{return FALSE;}	
		
	if( ($parameter eq "Malpha") && ($value!~/^[01]$/) )										#Malpha
	{return FALSE;}	

	if(  ($parameter eq "Small_Diff") && ( ($value!~/^([+]?)(?=\d|\.\d)\d*(\.\d)?([Ee]([+-]?\d+))?$/) || ($value>1.0) )  )	#Small_Diff
	{return FALSE;}	
	
	if( ($parameter eq "kappa") && ($value!~/^(\+?)(?=\d|\.\d)\d*(\.\d)?([Ee]([+-]?\d+))?$/) )	#kappa
	{return FALSE;}		

	if( ($parameter eq "alpha") && ($value!~/^(\+?)(?=\d|\.\d)\d*(\.\d)?([Ee]([+-]?\d+))?$/) )	#alpha
	{return FALSE;}	

	if( ($parameter eq "omega") && ($value!~/^(\+?)(?=\d|\.\d)\d*(\.\d)?([Ee]([+-]?\d+))?$/) )	#omega
	{return FALSE;}	
	
	if(  ($parameter eq "ncatG") && ( ($value!~/^(\d+)$/) || ($value==0) )  )					#ncatG
	{return FALSE;}		
	
	return TRUE;
}

sub util_write_file
{
	my($file, @string_data) = @_;
	
	#INTERNAL USE ONLY!
	if( !defined((caller(1))[3]) || (caller(1))[3]!~/^Control::([^:]+)$/ )
	{Carp::croak((caller(0))[3]." failed: can only be called by \'CodemlWrapper::Control\' subroutines, exiting;");}
	
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
	if( !defined((caller(1))[3]) || (caller(1))[3]!~/^Control::([^:]+)$/ )
	{Carp::croak((caller(0))[3]." failed: can only be called by \'CodemlWrapper::Control\' subroutines, exiting;");}
		
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

}

1
