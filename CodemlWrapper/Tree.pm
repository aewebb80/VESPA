#!/usr/bin/env perl
use strict;
use warnings;
use Carp qw(carp cluck croak confess);
#CodemlWrapper::Tree v1.01 by Tom Walsh (2012/03/28)

=head1 NAME

CodemlWrapper::Tree - Wrapper for a Codeml tree file.

=head1 SYNOPSIS

	use CodemlWrapper::Tree;
	
	#Create list of regexes.
	my @genes = ('Human', 'Chimp', 'Gorilla');
	
	#Create new CodemlWrapper::Tree object, load from file.
	my $tree = Tree->new('tree.nwk');	
		
	#Ensure tree has no branch lengths or labels.
	$tree->remove_attributes('BRANCH-LENGTHS', 'BRANCH-LABELS');
	
	#Ensure tree is unrooted.
	$tree->unroot_tree();
	
	#Label the tree foreground using gene names (in this example these are taken 
	#from species names). Each gene is set as a foreground gene and the common 
	#ancestor of each monophyletic group of foreground genes is labelled as foreground.
	$tree->set_foreground_genes(@genes);
	
	#Save the Codeml tree file to the current directory in Newick format. 
	$tree->save('tree', 'NEWICK');	

=head1 DESCRIPTION

Codeml is a component of Ziheng Yang's PAML package. Codeml performs selection 
analysis on gene sequences, allowing episodes of positive selection of a gene on 
evolutionary timescales to be inferred. For information on PAML see the PAML manual.

CodemlWrapper::Tree is a Perl wrapper module for the Codeml tree file used by Codeml 
(PAML4.4e). It helps with processing of Codeml tree files and is used mainly in 
the context of the CodemlWrapper::Job module to handle the Codeml tree file for each 
task in a Codeml job. A Codeml tree file is written to disk by calling L</save>, 
and can be read from disk using L</load>. 

Trees can be loaded in 3 similar file formats: 'PAML' tree format, 'PHYLIP' tree 
format and Newick tree format. Within the module, these are named 'PAML', 'PHYLIP' 
and 'NEWICK' respectively. The output format of a Codeml alignment is set using 
L</set_output_format>, and can be obtained using L</get_output_format>. Although 
PAML allows for more than one tree per file, for simplicity's sake CodemlWrapper::Tree 
only supports one tree per file. Newick format is the simplest of the three, and 
is the default output format. 

Codeml trees also allow labelling of particular lineages. This can be done by 
manually adding a '#1' to the branch (or branches) of interest, so that any branch
followed by '#1' is considered a 'foreground' lineage, while all other branches 
are considered 'background'. It can also be done automatically through CodemlWrapper::Tree 
by passing a list of gene names to L</set_foreground_genes>.

Most of the time Codeml analyses will be used with an unrooted tree, and in fact 
CodemlWrapper::Job unroots trees automatically. It's possible to test whether a tree is 
rooted using L</is_rooted>, and unroot the tree with L</unroot_tree>. By convention, 
a tree is unrooted if it has a multifurcation at the root. If a tree is bifurcated, 
CodemlWrapper::Tree unroots it by removing the branch whose clade has the larger number
of genes and reattaching its children to the root.

Comparison of trees can be performed using L</tree_match>. This takes as arguments 
a reference to the CodemlWrapper::Tree being compared to this one, and optional attribute
names. The two attributes that can be considered are 'BRANCH-LENGTHS' and 'BRANCH-LABELS'.

CodemlWrapper::Tree usually runs silently. To see full output, include 'Verbose' 
as an argument when creating a L</new> object.

=head1 OBJECT METHODS

=cut

package Tree;
{

use constant TRUE  => 1;
use constant FALSE => 0;
use constant MIN_GENES => 7;
use constant MODULE_NAME_LENGTH => 4;
use constant DEFAULT_OUTPUT_FORMAT => 'NEWICK';

my @supported_file_types = ('PAML',	'PHYLIP', 'NEWICK');
my @supported_attributes = ('BRANCH-LENGTHS', 'BRANCH-LABELS');

=head2 new

 Usage   : $obj = CodemlWrapper::Tree->new();
           $obj = CodemlWrapper::Tree->new($file);
           $obj = CodemlWrapper::Tree->new('Verbose');
           $obj = CodemlWrapper::Tree->new($file, 'Verbose');
 Function: Creates a new CodemlWrapper::Tree object. If a file name is given,  
           loads CodemlWrapper::Tree object from that file. If 'Verbose' is 
           specified, the CodemlWrapper::Tree object will have verbose output.
 Returns : $obj (CodemlWrapper::Tree object)
 Args    : $file (optional file name)

=cut

sub new
{
	my($self, @arguments) = @_;
	undef my %Tree;
	$self = \%Tree;
	
	bless( $self, 'Tree');
	
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

=head2 copy

 Usage   : $copy = $original->copy();
 Function: Creates a copy of a CodemlWrapper::Tree object identical to the original.
 Returns : $copy (CodemlWrapper::Tree object)
 Args    : none

=cut

sub copy
{
	my($self) = @_;
	my $sub = $self->get_nickname().substr((caller(0))[3], MODULE_NAME_LENGTH);
	
	undef my %Tree;
	my $copy = \%Tree;

	bless( $copy, 'Tree');

	$copy->{'INPUT-FORMAT'} = $self->{'INPUT-FORMAT'};
	$copy->{'OUTPUT-FORMAT'} = $self->{'OUTPUT-FORMAT'};
	$copy->{'FILE-STRING'} = $self->{'FILE-STRING'};
	
	my $tree_string = $self->get_tree_string('BRANCH-LABELS','BRANCH-LENGTHS');
	$copy->{'TREE'} = $copy->parse_tree($tree_string);

	$copy->{'FOREGROUND'} = $self->{'FOREGROUND'};
	
	print($sub.": tree copied...\n" ) if($self->{'VERBOSE'});
	
	return $copy;
}

=head2 load

 Usage   : $obj->load($file);
 Function: Loads CodemlWrapper::Tree object from given file.
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
	
	delete $self->{'INPUT-FORMAT'};
	delete $self->{'OUTPUT-FORMAT'};
	delete $self->{'FILE-STRING'};
	delete $self->{'FOREGROUND'};
	delete $self->{'TREE'};
	
	$self->{'OUTPUT-FORMAT'} = DEFAULT_OUTPUT_FORMAT;
	
	my $file_string = util_read_file($file);		
	$self->{'FILE-STRING'} = $file_string;		#THIS IS READ-ONLY FROM THIS POINT ON!!
	
	undef my $spec_leaf_tally;
	undef my $specified_tree_count;
	undef my $tree_data;
		
	if($file_string=~/^((\s*)(\d+)(\s+)(\d+)(\s*))\n((\(([^;]+);([^\(]*))+)$/) 
	{
		$self->{'INPUT-FORMAT'} = 'PAML';
		$spec_leaf_tally = $3;
		$specified_tree_count  = $5;
		$tree_data 		  = $7;			
	}
	elsif($file_string=~/^((\s*)(\d+)(\s*))\n((\(([^;]+);([^\(]*))+)$/)
	{
		$self->{'INPUT-FORMAT'} = 'PHYLIP';
		$specified_tree_count  = $3;
		$tree_data 		  = $5;			
	}
	elsif($file_string=~/^([^\w\(]*)((\(([^;]+);([^\(]*))+)$/)
	{
		$self->{'INPUT-FORMAT'} = 'NEWICK';
		$tree_data 		  = $2;		
	}
	else
	{
		Carp::croak($sub." failed: invalid Codeml tree file \'".$file."\', exiting;");			
	}		
			
	if( defined($specified_tree_count) && ($specified_tree_count != 1) )
	{Carp::croak($sub." failed: specified number of trees (".$specified_tree_count.") in \'".$file."\' is not supported (max=1), exiting;");}
	
	$tree_data=~s/\n//g;
	my @tree_array = ($tree_data=~/\([^;]+;/g);
	
	if(@tree_array==0)
	{Carp::croak($sub." failed: no valid tree found in \'".$file."\', exiting;");}
	
	if(@tree_array>1)
	{Carp::croak($sub." failed: number of trees (".@tree_array.") in \'".$file."\' is not supported (max=1)");}
	
	$self->{'TREE'} = $self->parse_tree($tree_array[0]);
		
	my $gene_count = $self->get_gene_count();
					
	if( defined($spec_leaf_tally) && ($gene_count!=$spec_leaf_tally) )
	{Carp::croak($sub." failed: specified number of genes (".$spec_leaf_tally.") in \'".$file."\' differs from apparent number (".$gene_count."), exiting;");}	

	if($gene_count < MIN_GENES)
	{Carp::croak($sub." failed: number of genes (".$gene_count.") in \'".$file."\' is below minimum number (".MIN_GENES."), exiting;");}	

	my @foreground_genes = sort {$a cmp $b} util_recurse_tree($self->{'TREE'}, \&util_node_get_foreground_genes);
	$self->{'FOREGROUND'} = [@foreground_genes] if(@foreground_genes > 0);

	print($sub.": Codeml tree with ".$gene_count." genes loaded from ".$self->{'INPUT-FORMAT'}." file \'".$file."\'...\n" ) if($self->{'VERBOSE'});		
}

=head2 save

 Usage   : $obj->save($file, $format);
 Function: Saves CodemlWrapper::Tree object to given file.  If a supported 
           format is specified, saves the file in that format.
 Returns : none
 Args    : $file (file name), $format (optional string)

=cut

sub save
{	
	my($self, $file, $output_format) = @_;
	my $sub = $self->get_nickname().substr((caller(0))[3], MODULE_NAME_LENGTH);
		
	my $output_string = '';
	
	if(!defined($file))
	{Carp::croak($sub." failed: no file name given, exiting;");}
	
	if(!defined($self->{'TREE'}))
	{Carp::croak($sub." failed: no tree defined, exiting;");}
	
	$self->set_output_format($output_format) if(defined($output_format));	
	
	my $gene_count = $self->get_gene_count();
	
	if($self->{'OUTPUT-FORMAT'} eq 'PAML')
	{
		$output_string = $gene_count." 1\n".$self->get_tree_string('BRANCH-LENGTHS', 'BRANCH-LABELS')."\n";
	}
	elsif($self->{'OUTPUT-FORMAT'} eq 'PHYLIP')
	{
		$output_string = "1\n".$self->get_tree_string('BRANCH-LENGTHS', 'BRANCH-LABELS')."\n";
	}
	elsif($self->{'OUTPUT-FORMAT'} eq 'NEWICK')
	{		
		$output_string .= $self->get_tree_string('BRANCH-LENGTHS', 'BRANCH-LABELS')."\n";
	}
	else
	{
		Carp::croak($sub." failed: no output format defined, exiting;");
	}
	
	util_write_file($file, $output_string);
	
	print($sub.": Codeml tree with ".$gene_count." genes saved to ".$self->{'OUTPUT-FORMAT'}." file \'".$file."\'...\n" ) if($self->{'VERBOSE'});	
}

=head2 get_genes

 Usage   : @genes = $obj->get_genes();
 Function: Gets the genes in the tree.
 Returns : @genes (array of strings)
 Args    : none

=cut

#This subroutine relies on the fact that the leaves are sorted when the tree is being input. 
sub get_genes
{
	my($self) = @_;
	my $sub = $self->get_nickname().substr((caller(0))[3], MODULE_NAME_LENGTH);
	
	if(!defined($self->{'TREE'}))
	{Carp::croak($sub." failed: no tree defined, exiting;");}
	
	my @leaves = @{$self->{'TREE'}{'LEAVES'}};
	
	print($sub.": ".@leaves." genes retrieved...\n" ) if($self->{'VERBOSE'});	
	
	return @leaves;
}

=head2 get_gene_count

 Usage   : $gene_count = $obj->get_gene_count();
 Function: Gets the number of genes in the tree.
 Returns : $gene_count (number)
 Args    : none

=cut

sub get_gene_count
{
	my($self) = @_;
	my $sub = $self->get_nickname().substr((caller(0))[3], MODULE_NAME_LENGTH);
	
	if(!defined($self->{'TREE'}))
	{Carp::croak($sub." failed: no tree defined, exiting;");}
	
	my $gene_count = @{$self->{'TREE'}{'LEAVES'}};

	print($sub.": ".$gene_count." genes present...\n" ) if($self->{'VERBOSE'} && ( !defined((caller(1))[3]) || (substr((caller(1))[3], 0, 6) ne 'Tree::') ) );
	
	return $gene_count;
}

=head2 set_foreground_genes

 Usage   : $obj->set_foreground_genes(@genes);
 Function: Sets foreground branches of Codeml tree such that specified genes are 
           marked as foreground.
 Returns : none
 Args    : @genes (array of strings)

=cut

sub set_foreground_genes
{
	my($self, @genes) = @_;
	my $sub = $self->get_nickname().substr((caller(0))[3], MODULE_NAME_LENGTH);
	
	if(!defined($self->{'TREE'}))
	{Carp::croak($sub." failed: no tree defined, exiting;");}
	
	if(@genes==0)
	{Carp::croak($sub." failed: no genes given, exiting;");}
				
	util_recurse_tree($self->{'TREE'}, \&util_node_remove_attributes, ['BRANCH-LABELS']);
	
	my $foreground_gene_count = 0;
	foreach my $leaf_node (@{$self->{'LEAF-NODES'}})
	{
		if(grep {$leaf_node->{'NODE-NAME'} eq $_} @genes)
		{
			$leaf_node->{'BRANCH-LABELS'} = '#1';
			$foreground_gene_count++;
		}
	}
	
	if($foreground_gene_count != @genes)
	{Carp::croak($sub." failed: not all specified genes are in tree, exiting;");}		
			
	if($foreground_gene_count == @{$self->{'TREE'}{'LEAVES'}})
	{Carp::croak($sub." failed: can't label all genes as foreground, exiting;");}
	
	util_recurse_tree($self->{'TREE'}, \&util_node_coalesce_foreground_branch_labels);
	
	$self->{'FOREGROUND'} = [@genes];
	
	print($sub.": ".$foreground_gene_count." foreground genes set...\n" ) if($self->{'VERBOSE'});
}

=head2 get_foreground_genes

 Usage   : @foreground_genes = $obj->get_foreground_genes();
 Function: Gets the foreground genes of the tree. If the tree is unlabelled, an
           empty array is returned.
 Returns : @foreground_genes (array of strings)
 Args    : none

=cut

sub get_foreground_genes
{
	my($self) = @_;
	my $sub = $self->get_nickname().substr((caller(0))[3], MODULE_NAME_LENGTH);
	
	undef my @foreground_genes;
	
	if(!defined($self->{'TREE'}))
	{Carp::croak($sub." failed: no tree defined, exiting;");}
	
	@foreground_genes = util_recurse_tree($self->{'TREE'}, \&util_node_get_foreground_genes);
	
	@foreground_genes = sort {$a cmp $b} @foreground_genes;

	print($sub.": ".@foreground_genes." foreground genes retrieved...\n" ) if($self->{'VERBOSE'});
	
	return @foreground_genes;	
}

=head2 is_rooted

 Usage   : if($obj->is_rooted()) {...}
 Function: Tests if the tree is rooted (i.e. bifurcation at root).
 Returns : Boolean indicating if the tree is rooted
 Args    : none

=cut

sub is_rooted
{
	my($self) = @_;
	my $sub = $self->get_nickname().substr((caller(0))[3], MODULE_NAME_LENGTH);
	
	if(!defined($self->{'TREE'}))
	{Carp::croak($sub." failed: no tree defined, exiting;");}
	
	undef my $status;
	
	if(@{$self->{'TREE'}{'CHILDREN'}}==2)
	{$status = TRUE;}
	else
	{$status = FALSE;}

	print($sub.": ".($status ? "TRUE" : "FALSE")."...\n" ) if($self->{'VERBOSE'});
	
	return $status;
}

=head2 unroot_tree

 Usage   : $obj->unroot_tree();
 Function: Unroots the tree by forcing a multifurcation at the root. If the tree 
           is bifurcated at the root, the branch whose clade has the larger number
           of genes is removed and its children reattached at the root.
 Returns : none
 Args    : none

=cut

sub unroot_tree
{
	my($self) = @_;
	my $sub = $self->get_nickname().substr((caller(0))[3], MODULE_NAME_LENGTH);
	
	if(!defined($self->{'TREE'}))
	{Carp::croak($sub." failed: no tree defined, exiting;");}
			
	if(@{$self->{'TREE'}{'CHILDREN'}}==2)
	{
		undef my $chosen_index;
		undef my $chosen_child;
				
		if( @{$self->{'TREE'}{'CHILDREN'}[0]{'LEAVES'}} > @{$self->{'TREE'}{'CHILDREN'}[1]{'LEAVES'}})
		{$chosen_index = 0;}
		elsif(@{$self->{'TREE'}{'CHILDREN'}[0]{'LEAVES'}} < @{$self->{'TREE'}{'CHILDREN'}[1]{'LEAVES'}})
		{$chosen_index = 1;}
		elsif( join(' ', @{$self->{'TREE'}{'CHILDREN'}[0]{'LEAVES'}}) gt join(' ', @{$self->{'TREE'}{'CHILDREN'}[1]{'LEAVES'}}) )
		{$chosen_index = 0;}
		elsif( join(' ', @{$self->{'TREE'}{'CHILDREN'}[0]{'LEAVES'}}) lt join(' ', @{$self->{'TREE'}{'CHILDREN'}[1]{'LEAVES'}}) )
		{$chosen_index = 1;}
		
		if(!defined($chosen_index))
		{Carp::croak($sub." failed: can't find consistent way to unroot tree, exiting;");}
				
		$chosen_child = splice(@{$self->{'TREE'}{'CHILDREN'}}, $chosen_index, 1);
		
		undef my @grandchildren;
		foreach my $grandchild (@{$chosen_child->{'CHILDREN'}})
		{
			push(@grandchildren, $grandchild);			
			$grandchild->{"PARENT"} = $self->{'TREE'};
			
			if(defined($chosen_child->{'BRANCH-LABELS'}))
			{
				$grandchild->{'BRANCH-LABELS'} = (@{$grandchild->{'CHILDREN'}}>0) ? "\'#1\'" : "#1" ;	
			}
			
			$grandchild->{'BRANCH-LENGTHS'} += $chosen_child->{'BRANCH-LENGTHS'} if(defined($chosen_child->{'BRANCH-LENGTHS'}));
		}
		
		splice(@{$self->{'TREE'}{'CHILDREN'}}, $chosen_index, 0, @grandchildren);
		
		undef $chosen_child;
	}
	
	print($sub.": tree unrooted...\n" ) if($self->{'VERBOSE'});
}

=head2 remove_attributes

 Usage   : $obj->remove_attributes($attribute);
           $obj->remove_attributes(@attributes);
 Function: Removes the specified attribute(s) from the Codeml tree.
 Returns : none
 Args    : $attribute (string) OR @attributes (array of strings)

=cut

sub remove_attributes
{
	my ($self, @attributes) = @_;
	my $sub = $self->get_nickname().substr((caller(0))[3], MODULE_NAME_LENGTH);
	
	if(!defined($self->{'TREE'}))
	{Carp::croak($sub." failed: no tree defined, exiting;");}
	
	foreach my $attribute (@attributes)
	{
		if(!(grep {$attribute eq $_} @supported_attributes))
		{Carp::croak($sub." failed: tree attribute \"".$attribute."\" not supported, exiting;");}
	}
	
	util_recurse_tree($self->{'TREE'}, \&util_node_remove_attributes, [@attributes]);
	
	print($sub.": removed attribute(s) \'".join("\', \'", @attributes)."\'...\n" ) if($self->{'VERBOSE'});
}

=head2 get_clade_strings

 Usage   : @clades = $obj->get_clade_strings(@regexes);
 Function: Gets list of Newick strings, each of which represents a monophyletic 
           clade within the tree containing genes that match one or more of the 
           specified regexes.
 Returns : @clades (array of strings)
 Args    : @regexes (array of strings)

=cut

sub get_clade_strings
{
	my($self, @regexes) = @_;
	my $sub = $self->get_nickname().substr((caller(0))[3], MODULE_NAME_LENGTH);
	
	if(!defined($self->{'TREE'}))
	{Carp::croak($sub." failed: no tree defined, exiting;");}
	
	if(@regexes==0)
	{Carp::croak($sub." failed: no regexes given, exiting;");}
				
	my @escaped_regexes = map { util_escape_meta_characters($_) } @regexes;

	my $clade_strings = util_recurse_tree($self->{'TREE'}, \&util_node_get_clade_strings, \@escaped_regexes);
	
	my @clades = split(/;/, $clade_strings);

	print($sub.": ".@clades." matching clades found...\n" ) if($self->{'VERBOSE'});
	
	return @clades;
}

=head2 get_tree_string

 Usage   : $tree_string = $obj->get_tree_string(@attributes);
 Function: Gets a string representing the Codeml tree, including any 
           specified attributes (e.g. 'BRANCH-LENGTHS', 'BRANCH-LABELS'). 
 Returns : $tree_string (string)
 Args    : @attributes (optional array of attributes)

=cut

sub get_tree_string
{
	my ($self, @attributes) = @_;
	my $sub = $self->get_nickname().substr((caller(0))[3], MODULE_NAME_LENGTH);
	
	if(!defined($self->{'TREE'}))
	{Carp::croak($sub." failed: no tree defined, exiting;");}
	
	foreach my $attribute (@attributes)
	{
		if(!(grep {$attribute eq $_} @supported_attributes))
		{Carp::croak($sub." failed: tree attribute \"".$attribute."\" not supported, exiting;");}
	}

	print($sub.": tree string retrieved...\n" ) if($self->{'VERBOSE'} && ( !defined((caller(1))[3]) || (substr((caller(1))[3], 0, 6) ne 'Tree::') ) );
	
	return util_recurse_tree($self->{'TREE'}, \&util_node_get_string, \@attributes).";";
}

=head2 get_ladderised_tree_string

 Usage   : $ladderised_tree_string = $obj->get_ladderised_tree_string(@attributes);
 Function: Gets a ladderised string representing the Codeml tree, including any 
           specified attributes (e.g. 'BRANCH-LENGTHS', 'BRANCH-LABELS'). 
 Returns : $ladderised_tree_string (string)
 Args    : @attributes (optional array of attributes)

=cut

sub get_ladderised_tree_string
{
	my ($self, @attributes) = @_;
	my $sub = $self->get_nickname().substr((caller(0))[3], MODULE_NAME_LENGTH);
	
	if(!defined($self->{'TREE'}))
	{Carp::croak($sub." failed: no tree defined, exiting;");}
	
	foreach my $attribute (@attributes)
	{
		if(!(grep {$attribute eq $_} @supported_attributes))
		{Carp::croak($sub." failed: tree attribute \"".$attribute."\" not supported, exiting;");}
	}

	print($sub.": ladderised tree string retrieved...\n" ) if($self->{'VERBOSE'});

	return util_recurse_tree($self->{'TREE'}, \&util_node_get_string, \@attributes, \&util_ladderise_children).";";
}

=head2 get_file_string

 Usage   : $text = $obj->get_file_string();
 Function: Gets text of tree file from which this CodemlWrapper::Tree object was loaded.
 Returns : $text (string)
 Args    : none

=cut

sub get_file_string
{
	my($self) = @_;
	my $sub = $self->get_nickname().substr((caller(0))[3], MODULE_NAME_LENGTH);
			
	if(!defined($self->{'TREE'}))
	{Carp::croak($sub." failed: no tree defined, exiting;");}

	print($sub.": file string retrieved...\n" ) if($self->{'VERBOSE'});
	
	return $self->{'FILE-STRING'};
}

=head2 tree_match

 Usage   : if($obj->tree_match($other, @attributes)) {...}
 Function: Tests if another CodemlWrapper::Tree object is identical to this one with 
           respect to topology, in addition to the specified attributes 
           (e.g. 'BRANCH-LENGTHS', 'BRANCH-LABELS'). 
 Returns : Boolean indicating if other CodemlWrapper::Tree object is identical to this one
 Args    : $other (CodemlWrapper::Tree object), @attributes (optional array of strings)

=cut

sub tree_match
{
	my($self, $other, @attributes) = @_;
	my $sub = $self->get_nickname().substr((caller(0))[3], MODULE_NAME_LENGTH);
	
	if(!defined($self->{'TREE'}))
	{Carp::croak($sub." failed: no tree defined, exiting;");}
	
	if(!defined($other))
	{Carp::croak($sub." failed: no other \'CodemlWrapper::Tree\' object given, exiting;");}
	
	if(ref($other) ne "Tree")
	{Carp::croak($sub." failed: can only compare to a \'CodemlWrapper::Tree\' object, exiting;");}
	
	if(!defined($other->{'TREE'}))
	{Carp::croak($sub." failed: other tree not defined, exiting;");}
	
	foreach my $attribute (@attributes)
	{
		if(!(grep {$attribute eq $_} @supported_attributes))
		{Carp::croak($sub." failed: tree attribute \"".$attribute."\" not supported, exiting;");}
	}

	my $self_tree_string  = util_recurse_tree($self->{'TREE'}, \&util_node_get_string, \@attributes, \&util_ladderise_children);
	my $other_tree_string = util_recurse_tree($other->{'TREE'}, \&util_node_get_string, \@attributes, \&util_ladderise_children);
	my $status = ($self_tree_string eq $other_tree_string);
	
	print($sub.": ".($status ? "TRUE" : "FALSE")."...\n" ) if($self->{'VERBOSE'});
			
	return $status;
}

=head2 set_output_format

 Usage   : $obj->set_output_format($output_format);
 Function: Sets the tree format used for file output. The supported tree output 
           formats are 'PAML', 'PHYLIP' and 'NEWICK'.
 Returns : none
 Args    : $output_format (string)

=cut

sub set_output_format
{	
	my($self, $output_format) = @_;
	my $sub = $self->get_nickname().substr((caller(0))[3], MODULE_NAME_LENGTH);
	
	$output_format = uc($output_format);
	
	if(! grep { $_ eq $output_format } @supported_file_types )
	{Carp::croak($sub." failed: file type \"".$output_format."\" not supported, exiting;");}
	
	$self->{'OUTPUT-FORMAT'} = $output_format;	

	print($sub.": output format set to \'".$output_format."\'...\n" ) if($self->{'VERBOSE'} && ( !defined((caller(1))[3]) || (substr((caller(1))[3], 0, 6) ne 'Tree::') ) );
}

=head2 get_output_format

 Usage   : $format = $obj->get_output_format();
 Function: Gets the tree format used for file output.
 Returns : $format (string)
 Args    : none

=cut

sub get_output_format
{	
	my($self) = @_;
	my $sub = $self->get_nickname().substr((caller(0))[3], MODULE_NAME_LENGTH);

	print($sub.": \'".$self->{'OUTPUT-FORMAT'}."\'...\n" ) if($self->{'VERBOSE'} && ( !defined((caller(1))[3]) || (substr((caller(1))[3], 0, 6) ne 'Tree::') ) );

	return $self->{'OUTPUT-FORMAT'};		
}

sub DESTROY
{
	my($self) = @_;
			
	%$self = ();
}

sub get_nickname
{
	my ($self) = @_;

	#INTERNAL USE ONLY!
	if( !defined((caller(1))[3]) || (substr((caller(1))[3], 0, 6) ne 'Tree::') )
	{Carp::croak((caller(0))[3]." failed: can only be called by \'CodemlWrapper::Tree\' subroutines, exiting;");}
	
	my $nickname = 'Tree';

	if( defined($self->{'TREE'}) )
	{
		my $all = scalar(@{$self->{'TREE'}{'LEAVES'}});
		
		my $fg = defined($self->{'FOREGROUND'}) ? scalar(@{$self->{'FOREGROUND'}}).'/' : '';
		
		$nickname .= "(".$fg.$all.")";
	}

	return $nickname;
}

sub parse_tree
{
	my($self, $tree_string) = @_;
	
	#INTERNAL USE ONLY!
	if( !defined((caller(1))[3]) || (caller(1))[3]!~/^Tree::([^:]+)$/ )
	{Carp::croak((caller(0))[3]." failed: can only be called by \'CodemlWrapper::Tree\' subroutines, exiting;");}
	
	#Remove semi-colon from tree string.
	if($tree_string=~/(.+);/)
	{$tree_string = $1;}
	else
	{Carp::croak((caller(0))[3]." failed: no semi-colon present, exiting;");}
	
	undef my %tree_hash;
	my $root_ref = \%tree_hash;	
	$root_ref->{'STRING'} = $tree_string;

	my @active_list = $root_ref;
	while(@active_list>0)
	{	
		my $tree_ref = shift(@active_list);
		
		if($tree_ref->{'STRING'}=~/^\s*(\(((.)+,(.)+)\))?\s*([^\(\)\[\];:,\$#]+)?\s*(:((\d)+(\.(\d)+)?))?\s*(\'(#|\$)(0|1)\'|\"(#|\$)(0|1)\"|(#|\$)(0|1))?\s*$/)
		{
			my $subtree		 	= $1;		#Regex: (\(((.)+,(.)+)\))?
			my $node_name 		= $5;		#Regex: ([^\(\)\[\];:,\$#]+)?	
			my $branch_length 	= $7;		#Regex: (:((\d)+(\.(\d)+)?))?		
			my $node_label  	= $11;		#Regex: (\'(#|\$)(0|1)\'|\"(#|\$)(0|1)\"|(#|\$)(0|1))?			
			
			undef my @leaves;
			while($tree_ref->{'STRING'}=~/(^|\(|,)([^\(,\);]+)/g)	
			{
				my $leaf = $2;
				$leaf=~s/((:|#|\$)(.+))//g;
				push(@leaves, $leaf) if($leaf);
			}			
			$tree_ref->{'LEAVES'} = [sort(@leaves)];
			
			if(defined($node_name) && ($node_name ne ''))
			{
				$node_name=~s/(^\s+|\s+$)//g;
				$tree_ref->{'NODE-NAME'} = $node_name;
			}
								
			if(defined($branch_length) && ($branch_length ne ''))
			{$tree_ref->{'BRANCH-LENGTHS'} = $branch_length;}	
			
			if(defined($node_label) && ($node_label ne ''))
			{	
				$tree_ref->{'BRANCH-LABELS'} = $node_label;
			
				if( ($node_label ne "#1") && ($node_label ne "\'#1\'") )
				{Carp::croak((caller(0))[3]." failed: label \"".$node_label."\" is not supported in module CodemlWrapper::Tree, exiting;");}
			}	
			
			my @child_strings = Tree::get_child_strings($subtree);
			
			if(@child_strings > 1)
			{					
				foreach my $child_string (@child_strings)
				{
					push( @{$tree_ref->{'CHILDREN'}}, {} );
					$tree_ref->{'CHILDREN'}[-1]{'STRING'} = $child_string;
					$tree_ref->{'CHILDREN'}[-1]{'PARENT'} = $tree_ref;
					push(@active_list, $tree_ref->{'CHILDREN'}[-1]);
				}
			}
			elsif(@child_strings == 1)
			{
				Carp::croak((caller(0))[3]." failed: tree has non-furcating node (Too many parentheses? Missing comma?), exiting;");
			}
			else
			{
				push(@{$self->{'LEAF-NODES'}}, $tree_ref);
			}
		}
		else
		{
			Carp::croak((caller(0))[3]." failed: subtree appears to be invalid \"".$tree_ref->{'STRING'}."\"");
		}
		
		delete $tree_ref->{'STRING'};
	}
	
	return $root_ref;
}

sub get_child_strings
{
	my ($tree_string) = @_;
	undef my @child_strings;
	
	if(defined($tree_string) && ($tree_string ne ''))
	{	
		$tree_string=~s/(^\()|(\);?$)//g;
		
		my $level = 0;
		my ($s, $l) = (undef, undef);
		for( my $i=0 ; $i<length($tree_string) ; $i++ )
		{
			my $c = substr($tree_string, $i, 1);
								
			if($c eq "("){$level++;}					
			elsif($c eq ")" ){$level--;}
								
			$s = $i if(!defined($s));
													
			if( ($c eq ",") && ($level==0) )
			{$l = ($i-$s);}
			elsif($i==(length($tree_string)-1))
			{$l = (length($tree_string)-$s);}		
			
			if(defined($s) && defined($l))
			{	
				my $cs = substr($tree_string, $s, $l);
				
				#Trim leading and trailing whitespace.
				$cs=~s/(^\s+|\s+$)//g;
				
				push(@child_strings, $cs);
				($s, $l) = (undef) x 2;
			}
		}
	}
	
	return @child_strings;
}

#Tree recursion subroutine.
sub util_recurse_tree
{
	my($ref_tree, $ref_func, $ref_args, $ref_sort) = @_;
	
	#INTERNAL USE ONLY!
	if( !defined((caller(1))[3]) || (caller(1))[3]!~/^Tree::([^:]+)$/ )
	{Carp::croak((caller(0))[3]." failed: can only be called by \'CodemlWrapper::Tree\' subroutines, exiting;");}
	
	if(defined($ref_tree))
	{		
		undef my @children_results;
		
		if(defined($ref_tree->{'CHILDREN'}))
		{
			my @children = @{$ref_tree->{'CHILDREN'}};
			
			if(defined($ref_sort))
			{@children = $ref_sort->(@children);}
		
			foreach my $ref_child (@children)
			{
				my @child_results = util_recurse_tree($ref_child, $ref_func, $ref_args, $ref_sort);
				
				push(@children_results, @child_results);
			}
		}		
		
		return $ref_func->($ref_tree, $ref_args, @children_results);
	}
}

#Subroutine for recursion
sub util_node_get_string
{
	my($node, $ref_args, @child_strings) = @_;
	
	#INTERNAL USE ONLY!
	if( !defined((caller(1))[3]) || (caller(1))[3]!~/^Tree::([^:]+)$/ )
	{Carp::croak((caller(0))[3]." failed: can only be called by \'CodemlWrapper::Tree\' subroutines, exiting;");}
	
	my $string = '';
	undef my @visible_attributes;
				
	if(defined($ref_args))
	{@visible_attributes = @{$ref_args};}
					
	if(defined($node))
	{					
		if(@child_strings>0)
		{						
			$string .= "(".join(",", @child_strings).")";
		}
		
		$string .= $node->{'NODE-NAME'} if defined($node->{'NODE-NAME'});	
		
		if( defined($node->{'BRANCH-LENGTHS'}) && (grep { $_ eq 'BRANCH-LENGTHS' } @visible_attributes) )
		{$string .= ":".$node->{'BRANCH-LENGTHS'};}
		
		if( defined($node->{'BRANCH-LABELS'}) && (grep { $_ eq 'BRANCH-LABELS' } @visible_attributes) )
		{$string .= $node->{'BRANCH-LABELS'};}					
	}
	
	return $string;	
}

#Subroutine for recursion
sub util_node_get_clade_strings
{
	my($node, $ref_args, @child_clade_strings) = @_;
	
	#INTERNAL USE ONLY!
	if( !defined((caller(1))[3]) || (caller(1))[3]!~/^Tree::([^:]+)$/ )
	{Carp::croak((caller(0))[3]." failed: can only be called by \'CodemlWrapper::Tree\' subroutines, exiting;");}
	
	my $clade_string = '';
	undef my @regexes;
				
	if(defined($ref_args))
	{@regexes = @{$ref_args};}

	if(defined($node))
	{				
		my $mismatch = 0;
		foreach my $leaf (@{$node->{'LEAVES'}})
		{
			if( !(grep { $leaf=~/$_/  } @regexes) )
			{$mismatch = 1; last;}
		}		
		
		if($mismatch)
		{
			@child_clade_strings = grep { $_ } @child_clade_strings;
			
			if(@child_clade_strings>0)
			{
				$clade_string .= join(";", @child_clade_strings);
			}
		}
		else
		{
			if(@child_clade_strings>0)
			{						
				$clade_string .= "(".join(",", @child_clade_strings).")";
			}
			
			if(!defined($node->{'CHILDREN'}))
			{
				$clade_string .= $node->{'NODE-NAME'} if defined($node->{'NODE-NAME'});	
			}
		}
	}
	
	return $clade_string;	
}
	
#Subroutine for recursion
sub util_node_remove_attributes
{
	my($node, $ref_args) = @_;
	
	#INTERNAL USE ONLY!
	if( !defined((caller(1))[3]) || (caller(1))[3]!~/^Tree::([^:]+)$/ )
	{Carp::croak((caller(0))[3]." failed: can only be called by \'CodemlWrapper::Tree\' subroutines, exiting;");}
	
	undef my @unwanted_attributes;
	
	if(defined($ref_args))
	{@unwanted_attributes = @{$ref_args};}
	
	if( defined($node->{'BRANCH-LABELS'}) && (grep { $_ eq 'BRANCH-LABELS' } @unwanted_attributes) )
	{delete $node->{'BRANCH-LABELS'};}
	
	if( defined($node->{'BRANCH-LENGTHS'}) && (grep { $_ eq 'BRANCH-LENGTHS' } @unwanted_attributes) )
	{delete $node->{'BRANCH-LENGTHS'};}
}

#Subroutine for recursion
sub util_node_get_foreground_genes
{
	my($node, $ref_args, @child_foreground_genes) = @_;
	
	#INTERNAL USE ONLY!
	if( !defined((caller(1))[3]) || (caller(1))[3]!~/^Tree::([^:]+)$/ )
	{Carp::croak((caller(0))[3]." failed: can only be called by \'CodemlWrapper::Tree\' subroutines, exiting;");}
	
	undef my @foreground_genes;
	
	if(defined($node->{'BRANCH-LABELS'}) && ($node->{'BRANCH-LABELS'}=~/^(#1|'#1')$/) )
	{		
		push(@foreground_genes, @{$node->{'LEAVES'}});
	}
	else
	{
		push(@foreground_genes, @child_foreground_genes);
	}
	
	return @foreground_genes;
}

#Subroutine for recursion
sub util_node_coalesce_foreground_branch_labels
{
	my($node) = @_;
	
	#INTERNAL USE ONLY!
	if( !defined((caller(1))[3]) || (caller(1))[3]!~/^Tree::([^:]+)$/ )
	{Carp::croak((caller(0))[3]." failed: can only be called by \'CodemlWrapper::Tree\' subroutines, exiting;");}
	
	if(defined($node->{"PARENT"}) && defined($node->{'CHILDREN'}))
	{
		my @children = @{$node->{'CHILDREN'}};
		my @labelled = grep { defined($_->{'BRANCH-LABELS'}) && ($_->{'BRANCH-LABELS'}=~/^(#1|'#1')$/) } @children;
				
		if(@labelled==@children)
		{
			foreach my $child (@children)
			{delete $child->{'BRANCH-LABELS'};}
			
			$node->{'BRANCH-LABELS'} = "\'#1\'";
		}
	}
}

#Subroutine for sorting children
sub util_ladderise_children
{
	my(@children) = @_;
	
	#INTERNAL USE ONLY!
	if( !defined((caller(1))[3]) || (caller(1))[3]!~/^Tree::([^:]+)$/ )
	{Carp::croak((caller(0))[3]." failed: can only be called by \'CodemlWrapper::Tree\' subroutines, exiting;");}
	
	#This relies on the fact that the leaves are sorted when the tree is being input. 
	my @ladderised_children = sort { (@{$b->{'LEAVES'}} <=> @{$a->{'LEAVES'}}) || (join('', @{$a->{'LEAVES'}}) cmp join('', @{$b->{'LEAVES'}})) } @children;		
	
	return @ladderised_children;
}

sub util_write_file
{
	my($file, @string_data) = @_;
	
	#INTERNAL USE ONLY!
	if( !defined((caller(1))[3]) || (caller(1))[3]!~/^Tree::([^:]+)$/ )
	{Carp::croak((caller(0))[3]." failed: can only be called by \'CodemlWrapper::Tree\' subroutines, exiting;");}
	
	my $string = '';
	
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
	if( !defined((caller(1))[3]) || (caller(1))[3]!~/^Tree::([^:]+)$/ )
	{Carp::croak((caller(0))[3]." failed: can only be called by \'CodemlWrapper::Tree\' subroutines, exiting;");}
		
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

sub util_escape_meta_characters
{
	my($string) = @_;
	
	#INTERNAL USE ONLY!
	if( !defined((caller(1))[3]) || (caller(1))[3]!~/^Tree::([^:]+)$/ )
	{Carp::croak((caller(0))[3]." failed: can only be called by \'CodemlWrapper::Tree\' subroutines, exiting;");}
		
	my @meta_chars = ( '\^', '\$', '\.', '\+', '\?', '\*', '\{', '\}', '\(', '\)', '\/', '\|', '\[', '\]' );
	
	foreach my $char (@meta_chars)
	{$string =~ s/$char/$char/;}
		
	return $string;
}

}

1
