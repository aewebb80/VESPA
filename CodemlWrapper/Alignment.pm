#!/usr/bin/perl -w
use strict;
use Carp qw(carp cluck croak confess);
#CodemlWrapper::Alignment v1.01 by Tom Walsh (2012/07/12)

=head1 NAME

CodemlWrapper::Alignment - Wrapper for a Codeml nucleotide alignment file.

=head1 SYNOPSIS

	use CodemlWrapper::Alignment;
		
	#Create new CodemlWrapper::Alignment object, load from file.
	my $alignment = Alignment->new('align.fasta');	
		
	#Set gap tolerance to two less than the number of sequences.
	my $gap_tolerance = $alignment->get_sequence_count() - 2;
	
	#Get sites for which the number of sequences with a gap exceeds the gap tolerance.
	my @gap_sites = $alignment->get_gap_sites($gap_tolerance);
	
	#Remove gap sites.
	$alignment->remove_sites(@gap_sites);
	
	#Save the Codeml alignment file in sequential Phylip format to the current directory. 
	$alignment->save('align.phy', 'PHYLIP-SEQUENTIAL');	

=head1 DESCRIPTION

Codeml is a component of Ziheng Yang's PAML package. Codeml performs selection 
analysis on gene sequences, allowing episodes of positive selection of a gene on 
evolutionary timescales to be inferred. For information on PAML see the PAML manual.

CodemlWrapper::Alignment is a Perl wrapper module for the nucleotide alignment file used 
by Codeml (PAML4.4e). It helps with processing of Codeml nucleotide alignments 
and is used mainly in the context of the CodemlWrapper::Job module to handle the Codeml 
alignment file for each task in a Codeml job. A Codeml alignment file is written 
to disk by calling L</save>, and can be read from disk using L</load>. 

On loading from an alignment file, lowercase letters in each sequence are converted 
to uppercase, and stop codons are removed from the end of sequences, if present 
in all sequences. The format of an input file is inferred from its content. Supported 
formats include PHYLIP sequential, PHYLIP interleaved and FASTA format. Within the 
module, these are named 'PHYLIP-SEQUENTIAL', 'PHYLIP-INTERLEAVED' and 'FASTA' respectively. 
Alignments can be saved in HTML using the 'HTML' format. (NB: HTML files saved by 
CodemlWrapper::Alignment aren't really conforming to any sort of bioinformatics 
standard, and should only be used by the CodemlWrapper::Alignment module.) The 
default output format is PHYLIP sequential. The output format of a Codeml alignment 
is set using L</set_output_format>, and can be obtained using L</get_output_format>. 

As described in the PAML manual, the PHYLIP sequential and interleaved format used 
by PAML (and this module) differs slightly from the standard with respect to sequence 
headers. The standard PHYLIP format limits header length to 10 characters, but PAML 
and CodemlWrapper::Alignment support headers of up to 30 characters in length, where two 
or more consecutive spaces indicate the end of the header. 

Positions in the alignment are given in terms of codon sites, so for example 
a nucleotide alignment of length 15 has 5 codon sites, numbered 1-5. 

Ambig sites and gap sites can be identified using L</get_ambig_sites> and L</get_gap_sites>, 
respectively. A tolerance parameter can be set in either case, such that a tolerance 
of zero will identify all ambig/gap sites, a tolerance of one will identify all sites 
for which more than one sequence has an ambig/gap, and so on. If no tolerance is given, 
a tolerance of zero is used.

Sites can be removed using L</remove_sites>. If you do this, it's critical 
that you record the changes and review the resulting alignment. This can be done 
by saving the original alignment in 'HTML' format with sites masked instead 
of removed, using L</set_mask_sites>.

CodemlWrapper::Alignment handles site mappings between three different forms of 
the sequence: 'ALIGNMENT' (the full alignment sequence), 'UNGAPPED' (the alignment 
sequence without gaps) and 'UNMASKED' (the alignment sequence without masked sites).

Codeml outputs positive site positions in terms of the input Codeml alignment, but 
often we're more interested in site positions in the ungapped sequence of the alignment. 
This can be done by mapping sites from the 'ALIGNMENT' to 'UNGAPPED' forms of a sequence. 

If the Codeml alignment was taken from a longer alignment from which gap/ambig
sites had been removed, we might be interested in mapping positive site positions 
to that original alignment. This is possible if the original alignment was saved
in 'HTML' format with gap/ambig sites masked. The 'UNMASKED' form of that HTML 
alignment corresponds to the Codeml alignment, so sites given in terms of the 
Codeml alignment can be mapped to the original alignment by mapping from the 
'UNMASKED' to 'ALIGNMENT' forms of a sequence in the original HTML alignment.

HTML alignments can also be used to display positively selected sites: these can 
be set using L</set_pos_sites>.

CodemlWrapper::Alignment usually runs silently. To see full output, include 
'Verbose' as an argument when creating a L</new> object.

=head1 OBJECT METHODS

=cut

package Alignment;
{

use constant TRUE  => 1;
use constant FALSE => 0;
use constant MIN_SEQUENCES => 7;
use constant MODULE_NAME_LENGTH => 9;
use constant DEFAULT_OUTPUT_FORMAT => 'PHYLIP-SEQUENTIAL';

my @supported_file_types = ('FASTA', 'PHYLIP-INTERLEAVED', 'PHYLIP-SEQUENTIAL', 'HTML');

my @html_boilerplate_start = 
(
	'<html>',
	'<head>',
	'<style>',
	'<!--',
	'body {background-color:#FFFFFF;color:#000000}',
	'span.S {background-color:#FFFFFF}',
	'span.M {background-color:#000000}',
	'span.F {background-color:#FF5555}',
	'span.B {background-color:#55AAFF}',
	'span.N {background-color:#888888}',
	'-->',
	'</style>',
	'</head>',
	'<body>',
	'<pre>'
);

my @html_boilerplate_end = 
(
	'</pre>',
	'</body>',
	'</html>'
);

my %span = 
(
	'S' => '<span class="S">',
	'M' => '<span class="M">',
	'F' => '<span class="F">',
	'B' => '<span class="B">',
	'N' => '<span class="N">'
);

my $se = '</span>';

#Init standard codon hash.
my %standard_codons = 
( 	
	'GCT' => 'A', 'GCC' => 'A', 'GCA' => 'A', 'GCG' => 'A',						#ALANINE
	'TGT' => 'C', 'TGC' => 'C', 												#CYSTEINE
	'GAT' => 'D', 'GAC' => 'D', 												#ASPARTIC ACID
	'GAA' => 'E', 'GAG' => 'E',													#GLUTAMIC ACID
	'TTT' => 'F', 'TTC' => 'F', 												#PHENYLALANINE
	'GGT' => 'G', 'GGC' => 'G', 'GGA' => 'G', 'GGG' => 'G',						#GLYCINE
	'CAT' => 'H', 'CAC' => 'H',													#HISTIDINE
	'ATT' => 'I', 'ATC' => 'I', 'ATA' => 'I', 									#ISOLEUCINE
	'AAA' => 'K', 'AAG' => 'K',													#LYSINE
	'CTT' => 'L', 'CTC' => 'L', 'CTA' => 'L', 
	'CTG' => 'L', 'TTA' => 'L', 'TTG' => 'L',									#LEUCINE
	'ATG' => 'M',																#METHIONINE
	'AAT' => 'N', 'AAC' => 'N', 												#ASPARAGINE
	'CCT' => 'P', 'CCC' => 'P', 'CCA' => 'P', 'CCG' => 'P',						#PROLINE
	'CAA' => 'Q', 'CAG' => 'Q',													#GLUTAMINE
	'CGT' => 'R', 'CGC' => 'R', 'CGA' => 'R', 									#ARGININE
	'CGG' => 'R', 'AGA' => 'R', 'AGG' => 'R',
	'TCT' => 'S', 'TCC' => 'S', 'TCA' => 'S',  									#SERINE
	'TCG' => 'S', 'AGT' => 'S', 'AGC' => 'S',
	'ACT' => 'T', 'ACC' => 'T', 'ACA' => 'T', 'ACG' => 'T',						#THREONINE
	'GTT' => 'V', 'GTC' => 'V', 'GTA' => 'V', 'GTG' => 'V',						#VALINE
	'TGG' => 'W',																#TRYPTOPHAN
	'TAT' => 'Y', 'TAC' => 'Y', 												#TYROSINE
	'---' => '-'
);

#Init redundant codon hash.
my %redundant_codons = 
(
	#2-FOLD DEGENERATE
	'TGY' => 'C', 																#CYSTEINE
	'GAY' => 'D', 																#ASPARTIC ACID
	'GAR' => 'E', 																#GLUTAMIC ACID
	'TTY' => 'F', 																#PHENYLALANINE
	'CAY' => 'H', 																#HISTIDINE
	'AAR' => 'K', 																#LYSINE
	'AAY' => 'N', 																#ASPARAGINE
	'CAR' => 'Q', 																#GLUTAMINE
	'AGY' => 'S', 																#SERINE
	'TAY' => 'Y', 																#TYROSINE
		
	#3-FOLD DEGENERATE
	'ATH' => 'I', 'ATY' => 'I', 'ATM' => 'I', 'ATW' => 'I',						#ISOLEUCINE
	
	#4-FOLD DEGENERATE
	'GCN' => 'A', 'GCR' => 'A', 'GCY' => 'A', 'GCK' => 'A',   					#ALANINE
	'GCM' => 'A', 'GCS' => 'A', 'GCW' => 'A', 'GCB' => 'A',  
	'GCD' => 'A', 'GCH' => 'A', 'GCV' => 'A',	
	'GGN' => 'G', 'GGR' => 'G', 'GGY' => 'G', 'GGK' => 'G',   					#GLYCINE
	'GGM' => 'G', 'GGS' => 'G', 'GGW' => 'G', 'GGB' => 'G', 
	'GGD' => 'G', 'GGH' => 'G', 'GGV' => 'G', 	
	'CCN' => 'P', 'CCR' => 'P', 'CCY' => 'P', 'CCK' => 'P', 					#PROLINE
	'CCM' => 'P', 'CCS' => 'P', 'CCW' => 'P', 'CCB' => 'P', 
	'CCD' => 'P', 'CCH' => 'P', 'CCV' => 'P', 
	'TCN' => 'S', 'TCR' => 'S', 'TCY' => 'S', 'TCK' => 'S',   					#SERINE
	'TCM' => 'S', 'TCS' => 'S', 'TCW' => 'S', 'TCB' => 'S', 
	'TCD' => 'S', 'TCH' => 'S', 'TCV' => 'S',
	'ACN' => 'T', 'ACR' => 'T', 'ACY' => 'T', 'ACK' => 'T', 					#THREONINE
	'ACM' => 'T', 'ACS' => 'T', 'ACW' => 'T', 'ACB' => 'T', 
	'ACD' => 'T', 'ACH' => 'T', 'ACV' => 'T', 
	'GTN' => 'V', 'GTR' => 'V', 'GTY' => 'V', 'GTK' => 'V',  					#VALINE
	'GTM' => 'V', 'GTS' => 'V', 'GTW' => 'V', 'GTB' => 'V', 
	'GTD' => 'V', 'GTH' => 'V', 'GTV' => 'V', 
		
	#'6-FOLD' DEGENERATE
	'TTR' => 'L', 'YTA' => 'L', 'YTG' => 'L', 'CTN' => 'L', 'CTR' => 'L',   	#LEUCINE
	'CTY' => 'L', 'CTK' => 'L', 'CTM' => 'L', 'CTS' => 'L', 'CTW' => 'L',  
	'CTB' => 'L', 'CTD' => 'L', 'CTH' => 'L', 'CTV' => 'L',
	'AGR' => 'R', 'MGA' => 'R', 'MGG' => 'R', 'CGN' => 'R', 'CGR' => 'R',   	#ARGININE
	'CGY' => 'R', 'CGK' => 'R', 'CGM' => 'R', 'CGS' => 'R', 'CGW' => 'R', 
	'CGB' => 'R', 'CGD' => 'R', 'CGH' => 'R', 'CGV' => 'R'
);

=head2 new

 Usage   : $obj = CodemlWrapper::Alignment->new();
           $obj = CodemlWrapper::Alignment->new($file);
           $obj = CodemlWrapper::Alignment->new('Verbose');
           $obj = CodemlWrapper::Alignment->new($file, 'Verbose');
 Function: Creates a new CodemlWrapper::Alignment object. If a file name is specified,  
           loads CodemlWrapper::Alignment object from that file. If 'Verbose' is 
           specified, the CodemlWrapper::Alignment object will have verbose output.
 Returns : $obj (CodemlWrapper::Alignment object)
 Args    : $file (optional file name), 'Verbose' (optional string)

=cut

sub new
{
	my($self, @arguments) = @_;
	undef my %Alignment;
	$self = \%Alignment;
	
	bless($self, 'Alignment');
	
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
 Function: Creates a copy of a CodemlWrapper::Alignment object identical to the original.
 Returns : $copy (CodemlWrapper::Alignment object)
 Args    : none

=cut

sub copy
{
	my($self) = @_;
	my $sub = $self->get_nickname().substr((caller(0))[3], MODULE_NAME_LENGTH);
	
	undef my %Alignment;
	my $copy = \%Alignment;
	
	bless($copy, 'Alignment');

	$copy->{'VERBOSE'} = $self->{'VERBOSE'};
	
	$copy->{'INPUT-FORMAT'} = $self->{'INPUT-FORMAT'};
	$copy->{'OUTPUT-FORMAT'} = $self->{'OUTPUT-FORMAT'};
	$copy->{'FILE-STRING'} = $self->{'FILE-STRING'};
	
	%{$copy->{'INDX2HDR'}} 	 = map { $_ => $self->{'INDX2HDR'}{$_} }   keys(%{$self->{'INDX2HDR'}});
	%{$copy->{'HDR2INDX'}}   = map { $_ => $self->{'HDR2INDX'}{$_} }   keys(%{$self->{'HDR2INDX'}});
	%{$copy->{'SEQS'}} = map { $_ => $self->{'SEQS'}{$_} } keys(%{$self->{'SEQS'}});
	
	%{$copy->{'FGS'}} = map { $_ => $self->{'FGS'}{$_} }   keys(%{$self->{'FGS'}});
	$copy->{'FORMATTING'} = $self->{'FORMATTING'};
	
	$copy->update_mappings();

	print($sub.": alignment copied...\n" ) if($self->{'VERBOSE'});

	return $copy;
}

=head2 load

 Usage   : $obj->load($file);
 Function: Loads CodemlWrapper::Alignment object from the specified file.
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

	delete $self->{'FORMATTING'};
	delete $self->{'INPUT-FORMAT'};
	delete $self->{'OUTPUT-FORMAT'};
	delete $self->{'FILE-STRING'};
	delete $self->{'INDX2HDR'};
	delete $self->{'HDR2INDX'};
	delete $self->{'SEQS'};
	delete $self->{'FGS'};	

	$self->{'OUTPUT-FORMAT'} = DEFAULT_OUTPUT_FORMAT;
	
	my @file_array = util_read_file($file);
	@file_array = grep { $_!~/^\s*$/ } @file_array;	
			
	my $file_string = join("\n", @file_array);
	$self->{'FILE-STRING'} = $file_string;		#THIS IS READ-ONLY FROM THIS POINT ON!!
	
	my $entry_index = 1;
	my $stop_count = 0;
	
	my $first_line = $file_array[0];
	if($first_line=~/^>.+/)
	{
		$self->{'INPUT-FORMAT'} = 'FASTA';
		undef my $expected_length;
		
		while($file_string=~/>([^\n]+)\s*\n([^>]+)/g)
		{
			my $seqname  = $1;
			my $sequence = $2;
											
			if(defined($self->{'SEQS'}{$seqname}))
			{Carp::croak($sub." failed: name of entry ".$entry_index." (\'".$seqname."\') is duplicate of entry \'".$self->{'HDR2INDX'}{$seqname}."\', exiting;");}
			
			$sequence = util_trim_sequence($sequence);	
							
			if(!defined($expected_length))
			{
				$expected_length = length($sequence);
			}
			elsif(length($sequence) != $expected_length)
			{
				Carp::croak($sub." failed: FASTA sequences have different lengths, exiting;");	
			}
			
			if(length($sequence)==0)
			{Carp::croak($sub." failed: entry ".$entry_index." (\'".$seqname."\') has no sequence, exiting;");}
	
			if(length($sequence) % 3 != 0)
			{Carp::croak($sub." failed: sequence ".$entry_index." (\'".$seqname."\') has incomplete reading frame, exiting;");}
					
			if( (substr($sequence, -3, 3) eq "TGA") || (substr($sequence, -3, 3) eq "TAG") || (substr($sequence, -3, 3) eq "TAA") )
			{
				$sequence = substr($sequence, 0, -3);
				$stop_count++;
			}
			
			while($sequence=~/(TGA|TAG|TAA)/g)
			{
				if(length($`) % 3 == 0)
				{Carp::croak($sub." failed: sequence ".$entry_index." (\'".$seqname."\') has internal STOP codon(s), exiting;");}
			}
						
			$self->{'SEQS'}{$seqname} = $sequence;
			
			$self->{'FORMATTING'} = 'S' x length($sequence);
			
			$self->{'HDR2INDX'}{$seqname} = $entry_index;
			$self->{'INDX2HDR'}{$entry_index} = $seqname;
			$entry_index++;
		}
				
		if($entry_index==0)
		{
			Carp::croak($sub." failed: no valid FASTA sequences in \'".$file."\', exiting;");
		}	
	}
	elsif($first_line=~/(\s*)(\d+)(\s+)(\d+)(\s+)(I)(\s*)/)
	{
		$self->{'INPUT-FORMAT'} = 'PHYLIP-INTERLEAVED';
		my $spec_seq_tally = $2;
		my $spec_seq_length = $4;
		
		my @phylip_array = @file_array[1..(@file_array-1)];
								
		if(@phylip_array == 0)
		{Carp::croak($sub." failed: no data present in \'".$file."\', exiting;");}
		
		if(@phylip_array < $spec_seq_tally)
		{Carp::croak($sub." failed: specified number of sequences (".$spec_seq_tally.") exceeds number of lines of data (".@phylip_array."), exiting;");}
		
		if($spec_seq_tally == 0)
		{Carp::croak($sub." failed: zero sequences specified, exiting;");}
		
		if($spec_seq_length == 0)
		{Carp::croak($sub." failed: sequence length of zero specified, exiting;");}
		
		if($spec_seq_length % 3 != 0)
		{Carp::croak($sub." failed: specified sequence length (".$spec_seq_length.") has incomplete reading frame, exiting;");}
					
		for(my $i=0 ; $i<$spec_seq_tally ; $i++)
		{
			undef my $seqname;
			undef my $sequence;
							
			if($phylip_array[$i]=~/^(([^\,\:\#\(\)\$\=]*?)([^\,\:\#\(\)\$\=\s]+?)([^\,\:\#\(\)\$\=]*?))(\s\s)(\s*)(.*)/)
			{
				$seqname  = $1;
				$sequence = $7;
			}
			else
			{
				Carp::croak($sub." failed: no valid sequence name on line ".($i+1)." of file \'".$file."\', exiting;");
			}
			
			if(length($seqname)>30)
			{Carp::croak($sub." failed: sequence name \'".$seqname."\' is longer than the maximum 30 characters allowed in Codeml, exiting;");}
			
			if(defined($self->{'SEQS'}{$seqname}))
			{Carp::croak($sub." failed: name of entry ".$entry_index." (\'".$seqname."\') is duplicate of entry \'".$self->{'HDR2INDX'}{$seqname}."\', exiting;");}
			
			my $next_line = $spec_seq_tally + $i;
			while(defined($phylip_array[$next_line]))
			{
				$sequence .= $phylip_array[$next_line];
				$next_line += $spec_seq_tally;
			}
			
			$sequence = util_trim_sequence($sequence);
							
			if(length($sequence) != $spec_seq_length)
			{Carp::croak($sub." failed: specified length of sequence \'".$seqname."\' (".$spec_seq_length.") differs from apparent length (".length($sequence)."), exiting;");}
			
			if(length($sequence)==0)
			{Carp::croak($sub." failed: PHYLIP entry ".$entry_index." (\'".$seqname."\') has no sequence, exiting;");}
	
			if(length($sequence) % 3 != 0)
			{Carp::croak($sub." failed: PHYLIP entry ".$entry_index." (\'".$seqname."\') has incomplete reading frame, exiting;");}
								
			if( (substr($sequence, -3, 3) eq "TGA") || (substr($sequence, -3, 3) eq "TAG") || (substr($sequence, -3, 3) eq "TAA") )
			{
				$sequence = substr($sequence, 0, -3);
				$stop_count++;
			}
			
			while($sequence=~/(TGA|TAG|TAA)/g)
			{
				if(length($`) % 3 == 0)
				{Carp::croak($sub." failed: sequence ".$entry_index." (\'".$seqname."\') has internal STOP codon(s), exiting;");}
			}	
					
			$self->{'SEQS'}{$seqname} = $sequence;	
			
			$self->{'FORMATTING'} = 'S' x length($sequence);
							
			$self->{'HDR2INDX'}{$seqname} = $entry_index;
			$self->{'INDX2HDR'}{$entry_index} = $seqname;
			$entry_index++;
		}
		
		if(keys(%{$self->{'SEQS'}}) != $spec_seq_tally)
		{Carp::croak($sub." failed: specified number of sequences (".$spec_seq_tally.") differs from apparent number (".keys(%{$self->{'SEQS'}})."), exiting;");}		
	}
	elsif($first_line=~/(\s*)(\d+)(\s+)(\d+)(\s*)(S\s*)?/)
	{
		$self->{'INPUT-FORMAT'} = 'PHYLIP-SEQUENTIAL';
		my $spec_seq_tally = $2;
		my $spec_seq_length = $4;
		
		my @phylip_array = @file_array[1..(@file_array-1)];
											
		if(@phylip_array == 0)
		{Carp::croak($sub." failed: no data present in \'".$file."\', exiting;");}
		
		if(@phylip_array < $spec_seq_tally)
		{Carp::croak($sub." failed: specified number of sequences (".$spec_seq_tally.") exceeds number of lines of data (".@phylip_array.") in file \'".$file."\', exiting;");}
		
		if($spec_seq_tally == 0)
		{Carp::croak($sub." failed: zero sequences specified, exiting;");}
		
		if($spec_seq_length == 0)
		{Carp::croak($sub." failed: sequence length of zero specified, exiting;");}
		
		if($spec_seq_length % 3 != 0)
		{Carp::croak($sub." failed: specified sequence length (".$spec_seq_length.") has incomplete reading frame, exiting;");}
					
		my $next_line = 0;
		for(my $i=0 ; $i<$spec_seq_tally ; $i++)
		{
			undef my $seqname;
			undef my $sequence;
			if($phylip_array[$next_line]=~/^(([^\,\:\#\(\)\$\=]*?)([^\,\:\#\(\)\$\=\s]+?)([^\,\:\#\(\)\$\=]*?))(\s\s)(\s*)(.*)/)
			{
				$seqname  = $1;
				$sequence = $7;
				$next_line++;
			}
			else
			{
				Carp::croak($sub." failed: no valid sequence name on line ".($next_line+1)." of file \'".$file."\', exiting;");
			}
			
			if(length($seqname)>30)
			{Carp::croak($sub." failed: sequence name \'".$seqname."\' is longer than the maximum 30 characters allowed in Codeml, exiting;");}
			
			if(defined($self->{'SEQS'}{$seqname}))
			{Carp::croak($sub." failed: name of entry ".$entry_index." (\'".$seqname."\') is duplicate of entry \'".$self->{'HDR2INDX'}{$seqname}."\', exiting;");}
			
			$sequence = util_trim_sequence($sequence);	
							
			while( (length($sequence)<$spec_seq_length) && defined($phylip_array[$next_line]) )
			{
				my $seq_chunk = $phylip_array[$next_line];
									
				$seq_chunk = util_trim_sequence($seq_chunk);	
				
				$sequence .= $seq_chunk;
				
				$next_line++;
			}
									
			if(length($sequence) != $spec_seq_length)
			{Carp::croak($sub." failed: specified length of sequence \'".$seqname."\' (".$spec_seq_length.") differs from apparent length (".length($sequence)."), exiting;");}
			
			if(length($sequence)==0)
			{Carp::croak($sub." failed: PHYLIP entry ".$entry_index." (\'".$seqname."\') has no sequence, exiting;");}
	
			if(length($sequence) % 3 != 0)
			{Carp::croak($sub." failed: PHYLIP entry ".$entry_index." (\'".$seqname."\') has incomplete reading frame, exiting;");}
						
			if( (substr($sequence, -3, 3) eq "TGA") || (substr($sequence, -3, 3) eq "TAG") || (substr($sequence, -3, 3) eq "TAA") )
			{
				$sequence = substr($sequence, 0, -3);
				$stop_count++;
			}
			
			while($sequence=~/(TGA|TAG|TAA)/g)
			{
				if(length($`) % 3 == 0)
				{Carp::croak($sub." failed: sequence ".$entry_index." (\'".$seqname."\') has internal STOP codon(s), exiting;");}
			}
			
			$self->{'SEQS'}{$seqname} = $sequence;	
			
			$self->{'FORMATTING'} = 'S' x length($sequence);
							
			$self->{'HDR2INDX'}{$seqname} = $entry_index;
			$self->{'INDX2HDR'}{$entry_index} = $seqname;
			$entry_index++;
		}
		
		if(keys(%{$self->{'SEQS'}}) != $spec_seq_tally)
		{Carp::croak($sub." failed: specified number of sequences (".$spec_seq_tally.") differs from apparent number (".keys(%{$self->{'SEQS'}})."), exiting;");}
	}
	elsif($first_line=~/<html>/)
	{
		$self->{'INPUT-FORMAT'} = 'HTML';
		undef my $expected_length;
		
		my $i = 0;
		while($i < @html_boilerplate_start)
		{ 
			if($file_array[$i] ne $html_boilerplate_start[$i])
			{Carp::croak($sub." failed: invalid HTML alignment file \'".$file."\', exiting;");}
			
			$i++;
		}
		
		my $j = @file_array - 1;
		my $x = @file_array - @html_boilerplate_end;
		while ($j >= $x )
		{ 
			if($file_array[$j] ne $html_boilerplate_end[$j-$x])
			{Carp::croak($sub." failed: invalid HTML alignment file \'".$file."\', exiting;");}
			
			$j--;
		}
		
		my @html_array =  @file_array[$i..$j];
		
		if(@html_array == 0)
		{Carp::croak($sub." failed: no data present in \'".$file."\', exiting;");}
		
		undef my %seq_data;
		undef my @seq_list;
		foreach my $line (@html_array)
		{
			if($line=~/(\S+)\s+(.+)/)
			{
				my ($name, $part) = ($1, $2);
				
				if( defined($seq_data{$name}) )
				{
					$seq_data{$name} .= $part;
				}
				else
				{
					$seq_data{$name} = $part;
					push(@seq_list, $name);
				}
			}
			else
			{
				Carp::croak($sub." failed: invalid HTML alignment file \'".$file."\', exiting;");
			}
		}

		foreach my $name (@seq_list)
		{
			my $html_sequence = $seq_data{$name};
			my $formatting = '';
			my $sequence = '';
			
			while($html_sequence=~/<span class=\"([^\"]+)\">([^<]+)<\/span>/g)
			{
				my ($fmt, $seq) = ($1, $2);
				
				if( ($fmt eq 'F') || ($fmt eq 'B') )
				{
					$self->{'FGS'}{$name} = $fmt;
					$fmt = 'P';
				}
				
				$formatting .= $fmt x length($seq);
				$sequence .= $seq;
			}
		
			$sequence = util_trim_sequence($sequence);
			
			if(!defined($expected_length))
			{
				$expected_length = length($sequence);
			}
			elsif(length($sequence) != $expected_length)
			{
				Carp::croak($sub." failed: sequences have different lengths, exiting;");	
			}
			
			if(!defined($self->{'FORMATTING'}))
			{
				$self->{'FORMATTING'} = $formatting;
			}
			elsif($self->{'FORMATTING'} ne $formatting)
			{
				Carp::croak($sub." failed: sequences have different masked/positive sites, exiting;");		
			}
		
			if(length($sequence)==0)
			{Carp::croak($sub." failed: HTML entry ".$entry_index." (\'".$name."\') has no sequence, exiting;");}

			if(length($sequence) % 3 != 0)
			{Carp::croak($sub." failed: HTML entry ".$entry_index." (\'".$name."\') has incomplete reading frame, exiting;");}
						
			if( (substr($sequence, -3, 3) eq "TGA") || (substr($sequence, -3, 3) eq "TAG") || (substr($sequence, -3, 3) eq "TAA") )
			{
				$sequence = substr($sequence, 0, -3);
				$stop_count++;
			}	
			
			while($sequence=~/(TGA|TAG|TAA)/g)
			{
				if(length($`) % 3 == 0)
				{Carp::croak($sub." failed: sequence ".$entry_index." (\'".$name."\') has internal STOP codon(s), exiting;");}
			}
			
			$self->{'SEQS'}{$name} = $sequence;

			$self->{'HDR2INDX'}{$name} = $entry_index;
			$self->{'INDX2HDR'}{$entry_index} = $name;			
			
			$entry_index++;
		}
	}
	else
	{
		Carp::croak($sub." failed: invalid Codeml alignment file \'".$file."\', exiting;");		
	}	
			
	if($stop_count!=0 && $stop_count!=keys(%{$self->{'SEQS'}}))
	{
		Carp::croak($sub." failed: STOP codons in some, but not all sequences in \'".$file."\', exiting;");
	}
	
	$self->update_mappings();	
	
	#Check for all-gaps sites.	
	my $codon_count = $self->get_codon_count();
	my $sequence_count = $self->get_sequence_count();
	for(my $c=0, my $s=1 ; $c<$codon_count ; $c++, $s++)
	{
		my @alignment_site = map { substr($self->{'SEQS'}{$_}, ($c*3), 3) } keys(%{$self->{'SEQS'}});
				
		if( (grep { $_ eq "---" } @alignment_site) == $sequence_count )
		{Carp::croak($sub." failed: no sequence in site ".$s." of alignment in \'".$file."\', exiting;");}
	}
	
	if($sequence_count < MIN_SEQUENCES)
	{Carp::croak($sub." failed: number of sequences (".$sequence_count.") in \'".$file."\' is below minimum number (".MIN_SEQUENCES."), exiting;");}
	
	foreach my $seqname ($self->get_sequence_names())
	{
		if($seqname=~/^(FOREGROUND|ALIGNMENT)$/i)
		{Carp::croak($sub." failed: can't accept \'".$&."\' as a sequence name (clashes with \'Codeml::Alignment\' keyword), exiting;");}
	}
	
	print($sub.": Codeml alignment with ".$self->get_sequence_count()." sequences of ".$self->get_codon_count()." codons loaded from ".$self->{'INPUT-FORMAT'}." file \'".$file."\'...\n" ) if($self->{'VERBOSE'});	
}

=head2 save

 Usage   : $obj->save($file);
           $obj->save($file, $format);
 Function: Saves CodemlWrapper::Alignment object to the specified file. If a supported 
           format is specified, saves the file in that format.
 Returns : none
 Args    : $file (file name), $format (optional string)

=cut

sub save
{
	my ($self, $file, $output_format) = @_;
	my $sub = $self->get_nickname().substr((caller(0))[3], MODULE_NAME_LENGTH);
	undef my $output_string;
	
	if(!defined($file))
	{Carp::croak($sub." failed: no file name given, exiting;");}
		
	if(!defined($self->{'SEQS'}))
	{Carp::croak($sub." failed: no alignment defined, exiting;");}
	
	$self->set_output_format($output_format) if(defined($output_format));
	
	my @entry_indices = sort { $a <=> $b } keys(%{$self->{'INDX2HDR'}});
	my @entry_names   = map { $self->{'INDX2HDR'}{$_} } @entry_indices;
	
	if($self->{'OUTPUT-FORMAT'} eq 'FASTA')
	{
		foreach my $seqname (@entry_names)
		{
			my $sequence = $self->get_sequence($seqname, 'FASTA');	
			
			$output_string .= ">".$seqname."\n".$sequence;
		}	
	}
	elsif($self->{'OUTPUT-FORMAT'} eq 'PHYLIP-INTERLEAVED')
	{
		my $sequence_tally = @entry_names;
		my $sequence_lengths = length($self->{'SEQS'}{$entry_names[0]});
		my $sn_width = util_get_header_width(@entry_names);
	
		$output_string = " ".@entry_names." ".$sequence_lengths." I\n";
	
		undef my $seq_lines;
		my %split_sequences = ();		
		foreach my $seqname (@entry_names)
		{
			my $formatted_sequence = $self->get_sequence($seqname, 'PHYLIP-INTERLEAVED');
			my @sequence_array = split(/\n/, $formatted_sequence);	
			$split_sequences{$seqname} =  [@sequence_array];
			
			unless(defined($seq_lines))
			{$seq_lines = @sequence_array;}
		}
		
		for(my $lineI=0 ; $lineI<$seq_lines ; $lineI++)
		{
			foreach my $seqname (@entry_names)
			{
				if($lineI==0)
				{$output_string .= $seqname.(" " x ($sn_width - length($seqname)));}
				else
				{$output_string .= (" " x $sn_width);}
				
				$output_string .= $split_sequences{$seqname}[$lineI]."\n";
			}
			
			$output_string .= "\n";
		}				
	}
	elsif($self->{'OUTPUT-FORMAT'} eq 'PHYLIP-SEQUENTIAL')
	{
		my $sequence_tally = @entry_names;
		my $sequence_lengths = length($self->{'SEQS'}{$entry_names[0]});
		my $sn_width = util_get_header_width(@entry_names);
		
		$output_string = " ".@entry_names." ".$sequence_lengths."\n";
		
		undef my $seq_lines;
		my %split_sequences = ();		
		foreach my $seqname (@entry_names)
		{
			my $formatted_sequence = $self->get_sequence($seqname, 'PHYLIP-SEQUENTIAL');
			my @sequence_array = split( /\n/, $formatted_sequence);	
			$split_sequences{$seqname} =  [@sequence_array];
			
			unless(defined($seq_lines))
			{$seq_lines = @sequence_array;}
		}
	
		foreach my $seqname (@entry_names)
		{	
			for(my $lineS=0 ; $lineS<$seq_lines ; $lineS++)
			{
				
				if($lineS==0)
				{$output_string .= $seqname.(" " x ($sn_width - length($seqname)));}
				else
				{$output_string .= (" " x $sn_width);}
				
				$output_string .= $split_sequences{$seqname}[$lineS]."\n";
			}
			
			$output_string .= "\n";				
		}
	}
	elsif($self->{'OUTPUT-FORMAT'} eq 'HTML')
	{
		my $sequence_tally = @entry_names;
		my $sequence_lengths = length($self->{'SEQS'}{$entry_names[0]});
		my $sn_width = util_get_header_width(@entry_names);
	
		$output_string = join("\n", @html_boilerplate_start)."\n";
	
		undef my $seq_lines;
		undef my %split_sequences;		
		foreach my $seqname (@entry_names)
		{
			my $formatted_sequence = $self->get_sequence($seqname, 'HTML');
			my @sequence_array = split( /\n/, $formatted_sequence);	
			$split_sequences{$seqname} =  [@sequence_array];
			
			unless(defined($seq_lines))
			{$seq_lines = @sequence_array;}
		}
		
		for(my $lineI=0 ; $lineI<$seq_lines ; $lineI++)
		{
			foreach my $seqname (@entry_names)
			{
				$output_string .= $seqname.(" " x ($sn_width - length($seqname))).$split_sequences{$seqname}[$lineI]."\n";
			}
			
			$output_string .= "\n";
		}		
		
		$output_string .= join("\n", @html_boilerplate_end);
	}
	else
	{
		Carp::croak($sub." failed: no output format specified, exiting;");
	}
	
	util_write_file($file, $output_string);
	
	print($sub.": Codeml alignment with ".$self->get_sequence_count()." sequences of ".$self->get_codon_count()." codons saved to ".$self->{'OUTPUT-FORMAT'}." file  \'".$file."\'...\n" ) if($self->{'VERBOSE'});
}

=head2 set_foreground_sequences

 Usage   : $obj->set_foreground_sequences(@sequence_names);
 Function: Sets the sequences that are in the foreground of an associated 
           CodemlWrapper::Tree object, so that any positive sites in the alignment 
           will be highlighted in a different colour in foreground and background
           sequences. 
 Returns : none
 Args    : @sequence_names (array of strings)

=cut

sub set_foreground_sequences
{
	my($self, @foreground_seqnames) = @_;
	my $sub = $self->get_nickname().substr((caller(0))[3], MODULE_NAME_LENGTH);
	
	if(!defined($self->{'SEQS'}))
	{Carp::croak($sub." failed: no alignment defined, exiting;");}

	if( (grep { defined($self->{'SEQS'}{$_}) } @foreground_seqnames) != @foreground_seqnames )
	{Carp::croak($sub." failed: one or more sequences not in alignment, exiting;");}

	if(@foreground_seqnames > 0)
	{	
		if(@foreground_seqnames >= $self->get_sequence_count())
		{Carp::croak($sub." failed: too many foreground sequences specified, exiting;");}
	
		my %fgsn = map {$_ => 1} @foreground_seqnames;
		foreach my $name (keys(%{$self->{'SEQS'}}))
		{
			$self->{'FGS'}{$name} = defined($fgsn{$name}) ? 'F' : 'B';
		}
	}
	else
	{
		delete $self->{'FGS'};		
	}
	
	
	print($sub.": ".@foreground_seqnames." foreground sequences set...\n" ) if($self->{'VERBOSE'});
}

=head2 get_foreground_sequences

 Usage   : @sequence_names = $obj->get_foreground_sequences();
 Function: Gets a list of names of sequences that are considered to be in the foreground. 
 Returns : @sequence_names (array of strings)
 Args    : none

=cut

sub get_foreground_sequences
{
	my($self) = @_;
	my $sub = $self->get_nickname().substr((caller(0))[3], MODULE_NAME_LENGTH);
	
	undef my @foreground_seqnames;
	
	if(!defined($self->{'SEQS'}))
	{Carp::croak($sub." failed: no alignment defined, exiting;");}
	
	if(defined($self->{'FGS'}))
	{
		my @seqnames = sort { $self->{'HDR2INDX'}{$a} <=> $self->{'HDR2INDX'}{$b} } keys(%{$self->{'FGS'}});
		@foreground_seqnames = grep { $self->{'FGS'}{$_} eq 'F' } @seqnames;
	}

	print($sub.": ".@foreground_seqnames." foreground sequences retrieved...\n" ) if($self->{'VERBOSE'});

	return @foreground_seqnames;
}

=head2 get_sequence_names

 Usage   : @sequence_names = $obj->get_sequence_names();
 Function: Gets a list of all sequence names in the alignment, ordered as they 
           were in the original alignment file.
 Returns : @sequence_names (array of strings)
 Args    : none

=cut

sub get_sequence_names
{
	my($self) = @_;
	my $sub = $self->get_nickname().substr((caller(0))[3], MODULE_NAME_LENGTH);
	
	if(!defined($self->{'SEQS'}))
	{Carp::croak($sub." failed: no alignment defined, exiting;");}
	
	my @sequence_names = sort { $self->{'HDR2INDX'}{$a} <=> $self->{'HDR2INDX'}{$b} } keys(%{$self->{'SEQS'}});
	
	print($sub.": ".@sequence_names." sequence names retrieved...\n" ) if($self->{'VERBOSE'} && ( !defined((caller(1))[3]) || (substr((caller(1))[3], 0, 11) ne 'Alignment::') ) );
	
	return @sequence_names;
}

=head2 get_sequence_count

 Usage   : $sequence_count = $obj->get_sequence_count();
 Function: Gets the number of sequences in the alignment.
 Returns : $sequence_count (number)
 Args    : none

=cut

sub get_sequence_count
{
	my($self) = @_;
	my $sub = $self->get_nickname().substr((caller(0))[3], MODULE_NAME_LENGTH);
	
	if(!defined($self->{'SEQS'}))
	{Carp::croak($sub." failed: no alignment defined, exiting;");}
	
	my $sequence_count = keys(%{$self->{'SEQS'}});
	
	print($sub.": ".$sequence_count." sequences in alignment...\n" ) if($self->{'VERBOSE'} && ( !defined((caller(1))[3]) || (substr((caller(1))[3], 0, 11) ne 'Alignment::') ) );
	
	return $sequence_count;
}

=head2 get_codon_count

 Usage   : $codon_count = $obj->get_codon_count();
 Function: Gets the number of codon sites in the alignment.
 Returns : $codon_count (number)
 Args    : none

=cut

sub get_codon_count
{
	my($self) = @_;
	my $sub = $self->get_nickname().substr((caller(0))[3], MODULE_NAME_LENGTH);
	
	if(!defined($self->{'SEQS'}))
	{Carp::croak($sub." failed: no alignment defined, exiting;");}
	
	my @keys = keys(%{$self->{'SEQS'}});
	
	my $codon_count = length($self->{'SEQS'}{$keys[0]})/3;
	
	print($sub.": ".$codon_count." codons in alignment...\n" ) if($self->{'VERBOSE'} && ( !defined((caller(1))[3]) || (substr((caller(1))[3], 0, 11) ne 'Alignment::') ) );
	
	return $codon_count;
}

=head2 get_sequence

 Usage   : $sequence = $obj->get_sequence($name);
           $sequence = $obj->get_sequence($name, $format);
 Function: Gets the nucleotide alignment sequence with the specified name. If a 
           format is specified, gets the sequence in that format. Possible formats
           include PHYLIP sequential ('PHYLIP-SEQUENTIAL'), PHYLIP interleaved 
           ('PHYLIP-INTERLEAVED'), FASTA format ('FASTA') and HTML ('HTML'). In 
           'HTML', the sequence will be contained in a series of '<span>' html tags. 
           In all other formats, a formatted sequence is simply composed of lines 
           of up to 60 nucleotides.
 Returns : $sequence (string)
 Args    : $name (string), $format (optional string)

=cut

sub get_sequence
{
	my($self, $name, $format) = @_;
	my $sub = $self->get_nickname().substr((caller(0))[3], MODULE_NAME_LENGTH);

	if(!defined($self->{'SEQS'}))
	{Carp::croak($sub." failed: no alignment defined, exiting;");}
	
	if(!defined($name))
	{Carp::croak($sub." failed: no sequence name given, exiting;");}
		
	if(!defined($self->{'SEQS'}{$name}))
	{Carp::croak($sub." failed: sequence \'".$name."\' not in alignment, exiting;");}
	
	my $sequence = $self->{'SEQS'}{$name};
	
	if(defined($format))
	{
		my $formatted_sequence = '';
		
		if( ($format eq 'FASTA') || ($format eq 'PHYLIP-INTERLEAVED') || ($format eq 'PHYLIP-SEQUENTIAL') )
		{
			for(my $i=0 ; $i<length($sequence) ; $i+=60)
			{
				$formatted_sequence .= substr($sequence, $i, 60)."\n";
			}
		}
		elsif($format eq 'HTML')
		{		
			for(my $i=0 ; $i<length($sequence) ; $i+=60)
			{
				my $seq_line = substr($sequence, $i, 60);
				my $format_line = substr($self->{'FORMATTING'}, $i, 60);
				
				my $j = 0;
				while($format_line=~/((MMM)+|(PPP)+|(NNN)+|(SSS)+)/g)
				{
					my $l = length($&);
					my $t = substr($&, 0, 1);
					my $sc = substr($seq_line, $j, $l);
					
					if($t eq 'P')
					{
						$t = defined($self->{'FGS'}{$name}) ? $self->{'FGS'}{$name} : 'F' ;
					}
					
					$formatted_sequence .= $span{$t} . $sc . $se;
					$j += $l;
				}
				
				$formatted_sequence .= "\n";
			}
		}
		else
		{
			Carp::croak($sub." failed: unknown format \'".$format."\', exiting;");	
		}
		
		$sequence = $formatted_sequence;
	}
	
	print($sub.": sequence \'".$name."\' retrieved".(defined($format) ? " in \'".$format."\' format" : "")."...\n" ) if($self->{'VERBOSE'} && ( !defined((caller(1))[3]) || (substr((caller(1))[3], 0, 11) ne 'Alignment::') ) );

	return $sequence;
}

=head2 get_residue

 Usage   : $residue = $obj->get_residue($name, $site);
 Function: Gets amino acid residue, gap or ambiguity character corresponding to 
           the codon in the given site in the nucleotide alignment. Site position 
           can range from 1 to N, where N is the number of codons in the alignment.
 Returns : $residue (character)
 Args    : $name (string), $site (number)

=cut

sub get_residue
{
	my($self, $name, $site) = @_;
	my $sub = $self->get_nickname().substr((caller(0))[3], MODULE_NAME_LENGTH);
	
	undef my $residue;
	
	my $codon = $self->get_codon($name, $site);
	
	#If standard codon, translate it..
	if(defined($standard_codons{$codon}))			
	{
		$residue = $standard_codons{$codon};
	}
	#..otherwise, if translatable redundant codon, translate it..
	elsif(defined($redundant_codons{$codon}))		
	{
		$residue = $redundant_codons{$codon};
	}
	#..otherwise, just translate to ambiguity 'X'. 
	else											
	{
		$residue = 'X';
	}					

	print($sub.": residue \'".$residue."\' retrieved from site ".$site." of sequence \'".$name."\'...\n" ) if($self->{'VERBOSE'});
	
	return $residue;
}

=head2 get_codon

 Usage   : $codon = $obj->get_codon($name, $site);
 Function: Gets triplet of nucleotides, gaps or ambiguity characters in the given 
           site in the nucleotide alignment. Site position can range from 1 to N, 
           where N is the number of codons in the alignment.
 Returns : $codon (character)
 Args    : $name (string), $site (number)

=cut

sub get_codon
{
	my($self, $name, $site) = @_;
	my $sub = $self->get_nickname().substr((caller(0))[3], MODULE_NAME_LENGTH);
	
	if(!defined($self->{'SEQS'}))
	{Carp::croak($sub." failed: no alignment defined, exiting;");}

	if(!defined($name))
	{Carp::croak($sub." failed: sequence name not specified, exiting;");}	
	
	if(!defined($self->{'SEQS'}{$name}))
	{Carp::croak($sub." failed: sequence \'".$name."\' not in alignment, exiting;");}
	
	if(!defined($site))
	{Carp::croak($sub." failed: sequence site not specified, exiting;");}	
	
	if( ($site!~/^\d+$/) || ($site < 1) || ($site > $self->{'LEN'}{'ALIGNMENT'}{$name}) )
	{Carp::croak($sub." failed: invalid sequence site: \'".$site."\', exiting;");}	
	
	my $sequence = $self->{'SEQS'}{$name};
	my $codon = substr($sequence, ($site-1)*3, 3);

	print($sub.": codon \'".$codon."\' retrieved from site ".$site." of sequence \'".$name."\'...\n" ) if($self->{'VERBOSE'} && ( !defined((caller(1))[3]) || (substr((caller(1))[3], 0, 11) ne 'Alignment::') ) );
	
	return $codon;
}

=head2 map_sites

 Usage   : $mapped_site = $obj->map_sites($source, $destination, $name, $source_site);
           @mapped_sites = $obj->map_sites($source, $destination, $name, @source_sites);
 Function: For a specified sequence, maps site position(s) from one form of the  
           sequence to another. Each sequence has three different forms: 'ALIGNMENT' 
           (the full alignment sequence), 'UNGAPPED' (the alignment sequence without 
           gaps) and 'UNMASKED' (the alignment sequence without masked sites). 
           If a sequence has no gaps and no sites are masked, these three forms 
           are identical. Sites can be mapped from 'ALIGNMENT' to 'UNGAPPED' and 
           vice versa, as well as from 'ALIGNMENT' to 'UNMASKED' and vice versa.
           Input sites are accepted in the range 1 to S (where S is the number
           of sites in the source sequence form), and given in the range 1 to D
           (where D is the number of sites in the destination sequence form).
           Note that gap sites in the 'ALIGNMENT' sequence are not mapped to the 
           'UNGAPPED' sequence, nor are masked sites in the 'ALIGNMENT' sequence
           mapped to the 'UNMASKED' sequence. 
 Returns : $mapped_site (number) OR @mapped_sites (array of numbers)
 Args    : $source (string), $destination (string), $name (string), $source_site (number) OR
           $source (string), $destination (string), $name (string), @source_sites (array of numbers)

=cut

sub map_sites
{
	my($self, $source, $destination, $name, @source_sites) = @_;
	my $sub = $self->get_nickname().substr((caller(0))[3], MODULE_NAME_LENGTH);
	
	if(!defined($self->{'SEQS'}))
	{Carp::croak($sub." failed: no alignment defined, exiting;");}
	
	if( !defined($source) )
	{Carp::croak($sub." failed: no mapping source specified, exiting;");}
	
	if( !defined($destination) )
	{Carp::croak($sub." failed: no mapping destination specified, exiting;");}
	
	if( !defined($self->{'MAP'}{$source}{$destination}) )
	{Carp::croak($sub." failed: no mapping from \'".$source."\' to \'".$destination."\', exiting;");}
	
	if(!defined($name))
	{Carp::croak($sub." failed: alignment sequence name not specified, exiting;");}	
	
	if(!defined($self->{'SEQS'}{$name}))
	{Carp::croak($sub." failed: alignment sequence \'".$name."\' not in alignment, exiting;");}
	
	if(@source_sites==0)
	{Carp::croak($sub." failed: alignment sites not specified, exiting;");}	
	
	if(!wantarray && @source_sites>1)
	{Carp::croak($sub." failed: can't map multiple sites in scalar context, exiting;");}
	
	my @numerical_sites = grep { $_=~/^\d+$/} @source_sites;
	if(@source_sites!=@numerical_sites)
	{Carp::croak($sub." failed: non-numerical alignment site(s), exiting;");}	
	
	my @sites_within_range = grep { ($_ >= 1) && ($_ <= $self->{'LEN'}{$source}{$name}) } @source_sites;
	if(@sites_within_range != @source_sites)
	{Carp::croak($sub." failed: site(s) out of range, exiting;");}	

	undef my @sites;
	foreach my $s (@source_sites)
	{
		my $site = $self->{'MAP'}{$source}{$destination}{$name}{$s};
		push(@sites, $site) if(defined($site));
	}

	print($sub.": ".@source_sites." \'".$source."\' sites mapped to ".@sites." \'".$source."\' sites of sequence \'".$name."\'...\n" ) if($self->{'VERBOSE'});
	
	if(wantarray)
	{
		return @sites;
	}
	else
	{
		return $sites[0];
	}
}

=head2 set_pos_sites

 Usage   : $obj->set_pos_sites(@sites);
 Function: Sets pos sites in the alignment. When the alignment is saved as a HTML
           file, these sites are highlighted in red (or red and blue, if foreground
           sequences were set), while all other unmasked sites are given a grey 
           background. Setting pos sites clears any previous pos sites. Attempting 
           to set pos sites in a masked region will have no effect on the alignment
           and will produce a warning.
 Returns : none
 Args    : @sites (array of numbers)

=cut

sub set_pos_sites
{
	my($self, @sites) = @_;
	my $sub = $self->get_nickname().substr((caller(0))[3], MODULE_NAME_LENGTH);

	my $codon_count = $self->get_codon_count();
	my @indices = map { ($_-1)*3 } @sites;
	
	if(!defined($self->{'SEQS'}))
	{Carp::croak($sub." failed: no alignment defined, exiting;");}
	
	if(@sites==0)
	{Carp::croak($sub." failed: sites not specified, exiting;");}	
	
	my @numerical_sites = grep { $_=~/^\d+$/} @sites;
	if(@sites != @numerical_sites)
	{Carp::croak($sub." failed: non-numerical site(s), exiting;");}	
	
	my @sites_within_range = grep { ($_ >= 1) && ($_ <= $codon_count) } @sites;
	if(@sites_within_range != @sites)
	{Carp::croak($sub." failed: site(s) out of range, exiting;");}
	
	if($self->{'FORMATTING'}=~/PPP/)
	{
		$self->{'FORMATTING'}=~s/NNN/SSS/g;
		$self->{'FORMATTING'}=~s/PPP/SSS/g;
	}
	
	my $set_count = 0;
	foreach my $i (@indices)
	{
		my $codon = substr($self->{'FORMATTING'}, $i, 3);
		
		if($codon ne 'MMM')
		{
			substr($self->{'FORMATTING'}, $i, 3) = 'PPP';
			$set_count++;
		}
		else
		{
			Carp::carp($sub." can't set pos site (".$i.") in masked region...");
		}
		
	}
	
	$self->{'FORMATTING'}=~s/SSS/NNN/g;
	
	print($sub.": ".$set_count." of ".@sites." pos sites set...\n" ) if($self->{'VERBOSE'});
}

=head2 get_pos_sites

 Usage   : @pos_site_positions = $obj->get_pos_sites();
 Function: Gets a list of pos sites in the alignment. Pos site positions range 
           from 1 to N, where N is the alignment length.
 Returns : @pos_site_positions (array of numbers)
 Args    : none

=cut

sub get_pos_sites
{
	my($self) = @_;
	my $sub = $self->get_nickname().substr((caller(0))[3], MODULE_NAME_LENGTH);

	undef my @pos_site_positions;
	
	if(!defined($self->{'SEQS'}))
	{Carp::croak($sub." failed: no alignment defined, exiting;");}
	
	while($self->{'FORMATTING'}=~/PPP/g)
	{
		my $l = length($`);
		my $p = ($l/3) + 1;
		push(@pos_site_positions, $p);
	}
	
	print($sub.": ".@pos_site_positions." pos sites retrieved...\n" ) if($self->{'VERBOSE'});
	
	return @pos_site_positions;
}

=head2 set_mask_sites

 Usage   : $obj->set_mask_sites(@sites);
 Function: Sets mask sites in the alignment. When the alignment is saved as a HTML
           file, these sites are masked in a black background. This is useful for
           marking regions of an alignment that are deleted before further analysis
           is done. Setting mask sites clears any previous mask sites.
 Returns : none
 Args    : @sites (array of numbers)

=cut

sub set_mask_sites
{
	my($self, @sites) = @_;
	my $sub = $self->get_nickname().substr((caller(0))[3], MODULE_NAME_LENGTH);	

	my @indices = map { ($_-1)*3 } @sites;
	
	if(!defined($self->{'SEQS'}))
	{Carp::croak($sub." failed: no alignment defined, exiting;");}
	
	if(@sites==0)
	{Carp::croak($sub." failed: sites not specified, exiting;");}	
	
	my @numerical_sites = grep { $_=~/^\d+$/} @sites;
	if(@sites != @numerical_sites)
	{Carp::croak($sub." failed: non-numerical site(s), exiting;");}	
	
	my @sites_within_range = grep { ($_ >= 1) && ($_ <= $self->get_codon_count()) } @sites;
	if(@sites_within_range != @sites)
	{Carp::croak($sub." failed: site(s) out of range, exiting;");}

	if($self->{'FORMATTING'}=~/MMM/)
	{
		my $rep = ($self->{'FORMATTING'}=~/PPP/) ? 'NNN' : 'SSS';
		$self->{'FORMATTING'}=~s/MMM/$rep/g;
	}

	if(@sites == $self->get_codon_count())
	{Carp::croak($sub." failed: masking ".$self->get_codon_count()." sites leaves no alignment, exiting;");}	
		
	my $pos_sites_before = () = $self->{'FORMATTING'}=~/PPP/g;
	
	foreach my $i (@indices)
	{
		substr($self->{'FORMATTING'}, $i, 3) = 'MMM';
	}
	
	my $pos_sites_after = () = $self->{'FORMATTING'}=~/PPP/g;
	
	if( ($pos_sites_before > 0) && ($pos_sites_after == 0) )
	{$self->{'FORMATTING'}=~s/NNN/SSS/g;}
	
	$self->update_mappings();

	foreach my $name (keys(%{$self->{'SEQS'}}))
	{
		my @unmasked_codons = grep { (substr($self->{'FORMATTING'}, ($_-1)*3, 3) ne 'MMM') && (substr($self->{'SEQS'}{$name}, ($_-1)*3, 3) ne '---') } 1..$self->get_codon_count();
		
		if(@unmasked_codons == 0)
		{Carp::croak($sub." failed: masking ".@sites." sites leaves no sequence in \'".$name."\', exiting;");}
	}
	
	if($self->{'VERBOSE'})
	{
		my $verbose_string = ($pos_sites_after < $pos_sites_before) ? " (".($pos_sites_before - $pos_sites_after)." pos sites masked)" : "";
		print($sub.": ".@sites." mask sites set".$verbose_string."...\n" ) ;
	}
}

=head2 get_mask_sites

 Usage   : @mask_site_positions = $obj->get_mask_sites();
 Function: Gets a list of mask sites in the alignment. Mask site positions range 
           from 1 to N, where N is the alignment length.
 Returns : @mask_site_positions (array of numbers)
 Args    : none

=cut

sub get_mask_sites
{
	my($self) = @_;
	my $sub = $self->get_nickname().substr((caller(0))[3], MODULE_NAME_LENGTH);

	undef my @mask_site_positions;
	
	if(!defined($self->{'SEQS'}))
	{Carp::croak($sub." failed: no alignment defined, exiting;");}
	
	while($self->{'FORMATTING'}=~/MMM/g)
	{
		my $l = length($`);
		my $m = ($l/3) + 1;
		push(@mask_site_positions, $m);
	}
	
	print($sub.": ".@mask_site_positions." mask sites retrieved...\n" ) if($self->{'VERBOSE'});
	
	return @mask_site_positions;
}

=head2 get_ambig_sites

 Usage   : @ambig_site_positions = $obj->get_ambig_sites($ambig_tolerance);
 Function: Gets a list of sites for which the number of sequences with ambig 
           characters at that site exceeds the specified 'ambig tolerance'. 
           Ambig tolerance can range from 0 to S-2, where S is the number of 
           sequences. If no ambig tolerance is given, an ambig tolerance of 
           zero is used. Ambig site positions range from 1 to N, where N is 
           the alignment length.
 Returns : @ambig_site_positions (array of numbers)
 Args    : $ambig_tolerance (optional number)

=cut

sub get_ambig_sites
{
	my($self, $ambig_tolerance) = @_;
	my $sub = $self->get_nickname().substr((caller(0))[3], MODULE_NAME_LENGTH);
	
	$ambig_tolerance = 0 if(!defined($ambig_tolerance));
	undef my @ambig_site_positions;
	
	if(!defined($self->{'SEQS'}))
	{Carp::croak($sub." failed: no alignment defined, exiting;");}

	if($ambig_tolerance < 0)
	{Carp::croak($sub." failed: ambig tolerance can't be less than zero, exiting;");}
	
	if($ambig_tolerance > ($self->get_sequence_count()-2))
	{Carp::croak($sub." failed: ambig tolerance must be at least two less than number of sequences in the alignment, exiting;");}
	
	for(my $c=0, my $s=1 ; $c < $self->get_codon_count() ; $c++, $s++)
	{
		my @alignment_site = map { substr($self->{'SEQS'}{$_}, $c*3, 3) } keys(%{$self->{'SEQS'}});
		
		if( (grep { $_=~/[BDHKMNRSVWY\?]/ } @alignment_site) > $ambig_tolerance )
		{
			push(@ambig_site_positions, $s);
		}
	}

	print($sub.": ".@ambig_site_positions." sites found with ambiguous codons in ".($ambig_tolerance + 1)." or more sequences...\n" ) if($self->{'VERBOSE'});
	
	return @ambig_site_positions;
}

=head2 get_gap_sites

 Usage   : @gap_site_positions = $obj->get_gap_sites($gap_tolerance);
 Function: Gets a list of sites for which the number of sequences with gaps at 
           that site exceeds the specified 'gap tolerance'. Gap tolerance can  
           range from 0 to S-2, where S is the number of sequences. If no gap
           tolerance is given, a gap tolerance of zero is used. Gap site  
           positions range from 1 to N, where N is the alignment length.
 Returns : @gap_site_positions (array of numbers)
 Args    : $gap_tolerance (optional number)

=cut

sub get_gap_sites
{
	my($self, $gap_tolerance) = @_;
	my $sub = $self->get_nickname().substr((caller(0))[3], MODULE_NAME_LENGTH);

	$gap_tolerance = 0 if(!defined($gap_tolerance));	
	undef my @gap_site_positions;
	
	if(!defined($self->{'SEQS'}))
	{Carp::croak($sub." failed: no alignment defined, exiting;");}
	
	if($gap_tolerance < 0)
	{Carp::croak($sub." failed: gap tolerance can't be less than zero, exiting;");}
	
	if($gap_tolerance > ($self->get_sequence_count()-2))
	{Carp::croak($sub." failed: gap tolerance must be at least two less than number of sequences in the alignment, exiting;");}
	
	for(my $c=0, my $s=1 ; $c < $self->get_codon_count() ; $c++, $s++)
	{
		my @alignment_site = map { substr($self->{'SEQS'}{$_}, $c*3, 3) } keys(%{$self->{'SEQS'}});
					
		if( (grep { $_ eq "---" } @alignment_site) > $gap_tolerance )
		{
			push(@gap_site_positions, $s);
		}
	}

	print($sub.": ".@gap_site_positions." sites found with gap codons in ".($gap_tolerance + 1)." or more sequences...\n" ) if($self->{'VERBOSE'});
	
	return @gap_site_positions;
}

=head2 remove_sites

 Usage   : $obj->remove_sites(@sites);
 Function: Removes specified sites from the alignment. Site positions are accepted
           in the range 1 to N, where N is the length of the original alignment.
 Returns : none
 Args    : @sites (array of numbers)

=cut

sub remove_sites
{
	my($self, @sites) = @_;
	my $sub = $self->get_nickname().substr((caller(0))[3], MODULE_NAME_LENGTH);
	
	my @indices = map { ($_-1)*3 } @sites;
	
	if(!defined($self->{'SEQS'}))
	{Carp::croak($sub." failed: no alignment defined, exiting;");}
	
	if(@sites==0)
	{Carp::croak($sub." failed: sites not specified, exiting;");}	
	
	my @numerical_sites = grep { $_=~/^\d+$/} @sites;
	if(@sites != @numerical_sites)
	{Carp::croak($sub." failed: non-numerical site(s), exiting;");}	
	
	my @sites_within_range = grep { ($_ >= 1) && ($_ <= $self->get_codon_count()) } @sites;
	if(@sites_within_range != @sites)
	{Carp::croak($sub." failed: site(s) out of range, exiting;");}
	
	if(@sites == $self->get_codon_count())
	{Carp::croak($sub." failed: removing ".@sites." sites leaves no alignment, exiting;");}	
		
	foreach my $site_index (reverse(@indices))
	{
		foreach my $seq (keys(%{$self->{'SEQS'}}))
		{
			substr($self->{'SEQS'}{$seq}, $site_index, 3) = '';
		}
		
		substr($self->{'FORMATTING'}, $site_index, 3) = '';
	}
	
	if( ($self->{'FORMATTING'}=~/NNN/) && ($self->{'FORMATTING'}!~/PPP/) )
	{$self->{'FORMATTING'}=~s/NNN/SSS/g;}
	
	$self->update_mappings();

	foreach my $name (keys(%{$self->{'SEQS'}}))
	{
		my @unmasked_codons = grep { (substr($self->{'FORMATTING'}, ($_-1)*3, 3) ne 'MMM') && (substr($self->{'SEQS'}{$name}, ($_-1)*3, 3) ne '---') } 1..$self->get_codon_count();
		
		if(@unmasked_codons == 0)
		{Carp::croak($sub." failed: removing ".@sites." sites leaves no sequence in \'".$name."\', exiting;");}
	}
	
	print($sub.": ".@sites." sites removed, leaving ".$self->get_codon_count()." sites...\n" ) if($self->{'VERBOSE'});
}

=head2 get_file_string

 Usage   : $text = $obj->get_file_string();
 Function: Gets text of alignment file from which this CodemlWrapper::Alignment object was loaded.
 Returns : $text (string)
 Args    : none

=cut

sub get_file_string
{
	my($self) = @_;
	my $sub = $self->get_nickname().substr((caller(0))[3], MODULE_NAME_LENGTH);
			
	if(!defined($self->{'SEQS'}))
	{Carp::croak($sub." failed: no alignment defined, exiting;");}
	
	print($sub.":...\n" ) if($self->{'VERBOSE'});
	
	return $self->{'FILE-STRING'};
}

=head2 alignment_match

 Usage   : if($obj->alignment_match($other)) {...}
 Function: Tests if another CodemlWrapper::Alignment object is identical to this one. 
 Returns : Boolean indicating if other CodemlWrapper::Alignment object is identical to this one
 Args    : $other (CodemlWrapper::Alignment object)

=cut

sub alignment_match
{
	my($self, $other) = @_;
	my $sub = $self->get_nickname().substr((caller(0))[3], MODULE_NAME_LENGTH);
	
	my $identity = TRUE;
	
	if(!defined($self->{'SEQS'}))
	{Carp::croak($sub." failed: no alignment defined, exiting;");}
	
	if(!defined($other))
	{Carp::croak($sub." failed: no other \'CodemlWrapper::Alignment\' object given, exiting;");}
	
	if(ref($other) ne "Alignment")
	{Carp::croak($sub." failed: can only compare to a \'CodemlWrapper::Alignment\' object, exiting;");}
	
	if(!defined($other->{'SEQS'}))
	{Carp::croak($sub." failed: other alignment not defined, exiting;");}
	
	my @self_headers  = $self->get_sequence_names();
	my @other_headers = $other->get_sequence_names();
	
	if(@self_headers==@other_headers)
	{
		foreach my $h (@self_headers)
		{
			my $self_seq  = $self->get_sequence($h);
			my $other_seq = $other->get_sequence($h);
			
			if( ($self_seq ne $other_seq) || (!$self_seq) || (!$other_seq) )
			{
				$identity = FALSE;
			}
		}
	}
	else
	{
		$identity = FALSE;
	}
	
	print($sub.": ".($identity ? "TRUE" : "FALSE")."...\n" ) if($self->{'VERBOSE'});
	
	return $identity;
}

=head2 set_output_format

 Usage   : $obj->set_output_format($format);
 Function: Sets the alignment format used for file output. The supported alignment
           output formats are PHYLIP sequential ('PHYLIP-SEQUENTIAL'), PHYLIP interleaved 
           ('PHYLIP-INTERLEAVED'), FASTA format ('FASTA') and HTML ('HTML').
 Returns : none
 Args    : $format (string)

=cut

sub set_output_format
{	
	my($self, $output_format) = @_;
	my $sub = $self->get_nickname().substr((caller(0))[3], MODULE_NAME_LENGTH);
	
	$output_format = uc($output_format);
	
	unless( grep { $_ eq $output_format } @supported_file_types )
	{Carp::croak($sub." failed: \'".$output_format."\' is not a supported file type, exiting;");}
	
	$self->{'OUTPUT-FORMAT'} = $output_format;	
	
	print($sub.": output format set to \'".$output_format."\'...\n" ) if($self->{'VERBOSE'} && ( !defined((caller(1))[3]) || (substr((caller(1))[3], 0, 11) ne 'Alignment::') ) );
}

=head2 get_output_format

 Usage   : $format = $obj->get_output_format();
 Function: Gets the alignment format used for file output.
 Returns : $format (string)
 Args    : none

=cut

sub get_output_format
{	
	my($self) = @_;
	my $sub = $self->get_nickname().substr((caller(0))[3], MODULE_NAME_LENGTH);

	print($sub.": output format (\'".$self->{'OUTPUT-FORMAT'}."\') retrieved...\n" ) if($self->{'VERBOSE'});

	return $self->{'OUTPUT-FORMAT'};		
}

sub DESTROY
{
	my( $self ) = @_;
			
	undef %$self;
}

sub get_nickname
{
	my ($self) = @_;

	#INTERNAL USE ONLY!
	if( !defined((caller(1))[3]) || (substr((caller(1))[3], 0, 11) ne 'Alignment::') )
	{Carp::croak((caller(0))[3]." failed: can only be called by \'CodemlWrapper::Alignment\' subroutines, exiting;");}
	
	my $nickname = 'Alignment';

	if( defined($self->{'SEQS'}) )
	{
		my @keys = keys(%{$self->{'SEQS'}});
		my $sequence_count = scalar(@keys);
		my $codon_count = length($self->{'SEQS'}{$keys[0]})/3;
		$nickname .= "(".$sequence_count."x".$codon_count.")";
	}

	return $nickname;
}

sub update_mappings
{
	my ($self) = @_;
	
	#INTERNAL USE ONLY!
	if( !defined((caller(1))[3]) || (caller(1))[3]!~/^Alignment::([^:]+)$/ )
	{Carp::croak((caller(0))[3]." failed: can only be called by \'CodemlWrapper::Alignment\' subroutines, exiting;");}

	if(!defined($self->{'SEQS'}))
	{Carp::croak((caller(1))[3]."\' failed: no alignment defined, exiting;");}		

	foreach my $name (keys(%{$self->{"SEQS"}}))
	{
		my $seq = $self->{'FORMATTING'};
		my $len = length($seq)/3;
		
		for(my $x=1, my $y=1; $x<=$len ; $x++)
		{
			my $codon = substr($seq, ($x-1)*3, 3);
			
			if( $codon ne 'MMM' )
			{
				$self->{'MAP'}{'ALIGNMENT'}{'UNMASKED'}{$name}{$x} = $y;
				$self->{'MAP'}{'UNMASKED'}{'ALIGNMENT'}{$name}{$y} = $x;
				$y++;
			}
			
			if($x == $len)
			{
				$self->{'LEN'}{'ALIGNMENT'}{$name} = $x if( !defined($self->{'LEN'}{'ALIGNMENT'}{$name}) );
				$self->{'LEN'}{'UNMASKED'}{$name} = $y if( !defined($self->{'LEN'}{'UNMASKED'}{$name}) );
			}
		}
		

	}	
		
	foreach my $name (keys(%{$self->{"SEQS"}}))
	{
		my $seq = $self->{'SEQS'}{$name};
		my $len = length($seq)/3;
		
		for(my $x=1, my $y=1; $x<=$len ; $x++)
		{
			my $codon = substr($seq, ($x-1)*3, 3);
			
			if( $codon ne "---" )
			{
				$self->{'MAP'}{'ALIGNMENT'}{'UNGAPPED'}{$name}{$x} = $y;
				$self->{'MAP'}{'UNGAPPED'}{'ALIGNMENT'}{$name}{$y} = $x;
				$y++;
			}
			
			if($x == $len)
			{
				$self->{'LEN'}{'ALIGNMENT'}{$name} = $x if( !defined($self->{'LEN'}{'ALIGNMENT'}{$name}) );
				$self->{'LEN'}{'UNGAPPED'}{$name} = $y if( !defined($self->{'LEN'}{'UNGAPPED'}{$name}) );
			}
		}
		

	}
}

sub util_trim_sequence
{
	my($sequence) = @_;
	
	#INTERNAL USE ONLY!
	if( !defined((caller(1))[3]) || (caller(1))[3]!~/^Alignment::([^:]+)$/ )
	{Carp::croak((caller(0))[3]." failed: can only be called by \'CodemlWrapper::Alignment\' subroutines, exiting;");}
	
	my $trim_seq = $sequence;
	
	$trim_seq = uc($trim_seq);	
			
	$trim_seq=~s/[^ABCDGHKMNRSTUVWY\-\?\.]//g;
			
	$trim_seq=~s/U/T/g;		
	
	return $trim_seq;
}

sub util_get_header_width
{
	my(@headers) = @_;
	
	#INTERNAL USE ONLY!
	if( !defined((caller(1))[3]) || (caller(1))[3]!~/^Alignment::([^:]+)$/ )
	{Carp::croak((caller(0))[3]." failed: can only be called by \'CodemlWrapper::Alignment\' subroutines, exiting;");}
	
	my $header_width = 0;
	
	foreach my $h (@headers)
	{			
		if(length($h)>$header_width)
		{$header_width = length($h);}
	}
	$header_width += 2;
				
	return $header_width;
}

sub util_write_file
{
	my($file, @string_data) = @_;
	
	#INTERNAL USE ONLY!
	if( !defined((caller(1))[3]) || (caller(1))[3]!~/^Alignment::([^:]+)$/ )
	{Carp::croak((caller(0))[3]." failed: can only be called by \'CodemlWrapper::Alignment\' subroutines, exiting;");}
	
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
	if( !defined((caller(1))[3]) || (caller(1))[3]!~/^Alignment::([^:]+)$/ )
	{Carp::croak((caller(0))[3]." failed: can only be called by \'CodemlWrapper::Alignment\' subroutines, exiting;");}
		
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
