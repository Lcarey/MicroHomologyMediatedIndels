#!/usr/bin/perl
# convert the output of trf -ngs 
#  to a fasta file for use with blast
#  @Scer_1 230218bp [GAP=0; TRNA=4; CODING=101]
#7 28 8 2.8 8 100 0 22 36 63 0 0 0.95 CACACCCA CACACCCACACACCCACACACC CCACAC ACACCACACACCACACCACACCCACACACACACATCCTAACACTACCCTA
#20 49 12 2.5 12 100 0 30 40 60 0 0 0.97 CCACACACCACA CCACACACCACACCACACACCACACCACAC CCACACCACACCCACACAC CCACACACACACATCCTAACACTACCCTAACACAGCCCTAATCTAACCCT

# we need the header and the 15th column

use warnings; 
use strict;

my $chr = '' ; 

my $MINIMUM_ENTROPY_THRESHOLD = 1.5 ; # col 12
my $MINIMUM_REPEAT_LENGTH = 10 ; # col 2
my $MINIMUM_NUM_REPS_OF_SEQ = 2 ; # col 3
while(<>){
	chomp ; 
	if ( m/^@([\w-]+)\s/ ) {
		$chr = $1 ; 
	} elsif ( m/^@([\w-]+)$/ ) {
		$chr = $1 ; 
	} else {
		my @l = split ;
		my $seq = $l[14] . $l[15] . $l[16] ; 
		my $pos_start = $l[0] ;
		my $pos_end = $l[1] ; 
		my $entropy = $l[12] ; 
		next unless $entropy >= $MINIMUM_ENTROPY_THRESHOLD ; 
		next unless $l[2] >= $MINIMUM_REPEAT_LENGTH ; 
		next unless $l[3] >= $MINIMUM_NUM_REPS_OF_SEQ ; 
		print(">$chr.$pos_start.$pos_end\n$seq\n") ; # fasta output
	}
}

