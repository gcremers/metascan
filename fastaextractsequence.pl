#!/usr/bin/perl -w
use strict;

use Bio::SeqIO;

my $usage = "$0 inputfile motive outputfilename\n";

# check number of arguments
if (scalar @ARGV != 3) {
        die $usage; }


my $inputfile=$ARGV[0];
my $motiv=$ARGV[1];
my $outputfile=$ARGV[2];

my $seqin = Bio::SeqIO->new(-file => "$inputfile", -format => "fasta");

my $seqout = Bio::SeqIO->new(-file => ">>$outputfile.fasta", -format => "fasta");

while(my $seq = $seqin->next_seq)
  { 
  if($seq->seq =~ /$motiv/i) {
    $seqout->write_seq($seq);
  }
}

# for motive and/or conserved region searches, use brackets when stating the motive. Bases can be entered as higher or lower cases. Usage: $perl /scratch/scripts_common/fastaextractsequence.pl inputfile "motive" outputfile
