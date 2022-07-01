#!/usr/bin/perl -w
use strict;

use Bio::Seq;
use Bio::SeqIO;

my $usage = "$0 inputfile motive outputfile\n";

# check number of arguments
if (scalar @ARGV != 3) {
        die $usage; }


my $inputfile=$ARGV[0];
my $motiv=$ARGV[1];
my $outputfile=$ARGV[2];

my $seqin  = Bio::SeqIO->new(-format => 'fasta', -file => $inputfile);
my $seqout = Bio::SeqIO->new(-format => 'fasta', -file => ">>$outputfile.fasta");

while((my $seq = $seqin->next_seq())) {
 
  if(($seq->display_id =~ /$motiv/i)||($seq->desc =~ /$motiv/i)) {
    $seqout->write_seq($seq);
  }
}

# to iterate over multiple directories, go the root/base/home directory from where your subdirs are and run : for dir in */; do  cd $dir && /path/to/fastaextract_desc_id.pl *.all.faa methanol ../methanol &&  cd ..; done
#change methanol to whatever you want, this is your search term. The second one is the name of the fasta in the base directory

