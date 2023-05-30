#!/usr/bin/env perl

#    Metascan - Metabolic scan of (meta)genomes
#
#    Copyright (C) 2018- Geert Cremers
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.



use strict;
use warnings;
use Storable;
use File::Copy;
use Time::Piece;
use Time::Seconds;
use XML::Simple;
use Digest::MD5;
use List::Util qw(min max sum sum0);
use List::MoreUtils qw(uniq);
use Scalar::Util qw(openhandle);
use Data::Dumper;
use Bio::SearchIO;
use Bio::Root::Version;
use Bio::SeqIO;
use Bio::Seq;
use Bio::SeqFeature::Generic;
use Bio::Tools::GFF;
use Bio::Tools::GuessSeqFormat;
use FindBin;
use 5.010;
use Cwd qw(getcwd);
use Cwd qw(abs_path);
use lib "$FindBin::RealBin/../perl5"; # for bundled Perl modules
use File::Path qw(remove_tree);
use File::Basename;


#finding Prokka databases
my $prokpath = qx(which prokka);
chomp $prokpath; 
my $absFilePath = abs_path("$prokpath");
my ($dbloc)= split('/prokka', $absFilePath);

# Change these three paths to the location where your databases are stored: The auxillary files also go into the databasedir.
#my $databasedir="/path/to/metascan_databases";
#my $databasedir_blastnt="/path/to/blast/nt_v5";
#my $databasedir_blastn="/path/to/blast/16S_ribosomal_RNA";

my $databasedir="/vol/micro-databases/metascan";
my $databasedir_blastnt="/vol/micro-databases/blast/nt";
my $databasedir_blastn="/vol/micro-databases/blast16S/16S_ribosomal_RNA";

# When using a custom HMM profile, add: <code> #CYCLE <tab> cycle_name </code> to the first line of the hmm profile in  order to give it a cycle name    
# . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 

# list of exceptions to /product annotations
my @GOOD_PROD = (
 'methylmalonyl-CoA mutase, C-terminal domain [ec:5.4.99.2]',
 'methylmalonyl-CoA mutase, N-terminal domain [ec:5.4.99.2]',
 'Conserved virulence factor B',
 'Cypemycin N-terminal methyltransferase',
);
my %GOOD_PROD = (map { ($_=>1) } @GOOD_PROD);

# global variables
my @CMDLINE = ($0, @ARGV);
my $OPSYS = $^O;
my $BINDIR = "$FindBin::RealBin/../binaries/$OPSYS";
my $EXE = $FindBin::RealScript;
my $VERSION = "1.3";
my $AUTHOR = 'G. Cremers';
my $URL = 'gitlab.science.ru.nl/gcremers/metascan';
my $METASCAN_PMID = 'nan';
my $METASCAN_DOI = 'https://doi.org/10.3389/fbinf.2022.861505';
my $DBDIR = "$FindBin::RealBin/../db";
my $HYPO = 'hypothetical protein';
my $UNANN = 'unannotated protein';
my $MAXCONTIGIDLEN = 37;  # Genbank rule
my $SIGNALP_MAXSEQ = 10_000;  # maximum allowed input for signalp
my @LOG; # buffer up log lines before we have log file ready
my $dir = shift @ARGV // 'An empty entry';


my $KEGG_hmm="$databasedir/meta.nonkey";
my $KO_file="$databasedir/ko00000.keg";
my $KO_hydro="$databasedir/hydro.desc";
my %seq;
my @seq;
my @database;
my $workdir = getcwd();

my $PARALLELCMD = "parallel --pipepart --recend '//' --gnu --plain";
my $PARALLELCMD2= "parallel --gnu --plain";
my $starttime = localtime;
my $rnammer_mode = 'bac';
my $barrnap_mode = 'bac';

# K22896 K22897	K22898	K22899 are added to nitfix. Seperated from K02588
my $module_file="$databasedir/mod.cyc.ko.txt";
my $process_file="$databasedir/proc.ko.txt";

(-d $dir) or err("$dir is not a valid a directory. Please supply a directory containing (a)FASTA(s) on the command line");

# command line options 

my(@Options, $quiet, $debug, $kingdom, $force,$outdir, $prefix, $cpus, $gcode, $gffver, $locustag, $increment, $mincontiglen, $eval, $hmms, $centre, $rawproduct, $compliant, $listdb, $citation, $rnammer, $addgenes, $depth, $bothhmms, $checkmqi, $mapping, $nozero, $norrna, $nokegg, $prokka, $trna, $ncrna, $crispr, $gram, $sizeperc, $ekegg, $enokegg, $smalltrgt, $sizepercpart, $part_eval_out, $restore,,$cut_nc, $cut_tc, $threshold, $nt); 

setOptions();

# . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 
# table of tools we need/optional and min versions
my $BIDEC = '(\d+\.\d+)';  # pattern of NN.NN for versions that can be compared
my %tools = (
   'parallel' => {
     GETVER  => "parallel --version | grep '^GNU parallel 2'",
     REGEXP  => qr/GNU parallel (\d+)/,
     MINVER  => "20130422",
     NEEDED  => 1,
   },
   'prodigal' => {
     GETVER  => "prodigal -v 2>&1 | grep -i '^Prodigal V'",
     REGEXP  => qr/($BIDEC)/,
     MINVER  => "2.6",
     MAXVER  => "2.69",  # changed cmdline options in 2.70 git :-/
     NEEDED  => 1,
   },
   'hmmsearch' => {
     GETVER  => "hmmsearch -h | grep '^# HMMER'",
     REGEXP  => qr/HMMER\s+($BIDEC)/,
     MINVER  => "3.1",
     NEEDED  => 1,
   },
   'hmmpress' => {
     GETVER  => "hmmpress -h | grep '^# HMMER'",
     REGEXP  => qr/HMMER\s+($BIDEC)/,
     MINVER  => "3.1", #3.2 might cause problems in building databases because of similar DESC fields
 #    MAXVER  => "3.1",
     NEEDED  => 1,
   },
   'table2asn' => {
    GETVER  => "table2asn -version",
    REGEXP  => qr/table2asn:\s+($BIDEC)/,
    MINVER  => "1.27",
    NEEDED  => 1,
   },
   'rnammer' => {
     GETVER  => "rnammer -V 2>&1 | grep -i 'rnammer [0-9]'",
     REGEXP  => qr/($BIDEC)/,
     MINVER  => "1.2",
     NEEDED  => 0, # only if --rnammer used
   },
   'barrnap' => {
     GETVER  => "barrnap --version 2>&1",
     REGEXP  => qr/($BIDEC)/,  
     MINVER  => "0.4",  
     NEEDED  => 1, # remove if --norrna is used
   },
   'blastn' => {
     GETVER  => "blastn -version",
     REGEXP  => qr/blastn:\s+($BIDEC)/,
     MINVER  => "2.1",
     NEEDED  => 1,
   },
   'blastp' => {
     GETVER  => "blastp -version",
     REGEXP  => qr/blastp:\s+($BIDEC)/,
     MINVER  => "2.1",
     NEEDED  => 0, # only if --prokka is used
   },
   'checkm' => {
     GETVER  => "checkm -h |grep 'CheckM'",
     REGEXP  => qr/($BIDEC)/,
     MINVER  => "1",
     NEEDED  => 0, # only if --checkmqi used
   },
#  'makeblastdb' => {
#    GETVER  => "makeblastdb -version",
#    REGEXP  => qr/makeblastdb:\s+($BIDEC)/,
#    MINVER  => "2.2",
#    NEEDED  => 0,  # only if --proteins used /remove?
#  },
   'prokka' => {
     GETVER  => "prokka --version 2>&1",
     REGEXP  => qr/prokka\s+($BIDEC)/,
     MINVER  => "1.0",
     NEEDED  => 0,
   },
   'samtools' => {
     GETVER  => "samtools 2>&1 | grep 'Version'",
     REGEXP  => qr/Version:\s+($BIDEC)/,
     MINVER  => "1.1",  #This is actually 1.6, but since were at 1.13 now.
     NEEDED  => 0, # only if --mapping used
   },
   'bwa' => {
     GETVER  => "bwa 2>&1 | grep 'Version'",
     REGEXP  => qr/Version:\s+($BIDEC)/,
     MINVER  => "0.7",
     NEEDED  => 0, # only if --mapping used
   },
   'signalp' => {
     # this is so long-winded as -v changed meaning (3.0=version, 4.0=verbose !?)
    GETVER  => "if [ \"`signalp -version 2>&1 | grep -Eo '[0-9]+\.[0-9]+'`\" != \"\" ]; then echo `signalp -version 2>&1 | grep -Eo '[0-9]+\.[0-9]+'`; else signalp -v < /dev/null 2>&1 | egrep ',|# SignalP' | sed 's/^# SignalP-//'; fi",
   #  GETVER  => "signalp -v < /dev/null 2>&1 | egrep ',|# SignalP' | sed 's/^# SignalP-//'",
     REGEXP  => qr/^($BIDEC)/,
     MINVER  => "3.0",
     MAXVER  => "5.0",  # changed cmdline options in 6.0
     NEEDED  => 0,  # only if --prokka or --gram used
   },
   'aragorn' => {
     GETVER  => "aragorn -h 2>&1 | grep -i '^ARAGORN v'",
     REGEXP  => qr/($BIDEC)/,
     MINVER  => "1.2",
     NEEDED  => 0, #only if --prokka or --trna is used 
   }, 
   'minced' => {
     GETVER => "minced --version | sed -n '1p'",
     REGEXP => qr/minced\s+\d+\.(\d+\.\d+)/,
     MINVER => "1.6",
     NEEDED => 0,#only if --crispr or --prokka is used
    },
   'cmscan' => {
     GETVER  => "cmscan -h | grep '^# INFERNAL'",
     REGEXP  => qr/INFERNAL\s+($BIDEC)/,
     MINVER  => "1.1",
     NEEDED  => 0,  # only if --ncrna or --prokka is used 
   },
   'cmpress' => {
     GETVER  => "cmpress -h | grep '^# INFERNAL'",
     REGEXP  => qr/INFERNAL\s+($BIDEC)/,
     MINVER  => "1.1",
     NEEDED  => 0,  #only if --ncrna or --prokka is used 
   },

   # now just the standard unix tools we need
   'less'    => { NEEDED=>1 },
   'grep'    => { NEEDED=>1 },  # yes, we need this before we can test versions :-/
   'egrep'   => { NEEDED=>1 },
   'sed'     => { NEEDED=>1 },
   'find'    => { NEEDED=>1 },
);  
# . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 
# functions to check if tool is installed and correct version

# only check the dependencies that are used
if ($norrna){
   $tools{'barrnap'}{'NEEDED'}=0;
}
if ($rnammer){
   $tools{'rnammer'}{'NEEDED'}=1;
}
if ($prokka){
   $tools{'blastp'}{'NEEDED'}=1;
   $tools{'aragorn'}{'NEEDED'}=1;
   $tools{'minced'}{'NEEDED'}=1;
   $tools{'cmscan'}{'NEEDED'}=1;   
   $tools{'cmpress'}{'NEEDED'}=1;
   $tools{'prokka'}{'NEEDED'}=1;
}
if ($mapping){
   $tools{'samtools'}{'NEEDED'}=1;   
   $tools{'bwa'}{'NEEDED'}=1;
}
if ($ncrna){
   $tools{'cmscan'}{'NEEDED'}=1;   
   $tools{'cmpress'}{'NEEDED'}=1;
}
if ($checkmqi){
   $tools{'checkm'}{'NEEDED'}=1;
}   
if ($gram){
   $tools{'signalp'}{'NEEDED'}=1;
}   
if ($trna){
   $tools{'trna'}{'NEEDED'}=1;
}   
if ($crispr){
   $tools{'minced'}{'NEEDED'}=1;
}

sub check_tool {
   my($toolname) = @_;
   my $t = $tools{$toolname};
   my $fp = find_exe($toolname);
   err("Can't find required '$toolname' in your \$PATH") if !$fp and $t->{NEEDED};
   if ($fp) {
      $t->{HAVE} = $fp;
      msg("Looking for '$toolname' - found $fp");
      if ($t->{GETVER}) {
         my($s) = qx($t->{GETVER});
         if (defined $s) {
            $s =~ $t->{REGEXP};
            $t->{VERSION} = $1 if defined $1;
            msg("Determined $toolname version is $t->{VERSION}");
            if (defined $t->{MINVER} and $t->{VERSION} < $t->{MINVER}) {
               err("Metascan needs $toolname $t->{MINVER} or higher. Please upgrade and try again.");
                                                                       }
            if (defined $t->{MAXVER} and $t->{VERSION} > $t->{MAXVER}) {
               err("Metascan needs a version of $toolname between $t->{MINVER} and $t->{MAXVER}. Please downgrade and try again."); 
                                                                       }
         }
         else {
            err("Could not determine version of $toolname - please install version",
            $t->{MINVER}, "or higher");  # FIXME: or less <= MAXVER if given
         }
      }
   }
}

$ENV{"GREP_OPTIONS"} = '';
for my $toolname (sort keys %tools) {
    if ($tools{$toolname}{'NEEDED'} == 1) {
       check_tool($toolname);
    }
}

open my $bins_fh, '>', "$dir/depths.bins";

#getting the setup data from the files

print "Preparing setup data from files\n";
print "Kegg files\n";
my %keg;
open (KO, "<$KO_file") or die "Can't open $KO_file";
while (my $line = <KO>){
   chomp $line;
   if ($line =~ /^D/) {
      my ($D, $ko, $desc) = split(/\s+/, $line, 3);
      $keg{$ko}=$desc;
   }
}
close (KO);

print "Hydrogen files\n";
open (KO, "<$KO_hydro") or die "Can't open $KO_hydro";
while (my $line = <KO>){
   chomp $line;
   my ($ko, $desc) = split(/\t+/, $line, 2);
   $keg{$ko}=$desc;
}

opendir (DB, $databasedir) or die "Could not open '$databasedir' for reading '$!'\n";
my @hmmdatabases = grep /\.hmm$/i && -f, map "$databasedir/$_", readdir (DB);
closedir (DB);

#not sure whether to include kegg in the hmms only option. For now it does, as there is an option to opt out
if ($bothhmms){
   push @hmmdatabases, $bothhmms;}
if ($hmms){
   undef @hmmdatabases;
   push @hmmdatabases, $hmms;}
if (!$nokegg){
   push @hmmdatabases, $KEGG_hmm;}

if ($nt){
   $databasedir_blastn=$databasedir_blastnt;
}

my @combined_cycles;
my %cycles_krona;
my %combined_cycles;
my %datastructure;

foreach (@hmmdatabases){
   my $datab = fileparse("$_",  ".hmm");
   my %HoA_ACCKO;
   my @getkey;
   my @getvalue;
   my @getacc;
   my %cycles;
   my $rec ={};
   my $cycle_tmp;
   my @cycle_krona;
   open (KEYS, "<$_") or die "Can't open data"; # opens the file read only
   my $name=0;
   my $desc=0;
   my $konmbr=0;
   my $compo=0;
   my $string = '';
   while (my $line = <KEYS>){
      chomp $line;
      $string .=$line;
      if ($line =~ /NAME/) {
         my ($dud, $key) = split(/\s+/, $line, 2);
         push (@getkey,$key);
         $name++;
      }
      if ($line =~ /^ACC/) {
         my ($dud, $acc) = split(/KO\:/, $line, 3);
         push (@getacc,$acc);
         $konmbr++;
      }
      if ($line =~ /DESC/) {
         my ($dud2, $value) = split(/\s+/, $line, 2);
         push (@getvalue, $value);
         $desc++;
      }
      if ($line =~ /COMPO/){$compo++;} 
      if ($line =~ /CYCLE/){ 
         my ($dud3, $cycle) = split(/\s/, $line, 2);
         $cycle_tmp=$cycle;
         push (@cycle_krona, $cycle);
      }
   }
   my @array_cycles  = split("\/\/", $string );
   @array_cycles = grep { $_ ne "\n" } @array_cycles;
   push (@array_cycles, @combined_cycles);

   my $ddiff= $compo-$desc;  
   my $ndiff= $compo-$name;
   my $kdiff= $compo-$konmbr;
   $compo eq $konmbr or print "There are $kdiff entries in the HMM profile-file that lack a ACC field. Please add them in the metadata.\n";
   $compo eq $desc or print "There are $ddiff entries in the HMM profile-file that lack a DESC field. Please add them in the metadata.\n";
   $compo eq $name or die "There are $ndiff entries in the HMM profile-file that lack a NAME field. Please add them in the metadata.\n";
   close (KEYS);

   @cycles{@getkey} = @getvalue;
   @combined_cycles{@getkey} = @getvalue;
   $rec->{genes}    = {%cycles};
   foreach my $i (@getacc){
      $cycles_krona{$i}=$cycle_tmp;
   }
   push @{ $HoA_ACCKO{ $getacc[$_] } }, $getkey[$_] for 0 .. $#getacc; #create hash of arrays for counting overview file
   $datastructure{$datab}= {%HoA_ACCKO};
}


#creating buckets to store overall values for full overviewfile

my $checkm_phyl;
my %big_hash;
my $total_gsum = 0; my $total_osum = 0; my $total_odsum = 0; my $total_gdsum = 0;my $total_keggsum = 0;

open my $tot_ribo_ovw_fh , '>>', "$dir/ribosomal.ovw" or die $!; 

opendir (DH, $dir) or die "Could not open '$dir' for reading '$!'\n";
my @fastas = grep {$_ =~ /\.fasta$/ or $_ =~/\.fna$/ or  $_ =~/\.fa$/} readdir (DH);
msg("Metascan found", scalar(@fastas), "fasta files to analyse");
msg(".*..**...****************METASCAN***************...**..*.");
if (! @fastas) {
   msg("There are no fasta files in this directory. Allowed extentions: .fasta .fna .fa");
}

# Restore hash / data in case the process after days of running

if ($restore) {
   $force=1;
 # identify unanalysed fastas-> substract @fasta_analysed from @fastas
   open my $fastanal_fh , '<', "$dir/analyzedfastas.txt" or die $!;
   chomp(my @analyzedfastas = <$fastanal_fh>);
   close $fastanal_fh;
   my %in_bl = map {$_ => 1} @analyzedfastas;
   my @fastas2  = grep {not $in_bl{$_}} @fastas;
   @fastas = @fastas2;

   #restore data
   %big_hash = %{retrieve("file_hash.txt")};

   $total_osum = scalar(@analyzedfastas);

   open my $gensum_fh , '<', "$dir/gensum.txt" or die $!;
   chomp($total_gsum = <$gensum_fh>);
   close $gensum_fh;

   open my $gendepthsum_fh , '<', "$dir/gendepthsum.txt" or die $!;
   chomp($total_gdsum = <$gendepthsum_fh>);
   close $gendepthsum_fh;

   open my $orgdepthsum_fh , '<', "$dir/orgdepthsum.txt" or die $!;
   chomp($total_odsum = <$orgdepthsum_fh>);
   close $orgdepthsum_fh;

   open my $keggsum_fh , '<', "$dir/keggsum.txt" or die $!;
   chomp($total_keggsum = <$keggsum_fh>);
   close $keggsum_fh;
}

# Store hash in file -> main lines at the end of the loop: line ~2343

open my $fastanal_fh , '>>', "$dir/analyzedfastas.txt" or die $!;

if ($centre){
   if ($centre eq '--compliant'){
      err("Please rename yor command line '--centre --compliant' to '--centre ANYTHING --compliant', otherwise compliant will not run.");
   }
}
#if ($prokka and not $compliant) {
#         err("If you run --prokka without --compliant, you might run into trouble later in the line. Assuming you use --prokka because you are thinking of submitting the data, you might need to run Metascan again");
#      }

my $current_bin;

foreach my $bin (@fastas) {
   my $depth_value=0;
   # . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 
   # welcome message
   msg("This is $EXE $VERSION");
   msg("Written by $AUTHOR");
   msg("Homepage is $URL");
   msg("Local time is $starttime");
   msg("You are", $ENV{USER} || 'not telling me who you are!');
   msg("Operating system is $OPSYS");
   msg("You have enabled DEBUG mode. Temporary files will NOT be deleted.") if $debug;
   msg("Command: @CMDLINE");
   msg("Output E-value setting for HMM:\tWithout Kegg: $enokegg\tWith Kegg: $ekegg");
   msg("E-value cut off for RNA and small proteins: $eval");
   msg("Cut-off value for a protein to be considered small: $smalltrgt ");
   msg("Size range Query-Target length: $sizeperc %");
   msg("Size range Query-Target length partials: $sizepercpart %");
   msg("Working directory: $workdir");

   # . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 
   # check options

   if ($compliant) {
      msg("Enabling options to ensure Genbank/ENA/DDJB submission compliance.");
      $centre ||= 'Metascan';
      $mincontiglen = 200;
   }

   my $in = "$dir/$bin"; 
   $current_bin= $bin;

   # generate locustag
   $locustag ||= generate_locus_tag($in);
   # http://www.ncbi.nlm.nih.gov/genomes/static/Annotation_pipeline_README.txt
   $prefix ||= $locustag; # NCBI wants US format, ech.
   $outdir ||= $prefix;

   #create a reference file to trace back the bin to the folder/locus name
   open my $bin_id_fh , '>>', "$dir/bin.id";
   print {$bin_id_fh} "$locustag\t$in\n";
   close $bin_id_fh;

   if (-d "$dir/$outdir") {
      if ($force) { 
         msg("Re-using existing --outdir $dir/$outdir")
      }
      else {
         err("Folder '$outdir' already exists! Please change --outdir or use --force");
      }
   }
   else {    
      msg("Creating new output folder: $dir/$outdir");
      runcmd("mkdir -p \Q$dir/$outdir\E");
   }
   if (-d "$dir/$outdir/hydrogenases") {
      if ($force) { 
         msg("Re-using existing --outdir $dir/$outdir/hydrogenases")
      }
      else {
         err("Folder '$outdir/hydrogenases' already exists! Please change --outdir or use --force");
      }
   }
   else {    
      msg("Creating new output folder: $dir/$outdir/hydrogenases");
      runcmd("mkdir -p \Q$dir/$outdir/hydrogenases\E");
   }

   msg("Using filename prefix: $prefix.XXX");
   my $logfile = "$dir/$outdir/$prefix.log";
   msg("Writing log to: $logfile");
   open LOG, '>', $logfile or err("Can't open logfile");

   msg("Loading and checking input file: $in");
   my $fin = Bio::SeqIO->new(-file=>$in, -format=>'fasta');
   my $fout = Bio::SeqIO->new(-file=>">$dir/$outdir/$prefix.fna", -format=>'fasta');

   # . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 
   # read in sequences; remove small contigs; replace ambig with N
   my $ncontig = 0;
   my $contigprefix = $locustag || $prefix || $outdir || '';
   $contigprefix .= '_' if $contigprefix;
   my $total_bp = 0;
   while (my $seq = $fin->next_seq) {
      my $id =$seq->id;
      if ($seq->length < $mincontiglen) {
         msg("Skipping short (<$mincontiglen bp) contig: $id");
         next;
      }
      $ncontig++;
      if ($id =~ m/\|/) {
         msg("Changing illegal '|' to '_' in sequence name: $id");
         $id =~ s/\|/_/g;
      }
      if (exists $seq{$id}) {
         err("Uh oh! Sequence file '$in' contains duplicate sequence ID:", $seq->id);
      }
   # http://www.ncbi.nlm.nih.gov/genomes/static/Annotation_pipeline_README.txt
   # leave contigs names as-is unless they are in --compliant mode or want --centre set
      if ($centre) {  
      $id = sprintf "gnl|$centre|Contig%d", $ncontig;  #${contigprefix}%d
      }
      if (length($id) > $MAXCONTIGIDLEN) {
         msg("Contig ID must <= $MAXCONTIGIDLEN chars long: $id");
         err("Please rename your contigs OR try '--centre X --compliant' to generate clean contig names.");
      }
      my $s = $seq->seq;
      $s = uc($s);
      $s =~ s/[*-]//g;      # replace pads/gaps with nothing
      $s =~ s/[^ACTG]/N/g;  # replace wacky IUPAC with N

      $seq->id($id);
      $seq->seq($s);
      $seq->desc(undef);
      $fout->write_seq($seq);

      $seq{$id}{DNA} = $seq;
      push @seq, $id;  # this array it used to preserve original contig order
      $total_bp += $seq->length;
   }

   undef $fin;
   undef $fout; 

   msg("Wrote $ncontig contigs totalling $total_bp bp.");
   if ($ncontig < 1){
      nofasta("FASTA file '$in' contains no suitable sequence entries");  
      next;
   } 
   $total_osum++; #fasta is approved from here on and thus counted


   #Declining --prokka for large bins (unbinned or small contigs)  
   my $toobig;
   if ($total_bp < 15000000){
       $toobig=0;
   }
   elsif ($total_bp >= 15000000){
       $toobig=1;
       msg("Bin too large: Not using --prokka")
   }

   #$centre or err("You must set --centre or the NCBI tbl2asn tool won't work properly, sorry.");
   $eval >= 0 or err("Invalid --evalue, must be >= 0");
   $increment >= 1 or err("Invalid --increment, must be >= 1");
   msg("Setting HMMER_NCPU=1");
   $ENV{HMMER_NCPU} = 1;
   $gcode ||= 11;
   msg("Using genetic code table $gcode.");
   ($gcode < 1 or $gcode > 25) and err("Invalid genetic code, must be 1..25");
             
   # check BioPerl version
   my $minbpver = "1.006002"; # for Bio::SearchIO::hmmer3
   my $bpver = $Bio::Root::Version::VERSION;
   my ($bpv) = $bpver =~ qr/^(\d+\.\d+)/;     
   msg("You have BioPerl $bpver");
   err("Please install BioPerl $minbpver or higher") if $bpv < $minbpver;

   # . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 
   # Determine CPU cores available

   my $num_cores = num_cpu();
   msg("System has $num_cores cores.");
   if (!defined $cpus or $cpus < 0) {
      $cpus = 1;
   }
   elsif ($cpus == 0) {
      $cpus = $num_cores;
   }
   elsif ($cpus > $num_cores) {
      msg("Option --cpu asked for $cpus cores, but system only has $num_cores");
      $cpus = $num_cores;
   }
   msg("Will use maximum of $cpus cores.");

   # . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 
   # set up options based on --mode

#not sure whether to keep this


   if ($kingdom =~ m/bac|prok/i) {
     $kingdom = 'Bacteria';
     $gcode ||= 11;
     $rnammer_mode = 'bac';
     $barrnap_mode = 'bac';
   }
   elsif ($kingdom =~ m/arch/i) {
     $kingdom = 'Archaea';  
     $gcode ||= 11;
     $rnammer_mode = 'arc';
     $barrnap_mode = 'arc';
   }
   else {
     err("Can't parse --mode '$kingdom'. Choose from: Bacteria Archaea");
        }
   msg("Annotating as >>> $kingdom <<<");


   # . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 
   # check if --setupdb has been run for BLAST


   # . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 
   # add our included binaries to the END of the PATH

   add_bundle_to_path();

   if ($rnammer and $tools{'rnammer'}->{HAVE}) {
      msg("Will use RNAmmer instead of Barrnap for rRNA prediction");
      $tools{'barrnap'}->{HAVE} = 0;
                                              }
   # . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 
   # check if optional tools are installed if the option was enabled
   
 # if ($checkmqi and ! $tools{'checkm'}->{HAVE}) {
  #    err("You need to install 'checkm' to use the --checkmqi option.");
   #                                              }
   if ($mapping and ! $tools{'bwa'}->{HAVE}) {
      err("You need to install BWA to use the --mapping option.");
                                             }
   if ($mapping and ! $tools{'metabat2'}->{HAVE}) {
      err("You need to install Metabat2 to use the --mapping option.");
                                                  }
   if ($mapping and ! $tools{'samtools'}->{HAVE}) {
      err("You need to install Samtools to use the --mapping option.");
                                                  }
   if ($prokka || $ncrna and ! $tools{'cmscan'}->{HAVE}) {
      err("You need to install 'cmscan' to use the --prokka or --ncrna option.");
                                                         }
   if ($prokka || $ncrna and ! $tools{'cmpress'}->{HAVE}) {
      err("You need to install 'cmpress' to use the --prokka or --ncrna option.");
                                                          }
   if ($prokka || $trna and ! $tools{'aragorn'}->{HAVE}) {
      err("You need to install 'aragorn' to use the --prokka or --trna option.");
                                                         }
   if ($gram and ! $tools{'signalp'}->{HAVE}) {
      err("You need to install 'signalp' to use the --gram option.");
                                              }
   if ($prokka || $crispr and ! $tools{'minced'}->{HAVE}) {
      err("You need to install 'minced' to use the --prokka or --crispr option.");
                                              }

   # . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 
   # tRNA + tmRNA
   if ($toobig == 0){
      if ($prokka or $trna){
         msg("Using --prokka or --trna : Predicting tRNAs and tmRNAs");
         my $aragorn_opt = '';
         # -l : Assume that each sequence has a linear topology. Search does not wrap
         # -w : batch mode
         my $cmd = "aragorn -l -gc$gcode $aragorn_opt -w \Q$dir/$outdir/$prefix.fna\E"; # -t/-m 
         msg("Running: $cmd");
         my $num_trna=0;
         open TRNA, '-|', $cmd;
         my $sid;
         while (<TRNA>) {
            chomp;
            if (m/^>(\S+)/) {
               $sid = $1;
               next;
            }
            my @x = split m/\s+/;
            next unless @x == 5 and $x[0] =~ m/^\d+$/;
            if ($x[1] =~ m/\?/) {
                  msg("tRNA $x[2] is a pseudo/wacky gene - skipping.");
                  next;
            }   
            msg("@x");
            # in linear mode (-l) we can get -ve coordinates
            $x[2] =~ m/(c)?\[-?(\d+),(\d+)\]/;
            my($revcom, $start, $end) = ($1,$2,$3);
            # bug fix for aragorn when revcom trna ends at start of contig!
            #  if (defined $revcom and $start > $end) { 
            #    msg("Activating kludge for Aragorn bug for tRNA end at contig start");
            #    $start = 1;
            #  }
            if ($start > $end) {
               msg("tRNA $x[2] has start($start) > end ($end) - skipping.");
               next;
            }
            # correct strange coordinates in -l mode
            $start = max( $start, 1 );
            $end = min( $end, $seq{$sid}{DNA}->length );

            if (abs($end-$start) > 500) {
               msg("tRNA/tmRNA $x[2] is too big (>500bp) - skipping.");
               next;
            }
            # end kludge
            $num_trna++;
 
            my $ftype = 'tRNA';
            my $product = $x[1].$x[4];
            my @gene = ();
            if ($x[1] =~ m/^(tmRNA)/) {
               $ftype = $1;
               $product = "transfer-messenger RNA, SsrA";
               @gene = ('gene' => 'ssrA')
            }

            my $tool = "Aragorn:".$tools{aragorn}->{VERSION};
            push @{$seq{$sid}{FEATURE}}, Bio::SeqFeature::Generic->new( 
               -primary    => $ftype, # tRNA or tmRNA
               -seq_id     => $sid,
               -source     => $tool,
               -start      => $start,
               -end        => $end,
               -strand     => (defined $revcom ? -1 : +1),
               -score      => undef,
               -frame      => 0,
               -tag        => {
                 'product' => $product,
                 'inference' => "COORDINATES:profile: $tool",
                 @gene,       }
            );
         }
         msg("Found $num_trna tRNAs");
      }
   }
   elsif ($toobig == 1){ 
      msg("Bin too large; ignoring --prokka or --trna : Not predicting tRNAs and tmRNAs");
   }

   # . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 
   # rRNA
   if (!$norrna){
   msg("Predicting Ribosomal RNAs");
      if ($tools{'barrnap'}->{HAVE}) {
         msg("Running Barrnap with $cpus threads");
         my $num_rrna=0;
         open my $BARRNAP, '-|', "barrnap --kingdom $barrnap_mode --threads $cpus --quiet \Q$dir/$outdir/$prefix.fna\E";
         my $gff = Bio::Tools::GFF->new(-fh => $BARRNAP, -gff_version => 3);
         while (my $feat = $gff->next_feature) {
            $feat->remove_tag('Name'); # only want /product
            push @{$seq{$feat->seq_id}{FEATURE}}, $feat;
            $num_rrna++;
            msg("$num_rrna", $feat->seq_id, $feat->start, $feat->get_tag_values('product'));
         }
      undef $gff;
      msg("Found $num_rrna rRNAs");
   }
   elsif ($tools{'rnammer'}->{HAVE}) {
      msg("Running RNAmmer");
      my $rnammerfn = "$dir/$outdir/rnammer.xml";
      my $num_rrna = 0;
      my $rnammer_opt = $cpus != 1 ? "-multi" : "";
      runcmd("rnammer -S $rnammer_mode $rnammer_opt -xml \Q$rnammerfn\E \Q$dir/$outdir/$prefix.fna\E");
      my $xml = XML::Simple->new(ForceArray => 1);
      my $data = $xml->XMLin($rnammerfn);
      for my $entry (@{$data->{entries}[0]->{entry}}) {
         my $sid = $entry->{sequenceEntry}[0];
         next unless exists $seq{$sid};
         my $desc = $entry->{mol}[0];
         $desc =~ s/s_r/S ribosomal /i; # make it English '23S_rRNA => 23S ribosomal RNA'
         $num_rrna++;
         my $tool = "RNAmmer:".$tools{rnammer}->{VERSION};
         push @{$seq{$sid}{FEATURE}}, Bio::SeqFeature::Generic->new( 
            -primary    => 'rRNA',
            -seq_id     => $sid,
            -source     => $tool, # $data->{predictor}[0]
            -start      => $entry->{start}[0],
            -end        => $entry->{stop}[0],
            -strand     => $entry->{direction}[0],
            -score      => undef, # $entry->{score}[0],
            -frame      => 0,
            -tag        => {
              'product' => $desc,
              'inference' => "COORDINATES:profile: $tool",  # FIXME version ?
                           }
         );
         msg(join "\t", $num_rrna, $desc, $sid, $entry->{start}[0], $entry->{stop}[0], $entry->{direction}[0]);
      }
      delfile($rnammerfn);
      msg("Found $num_rrna rRNAs");
   }
   else {
      msg("You need either Barrnap or RNAmmer installed to predict rRNAs!");
        }
   }
   else {
      msg("Disabling rRNA search: --kingdom=$kingdom or --norrna=$norrna");
   }
   # . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 
   # ncRNA via Rfam + Infernal
   if ($toobig == 0){
      if ($prokka or $ncrna) {
         msg("Using --prokka or --ncrna : Predicting ncRNAs");
         my $cmdb= "$dbloc/prokka/db/cm/Bacteria";
         msg("Preparing HMMER annotation source");
         -r "$cmdb.i1m" or err("Your CM is not indexed, please run: cmpress $cmdb #*Make sure to use the full path");
         if (-r "$cmdb.i1m") {
            msg("Scanning for ncRNAs... please be patient.");
            my $num_ncrna = 0;
            my $tool = "Infernal:".$tools{'cmscan'}->{VERSION};
            my $icpu = $cpus || 1;
            my $cmd = "cmscan --rfam --cpu $icpu -E $eval --tblout /dev/stdout -o /dev/null --noali $cmdb \Q$dir/$outdir/$prefix.fna\E";
            msg("Running: $cmd");
            open INFERNAL, '-|', $cmd;
            while (<INFERNAL>) {
               next if m/^#/;       # ignore comments
               my @x = split ' ';   # magic Perl whitespace splitter
               #    msg("DEBUG: ", join("~~~", @x) );
               next unless @x > 9;  # avoid incorrect lines
               next unless defined $x[1] and $x[1] =~ m/^RF\d/;
               my $sid = $x[2];
               next unless exists $seq{$sid} && $x[0] !~ m/tRNA/ or m/rRNA/ or m/tmRNA/;
               push @{$seq{$sid}{FEATURE}}, Bio::SeqFeature::Generic->new( 
                     -primary    => 'ncRNA',
                     -seq_id     => $sid,
                     -source     => $tool,
                     -start      => min($x[7], $x[8]),
                     -end        => max($x[7], $x[8]),
                     -strand     => ($x[9] eq '-' ? -1 : +1),
                     -score      => undef,  # possibly x[16] but had problems here with '!'
                     -frame      => 0,
                     -tag        => {
                        'product' => $x[0], #trna rrna
                        'inference' => "COORDINATES:profile: $tool",
                                    }
               );
               $num_ncrna++;    
               msg("$num_ncrna ncRNA $x[0] $sid $x[7]..$x[8]");
            } 
            msg("Found $num_ncrna ncRNAs.");
         }
         else {
            msg("Disabling ncRNA search, can't find $cmdb index file.");
         }
      }
   }
   elsif ($toobig == 1){ 
      msg("Bin too large; ignoring --prokka or --ncrna : Not predicting ncRNAs");
   }

   # . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 
   # Tally all the RNA features __ which we want to exclude overlaps with CDS __

   my @allrna;
   for my $sid (@seq) {
      push @allrna, (grep { $_->primary_tag =~ m/[tr]RNA/ } @{ $seq{$sid}{FEATURE} });
   }
   msg("Total of", scalar(@allrna), "tRNA + rRNA features");
   # . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
   # CRISPRs
   if ($toobig == 0){
      if ($prokka or $crispr) {
         msg("Using --prokka or --crispr : Predicting CRISPRS");
         if ($tools{'minced'}->{HAVE}) {
            msg("Searching for CRISPR repeats");
            my $num_crispr=0;
            open my $MINCED, '-|', "minced -gff \Q$dir/$outdir/$prefix.fna\E";
            my $gff = Bio::Tools::GFF->new(-fh => $MINCED, -gff_version => 3);
            while (my $feat = $gff->next_feature) {
               # format it properly for NCBI
               $feat->primary_tag("repeat_region");
               $feat->remove_tag('ID');
               $feat->add_tag_value('rpt_family', 'CRISPR');
               push @{$seq{$feat->seq_id}{FEATURE}}, $feat;
               # there should be no CDS features overlapping with CRISPRs, but prodigal
               # will occasionally create ORFs in these regions. Got to stop that.
               push @allrna, $feat;
               $num_crispr++;
               msg("CRISPR $num_crispr", $feat->seq_id, $feat->start, "with", $feat->score, "spacers");
            }
            msg("Found $num_crispr CRISPRs");
         }
      }
   }
   elsif ($toobig == 1){ 
      msg("Bin too large; ignoring --prokka or --crispr : Not predicting CRISPRS");
   }

   # . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 
   # CDS
   msg("Predicting coding sequences");
   my $prodigal_mode = ($total_bp >= 500000 && $total_bp <= 14000000) ? 'single' : 'meta';
   msg("Contigs total $total_bp bp, so using $prodigal_mode mode");
   my $num_cds=0;
   my $cmd0 = "prodigal -i \Q$dir/$outdir/$prefix.fna\E -m -g $gcode -p $prodigal_mode -f sco -q"; #removed -c because of fragmented genomes
   msg("Running: $cmd0");
   open my $PRODIGAL, '-|', $cmd0;
   my $sid;
   #create a reference file to trace back the bin to the Prodigal metadata;
   open my $prodigal_id_fh , '>>', "$dir/prodigal.txt" or die "Could not open file\n";
   print {$prodigal_id_fh} "$in\n";
   while (<$PRODIGAL>) {
      chomp;
      if (m/seqhdr="([^\s\"]+)"/) {  
      $sid = $1;
      next;
      }
      elsif (m/^>\d+_(\d+)_(\d+)_([+-])$/) {   
         my $tool = "Prodigal:".$tools{prodigal}->{VERSION}; # FIXME: why inner loop?
         my $cds = Bio::SeqFeature::Generic->new(
            -primary    => 'CDS',
            -seq_id     => $sid,
            -source     => $tool,
            -start      => $1,
            -end        => $2,
            -strand     => ($3 eq '+' ? +1 : -1),
            -score      => undef,
            -frame      => 0,
            -tag        => {
              'inference' => "ab initio prediction: $tool", 
                           }
         );
         my $nulen = $cds->length; my $aalen = ($cds->length)/3;
         $cds->add_tag_value('NUlength', $nulen);
         $cds->add_tag_value('AAlength', $aalen);
         
         my $overlap;
         for my $rna (@allrna) {
         # same contig, overlapping (could check same strand too? not sure)
            if ($rna->seq_id eq $sid and $cds->overlaps($rna)) { 
               $overlap = $rna;
               last;
            }	
         }
         # mark partial genes on the ends of the contigs
         if ($cds->start <= 3) { $cds->add_tag_value('part', '1');
         }
         elsif ($cds->end >= $seq{$cds->seq_id}{DNA}->length - 2)  { $cds->add_tag_value('part', '1');
         }
         else {$cds->add_tag_value('part', '0');}         
         # mitochondria are highly packed, so don't exclude as CDS/tRNA often overlap.
         if ($overlap) {
         my $type = $overlap->primary_tag;
            msg("Excluding CDS which overlaps existing RNA ($type) at $sid:$1..$2 on $3 strand");
         }
         else {
            $num_cds++;
            push @{$seq{$sid}{FEATURE}}, $cds;
            ## BUG James Doonan - ensure no odd features extending beyond contig
            if ($cds->end > $seq{$cds->seq_id}{DNA}->length )  { 
               err("CDS end", $cds->end, "is beyond length", $seq{$sid}{DNA}->length, "in contig $sid") 
            }
         }
      }   
   }
   close $prodigal_id_fh; 
   msg("Found $num_cds CDS");
   # . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 
   # Connect features to their parent sequences

   msg("Connecting features back to sequences");
   for my $sid (@seq) {
      for my $f (@{ $seq{$sid}{FEATURE} }) {
         $f->attach_seq( $seq{$sid}{DNA} );
      }
   }
   # . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 
   # Find signal peptide leader sequences

   if ($tools{signalp}->{HAVE}) {
      my $sigpver = substr $tools{signalp}{VERSION}, 0, 1;  # first char, expect 3 or 4 or 5
      if ($kingdom eq 'Bacteria' and $sigpver==3 || $sigpver==4 || $sigpver==5) {
         if ($gram) {
            $gram = $gram =~ m/\+|[posl]/i ? 'gram+' : 'gram-';
            msg("Looking for signal peptides at start of predicted proteins");
            msg("Treating $kingdom as $gram");
            my $spoutfn = "$outdir/signalp.faa";
            open my $spoutfh, '>', $spoutfn;
            my $spout = Bio::SeqIO->new(-fh=>$spoutfh, -format=>'fasta');
            my %cds;
            my $count=0;
            for my $sid (@seq) {
               for my $f (@{ $seq{$sid}{FEATURE} }) {
                  next unless $f->primary_tag eq 'CDS';
                  $cds{++$count} = $f;
                  my $seq = $f->seq->translate(-codontable_id=>$gcode, -complete=>0);
                  $seq->display_id($count);
                  $spout->write_seq($seq);
               }
            }
            if ($count > $SIGNALP_MAXSEQ) {
               msg("Skipping signalp because it can not handle >$SIGNALP_MAXSEQ sequences.");
            }
            else {
               my $opts = $sigpver==3 ? '-m hmm' : '';
               my $cmd = "signalp -t $gram -f short $opts \Q$spoutfn\E 2> /dev/null";
               msg("Running: $cmd");
               my $tool = "SignalP:".$tools{signalp}->{VERSION};
               my $num_sigpep = 0;
               open SIGNALP, '-|', $cmd;
               while (<SIGNALP>) {
                  my @x = split m/\s+/;
                  if ($sigpver == 3) {
                     next unless @x == 7 and $x[6] eq 'Y'; # has sig_pep
                     my $parent = $cds{ $x[0] };
                     my $prob = $x[5];
                     my $cleave = $x[3];
                     my $start = $parent->strand > 0 ? $parent->start : $parent->end;
                     # need to convert to DNA coordinates
                     my $end = $start + $parent->strand * ($cleave*3 - 1);
                     my $sigpep = Bio::SeqFeature::Generic->new( 
                        -seq_id     => $parent->seq_id,
                        -source_tag => $tool,
                        -primary    => 'sig_peptide',
                        -start      => min($start, $end),
                        -end        => max($start, $end),
                        -strand     => $parent->strand,
                        -frame      => 0,    # PHASE: compulsory for peptides, can't be '.'
                        -tag        => {
                        # 'ID' => $ID,
                        # 'Parent' => $x[0],  # don't have proper IDs yet....
                          'product' => "putative signal peptide", 
                          'inference' => "ab initio prediction: $tool", 
                          'note' => "predicted cleavage at residue $x[3] with probability $prob",
                                       }
                     );
                     push @{$seq{$parent->seq_id}{FEATURE}}, $sigpep;
                     $num_sigpep++;
                  }
                  else {
                     # msg("sigp$sigpver: @x");
                     next unless @x==12 and $x[9] eq 'Y'; # has sig_pep
                     my $parent = $cds{ $x[0] };
                     my $cleave = $x[2];
                     my $start = $parent->strand > 0 ? $parent->start : $parent->end;
                     # need to convert to DNA coordinates 
                     my $end = $start + $parent->strand * ($cleave*3 - 1);
                     my $sigpep = Bio::SeqFeature::Generic->new( 
                        -seq_id     => $parent->seq_id,
                        -source_tag => $tool,
                        -primary    => 'sig_peptide',
                        -start      => min($start, $end),
                        -end        => max($start, $end),
                        -strand     => $parent->strand,
                        -frame      => 0,    # PHASE: compulsory for peptides, can't be '.'
                        -tag        => {
                        # 'ID' => $ID,
                        # 'Parent' => $x[0],  # don't have proper IDs yet....
                          'product' => "putative signal peptide", 
                          'inference' => "ab initio prediction: $tool", 
                          'note' => "predicted cleavage at residue $x[2]",
                                       }
                     );
                     push @{$seq{$parent->seq_id}{FEATURE}}, $sigpep;
                     $num_sigpep++;
                  }	
               }
               msg("Found $num_sigpep signal peptides");
            }
            delfile($spoutfn);
         }
         else {
            msg("Option --gram not specified, will NOT check for signal peptides.");
         }
      }
   }
   # . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 
   # Annotate CDS:primary data source is a set of HMM profiles of 126 metabolic proteins

    # set options for stringency

   my $threshold = "-E %e";
   if ($cut_nc) {$threshold = '--cut_nc'};
   if ($cut_tc) {$threshold = '--cut_tc'};

   # these should accept hmm on STDIN and write report to STDOUT
   my $HMMER3CMD = "hmmsearch $threshold   --noali --notextw --acc --cpu 1 - $dir/mytargets.faa";

   foreach my $i (@hmmdatabases) {
      msg("Preparing HMMER annotation source");
      -r "$i.h3i" or err("Your HMM is not indexed, please run: hmmpress $i #*Make sure to use the full path");
      my $src = $i;
      $src =~ s{^.*/}{};
      $src =~ s/.hmm$//;
      msg("Using /inference source as '$src'");
      push @database, {
         DB  => $i,
         SRC => $compliant ? "" : "protein motif:$src:",
         FMT => 'hmmer',
         CMD => $HMMER3CMD,
         ID =>  $src,
         VERSION => 3,
         STRNGT => $threshold, #prefilter
      };
   }
   # we write out all the CDS which haven't been annotated yet and then search them
   my $faa_name = "$dir/mytargets.faa";
   open my $faa, '>', $faa_name;
   my %cds;
   my $count=0;
   for my $sid (@seq) {
      for my $f (@{ $seq{$sid}{FEATURE} }) {
         next unless $f->primary_tag eq 'CDS';
         next if $f->has_tag('product');
	 $cds{++$count} = $f;
	 print $faa ">$count;;$f->{'_gsf_tag_hash'}->{'AAlength'}->[0];;$f->{'_gsf_tag_hash'}->{'part'}->[0]\n",
         $f->seq->translate(-codontable_id=>$gcode, -complete=>0)->seq,"\n";
      }
   }
   close $faa;
   msg("There are still $count unannotated CDS left (started with $num_cds)");
   
   # for each sequence/profile database in order,
   for my $db (@database) {
      # create a unqiue output name so we can save them in --debug mode
      my $outname = $db->{DB};
      $outname =~ s{^.*/}{};
      next if $count <= 0;
      msg("Will use ", $db->{FMT}, "to search against", $db->{DB}, "with $cpus CPUs");
       my $cmd = $db->{CMD};
   #     $cmd =~ s/%i/{}/g;
   #     $cmd =~ s/%o/{}.out/g;
         $cmd =~ s/%e/$eval/g;
         $cmd =~ s,%d,$db->{DB},g;
         $cmd =~ s/$threshold/$db->{STRNGT}/g;

      msg("Annotating with >>> $db->{STRNGT} <<<");

      my $paropts = $cpus > 0 ? " -j $cpus" : "";
      my $hmm_name = "$dir/$outdir/$outname.".$db->{FMT};

      if(-s $db->{DB} > 150000000 ){
         runcmd("${PARALLELCMD}$paropts -a $db->{DB} $cmd > \Q$hmm_name\E ");}#2> /dev/null $KEGG_hmm
      else {
         runcmd("cat $db->{DB} | $cmd > \Q$hmm_name\E 2> /dev/null");} 
  
      my $dbid = $db->{ID};

      my $eval_out =$ekegg;

      if ($nokegg){
      $eval_out =$enokegg;
      }
     
      open my $hmm_fh, '>>', "$dir/$outdir/$prefix.$dbid.hmm.tbl";
      my $bls = Bio::SearchIO->new(-file=>$hmm_name, -format=>$db->{FMT});
      while (my $res = $bls->next_result) {
         while(my $hit = $res->next_hit){ 
            my $name=$hit->name;
            my ($CDS_id,$tlength,$partial) = split( m/;;/, $name, 3);
            my $DESC=$res->{'_querydesc'};
            my $NAME=$res->{'_queryname'};
            my $ACC=$res->{'_queryacc'};
            my $bit=$hit->score;
            my $evalue= $hit->significance;   
            my $qlength=$res->{'_querylength'}; my $Qlow=$qlength/100*(100-$sizeperc); my $Qhigh=$qlength/100*(100+$sizeperc);
                                                my $Qpartlow=$qlength/100*(100-$sizepercpart); my $Qparthigh=$qlength/100*(100+$sizepercpart);
            if ($partial == 0) {
               if($tlength > $Qlow && $tlength < $Qhigh){
                  if ($qlength <= $smalltrgt){ #qlentgh of tlength; not sure yet
                     print {$hmm_fh} "$CDS_id\t$NAME\t$bit\t$evalue\t$dbid\t$prefix\t$ACC\t$DESC\n";
                  }
                  if (($qlength >= $smalltrgt) and ($evalue <= $eval_out)){
                     print {$hmm_fh} "$CDS_id\t$NAME\t$bit\t$evalue\t$dbid\t$prefix\t$ACC\t$DESC\n";
                  }
               }
            }
            else {
               if($qlength > $Qpartlow && $tlength < $Qparthigh){
                  if ($qlength <= $smalltrgt){
                     print {$hmm_fh} "$CDS_id\t$NAME\t$bit\t$evalue\t$dbid\t$prefix\t$ACC\t$DESC\n";
                  }
                  if (($qlength >= $smalltrgt) and ($evalue <= $part_eval_out)){
                     print {$hmm_fh} "$CDS_id\t$NAME\t$bit\t$evalue\t$dbid\t$prefix\t$ACC\t$DESC\n";
                  }
               }
            }
         }
      }
      close $hmm_fh;
      runcmd("cat \Q$dir/$outdir/$prefix.$dbid.hmm.tbl\E >>\Q$dir/$outdir/$prefix.total.tbl\E");#adds up after repeated use, even when using force
      delfile($hmm_name,  "$dir/$outdir/$prefix.$dbid.hmm.tbl"); #<- hmm output
   }
   
   runcmd("sort -nr -t '\t' -k1,1 -k3,3 \Q$dir/$outdir/$prefix.total.tbl\E >\Q$dir/$outdir/$prefix.total.sort.tbl\E");
   delfile ("$dir/$outdir/$prefix.total.tbl");
   runcmd("awk '!x[\$1]++' \Q$dir/$outdir/$prefix.total.sort.tbl\E > \Q$dir/$outdir/$prefix.total.uniq.tbl\E");
  
   my $num_cleaned=0;
   open (UNIQ,  "<", "$dir/$outdir/$prefix.total.uniq.tbl")  or die "Could not open file\n";
   while ( my $line = <UNIQ> ) {
      if ($line !~ /^#/){ 
         chomp $line;
         my ($pid,$qname,$score,$eva,$cycle,$dirname,$acc,$prod) = split(/\t/, $line, 8); 
         my ($gene,$EC) = ('na', 'na');
         if ($prod =~ m/;/){
            ($gene,$prod) = split( m/;\s/, $prod, 2);
         }
         if ($prod =~ m/\[ec.*\]/) {
            ($EC)= $prod =~ /\[(.*?)\]/;
         }
         my @ko = $acc =~ m/K\d+/g;
         foreach my $i(@ko){
            if (defined $cds{$pid}){
               $cds{$pid}->add_tag_value('KO', $i);
            } 
         }
         if ($acc =~ m/hydro/g){
            $cds{$pid}->add_tag_value('KO', 'hydro');
         }
         if ($acc =~ m/UNKNOWN/g){
            $cds{$pid}->add_tag_value('KO', 'unknown');
         }
         my $cleanprod = $rawproduct ? $prod : cleanup_product($prod);
         if ($cleanprod ne $prod) {
            msg("Modify product: $prod => $cleanprod");
            # we remove any special /gene or /EC if the /product is 'hypothetical protein' !
            if ($cleanprod eq $HYPO) {
               if (defined $cds{$pid}){
                  $cds{$pid}->add_tag_value('note', $prod);
   	          $cds{$pid}->remove_tag('gene') if $cds{$pid}->has_tag('gene');
                  $cds{$pid}->add_tag_value('gene', 'hypo');
   	       } 
            }
   	    $num_cleaned++;
         } 
         if (defined $cds{$pid}){
            if ($cds{$pid}->{'_gsf_tag_hash'}->{'part'}->[0] == '0'){
            $cds{$pid}->add_tag_value('product', $cleanprod);
            }
            elsif ($cds{$pid}->{'_gsf_tag_hash'}->{'part'}->[0] == '1'){
            $cds{$pid}->add_tag_value('product', 'partial '.$cleanprod);
            }
            #$cds{$pid}->add_tag_value('product', $cleanprod);
            $cds{$pid}->add_tag_value('EC_number', $EC) if $EC;
            $cds{$pid}->add_tag_value('gene', $gene) if $gene;
            $cds{$pid}->add_tag_value('cycle', $cycle);
            $cds{$pid}->add_tag_value('inference', "$cycle\: ".$prod) if $cycle;
            $cds{$pid}->add_tag_value('eva', $eva);
            $cds{$pid}->add_tag_value('score', $score);
            $cds{$pid}->add_tag_value('key', $qname); 
            $cds{$pid}->add_tag_value('method', 'hmmer'); 
         }
      }
   }
   msg("Cleaned $num_cleaned /product names") if $num_cleaned > 0;
   close UNIQ;
   delfile($faa_name);

   #blast remainig CDS

   msg("BLASTing remaining proteins");

   my $BLASTPCMD = "blastp -query - -db %d -evalue %e -num_threads 1 -num_descriptions 1 -num_alignments 1 -seg no";

   my @blastdatabase = ({
      DB  => "$dbloc/prokka/db/kingdom/Bacteria/sprot",
      SRC => 'similar to AA sequence:UniProtKB:',
      FMT => 'blast',
      CMD => $BLASTPCMD,},
   );
   if ($toobig == 0){
      if ($prokka){
         for my $db (@blastdatabase){
         # we write out all the CDS which haven't been annotated yet and then search them
            my $outname = $db->{DB};
            $outname =~ s{^.*/}{};
            my $faa_name = "$dir/mytargets.faa";
            open my $faa, '>', $faa_name;
            my %cds;
            my $count=0;
            for my $sid (@seq) {
               for my $f (@{ $seq{$sid}{FEATURE} }) {
	          next unless $f->primary_tag eq 'CDS';
                  next if $f->has_tag('product');
	          $cds{++$count} = $f;
	          print $faa ">$count\n",
                  $f->seq->translate(-codontable_id=>$gcode, -complete=>0)->seq,"\n";
               }
            }
            close $faa;
     
            next if $count <= 0;
            msg("There are still $count unannotated CDS left (started with $num_cds)");
            msg("Will use ", $db->{FMT}, "to search against", $db->{DB}, "with $cpus CPUs");
            msg("Annotating CDS, please be patient.");

            my $cmd = $db->{CMD}; # get the right database!!!
          #    $cmd =~ s/%i/{}/g;
          #    $cmd =~ s/%o/{}.out/g;
               $cmd =~ s/%e/$eval/g;
               $cmd =~ s,%d,$db->{DB},g;

          #
          # **** PARALLEL RUN! ****
          #
            my $faa_bytes = -s $faa_name;
            my $bsize = int($faa_bytes / $cpus / 2);
            my $paropts = $cpus > 0 ? " -j $cpus" : "";
   
            my $bls_name = "$dir/$outdir/$outname.".$db->{FMT}; 
            runcmd(
               "cat \Q$faa_name\E | ${PARALLELCMD2}$paropts --block $bsize --recstart '>' --pipe ". 
               "$cmd > \Q$bls_name\E");
   
            my $num_cleaned=0;
            my $bls = Bio::SearchIO->new(-file=>$bls_name, -format=>$db->{FMT}, -version=>$db->{VERSION});
            while (my $res = $bls->next_result) {
               my $hit = $res->next_hit or next;
               my($pid,$prod,$gene,$EC) = ($res->query_name, $hit->description, '', '');
               if ($prod =~ m/~~~/) {
                  ($EC,$gene,$prod) = split m/~~~/, $prod;
                  $EC =~ s/n\d+/-/g; # collapse transitionary EC numbers
               }
               my $cleanprod = $rawproduct ? $prod : cleanup_product($prod);      
               if ($cleanprod ne $prod) {
                  msg("Modify product: $prod => $cleanprod");
               # we remove any special /gene or /EC if the /product is 'hypothetical protein' !
                  if ($cleanprod eq $HYPO) {
   	             $cds{$pid}->add_tag_value('note', $prod);
  	             $cds{$pid}->remove_tag('gene') if $cds{$pid}->has_tag('gene');
	             $cds{$pid}->remove_tag('EC_number') if $cds{$pid}->has_tag('EC_number');
	             $EC = $gene = undef;
	          }
	          $num_cleaned++;
               }
               if ($cds{$pid}->{'_gsf_tag_hash'}->{'part'}->[0] == '0'){
               $cds{$pid}->add_tag_value('product', $cleanprod);
               }
               elsif ($cds{$pid}->{'_gsf_tag_hash'}->{'part'}->[0] == '1'){
               $cds{$pid}->add_tag_value('product', 'partial '.$cleanprod);
               }
               $cds{$pid}->add_tag_value('method', $db->{FMT});
               $cds{$pid}->add_tag_value('EC_number', $EC) if $EC;
               $cds{$pid}->add_tag_value('gene', $gene) if $gene;
               # only add /inference if this DB has a proper SRC to atrribute to
               $cds{$pid}->add_tag_value('inference', "$db->{SRC}\: ".$hit->name) if $db->{SRC};
            }
            msg("Cleaned $num_cleaned /product names") if $num_cleaned > 0;
            delfile($bls_name, $faa_name);
         }
      }
   }
   # . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 
   # Label unannotated proteins as 'hypothetical protein'

   my $num_hypo=0;
   for my $sid (@seq) {
      for my $f ( @{ $seq{$sid}{FEATURE} }) { 
         if ($f->primary_tag eq 'CDS' and not $f->has_tag('product')) {
         next if $f->{'_gsf_tag_hash'}->{'part'}->[0] == '1';
         $f->add_tag_value('product', $HYPO);
         $f->add_tag_value('method', 'non');#needed in line 1377 to prevent hypoths to go into the mix  
         $num_hypo++;
         }
      }
   }
   msg("Labelling remaining $num_hypo proteins as '$HYPO'") if $num_hypo > 0;
   # . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 
   # Add locus_tags and protein_id[CDS only] (and Parent genes if asked)

   msg("Adding /locus_tag identifiers");
   my $num_lt=0;
   for my $sid (@seq) {
      for my $f ( sort { $a->start <=> $b->start } @{ $seq{$sid}{FEATURE} }) {
         next unless $f->primary_tag =~ m/CDS|RNA/;
         $num_lt++;
         my $ID = sprintf("${locustag}_%05d", $num_lt * $increment);
         $f->add_tag_value('ID', $ID);
         $f->add_tag_value('locus_tag', $ID);
         # it seems CDS features _must_ have a valid /protein_id to be output by tbl2asn into .gbk
         if ($centre and $f->primary_tag eq 'CDS') {
            $f->add_tag_value('protein_id', "gnl|$centre|$ID")
         }
         if ($prokka or $addgenes) {
            # make a 'sister' gene feature for the CDS feature
            # (ideally it would encompass the UTRs/promoters as well, but we don't know them)
            my $gene_id = "${ID}_gene";
            my $g = Bio::SeqFeature::Generic->new(
               -primary    => 'gene',
               -seq_id     => $f->seq_id,
               -start      => $f->start,
               -end        => $f->end,
               -strand     => $f->strand,
               -source_tag => $EXE,
               -tag        => { 
                 'locus_tag' => $ID, 
                 'ID'        => $gene_id, # Add suffix to ID for GFF output
                              },
            );
            # Make a Parent tag from the CDS to the gene
            $f->add_tag_value('Parent', $gene_id);
            # copy the /gene from the CDS
            if (my $gENE = TAG($f, 'gene')) {
               $g->add_tag_value('gene', $gENE);
            }
            push @{ $seq{$sid}{FEATURE} }, $g;
         }
      }
   }
   msg("Assigned $num_lt locus_tags to CDS and RNA features.");

   # . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 
   # Write it all out!

   msg("Writing outputs to $dir/$outdir/");
   open my $gff_fh, '>', "$dir/$outdir/$prefix.gff";
   open my $ovw_tsv_fh, '>', "$dir/$outdir/$prefix.tsv";
   open my $kegg_fh, '>', "$dir/$outdir/$prefix.kegg";
   open my $hmm_tot_fh, '>>', "$dir/metagenome.tsv";
   open my $ovw_fh, '>', "$dir/$outdir/$prefix.ovw";
   open my $tbl_fh, '>', "$dir/$outdir/$prefix.tabel";
   my $faa_fh = Bio::SeqIO->new(-file=>">$dir/$outdir/$prefix.hmm.faa", -format=>'fasta');
   my $ffn_fh = Bio::SeqIO->new(-file=>">$dir/$outdir/$prefix.hmm.ffn", -format=>'fasta');
   my $f16_fh = Bio::SeqIO->new(-file=>">$dir/$outdir/$prefix.f16", -format=>'fasta');
   my $fsa_fh = Bio::SeqIO->new(-file=>">$dir/$outdir/$prefix.fsa", -format=>'fasta');
   my $f16fn_fh = Bio::SeqIO->new(-file=>">$dir/$outdir/$prefix.fall", -format=>'fasta');
   my $faaALL_fh = Bio::SeqIO->new(-file=>">$dir/$outdir/$prefix.all.faa", -format=>'fasta');
   my $ffnALL_fh = Bio::SeqIO->new(-file=>">$dir/$outdir/$prefix.all.ffn", -format=>'fasta');
   my $gff_factory = Bio::Tools::GFF->new(-gff_version=>$gffver);
   my $hydrofaa_fh = Bio::SeqIO->new(-file=>">$dir/$outdir/hydrogenases/hydrogenases.faa", -format=>'fasta');
   my $hydroffn_fh = Bio::SeqIO->new(-file=>">$dir/$outdir/hydrogenases/hydrogenases.ffn", -format=>'fasta');

   print {$hmm_tot_fh} "FASTA file\tLocus_tag\tContig\tHMM profile\tBit-score\tE-value\tCycle\tBinfolder\tGene\tEC number\tKO number\tDescription\n";
   print {$gff_fh} "##gff-version $gffver\n";
   for my $id (@seq) {
      print $gff_fh "##sequence-region $id 1 ", $seq{$id}{DNA}->length, "\n";
   }

   my $fsa_desc = "[gcode=$gcode]";
   my %locacc;
   for my $sid (@seq) {
      my $ctg = $seq{$sid}{DNA};
      my $contig_end = $seq{$sid}{DNA}->length;
      $ctg->desc($fsa_desc);
      $fsa_fh->write_seq($ctg);
      $ctg->desc(undef);
      print {$tbl_fh} ">Feature $sid\n";
      for my $f ( sort { $a->start <=> $b->start } @{ $seq{$sid}{FEATURE} }) {
         next if ($f->primary_tag eq 'CDS' and not $f->has_tag('product')); 
         # Add a GFF "Name" tag if we have /gene (better than "ID" /locus_tag)
         if (my $name = TAG($f, 'gene')) {
            $f->add_tag_value('Name', $name);
         }
         # Make sure we have valid frames/phases (GFF column 8)
         $f->frame( $f->primary_tag eq 'CDS' ? 0 : '.' );
    
         print {$gff_fh} $f->gff_string($gff_factory),"\n";
         my ($L,$R) = ($f->strand >= 0) ? ($f->start,$f->end) : ($f->end,$f->start);

         # add > or < to partial translations for gbk file and subsequent submission
         # mark partial genes on the ends of the contigs

         if (($L <= 3 or $L >= $seq{$sid}{DNA}->length - 2) and ($R >= $seq{$sid}{DNA}->length - 2 or $R <= 3)) { print {$tbl_fh} "<$L\t>$R\t",$f->primary_tag,"\n";} #check if this works as it should

         if ($L <= 3 or $L >= $seq{$sid}{DNA}->length - 2) { print {$tbl_fh} "<$L\t$R\t",$f->primary_tag,"\n";
         }
         elsif ($R >= $seq{$sid}{DNA}->length - 2 or $R <= 3)  { print {$tbl_fh} "$L\t>$R\t",$f->primary_tag,"\n";}
         else { print {$tbl_fh} "$L\t$R\t",$f->primary_tag,"\n";
         } 

         for my $tag ($f->get_all_tags) {
            next if $tag =~ m/^[A-Z]/ and $tag !~ m/EC_number/i ; # remove GFF specific tags (start with uppercase letter)
            for my $value ($f->get_tag_values($tag)) {
               print {$tbl_fh} "\t\t\t$tag\t$value\n";
            }
         }
         #overview file
         if ($f->primary_tag eq 'CDS') {
            my $contig = $f ->{'_gsf_seq'}->{'primary_id'};
            my $qname = $f->{'_gsf_tag_hash'}->{'key'}->[0]; #not from blast
            my $score = $f->{'_gsf_tag_hash'}->{'score'}->[0];#not from blast
            my $eva = $f->{'_gsf_tag_hash'}->{'eva'}->[0]; #not from blast
            my $cycle = $f->{'_gsf_tag_hash'}->{'cycle'}->[0]; #not from blast
            my $gene = $f->{'_gsf_tag_hash'}->{'gene'}->[0];
            my $EC = $f->{'_gsf_tag_hash'}->{'EC_number'}->[0];
            my $acc = $f->{'_gsf_tag_hash'}->{'KO'}->[0]; #not from blast
            my $prod= $f->{'_gsf_tag_hash'}->{'product'}->[0];
            my $locus = $f->{'_gsf_tag_hash'}->{'ID'}->[0];
            if ($f->{'_gsf_tag_hash'}->{'method'}->[0] eq 'hmmer') {
               print {$hmm_tot_fh} "$in\t$locus\t$contig\t$qname\t$score\t$eva\t$cycle\t$locustag\t$gene\t$EC\t$acc\t$prod\n";
            }
            if (defined $acc) {$locacc{$locus}=$acc;}
         }
         #Kegg file
         for my $tag ($f->get_all_tags) {
            if ($tag =~ m/KO/){
               my $locus= $f->{'_gsf_tag_hash'}->{'locus_tag'}->[0];
               my $genename= $f->{'_gsf_tag_hash'}->{'Name'}->[0];
               for my $value ($f->get_tag_values($tag)){
                  if ($f->has_tag('Name')){
                     my ($genenam, $moregenes)= split m/,/, $genename;
                     print {$kegg_fh} "$genenam\t$value\n";}
                  else {
                     print {$kegg_fh} "$locus\t$value\n";}
               }
            }
         }
         my $p = $seq{$sid}{DNA}->trunc($f->location);
         $p->display_id( TAG($f, 'locus_tag') );
         $p->desc( TAG($f, 'product') ) if $f->has_tag('product');

         if ($f->primary_tag eq 'CDS' and $f->{'_gsf_tag_hash'}->{'method'}->[0] !~ m/non|blast/) {
            $faa_fh->write_seq($p->translate(-codontable_id=>$gcode, -complete=>0) ); 
            if ($f->{'_gsf_tag_hash'}->{'KO'}->[0] eq 'hydro' ){
               $hydrofaa_fh->write_seq($p->translate(-codontable_id=>$gcode, -complete=>0) );
            }
         }

         if ($f->primary_tag eq 'CDS') {
            $faaALL_fh->write_seq($p->translate(-codontable_id=>$gcode, -complete=>0) ); 
         }

         if ($f->primary_tag =~ m/^(rRNA)$/) {
            $f16_fh->write_seq($p);
         }

         if ($f->primary_tag eq 'CDS' and $f->{'_gsf_tag_hash'}->{'method'}->[0] !~ m/non|blast/) {
            $ffn_fh->write_seq($p); 
               if ($f->{'_gsf_tag_hash'}->{'KO'}->[0] eq 'hydro' ){
                  $hydroffn_fh->write_seq($p); 
               }                     
         }

         if ($f->primary_tag eq 'CDS') {
            $ffnALL_fh->write_seq($p); 
         }

         if ($f->primary_tag eq 'CDS' or $f->primary_tag =~ m/^(rRNA)$/) {
            $f16fn_fh->write_seq($p); 
         }
      }
   } 

   if (@seq) {
      print $gff_fh "##FASTA\n";
      my $seqio = Bio::SeqIO->new(-fh=>$gff_fh, -format=>'fasta');
      for my $sid (@seq) {
         $seqio->write_seq($seq{$sid}{DNA});
      }
   }
   undef $gff_fh;
  
   #getting the value of the bincoverage from an overview file (format: binname\tcoverage\tcheckmlineage)
   my ($bin_base, $x)= split(/.([^.]+)$/, $bin);

   if ($depth) {
      open(DEPTH, "<$depth") or die "Unable to open: $depth\n";
      while(<DEPTH>) { 
         if(grep(/$bin_base/, $_)) { #prefix
            my @tsv_split = split /\t/, $_;
            $depth_value=$tsv_split[1];
            $checkm_phyl=$tsv_split[2];
         }
      }
      close(DEPTH);
   }

   # Binmate pipeline uses Paired End reads only  
   if ($mapping) {
      if ($toobig==0){ #using bin coverage for depth calculation 
         if (length( $depth // '') && $depth_value!=0) {msg("Abundance data already available");}
         else {
            opendir (MAPPING, $mapping) or die "Could not open '$mapping' for reading '$!'\n";
            my @fastq = grep {$_ =~ /fastq/i or $_ =~/fq/i} readdir (MAPPING);
            foreach (@fastq) {
               msg("Mapping Fastq reads to contigs");
               map_sub($bin);
            }
            closedir(MAPPING);
            if (glob("$dir/$outdir/map/$prefix.*.sorted.bam.bai")) {
               my @checkm_tot;
               my $num_con=0;
               my $loopcount=0;
               foreach (@fastq) {
                  $loopcount ++;
                  runcmd("checkm coverage -x $x --all_reads --min_align 0.95 -q -t $cpus $dir/$outdir/map/  $dir/$outdir/map/$prefix.$_.checkm_cov.tsv $dir/$outdir/map/$prefix.$_.sorted.bam");          
                  delfile("$dir/$outdir/map/$prefix.$_.sorted.bam.bai", "$dir/$outdir/map/$prefix.$_.sorted.bam");
                  open(CHECKM, "<$dir/$outdir/map/$prefix.$_.checkm_cov.tsv") or die "Unable to open file\n";
                  while (<CHECKM>)  {
                     next if $. == 1;
                     my @checkm_cov= split (/\t/, "$_");
                     $num_con ++;
                     push @checkm_tot, $checkm_cov[4]; 
                  }
                  close(CHECKM);
               }
               my $contig_count= $num_con/$loopcount;
               my $checkm_tot_sum = sum0 @checkm_tot;
               $depth_value= $checkm_tot_sum /$contig_count;
            }
         }
      }
   }
   if ($mapping) {
      if ($toobig==1){ #using gene coverage for depth calculation
         if (length( $depth // '') && $depth_value!=0)  {msg("Abundance data already available");}
         else {
            opendir (MAPPING, $mapping) or die "Could not open '$mapping' for reading '$!'\n";
            my @fastq = grep {$_ =~ /fastq/i or $_ =~/fq/i} readdir (MAPPING);
            system("cp $dir/$outdir/$prefix.hmm.ffn $dir/");
            foreach (@fastq) {
               msg("Mapping Fastq reads to genes");
               map_sub("$dir/$prefix.hmm.ffn");
            }
            closedir(MAPPING);
            if (glob("$dir/$outdir/map/$prefix.*.sorted.bam.bai")) {
               my @checkm_tot;
               foreach (@fastq) {
                  runcmd("checkm coverage -x ffn --all_reads --min_align 0.95 -q -t $cpus $dir/$outdir/map/  $dir/$outdir/map/$prefix.$_.checkm_cov.tsv $dir/$outdir/map/$prefix.$_.sorted.bam");          
                  delfile("$dir/$outdir/map/$prefix.$_.sorted.bam.bai", "$dir/$outdir/map/$prefix.$_.sorted.bam");
                  my %checkm_cov;
                  open(CHECKM, "<$dir/$outdir/map/$prefix.$_.checkm_cov.tsv") or die "Unable to open file\n";
                  while (<CHECKM>)  {
                     chomp $_;
                     next if $. == 1;
                     my @line= split (/\t/, "$_");
                     $checkm_cov{$line[0]} += $line[4]; 
                  }
                  close(CHECKM);
                  foreach my $cyc (sort keys %datastructure){
                     my $selector = $cyc;
                     foreach my $ko (sort keys %{ $datastructure {$cyc} }) {
                        my @loci;
                        @loci = grep { $locacc{$_} eq $ko } keys %locacc;
                        for my $loc (@loci){
                           $big_hash{$bin}{$cyc}{$ko}->{gendepth} += $checkm_cov{$loc};
                           unless ($selector =~ 'meta.nonkey'){
                              $total_gdsum += $big_hash{$bin}{$cyc}{$ko}{'gendepth'};
                           }
                        }
                     }      
                  }
               }
            }
            delfile("$dir/$prefix.hmm.ffn")
         }
      }
   } 
   elsif (!$mapping){
      if ($toobig==1){
         foreach my $cyc (sort keys %datastructure){
            foreach my $ko (sort keys %{ $datastructure {$cyc} }) {
               $big_hash{$bin}{$cyc}{$ko}->{gendepth}=0;
               $big_hash{$bin}{$cyc}{$ko}->{orgdepth}=0;
            }      
         }
      }
   }
   # CheckM Bin Quality and Identification
   if ($checkmqi){
      if (length( $depth // '' )) {msg("Checkm data already available");}
      else {
         if (-s -d -e "$dir/$outdir/out/") {
            msg("$dir/$outdir/out/ is not empty. Deleting contents from the folder");
            remove_tree("$dir/$outdir/out/", {keep_root => 1});
            msg("Removed contents from $dir/$outdir/out/");
         }
         # Place bins in the reference genome tree
         system('source activate checkm');
         runcmd("checkm tree -x fsa -q -t $cpus --pplacer_threads $cpus  $dir/$outdir/ $dir/$outdir/out/");
         # Assess phylogenetic markers found in each bin
         runcmd("checkm tree_qa -o 2 --tab_table -q $dir/$outdir/out/ > $dir/$outdir/out/$prefix.phylo_markers.tsv");
         # Infer lineage-specific marker sets for each bin
         runcmd("checkm lineage_set -q $dir/$outdir/out/ $dir/$outdir/out/markers.mf");
         # Identify marker genes in bins
         runcmd("checkm analyze -q -t $cpus $dir/$outdir/out/markers.mf -x fsa $dir/$outdir/ $dir/$outdir/out/");
         # Assess bins for contamination and completeness
         runcmd("checkm qa -o 2 -t $cpus -q --tab_table  $dir/$outdir/out/markers.mf $dir/$outdir/out/ > $dir/$outdir/out/$prefix.quality.tsv");
         # Bar plot of bin completeness, contamination, and strain heterogeneity
         runcmd("checkm bin_qa_plot --image_type pdf --dpi 300 -q -x fsa $dir/$outdir/out/ $dir/$outdir/ $dir/$outdir/");
         runcmd("mv $dir/$outdir/bin_qa_plot.pdf $dir/$outdir/$prefix.bin_qa_plot.pdf");
         system('source deactivate');
         open(PHYLO, "<$dir/$outdir/out/$prefix.phylo_markers.tsv") or die "Unable to open file\n";
         while(<PHYLO>) {
             if(grep(/$prefix/, $_)) {
                my @tsv_split = split /\t/, $_;
                $checkm_phyl=$tsv_split[4];
             }
         }
         close(PHYLO);
      } 
   }
            
   #..............................................................
   # 16S blast
   msg("BLAST-ing 16/5/28S");
   #reading fasta files into a hash
   my %blast16;
   my $header;

   open (NAMES, "<$dir/$outdir/$prefix.f16") or die "can not open fasta";
   while (my $line = <NAMES>){
      chomp $line;
      $line =~ s/\r//;
      my $first_char = substr($line,0,1);
      if ($first_char eq ">"){
         $line =~ s/ /_/;   
         $line =~ s/,//;
         $header = substr($line, 1);
      }
      else {
         $blast16{$header} .= $line;}
   }
   close NAMES;

   my $blastnsettings ="\Q6  sseqid length bitscore pident sblastname stitle\E";
   my $BLASTNCMD = "blastn -query - -db %d -evalue %e -num_threads $cpus -max_target_seqs 1 -max_hsps 1 -outfmt $blastnsettings";

   my @blastdb = ({
   DB  => $databasedir_blastn,
   SRC => 'similar to Nt database NCBI:',
   FMT => 'blastn',
   CMD => $BLASTNCMD,
       },
   );

   for my $bl_db (@blastdb) {
      # create a unqiue output name so we can save them in --debug mode
      my $outname = $bl_db->{DB};
      $outname =~ s{^.*/}{};
   
      msg("Will use", $bl_db->{FMT}, "to search against", $bl_db->{DB}, "with $cpus CPUs");

      my $cmd2 = $bl_db->{CMD};
     #   $cmd2 =~ s/%i/{}/g;
     #   $cmd2 =~ s/%o/{}.out/g;
         $cmd2 =~ s/%e/$eval/g;
         $cmd2 =~ s,%d,$bl_db->{DB},g;
 
      print $ovw_fh ("\n\t****Ribosomal identification****\n\n");
      print $tot_ribo_ovw_fh ("\n\t****Ribosomal identification****\t$prefix\n");
      printf $ovw_fh ("%-39s %-31s %-7s %-7s %-7s %-15s %-90s\n", "Query", "Sequence ID", "Length", "Bit", "\% Match", "Subject Name", "Subject Title");
      printf $tot_ribo_ovw_fh ("%-39s %-31s %-7s %-7s %-7s %-15s %-90s\n", "Query", "Sequence ID", "Length", "Bit", "\% Match", "Subject Name", "Subject Title");

      my $bln_name = "$dir/$outdir/$outname.blastn";
      my $temp=undef;
      for $header ( sort keys %blast16){
         print $ovw_fh "$header:\t";
         print $tot_ribo_ovw_fh "$header;\t";
    
         open my $temp, '>', 'temp.txt' or die $!;
         print $temp ($blast16{$header});
         close $temp;
         undef $temp;
         runcmd("cat \Qtemp.txt\E | $cmd2 > \Q$bln_name\E");
             
         open(FILE, "$dir/$outdir/$outname.blastn") or die "Error: no file found.";
         my $output = do {local $/; <FILE> };
         print {$ovw_fh} "$output\n";
         print {$tot_ribo_ovw_fh} "$output\n";  
      }
 
      print $ovw_fh ("\n\t***Checkm RP identification***\n");
      if ($depth or $checkmqi) {print $ovw_fh("\n$checkm_phyl\n");
         printf $tot_ribo_ovw_fh ("\nCheckm: $checkm_phyl\n");
      }
      delfile("$dir/$outdir/$outname.blastn", 'temp.txt')
   }
   # . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 
   # Output a general .txt file with statistics about the annotation
   msg("Generating annotation statistics file");
   open my $txtFh, '>', "$dir/$outdir/$prefix.txt";
   select $txtFh;
   printf "contigs: %d\n", scalar(@seq); 
   printf "bases: %d\n", sum( map { $seq{$_}{DNA}->length } @seq );
   my %count;
   for my $sid (@seq) {
      for my $f (@{ $seq{$sid}{FEATURE} }) {
      $count{ $f->primary_tag }++;
      }
   }
   for my $ft (keys %count) {
      printf "%s: %d\n", $ft, $count{$ft};
   }
   select STDOUT;
   close $txtFh;
   undef $txtFh;
   #......................................................................
   # trying to parse the data into a neat hash

   msg("Generating overview files");
   # now we have to count the occurences of genes(key) in the hmm tab file
   my $count_cycle=0; #excluding kegg
   my $count_kegg=0; #including kegg
   print {$bins_fh} "$bin\t$depth_value\n";
   $total_odsum += $depth_value;
   foreach my $cyc (sort keys %datastructure){ 
      my $selector = $cyc;
      foreach my $ko (sort keys %{ $datastructure {$cyc} }) { #pmo/amoA is used here twice in the count. methane and nitrogen
         my $ocount=0;
         my $count=0;         
         open my $names_fh, '<', "$dir/$outdir/$prefix.total.uniq.tbl";
         while (my $line = <$names_fh>) {
            chomp $line;
            if ($line =~ m/KO\:$ko\b/) {
               $count++;
               $ocount=1;
               $count_kegg++; #one bin
               $total_keggsum++; #all bins
               unless ($selector =~ 'meta.nonkey'){
                  $count_cycle++; #one bin
                  $total_gsum++; #all bins
               }
            }      
         }
         $big_hash{$bin}{$cyc}{$ko}->{gene} = $count;
         $big_hash{$bin}{$cyc}{$ko}->{organism}= $ocount;
         $big_hash{$bin}{$cyc}{$ko}->{depth}= $depth_value;

         my @ko = $ko=~ m/K\d+|hydro\w+|UNKNOWN\d+/g;
         my @kodesc;
         foreach my $i(@ko){
             push @kodesc, $keg{$i};
         }
         my $d = join ",", @kodesc;
         $big_hash{$bin}{$cyc}{$ko}->{description} = $d;

         if ($toobig==0){
            if (length( $depth_value // '') && $depth_value!=0){
               $big_hash{$bin}{$cyc}{$ko}->{orgdepth}= ($depth_value*$ocount);
               $big_hash{$bin}{$cyc}{$ko}->{gendepth}= ($depth_value*$count);
               unless ($selector =~ 'meta.nonkey'){
                  $total_gdsum += $big_hash{$bin}{$cyc}{$ko}{'gendepth'};
               }
            }
            else { 
               $big_hash{$bin}{$cyc}{$ko}->{orgdepth}=0;
               $big_hash{$bin}{$cyc}{$ko}->{gendepth}=0;
            }
         } 
      # gene-centric mapping!!!!!!!!!!!
      }
   }
   #printing out the data

   print $ovw_fh ("\nMetabolic potential\n");
   if ($count_cycle == 0) {print {$ovw_fh} "No key genes found\n";}
   else {print {$ovw_fh} "Total number of key genes: $count_cycle\n";}
   if ($count_kegg == 0) {print {$ovw_fh} "No metabolic genes found\n";}
   else {print {$ovw_fh} "Total number of metabolic genes including kegg: $count_kegg\nNOTE: pmoA/amoA is counted twice in this number\n\n";} 
   if (length( $depth_value // '' )&& $depth_value!=0) {printf $ovw_fh ("%-15s %-5.4f\n\n", "Depth of bin:", $depth_value);}
   else {print {$ovw_fh} "Depth of bin: Not applicable\n";}
   foreach  $bin (keys %big_hash){
      next unless $bin eq $current_bin;
      print {$ovw_fh} "\nbin: $bin\n\n";
      print {$ovw_tsv_fh} "\nbin: $bin\n\n";
      foreach my $cyc (sort keys %{$big_hash{ $bin } } ) {
         my $selector = $cyc;
         unless ($selector =~ 'meta.nonkey'){
            print {$ovw_fh} "\ncycle: $cyc\n\n";
            printf $ovw_fh ("%-15s %-10s %-10s %-60s\n", "Gene", "N#", "Percentage", "Full description");
         }
         print {$ovw_tsv_fh} "\ncycle: $cyc\n\n";
         print {$ovw_tsv_fh} "Gene\tN#\tFull description\n";
         foreach my $ko ( reverse sort {$big_hash{$bin}{$cyc}->{$a}{gene} <=> $big_hash{$bin}{$cyc}->{$b}{gene}} keys %{$big_hash{ $bin }{ $cyc } } ) {
            my $perc;
            if ($count_cycle >= 1) {
               $perc = ($big_hash{$bin}{$cyc}{$ko}{'gene'}/$count_cycle*100); #hydro
            }
            else {$perc=0;}
            if ($nozero) {
               if ($big_hash{$bin}{$cyc}{$ko}->{gene} >= 1) {
                  unless ($selector =~ 'meta.nonkey'){
                     printf $ovw_fh ("%-15s %-10s %-10.2f %-60s \n",$ko,$big_hash{$bin}{$cyc}{$ko}{'gene'},$perc,$big_hash{$bin}{$cyc}{$ko}{'description'});
                  }
                  printf $ovw_tsv_fh ("%s\t%s\t%s\n", $ko,$big_hash{$bin}{$cyc}{$ko}{'gene'},$big_hash{$bin}{$cyc}{$ko}{'description'});
               }
            }
            else {
               unless ($selector =~ 'meta.nonkey'){
                  printf $ovw_fh ("%-15s %-10s %-10.2f %-60s \n",$ko,$big_hash{$bin}{$cyc}{$ko}{'gene'},$perc,$big_hash{$bin}{$cyc}{$ko}{'description'});
               }
               printf $ovw_tsv_fh ("%s\t%s\t%s\n",$ko,$big_hash{$bin}{$cyc}{$ko}{'gene'},$big_hash{$bin}{$cyc}{$ko}{'description'});
            }
         }
      }
   } 
   
   #
   # . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 
   # Use tbl2asn tool to make .gbk and .sqn for us
   # NOTE: the -A and -N  fake accession and version needed for Mauve!
   # BUT: -A xxx has a bug where it uses xxx as the LOCUS (contig) ID for 1st contig
   # SO: *sigh*

    # tbl2asn should not be used for unbinned or small contigs. This takes up to days to do this task and it is not nesecairy as noone is going to submit that data to genbank anyway. Therefor tbl2asn is limited to 15 MB as this is the current maximum size of a bacterial genome plus a bit. Files bigger than this will logically be either unbinnned or small deleted fragments.

   unless ($total_bp >= 15000000){
      my $tbl2asn_opt = @seq > 10_000 ? '-M b' : '-M n';  # Issue 93 - big assemblies

      msg("Generating Genbank and Sequin files");
      runcmd(
        "table2asn -V b -a s -verbose -U -N 1 -c befw -W -locus-tag-prefix $locustag -i \Q$dir/$outdir/$prefix.fsa\E ");  

#      runcmd(
#      "tbl2asn -V b -a r10k -l paired-ends $tbl2asn_opt -N 1 -y 'Annotated using $EXE $VERSION from $URL'".
#      " -Z \Q$dir/$outdir/$prefix.err\E -i \Q$dir/$outdir/$prefix.fsa\E 2> /dev/null"
#      );
      delfile("$dir/$outdir/errorsummary.val");
      delfile( map { "$dir/$outdir/$prefix.$_" } qw(dr fixedproducts ecn val) );

      msg("Repairing broken .GBK output that tbl2asn produces...");
      runcmd("sed 's/COORDINATES: profile/COORDINATES:profile/' < \Q$dir/$outdir/$prefix.gbf\E > \Q$dir/$outdir/$prefix.gbk\E");
      delfile("$dir/$outdir/$prefix.gbf");
   }
   my $genbank  = Bio::SeqIO->newFh(-file => "$dir/$outdir/$prefix.gbk", -format => 'genbank' ); 
   my $embl = Bio::SeqIO->newFh(-file => ">$dir/$outdir/$prefix.embl", -format => 'EMBL' ); 
    #note: you might want to quote -format to keep older perl's from complaining.
    # when you get a BAD LOCUS NAME error; use --compliant, then the name is too long
    print $embl $_ while <$genbank>;
   msg("when you get a BAD LOCUS message error , the contig name is too long. Use --compliant");
   # Some final log output
   msg("Output files:");
   foreach (qx(find \Q$dir/$outdir\E -type f -name "$prefix.*")) {
      chomp;
      msg($_);
   }
   

   msg("Annotation finished successfully.");
   my $endtime = localtime;
   my $walltime = $endtime - $starttime;
   my $pretty = sprintf "%.2f minutes", $walltime->minutes;
   msg("Walltime used: $pretty");
   msg("If you use this result please cite the Metascan paper:");
   msg("This script is based on: Seemann T (2014) Prokka: rapid prokaryotic genome annotation. Bioinformatics. 30(14):2068-9.");
   msg("Type 'metascan . --citation' for more details.");
   msg("************************************");
  
   undef $outdir;
   undef $prefix;
   undef $locustag;
   undef %seq;
   undef @seq;
   undef $sid;
   undef $checkm_phyl;
   undef @database;
   undef $count_cycle;
   close LOG;

   #____________________________________________________________________________________________________-

   sub cleanup_product {
      my $p = shift; 

      # check the whitelist first
      return $p if exists $GOOD_PROD{$p};

      return $HYPO if $p =~ m/DUF\d|UPF\d|conserved|domain of unknown|\b[CN].term|paralog/i; #[CN].term , return pos
   #  return $HYPO if $p =~ m/^[A-Z]+$/;
      return $HYPO if $p !~ m/[a-z]/;

      $p=~ s/\bhomolog( \d)?\b//g;

      $p =~ s/^arCOG\d+\s+//;
      $p =~ s/\((EC|COG).*?\)//;
      $p =~ s/\s+\w+\d{4,}c?//; # remove possible locus tags
      $p =~ s/ and (inactivated|related) \w+//;
      $p =~ s/,\s*family$//;

   #  $p =~ s/\b([A-Z][a-z]{3,5})/\l$1/g;  # lower the case for protein desc initials
   #  $p =~ s/(rossman|willebrand)/\u$1/; # exception to the rule

      $p =~ s/^(potential|possible|probable|predicted|uncharacteri.ed)/putative/i; #potential als in potentiaal , probable,uncharacteri.ed, putative
      if ($p =~ m/(domain|family|binding|fold|like|motif|repeat)\s*$/i and $p !~ m/,/)  { # domain, family, binding, motif, repeat
         $p .= " protein";
      }
      return $p;
   }

   #----------------------------------------------------------------------

   sub TAG {
      my($f, $tag) = @_;
      return unless $f->has_tag($tag);
      return ($f->get_tag_values($tag))[0];
   }

   #----------------------------------------------------------------------

#   sub num_digits {
#      my $n = shift;
#      return max( 2, 1+int(log($n)/log(10)) );
#   }

   #----------------------------------------------------------------------

   sub num_cpu {
      if ( $^O =~ m/linux/i ) {
        my($num) = qx(grep -c ^processor /proc/cpuinfo);
        return $1 if $num =~ m/^(\d+)/;
      }
      elsif ( $^O =~ m/darwin/i ) {
         my($num) = qx(system_profiler SPHardwareDataType | grep Cores);
         return $1 if $num =~ /.*Cores: (\d+)/;
      }
      return 1;
   }

   #----------------------------------------------------------------------

   sub find_exe {
      my($bin) = shift;
      for my $dir (File::Spec->path) {
         my $exe = File::Spec->catfile($dir, $bin);
         return $exe if -x $exe; 
      }
      return;
   }

   #----------------------------------------------------------------------

   sub msg {
      my $t = localtime;
      my $line = "[".$t->hms."] @_\n";
      print STDERR $line unless $quiet;
      if (openhandle(\*LOG)) {
         # write out any buffered log lines
         if (@LOG) {        
            print LOG @LOG;
            @LOG=();
         }
         # write out the current log line
         print LOG $line;
      }
      else {
         # buffer this log line for later
         push @LOG, $line;  
      }
   }

   #----------------------------------------------------------------------

   sub err {
      $quiet=0;
      msg(@_);
      exit(2);
   }

   #----------------------------------------------------------------------

   sub nofasta {
      $quiet=0;
      msg(@_);
      msg("************************************");
      undef $outdir;
      undef $prefix;
      undef $locustag;
      undef %seq;
      undef @seq;
      undef $sid;
      undef $checkm_phyl;
      undef @database;
      undef $count_cycle;
      close LOG;
      return;
   }

   #----------------------------------------------------------------------



   sub runcmd {
      msg("Running:", @_);
      system(@_)==0 or err("Could not run command:", @_);
   }

   #----------------------------------------------------------------------

   sub delfile {
      for my $file (@_) {
         if ($debug) {
            msg("In --debug mode, saving temporary file:", $file);
         }
         else {
            msg("Deleting unwanted file:", $file);
            unlink $file;
         }
      }
   }

   #----------------------------------------------------------------------

   sub version {
      print STDERR "$EXE $VERSION\n";
      exit;
   }

   #----------------------------------------------------------------------

   sub show_citation {
      print STDERR << "EOCITE";
  
      This script is based on Prokka:

      Seemann T, "Prokka: Rapid Prokaryotic Genome Annotation", 
      Bioinformatics, 2014 Jul 15;30(14):2068-9.

      To cite Metascan:

      PMID:$METASCAN_PMID
      doi:$METASCAN_DOI
      http://www.ncbi.nlm.nih.gov/pubmed/$METASCAN_PMID
    
      Thank you.
EOCITE

   exit;
   }

   #----------------------------------------------------------------------

   sub add_bundle_to_path {
      for my $dir ($BINDIR, "$BINDIR/../common", $FindBin::RealBin) {
         if (-d $dir) {
            msg("Appending to PATH: $dir");
            $ENV{PATH} .= ":$dir";
         }
      }
   }

   #----------------------------------------------------------------------

   sub kingdoms {
      return map { m{kingdom/(\w+?)/}; $1 } glob("$dbloc/prokka/db/kingdom/*/*.pin");
   }

   sub genera {
      return map { m{([^/]+)(\.\d+)?\.pin$}; $1 } glob("$dbloc/prokka/db/genus/*.pin");
   }

   sub hmms {
      return map { m{([^/]+)\.hmm\.h3m$}; $1 } glob("$dbloc/prokka/db/hmm/*.h3m");
   }

   sub cms {
      return map { m{([^/]+)\.i1m$}; $1 } glob("$dbloc/prokka/db/cm/*.i1m");
   }

   #----------------------------------------------------------------------

   sub list_db {
      msg( "Looking for Prokka databases in: $dbloc/prokka/db" );
      msg( "* Kingdoms:", kingdoms() );
      msg( "* Genera:", genera() );
      msg( "* HMMs:", hmms() );
     msg( "* CMs:", cms() );
      exit(0);
   }


   #----------------------------------------------------------------------

   sub generate_locus_tag {
      my($fname) = @_;
      msg("Generating locus_tag from '$fname' contents.");
      open my $fh, '<', $fname;
      my $md5 = Digest::MD5->new;
      $md5->addfile($fh);
      my $cs = $md5->hexdigest;
      close $fh;
      my $lt = '';
      for my $i (0 .. 7) {
         my $c = uc(substr($cs,$i,1));
         $c = chr( ord('F') + 1 + $c ) if $c =~ m/^\d$/;
         $lt .= $c;
      }
      msg("Setting --locustag ${lt} from MD5 $cs");
      return $lt;
   }
   #----------------------------------------------------------------------

   # Option setting routines

   sub setOptions {
      use Getopt::Long;

      @Options = (

        "\n Warning, when using the 'centre' option, it NEEDS input on the command line. Metascan will run without one, but it will then not do what you expect it to do.\n And since Metascan can run quite a while it would be an unfortunate waste of time and effort to find out it ignored the next option on the command line.\n This goes for every entry needing a string [X] input, but for centre it may not always be obvious it needs one",
        "\nGeneral:",
        {OPT=>"help",           VAR=>\&usage,                              DESC=>"This help"},
        {OPT=>"version",        VAR=>\&version,                            DESC=>"Print version and exit"},
        {OPT=>"citation",       VAR=>\&show_citation,                      DESC=>"Print citation for referencing Metascan"},
        {OPT=>"quiet!",         VAR=>\$quiet,         DEFAULT=>0,          DESC=>"No screen output"},
        {OPT=>"debug!",         VAR=>\$debug,         DEFAULT=>0,          DESC=>"Debug mode: keep all temporary files"},
        {OPT=>"restore!",       VAR=>\$restore,       DEFAULT=>0,          DESC=>"Restore data and restart from breaking point"},
       
        "\nSetup:",
        {OPT=>"listdb",         VAR=>\&list_db,                            DESC=>"List all configured Prokka databases"},
#        {OPT=>"setupdb",        VAR=>\&setup_db,                           DESC=>"Index all installed databases. Manually please until thoroughly checked"},
#        {OPT=>"cleandb",        VAR=>\&clean_db,                           DESC=>"Remove all database indices"},
#        {OPT=>"depends",        VAR=>\&list_depends,                       DESC=>"List all software dependencies"},

        "\nOutputs:",
        #{OPT=>"outdir=s",      VAR=>\$outdir,        DEFAULT=>'',         DESC=>"Output folder [auto]"},
        {OPT=>"force!",         VAR=>\$force,         DEFAULT=>0,          DESC=>"Force overwriting existing output folder"},
        {OPT=>"prefix=s",       VAR=>\$prefix,        DEFAULT=>'',         DESC=>"Filename output prefix [auto]"},
        {OPT=>"locustag=s",     VAR=>\$locustag,      DEFAULT=>'',         DESC=>"Locus tag prefix [auto]"},
        {OPT=>"increment=i",    VAR=>\$increment,     DEFAULT=>1,          DESC=>"Locus tag counter increment"},
        {OPT=>"gffver=i",       VAR=>\$gffver,        DEFAULT=>3,          DESC=>"GFF version"},
        {OPT=>"addgenes!",      VAR=>\$addgenes,      DEFAULT=>'',         DESC=>"Add 'gene' feature to each 'CDS' "},
        {OPT=>"compliant!",     VAR=>\$compliant,     DEFAULT=>0,          DESC=>"Force Genbank/ENA/DDJB compliance: --addgenes --mincontiglen 200 --centre XXX"},
        {OPT=>"centre=s",       VAR=>\$centre,        DEFAULT=>'',         DESC=>"Sequencing centre ID. This option NEEDS command line input!"},
        {OPT=>"nozero!",        VAR=>\$nozero,        DEFAULT=>0,          DESC=>"Do not list negative hits in output"},

        "\nAnnotations:",
        {OPT=>"kingdom=s",      VAR=>\$kingdom,       DEFAULT=>'Bacteria', DESC=>"rRNA mode: Bacteria Archaea" },
        {OPT=>"gcode=i",        VAR=>\$gcode,         DEFAULT=>0,          DESC=>"Genetic code / Translation table (set if --kingdom is set)"},
        {OPT=>"hmms=s",         VAR=>\$hmms,          DEFAULT=>'',         DESC=>"User supplied HMM to annotate from"},
        {OPT=>"bothhmms=s",     VAR=>\$bothhmms,      DEFAULT=>'',         DESC=>"Both User supplied and metabolic HMM to annotate from"},
        {OPT=>"rawproduct!",    VAR=>\$rawproduct,    DEFAULT=>0,          DESC=>"Do not clean up /product annotation"},
        {OPT=>"norrna!",        VAR=>\$norrna,        DEFAULT=>0,          DESC=>"Don't run rRNA search"},
        {OPT=>"nokegg!",        VAR=>\$nokegg,        DEFAULT=>0,          DESC=>"Don't run Kegg database for KEGG"},
        {OPT=>"gram=s",         VAR=>\$gram,          DEFAULT=>'',         DESC=>"Gram: -/neg +/pos"},
 
        "\nComputation:",
        {OPT=>"cut_nc!",       VAR=>\$cut_nc,         DEFAULT=>'',         DESC=>"Use cut-nc instead of evalue for HMM threshold. WILL NOT RESULT IN OVERVIEW FILES (yet)"},
        {OPT=>"cut_tc!",       VAR=>\$cut_tc,         DEFAULT=>'',         DESC=>"Use cut-tc instead of evalue for HMM threshold. WILL NOT RESULT IN OVERVIEW FILES (yet)"},
        {OPT=>"cpus=i",         VAR=>\$cpus,          DEFAULT=>8,          DESC=>"Number of CPUs to use [0=all]"},
        {OPT=>"mincontiglen=i", VAR=>\$mincontiglen,  DEFAULT=>1,          DESC=>"Minimum contig size [NCBI needs 200]"},
        {OPT=>"rnammer!",       VAR=>\$rnammer,       DEFAULT=>0,          DESC=>"Prefer RNAmmer over Barrnap for rRNA prediction"},
        {OPT=>"evalue=f",       VAR=>\$eval,          DEFAULT=>1E-06,      DESC=>"Similarity e-value cut-off for RNA and small proteins"},
        {OPT=>"e-kegg=f",       VAR=>\$ekegg,         DEFAULT=>1E-50,      DESC=>"E-value cut-off for big proteins with the use of the Kegg database"},
        {OPT=>"e-nokegg=f",     VAR=>\$enokegg,       DEFAULT=>1E-100,     DESC=>"E-value cut-off for big proteins for key genes only"},
        {OPT=>"e-partial=f",    VAR=>\$part_eval_out, DEFAULT=>1E-80,      DESC=>"E-value cut-off for partials genes"},
        {OPT=>"size=f",         VAR=>\$sizeperc,      DEFAULT=>20,         DESC=>"Size range Query-Target length; 20 => 80-120%"},
        {OPT=>"size-part=f",    VAR=>\$sizepercpart,  DEFAULT=>30,         DESC=>"Size range Query-Target length for partial genes; 30 => 70-130%"},
        {OPT=>"smalltrgt=f",    VAR=>\$smalltrgt,     DEFAULT=>200,        DESC=>"M--centre XXXax value for a target sequence to be considered a small protein. Lower is more stringent"},
  
        "\nAdditional Options:",
        {OPT=>"depth=s",        VAR=>\$depth,         DEFAULT=>'',         DESC=>"Include Depth of Genes. Use the Binmate TSV overview file"},
        {OPT=>"checkm!",        VAR=>\$checkmqi,      DEFAULT=>'',         DESC=>"CheckM for Bin Quality Control and Identification"},
        {OPT=>"mapping=s",      VAR=>\$mapping,       DEFAULT=>'',         DESC=>"Map reads to genes. Requires dir name containing FASTQ file(s)"},
        {OPT=>"nt!",            VAR=>\$nt,            DEFAULT=>0,          DESC=>"Run nt databse for 16S BLAST"}, 

        "\nProkka annotation options",
        {OPT=>"prokka!",        VAR=>\$prokka,        DEFAULT=>'',         DESC=>"Use all prokka options; Combine with --compliant when submition is an option"},
        {OPT=>"trna!",          VAR=>\$trna,          DEFAULT=>'',         DESC=>"Search for tRNA and tmRNA"},
        {OPT=>"ncrna!",         VAR=>\$ncrna,         DEFAULT=>'',         DESC=>"Search for ncRNA"},
        {OPT=>"crispr!",        VAR=>\$crispr,        DEFAULT=>'',         DESC=>"Search for CRISPRS"},
      );
  
      &GetOptions(map {$_->{OPT}, $_->{VAR}} grep { ref } @Options) || usage();

      # Now setup default values.
      foreach (@Options) {
         if (ref $_ && defined($_->{DEFAULT}) && !defined(${$_->{VAR}})) {
            ${$_->{VAR}} = $_->{DEFAULT};
         }
      }
   }
   #----------------------------------------------------------------------
   # usage

   sub usage {
      print STDERR 
      "Name:\n  ", ucfirst($EXE), " $VERSION by $AUTHOR\n",
      "Synopsis:\n  rapid prokarotic gene identification\n",
      "Usage:\n  $EXE <dir> [options] \n";
      foreach (@Options) {
         if (ref) {
            my $def = defined($_->{DEFAULT}) ? " (default '$_->{DEFAULT}')" : "";
            $def = ($def ? ' (default OFF)' : '(default ON)') if $_->{OPT} =~ m/!$/;
            my $opt = $_->{OPT};
            $opt =~ s/!$//; 
            $opt =~ s/=s$/ [X]/; 
            $opt =~ s/=i$/ [N]/;
            $opt =~ s/=f$/ [n.n]/;
            printf STDERR "  --%-15s %s%s\n", $opt, $_->{DESC}, $def;
         }
         else {
            print STDERR "$_\n";
         }      
      }
      exit(1);
   }
   #-------------------------------------------------------------------------
   # Mapping subroutine $bin #-k 49 -B 10 -O 7 -E 6 -w 2 -T 50
   sub map_sub {
      if (! -d "$dir/$outdir/map" or ! -e "$dir/$outdir/map/@_"){
         runcmd("mkdir -p $dir/$outdir/map/ && cp $dir/@_ $dir/$outdir/map/");
      }
      if (! -e "$dir/$outdir/map/@_.amb"){
         runcmd("bwa index $dir/$outdir/map/@_");
      }
      if (-s "$dir/$outdir/map/@_"){ 
         runcmd("bwa mem -M -t $cpus  $dir/$outdir/map/@_ $mapping/$_ > $dir/$outdir/map/$prefix.$_.sam");
         runcmd("samtools view -@ $cpus -h -bt $dir/$outdir/map/@_ -o $dir/$outdir/map/$prefix.$_.unsorted.bam $dir/$outdir/map/$prefix.$_.sam");
         runcmd("samtools sort -@ $cpus $dir/$outdir/map/$prefix.$_.unsorted.bam -o $dir/$outdir/map/$prefix.$_.sorted.bam");
         runcmd("samtools index -b -@ $cpus $dir/$outdir/map/$prefix.$_.sorted.bam");
         delfile ("$dir/$outdir/map/$prefix.$_.unsorted.bam", "$dir/$outdir/map/$prefix.$_.sam")
      }
   }
   #store hash and data in a file:
   store(\%big_hash, "$dir/file_hash.txt");

   open my $gensum_fh , '>', "$dir/gensum.txt" or die $!;
   print {$gensum_fh} "$total_gsum\n";
   close $gensum_fh;

   open my $gendepthsum_fh , '>', "$dir/gendepthsum.txt" or die $!;
   print {$gendepthsum_fh} "$total_gdsum\n";
   close $gendepthsum_fh;
  
   open my $orgdepthsum_fh , '>', "$dir/orgdepthsum.txt" or die $!;
   print {$orgdepthsum_fh} "$total_odsum\n";

   close $orgdepthsum_fh;

   open my $keggsum_fh , '>', "$dir/keggsum.txt" or die $!;
   print {$keggsum_fh} "$total_keggsum\n";
   close $keggsum_fh;

   print {$fastanal_fh} "$current_bin\n"; #organisms/bins

} 

closedir (DH); #END MAINLOOP
close $fastanal_fh;
#----------------------------------------------------------------------

#sub setup_db {

#   add_bundle_to_path();
#   clean_db();

#   check_tool('makeblastdb');
#   for my $sprot (<$DBDIR/kingdom/*/sprot>) {
#      msg("Making kingdom BLASTP database: $sprot");
#      runcmd("makeblastdb -hash_index -dbtype prot -in \Q$sprot\E -logfile /dev/null");
#   }
#   for my $genus (<$DBDIR/genus/*>) {
#      msg("Making genus BLASTP database: $genus");
#      runcmd("makeblastdb -hash_index -dbtype prot -in \Q$genus\E -logfile /dev/null");
#   }

#   check_tool('hmmpress');
#   for my $hmm (<$DBDIR/hmm/*.hmm>) {
#      msg("Pressing HMM database: $hmm");
#      runcmd("hmmpress \Q$hmm\E");
#   }

#   check_tool('cmpress');
#      for my $cm (<$DBDIR/cm/{Viruses,Bacteria}>) {
#      msg("Pressing CM database: $cm");    
#      runcmd("cmpress \Q$cm\E");
#   }
 
#   list_db();
#}

#----------------------------------------------------------------------

#sub list_depends {
#   for my $t (sort keys %tools) {
#      my $ver = $tools{$t}{MINVER} || '0';
#      my $opt = $tools{$t}{NEEDED} ? 'compulsory' : 'optional';
#      my $inc = -x "$BINDIR/$t" || -x "$BINDIR/../common/$t" || -x "$BINDIR/../../bin/$t" 
#           ? 'bundled' 
#           : 'not bundled';
#      print "$t >= $ver ($opt, $inc)\n";
#   }
#   exit(0);
#}

#----------------------------------------------------------------------

#sub clean_db {
#   msg("Cleaning databases in $DBDIR");
#   delfile( <$DBDIR/kingdom/*/sprot.p??> );
#   delfile( <$DBDIR/genus/*.p??> );
#   delfile( <$DBDIR/hmm/*.h3?> );
#   delfile( <$DBDIR/cm/*.i1?> );
#   msg("Cleaning complete.");
#}

#______________________________________________________________________________________________________________

# Creating overview file of the total samples.

#tally up all genes, organisms and depth

my %ovw_results;
my %ko_results;

for my $bin (keys %big_hash) {
   for my $cyc (keys %{$big_hash{ $bin }}) {
      for my $ko(keys %{$big_hash{ $bin }{ $cyc }} ) {
         $ovw_results{$cyc}{$ko}->{description}= $big_hash{ $bin }{ $cyc }{ $ko }{'description'};
         $ko_results{$ko}->{description}= $big_hash{ $bin }{ $cyc }{ $ko }{'description'};
#---------------------------------------------------------------#
         $ovw_results{$cyc}{$ko}->{totg} += $big_hash{$bin}{$cyc}{$ko}{'gene'};
         $ko_results{$ko}->{totg} += $big_hash{$bin}{$cyc}{$ko}{'gene'};
         if ($total_gsum != 0){
            $ovw_results{$cyc}{$ko}->{percg} = ($ovw_results{$cyc}{$ko}{'totg'}/$total_gsum*100);
            $ko_results{$ko}->{percg} = ($ovw_results{$cyc}{$ko}{'totg'}/$total_gsum*100);
         }
         else {$ovw_results{$cyc}{$ko}->{percg}='0';
         $ko_results{$ko}->{percg}='0';}
#---------------------------------------------------------------#
         $ovw_results{$cyc}{$ko}->{toto} += $big_hash{$bin}{$cyc}{$ko}{'organism'};
         $ko_results{$ko}->{toto} += $big_hash{$bin}{$cyc}{$ko}{'organism'};
         if ($total_osum != 0){
            $ovw_results{$cyc}{$ko}->{perco} =($ovw_results{$cyc}{$ko}{'toto'}/$total_osum*100);
            $ko_results{$ko}->{perco} =($ovw_results{$cyc}{$ko}{'toto'}/$total_osum*100);
         }
         else {$ovw_results{$cyc}{$ko}->{perco}='0';
         $ko_results{$ko}->{perco}='0';}
#---------------------------------------------------------------#
         $ovw_results{$cyc}{$ko}->{totod} += $big_hash{$bin}{$cyc}{$ko}{'orgdepth'};
         $ko_results{$ko}->{totod} += $big_hash{$bin}{$cyc}{$ko}{'orgdepth'};
         if ($total_odsum != 0){
            $ovw_results{$cyc}{$ko}->{percod}=($ovw_results{$cyc}{$ko}{'totod'}/$total_odsum*100); #depth per organism
            $ko_results{$ko}->{percod}=($ovw_results{$cyc}{$ko}{'totod'}/$total_odsum*100); #depth per organism
         }
         else {$ovw_results{$cyc}{$ko}->{percod}='0';
         $ko_results{$ko}->{percod}='0';}
#---------------------------------------------------------------#
         $ovw_results{$cyc}{$ko}->{totgd} += $big_hash{$bin}{$cyc}{$ko}{'gendepth'};
         $ko_results{$ko}->{totgd} += $big_hash{$bin}{$cyc}{$ko}{'gendepth'};
         if ($total_gdsum != 0){
            $ovw_results{$cyc}{$ko}->{percgd}=($ovw_results{$cyc}{$ko}{'totgd'}/$total_gdsum*100); #depth per gene
            $ko_results{$ko}->{percgd}=($ovw_results{$cyc}{$ko}{'totgd'}/$total_gdsum*100); #depth per gene
         }
         else {$ovw_results{$cyc}{$ko}->{percgd}='0';
         $ko_results{$ko}->{percgd}='0';}
      }
   }
}

my %modulesKO;
my @c;
open (MKO, "<$module_file") or die "Can't open $module_file";
while (my $line = <MKO>){
   chomp $line;
   my ($c, $module, $ko_mod) = split(/\t/, $line, 3);
   my @ko_mod = split "\t", $ko_mod;
   for my $ko (@ko_mod){
       $modulesKO{$c}{$module}->{$ko}=$ko_results{$ko};
   }
}
close (MKO);

my %processKO;
my @p;
open (PKO, "<$process_file") or die "Can't open $process_file";
while (my $line = <PKO>){
   chomp $line;
   my ($p, $proc, $ko_proc) = split(/\t/, $line, 3);  
   my @ko_proc = split "\t", $ko_proc;
   for my $ko (@ko_proc){
       $processKO{$p}{$proc}->{$ko}=$ko_results{$ko};
   }
}
close (PKO);

open my $krona_mod_pd_fh , '>', "$dir/krona.mod.gd.tsv" or die $!;
open my $krona_mod_g_fh , '>', "$dir/krona.mod.g.tsv" or die $!;
open my $krona_mod_o_fh , '>', "$dir/krona.mod.o.tsv" or die $!;
open my $krona_mod_d_fh , '>', "$dir/krona.mod.od.tsv" or die $!;
open my $krona_proc_pd_fh , '>', "$dir/krona.proc.gd.tsv" or die $!;
open my $krona_proc_g_fh , '>', "$dir/krona.proc.g.tsv" or die $!;
open my $krona_proc_o_fh , '>', "$dir/krona.proc.o.tsv" or die $!;
open my $krona_proc_d_fh , '>', "$dir/krona.proc.od.tsv" or die $!;
open my $krona_pd_fh , '>', "$dir/krona.gd.tsv" or die $!;
open my $krona_g_fh , '>', "$dir/krona.g.tsv" or die $!;
open my $krona_o_fh , '>', "$dir/krona.o.tsv" or die $!;
open my $krona_d_fh , '>', "$dir/krona.od.tsv" or die $!;



#printing out the data
open my $totovw_fh , '>', "$dir/total.ovw" or die $!;
open my $totovw_tsv_fh , '>', "$dir/total.tsv" or die $!;

print {$totovw_fh} "Metabolic Overview\n\n";
print {$totovw_tsv_fh} "Metabolic Overview\n\n";

if ($total_gsum == 0 and $total_keggsum == 0) {
   print {$totovw_fh} "No metabolic genes found\n";
   print {$totovw_tsv_fh} "No metabolic genes found\n";}

else {
   print {$totovw_fh} "Total key genes: $total_gsum\n";
   if (!$nokegg){print {$totovw_fh} "Total metabolic genes including kegg : $total_keggsum\nNOTE: pmoA/amoA is counted twice in this number\n\n";}
   print {$totovw_tsv_fh} "Total key genes: $total_gsum\n";}
   if (!$nokegg){print {$totovw_tsv_fh} "Total metabolic genes including kegg: $total_keggsum\nNOTE: pmoA/amoA is counted twice in this number\n\n";}
if ($total_osum == 0) {
   print {$totovw_fh} "No organisms found\n";
   print {$totovw_tsv_fh} "No organisms found\n";}

else {
   print {$totovw_fh} "Total number of organisms: $total_osum\n"; 
   print {$totovw_tsv_fh} "Total number of organisms in cycles: $total_osum\n";} 

if (length( $total_odsum // '' )&& $total_odsum!=0 || length( $total_gdsum // '' )&& $total_gdsum!=0) {
   printf $totovw_fh ("%-15s %-5.2f\n", "Depth of bin at organism level:", $total_odsum);
   printf $totovw_tsv_fh ("%-15s %-5.2f\n", "Depth of bin at organism level:", $total_odsum);
   printf $totovw_fh ("%-15s %-5.2f\n", "Depth of bin at gene level:", $total_gdsum);
   printf $totovw_tsv_fh ("%-15s %-5.2f\n", "Depth of bin at gene level:", $total_gdsum);}

else {
   print {$totovw_fh} "Depth of bin: Not applicable\n";
   print {$totovw_tsv_fh} "Depth of bin: Not applicable\n";}

for my $cyc (sort keys %ovw_results){
   my $selector = $cyc;
   unless ($selector =~ 'meta.nonkey'){
      print {$totovw_fh} "\n$cyc\n\n";
      printf $totovw_fh ("%-10s %-10s %-10s %-10s %-10s %-10s %-10s %-10s %-10s %-60s\n", "Gene", "N# genes", '%genes', "N# org",'%Org' ,"O-Depth", '%O-Depth',"G-Depth",'%G-Depth', "Full description");
   }
   print {$totovw_tsv_fh} "\n$cyc\n\n";
   printf $totovw_tsv_fh  ("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", "Gene", "N# genes", '%genes', "N# org",'%Org' ,"O-Depth", '%O-Depth',"G-Depth",'%G-Depth', "Full description");
   for my $ko (reverse sort {$ovw_results{$cyc}->{$a}{totg} <=> $ovw_results{$cyc}->{$b}{totg}} keys %{$ovw_results { $cyc }} ) {
      if ($total_gsum != 0){
         if ($nozero) { 
            if ($ovw_results{$cyc}{$ko}{'totg'} >= 1) {
               unless ($selector =~ 'meta.nonkey'){
                  printf $totovw_fh ("%-10s %-10s %-10.2f %-10s %-10.2f %-10.2f %-10.2f %-10.2f %-10.2f %-60s \n",$ko, $ovw_results{$cyc}{$ko}{'totg'}, $ovw_results{$cyc}{$ko}{'percg'}, $ovw_results{$cyc}{$ko}{'toto'}, $ovw_results{$cyc}{$ko}{'perco'}, $ovw_results{$cyc}{$ko}{'totod'}, $ovw_results{$cyc}{$ko}{'percod'}, $ovw_results{$cyc}{$ko}{'totgd'}, $ovw_results{$cyc}{$ko}{'percgd'},  $ovw_results{$cyc}{$ko}{'description'});
                  printf $totovw_tsv_fh ("%s\t%s\t%.2f\t%s\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%s\n",$ko, $ovw_results{$cyc}{$ko}{'totg'}, $ovw_results{$cyc}{$ko}{'percg'}, $ovw_results{$cyc}{$ko}{'toto'}, $ovw_results{$cyc}{$ko}{'perco'}, $ovw_results{$cyc}{$ko}{'totod'}, $ovw_results{$cyc}{$ko}{'percod'}, $ovw_results{$cyc}{$ko}{'totgd'}, $ovw_results{$cyc}{$ko}{'percgd'},  $ovw_results{$cyc}{$ko}{'description'});
                  print {$krona_pd_fh} "$ovw_results{$cyc}{$ko}{'totgd'}\t$cyc\t$ko\t$ovw_results{$cyc}{$ko}{'description'}\n";
                  print {$krona_g_fh} "$ovw_results{$cyc}{$ko}{'totg'}\t$cyc\t$ko\t$ovw_results{$cyc}{$ko}{'description'}\n";
                  print {$krona_o_fh} "$ovw_results{$cyc}{$ko}{'toto'}\t$cyc\t$ko\t$ovw_results{$cyc}{$ko}{'description'}\n";
                  print {$krona_d_fh} "$ovw_results{$cyc}{$ko}{'totod'}\t$cyc\t$ko\t$ovw_results{$cyc}{$ko}{'description'}\n";
               }
            }
         }
         else {
            unless ($selector =~ 'meta.nonkey'){
               printf $totovw_fh ("%-10s %-10s %-10.2f %-10s %-10.2f %-10.2f %-10.2f %-10.2f %-10.2f  %-60s \n",,$ko, $ovw_results{$cyc}{$ko}{'totg'}, $ovw_results{$cyc}{$ko}{'percg'}, $ovw_results{$cyc}{$ko}{'toto'}, $ovw_results{$cyc}{$ko}{'perco'}, $ovw_results{$cyc}{$ko}{'totod'}, $ovw_results{$cyc}{$ko}{'percod'}, $ovw_results{$cyc}{$ko}{'totgd'}, $ovw_results{$cyc}{$ko}{'percgd'},  $ovw_results{$cyc}{$ko}{'description'});
               printf $totovw_tsv_fh ("%s\t%s\t%.2f\t%s\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%s\n",$ko, $ovw_results{$cyc}{$ko}{'totg'}, $ovw_results{$cyc}{$ko}{'percg'}, $ovw_results{$cyc}{$ko}{'toto'}, $ovw_results{ $cyc}{$ko}{'perco'}, $ovw_results{$cyc}{$ko}{'totod'}, $ovw_results{$cyc}{$ko}{'percod'}, $ovw_results{$cyc}{$ko}{'totgd'}, $ovw_results{$cyc}{$ko}{'percgd'},  $ovw_results{$cyc}{$ko}{'description'});
               print {$krona_pd_fh} "$ovw_results{$cyc}{$ko}{'totgd'}\t$cyc\t$ko\t$ovw_results{$cyc}{$ko}{'description'}\n";
               print {$krona_g_fh} "$ovw_results{$cyc}{$ko}{'totg'}\t$cyc\t$ko\t$ovw_results{$cyc}{$ko}{'description'}\n";
               print {$krona_o_fh} "$ovw_results{$cyc}{$ko}{'toto'}\t$cyc\t$ko\t$ovw_results{$cyc}{$ko}{'description'}\n";
               print {$krona_d_fh} "$ovw_results{$cyc}{$ko}{'totod'}\t$cyc\t$ko\t$ovw_results{$cyc}{$ko}{'description'}\n";
            }
         }
      }
   }
}
###
open my $totmod_fh , '>', "$dir/mod.tsv" or die $!;

for my $cyc (sort keys %modulesKO){ 
   my $selector = $cyc;
   unless ($selector =~ 'meta.nonkey'){
      print {$totmod_fh} "$cyc\t\t";
      printf $totmod_fh ("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", "Gene", "N# genes", '%genes', "N# org",'%Org' ,"O-Depth", '%O-Depth',"G-Depth",'%G-Depth', "Full description");
   }
   for my $module  (keys %{$modulesKO { $cyc }} ) {
      print {$totmod_fh} "\t$module\n";
      for my $ko (keys %{$modulesKO { $cyc }{ $module }} ) {
         if ($total_gsum != 0){
            if ($nozero) { 
               if ($modulesKO{$cyc}{$module}{$ko}{'totg'} >= 1) {
                  unless ($selector =~ 'meta.nonkey'){
                     printf $totmod_fh ("\t\t%s\t%s\t%.2f\t%s\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%s\n",$ko, $modulesKO{$cyc}{$module}{$ko}{'totg'}, $modulesKO{$cyc}{$module}{$ko}{'percg'}, $modulesKO{$cyc}{$module}{$ko}{'toto'}, $modulesKO{$cyc}{$module}{$ko}{'perco'}, $modulesKO{$cyc}{$module}{$ko}{'totod'}, $modulesKO{$cyc}{$module}{$ko}{'percod'}, $modulesKO{$cyc}{$module}{$ko}{'totgd'}, $modulesKO{$cyc}{$module}{$ko}{'percgd'},  $modulesKO{$cyc}{$module}{$ko}{'description'});
                     print {$krona_mod_pd_fh} "$modulesKO{$cyc}{$module}{$ko}{'totgd'}\t$cyc\t$module\t$modulesKO{$cyc}{$module}{$ko}{'description'}\n";
                     print {$krona_mod_g_fh} "$modulesKO{$cyc}{$module}{$ko}{'totg'}\t$cyc\t$module\t$modulesKO{$cyc}{$module}{$ko}{'description'}\n";
                     print {$krona_mod_o_fh} "$modulesKO{$cyc}{$module}{$ko}{'toto'}\t$cyc\t$module\t$modulesKO{$cyc}{$module}{$ko}{'description'}\n";
                     print {$krona_mod_d_fh} "$modulesKO{$cyc}{$module}{$ko}{'totod'}\t$cyc\t$module\t$modulesKO{$cyc}{$module}{$ko}{'description'}\n";
                  }
               }
            }
            else {
               unless ($selector =~ 'meta.nonkey'){
                  printf $totmod_fh ("\t\t%s\t%s\t%.2f\t%s\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%s\n",$ko, $modulesKO{$cyc}{$module}{$ko}{'totg'}, $modulesKO{$cyc}{$module}{$ko}{'percg'}, $modulesKO{$cyc}{$module}{$ko}{'toto'}, $modulesKO{$cyc}{$module}{$ko}{'perco'}, $modulesKO{$cyc}{$module}{$ko}{'totod'}, $modulesKO{$cyc}{$module}{$ko}{'percod'}, $modulesKO{$cyc}{$module}{$ko}{'totgd'}, $modulesKO{$cyc}{$module}{$ko}{'percgd'},  $modulesKO{$cyc}{$module}{$ko}{'description'});
                  print {$krona_mod_pd_fh} "$modulesKO{$cyc}{$module}{$ko}{'totgd'}\t$cyc\t$module\t$modulesKO{$cyc}{$module}{$ko}{'description'}\n";
                  print {$krona_mod_g_fh} "$modulesKO{$cyc}{$module}{$ko}{'totg'}\t$cyc\t$module\t$modulesKO{$cyc}{$module}{$ko}{'description'}\n";
                  print {$krona_mod_o_fh} "$modulesKO{$cyc}{$module}{$ko}{'toto'}\t$cyc\t$module\t$modulesKO{$cyc}{$module}{$ko}{'description'}\n";
                  print {$krona_mod_d_fh} "$modulesKO{$cyc}{$module}{$ko}{'totod'}\t$cyc\t$module\t$modulesKO{$cyc}{$module}{$ko}{'description'}\n";
               }
            }
         }
      }
   }
   print {$totmod_fh} "\n";#here would come the sum of each of all the modules
}

###test

open my $totproc_tsv_fh , '>', "$dir/proc.tsv" or die $!;

for my $cyc (sort keys %processKO){
   my $selector = $cyc;
   unless ($selector =~ 'meta.nonkey'){
      print {$totproc_tsv_fh} "$cyc\t\t";
      printf $totproc_tsv_fh  ("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", "Gene", "N# genes", '%genes', "N# org",'%Org' ,"O-Depth", '%O-Depth',"G-Depth",'%G-Depth', "Full description");
   }
   for my $proc  (sort keys %{$processKO { $cyc }} ) {
      print {$totproc_tsv_fh} "\t$proc\n";
      for my $ko (keys %{$processKO { $cyc }{ $proc }} ) {
         if ($total_gsum != 0){
            if ($nozero) { 
               if ($processKO{$cyc}{$proc}{$ko}{'totg'} >= 1) {
                  unless ($selector =~ 'meta.nonkey'){
                     printf $totproc_tsv_fh ("\t\t%s\t%s\t%.2f\t%s\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%s\n",$ko, $processKO{$cyc}{$proc}{$ko}{'totg'}, $processKO{$cyc}{$proc}{$ko}{'percg'}, $processKO{$cyc}{$proc}{$ko}{'toto'}, $processKO{$cyc}{$proc}{$ko}{'perco'}, $processKO{$cyc}{$proc}{$ko}{'totod'}, $processKO{$cyc}{$proc}{$ko}{'percod'}, $processKO{$cyc}{$proc}{$ko}{'totgd'}, $processKO{$cyc}{$proc}{$ko}{'percgd'},  $processKO{$cyc}{$proc}{$ko}{'description'});
                     print {$krona_proc_pd_fh} "$processKO{$cyc}{$proc}{$ko}{'totgd'}\t$cyc\t$proc\t$ko\t$processKO{$cyc}{$proc}{$ko}{'description'}\n";
                     print {$krona_proc_g_fh} "$processKO{$cyc}{$proc}{$ko}{'totg'}\t$cyc\t$proc\t$ko\t$processKO{$cyc}{$proc}{$ko}{'description'}\n";
                     print {$krona_proc_o_fh} "$processKO{$cyc}{$proc}{$ko}{'toto'}\t$cyc\t$proc\t$ko\t$processKO{$cyc}{$proc}{$ko}{'description'}\n";
                     print {$krona_proc_d_fh} "$processKO{$cyc}{$proc}{$ko}{'totod'}\t$cyc\t$proc\t$ko\t$processKO{$cyc}{$proc}{$ko}{'description'}\n";
                  }
               }
            }
            else {
               unless ($selector =~ 'meta.nonkey'){
                  printf $totproc_tsv_fh ("\t\t%s\t%s\t%.2f\t%s\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%s\n",$ko, $processKO{$cyc}{$proc}{$ko}{'totg'}, $processKO{$cyc}{$proc}{$ko}{'percg'}, $processKO{$cyc}{$proc}{$ko}{'toto'}, $processKO{$cyc}{$proc}{$ko}{'perco'}, $processKO{$cyc}{$proc}{$ko}{'totod'}, $processKO{$cyc}{$proc}{$ko}{'percod'}, $processKO{$cyc}{$proc}{$ko}{'totgd'}, $processKO{$cyc}{$proc}{$ko}{'percgd'},  $processKO{$cyc}{$proc}{$ko}{'description'});
                  print {$krona_proc_pd_fh} "$processKO{$cyc}{$proc}{$ko}{'totgd'}\t$cyc\t$proc\t$ko\t$processKO{$cyc}{$proc}{$ko}{'description'}\n";
                  print {$krona_proc_g_fh} "$processKO{$cyc}{$proc}{$ko}{'totg'}\t$cyc\t$proc\t$ko\t$processKO{$cyc}{$proc}{$ko}{'description'}\n";
                  print {$krona_proc_o_fh} "$processKO{$cyc}{$proc}{$ko}{'toto'}\t$cyc\t$proc\t$ko\t$processKO{$cyc}{$proc}{$ko}{'description'}\n";
                  print {$krona_proc_d_fh} "$processKO{$cyc}{$proc}{$ko}{'totod'}\t$cyc\t$proc\t$ko\t$processKO{$cyc}{$proc}{$ko}{'description'}\n";
               }
            }
         }
      }
   }
   print {$totproc_tsv_fh} "\n";#here would come the sum of each of all the processes
}

delfile ("$dir/analyzedfastas.txt", "$dir/gensum.txt", "$dir/gendepthsum.txt", "$dir/orgdepthsum.txt", "$dir/keggsum.txt", "$dir/file_hash.txt"); #remove this line and the restore option can be used as a append option

#creating Krona files:

runcmd("ktImportText \Q$dir/krona.g.tsv,Genes\E \Q$dir/krona.gd.tsv,Gene Depth\E \Q$dir/krona.o.tsv,Organisms\E \Q$dir/krona.od.tsv,Organism Depth\E \Q$dir/krona.mod.g.tsv,Modules Genes\E \Q$dir/krona.mod.gd.tsv,Modules Gene Depth\E \Q$dir/krona.mod.o.tsv,Modules Organisms\E \Q$dir/krona.mod.od.tsv,Modules Organism Depth\E \Q$dir/krona.proc.g.tsv,Process Genes\E \Q$dir/krona.proc.gd.tsv,Process Gene Depth\E \Q$dir/krona.proc.o.tsv,Process Organisms\E \Q$dir/krona.proc.od.tsv,Process Organism Depth\E -o $dir/krona.html");

delfile ("$dir/krona.mod.g.tsv", "$dir/krona.mod.o.tsv", "$dir/krona.mod.od.tsv", "$dir/krona.mod.gd.tsv", "$dir/krona.proc.g.tsv", "$dir/krona.proc.o.tsv", "$dir/krona.proc.od.tsv", "$dir/krona.proc.gd.tsv", "$dir/krona.o.tsv", "$dir/krona.od.tsv", "$dir/krona.gd.tsv", "$dir/krona.g.tsv");

# to iterate over multiple directories, go the root/base/home directory from where your subdirs are and run : for dir in */; do  cd $dir && /path/fastaextract_desc_id.pl *.all.faa methanol ../methanol &&  cd ..; done
#change methanol to whatever you want, this is your search term. The second one is the name of the fasta in the base directory

__DATA__
