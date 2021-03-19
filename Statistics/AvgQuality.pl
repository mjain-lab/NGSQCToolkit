#! /usr/bin/perl

use strict;
use warnings;
use List::Util qw(sum min max);
use Getopt::Long;
use File::Basename;

# Parameter variables
my $file;
my $helpAsked;
my $outFile = "";

GetOptions(
			"i=s" => \$file,
			"h|help" => \$helpAsked,
			"o|outputFile=s" => \$outFile,
		  );
if(defined($helpAsked)) {
	prtUsage();
	exit;
}
if(!defined($file)) {
	prtError("No input files are provided");
}

my ($fileName, $filePath) = fileparse($file);
$outFile = $file . "_qual_stat" if($outFile eq "");

open(I, "<$file") or die "Can not open file: $file\n";
open(O, ">$outFile") or die "Can not open file: $outFile\n";


my $prevFastaSeqId = "";
my $fastaSeqId = "";
my $fastaSeq = "";
my $seqCount = 0;
my $ttlQual = 0;

while(my $line = <I>) {
	chomp $line;
	if($line =~ /^>/) {
		$line =~ s/^>//;
		$prevFastaSeqId = $fastaSeqId;
		$fastaSeqId = $line;
		if($fastaSeq ne "") {
			prtQuality($prevFastaSeqId, $fastaSeq);
		}
		$fastaSeq = "";
	}
	else {
		$fastaSeq .= $line . " ";
	}
}
if($fastaSeq ne "") {
	$prevFastaSeqId = $fastaSeqId;
	prtQuality($prevFastaSeqId, $fastaSeq);
}

printf O "Final quality average: %0.2f\n", $ttlQual/$seqCount;
close(O);
close(I);

exit;


sub prtQuality {
	my $id = $_[0];
	my $qualStr = $_[1];
	chop $qualStr;
	$seqCount++;
	my @qVal = split(/\s+/, $qualStr);
	my $sum = sum(@qVal);
	my $qual = sprintf "%0.2f", $sum/(scalar @qVal);
	$ttlQual += $qual;
	print O "$id\t$qual\n";	
}


sub prtHelp {
	print "\n$0 options:\n\n";
	print "### Input quality (FASTA) (Required)\n";
	print "  -i <Quality file>\n";
	print "    Quality file in FASTA format\n";
	print "\n";
	print "### Other options [Optional]\n";
	print "  -h | -help\n";
	print "    Prints this help\n";
	print "  -o | -outputFile <Output file name>\n";
	print "    Output will be stored in the given file\n";
	print "    default: By default, quality statistics file will be stored where the input file is\n";
	print "\n";
}

sub prtError {
	my $msg = $_[0];
	print STDERR "+======================================================================+\n";
	printf STDERR "|%-70s|\n", "  Error:";
	printf STDERR "|%-70s|\n", "       $msg";
	print STDERR "+======================================================================+\n";
	prtUsage();
	exit;
}

sub prtUsage {
	print "\nUsage: perl $0 <options>\n";
	prtHelp();
}
