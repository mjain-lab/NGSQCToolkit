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
$outFile = $file . "_illu.fq" if($outFile eq "");

open(I, "<$file") or die "Can not open file: $file\n";
open(OF, ">$outFile") or die "Can not open file: $outFile\n";

while(my $line = <I>) {
	chomp($line);
	my $id = $line;
	$id =~ s/^\@//;
	print OF "\@$id\n";
	my $seq = <I>;
	print OF "$seq";
	<I>;
	print OF "+\n";
	my $qualLine = <I>;
	chomp($qualLine);
	print OF &IlluToSanger($qualLine), "\n";
}

exit;

sub IlluToSanger {
	my $qualLine = $_[0];
	my @ASCII = unpack("C*", $qualLine);
	my @newASCII = ();
	foreach my $val (@ASCII) {
		my $Q = sprintf "%0.0f", 10 * log(1 + 10 ** (($val - 64) / 10.0)) / log(10);
		push(@newASCII, $Q+64);
	}
	my $retQualLine = pack("C*", @newASCII);
	return $retQualLine;
}

sub prtHelp {
	print "\n$0 options:\n\n";
	print "### Input reads (FASTQ) (Required)\n";
	print "  -i <Solexa FASTQ read file>\n";
	print "    Read file in Solexa FASTQ format\n";
	print "\n";
	print "### Other options [Optional]\n";
	print "  -h | -help\n";
	print "    Prints this help\n";
	print "  -o | -outputFile <Output file name>\n";
	print "    Output will be stored in the given file\n";
	print "    default: By default, file will be stored where the input file is\n";
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
