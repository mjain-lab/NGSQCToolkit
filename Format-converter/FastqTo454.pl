#! /usr/bin/perl

use strict;
use warnings;
use List::Util qw(sum min max);
use Getopt::Long;
use File::Basename;

# Parameter variables
my $file;
my $helpAsked;
my $outFolder = "";
my $subVal;

my $seqFormat = "a";		# 1: Sanger; 2: Solexa; 3: Illumina 1.3+; 4: Illumina 1.5+;

GetOptions(
			"i=s" => \$file,
			"h|help" => \$helpAsked,
			"o|outputFolder=s" => \$outFolder,
			"v|fastqVariant=s" => \$seqFormat,
		  );
if(defined($helpAsked)) {
	prtUsage();
	exit;
}
if(!defined($file)) {
	prtError("No input files are provided");
}

my ($fileName, $filePath) = fileparse($file);
$outFolder = $filePath if($outFolder eq "");
$outFolder .= "/" if($outFolder !~ /\/$/);
if(! -e $outFolder) {
	mkdir($outFolder) or die "Can not create output folder: $outFolder\n";
}

if($seqFormat =~ /a/i) {
	print "Checking FASTQ format: File $file...\n";
	my $nLines = checkFastQFormat($file, 1);
}
if($seqFormat == 1) {
	$subVal = 33;
	print "Input FASTQ file format: Sanger\n";
}
if($seqFormat == 2) {
	$subVal = 64;
	print "Input FASTQ file format: Solexa\n";
}
if($seqFormat == 3) {
	$subVal = 64;
	print "Input FASTQ file format: Illumina 1.3+\n";
}
if($seqFormat == 4) {
	$subVal = 64;
	print "Input FASTQ file format: Illumina 1.5+\n";
}
if($seqFormat == 5) {
	$subVal = 33;
	print "Input FASTQ file format: Illumina 1.8+\n";
}


my $outFnaFile = $outFolder . $fileName . "_fna";
my $outQualFile = $outFolder . $fileName . "_qual";

open(I, "<$file") or die "Can not open file: $file\n";
open(OF, ">$outFnaFile") or die "Can not open file: $outFnaFile\n";
open(OQ, ">$outQualFile") or die "Can not open file: $outQualFile\n";


while(my $line = <I>) {
	chomp($line);
	my $id = $line;
	$id =~ s/^\@//;
	print OF ">$id\n";
	my $seq = <I>;
	chomp $seq;
	print OF formatSeq($seq), "\n";
	<I>;
	print OQ ">$id\n";
	my $qualLine = <I>;
	chomp($qualLine);
	print OQ &IlluToPhred($qualLine), "\n";
}

exit;

sub IlluToPhred {
	my $qualLine = $_[0];
	my $retQualLine = "";
	my @ASCII = unpack("C*", $qualLine);
	my @newASCII = ();
	foreach my $val (@ASCII) {
		my $newVal = $val - $subVal;
		$retQualLine .= $newVal . " ";
	}
	chop $retQualLine;
	return formatQualSeq($retQualLine);
}

sub prtHelp {
	print "\n$0 options:\n\n";
	print "### Input reads (FASTQ) (Required)\n";
	print "  -i <Illumina FASTQ read file>\n";
	print "    Read file in Illumina FASTQ format\n";
	print "\n";
	print "### Other options [Optional]\n";
	print "  -h | -help\n";
	print "    Prints this help\n";
	print "  -o | -outputFolder <Output folder name>\n";
	print "    Output will be stored in the given folder\n";
	print "    default: By default, files will be stored where the input file is\n";
	print "  -v | -fastqVariant <FASTQ variant>\n";
	print "    FASTQ variants:\n";
	print "      1 = Sanger (Phred+33, 33 to 73)\n";
	print "      2 = Solexa (Phred+64, 59 to 104)\n";
	print "      3 = Illumina (1.3+) (Phred+64, 64 to 104)\n";
	print "      4 = Illumina (1.5+) (Phred+64, 66 to 104)\n";
	print "      5 = Illumina (1.8+) (Phred+33, 33 to 74)\n";
	print "      A = Automatic detection of FASTQ variant\n";
	print "    default: \"A\"\n";
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

sub formatSeq {
	my $seq = $_[0];
    my $newSeq = "";
    my $ch = 60;
    my $len = length $seq;
    for(my $i=0; $i<$len; $i+=$ch) {
        $newSeq .= substr($seq, $i, $ch) . "\n";
    }
    chomp($newSeq);         # To remove \n at the end of the whole sequence..
    return $newSeq;
}

sub formatQualSeq {
	my $qualSeq = $_[0];
	my $fQSeq = "";
	my $ch = 60;
	my $valCount = 0;
	my @arr = split(/\s+/, $qualSeq);
	for(my $i=0; $i<@arr; $i++) {
		$valCount++;
		if($valCount % $ch == 0) {
			$fQSeq .= $arr[$i] . "\n";
		}
		else {
			$fQSeq .= $arr[$i] . " ";
		}
	}
	$fQSeq =~ s/\s+$//;
	return $fQSeq;
}

sub checkFastQFormat {				# Takes FASTQ file as an input and if the format is incorrect it will print error and exit, otherwise it will return the number of lines in the file.
	my $file = $_[0];
	my $isVariantIdntfcntOn = $_[1];
	my $lines = 0;
	open(F, "<$file") or die "Can not open file $file\n";
	my $counter = 0;
	my $minVal = 1000;
	my $maxVal = 0;
	while(my $line = <F>) {
		$lines++;
		$counter++;
		next if($line =~ /^\n$/);
		if($counter == 1 && $line !~ /^\@/) {
			prtErrorExit("Invalid FASTQ file format.\n\t\tFile: $file");
		}
		if($counter == 3 && $line !~ /^\+/) {
			prtErrorExit("Invalid FASTQ file format.\n\t\tFile: $file");
		}
		if($counter == 4 && $lines < 1000000) {
			chomp $line;
			my @ASCII = unpack("C*", $line);
			$minVal = min(min(@ASCII), $minVal);
			$maxVal = max(max(@ASCII), $maxVal);
		}
		if($counter == 4) {
			$counter = 0;
		}
	}
	close(F);
	my $tseqFormat = 0;
	if($minVal >= 33 && $minVal <= 73 && $maxVal >= 33 && $maxVal <= 73) {
		$tseqFormat = 1;
	}
	elsif($minVal >= 66 && $minVal <= 105 && $maxVal >= 66 && $maxVal <= 105) {
		$tseqFormat = 4;			# Illumina 1.5+
	}
	elsif($minVal >= 64 && $minVal <= 105 && $maxVal >= 64 && $maxVal <= 105) {
		$tseqFormat = 3;			# Illumina 1.3+
	}
	elsif($minVal >= 59 && $minVal <= 105 && $maxVal >= 59 && $maxVal <= 105) {
		$tseqFormat = 2;			# Solexa
	}
	elsif($minVal >= 33 && $minVal <= 74 && $maxVal >= 33 && $maxVal <= 74) {
		$tseqFormat = 5;			# Illumina 1.8+
	}
	if($isVariantIdntfcntOn) {
		$seqFormat = $tseqFormat;
	}
	else {
		if($tseqFormat != $seqFormat) {
			print STDERR "Warning: It seems the specified variant of FASTQ doesn't match the quality values in input FASTQ files.\n";
		}
	}
	return $lines;
}

