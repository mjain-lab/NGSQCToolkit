#! /usr/bin/perl

use strict;
use warnings;
use Getopt::Long;
use File::Basename;

# Parameter variables
my @files;
my $helpAsked;
my $homoPolyLen	= 8;
my $outFolder = "";
my $minReadLen = 100;

GetOptions(
			"i=s{1,2}" => \@files,
			"h|help" => \$helpAsked,
			"l|minReadLen=i" => \$minReadLen,
			"n|homoPolyLen=i" => \$homoPolyLen,
 			"o|outputFolder=s" => \$outFolder,
		  );
if(defined($helpAsked)) {
	prtUsage();
	exit;
}
if(@files == 0) {
	prtError("No input files are provided");
}

# Variables
my $seqFile = $files[0];
my $qualFile;
$qualFile = $files[1] if(scalar(@files) == 2);
my $outSeqFile;
my $outQualFile;
my $prevFastaSeqId = "";
my $fastaSeqId = "";
my $fastaSeq = "";
my $qualSeqId = "";
my $qualSeq = "";
my $seqCount = 0;
my $trimCount = 0;
my $trashCount = 0;

	my ($seqFileName, $filePath) = fileparse($seqFile);
	my ($qualFileName) = fileparse($qualFile) if(scalar(@files) == 2);
	$outFolder = $filePath if($outFolder eq "");
	$outFolder .= "/" if($outFolder !~ /\/$/);
	if(! -e $outFolder) {
		mkdir($outFolder) or die "Can not create output folder: $outFolder\n";
	}
	$outSeqFile = $outFolder . $seqFileName . "_trimmed";
	$outQualFile = $outFolder . $qualFileName . "_trimmed" if(scalar(@files) == 2);

	open(I, "<$seqFile") or die "Can not open file: $seqFile\n"; 
	open(Q, "<$qualFile") or die "Can not open file: $qualFile\n" if(scalar(@files) == 2);
	open(OI, ">$outSeqFile") or die "Can not open file: $outSeqFile\n"; 
	open(OQ, ">$outQualFile") or die "Can not open file: $outQualFile\n" if(scalar(@files) == 2);

if(scalar(@files) == 2) { 
	while(my $line = <I>) {
		chomp $line;
		my $qualLine = <Q>;
		chomp($qualLine);
		if($line =~ /^>/) {
			$seqCount++;
			$prevFastaSeqId = $fastaSeqId;
			$fastaSeqId = $line;
			$qualSeqId = $qualLine;
			if($fastaSeqId ne $qualSeqId) {
				print STDERR "Error: Read Id doesn't match in sequence and quality file for read number $seqCount in sequence file.\n";
				exit(-1);
			}
			if($fastaSeq ne "") {
				processSeq();
			}
			$fastaSeq = "";
			$qualSeq = "";
		}
		else {
			$fastaSeq .= $line;
			$qualSeq .= $qualLine . " ";
		}
	}
	if($fastaSeq ne "") {
		$prevFastaSeqId = $fastaSeqId;
		processSeq();
	}
}
else {
	while(my $line = <I>) {
		chomp $line;
		if($line =~ /^>/) {
			$seqCount++;
			$prevFastaSeqId = $fastaSeqId;
			$fastaSeqId = $line;
			if($fastaSeq ne "") {
				processSeq();
			}
			$fastaSeq = "";
		}
		else {
			$fastaSeq .= $line;
		}
	}
	if($fastaSeq ne "") {
		$prevFastaSeqId = $fastaSeqId;
		processSeq();
	}
}

print "Number of reads/sequences trashed with length < $minReadLen: $trashCount\n";
print "Number of reads/sequences trimmed containing homopolymer: $trimCount\n";
print "Trimmed read/sequence file: $outSeqFile\n";
print "Trimmed quality file: $outQualFile\n" if(scalar(@files) == 2);
exit;


sub processSeq {
	if(length $fastaSeq < $minReadLen) {
		$trashCount++;
		return;
	}
	if($homoPolyLen != 0) {
		if(hasPolyChar(\$fastaSeq)) {
			$trimCount++;
			$qualSeq = trimQualSeq($qualSeq, length $fastaSeq) if(scalar(@files) == 2);
		}
	}
	if(length $fastaSeq < $minReadLen) {
		$trashCount++;
		return;
	}
	print OI "$prevFastaSeqId\n";
	print OI formatSeq($fastaSeq), "\n";
	print OQ "$prevFastaSeqId\n" if(scalar(@files) == 2);
	print OQ formatQualSeq($qualSeq), "\n" if(scalar(@files) == 2);
}

sub hasPolyChar {
	my $seqRef = $_[0];
	my $flag = 0;
	if($$seqRef =~ s/(A{$homoPolyLen,}).*//i) {
		$flag = 1;
	}
	if($$seqRef =~ s/(T{$homoPolyLen,}).*//i) {
		$flag = 1;
	}
	if($$seqRef =~ s/(G{$homoPolyLen,}).*//i) {
		$flag = 1;
	}
	if($$seqRef =~ s/(C{$homoPolyLen,}).*//i) {
		$flag = 1;
	}
	return $flag;
}

sub trimQualSeq {
	my $qualSeq = $_[0];
	my $seqLen = $_[1];
	$qualSeq =~ /^((\d{1,2}\s+){$seqLen})/;
	my $trimmedQualSeq = $1;
	$trimmedQualSeq =~ s/\s+$//;
	return $trimmedQualSeq;
}

sub formatSeq {
	my $seq = $_[0];
    my $newSeq = "";
    my $ch = 60;
    for(my $i=0; $i<length $seq; $i+=$ch) {
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
	my @arr = split(" ", $qualSeq);
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

sub prtHelp {
	print "\n$0 options:\n\n";
	print "### Input reads/sequences (FASTA format; .fna and .qual files) (Required)\n";
	print "  -i <Read/Sequence file> [Quality file (optional)]\n";
	print "    Read/Sequence and quality file in FASTA format\n";
	print "\n";
	print "### Other options [Optional]\n";
	print "  -h | -help\n";
	print "    Prints this help\n";
	print "  -l | -minReadLen <Integer>\n";
	print "    Minimum length of a read/sequence to be retained in output\n";
	print "    default: 100\n";
	print "  -n | -homoPolyLen <Integer>\n";
	print "    Minimum length of the homopolymer to be trimmed\n";
	print "      For eg.: -n 8, will trim the right end of read/sequence from the homopolymer of at least 8 bases long\n";
	print "      Note:- use -n 0 to skip homopolymer trimming (for only length filtering)\n";
	print "    default: 8\n";
	print "  -o | -outputFolder <Output folder name/path>\n";
	print "    Output will be stored in the given folder\n";
	print "    default: By default, files will be stored where the input files are\n";
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
