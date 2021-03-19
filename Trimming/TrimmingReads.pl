#! /usr/bin/perl

use strict;
use warnings;
use Getopt::Long;
use File::Basename;
use List::Util qw(sum min max);

my $seqFormat = "a";		# 1: Sanger; 2: Solexa; 3: Illumina 1.3+; 4: Illumina 1.5+;
my $subVal;
my $subVal2;

# Parameter variables
my $file;
my $file2;
my $helpAsked;
my $rTrimBases = 0;
my $lTrimBases = 0;
my $qCutOff = 0;
my $lenCutOff = -1;
my $outFile = "";
my $isQualTrimming = 1;

GetOptions(
			"i=s" => \$file,
			"irev=s" => \$file2,
			"h|help" => \$helpAsked,
			"l|leftTrimBases=i" => \$lTrimBases,
			"o|outputFile=s" => \$outFile,
			"r|rightTrimBases=i" => \$rTrimBases,
			"q|qualCutOff=i" => \$qCutOff,
			"n|lenCutOff=i" => \$lenCutOff,
		  );
if(defined($helpAsked)) {
	prtUsage();
	exit;
}
if(!defined($file)) {
	prtError("No input files are provided");
}

if($file2) {
	$outFile = $file . "_trimmed";
	my $outFile2 = $file2 . "_trimmed";
	open(I1, "$file") or die "Can not open file $file\n";
	open(I2, "$file2") or die "Can not open file $file2\n";
	open(O1, ">$outFile") or die "Can not create file $outFile\n";
	open(O2, ">$outFile2") or die "Can not create file $outFile2\n";
	my $tmpLine = <I1>;
	close(I1);
	if($tmpLine =~ /^@/) {
		print "Input read/sequence format: FASTQ (Paired-end)\n";
		if($lTrimBases == 0 && $rTrimBases == 0 && $qCutOff == 0 && $lenCutOff == -1) {
			print "Trimming parameters are not set.\nNothing to do.\nExiting...\n";
			unlink($outFile);
			unlink($outFile2);
			exit;
		}
		print "Checking FASTQ variant: File $file...\n";
		my $nLines = checkFastQFormat($file, 1);
		$subVal = getSubVal($seqFormat);

		print "Checking FASTQ variant: File $file2...\n";
		my $nLines2 = checkFastQFormat($file2, 1);
		$subVal2 = getSubVal($seqFormat);

		if($nLines != $nLines2) {
			prtErrorExit("Number of reads in paired-end data files are not same.\n\t\tFiles: $file, $file2");
		}
		
		if($subVal != $subVal2) {
			prtErrorExit("FASTQ variant of paired-end data files are not same.\n\t\tFiles: $file, $file2");
		}

		if($lTrimBases != 0 || $rTrimBases != 0) {
			print "Trimming $lTrimBases bases from left end and $rTrimBases bases from right end";
			$isQualTrimming = 0;
		}
		else {
			print "Trimming based on PHRED quality score (< $qCutOff)";
		}
		print " followed by length filtering (< $lenCutOff bp)\n";

		open(I1, "$file") or die "Can not open file $file\n";
		my $c = 0;
		my $currId1 = "";
		my $currId2 = "";
		my $currSeq = "";
		my $currQual = "";
		my $curr2Id1 = "";
		my $curr2Id2 = "";
		my $currSeq2 = "";
		my $currQual2 = "";
		
		
		while(my $line = <I1>) {
			my $line2 = <I2>;
			chomp $line;
			chomp $line2;
	        $c++;
	        if($c == 5) {
	                $c = 1;
	        }
	        if($c == 1) {
	        	$currId1 = $line;
	        	$curr2Id1 = $line2;
	        }
	        if($c == 3) {
	        	$currId2 = $line;
	        	$curr2Id2 = $line2;
	        }
	        if($isQualTrimming == 0) {
		        if($c == 2) {
		        	$currSeq = trimSeq($line);
		        	$currSeq2 = trimSeq($line2);
		        }
		        elsif($c == 4) {
		        	$currQual = trimSeq($line);
		        	$currQual2 = trimSeq($line2);
		        	print O1 "$currId1\n$currSeq\n$currId2\n$currQual\n" if((length $currSeq >= $lenCutOff) && (length $currSeq2 >= $lenCutOff));
		        	print O2 "$curr2Id1\n$currSeq2\n$curr2Id2\n$currQual2\n" if((length $currSeq >= $lenCutOff) && (length $currSeq2 >= $lenCutOff));
		        }
	        }
	        else {
		        if($c == 4) {
					$currQual = trimSeq4Qual($line);
					$currQual2 = trimSeq4Qual($line2);
					if(defined($currQual) && defined($currQual2)) {
						my $len = length $currQual;
						my $len2 = length $currQual2;
						$currSeq =~ /^(.{$len})/;
						$currSeq = $1;
						$currSeq2 =~ /^(.{$len2})/;
						$currSeq2 = $1;
			        	print O1 "$currId1\n$currSeq\n$currId2\n$currQual\n" if((length $currSeq >= $lenCutOff) && (length $currSeq2 >= $lenCutOff));
			        	print O2 "$curr2Id1\n$currSeq2\n$curr2Id2\n$currQual2\n" if((length $currSeq >= $lenCutOff) && (length $currSeq2 >= $lenCutOff));
					}
		        }
		        elsif($c == 2) {
		        	$currSeq = $line;
		        	$currSeq2 = $line2;
		        }
	        }
		}
		print "Trimmed files are generated: $outFile $outFile2\n";
	}
	else {
		print "Error:::\n\tPaired-end sequeneing data need to be in FASTQ format\n";
		exit;
	}
	close(O2);
	close(O1);
	close(I2);
	close(I1);
}
else {
	$outFile = $file . "_trimmed" if($outFile eq "");
	
	open(I, "$file") or die "Can not open file $file\n";
	open(O, ">$outFile") or die "Can not create file $outFile\n";
	my $tmpLine = <I>;
	close(I);
	if($tmpLine =~ /^@/) {
		print "Input read/sequence format: FASTQ\n";
		if($lTrimBases == 0 && $rTrimBases == 0 && $qCutOff == 0 && $lenCutOff == -1) {
			print "Trimming parameters are not set.\nNothing to do.\nExiting...\n";
			unlink($outFile);
			exit;
		}
		print "Checking FASTQ variant: File $file...\n";
		my $nLines = checkFastQFormat($file, 1);
		$subVal = getSubVal($seqFormat);


		if($lTrimBases != 0 || $rTrimBases != 0) {
			print "Trimming $lTrimBases bases from left end and $rTrimBases bases from right end";
			$isQualTrimming = 0;
		}
		else {
			print "Trimming based on PHRED quality score (< $qCutOff)";
		}
		print " followed by length filtering (< $lenCutOff bp)\n";
	
		open(I, "$file") or die "Can not open file $file\n";	
		my $c = 0;
		my $currId1 = "";
		my $currId2 = "";
		my $currSeq = "";
		my $currQual = "";
		
		
		while(my $line = <I>) {
			chomp $line;
	        $c++;
	        if($c == 5) {
	                $c = 1;
	        }
	        if($c == 1) {
	        	$currId1 = $line;
	        }
	        if($c == 3) {
	        	$currId2 = $line;
	        }
	        if($isQualTrimming == 0) {
		        if($c == 2) {
		        	$currSeq = trimSeq($line);
		        }
		        elsif($c == 4) {
		        	$currQual = trimSeq($line);
		        	print O "$currId1\n$currSeq\n$currId2\n$currQual\n" if(length $currSeq >= $lenCutOff);
		        }
	        }
	        else {
		        if($c == 4) {
					$currQual = trimSeq4Qual($line);
					if(defined($currQual)) {
						my $len = length $currQual;
						$currSeq =~ /^(.{$len})/;
						$currSeq = $1;
			        	print O "$currId1\n$currSeq\n$currId2\n$currQual\n" if(length $currSeq >= $lenCutOff);
					}
		        }
		        elsif($c == 2) {
		        	$currSeq = $line;
		        }
	        }
		}
		print "Trimmed file is generated: $outFile\n";
	}
	else {
		print "Input read/sequence format: FASTA\n";
		if($qCutOff != 0) {
			print "Warning: Quality trimming can not be performed for FASTA files\n";
			$qCutOff = 0;
		}
		if($lTrimBases == 0 && $rTrimBases == 0 && $lenCutOff == -1) {
			print "Trimming parameters are not set.\nNothing to do.\nExiting...\n";
			unlink($outFile);
			exit;
		}
		print "Trimming $lTrimBases bases from left end and $rTrimBases bases from right end followed by length filtering (< $lenCutOff bp)\n";
		open(I, "$file") or die "Can not open file $file\n";
		my $prevFastaSeqId = "";
		my $fastaSeqId = "";
		my $fastaSeq = "";
		
		while(my $line = <I>) {
			chomp $line;
			if($line =~ /^>/) {
				$prevFastaSeqId = $fastaSeqId;
				$fastaSeqId = $line;
				if($fastaSeq ne "") {
					my $outSeq = trimSeq($fastaSeq);
					print O "$prevFastaSeqId\n", formatSeq($outSeq), "\n" if(length $outSeq >= $lenCutOff);
				}
				$fastaSeq = "";
			}
			else {
				$fastaSeq .= $line;
			}
		}
		if($fastaSeq ne "") {
			$prevFastaSeqId = $fastaSeqId;
			my $outSeq = trimSeq($fastaSeq);
			print O "$prevFastaSeqId\n", formatSeq($outSeq), "\n" if(length $outSeq >= $lenCutOff);
		}
		print "Trimmed read/sequence file is generated: $outFile\n";
	}
	close(O);
	close(I);
}


sub trimSeq {
	my $seq = $_[0];
	$seq =~ /^.{$lTrimBases}(.+).{$rTrimBases}$/;
	return $1;
}

sub trimSeq4Qual {
	my $qual = $_[0];
	my @ASCII = unpack("C*", $qual);
	my $trimCount = 0;
	for(my $i=@ASCII; $i>0; $i--) {
		my $val = $ASCII[$i-1];
		$val -= $subVal;
		if($val < $qCutOff) {
			$trimCount++;
		}
		else {
			last;
		}
	}
	$qual =~ /^(.+).{$trimCount}$/;
	return $1;
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

sub prtHelp {
	print "\n$0 options:\n\n";
	print "### Input reads/sequences (FASTQ/FASTA) (Required)\n";
	print "  -i <Forward read/sequence file>\n";
	print "    File containing reads/sequences in either FASTQ or FASTA format\n";
	print "\n";
	print "### Input reads/sequences (FASTQ) [Optional]\n";
	print "  -irev <Reverse read/sequence file of paired-end data>\n";
	print "    File containing reverse reads/sequences of paired-end data in FASTQ format\n";
	print "\n";
	print "### Other options [Optional]\n";
	print "  -h | -help\n";
	print "    Prints this help\n";
	print "--------------------------------- Trimming Options ---------------------------------\n";
	print "  -l | -leftTrimBases <Integer>\n";
	print "    Number of bases to be trimmed from left end (5' end)\n";
	print "    default: 0\n";
	print "  -r | -rightTrimBases <Integer>\n";
	print "    Number of bases to be trimmed from right end (3' end)\n";
	print "    default: 0\n";
	print "  -q | -qualCutOff <Integer> (Only for FASTQ files)\n";
	print "    Cut-off PHRED quality score for trimming reads from right end (3' end)\n";
	print "      For eg.: -q 20, will trim bases having PHRED quality score less than 20 at 3' end of the read\n";
	print "      Note: Quality trimming can be performed only if -l and -r are not used\n";
	print "    default: 0 (i.e. quality trimming is OFF)\n";
	print "  -n | -lenCutOff <Integer>\n";
	print "    Read length cut-off\n";
	print "    Reads shorter than given length will be discarded\n";
	print "    default: -1 (i.e. length filtering is OFF)\n";
	print "--------------------------------- Output Options ---------------------------------\n";
	print "  -o | -outputFile <Output file name>\n";
	print "    Output will be stored in the given file\n";
	print "    default: By default, output file will be stored where the input file is\n";
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

sub prtErrorExit {
	my $errmsg = $_[0];
	print STDERR "Error:\t", $errmsg, "\n";
	exit;
}

sub prtUsage {
	print "\nUsage: perl $0 <options>\n";
	prtHelp();
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

sub getSubVal {
	my $seqFormat = $_[0];
	my $subVal = 0;
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
	return $subVal;
}