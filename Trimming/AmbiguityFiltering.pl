#! /usr/bin/perl

use strict;
use warnings;
use Getopt::Long;
use File::Basename;
use List::Util qw(sum min max);

# Parameter variables
my $file;
my $file2;
my $helpAsked;
my $numNBases = -1;
my $perNBases = -1;
my $end5Trim;
my $end3Trim;
my $lenCutOff = -1;
my $outFile = "";
my $processingFlag = 0; # 1: N count filter, 2: N percent filter, 3: End N trim

GetOptions(
			"i=s" => \$file,
			"irev=s" => \$file2,
			"h|help" => \$helpAsked,
			"c|countN=i" => \$numNBases,
			"o|outputFile=s" => \$outFile,
			"p|percentN=i" => \$perNBases,
			"t5|trim5EndN" => \$end5Trim,
			"t3|trim3EndN" => \$end3Trim,
			"n|lenCutOff=i" => \$lenCutOff,
		  );
if(defined($helpAsked)) {
	prtUsage();
	exit;
}
if(!defined($file)) {
	prtError("No input files are provided");
}

### Validating input filtering options
if($numNBases == -1 && $perNBases == -1 && !defined($end5Trim) && !defined($end3Trim) && $lenCutOff == -1) {
	print "Filtering or trimming parameters are not set.\nNothing to do.\nExiting...\n";
	exit;
}
if($numNBases != -1) {
	$processingFlag = 1;
	print "Filtering out reads containing ambiguous base count > $numNBases\n";
}
elsif($perNBases != -1) {
	$processingFlag = 2;
	print "Filtering out reads containing ambiguous base percentage > $perNBases%\n";
}
else {
	$processingFlag = 3;
	print "Trimming ambiguous bases from end(s) of reads followed by length filtering (< $lenCutOff bp)\n";
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
		print "Checking FASTQ format: File $file...\n";
		my $nLines = checkFastQFile($file, 1);

		print "Checking FASTQ format: File $file2...\n";
		my $nLines2 = checkFastQFile($file2, 1);

		if($nLines != $nLines2) {
			prtErrorExit("Number of reads in paired-end data files are not same.\n\t\tFiles: $file, $file2");
		}

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
	        if($c == 2) {
		        $currSeq = $line;
		        $currSeq2 = $line2;	        	
	        }
	        if($c == 4) {
		        $currQual = $line;
		        $currQual2 = $line2;
		        if($processingFlag == 1) {
		        	my $nC1 = getNCount($currSeq);
		        	my $nC2 = getNCount($currSeq2);
		        	if($nC1 <= $numNBases && $nC2 <= $numNBases) {
			        	print O1 "$currId1\n$currSeq\n$currId2\n$currQual\n" if((length $currSeq >= $lenCutOff) && (length $currSeq2 >= $lenCutOff));
			        	print O2 "$curr2Id1\n$currSeq2\n$curr2Id2\n$currQual2\n" if((length $currSeq >= $lenCutOff) && (length $currSeq2 >= $lenCutOff));
		        	}
		        }
		        elsif($processingFlag == 2) {
		        	my $nP1 = getNPercent($currSeq);
		        	my $nP2 = getNPercent($currSeq2);
		        	if($nP1 <= $perNBases && $nP2 <= $perNBases) {
			        	print O1 "$currId1\n$currSeq\n$currId2\n$currQual\n" if((length $currSeq >= $lenCutOff) && (length $currSeq2 >= $lenCutOff));
			        	print O2 "$curr2Id1\n$currSeq2\n$curr2Id2\n$currQual2\n" if((length $currSeq >= $lenCutOff) && (length $currSeq2 >= $lenCutOff));
		        	}
		        }
		        elsif($processingFlag == 3) {
		        	($currSeq, $currQual) = trimNsAndLenFilter($currSeq, $currQual);
		        	($currSeq2, $currQual2) = trimNsAndLenFilter($currSeq2, $currQual2);
		        	if($currSeq ne "-1" && $currSeq2 ne "-1") {
			        	print O1 "$currId1\n$currSeq\n$currId2\n$currQual\n" if((length $currSeq >= $lenCutOff) && (length $currSeq2 >= $lenCutOff));
			        	print O2 "$curr2Id1\n$currSeq2\n$curr2Id2\n$currQual2\n" if((length $currSeq >= $lenCutOff) && (length $currSeq2 >= $lenCutOff));
		        	}
		        }
	        }
		}
		print "Filtered files are generated: $outFile $outFile2\n";
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
		print "Checking FASTQ variant: File $file...\n";
		my $nLines = checkFastQFile($file, 1);
	
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
	        if($c == 2) {
		        $currSeq = $line;
	        }
	        if($c == 4) {
		        $currQual = $line;
		        if($processingFlag == 1) {
		        	my $nC1 = getNCount($currSeq);
		        	if($nC1 <= $numNBases) {
		        		print O "$currId1\n$currSeq\n$currId2\n$currQual\n" if(length $currSeq >= $lenCutOff);
		        	}
		        }
		        elsif($processingFlag == 2) {
		        	my $nP1 = getNPercent($currSeq);
		        	if($nP1 <= $perNBases) {
			        	print O "$currId1\n$currSeq\n$currId2\n$currQual\n" if(length $currSeq >= $lenCutOff);
		        	}
		        }
		        elsif($processingFlag == 3) {
		        	($currSeq, $currQual) = trimNsAndLenFilter($currSeq, $currQual);
		        	if($currSeq ne "-1") {
			        	print O "$currId1\n$currSeq\n$currId2\n$currQual\n" if(length $currSeq >= $lenCutOff);
		        	}
		        }
	        }
		}
		print "Filtered file is generated: $outFile\n";
	}
	else {
		print "Input read/sequence format: FASTA\n";

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
					processFastaSeq($prevFastaSeqId, $fastaSeq);
				}
				$fastaSeq = "";
			}
			else {
				$fastaSeq .= $line;
			}
		}
		if($fastaSeq ne "") {
			$prevFastaSeqId = $fastaSeqId;
			processFastaSeq($prevFastaSeqId, $fastaSeq);
		}
		print "Filtered file is generated: $outFile\n";
	}
	close(O);
	close(I);
}

sub processFastaSeq {
	my ($prevFastaSeqId, $fastaSeq) = @_;
	if($processingFlag == 1) {
		my $nC1 = getNCount($fastaSeq);
		if($nC1 <= $numNBases) {
			print O "$prevFastaSeqId\n", formatSeq($fastaSeq), "\n" if(length $fastaSeq >= $lenCutOff);
		}
	}
	elsif($processingFlag == 2) {
		my $nP1 = getNPercent($fastaSeq);
		if($nP1 <= $perNBases) {
			print O "$prevFastaSeqId\n", formatSeq($fastaSeq), "\n" if(length $fastaSeq >= $lenCutOff);
		}
	}
	elsif($processingFlag == 3) {
		($fastaSeq) = trimNsAndLenFilter($fastaSeq);
		if($fastaSeq ne "-1") {
			print O "$prevFastaSeqId\n", formatSeq($fastaSeq), "\n" if(length $fastaSeq >= $lenCutOff);
		}
	}	
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
	print "  -c | -countN <Integer>\n";
	print "    Maximum number of allowed ambiguous bases\n";
	print "    default: 0\n";
	print "  -p | -percentN <Integer>\n";
	print "    Maximum percentage of allowed ambiguous bases\n";
	print "    default: 0\n";
	print "  -t5 | -trim5EndN\n";
	print "    Trim ambiguous bases from 5' end of the sequence\n";
	print "    default: off\n";
	print "  -t3 | -trim3EndN\n";
	print "    Trim ambiguous bases from 3' end of the sequence\n";
	print "    default: off\n";
	print "  -n | -lenCutOff <Integer>\n";
	print "    Sequence length cut-off\n";
	print "    Sequences shorter than given length will be discarded\n";
	print "    default: -1 (i.e. length filtering is OFF)\n";
	print " NOTE: filtering can be performed using any one of (-c), (-p) and (-t5 and/or -t3) switches at a time\n";
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

sub checkFastQFile {				# Takes FASTQ file as an input and if the format is incorrect it will print error and exit, otherwise it will return the number of lines in the file.
	my $file = $_[0];
	my $lines = 0;
	open(F, "<$file") or die "Can not open file $file\n";
	my $counter = 0;
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
		if($counter == 4) {
			$counter = 0;
		}
	}
	close(F);
	return $lines;
}

sub getNCount {
	my ($seq) = @_;
	my $len = length $seq;
	my $c = 0;
	while($seq =~ /[^ATGC]/ig){$c++;}
	return $c;
}

sub getNPercent {
	my ($seq) = @_;
	my $len = length $seq;
	my $c = 0;
	while($seq =~ /[^ATGC]/ig){$c++;}
	return $c/$len*100;
}

sub trimNsAndLenFilter {
	my ($seq, $qual) = @_;
	if($end3Trim) {
		if($seq =~ s/([^ATGC]+)$//) {
			if($qual) {
				my $nCount = length $1;
				$qual =~ s/[^\n]{$nCount}$//;
			}
		}
	}
	if($end5Trim) {
		if($seq =~ s/^([^ATGC]+)//) {
			if($qual) {
				my $nCount = length $1;
				$qual =~ s/^[^\n]{$nCount}//;
			}
		}
	}
	if($seq) {
		if((length $seq) >= $lenCutOff) {
			return ($seq, $qual) if($qual);
			return ($seq) if(!$qual);
		}
		else {
			return -1;
		}
	}
	else {
		return -1;
	}
}
