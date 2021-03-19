#! /usr/bin/perl

use File::Basename;
#BEGIN {
#	my ($tmp, $path) = fileparse($0);
#	push ( @INC,"$path/lib");
#	#use lib "$path";
#}
use strict;
use warnings;
use Getopt::Long;
use List::Util qw(sum min max);
use Cwd qw(abs_path);
use IO::Zlib;
use FindBin qw($RealBin);
use lib "$RealBin/lib";
require "html.pl";

eval {
	require Parallel::ForkManager;
	require String::Approx;
	require GD::Graph::linespoints;
	require GD::Graph::bars;
	require GD::Graph::pie;
	require GD::Text::Wrap;
};

my $isGDMod = 1;

if($@) {
	my $errorText = join("", $@);
	if($errorText =~ /Parallel/) {
		print "Error:\n\tCan not find 'lib' folder with this perl program\n"; #module 'Parallel::ForkManager'\n";
		print "\tCopy the 'lib' folder, provided with the toolkit, to the directory where this perl program is and try again\n\n";
		exit;
	}
	elsif($errorText =~ /GD\/Graph\/linespoints/) {
		print STDERR "Warning:\n\tCan not find module 'GD::Graph'\n";
		print STDERR "\tGraphs for statistics will not be produced. \n\t\t\tOR \n\tInstall GD::Graph module and try again.\n\n";
		$isGDMod = 0;
	}
	elsif($errorText =~ /String\/Approx/) {
		print STDERR "Error:\n\tCan not find module 'String::Approx'\n";
		print STDERR "\tInstall it and try again\n\n";
		exit;
	}
}


# Stat variables.
my @totalBases = (0, 0);
my @totalHQBases = (0, 0);	### This is for HQ reads, does not include NO_VEC stat.
my @totalBasesAfterHQ = (0, 0);
my @totalHQBasesAfterHQ = (0, 0);
my @totalBasesFinal = (0, 0);
my @totalHQBasesFinal = (0, 0);
my @totalReadsAfterHQ = (0, 0);	### This is for HQ reads, does not include NO_VEC stat.
my @totalReads = (0, 0);
my @totalValidReadsNoPriAda = (0, 0);
my @totalValidReadsWithPriAda = (0, 0);
my @minLen = (1000, 1000, 1000, 1000);
my @maxLen = (0, 0, 0, 0);
my @positionSpecificBaseCount = ();
my @positionSpecificBaseCountHQ = ();	#### This is for final output reads only... Which may include NO_VEC according to the user input.
my @positionSpecificBaseCountWithRanges = ();
my @positionSpecificBaseCountHQWithRanges = ();
my @totalReadsFinal = ();
my @fileName = ();
my @outFileName = ();
my @readsWithN = (0, 0, 0, 0);
my @totalNs = (0, 0, 0, 0);
my @totalTrimmedReads = (0, 0);
my @priAdaLibNames = ("Genomic DNA/Chip-seq Library", "Paired End DNA Library", "DpnII gene expression Library", "NlaIII gene expression Library", "Small RNA Library", "Multiplexing DNA Library");
my $nLines = 0;
my $seqFormat = 0;		# 1: Sanger; 2: Solexa; 3: Illumina 1.3+; 4: Illumina 1.5+;
my $subVal = 0;			# 33: Sanger; 64: Illumina
my @qualDistribRaw = ();
my @qualDistribFinal = ();
my $qualDistribInterval = 1;
my @qualLabel = ();
my @gcDistribRaw = ();
my @gcDistribFinal = ();
my $gcDistribInterval = 5;
my @gcLabel = ();
my @baseCountRaw = ();
my @baseCountFinal = ();
my @charCountRaw = ();
my @charCountFinal = ();

#my @monoRepeat = ();				### Poly A, Poly T, Poly G, Poly C
#my @diRepeat = ();
#my @triRepeat = ();
#my @tetraRepeat = ();

my $font_spec = getFilePath($0) . "lib/Fonts/Dustismo_Sans.ttf";
my $f = getFilePath($0) . "lib/Fonts/LucidaSansDemiBold.ttf";


# Misc variables.
my $isPairedEnd = 0;
my $outFolder = "";
my $substrlen = 20;						# For removePriAda
my $mismLim = 1;						# For removePriAda
my $indOfAnalysis = 0;

# Parameter variables.
my @peFiles = ();
my @seFiles = ();
my @allFiles = ();
my $noOfInp4PE = 4;
my $noOfInp4SE = 3;
my $isOnlyStat;
my $priAdaLib = "";
my $cutOffReadLen4HQ = 70;
my $cutOffPhScore = 20;
#my $trimAfterUnknownCall;
my $noOfProcesses = 1;
my $helpAsked;
my $statOutFmt = 1;			# 1: Text format; 2: Tab-delimited format.
my $priAdaFile;
my @usrDefinedPriAda = ();
my $outputDataFmt = "t";		# t/T: Text; g/G: Gzip.
GetOptions(
			"pe=s{$noOfInp4PE}" => \@peFiles,
			"se=s{$noOfInp4SE}" => \@seFiles,
			"h|help" => \$helpAsked,
			"l|cutOffReadLen4HQ=f" => \$cutOffReadLen4HQ,
			"o|outputFolder=s" => \$outFolder,
			"z|outputDataCompression=s" => \$outputDataFmt,
			"t|statOutFmt=i" => \$statOutFmt,
			"onlyStat" => \$isOnlyStat,
			"p|processes=i" => \$noOfProcesses,
			"s|cutOffQualScore=i" => \$cutOffPhScore,
		  );
if($helpAsked) {
	prtUsage();
	exit;
}
if(@peFiles == 0 && @seFiles == 0) {
	prtError("No input files are provided");
}
# Validating inputs
my @tempFiles = ();
prtError("Missing inputs for paired-end files") if((scalar @peFiles)%$noOfInp4PE != 0);
for(my $i=0; $i<@peFiles; $i+=$noOfInp4PE) {
	my $str = "$peFiles[$i] $peFiles[$i+1] $peFiles[$i+2] $peFiles[$i+3]";
	if($peFiles[$i+2] =~ /^-/) {
		prtError("Missing inputs for paired-end files: at '-pe $str'")
	}
	if($peFiles[$i+2] =~ /^\d$/) {
		if($peFiles[$i+2] < 1 || $peFiles[$i+2] > 6) {
			prtError("Incorrect option for Primer/Adaptor library: at '-pe $str'");
		}
	}
	if($peFiles[$i+3] =~ /^-/) {
		prtError("Missing inputs for paired-end files: at '-pe $str'")
	}
	if($peFiles[$i+3] !~ /\d/ && $peFiles[$i+3] !~ /a/i) {
		prtError("Incorrect option for FASTQ variant: at '-pe $str'")
	}
	if($peFiles[$i+3] !~ /a/i) {
		if($peFiles[$i+3] < 1 || $peFiles[$i+3] > 5) {
			prtError("Incorrect option for FASTQ variant: at '-pe $str'");
		}
	}
	push(@tempFiles, $str);
}
@peFiles = ();
@peFiles = @tempFiles;
@tempFiles = ();
prtError("Missing inputs for single-end files") if((scalar @seFiles)%$noOfInp4SE != 0);
for(my $i=0; $i<@seFiles; $i+=$noOfInp4SE) {
	my $str = "$seFiles[$i] $seFiles[$i+1] $seFiles[$i+2]";
	if($seFiles[$i+1] =~ /^-/) {
		prtError("Missing inputs for single-end files: at '-se $str'")
	}
	if($seFiles[$i+1] =~ /^\d$/i) {
		if($seFiles[$i+1] < 1 || $seFiles[$i+1] > 6) {
			prtError("Incorrect option for Primer/Adaptor library: at '-se $str'");
		}
	}
	if($seFiles[$i+2] =~ /^-/) {
		prtError("Missing inputs for single-end files: at '-se $str'")
	}
	if($seFiles[$i+2] !~ /\d/ && $seFiles[$i+2] !~ /a/i) {
		prtError("Incorrect option for FASTQ variant: at '-se $str'")
	}
	if($seFiles[$i+2] !~ /a/i) {
		if($seFiles[$i+2] < 1 || $seFiles[$i+2] > 5) {
			prtError("Incorrect option for FASTQ variant: at '-se $str'");
		}
	}
	push(@tempFiles, $str);
}
@seFiles = ();
@seFiles = @tempFiles;
@tempFiles = ();
if($cutOffReadLen4HQ < 0 || $cutOffReadLen4HQ > 100) {
	prtError("Incorrect value for -l|cutOffReadLen4HQ option: at '-l $cutOffReadLen4HQ'");
}
if($statOutFmt < 1 || $statOutFmt > 2) {
	prtError("Incorrect value for -statOutFmt: at '-statOutFmt $statOutFmt'");
}
if($outputDataFmt !~ /^[tg]$/i) {
	prtError("Incorrect value for -f|outputDataFmt option: at '-f $outputDataFmt'");
}

my $pm = new Parallel::ForkManager($noOfProcesses);

@allFiles = (@peFiles, @seFiles);
my $pid;

foreach my $file (@allFiles) {
	@totalBases = (0, 0);
	@totalHQBases = (0, 0);
	@totalBasesAfterHQ = (0, 0);
	@totalHQBasesAfterHQ = (0, 0);
	@totalBasesFinal = (0, 0);
	@totalHQBasesFinal = (0, 0);
	@totalReadsAfterHQ = (0, 0);
	@totalReads = (0, 0);
	@totalValidReadsNoPriAda = (0, 0);
	@totalValidReadsWithPriAda = (0, 0);
	@minLen = (1000, 1000, 1000, 1000);
	@maxLen = (0, 0, 0, 0);
	@positionSpecificBaseCount = ();
	@positionSpecificBaseCountHQ = ();
	@positionSpecificBaseCountWithRanges = ();
	@positionSpecificBaseCountHQWithRanges = ();
	@totalReadsFinal = ();
	@fileName = ();
 	@outFileName = ();
	@readsWithN = (0, 0, 0, 0);
	@totalNs = (0, 0, 0, 0);
	@totalTrimmedReads = (0, 0);
	@qualDistribRaw = ();
	@qualDistribFinal = ();
	@qualLabel = ();
	@gcDistribRaw = ();
	@gcDistribFinal = ();
	@gcLabel = ();
	@charCountRaw = ();
	@charCountFinal = ();
	$priAdaFile = "";
	@usrDefinedPriAda = ();



	$file =~ s/\\([A-Za-z_\.])/\/$1/g;		# To remove '\' from the path of windows file
	$isPairedEnd = 0;
	my @inpData = split(/\s+/, $file);
	if($inpData[$#inpData-1] =~ /^n$/i) {
		undef $priAdaLib;
	}
	elsif($inpData[$#inpData-1] =~ /^\d$/) {
		$priAdaLib = $inpData[$#inpData-1] - 1;
	}
	else {
		$priAdaLib = "u";
		$priAdaFile = $inpData[$#inpData-1];
		open(PRIADA, "<$priAdaFile") or die "Can not open the user-defined primer/adapter file: $priAdaFile\n";
		@usrDefinedPriAda = <PRIADA>;
		for(my $i=0; $i<=$#usrDefinedPriAda; $i++) {
			$usrDefinedPriAda[$i] =~ s/\s+//g;
		}
	}
	$seqFormat = $inpData[$#inpData];
	
	$indOfAnalysis++;
$pid = $pm->start and next;
	print "Analysis has been started for \"$file\": Index: $indOfAnalysis\n";
	my $statFile = "";
	my $outFile1;
	my $outFile2;
	my $unPaired;
	my $outFile;
	if((scalar @inpData) == $noOfInp4PE) {
		@fileName = ($inpData[0], $inpData[1]);
		$outFolder = getFilePath($fileName[0]) . "IlluQC_Filtered_files" if($outFolder eq "");
		$outFolder .= "/" if($outFolder !~ /\/$/);
		if(! -e $outFolder) {
				mkdir($outFolder) or die "Can not create output folder: $outFolder\n";
		}
		$outFile1 = $outFolder . getFileName($fileName[0]) . "_filtered";
		$outFile2 = $outFolder . getFileName($fileName[1]) . "_filtered";
		$outFile1 .= ".gz" if($outputDataFmt =~ /g/i);
		$outFile2 .= ".gz" if($outputDataFmt =~ /g/i);
		$statFile = $outFolder . getFileName($fileName[0]) . "_" . getFileName($fileName[1]) . "_stat";
		$outFileName[0] = $outFile1;
		$outFileName[1] = $outFile2;
		if($seqFormat =~ /a/i) {
			print "$indOfAnalysis: Checking FASTQ format: File $fileName[0]...\n";
			$nLines = checkFastQFormat($fileName[0], 1);
			print "$indOfAnalysis: Checking FASTQ format: File $fileName[1]...\n";
			if($nLines != checkFastQFormat($fileName[1], 1)) {
				prtErrorExit("Number of reads in paired end files are not same.\n\t\tFiles: $fileName[0], $fileName[1]");
			}
			if($seqFormat == 1) {
				$subVal = 33;
				print "$indOfAnalysis: Input FASTQ file format: Sanger\n";
			}
			if($seqFormat == 2) {
				$subVal = 64;
				print "$indOfAnalysis: Input FASTQ file format: Solexa\n";
			}
			if($seqFormat == 3) {
				$subVal = 64;
				print "$indOfAnalysis: Input FASTQ file format: Illumina 1.3+\n";
			}
			if($seqFormat == 4) {
				$subVal = 64;
				print "$indOfAnalysis: Input FASTQ file format: Illumina 1.5+\n";
			}
			if($seqFormat == 5) {
				$subVal = 33;
				print "$indOfAnalysis: Input FASTQ file format: Illumina 1.8+\n";
			}
		}
		else {
			$nLines = checkFastQFormat($fileName[0], 0);
			if($nLines != checkFastQFormat($fileName[1], 0)) {
				prtErrorExit("Number of reads in paired end files are not same.\n\t\tFiles: $fileName[0], $fileName[1]");
			}
			if($seqFormat == 1 || $seqFormat == 5) {
				$subVal = 33;
			}
			else {
				$subVal = 64;
			}
		}
		print "$indOfAnalysis: Processing input files...\n";
		$unPaired = getFilePath($outFile1) . getFileName($fileName[0]) . "_" . getFileName($fileName[1]) . "_unPaired_HQReads";
		$unPaired .= ".gz" if($outputDataFmt =~ /g/i);
		my $t = sprintf("%0.0f", $nLines/4);
		print "$indOfAnalysis: Number of reads processed: " . "0/$t (0\%)...\n";
		processPairedEndFiles($fileName[0], $fileName[1], $outFile1, $outFile2, $unPaired);
		print "$indOfAnalysis: Number of reads processed: " . "$totalReads[0]/$totalReads[0] (100\%)...\n";
		if(!defined($isOnlyStat)) {
		}
		$isPairedEnd = 1;
	}
	else {
		$fileName[0] = $inpData[0]; #$arg;
		$outFolder = getFilePath($fileName[0]) . "IlluQC_Filtered_files" if($outFolder eq "");
		$outFolder .= "/" if($outFolder !~ /\/$/);
		if(! -e $outFolder) {
			if(!defined($isOnlyStat)) {
				mkdir($outFolder) or die "Can not create output folder: $outFolder\n";
			}
		}
		$outFile = $outFolder . getFileName($fileName[0]) . "_filtered";
		$outFile .= ".gz" if($outputDataFmt =~ /g/i);
		$outFileName[0] = $outFile;
		$statFile = $outFolder . getFileName($fileName[0]) . "_stat";
		if($seqFormat =~ /a/i) {
			print "$indOfAnalysis: Checking FASTQ format: File $fileName[0]...\n";
			$nLines = checkFastQFormat($fileName[0], 1);
			if($seqFormat == 1) {
				$subVal = 33;
				print "$indOfAnalysis: Input FASTQ file format: Sanger\n";
			}
			if($seqFormat == 2) {
				$subVal = 64;
				print "$indOfAnalysis: Input FASTQ file format: Solexa\n";
			}
			if($seqFormat == 3) {
				$subVal = 64;
				print "$indOfAnalysis: Input FASTQ file format: Illumina 1.3+\n";
			}
			if($seqFormat == 4) {
				$subVal = 64;
				print "$indOfAnalysis: Input FASTQ file format: Illumina 1.5+\n";
			}
			if($seqFormat == 5) {
				$subVal = 33;
				print "$indOfAnalysis: Input FASTQ file format: Illumina 1.8+\n";
			}
		}
		else {
			$nLines = checkFastQFormat($fileName[0], 0);
			if($seqFormat == 1 || $seqFormat == 5) {
				$subVal = 33;
			}
			else {
				$subVal = 64;
			}
		}
		print "$indOfAnalysis: Processing input files...\n";
		my $t = sprintf("%0.0f", $nLines/4);
		print "$indOfAnalysis: Number of reads processed: " . "0/$t (0\%)...\n";
		processSingleEndFiles($fileName[0], $outFile);
		print "$indOfAnalysis: Number of reads processed: " . "$totalReads[0]/$totalReads[0] (100\%)...\n";
		if(!defined($isOnlyStat)) {
		}
		$isPairedEnd = 0;
	}
	
	print "$indOfAnalysis: Analysis completed\n";
	
	print "$indOfAnalysis: Printing Statistics...\n";

	my $qualDistF1 = getFileName($fileName[0])."_qualDistribution.png";
	my $qualDistF2 = getFileName($fileName[1])."_qualDistribution.png" if($isPairedEnd);
	my $sumPieF;
	$sumPieF = getFileName($fileName[0]). "_summary.png";
	$sumPieF = getFileName($fileName[0]). "_" . getFileName($fileName[1]) ."_summary.png" if($isPairedEnd);
	my $gcDistF1 = getFileName($fileName[0])."_gcDistribution.png";
	my $gcDistF2 = getFileName($fileName[1])."_gcDistribution.png" if($isPairedEnd);
	my $baseCntF1 = getFileName($fileName[0])."_baseCompostion.png";
	my $baseCntF2 = getFileName($fileName[1])."_baseCompostion.png" if($isPairedEnd);
	my $avgQF1 = getFileName($fileName[0]) . "_avgQual.png";
	my $avgQF2 = getFileName($fileName[1]) . "_avgQual.png" if($isPairedEnd);
	my $QRangeRawF1 = getFileName($fileName[0]) . "_QualRangePerBase.png";
	my $QRangeFilteredF1 = getFileName($outFileName[0]) . "_QualRangePerBase.png" if(!defined($isOnlyStat));
	my $QRangeF1 = "$QRangeRawF1";
	$QRangeF1 .= ":::$QRangeFilteredF1" if(!defined($isOnlyStat));
	my $QRangeF2;
	if($isPairedEnd) {
		my $QRangeRawF2 = getFileName($fileName[1]) . "_QualRangePerBase.png";
		my $QRangeFilteredF2 = getFileName($outFileName[1]) . "_QualRangePerBase.png" if(!defined($isOnlyStat));
		$QRangeF2 = "$QRangeRawF2";
		$QRangeF2 .= ":::$QRangeFilteredF2" if(!defined($isOnlyStat));
	}

	my $c=0;
	foreach my $ref (@qualDistribRaw) {
		my $str = "";
		foreach my $val (@{$ref}) {
			if($c == 0) {
				$str = "0";
				$str .= "-$qualDistribInterval" if($qualDistribInterval>1);
			}
			else {
				$str = $qualDistribInterval*$c;
				$str .=  "-" . $qualDistribInterval*($c+1) if($qualDistribInterval>1);
			}
			$c++;
			push(@qualLabel, $str);
		}
		last;
	}
	my @file1 = (\@qualLabel, $qualDistribRaw[0]);
	my @file2 = (\@qualLabel, $qualDistribRaw[1]) if($isPairedEnd);
	if(!$isOnlyStat) {
		push(@file1, $qualDistribFinal[0]);
		push(@file2, $qualDistribFinal[1]) if($isPairedEnd);;
	}
	if($isGDMod) {
		drawQualDist(\@file1, $outFolder.$qualDistF1, getFileName($fileName[0]), 650, 350);
		drawQualDist(\@file2, $outFolder.$qualDistF2, getFileName($fileName[1]), 650, 300) if($isPairedEnd);
	}
	
	my $readsWPriAda = $totalReadsAfterHQ[0] - $totalReadsFinal[0];		# For Paired end, different number of contaminated sequences will be filtered in both the files. And we have to report total reads contaminated including both end files.
	my $readsLowQual = $totalReads[0] - $readsWPriAda - $totalReadsFinal[0];

	@file1 = (["", "", ""], [$readsWPriAda, $totalReadsFinal[0], $readsLowQual]);
	if($isGDMod) {
		drawSummaryPie(\@file1, $outFolder.$sumPieF, 500, 350);
	}

	$c=0;
	foreach my $ref (@gcDistribRaw) {
		foreach my $val (@{$ref}) {
			my $str = "";
			if($c == 0) {
				$str = "0-$gcDistribInterval";
			}
			else {
				$str = $gcDistribInterval*$c . "-" . $gcDistribInterval*($c+1);
			}
			$c++;
			push(@gcLabel, $str);
		}
		last;
	}

	@file1 = (\@gcLabel, $gcDistribRaw[0]);
	@file2 = (\@gcLabel, $gcDistribRaw[1]) if($isPairedEnd);
	if(!$isOnlyStat) {
		push(@file1, $gcDistribFinal[0]);
		push(@file2, $gcDistribFinal[1]) if($isPairedEnd);
	}

	if($isGDMod) {
		drawGCDist(\@file1, $outFolder.$gcDistF1, getFileName($fileName[0]), 550, 350);
		drawGCDist(\@file2, $outFolder.$gcDistF2, getFileName($fileName[1]), 550, 350) if($isPairedEnd);
	}

	@file1 = (["A", "T", "G", "C", "Non-ATGC"], $charCountRaw[0]);
	@file2 = (["A", "T", "G", "C", "Non-ATGC"], $charCountRaw[1]) if($isPairedEnd);
	if(!$isOnlyStat) {
		@file1 = (["A", "T", "G", "C", "Non-ATGC"], $charCountRaw[0], $charCountFinal[0]);
		@file2 = (["A", "T", "G", "C", "Non-ATGC"], $charCountRaw[1], $charCountFinal[1]) if($isPairedEnd);
	}
	if($isGDMod) {
		drawBaseComp(\@file1, $outFolder.$baseCntF1, getFileName($fileName[0]), 500, 300);
		drawBaseComp(\@file2, $outFolder.$baseCntF2, getFileName($fileName[1]), 500, 300) if($isPairedEnd);
	}

	open(STAT, ">$statFile") or die "Can not create statistics file $statFile\n";
	printStat(*STAT) if($statOutFmt == 1);
	printStatTab(*STAT) if($statOutFmt == 2);
	close(STAT);
	
	my $iFol = getFilePath(abs_path($fileName[0]));
	my $oFol = abs_path($outFolder) . "/";
	my $inpFs = getFileName($fileName[0]);
	my $seqFormatName;
	$inpFs .= ":::::" . getFileName($fileName[1]) if($isPairedEnd);
	my $htF = $oFol . "output_" . getFileName($fileName[0]);
	$htF .= "_" . getFileName($fileName[1]) if($isPairedEnd);
	$htF .= ".html";
	if($seqFormat == 1) {
		$seqFormatName = "Sanger";
	}
	elsif($seqFormat == 2) {
		$seqFormatName = "Solexa";
	}
	elsif($seqFormat == 3) {
		$seqFormatName = "Illumina 1.3+";
	}
	elsif($seqFormat == 4) {
		$seqFormatName = "Illumina 1.5+";
	}
	my @fileNames4HTML;
	@fileNames4HTML = ($outFile, $avgQF1, $baseCntF1, $gcDistF1, $qualDistF1, $sumPieF, $QRangeF1);
	@fileNames4HTML = ($outFile1, $outFile2, $unPaired, $avgQF1, $avgQF2, $baseCntF1, $baseCntF2, $gcDistF1, $gcDistF2, $qualDistF1, $qualDistF2, $sumPieF, $QRangeF1, $QRangeF2) if($isPairedEnd);
	htmlPrint(getFilePath(abs_path($0)), getFileName($0), $htF, $iFol, $isPairedEnd, $isOnlyStat, $inpFs, $seqFormatName, $statFile, $oFol, \@fileNames4HTML);
	$pm->finish;
}
$pm->wait_all_children;

print "================================================================\n";
print "Processing has been finished\n";
print "Output files are generated in $outFolder\n" if($outFolder ne "");
print "Output files are generated in the folder of input files\n" if($outFolder eq "");
print "================================================================\n";

exit;

sub openFileGetHandle {
	my ($file, $rOrw) = @_;
	my $fh;
	if($file =~ /\.gz$/i) {
		$fh = new IO::Zlib;
		$fh->open("$file", "rb") or die "Can not open file $file" if($rOrw eq "r");
		$fh->open("$file", "wb") or die "Can not create file $file" if($rOrw eq "w");
	}
	else {
		open($fh, "<$file") or die "Can not open file $file" if($rOrw eq "r");
		open($fh, ">$file") or die "Can not create file $file" if($rOrw eq "w");
	}
	return $fh;
}

sub processPairedEndFiles {
	my $file1 = $_[0];
	my $file2 = $_[1];
	my $outFile1 = $_[2];
	my $outFile2 = $_[3];
	my $unPaired = $_[4];
	$totalReads[0] = sprintf("%0.0f", $nLines/4);
	$totalReads[1] = sprintf("%0.0f", $nLines/4);
	
	my $fH1 = openFileGetHandle($file1, "r");
	*F1 = $fH1;
	my $fH2 = openFileGetHandle($file2, "r");
	*F2 = $fH2;

	if(!defined($isOnlyStat)) {
		my $ofH1 = openFileGetHandle($outFile1, "w");
		*OF1 = $ofH1;
		my $ofH2 = openFileGetHandle($outFile2, "w");
		*OF2 = $ofH2;
		my $ofupH = openFileGetHandle($unPaired, "w");
		*OFUP = $ofupH;
	}
	
	my $isEOF = 1;
	if($nLines/4 > 0) {
		$isEOF = 0;
	}
	
	my $lineCount = 0;
	
	while(!$isEOF) {
		my @fRead = ();
		my @rRead = ();
		for(my $i=0; $i<4; $i++) {
			$fRead[$i] = <F1>;
			$rRead[$i] = <F2>;
		}
		last if($fRead[0]=~ /^\n$/);
		last if($rRead[0]=~ /^\n$/);
		chomp(my $fQualLine = $fRead[3]);
		chomp(my $rQualLine = $rRead[3]);
		chomp(my $fSeqLine = $fRead[1]);
		chomp(my $rSeqLine = $rRead[1]);
		my $fNs = getNoOfNs($fSeqLine);
		my $rNs = getNoOfNs($rSeqLine);
		$totalNs[0] += $fNs;
		$totalNs[1] += $rNs;
		if($fNs) {
			$readsWithN[0]++;
		}
		if($rNs) {
			$readsWithN[1]++;
		}
		
		my @qualArr = ();
		my $isFReadOfHQ = isReadOfHQ($fQualLine, 0, \@qualArr);
		my $isRReadOfHQ = isReadOfHQ($rQualLine, 1, \@qualArr);
		my $fSeqLineLen = length $fSeqLine;
		my $rSeqLineLen = length $rSeqLine;
		my @ASCII = unpack("C*", $fQualLine);
		my $fAvgQual = sprintf "%.0f", (sum(@ASCII)/$fSeqLineLen);
		my @rASCII = unpack("C*", $rQualLine);
		my $rAvgQual = sprintf "%.0f", (sum(@rASCII)/$rSeqLineLen);
		$fAvgQual -= $subVal;
		$rAvgQual -= $subVal;
		$qualDistribRaw[0][getIndex($fAvgQual,$qualDistribInterval)]++;
		$qualDistribRaw[1][getIndex($rAvgQual,$qualDistribInterval)]++;
		my $fAs = $fSeqLine =~ s/A/A/gi;
		my $fTs = $fSeqLine =~ s/T/T/gi;
		my $fGs = $fSeqLine =~ s/G/G/gi;
		my $fCs = $fSeqLine =~ s/C/C/gi;
		my $fgcPercent = ($fGs + $fCs)/$fSeqLineLen*100; 
		$charCountRaw[0][0] += $fAs;
		$charCountRaw[0][1] += $fTs;
		$charCountRaw[0][2] += $fGs;
		$charCountRaw[0][3] += $fCs;
		$charCountRaw[0][4] += $fNs;
		my $rAs = $rSeqLine =~ s/A/A/gi;
		my $rTs = $rSeqLine =~ s/T/T/gi;
		my $rGs = $rSeqLine =~ s/G/G/gi;
		my $rCs = $rSeqLine =~ s/C/C/gi;
		my $rgcPercent = ($rGs + $rCs)/$rSeqLineLen*100;
		$charCountRaw[1][0] += $rAs;
		$charCountRaw[1][1] += $rTs;
		$charCountRaw[1][2] += $rGs;
		$charCountRaw[1][3] += $rCs;
		$charCountRaw[1][4] += $rNs; 
		$gcDistribRaw[0][getIndex($fgcPercent,$gcDistribInterval)]++;
		$gcDistribRaw[1][getIndex($rgcPercent,$gcDistribInterval)]++;
		if($isFReadOfHQ && $isRReadOfHQ) {
			$totalReadsAfterHQ[0]++;
			$totalReadsAfterHQ[1]++;
			$totalBasesAfterHQ[0] += $fSeqLineLen;
			$totalBasesAfterHQ[1] += $rSeqLineLen;
			$totalHQBasesAfterHQ[0] += $isFReadOfHQ;
			$totalHQBasesAfterHQ[1] += $isRReadOfHQ;
			if(defined($priAdaLib)) {
				my $isFWOPriAda = isWOPriAda($fSeqLine, 0, 1);
				my $isRWOPriAda = isWOPriAda($rSeqLine, 1, 1);
				if($isFWOPriAda && $isRWOPriAda) {
					$totalReadsFinal[0]++;
					$totalReadsFinal[1]++;
					$totalBasesFinal[0] += $fSeqLineLen;
					$totalBasesFinal[1] += $rSeqLineLen;
					$totalHQBasesFinal[0] += $isFReadOfHQ;
					$totalHQBasesFinal[1] += $isRReadOfHQ;
					$minLen[2] = $fSeqLineLen if($minLen[2] > $fSeqLineLen);
					$maxLen[2] = $fSeqLineLen if($maxLen[2] < $fSeqLineLen);
					$minLen[3] = $rSeqLineLen if($minLen[3] > $rSeqLineLen);
					$maxLen[3] = $rSeqLineLen if($maxLen[3] < $rSeqLineLen);
					$totalNs[2] += $fNs;
					$totalNs[3] += $rNs;
					if($fNs) {
						$readsWithN[2]++;
					}
					if($rNs) {
						$readsWithN[3]++;
					}
					for(my $x=0; $x<@qualArr; $x++) {
						my @row = @{$qualArr[$x]};
						for(my $y=0; $y<@row; $y++) {
							$positionSpecificBaseCountHQ[$x][$y] += $qualArr[$x][$y];
							my $ind = int($qualArr[$x][$y]/10);
							$ind-- if($qualArr[$x][$y]%10 == 0 && $qualArr[$x][$y] != 0);
							$positionSpecificBaseCountHQWithRanges[$x][$y][$ind]++;
						}
					}
					if(!defined($isOnlyStat)) {
						$qualDistribFinal[0][getIndex($fAvgQual,$qualDistribInterval)]++;
						$qualDistribFinal[1][getIndex($rAvgQual,$qualDistribInterval)]++;
						$gcDistribFinal[0][getIndex($fgcPercent,$gcDistribInterval)]++;
						$gcDistribFinal[1][getIndex($rgcPercent,$gcDistribInterval)]++;
						$charCountFinal[0][0] += $fAs;
						$charCountFinal[0][1] += $fTs;
						$charCountFinal[0][2] += $fGs;
						$charCountFinal[0][3] += $fCs;
						$charCountFinal[0][4] += $fNs;
						$charCountFinal[1][0] += $rAs;
						$charCountFinal[1][1] += $rTs;
						$charCountFinal[1][2] += $rGs;
						$charCountFinal[1][3] += $rCs;
						$charCountFinal[1][4] += $rNs; 
						print OF1 @fRead;
						print OF2 @rRead;
					}
				}
				else {
					if(!defined($isOnlyStat)) {
						if($isFWOPriAda) {
							print OFUP @fRead;
						}
						elsif($isRWOPriAda) {
							print OFUP @rRead;
						}
					}
				}
			}
			else {
				$totalReadsFinal[0]++;
				$totalReadsFinal[1]++;
				$totalBasesFinal[0] += $fSeqLineLen;
				$totalBasesFinal[1] += $rSeqLineLen;
				$totalHQBasesFinal[0] += $isFReadOfHQ;
				$totalHQBasesFinal[1] += $isRReadOfHQ;
				$minLen[2] = $fSeqLineLen if($minLen[2] > $fSeqLineLen);
				$maxLen[2] = $fSeqLineLen if($maxLen[2] < $fSeqLineLen);
				$minLen[3] = $rSeqLineLen if($minLen[3] > $rSeqLineLen);
				$maxLen[3] = $rSeqLineLen if($maxLen[3] < $rSeqLineLen);
				$totalNs[2] += $fNs;
				$totalNs[3] += $rNs;
				if($fNs) {
					$readsWithN[2]++;
				}
				if($rNs) {
					$readsWithN[3]++;
				}
				for(my $x=0; $x<@qualArr; $x++) {
					my @row = @{$qualArr[$x]};
					for(my $y=0; $y<@row; $y++) {
						$positionSpecificBaseCountHQ[$x][$y] += $qualArr[$x][$y];
						my $ind = int($qualArr[$x][$y]/10);
						$ind-- if($qualArr[$x][$y]%10 == 0 && $qualArr[$x][$y] != 0);
						$positionSpecificBaseCountHQWithRanges[$x][$y][$ind]++;
					}
				}				
				if(!defined($isOnlyStat)) {
					$qualDistribFinal[0][getIndex($fAvgQual,$qualDistribInterval)]++;
					$qualDistribFinal[1][getIndex($rAvgQual,$qualDistribInterval)]++;
					$gcDistribFinal[0][getIndex($fgcPercent,$gcDistribInterval)]++;
					$gcDistribFinal[1][getIndex($rgcPercent,$gcDistribInterval)]++;
					$charCountFinal[0][0] += $fAs;
					$charCountFinal[0][1] += $fTs;
					$charCountFinal[0][2] += $fGs;
					$charCountFinal[0][3] += $fCs;
					$charCountFinal[0][4] += $fNs;
					$charCountFinal[1][0] += $rAs;
					$charCountFinal[1][1] += $rTs;
					$charCountFinal[1][2] += $rGs;
					$charCountFinal[1][3] += $rCs;
					$charCountFinal[1][4] += $rNs; 
					print OF1 @fRead;
					print OF2 @rRead;
				}
			}
		}
		else {
			if(!defined($isOnlyStat)) {
				if($isFReadOfHQ) {
					my $isFWOPriAda = 1;
					$isFWOPriAda = isWOPriAda($fSeqLine, 0, 0) if(defined($priAdaLib));
					if($isFWOPriAda) {
						print OFUP @fRead;
					}
				}
				elsif($isRReadOfHQ) {
					my $isRWOPriAda = 1;
					$isRWOPriAda = isWOPriAda($rSeqLine, 1, 0) if(defined($priAdaLib));
					if($isRWOPriAda) {
						print OFUP @rRead;
					}
				}
			}
		}
		$lineCount += 4;
		if($lineCount >= $nLines) {
			$isEOF = 1;
		}
		if($lineCount % (100000*4) == 0) {
			my $tmpP = sprintf "%0.0f", ($lineCount/4/$totalReads[0]*100);
			print "$indOfAnalysis: Number of reads processed: " . $lineCount/4 . "/$totalReads[0] ($tmpP\%)...\n";
		}
	}
	close(OFUP);
	close(OF2);
	close(OF1);
	close(F2);
	close(F1);
}

sub processSingleEndFiles {
	my $file = $_[0];
	my $outFile = $_[1];
	$totalReads[0] = sprintf("%0.0f", $nLines/4);
	
	my $fH = openFileGetHandle($file, "r");
	*F = $fH;

	if(!defined($isOnlyStat)) {
		my $ofH = openFileGetHandle($outFile, "w");
		*OF = $ofH;
	}
	
	my $isEOF = 1;
	if($nLines/4 > 0) {
		$isEOF = 0;
	}
	
	my $lineCount = 0;
	
	while(!$isEOF) {
		my @fRead = ();
		for(my $i=0; $i<4; $i++) {
			$fRead[$i] = <F>;
		}
		last if($fRead[0]=~ /^\n$/);
		chomp(my $fQualLine = $fRead[3]);
		chomp(my $fSeqLine = $fRead[1]);
		my $fNs = getNoOfNs($fSeqLine);
		$totalNs[0] += $fNs;
		if($fNs) {
			$readsWithN[0]++;
		}

		my @qualArr = ();
		my $isFReadOfHQ = isReadOfHQ($fQualLine, 0, \@qualArr);
		my $fSeqLineLen = length $fSeqLine;
		my @ASCII = unpack("C*", $fQualLine);
		my $fAvgQual = sprintf "%.0f", (sum(@ASCII)/$fSeqLineLen);
		$fAvgQual -= $subVal;
		$qualDistribRaw[0][getIndex($fAvgQual,$qualDistribInterval)]++;
		my $fAs = $fSeqLine =~ s/A/A/gi;
		my $fTs = $fSeqLine =~ s/T/T/gi;
		my $fGs = $fSeqLine =~ s/G/G/gi;
		my $fCs = $fSeqLine =~ s/C/C/gi;
		my $fgcPercent = ($fGs + $fCs)/$fSeqLineLen*100; 
		$charCountRaw[0][0] += $fAs;
		$charCountRaw[0][1] += $fTs;
		$charCountRaw[0][2] += $fGs;
		$charCountRaw[0][3] += $fCs;
		$charCountRaw[0][4] += $fNs;
		$gcDistribRaw[0][getIndex($fgcPercent,$gcDistribInterval)]++;
		if($isFReadOfHQ) {
			$totalReadsAfterHQ[0]++;
			$totalBasesAfterHQ[0] += $fSeqLineLen;
			$totalHQBasesAfterHQ[0] += $isFReadOfHQ;
			if(defined($priAdaLib)) {
				if(isWOPriAda($fSeqLine, 0, 1)) {
					$totalReadsFinal[0]++;
					$totalBasesFinal[0] += $fSeqLineLen;
					$totalHQBasesFinal[0] += $isFReadOfHQ;
					$minLen[1] = $fSeqLineLen if($minLen[1] > $fSeqLineLen);
					$maxLen[1] = $fSeqLineLen if($maxLen[1] < $fSeqLineLen);
					$totalNs[1] += $fNs;
					if($fNs) {
						$readsWithN[1]++;
					}
					for(my $x=0; $x<@qualArr; $x++) {
						my @row = @{$qualArr[$x]};
						for(my $y=0; $y<@row; $y++) {
							$positionSpecificBaseCountHQ[$x][$y] += $qualArr[$x][$y];
							my $ind = int($qualArr[$x][$y]/10);
							$ind-- if($qualArr[$x][$y]%10 == 0 && $qualArr[$x][$y] != 0);
							$positionSpecificBaseCountHQWithRanges[$x][$y][$ind]++;
						}
					}				
					if(!defined($isOnlyStat)) {
						$qualDistribFinal[0][getIndex($fAvgQual,$qualDistribInterval)]++;
						$gcDistribFinal[0][getIndex($fgcPercent,$gcDistribInterval)]++;
						$charCountFinal[0][0] += $fAs;
						$charCountFinal[0][1] += $fTs;
						$charCountFinal[0][2] += $fGs;
						$charCountFinal[0][3] += $fCs;
						$charCountFinal[0][4] += $fNs;
						print OF @fRead;
					}
				}
			}
			else {
				$totalReadsFinal[0]++;
				$totalBasesFinal[0] += $fSeqLineLen;
				$totalHQBasesFinal[0] += $isFReadOfHQ;
				$minLen[1] = $fSeqLineLen if($minLen[1] > $fSeqLineLen);
				$maxLen[1] = $fSeqLineLen if($maxLen[1] < $fSeqLineLen);
				$totalNs[1] += $fNs;
				if($fNs) {
					$readsWithN[1]++;
				}
				for(my $x=0; $x<@qualArr; $x++) {
					my @row = @{$qualArr[$x]};
					for(my $y=0; $y<@row; $y++) {
						$positionSpecificBaseCountHQ[$x][$y] += $qualArr[$x][$y];
						my $ind = int($qualArr[$x][$y]/10);
						$ind-- if($qualArr[$x][$y]%10 == 0 && $qualArr[$x][$y] != 0);
						$positionSpecificBaseCountHQWithRanges[$x][$y][$ind]++;
					}
				}
				if(!defined($isOnlyStat)) {
					$qualDistribFinal[0][getIndex($fAvgQual,$qualDistribInterval)]++;
					$gcDistribFinal[0][getIndex($fgcPercent,$gcDistribInterval)]++;
					$charCountFinal[0][0] += $fAs;
					$charCountFinal[0][1] += $fTs;
					$charCountFinal[0][2] += $fGs;
					$charCountFinal[0][3] += $fCs;
					$charCountFinal[0][4] += $fNs;
					print OF @fRead;
				}
			}
		}
		else {
		}
		$lineCount += 4;
		if($lineCount >= $nLines) {
			$isEOF = 1;
		}
		if($lineCount % (100000*4) == 0) {
			my $tmpP = sprintf "%0.0f", ($lineCount/4/$totalReads[0]*100);
			print "$indOfAnalysis: Number of reads processed: " . $lineCount/4 . "/$totalReads[0] ($tmpP\%)...\n";
		}
	}
	
	close(OF);
	close(F);
}

sub checkFastQFormat {				# Takes FASTQ file as an input and if the format is incorrect it will print error and exit, otherwise it will return the number of lines in the file.
	my $file = $_[0];
	my $isVariantIdntfcntOn = $_[1];
	my $lines = 0;
	my $fH = openFileGetHandle($file, "r");
	*F = $fH;
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

sub getFilePath {
	my $name = $_[0];
	my $path = "";
	if($name =~ /\//) {
		$name =~ /(.+)\//;
		$path = $1 . "/";
	}
	else {
		$path = "./";
	}
	return $path;
}

sub getFileName {	# This sub takes a path of a file and returns just its name after separating the path from it.
	my $path = $_[0];
	my $name = "";
	$path =~ /([^\/]+)$/;
	$name = $1;
	return $name;	
}

sub prtErrorExit {
	my $errmsg = $_[0];
	print STDERR "Error:\t", $errmsg, "\n";
	exit;
}

sub isReadOfHQ {	# Criteria for HQ is greater than or equal to 70% of bases have phred score >= 20
	my $read = $_[0];
	my $v0Or1 = $_[1];		# 0 will be for forward reads and 1 for reverse reads.
	my $arrRef = $_[2];
	my $readLen = length $read;
	$minLen[$v0Or1] = $readLen if($minLen[$v0Or1] > $readLen);
	$maxLen[$v0Or1] = $readLen if($maxLen[$v0Or1] < $readLen);
	my $cutOffLen = sprintf("%0.0f", $readLen * $cutOffReadLen4HQ / 100);	# 70% length of read length is calculated.
	my $validBaseCount = 0;
	my @ASCII = unpack("C*", $read);
	my $c = 0;
	foreach my $val (@ASCII) {
		$val -= $subVal;
		$positionSpecificBaseCount[$v0Or1][$c] += $val;
		my $ind = int($val/10);
		$ind-- if($val%10 == 0 && $val != 0);
		$positionSpecificBaseCountWithRanges[$v0Or1][$c][$ind]++;
		$$arrRef[$v0Or1][$c] = $val;
		if($val >= $cutOffPhScore) {
			$validBaseCount++;
		}
		$c++;
	}
	$totalBases[$v0Or1] += $readLen;
	$totalHQBases[$v0Or1] += $validBaseCount;
	if($validBaseCount >= $cutOffLen) {
		return $validBaseCount;				# Return true.
	}
	else {
		return 0;				# Return false.
	}
}


sub qualGraph {				### Use this just for final graph generation...
	my $file = $_[0];
	my $v0Or1 = $_[1];		# 0 will be for forward reads and 1 for reverse reads.
	my $arrRef = $_[2];
	my $flag = 0;
	my $fH = openFileGetHandle($file, "r");
	*F = $fH;
	while(my $read = <F>) {
		chomp($read);
		$flag++;
		if($flag%4 == 0) {			# To obtain the quality value line.
			my @ASCII = ();
			@ASCII = unpack("C*", $read);
			my $c=0;
			foreach my $val (@ASCII) {
				$val -= $subVal;
				$$arrRef[$v0Or1][$c] += $val;
				$c++;
			}
		}
	}
	close(F);
}







sub isWOPriAda {
	my $seq = $_[0];
	my $v0Or1 = $_[1];
	my $isCountStatOn = $_[2];
	chomp($seq);

	my @arrGenomic = (
		"GATCGGAAGAGCTCGTATGCCGTCTTCTGCTTG",
		"ACACTCTTTCCCTACACGACGCTCTTCCGATCT",
		"AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT",
		"CAAGCAGAAGACGGCATACGAGCTCTTCCGATCT",
		"ACACTCTTTCCCTACACGACGCTCTTCCGATCT"
	);

	my @arrPE = (
		"GATCGGAAGAGCGGTTCAGCAGGAATGCCGAG",
		"ACACTCTTTCCCTACACGACGCTCTTCCGATCT",
		"AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT",
		"CAAGCAGAAGACGGCATACGAGATCGGTCTCGGCATTCCTGCTGAACCGCTCTTCCGATCT",
		"ACACTCTTTCCCTACACGACGCTCTTCCGATCT",
		"CGGTCTCGGCATTCCTGCTGAACCGCTCTTCCGATCT"
	);

	my @arrDpnII = (
		"GATCGTCGGACTGTAGAACTCTGAAC",
		"ACAGGTTCAGAGTTCTACAGTCCGAC",
		"CAAGCAGAAGACGGCATACGANN",
		"TCGTATGCCGTCTTCTGCTTG",
		"CAAGCAGAAGACGGCATACGA",
		"AATGATACGGCGACCACCGACAGGTTCAGAGTTCTACAGTCCGA",
		"CGACAGGTTCAGAGTTCTACAGTCCGACGATC"
	);

	my @arrNlaIII = (
		"TCGGACTGTAGAACTCTGAAC",
		"ACAGGTTCAGAGTTCTACAGTCCGACATG",
		"CAAGCAGAAGACGGCATACGANN",
		"TCGTATGCCGTCTTCTGCTTG",
		"CAAGCAGAAGACGGCATACGA",
		"AATGATACGGCGACCACCGACAGGTTCAGAGTTCTACAGTCCGA",
		"CCGACAGGTTCAGAGTTCTACAGTCCGACATG"
	);

	my @arrsmRNA = (
		"GTTCAGAGTTCTACAGTCCGACGATC",
		"TCGTATGCCGTCTTCTGCTTGT",
		"CAAGCAGAAGACGGCATACGA",
		"CAAGCAGAAGACGGCATACGA",
		"AATGATACGGCGACCACCGACAGGTTCAGAGTTCTACAGTCCGA",
		"CGACAGGTTCAGAGTTCTACAGTCCGACGATC"
	);

	my @arrmulPlex = (
		"GATCGGAAGAGCACACGTCT",
		"ACACTCTTTCCCTACACGACGCTCTTCCGATCT",
		"AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT",
		"GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT",
		"ACACTCTTTCCCTACACGACGCTCTTCCGATCT",
		"GATCGGAAGAGCACACGTCTGAACTCCAGTCAC",
		"GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT",
		"CAAGCAGAAGACGGCATACGAGATCGTGATGTGACTGGAGTTC",
		"CAAGCAGAAGACGGCATACGAGATACATCGGTGACTGGAGTTC",
		"CAAGCAGAAGACGGCATACGAGATGCCTAAGTGACTGGAGTTC",
		"CAAGCAGAAGACGGCATACGAGATTGGTCAGTGACTGGAGTTC",
		"CAAGCAGAAGACGGCATACGAGATCACTGTGTGACTGGAGTTC",
		"CAAGCAGAAGACGGCATACGAGATATTGGCGTGACTGGAGTTC",
		"CAAGCAGAAGACGGCATACGAGATGATCTGGTGACTGGAGTTC",
		"CAAGCAGAAGACGGCATACGAGATTCAAGTGTGACTGGAGTTC",
		"CAAGCAGAAGACGGCATACGAGATCTGATCGTGACTGGAGTTC",
		"CAAGCAGAAGACGGCATACGAGATAAGCTAGTGACTGGAGTTC",
		"CAAGCAGAAGACGGCATACGAGATGTAGCCGTGACTGGAGTTC",
		"CAAGCAGAAGACGGCATACGAGATTACAAGGTGACTGGAGTTC"
	);
	
	my @priAdas = (\@arrGenomic, \@arrPE, \@arrDpnII, \@arrNlaIII, \@arrsmRNA, \@arrmulPlex);
	my %checkedPriStr = ();	# The 20 bp from start and end are stored in this hash as key. So that next time when another pri/ada seq

	my @priAdaSeqs = ();
	if($priAdaLib eq "u") {
		@priAdaSeqs = @usrDefinedPriAda;
	}
	else {
		@priAdaSeqs = @{$priAdas[$priAdaLib]};
	}
	my @stat = ();
	my $priInd = 0;
	
	my $isMatched = 0;
	foreach my $priAda (@priAdaSeqs) {
		if(findSeq($priAda, $seq, \%checkedPriStr)) {
			$isMatched = 1;
			last;
		}
	}
	
	if($isMatched) {
		$totalValidReadsWithPriAda[$v0Or1]++ if($isCountStatOn);
		return 0;
	}
	else {
		$totalValidReadsNoPriAda[$v0Or1]++ if($isCountStatOn);
		return 1;
	}
}

sub findSeq {
	my $pri = $_[0];
	my $seq = $_[1];
	my $hashRef = $_[2];
	my $spri = substr($pri, 0, $substrlen);
	my $tmpInd = (length $pri) - $substrlen;
	$tmpInd = 0 if($tmpInd < 0);
	my $epri = substr($pri, $tmpInd, $substrlen);
	my $ans;
	if(!defined($$hashRef{$spri})) {
		my @catches = String::Approx::amatch($spri, ['I0 D0 S1'], $seq);
		if(@catches != 0) {
			return 1;
		}
		$$hashRef{$spri} = 1;
	}
	if(!defined($$hashRef{$epri})) {
		my @catches = String::Approx::amatch($epri, ['I0 D0 S1'], $seq);
		if(@catches != 0) {
			return 1;
		}
		$$hashRef{$epri} = 1;
	}
	return 0;
}


sub prtHelp {
	print "\n$0 options:\n\n";
	print "### Input reads (FASTQ) options (Atleast one option is required)\n";
	print "  -pe <Forward reads file> <Reverse reads file> <Primer/Adaptor library> <FASTQ variant>\n";
	print "    Paired-end read files (FASTQ) with primer/adaptor library and FASTQ variant\n";
	print "    User may choose from the provided primer/adaptor library or can give a file containing primer/adaptor sequences, one per line\n";
	print "    Multiple libraries can be given using multiple '-pe' options\n";
	print "      For eg.: -pe r1.fq r2.fq 3 1 -pe t1.fq t2.fq 2 A\n\n";
	print "  -se <Reads file> <Primer/Adaptor library> <FASTQ variant>\n";
	print "    Single-end read file (FASTQ) with primer/adaptor library and FASTQ variant\n";
	print "    Multiple libraries can be given using multiple '-se' options\n";
	print "      For eg.: -se r1.fq 3 2 -se t2.fq 2 2\n\n";
	print "    Primer/Adaptor libraries:\n";
	my $c = 1;
	foreach my $lib (@priAdaLibNames) {
		print "      $c = $lib\n";
		$c++;
	}
	print "      N = Do not filter for Primer/Adaptor\n";
	print "      <File> = File for user defined primer/adaptor sequences, one per line\n";
	print "\n";
	print "    FASTQ variants:\n";
	print "      1 = Sanger (Phred+33, 33 to 73)\n";
	print "      2 = Solexa (Phred+64, 59 to 104)\n";
	print "      3 = Illumina (1.3+) (Phred+64, 64 to 104)\n";
	print "      4 = Illumina (1.5+) (Phred+64, 66 to 104)\n";
	print "      5 = Illumina (1.8+) (Phred+33, 33 to 74)\n";
	print "      A = Automatic detection of FASTQ variant\n";
	print "\n";
	print "### Other options [Optional]\n";
	print "  -h | -help\n";
	print "    Prints this help\n";
	print "--------------------------------- QC Options ---------------------------------\n";
	print "  -l | -cutOffReadLen4HQ <Real number, 0 to 100>\n";
	print "    The cut-off value for percentage of read length that should be of given quality\n";
	print "    default: 70\n";
	print "  -s | -cutOffQualScore <Integer, 0 to 40>\n";
	print "    The cut-off value for PHRED quality score for high-quality filtering\n";
	print "    default: 20\n";
	print "----------------------------- Processing Options -----------------------------\n";
	print "  -p | -processes <Integer>\n";
	print "    Number of processes to be used\n";
	print "    default: 1\n";
	print "  -onlyStat\n";
	print "    Outputs only statistics without filtered data output\n";
	print "------------------------------- Output Options -------------------------------\n";
	print "  -t | -statOutFmt <Integer>\n";
	print "    Output format for statistics\n";
	print "    Formats:\n";
	print "      1 = formatted text\n";
	print "      2 = tab delimited\n";
	print "    default: 1\n";
	print "  -o | -outputFolder <Output folder name/path>\n";
	print "    Output will be stored in the given folder\n";
	print "    default: By default, output folder (IlluQC_Filtered_files) will be generated where the input files are\n";
	print "  -z | -outputDataCompression <Character>\n";
	print "    Output format for HQ filtered data\n";
	print "    Formats:\n";
	print "      t = text FASTQ files\n";
	print "      g = gzip compressed files\n";
	print "    default: t\n";
	print "\n";
}

sub prtUsage {
	print "\nUsage: perl $0 <options>\n";
	prtHelp();
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


sub getNoOfNs {			# This takes sequence and returns the number of N/. (unknown base call).
	my $seq = $_[0];
	my $count = 0;
	while($seq =~ /[N\.]/g) {
		$count++;
	}
	return $count;
}

sub getIndex {
	my $up = $_[0];
	my $down = $_[1];
	my $inp = $up/$down;
	return (sprintf "%0.0f", $up) if($down == 1);
	my $index = int((sprintf "%0.2f", $inp)+0.99)-1;
	$index = 0 if($index < 0);
	return $index;
}


sub drawBaseComp {
	my $dataRef = $_[0];
	my $fileNameWPath = $_[1];
	my $fileName = $_[2];
	my $width = $_[3];
	my $height = $_[4];
	
	my $mygraph = GD::Graph::bars->new($width, $height);
	
	$mygraph->set( 
		y_label => 'Count',
		y_min_value => 0,
		box_axis => 0,
		line_width => 3,
		transparent => 0,
		dclrs => [ qw(lred dgreen) ],
		legend_placement	=> 'BR',
		x_label_position	=> 1/2,
		long_ticks			=> 1,
		fgclr               => '#dddddd',
		l_margin 			=> 60,
		r_margin 			=> 60,
		b_margin			=> 50,
		t_margin			=> 50,
		show_values         => 1,
		bar_spacing         => 1,
		values_vertical 	=> 1,
	) or warn $mygraph->error;
	
	if(!defined($isOnlyStat)) {
		$mygraph->set_legend( $fileName, $fileName."_filtered");
	}
	else {
		$mygraph->set_legend( $fileName);
	}

    $mygraph->set_y_label_font($font_spec, 12);
    $mygraph->set_x_label_font($font_spec, 12);
    $mygraph->set_y_axis_font($font_spec, 10);
    $mygraph->set_x_axis_font($font_spec, 8);
    $mygraph->set_title_font($f, 11);
    $mygraph->set_legend_font($f, 8);
    $mygraph->set_values_font($f, 6);

	my $myImage = $mygraph->plot($dataRef);

	my $black = $myImage->colorAllocate(0,0,0);			# To set the color for the next time printing on the image.
	my $lred = $myImage->colorAllocate(255,0,0);	
	my $dgreen = $myImage->colorAllocate(0,127,0);
	my $dblue = $myImage->colorAllocate(0,0,127);
	
	my $sum1 = sum(@{$$dataRef[1]});
	my $sum2 = sum(@{$$dataRef[2]});
	
	my $wrapbox = GD::Text::Wrap->new( $myImage,
	        line_space  => 4,
	        text        => "Base composition for $fileName",
	        color		=> $dblue,
	);
	
	$wrapbox->set(align => 'center', width => $width);
	$wrapbox->set_font($f, 11);
	$wrapbox->draw(0,0);

	$wrapbox = GD::Text::Wrap->new( $myImage,
	        line_space  => 4,
	        text        => "A (" . (sprintf "%0.2f", @{$$dataRef[1]}[0]/$sum1*100) . "\%)",
	        color		=> $black,
	);
	
	$wrapbox->set(align => 'left', width => 300);
	$wrapbox->set_font($f, 8);
	$wrapbox->draw($width/2-220,$height-35);

	$wrapbox = GD::Text::Wrap->new( $myImage,
	        line_space  => 4,
	        text        => "T (" . (sprintf "%0.2f", @{$$dataRef[1]}[1]/$sum1*100) . "\%)",
	        color		=> $black,
	);
	
	$wrapbox->set(align => 'left', width => 300);
	$wrapbox->set_font($f, 8);
	$wrapbox->draw($width/2-130,$height-35);

	$wrapbox = GD::Text::Wrap->new( $myImage,
	        line_space  => 4,
	        text        => "G (" . (sprintf "%0.2f", @{$$dataRef[1]}[2]/$sum1*100) . "\%)",
	        color		=> $black,
	);
	
	$wrapbox->set(align => 'left', width => 300);
	$wrapbox->set_font($f, 8);
	$wrapbox->draw($width/2-40,$height-35);

	$wrapbox = GD::Text::Wrap->new( $myImage,
	        line_space  => 4,
	        text        => "C (" . (sprintf "%0.2f", @{$$dataRef[1]}[3]/$sum1*100) . "\%)",
	        color		=> $black,
	);
	
	$wrapbox->set(align => 'left', width => 300);
	$wrapbox->set_font($f, 8);
	$wrapbox->draw($width/2+50,$height-35);

	$wrapbox = GD::Text::Wrap->new( $myImage,
	        line_space  => 4,
	        text        => "Non-ATGC (" . (sprintf "%0.2f", @{$$dataRef[1]}[4]/$sum1*100) . "\%)",
	        color		=> $black,
	);
	
	$wrapbox->set(align => 'left', width => 300);
	$wrapbox->set_font($f, 8);
	$wrapbox->draw($width/2+140,$height-35);


	my $startRectX = $width/2-230;
	my $startRectY = $height-35;
	$myImage->filledRectangle($startRectX,$startRectY,$startRectX+8,$startRectY+8,$lred);
	$myImage->rectangle($startRectX,$startRectY,$startRectX+8,$startRectY+8,$black);

	$startRectX = $width/2-140;
	$myImage->filledRectangle($startRectX,$startRectY,$startRectX+8,$startRectY+8,$lred);
	$myImage->rectangle($startRectX,$startRectY,$startRectX+8,$startRectY+8,$black);

	$startRectX = $width/2-50;
	$myImage->filledRectangle($startRectX,$startRectY,$startRectX+8,$startRectY+8,$lred);
	$myImage->rectangle($startRectX,$startRectY,$startRectX+8,$startRectY+8,$black);

	$startRectX = $width/2+40;
	$myImage->filledRectangle($startRectX,$startRectY,$startRectX+8,$startRectY+8,$lred);
	$myImage->rectangle($startRectX,$startRectY,$startRectX+8,$startRectY+8,$black);

	$startRectX = $width/2+130;
	$myImage->filledRectangle($startRectX,$startRectY,$startRectX+8,$startRectY+8,$lred);
	$myImage->rectangle($startRectX,$startRectY,$startRectX+8,$startRectY+8,$black);


	if(!$isOnlyStat) {
		$wrapbox = GD::Text::Wrap->new( $myImage,
		        line_space  => 4,
		        text        => "A (" . (sprintf "%0.2f", @{$$dataRef[2]}[0]/$sum2*100) . "\%)",
		        color		=> $black,
		);
		
		$wrapbox->set(align => 'left', width => 300);
		$wrapbox->set_font($f, 8);
		$wrapbox->draw($width/2-220,$height-20);
	
		$wrapbox = GD::Text::Wrap->new( $myImage,
		        line_space  => 4,
		        text        => "T (" . (sprintf "%0.2f", @{$$dataRef[2]}[1]/$sum2*100) . "\%)",
		        color		=> $black,
		);
		
		$wrapbox->set(align => 'left', width => 300);
		$wrapbox->set_font($f, 8);
		$wrapbox->draw($width/2-130,$height-20);
	
		$wrapbox = GD::Text::Wrap->new( $myImage,
		        line_space  => 4,
		        text        => "G (" . (sprintf "%0.2f", @{$$dataRef[2]}[2]/$sum2*100) . "\%)",
		        color		=> $black,
		);
		
		$wrapbox->set(align => 'left', width => 300);
		$wrapbox->set_font($f, 8);
		$wrapbox->draw($width/2-40,$height-20);
	
		$wrapbox = GD::Text::Wrap->new( $myImage,
		        line_space  => 4,
		        text        => "C (" . (sprintf "%0.2f", @{$$dataRef[2]}[3]/$sum2*100) . "\%)",
		        color		=> $black,
		);
		
		$wrapbox->set(align => 'left', width => 300);
		$wrapbox->set_font($f, 8);
		$wrapbox->draw($width/2+50,$height-20);
	
		$wrapbox = GD::Text::Wrap->new( $myImage,
		        line_space  => 4,
		        text        => "Non-ATGC (" . (sprintf "%0.2f", @{$$dataRef[2]}[4]/$sum2*100) . "\%)",
		        color		=> $black,
		);
		
		$wrapbox->set(align => 'left', width => 300);
		$wrapbox->set_font($f, 8);
		$wrapbox->draw($width/2+140,$height-20);
	

	
		$startRectX = $width/2-230;
		$startRectY = $height-20;
		$myImage->filledRectangle($startRectX,$startRectY,$startRectX+8,$startRectY+8,$dgreen);
		$myImage->rectangle($startRectX,$startRectY,$startRectX+8,$startRectY+8,$black);
	
		$startRectX = $width/2-140;
		$myImage->filledRectangle($startRectX,$startRectY,$startRectX+8,$startRectY+8,$dgreen);
		$myImage->rectangle($startRectX,$startRectY,$startRectX+8,$startRectY+8,$black);
	
		$startRectX = $width/2-50;
		$myImage->filledRectangle($startRectX,$startRectY,$startRectX+8,$startRectY+8,$dgreen);
		$myImage->rectangle($startRectX,$startRectY,$startRectX+8,$startRectY+8,$black);
	
		$startRectX = $width/2+40;
		$myImage->filledRectangle($startRectX,$startRectY,$startRectX+8,$startRectY+8,$dgreen);
		$myImage->rectangle($startRectX,$startRectY,$startRectX+8,$startRectY+8,$black);
	
		$startRectX = $width/2+130;
		$myImage->filledRectangle($startRectX,$startRectY,$startRectX+8,$startRectY+8,$dgreen);
		$myImage->rectangle($startRectX,$startRectY,$startRectX+8,$startRectY+8,$black);
	}





	open(I, ">$fileNameWPath") or print STDERR "Error:\n\tCan not create image file: $fileNameWPath\n";
	binmode I;
	print I $myImage->png;
	close(I);
}


sub drawGCDist {
	my $dataRef = $_[0];
	my $fileNameWPath = $_[1];
	my $fileName = $_[2];
	my $width = $_[3];
	my $height = $_[4];
	
	my $mygraph = GD::Graph::linespoints->new($width, $height);
	
	$mygraph->set( 
		x_label => '% GC content',
		y_label => 'Number of reads',
		title => "GC content distribution for $fileName",
		y_min_value => 0,
		box_axis => 0,
		line_width => 3,
		transparent => 0,
		markers				=> [1],
		marker_size			=> 3,
		dclrs => [ qw(lred dgreen) ],
		legend_placement	=> 'BR',
		x_label_position	=> 1/2,
		x_labels_vertical	=> 1,
		long_ticks			=> 1,
		fgclr               => '#dddddd',
	) or warn $mygraph->error;
	
	if(!defined($isOnlyStat)) {
		$mygraph->set_legend( $fileName, $fileName."_filtered");
	}
	else {
		$mygraph->set_legend( $fileName, "a");
	}

    $mygraph->set_y_label_font($font_spec, 12);
    $mygraph->set_x_label_font($font_spec, 12);
    $mygraph->set_y_axis_font($font_spec, 10);
    $mygraph->set_x_axis_font($font_spec, 8);
    $mygraph->set_title_font($f, 11);
    $mygraph->set_legend_font($f, 8);

	my $myImage = $mygraph->plot($dataRef);
	open(I, ">$fileNameWPath") or print STDERR "Error:\n\tCan not create image file: $fileNameWPath\n";
	binmode I;
	print I $myImage->png;
	close(I);
}


sub drawSummaryPie {
	my $dataRef = $_[0];
	my $fileNameWPath = $_[1];
	my $width = $_[2];
	my $height = $_[3];
	my $mygraph = new GD::Graph::pie($width, $height);
	
	$mygraph->set( 
		title => "Summary of quality check and filtering",
		axislabelclr => 'black',
		pie_height => 40,
	
		l_margin => 15,
		r_margin => 15,
		b_margin => 50,
		start_angle => 45,
		dclrs => [ qw(lyellow lgreen lred) ],
		transparent => 0,
	) or warn $mygraph->error;
	
    $mygraph->set_label_font($f, 8);
    $mygraph->set_value_font(['verdana', 'arial'],14);
    $mygraph->set_title_font($f, 11);

	my $myImage = $mygraph->plot($dataRef);
	
	my $black = $myImage->colorAllocate(0,0,0);			# To set the color for the next time printing on the image.
	my $red = $myImage->colorAllocate(255,0,0);	
	my $yellow = $myImage->colorAllocate(255,255,0);	
	my $green = $myImage->colorAllocate(0,255,0);
	
	my $sum = sum(@{$$dataRef[1]});
	
	my $wrapbox = GD::Text::Wrap->new( $myImage,
	        line_space  => 4,
	        text        => (sprintf "Primer/Adaptor contaminated reads (%0.02f", @{$$dataRef[1]}[0]/$sum*100) . "\%)",
	        color		=> $black,
	);
	
	$wrapbox->set(align => 'left', width => 300);
	$wrapbox->set_font($f, 8);
	$wrapbox->draw($width/2-100,$height-45);

	$wrapbox = GD::Text::Wrap->new( $myImage,
	        line_space  => 4,
	        text        => (sprintf "High quality filtered reads (%0.02f", @{$$dataRef[1]}[1]/$sum*100) . "\%)",
	        color		=> $black,
	);
	
	$wrapbox->set(align => 'left', width => 300);
	$wrapbox->set_font($f, 8);
	$wrapbox->draw($width/2-100,$height-30);

	$wrapbox = GD::Text::Wrap->new( $myImage,
	        line_space  => 4,
	        text        => (sprintf "Low quality reads (%0.02f", @{$$dataRef[1]}[2]/$sum*100) . "\%)",
	        color		=> $black,
	);
	
	$wrapbox->set(align => 'left', width => 300);
	$wrapbox->set_font($f, 8);
	$wrapbox->draw($width/2-100,$height-15);

	my $startRectX = $width/2-120;
	my $startRectY = $height-45;
	$myImage->filledRectangle($startRectX,$startRectY,$startRectX+8,$startRectY+8,$yellow);
	$myImage->rectangle($startRectX,$startRectY,$startRectX+8,$startRectY+8,$black);

	$startRectY += 15;
	$myImage->filledRectangle($startRectX,$startRectY,$startRectX+8,$startRectY+8,$green);
	$myImage->rectangle($startRectX,$startRectY,$startRectX+8,$startRectY+8,$black);

	$startRectY += 15;
	$myImage->filledRectangle($startRectX,$startRectY,$startRectX+8,$startRectY+8,$red);
	$myImage->rectangle($startRectX,$startRectY,$startRectX+8,$startRectY+8,$black);

	open(I, ">$fileNameWPath") or print STDERR "Error:\n\tCan not create image file: $fileNameWPath\n";
	binmode I;
	print I $myImage->png;
	close(I);
	
}

sub drawQualDist {
	my $dataRef = $_[0];
	my $fileNameWPath = $_[1];
	my $fileName = $_[2];
	my $width = $_[3];
	my $height = $_[4];
	
	my $mygraph = GD::Graph::bars->new($width, $height);
	
	$mygraph->set( 
		x_label => 'Average phred quality score',
		y_label => 'Number of reads',
		title => "Quality distribution for $fileName",
		y_min_value => 0,
		box_axis => 0,
		line_width => 3,
		transparent => 0,
		dclrs => [ qw(lred dgreen) ],
		legend_placement	=> 'BR',
		x_label_position	=> 1/2,
		long_ticks			=> 1,
		fgclr               => '#dddddd',
		bar_spacing         => 1,
	) or warn $mygraph->error;
	
	if(!defined($isOnlyStat)) {
		$mygraph->set_legend( $fileName, $fileName."_filtered");
	}
	else {
		$mygraph->set_legend( $fileName);
	}

    $mygraph->set_y_label_font($font_spec, 12);
    $mygraph->set_x_label_font($font_spec, 12);
    $mygraph->set_y_axis_font($font_spec, 10);
    $mygraph->set_x_axis_font($font_spec, 8);
    $mygraph->set_title_font($f, 11);
    $mygraph->set_legend_font($f, 8);

	my $myImage = $mygraph->plot($dataRef);
	open(I, ">$fileNameWPath") or print STDERR "Error:\n\tCan not create image file: $fileNameWPath\n";
	binmode I;
	print I $myImage->png;
	close(I);
}


sub drawGraph {
	my @data = @{$_[0]};
	my $fileNameWPath = $_[1];
	my $fileName = $_[2];
	open(I, ">$fileNameWPath") or print STDERR "Error:\n\tCan not create image file: $fileNameWPath\n";
	binmode I;
	my $y_min = 0;
	my $y_max = 0;
	for(my $i=1; $i<@data; $i++) {
		$y_max = max($y_max, max(@{$data[$i]}));
	}
	$y_max = (sprintf "%0.0f",($y_max/5)) * 5 + 5;
	my $height = sprintf "%0.0f", $y_max * 300 / 45;
	my $width = sprintf "%0.0f", scalar @{$data[0]} * 600 / 75;
	my $mygraph = GD::Graph::linespoints->new($width, $height);
	$mygraph->set(
		x_label 			=> 'Base position',
		y_label				=> 'Average quality score',
		title				=> $fileName,
		y_min_value			=> $y_min,
		y_max_value			=> $y_max,
		x_label_skip		=> 2,
		y_tick_number		=> $y_max/5,
		y_label_skip		=> 1,
		markers				=> [7],
		marker_size			=> 3,
		long_ticks			=> 1,
		line_width			=> 2,
		dclrs				=> [ qw(lred dgreen) ],
		legend_placement	=> 'BR',
		x_label_position	=> 1/2,
		x_labels_vertical	=> 1,
		transparent			=> 0,
		r_margin			=> 10,
	fgclr               => '#dddddd',
	accentclr           => 'yellow',
	) or warn $mygraph->error;

	if(!defined($isOnlyStat)) {
		$mygraph->set_legend( $fileName, $fileName."_filtered");
	}
	else {
		$mygraph->set_legend( $fileName);
	}


    $mygraph->set_y_label_font($font_spec, 12);
    $mygraph->set_x_label_font($font_spec, 12);
    $mygraph->set_y_axis_font($font_spec, 10);
    $mygraph->set_x_axis_font($font_spec, 8);
    $mygraph->set_title_font($f, 11);
    $mygraph->set_legend_font($f, 8);
	my $myimage = $mygraph->plot(\@data) or die $mygraph->error;

	print I $myimage->png;
	close(I);
}

sub drawRangeGraph {
	my @data = @{$_[0]};
	my $fileNameWPath = $_[1];
	my $fileName = $_[2];
	my $height = 350;
	open(I, ">$fileNameWPath") or print STDERR "Error:\n\tCan not create image file: $fileNameWPath\n";
	binmode I;
	my $y_min = 0;
	my $y_max = 0;
	for(my $i=1; $i<@data; $i++) {
		$y_max = max($y_max, max(@{$data[$i]}));
	}
	$y_max = 100; #(sprintf "%0.0f",($y_max/5)) * 5 + 5;
	my $width = sprintf "%0.0f", scalar @{$data[0]} * 700 / 75;
	$width = max(700, $width);
	my $mygraph = GD::Graph::linespoints->new($width, $height);
	$mygraph->set(
		x_label 			=> 'Base position',
		y_label				=> 'Read count (%)',
		title				=> "Read count (%) per base for different quality score ranges for $fileName",
		y_min_value			=> $y_min,
		y_max_value			=> $y_max,
		x_label_skip		=> 2,
		y_tick_number		=> $y_max/5,
		y_label_skip		=> 1,
		markers				=> [7],
		marker_size			=> 3,
		long_ticks			=> 1,
		line_width			=> 2,
		dclrs				=> [ qw(lred dgreen lyellow blue) ],
		legend_placement	=> 'BR',
		x_label_position	=> 1/2,
		x_labels_vertical	=> 1,
		transparent			=> 0,
		r_margin			=> 10,
	fgclr               => '#dddddd',
	accentclr           => 'yellow',
	) or warn $mygraph->error;

		$mygraph->set_legend( "0-10", "11-20", "21-30", "31-40");


    $mygraph->set_y_label_font($font_spec, 12);
    $mygraph->set_x_label_font($font_spec, 12);
    $mygraph->set_y_axis_font($font_spec, 10);
    $mygraph->set_x_axis_font($font_spec, 8);
    $mygraph->set_title_font($f, 11);
    $mygraph->set_legend_font($f, 8);
	my $myimage = $mygraph->plot(\@data) or die $mygraph->error;

	print I $myimage->png;
	close(I);
}

sub prepareData4RangeGraph {
	my $STAT = $_[0];
	print $STAT "Read count (%) per base for different quality score ranges\n\n";
	my $c = 0;
	my @rangeGraphData = ();
	foreach my $arr (@positionSpecificBaseCountWithRanges) {
		my $arrFiltered = $positionSpecificBaseCountHQWithRanges[$c];
		print $STAT "\t", getFileName($fileName[$c]);
		print $STAT "\t\t\t\t\t", getFileName($outFileName[$c]) if(!defined($isOnlyStat));
		print $STAT "\n";
		print $STAT "Ranges\t0-10\t11-20\t21-30\t31-40";
		print $STAT "\t\t0-10\t11-20\t21-30\t31-40" if(!defined($isOnlyStat));
		print $STAT "\nBase\n";
		my $basePos = 1;
		foreach my $valArr (@$arr) {
			my $valArrF = @$arrFiltered[$basePos-1];
			@$valArr[0] = 0 if(! @$valArr[0]);
			@$valArr[1] = 0 if(! @$valArr[1]);
			@$valArr[2] = 0 if(! @$valArr[2]);
			@$valArr[3] = 0 if(! @$valArr[3]);
			my $total = @$valArr[0] + @$valArr[1] + @$valArr[2] + @$valArr[3];
			my $val1 = sprintf "%0.2f", @$valArr[0]/$total*100;
			my $val2 = sprintf "%0.2f", @$valArr[1]/$total*100;
			my $val3 = sprintf "%0.2f", @$valArr[2]/$total*100;
			my $val4 = sprintf "%0.2f", @$valArr[3]/$total*100;
			$rangeGraphData[$c*2][0][$basePos-1] = $basePos;
			$rangeGraphData[$c*2][1][$basePos-1] = $val1;
			$rangeGraphData[$c*2][2][$basePos-1] = $val2;
			$rangeGraphData[$c*2][3][$basePos-1] = $val3;
			$rangeGraphData[$c*2][4][$basePos-1] = $val4;
			print $STAT "$basePos\t$val1\t$val2\t$val3\t$val4";
			if(! defined($isOnlyStat)) {
				@$valArrF[0] = 0 if(! @$valArrF[0]);
				@$valArrF[1] = 0 if(! @$valArrF[1]);
				@$valArrF[2] = 0 if(! @$valArrF[2]);
				@$valArrF[3] = 0 if(! @$valArrF[3]);
				my $totalF = @$valArrF[0] + @$valArrF[1] + @$valArrF[2] + @$valArrF[3];
				my $valF1 = sprintf "%0.2f", @$valArrF[0]/$totalF*100;
				my $valF2 = sprintf "%0.2f", @$valArrF[1]/$totalF*100;
				my $valF3 = sprintf "%0.2f", @$valArrF[2]/$totalF*100;
				my $valF4 = sprintf "%0.2f", @$valArrF[3]/$totalF*100;
				$rangeGraphData[$c*2+1][0][$basePos-1] = $basePos;
				$rangeGraphData[$c*2+1][1][$basePos-1] = $valF1;
				$rangeGraphData[$c*2+1][2][$basePos-1] = $valF2;
				$rangeGraphData[$c*2+1][3][$basePos-1] = $valF3;
				$rangeGraphData[$c*2+1][4][$basePos-1] = $valF4;
				print $STAT "\t\t$valF1\t$valF2\t$valF3\t$valF4";
			}
			print $STAT "\n";
			$basePos++;
		}
		print $STAT "\n\n";
		$c++;
	}
	print $STAT "\n\n";
	if($isGDMod) {
		drawRangeGraph($rangeGraphData[0], $outFolder.getFileName($fileName[0])."_QualRangePerBase.png", getFileName($fileName[0]));
		drawRangeGraph($rangeGraphData[1], $outFolder.getFileName($outFileName[0])."_QualRangePerBase.png", getFileName($outFileName[0])) if(!defined($isOnlyStat));
		if($isPairedEnd) {
			drawRangeGraph($rangeGraphData[2], $outFolder.getFileName($fileName[1])."_QualRangePerBase.png", getFileName($fileName[1]));
			drawRangeGraph($rangeGraphData[3], $outFolder.getFileName($outFileName[1])."_QualRangePerBase.png", getFileName($outFileName[1])) if(!defined($isOnlyStat));
		}
	}
}

sub printStat {
	my $STAT = $_[0];
	my $tmpPer;
	my $inde = " " x 1;
	print $STAT "Parameters\n";
	my @graphData1 = ();
	my @graphData2 = ();
	if($isPairedEnd) {
		printf $STAT "$inde %-30s %s\n", "Library type", "Paired-end";
		printf $STAT "$inde %-30s %s  %s\n", "Input files", $fileName[0], $fileName[1];
		printf $STAT "$inde %-30s %s\n", "Primer/Adaptor library", defined($priAdaLib)?(($priAdaLib ne "u")?$priAdaLibNames[$priAdaLib]:"User defined ($priAdaFile)"):"NA";
		printf $STAT "$inde %-30s %s\n", "Cut-off read length for HQ", $cutOffReadLen4HQ."%";
		printf $STAT "$inde %-30s %s\n", "Cut-off quality score", $cutOffPhScore;
		printf $STAT "$inde %-30s %s\n", "Only statistics", defined($isOnlyStat)?"On":"Off";
		printf $STAT "$inde %-30s %s\n", "Number of processes", $noOfProcesses;

		print $STAT "\n\n";

		print $STAT "QC statistics\n";
		printf $STAT "$inde %-50s %-20s  %s\n", "File name", getFileName($fileName[0]), getFileName($fileName[1]);
		printf $STAT "$inde %-50s %-20d %d\n", "Total number of reads", $totalReads[0], $totalReads[1];
		printf $STAT "$inde %-50s %-20d %d\n", "Total number of HQ reads", $totalReadsAfterHQ[0], $totalReadsAfterHQ[1];
		$tmpPer = sprintf "%0.2f", $totalReadsAfterHQ[0]/$totalReads[0]*100;
		printf $STAT "$inde %-50s %-20s %0.2f%s\n", "Percentage of HQ reads", $tmpPer."%", $totalReadsAfterHQ[1]/$totalReads[1]*100, "%";
		printf $STAT "$inde %-50s %-20d %d\n", "Total number of bases", $totalBases[0], $totalBases[1];
		printf $STAT "$inde %-50s %-20d %d\n", "Total number of bases in HQ reads", $totalBasesAfterHQ[0], $totalBasesAfterHQ[1];
		printf $STAT "$inde %-50s %-20d %d\n", "Total number of HQ bases in HQ reads", $totalHQBasesAfterHQ[0], $totalHQBasesAfterHQ[1];
		$tmpPer = sprintf "%0.2f", $totalHQBasesAfterHQ[0]/$totalBasesAfterHQ[0]*100;
		printf $STAT "$inde %-50s %-20s %0.2f%s\n", "Percentage of HQ bases in HQ reads", $tmpPer."%" , $totalHQBasesAfterHQ[1]/$totalBasesAfterHQ[1]*100, "%";
		if(defined($priAdaLib)) {
			printf $STAT "$inde %-50s %-20d %d\n", "Number of Primer/Adaptor contaminated HQ reads", $totalValidReadsWithPriAda[0], $totalValidReadsWithPriAda[1];
		}
		else {
			printf $STAT "$inde %-50s %-20s %s\n", "Number of Primer/Adaptor contaminated HQ reads", "NA", "NA";
		}
		printf $STAT "$inde %-50s %-20d %d\n", "Total number of HQ filtered reads", $totalReadsFinal[0], $totalReadsFinal[1];
		$tmpPer = sprintf "%0.2f", $totalReadsFinal[0]/$totalReads[0]*100;
		printf $STAT "$inde %-50s %-20s %0.2f%s\n", "Percentage of HQ filtered reads", $tmpPer."%", $totalReadsFinal[1]/$totalReads[1]*100, "%";

		print $STAT "\n\n";

		print $STAT "Detailed QC statistics\n";
		my @arr = (
			["File name", getFileName($fileName[0]), getFileName($fileName[1]), getFileName($outFileName[0]), getFileName($outFileName[1])],
			["Minimum read length", $minLen[0], $minLen[1], $minLen[2], $minLen[3]],
			["Maximum read length", $maxLen[0], $maxLen[1], $maxLen[2], $maxLen[3]],
			["Average read length", (sprintf "%0.2f", $totalBases[0]/$totalReads[0]), (sprintf "%0.2f", $totalBases[1]/$totalReads[1]), (sprintf "%0.2f", $totalBasesFinal[0]/$totalReadsFinal[0]), (sprintf "%0.2f", $totalBasesFinal[1]/$totalReadsFinal[1])],
			["Total number of reads", $totalReads[0], $totalReads[1], $totalReadsFinal[0], $totalReadsFinal[1]],
			["Total number of reads with non-ATGC bases", $readsWithN[0], $readsWithN[1], $readsWithN[2], $readsWithN[3]],
			["Percentage of reads with non-ATGC bases", (sprintf "%0.2f", $readsWithN[0]/$totalReads[0]*100)."%", (sprintf "%0.2f", $readsWithN[1]/$totalReads[1]*100)."%", (sprintf "%0.2f", $readsWithN[2]/$totalReadsFinal[0]*100)."%", (sprintf "%0.2f", $readsWithN[3]/$totalReadsFinal[1]*100)."%"],
			["Total number of bases", $totalBases[0], $totalBases[1], $totalBasesFinal[0], $totalBasesFinal[1]],
			["Total number of HQ bases", $totalHQBases[0], $totalHQBases[1], $totalHQBasesFinal[0], $totalHQBasesFinal[1]],
			["Percentage of HQ bases", (sprintf "%0.2f", $totalHQBases[0]/$totalBases[0]*100)."%", (sprintf "%0.2f", $totalHQBases[1]/$totalBases[1]*100)."%", (sprintf "%0.2f", $totalHQBasesFinal[0]/$totalBasesFinal[0]*100)."%", (sprintf "%0.2f", $totalHQBasesFinal[1]/$totalBasesFinal[1]*100)."%"],
			["Total number of non-ATGC bases", $totalNs[0], $totalNs[1], $totalNs[2], $totalNs[3]],
			["Percentage of non-ATGC bases", (sprintf "%0.2f", $totalNs[0]/$totalBases[0]*100)."%", (sprintf "%0.2f", $totalNs[1]/$totalBases[1]*100)."%", (sprintf "%0.2f", $totalNs[2]/$totalBasesFinal[0]*100)."%", (sprintf "%0.2f", $totalNs[3]/$totalBasesFinal[1]*100)."%"],
		);
		for(my $i=0; $i<@arr; $i++) {
			if(!defined($isOnlyStat)) {
				printf $STAT "$inde %-45s  %-20s  %-20s  %-20s  %s\n", $arr[$i][0], $arr[$i][1], $arr[$i][2], $arr[$i][3], $arr[$i][4];
			}
			else {
				printf $STAT "$inde %-45s  %-20s  %s\n", $arr[$i][0], $arr[$i][1], $arr[$i][2];				
			}
		}
		
		print $STAT "\n\n";
		
		print $STAT "Average quality score at each base position of input reads\n\n";
		my $c = 0;
		@graphData1 = ();
		@graphData2 = ();
		foreach my $arr (@positionSpecificBaseCount) {
			print $STAT getFileName($fileName[$c]), "\n";
			my $basePos = 1;
			foreach my $val (@$arr) {
				my $outVal = sprintf "%0.2f", $val/$totalReads[$c];
				if($c == 0) {
					$graphData1[0][$basePos-1] = $basePos;
					$graphData1[1][$basePos-1] = $outVal;
				}
				else {
					$graphData2[0][$basePos-1] = $basePos;
					$graphData2[1][$basePos-1] = $outVal;
				}
				print $STAT $outVal, "\t";
				$basePos++;
			}
			print $STAT "\n\n";
			$c++;
		}
		print $STAT "\n\n";
		if(!defined($isOnlyStat)) {
			print $STAT "Average quality score at each base position of filtered reads\n\n";
			$c = 0;
			foreach my $arr (@positionSpecificBaseCountHQ) {
				print $STAT getFileName($outFileName[$c]), "\n";
				my $basePos = 1;
				foreach my $val (@$arr) {
					my $outVal = sprintf "%0.2f", $val/$totalReadsFinal[$c];
					if($c == 0) {
						$graphData1[2][$basePos-1] = $outVal;
					}
					else {
						$graphData2[2][$basePos-1] = $outVal;
					}
					print $STAT $outVal, "\t";
					$basePos++;
				}
				$c++;
				print $STAT "\n\n";
			}
		}
		print $STAT "\n\n";
		prepareData4RangeGraph($STAT);
		if($isGDMod) {
			drawGraph(\@graphData1, $outFolder.getFileName($fileName[0])."_avgQual.png", getFileName($fileName[0]));
			drawGraph(\@graphData2, $outFolder.getFileName($fileName[1])."_avgQual.png", getFileName($fileName[1]));
		}
	}
	else {
		printf $STAT "$inde %-30s %s\n", "Library type", "Single-end";
		printf $STAT "$inde %-30s %s\n", "Input file", $fileName[0];
		printf $STAT "$inde %-30s %s\n", "Primer/Adaptor library", defined($priAdaLib)?(($priAdaLib ne "u")?$priAdaLibNames[$priAdaLib]:"User defined ($priAdaFile)"):"NA";
		printf $STAT "$inde %-30s %s\n", "Cut-off read length for HQ", $cutOffReadLen4HQ."%";
		printf $STAT "$inde %-30s %s\n", "Cut-off quality score", $cutOffPhScore;
		printf $STAT "$inde %-30s %s\n", "Only statistics", defined($isOnlyStat)?"On":"Off";
		printf $STAT "$inde %-30s %s\n", "Number of processes", $noOfProcesses;

		print $STAT "\n\n";

		print $STAT "QC statistics\n";
		printf $STAT "$inde %-50s %s\n", "File name", getFileName($fileName[0]);
		printf $STAT "$inde %-50s %d\n", "Total number of reads", $totalReads[0];
		printf $STAT "$inde %-50s %d\n", "Total number of HQ reads", $totalReadsAfterHQ[0];
		$tmpPer = sprintf "%0.2f", $totalReadsAfterHQ[0]/$totalReads[0]*100;
		printf $STAT "$inde %-50s %s\n", "Percentage of HQ reads", $tmpPer."%";
		printf $STAT "$inde %-50s %d\n", "Total number of bases", $totalBases[0];
		printf $STAT "$inde %-50s %d\n", "Total number of bases in HQ reads", $totalBasesAfterHQ[0];
		printf $STAT "$inde %-50s %d\n", "Total number of HQ bases in HQ reads", $totalHQBasesAfterHQ[0];
		$tmpPer = sprintf "%0.2f", $totalHQBasesAfterHQ[0]/$totalBasesAfterHQ[0]*100;
		printf $STAT "$inde %-50s %s\n", "Percentage of HQ bases in HQ reads", $tmpPer."%";
		if(defined($priAdaLib)) {
			printf $STAT "$inde %-50s %d\n", "Number of Primer/Adaptor contaminated HQ reads", $totalValidReadsWithPriAda[0];
		}
		else {
			printf $STAT "$inde %-50s %s\n", "Number of Primer/Adaptor contaminated HQ reads", "NA";
		}
		printf $STAT "$inde %-50s %d\n", "Total number of HQ filtered reads", $totalReadsFinal[0];
		$tmpPer = sprintf "%0.2f", $totalReadsFinal[0]/$totalReads[0]*100;
		printf $STAT "$inde %-50s %s\n", "Percentage of HQ filtered reads", $tmpPer."%";

		print $STAT "\n\n";

		print $STAT "Detailed QC statistics\n";
		my @arr = (
			["File name", getFileName($fileName[0]), getFileName($outFileName[0])],
			["Minimum read length", $minLen[0], $minLen[1]],
			["Maximum read length", $maxLen[0], $maxLen[1]],
			["Average read length", (sprintf "%0.2f", $totalBases[0]/$totalReads[0]), (sprintf "%0.2f", $totalBasesFinal[0]/$totalReadsFinal[0])],
			["Total number of reads", $totalReads[0], $totalReadsFinal[0]],
			["Total number of reads with non-ATGC bases", $readsWithN[0], $readsWithN[1]],
			["Percentage of reads with non-ATGC bases", (sprintf "%0.2f", $readsWithN[0]/$totalReads[0]*100)."%", (sprintf "%0.2f", $readsWithN[1]/$totalReadsFinal[0]*100)."%"],
			["Total number of bases", $totalBases[0], $totalBasesFinal[0]],
			["Total number of HQ bases", $totalHQBases[0], $totalHQBasesFinal[0]],
			["Percentage of HQ bases", (sprintf "%0.2f", $totalHQBases[0]/$totalBases[0]*100)."%", (sprintf "%0.2f", $totalHQBasesFinal[0]/$totalBasesFinal[0]*100)."%"],
			["Total number of non-ATGC bases", $totalNs[0], $totalNs[1]],
			["Percentage of non-ATGC bases", (sprintf "%0.2f", $totalNs[0]/$totalBases[0]*100)."%", (sprintf "%0.2f", $totalNs[1]/$totalBasesFinal[0]*100)."%"],
		);
		for(my $i=0; $i<@arr; $i++) {
			if(!defined($isOnlyStat)) {
				printf $STAT "$inde %-50s  %-20s  %s\n", $arr[$i][0], $arr[$i][1], $arr[$i][2];
			}
			else {
				printf $STAT "$inde %-50s  %s\n", $arr[$i][0], $arr[$i][1];				
			}
		}
		
		print $STAT "\n\n";

		@graphData1 = ();
		print $STAT "Average quality score at each base position of input reads\n\n";
		my $c = 0;
		foreach my $arr (@positionSpecificBaseCount) {
			print $STAT getFileName($fileName[$c]), "\n";
			my $basePos = 1;
			foreach my $val (@$arr) {
				my $outVal = sprintf "%0.2f", $val/$totalReads[$c];
				$graphData1[0][$basePos-1] = $basePos;
				$graphData1[1][$basePos-1] = $outVal;
				print $STAT $outVal, "\t";
				$basePos++;
			}
			$c++;
			print $STAT "\n";
		}
		print $STAT "\n\n";
		$c = 0;
		if(!defined($isOnlyStat)) {
			print $STAT "Average quality score at each base position of filtered reads\n\n";
			foreach my $arr (@positionSpecificBaseCountHQ) {
				print $STAT getFileName($fileName[$c]), "\n";
				my $basePos = 1;
				foreach my $val (@$arr) {
					my $outVal = sprintf "%0.2f", $val/$totalReadsFinal[$c];
					$graphData1[2][$basePos-1] = $outVal;
					print $STAT $outVal, "\t";
					$basePos++;
				}
				$c++;
				print $STAT "\n";
			}
		}
		print $STAT "\n\n";
		prepareData4RangeGraph($STAT);
		if($isGDMod) {
			drawGraph(\@graphData1, $outFolder.getFileName($fileName[0])."_avgQual.png", getFileName($fileName[0]));
		}
	}
}

sub printStatTab {
	my $STAT = $_[0];
	my $tmpPer;
	my $inde = "\t";
	print $STAT "Parameters\n";
	my @graphData1 = ();
	my @graphData2 = ();
	if($isPairedEnd) {
		printf $STAT "\t%s\t%s\n", "Library type", "Paired-end";
		printf $STAT "\t%s\t%s\t%s\n", "Input files", $fileName[0], $fileName[1];
		printf $STAT "\t%s\t%s\n", "Primer/Adaptor library", defined($priAdaLib)?(($priAdaLib ne "u")?$priAdaLibNames[$priAdaLib]:"User defined ($priAdaFile)"):"NA";
		printf $STAT "\t%s\t%s\n", "Cut-off read length for HQ", $cutOffReadLen4HQ."%";
		printf $STAT "\t%s\t%s\n", "Cut-off quality score", $cutOffPhScore;
		printf $STAT "\t%s\t%s\n", "Only statistics", defined($isOnlyStat)?"On":"Off";
		printf $STAT "\t%s\t%s\n", "Number of processes", $noOfProcesses;

		print $STAT "\n\n";

		print $STAT "QC statistics\n";
		printf $STAT "\t%s\t%s\t%s\n", "File name", getFileName($fileName[0]), getFileName($fileName[1]);
		printf $STAT "\t%s\t%d\t%d\n", "Total number of reads", $totalReads[0], $totalReads[1];
		printf $STAT "\t%s\t%d\t%d\n", "Total number of HQ reads", $totalReadsAfterHQ[0], $totalReadsAfterHQ[1];
		$tmpPer = sprintf "%0.2f", $totalReadsAfterHQ[0]/$totalReads[0]*100;
		printf $STAT "\t%s\t%s\t%0.2f%s\n", "Percentage of HQ reads", $tmpPer."%", $totalReadsAfterHQ[1]/$totalReads[1]*100, "%";
		printf $STAT "\t%s\t%d\t%d\n", "Total number of bases", $totalBases[0], $totalBases[1];
		printf $STAT "\t%s\t%d\t%d\n", "Total number of bases in HQ reads", $totalBasesAfterHQ[0], $totalBasesAfterHQ[1];
		printf $STAT "\t%s\t%d\t%d\n", "Total number of HQ bases in HQ reads", $totalHQBasesAfterHQ[0], $totalHQBasesAfterHQ[1];
		$tmpPer = sprintf "%0.2f", $totalHQBasesAfterHQ[0]/$totalBasesAfterHQ[0]*100;
		printf $STAT "\t%s\t%s\t%0.2f%s\n", "Percentage of HQ bases in HQ reads", $tmpPer."%" , $totalHQBasesAfterHQ[1]/$totalBasesAfterHQ[1]*100, "%";
		if(defined($priAdaLib)) {
			printf $STAT "\t%s\t%d\t%d\n", "Number of Primer/Adaptor contaminated HQ reads", $totalValidReadsWithPriAda[0], $totalValidReadsWithPriAda[1];
		}
		else {
			printf $STAT "\t%s\t%s\t%s\n", "Number of Primer/Adaptor contaminated HQ reads", "NA", "NA";
		}
		printf $STAT "\t%s\t%d\t%d\n", "Total number of HQ filtered reads", $totalReadsFinal[0], $totalReadsFinal[1];
		$tmpPer = sprintf "%0.2f", $totalReadsFinal[0]/$totalReads[0]*100;
		printf $STAT "\t%s\t%s\t%0.2f%s\n", "Percentage of HQ filtered reads", $tmpPer."%", $totalReadsFinal[1]/$totalReads[1]*100, "%";

		print $STAT "\n\n";

		print $STAT "Detailed QC statistics\n";
		my @arr = (
			["File name", getFileName($fileName[0]), getFileName($fileName[1]), getFileName($outFileName[0]), getFileName($outFileName[1])],
			["Minimum read length", $minLen[0], $minLen[1], $minLen[2], $minLen[3]],
			["Maximum read length", $maxLen[0], $maxLen[1], $maxLen[2], $maxLen[3]],
			["Average read length", (sprintf "%0.2f", $totalBases[0]/$totalReads[0]), (sprintf "%0.2f", $totalBases[1]/$totalReads[1]), (sprintf "%0.2f", $totalBasesFinal[0]/$totalReadsFinal[0]), (sprintf "%0.2f", $totalBasesFinal[1]/$totalReadsFinal[1])],
			["Total number of reads", $totalReads[0], $totalReads[1], $totalReadsFinal[0], $totalReadsFinal[1]],
			["Total number of reads with non-ATGC bases", $readsWithN[0], $readsWithN[1], $readsWithN[2], $readsWithN[3]],
			["Percentage of reads with non-ATGC bases", (sprintf "%0.2f", $readsWithN[0]/$totalReads[0]*100)."%", (sprintf "%0.2f", $readsWithN[1]/$totalReads[1]*100)."%", (sprintf "%0.2f", $readsWithN[2]/$totalReadsFinal[0]*100)."%", (sprintf "%0.2f", $readsWithN[3]/$totalReadsFinal[1]*100)."%"],
			["Total number of bases", $totalBases[0], $totalBases[1], $totalBasesFinal[0], $totalBasesFinal[1]],
			["Total number of HQ bases", $totalHQBases[0], $totalHQBases[1], $totalHQBasesFinal[0], $totalHQBasesFinal[1]],
			["Percentage of HQ bases", (sprintf "%0.2f", $totalHQBases[0]/$totalBases[0]*100)."%", (sprintf "%0.2f", $totalHQBases[1]/$totalBases[1]*100)."%", (sprintf "%0.2f", $totalHQBasesFinal[0]/$totalBasesFinal[0]*100)."%", (sprintf "%0.2f", $totalHQBasesFinal[1]/$totalBasesFinal[1]*100)."%"],
			["Total number of non-ATGC bases", $totalNs[0], $totalNs[1], $totalNs[2], $totalNs[3]],
			["Percentage of non-ATGC bases", (sprintf "%0.2f", $totalNs[0]/$totalBases[0]*100)."%", (sprintf "%0.2f", $totalNs[1]/$totalBases[1]*100)."%", (sprintf "%0.2f", $totalNs[2]/$totalBasesFinal[0]*100)."%", (sprintf "%0.2f", $totalNs[3]/$totalBasesFinal[1]*100)."%"],
		);
		for(my $i=0; $i<@arr; $i++) {
			if(!defined($isOnlyStat)) {
				printf $STAT "\t%s\t%s\t%s\t%s\t%s\n", $arr[$i][0], $arr[$i][1], $arr[$i][2], $arr[$i][3], $arr[$i][4];
			}
			else {
				printf $STAT "\t%s\t%s\t%s\n", $arr[$i][0], $arr[$i][1], $arr[$i][2];				
			}
		}
		
		print $STAT "\n\n";
		
		print $STAT "Average quality score at each base position of input reads\n\n";
		my $c = 0;
		@graphData1 = ();
		@graphData2 = ();
		foreach my $arr (@positionSpecificBaseCount) {
			print $STAT getFileName($fileName[$c]), "\n";
			my $basePos = 1;
			foreach my $val (@$arr) {
				my $outVal = sprintf "%0.2f", $val/$totalReads[$c];
				if($c == 0) {
					$graphData1[0][$basePos-1] = $basePos;
					$graphData1[1][$basePos-1] = $outVal;
				}
				else {
					$graphData2[0][$basePos-1] = $basePos;
					$graphData2[1][$basePos-1] = $outVal;
				}
				print $STAT $outVal, "\t";
				$basePos++;
			}
			print $STAT "\n\n";
			$c++;
		}
		print $STAT "\n\n";
		if(!defined($isOnlyStat)) {
			print $STAT "Average quality score at each base position of filtered reads\n\n";
			$c = 0;
			foreach my $arr (@positionSpecificBaseCountHQ) {
				print $STAT getFileName($outFileName[$c]), "\n";
				my $basePos = 1;
				foreach my $val (@$arr) {
					my $outVal = sprintf "%0.2f", $val/$totalReadsFinal[$c];
					if($c == 0) {
						$graphData1[2][$basePos-1] = $outVal;
					}
					else {
						$graphData2[2][$basePos-1] = $outVal;
					}
					print $STAT $outVal, "\t";
					$basePos++;
				}
				$c++;
				print $STAT "\n\n";
			}
		}
		print $STAT "\n\n";
		prepareData4RangeGraph($STAT);
		if($isGDMod) {
			drawGraph(\@graphData1, $outFolder.getFileName($fileName[0])."_avgQual.png", getFileName($fileName[0]));
			drawGraph(\@graphData2, $outFolder.getFileName($fileName[1])."_avgQual.png", getFileName($fileName[1]));
		}
	}
	else {
		printf $STAT "\t%s\t%s\n", "Library type", "Single-end";
		printf $STAT "\t%s\t%s\n", "Input file", $fileName[0];
		printf $STAT "\t%s\t%s\n", "Primer/Adaptor library", defined($priAdaLib)?(($priAdaLib ne "u")?$priAdaLibNames[$priAdaLib]:"User defined ($priAdaFile)"):"NA";
		printf $STAT "\t%s\t%s\n", "Cut-off read length for HQ", $cutOffReadLen4HQ."%";
		printf $STAT "\t%s\t%s\n", "Cut-off quality score", $cutOffPhScore;
		printf $STAT "\t%s\t%s\n", "Only statistics", defined($isOnlyStat)?"On":"Off";
		printf $STAT "\t%s\t%s\n", "Number of processes", $noOfProcesses;

		print $STAT "\n\n";

		print $STAT "QC statistics\n";
		printf $STAT "\t%s\t%s\n", "File name", getFileName($fileName[0]);
		printf $STAT "\t%s\t%d\n", "Total number of reads", $totalReads[0];
		printf $STAT "\t%s\t%d\n", "Total number of HQ reads", $totalReadsAfterHQ[0];
		$tmpPer = sprintf "%0.2f", $totalReadsAfterHQ[0]/$totalReads[0]*100;
		printf $STAT "\t%s\t%s\n", "Percentage of HQ reads", $tmpPer."%";
		printf $STAT "\t%s\t%d\n", "Total number of bases", $totalBases[0];
		printf $STAT "\t%s\t%d\n", "Total number of bases in HQ reads", $totalBasesAfterHQ[0];
		printf $STAT "\t%s\t%d\n", "Total number of HQ bases in HQ reads", $totalHQBasesAfterHQ[0];
		$tmpPer = sprintf "%0.2f", $totalHQBasesAfterHQ[0]/$totalBasesAfterHQ[0]*100;
		printf $STAT "\t%s\t%s\n", "Percentage of HQ bases in HQ reads", $tmpPer."%";
		if(defined($priAdaLib)) {
			printf $STAT "\t%s\t%d\n", "Number of Primer/Adaptor contaminated HQ reads", $totalValidReadsWithPriAda[0];
		}
		else {
			printf $STAT "\t%s\t%s\n", "Number of Primer/Adaptor contaminated HQ reads", "NA";
		}
		printf $STAT "\t%s\t%d\n", "Total number of HQ filtered reads", $totalReadsFinal[0];
		$tmpPer = sprintf "%0.2f", $totalReadsFinal[0]/$totalReads[0]*100;
		printf $STAT "\t%s\t%s\n", "Percentage of HQ filtered reads", $tmpPer."%";

		print $STAT "\n\n";

		print $STAT "Detailed QC statistics\n";
		my @arr = (
			["File name", getFileName($fileName[0]), getFileName($outFileName[0])],
			["Minimum read length", $minLen[0], $minLen[1]],
			["Maximum read length", $maxLen[0], $maxLen[1]],
			["Average read length", (sprintf "%0.2f", $totalBases[0]/$totalReads[0]), (sprintf "%0.2f", $totalBasesFinal[0]/$totalReadsFinal[0])],
			["Total number of reads", $totalReads[0], $totalReadsFinal[0]],
			["Total number of reads with non-ATGC bases", $readsWithN[0], $readsWithN[1]],
			["Percentage of reads with non-ATGC bases", (sprintf "%0.2f", $readsWithN[0]/$totalReads[0]*100)."%", (sprintf "%0.2f", $readsWithN[1]/$totalReadsFinal[0]*100)."%"],
			["Total number of bases", $totalBases[0], $totalBasesFinal[0]],
			["Total number of HQ bases", $totalHQBases[0], $totalHQBasesFinal[0]],
			["Percentage of HQ bases", (sprintf "%0.2f", $totalHQBases[0]/$totalBases[0]*100)."%", (sprintf "%0.2f", $totalHQBasesFinal[0]/$totalBasesFinal[0]*100)."%"],
			["Total number of non-ATGC bases", $totalNs[0], $totalNs[1]],
			["Percentage of non-ATGC bases", (sprintf "%0.2f", $totalNs[0]/$totalBases[0]*100)."%", (sprintf "%0.2f", $totalNs[1]/$totalBasesFinal[0]*100)."%"],
		);
		for(my $i=0; $i<@arr; $i++) {
			if(!defined($isOnlyStat)) {
				printf $STAT "\t%s\t%s\t%s\n", $arr[$i][0], $arr[$i][1], $arr[$i][2];
			}
			else {
				printf $STAT "\t%s\t%s\n", $arr[$i][0], $arr[$i][1];				
			}
		}
		
		print $STAT "\n\n";


		@graphData1 = ();
		print $STAT "Average quality score at each base position of input reads\n\n";
		my $c = 0;
		foreach my $arr (@positionSpecificBaseCount) {
			print $STAT getFileName($fileName[$c]), "\n";
			my $basePos = 1;
			foreach my $val (@$arr) {
				my $outVal = sprintf "%0.2f", $val/$totalReads[$c];
				$graphData1[0][$basePos-1] = $basePos;
				$graphData1[1][$basePos-1] = $outVal;
				print $STAT $outVal, "\t";
				$basePos++;
			}
			$c++;
			print $STAT "\n";
		}
		print $STAT "\n\n";
		$c = 0;
		if(!defined($isOnlyStat)) {
			print $STAT "Average quality score at each base position of filtered reads\n\n";
			foreach my $arr (@positionSpecificBaseCountHQ) {
				print $STAT getFileName($fileName[$c]), "\n";
				my $basePos = 1;
				foreach my $val (@$arr) {
					my $outVal = sprintf "%0.2f", $val/$totalReadsFinal[$c];
					$graphData1[2][$basePos-1] = $outVal;
					print $STAT $outVal, "\t";
					$basePos++;
				}
				$c++;
				print $STAT "\n";
			}
		}
		print $STAT "\n\n";
		prepareData4RangeGraph($STAT);
		if($isGDMod) {
			drawGraph(\@graphData1, $outFolder.getFileName($fileName[0])."_avgQual.png", getFileName($fileName[0]));
		}
	}
}
