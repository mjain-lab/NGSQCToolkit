#! /usr/bin/perl

use File::Basename;
#BEGIN {
#	my ($tmp, $path) = fileparse($0);
#	push ( @INC,"$path/lib");
#	#use lib "$path";
#}
use strict;
use warnings;
use List::Util qw(sum min max);
use Getopt::Long;
use Cwd qw(abs_path);
use IO::Zlib;
use FindBin qw($RealBin);
use lib "$RealBin/lib";
require "454PEhtml.pl";

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


# Setting parameters
my $lowestValidLen = 50;
my @files = ();
my $noOfInp = 3;
my $helpAsked;
my $cutOffReadLen4HQ = 70;
my $cutOffPhScore = 20;
my $outFolder = "";
my $isOnlyStat;
my $statOutFmt = 1;
my $noOfProcesses = 1;
my $homoPolyLen = 0;
my $priAdaLib;
my $isLenFilterOn = 1;
my @priAdaLibNames = ("Rapid Library (Standard)", "Paired End Library", "Amplicon PE Library", "Small RNA Library");
my $priAdaFile;
my @usrDefinedPriAda = ();
my $outputDataFmt = "t";		# t/T: Text; g/G: Gzip.
my $linker;
my $linkerF = "GTTGGAACCGAAAGGGTTTGAATTCAAACCCTTTCGGTTCCAAC";
my $linkerR;

GetOptions(
			"i=s{$noOfInp}" => \@files,
			"h|help" => \$helpAsked,
			"l|cutOffReadLen4HQ=f" => \$cutOffReadLen4HQ,
			"n|homoPolyLen=i" => \$homoPolyLen,
			"o|outputFolder=s" => \$outFolder,
			"z|outputDataCompression=s" => \$outputDataFmt,
			"t|statOutFmt=i" => \$statOutFmt,
			"onlyStat" => \$isOnlyStat,
			"p|processes=i" => \$noOfProcesses,
			"s|cutOffQualScore=i" => \$cutOffPhScore,
			"m|minLen=i" => \$lowestValidLen,
			"f|lenFilter=s" => \$isLenFilterOn,
			"linker=s" => \$linkerF,
		  );
if(defined($helpAsked)) {
	prtUsage();
	exit;
}
if(@files == 0) {
	prtError("No input files are provided");
}
my @tempFiles = ();
prtError("Missing inputs for option -i") if((scalar @files)%$noOfInp != 0);
for(my $i=0; $i<@files; $i+=$noOfInp) {
	my $str = "$files[$i] $files[$i+1] $files[$i+2]";
	if($files[$i+2] =~ /^-/) {
		prtError("Missing inputs for option -i: at '-i $str'")
	}
	if($files[$i+2] =~ /^\d$/) {
		if($files[$i+2] < 1 || $files[$i+2] > 4) {
			prtError("Incorrect option for Primer/Adaptor library: at '-i $str'");
		}
	}
	push(@tempFiles, $str);
}
@files = ();
@files = @tempFiles;
if($cutOffReadLen4HQ < 0 || $cutOffReadLen4HQ > 100) {
	prtError("Incorrect value for -l|cutOffReadLen4HQ option: at '-l $cutOffReadLen4HQ'");
}
if($cutOffPhScore < 0 || $cutOffPhScore > 40) {
	prtError("Incorrect value for -s|cutOffPhScore option: at '-s $cutOffPhScore'");
}
if($statOutFmt < 1 || $statOutFmt > 2) {
	prtError("Incorrect value for -statOutFmt: at '-statOutFmt $statOutFmt'");
}
if($isLenFilterOn =~ /^N/i) {
	$isLenFilterOn = 0;
}
else {
	$isLenFilterOn = 1;
}
if($outputDataFmt !~ /^[tg]$/i) {
	prtError("Incorrect value for -f|outputDataFmt option: at '-f $outputDataFmt'");
}

my $tmpLinker = $linkerF;
$tmpLinker =~ tr/ATGCatgc/TACGtacg/;
$tmpLinker = reverse $tmpLinker;
if($linkerF ne $tmpLinker) {
	$linkerR = $tmpLinker;
}

my $pm = new Parallel::ForkManager($noOfProcesses);


my $trimCount = 0;
my @trimCountPE = (0, 0, 0);
my $seqCount = 0;
my $seqCountPE = 0;
my $ttlSeqCount = 0;
my $lt100 = 0;
my @lt100PE = (0, 0, 0);
my $hQCount = 0;
my @hQCountPE = (0, 0, 0);
my $lQCount = 0;
my @lQCountPE = (0, 0, 0);
my $maxRawLen = 0;
my $minRawLen = 1000000000000;
my @maxRawLenPE = (0, 0);
my @minRawLenPE = (1000000000000, 1000000000000);
my $avgRawLen = 0;
my $maxHQLen = 0;
my $minHQLen = 1000000000000;
my $avgHQLen = 0;
my @rawLen = ();
my @rawLenPE = ();
my @hQLen = ();
my $totalBases = 0;
my $totalHQBases = 0;
my @totalBasesPE = (0, 0, 0);
my @totalHQBasesPE = (0, 0, 0);
my $totalBasesUPOri = 0;
my $totalBasesAfterHQ = 0;
my @totalBasesAfterHQPE = (0, 0, 0);
my $totalHQBasesAfterHQ = 0;
my @totalHQBasesAfterHQPE = (0, 0, 0);
my $totalBasesFinal = 0;
my $totalHQBasesFinal = 0;
my $totalReadsFinal = 0;
my @totalReadsFinalPE = (0, 0, 0);
my @totalReadsFinalUP = (0, 0);
my $avgQual = 0;
my @avgQualPE = (0, 0);
my $avgQualFinal = 0;
my $totalValidReadsWithPriAda = 0;
my @totalValidReadsWithPriAdaPE = (0, 0, 0);
my $totalValidReadsNoPriAda = 0;
my $substrlen = 20;						# For removePriAda
my $mismLim = 1;						# For removePriAda

my $fastaSeqId = "";
my $fastaSeq = "";
my $qualSeqId = "";
my $qualSeq = "";
my $prevFastaSeqId = "";
my $indOfAnalysis = 0;

my @lenDistrib = ();
my $lenInterval = 40;
my @qualDistrib = ();
my $qualInterval = 1;
my @gcDistrib = ();
my $gcInterval = 5;
my @charCount = ();
my @charCountPE = ();

my $font_spec = getFilePath($0) . "lib/Fonts/Dustismo_Sans.ttf";
my $f = getFilePath($0) . "lib/Fonts/LucidaSansDemiBold.ttf";


foreach my $inpData (@files) {
	$indOfAnalysis++;
my $pid = $pm->start and next;
	$inpData =~ s/\\([A-Za-z_\.])/\/$1/g;		# To remove '\' from the path of windows file
	my @iData = split(" ", $inpData);
	my $seqFile = $iData[0];
	my $qualFile = $iData[1];
	$priAdaLib = $iData[2];
	print "Analysis has been started for \"$seqFile\": Index: $indOfAnalysis\n";
	if($priAdaLib =~ /^n$/i) {
		undef $priAdaLib;
	}
	elsif($priAdaLib =~ /^\d$/) {
		$priAdaLib = $priAdaLib - 1;
	}
	else {
		$priAdaFile = $priAdaLib;
		$priAdaLib = "u";
		open(PRIADA, "<$priAdaFile") or die "Can not open the user-defined primer/adapter file: $priAdaFile\n";
		@usrDefinedPriAda = <PRIADA>;
		for(my $i=0; $i<$#usrDefinedPriAda; $i++) {
			$usrDefinedPriAda[$i] =~ s/\s+//g;
		}
	}
	my ($seqFileName, $filePath) = fileparse($seqFile);
	my ($qualFileName) = fileparse($qualFile);
	$outFolder = $filePath . "454QC_Filtered_files" if($outFolder eq "");
	$outFolder .= "/" if($outFolder !~ /\/$/);
	if(! -e $outFolder) {
			mkdir($outFolder) or die "Can not create output folder: $outFolder\n";
	}
	my $outSeqFile = $outFolder . $seqFileName . "_filtered";
	my $outQualFile = $outFolder . $qualFileName . "_filtered";
	$outSeqFile .= ".gz" if($outputDataFmt =~ /g/i);
	$outQualFile .= ".gz" if($outputDataFmt =~ /g/i);
	my $statFile = $outFolder . $seqFileName . "_stat";

	my $iH;
	openFileGetHandle($seqFile, "r", \$iH);
	*I = $iH;
	my $qH;
	openFileGetHandle($qualFile, "r", \$qH);
	*Q = $qH;
	if(!defined($isOnlyStat)) {
		my $oiH;
		openFileGetHandle($outSeqFile, "w", \$oiH);
		*OI = $oiH;
		my $oqH;
		openFileGetHandle($outQualFile, "w", \$oqH);
		*OQ = $oqH;
	}
	open(STAT, ">$statFile") or die "Can not open file: $statFile\n";
	while(my $line = <I>) {
		$ttlSeqCount++ if($line =~ /^>/);
	}
	close(I);
	print "$indOfAnalysis: Number of reads processed: " . "0/$ttlSeqCount (0\%)...\n";
	undef($iH);
	openFileGetHandle($seqFile, "r", \$iH);
	*I = $iH;
	
	while(my $line = <I>) {
		chomp $line;
		my $qualLine = <Q>;
		chomp($qualLine);
		if($line =~ /^>/) {
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
	
	print "$indOfAnalysis: Number of reads processed: " . "$ttlSeqCount/$ttlSeqCount (100\%)...\n";
	print "$indOfAnalysis: Analysis completed\n";

	print "$indOfAnalysis: Printing Statistics...\n";

	if($statOutFmt == 1) {
		my $inde = " " x 1;
		my ($tmpPer, $tmpPer2, $tmpPer3);
		printf STAT "Parameters\n";
		printf STAT "$inde %-40s %s  %s\n", "Input files ", $seqFile, $qualFile;
		printf STAT "$inde %-40s %s\n", "Primer/Adaptor library", defined($priAdaLib)?(($priAdaLib ne "u")?$priAdaLibNames[$priAdaLib]:"User defined ($priAdaFile)"):"NA";
		printf STAT "$inde %-40s %s\n", "Linker sequence", (($linkerR)?"(+)5'$linkerF 3'/ (-)5'$linkerR 3'":"(+)5'$linkerF 3' / (-)5'$linkerF 3'");
		printf STAT "$inde %-40s %s\n", "Homopolymer trimming", "Off" if($homoPolyLen == 0);
		printf STAT "$inde %-40s %s\n", "Homopolymer trimming", "On" if($homoPolyLen != 0);
		printf STAT "$inde %-40s %s\n", "Length of the homopolymer to be removed", $homoPolyLen if($homoPolyLen != 0);
		printf STAT "$inde %-40s %s\n", "Length filter", ($isLenFilterOn)?"On":"Off";
		printf STAT "$inde %-40s %s\n", "Cut-off for minimum read length", $lowestValidLen if($isLenFilterOn);
		printf STAT "$inde %-40s %s\n", "Cut-off read length for HQ", $cutOffReadLen4HQ."%";
		printf STAT "$inde %-40s %s\n", "Cut-off quality score", $cutOffPhScore;
		printf STAT "$inde %-40s %s\n", "Only statistics", defined($isOnlyStat)?"On":"Off";
		printf STAT "$inde %-40s %s\n", "Number of processes", $noOfProcesses;
	
		print STAT "\n\n";
	
		print STAT "QC statistics\n";
		printf STAT "$inde %-70s %s\n", "File name", $seqFileName;
		printf STAT "$inde %-70s %d\n", "Total number of reads", $seqCount+$seqCountPE;
		print STAT "\n";
		printf STAT "$inde %-70s %-13s %s\n", "QC analysis of Paired reads:", "Paired", "(Read1 / Read2)";
		printf STAT "$inde %-70s %-13d (%d / %d)\n", "Total number of Paired reads", $seqCountPE, $seqCountPE, $seqCountPE;
		printf STAT "$inde %-70s %-13d (%d / %d)\n", "Total number of trimmed reads containing homopolymer", $trimCountPE[2], $trimCountPE[0], $trimCountPE[1];
		printf STAT "$inde %-70s %-13d (%d / %d)\n", "Total number of trashed reads (<$lowestValidLen bp in length after trimming)", $lt100PE[2], $lt100PE[0], $lt100PE[1];
		printf STAT "$inde %-70s %-13d (%d / %d)\n", "Total number of low quality reads (excluding <$lowestValidLen reads)", $lQCountPE[2], $lQCountPE[0], $lQCountPE[1];
		printf STAT "$inde %-70s %-13d (%d / %d)\n", "Total number of HQ reads", $hQCountPE[2], $hQCountPE[0], $hQCountPE[1];
		$tmpPer = ($seqCountPE)?(sprintf "%0.2f", $hQCountPE[2]/$seqCountPE*100):"0";
		$tmpPer2 = ($seqCountPE)?(sprintf "%0.2f", $hQCountPE[0]/$seqCountPE*100):"0";
		$tmpPer3 = ($seqCountPE)?(sprintf "%0.2f", $hQCountPE[1]/$seqCountPE*100):"0";
		printf STAT "$inde %-70s %-13s (%s / %s)\n", "Percentage of HQ reads", $tmpPer."%", $tmpPer2."%", $tmpPer3."%";
		printf STAT "$inde %-70s %-13d (%d / %d)\n", "Total number of bases", $totalBasesPE[2], $totalBasesPE[0], $totalBasesPE[1];
		printf STAT "$inde %-70s %-13d (%d / %d)\n", "Total number of bases in HQ reads", $totalBasesAfterHQPE[2], $totalBasesAfterHQPE[0], $totalBasesAfterHQPE[1];
		printf STAT "$inde %-70s %-13d (%d / %d)\n", "Total number of HQ bases in HQ reads", $totalHQBasesAfterHQPE[2], $totalHQBasesAfterHQPE[0], $totalHQBasesAfterHQPE[1];
		$tmpPer = ($totalBasesAfterHQPE[2])?(sprintf "%0.2f", $totalHQBasesAfterHQPE[2]/$totalBasesAfterHQPE[2]*100):"0";
		$tmpPer2 = ($totalBasesAfterHQPE[2])?(sprintf "%0.2f", $totalHQBasesAfterHQPE[0]/$totalBasesAfterHQPE[0]*100):"0";
		$tmpPer3 = ($totalBasesAfterHQPE[2])?(sprintf "%0.2f", $totalHQBasesAfterHQPE[1]/$totalBasesAfterHQPE[1]*100):"0";
		printf STAT "$inde %-70s %-13s (%s / %s)\n", "Percentage of HQ bases in HQ reads", $tmpPer."%", $tmpPer2."%", $tmpPer3."%";
		if(defined($priAdaLib)) {
			printf STAT "$inde %-70s %-13d (%d / %d)\n", "Number of Primer/Adaptor trimmed reads", $totalValidReadsWithPriAdaPE[2], $totalValidReadsWithPriAdaPE[0], $totalValidReadsWithPriAdaPE[1];
		}
		else {
			printf STAT "$inde %-70s %-13s (%s / %s)\n", "Number of Primer/Adaptor trimmed reads", "NA", "NA", "NA";
		}
		printf STAT "$inde %-70s %-13d (%d / %d)\n", "Total number of HQ filtered reads", $totalReadsFinalPE[2], $totalReadsFinalPE[0], $totalReadsFinalPE[1];
		$tmpPer = ($seqCountPE)?(sprintf "%0.2f", $totalReadsFinalPE[2]/$seqCountPE*100):"0";
		$tmpPer2 = ($seqCountPE)?(sprintf "%0.2f", $totalReadsFinalPE[0]/$seqCountPE*100):"0";
		$tmpPer3 = ($seqCountPE)?(sprintf "%0.2f", $totalReadsFinalPE[1]/$seqCountPE*100):"0";
		printf STAT "$inde %-70s %-13s (%s / %s)\n", "Percentage of HQ filtered reads", $tmpPer."%", $tmpPer2."%", $tmpPer3."%";
		printf STAT "$inde %-70s %d\n", "Total number of HQ filtered reads (Unpaired)", $totalReadsFinalUP[0];
		$tmpPer = ($seqCountPE)?(sprintf "%0.2f", $totalReadsFinalUP[0]/$seqCountPE*100):"0";
		printf STAT "$inde %-70s %s\n", "Percentage of HQ filtered reads (Unpaired)", $tmpPer."%";

		print STAT "\n";
		printf STAT "$inde %-70s %-13s\n", "QC analysis of Unpaired (UPOri*) reads:", "Total";
		printf STAT "$inde %-70s %-13d\n", "Total number of Unpaired reads", $seqCount;
		printf STAT "$inde %-70s %-13d\n", "Total number of trimmed reads containing homopolymer", $trimCount;
		printf STAT "$inde %-70s %-13d\n", "Total number of trashed reads (<$lowestValidLen bp in length after trimming)", $lt100;
		printf STAT "$inde %-70s %-13d\n", "Total number of low quality reads (excluding <$lowestValidLen reads)", $lQCount;
		printf STAT "$inde %-70s %-13d\n", "Total number of HQ reads", $hQCount;
		$tmpPer = sprintf "%0.2f", $hQCount/$seqCount*100;
		printf STAT "$inde %-70s %-13s\n", "Percentage of HQ reads", $tmpPer."%";
		printf STAT "$inde %-70s %-13d\n", "Total number of bases", $totalBasesUPOri;
		printf STAT "$inde %-70s %-13d\n", "Total number of bases in HQ reads", $totalBasesAfterHQ;
		printf STAT "$inde %-70s %-13d\n", "Total number of HQ bases in HQ reads", $totalHQBasesAfterHQ;
		$tmpPer = sprintf "%0.2f", $totalHQBasesAfterHQ/$totalBasesAfterHQ*100;
		printf STAT "$inde %-70s %-13s\n", "Percentage of HQ bases in HQ reads", $tmpPer."%";
		if(defined($priAdaLib)) {
			printf STAT "$inde %-70s %-13d\n", "Number of Primer/Adaptor trimmed reads", $totalValidReadsWithPriAda;
		}
		else {
			printf STAT "$inde %-70s %-13s\n", "Number of Primer/Adaptor trimmed reads", "NA";
		}
		printf STAT "$inde %-70s %-13d\n", "Total number of HQ filtered reads (Unpaired)", $totalReadsFinal;
		$tmpPer = sprintf "%0.2f", $totalReadsFinal/$seqCount*100;
		printf STAT "$inde %-70s %-13s\n", "Percentage of HQ filtered reads (Unpaired)", $tmpPer."%";
		print STAT ("-"x100)."\n";
		$tmpPer = sprintf "%0.2f", $totalReadsFinalPE[2]/($seqCount+$seqCountPE)*100;
		$tmpPer2 = sprintf "%0.2f", $totalReadsFinalUP[0]/($seqCount+$seqCountPE)*100;
		$tmpPer3 = sprintf "%0.2f", $totalReadsFinal/($seqCount+$seqCountPE)*100;
		printf STAT "$inde %-70s %-13d (%s)\n", "Number of HQ filtered reads (Paired)", $totalReadsFinalPE[2], $tmpPer."%";
		printf STAT "$inde %-70s %-13d (%s)\n", "Number of HQ filtered reads (UPPair*)", $totalReadsFinalUP[0], $tmpPer2."%";
		printf STAT "$inde %-70s %-13d (%s)\n", "Number of HQ filtered reads (UPOri*)", $totalReadsFinal, $tmpPer3."%";
		print STAT ("-"x100)."\n";
		$tmpPer = sprintf "%0.2f", ($totalReadsFinalPE[2]+$totalReadsFinalUP[0]+$totalReadsFinal)/($seqCount+$seqCountPE)*100;
		printf STAT "$inde %-70s %-13d (%s)\n", "Total number of HQ filtered reads", ($totalReadsFinalPE[2]+$totalReadsFinalUP[0]+$totalReadsFinal), $tmpPer."%";
		print STAT "\n";
		print STAT "* UPOri: Unpaired reads (i.e. reads without linker sequence) found in the input file\n";
		print STAT "  UPPair: One of the paired reads which passed QC\n";

		print STAT "\n\n";

######### Adding the statistics of PE, UPPair and UPOri reads
		$seqCount = $seqCount + $seqCountPE;
		$totalReadsFinal = $totalReadsFinalPE[2]+$totalReadsFinalUP[0]+$totalReadsFinal;
######### Done adding
		print STAT "Detailed QC statistics\n";
		my @arr = (
			["File name", $seqFileName, (fileparse($outSeqFile))[0]],
			["Total number of reads", $seqCount, $totalReadsFinal],
			["Minimum read length", $minRawLen, $minHQLen],
			["Maximum read length", $maxRawLen, $maxHQLen],
			["Average read length", (sprintf "%0.2f", $totalBases/$seqCount), (sprintf "%0.2f", $totalBasesFinal/$totalReadsFinal)],
			["Median read length", calcMedian(@rawLen), calcMedian(@hQLen)],
			["N25 length", calcN50(\@rawLen, 25), calcN50(\@hQLen, 25)],
			["N50 length", calcN50(\@rawLen, 50), calcN50(\@hQLen, 50)],
			["N75 length", calcN50(\@rawLen, 75), calcN50(\@hQLen, 75)],
			["N90 length", calcN50(\@rawLen, 90), calcN50(\@hQLen, 90)],
			["N95 length", calcN50(\@rawLen, 95), calcN50(\@hQLen, 95)],
			["Total number of bases", $totalBases, $totalBasesFinal],
			["Total number of HQ bases", $totalHQBases, $totalHQBasesFinal],
			["Percentage of HQ bases", (sprintf "%0.2f", $totalHQBases/$totalBases*100)."%", (sprintf "%0.2f", $totalHQBasesFinal/$totalBasesFinal*100)."%"],
			["Average quality score (Overall)", (sprintf "%0.2f", $avgQual/$seqCount), (sprintf "%0.2f", $avgQualFinal/$totalReadsFinal)],
		);
		for(my $i=0; $i<@arr; $i++) {
			if(!defined($isOnlyStat)) {
				printf STAT "$inde %-50s  %-20s  %s\n", $arr[$i][0], $arr[$i][1], $arr[$i][2];
			}
			else {
				printf STAT "$inde %-50s  %s\n", $arr[$i][0], $arr[$i][1];				
			}
		}
	
		print STAT "\n\n";
######### Subtracting the values back
		$seqCount = $seqCount - $seqCountPE;
		$totalReadsFinal = $totalReadsFinal-$totalReadsFinalPE[2]-$totalReadsFinalUP[0];
######### Done
	}
	elsif($statOutFmt == 2) {
		my ($tmpPer, $tmpPer2, $tmpPer3);
		printf STAT "Parameters\n";
		printf STAT "\t%s\t%s\t%s\n", "Input files ", $seqFile, $qualFile;
		printf STAT "\t%s\t%s\n", "Primer/Adaptor library", defined($priAdaLib)?(($priAdaLib ne "u")?$priAdaLibNames[$priAdaLib]:"User defined ($priAdaFile)"):"NA";
		printf STAT "\t%s\t%s\n", "Linker sequence", (($linkerR)?"(+)5'$linkerF 3'/ (-)5'$linkerR 3'":"(+)5'$linkerF 3' / (-)5'$linkerF 3'");
		printf STAT "\t%s\t%s\n", "Homopolymer trimming", "Off" if($homoPolyLen == 0);
		printf STAT "\t%s\t%s\n", "Homopolymer trimming", "On" if($homoPolyLen != 0);
		printf STAT "\t%s\t%s\n", "Length of the homopolymer to be removed", $homoPolyLen if($homoPolyLen != 0);
		printf STAT "\t%s\t%s\n", "Length filter", ($isLenFilterOn)?"On":"Off";
		printf STAT "\t%s\t%s\n", "Cut-off for minimum read length", $lowestValidLen if($isLenFilterOn);
		printf STAT "\t%s\t%s\n", "Cut-off read length for HQ", $cutOffReadLen4HQ."%";
		printf STAT "\t%s\t%s\n", "Cut-off quality score", $cutOffPhScore;
		printf STAT "\t%s\t%s\n", "Only statistics", defined($isOnlyStat)?"On":"Off";
		printf STAT "\t%s\t%s\n", "Number of processes", $noOfProcesses;
	
		print STAT "\n\n";
	
		print STAT "QC statistics\n";
		printf STAT "\t%s\t%s\n", "File name", $seqFileName;
		printf STAT "\t%s\t%d\n", "Total number of reads", $seqCount+$seqCountPE;
		print STAT "\n";
		printf STAT "\t%s\t%s\t%s\n", "QC analysis of Paired reads:", "Paired", "(Read1 / Read2)";
		printf STAT "\t%s\t%d\t(%d / %d)\n", "Total number of Paired reads", $seqCountPE, $seqCountPE, $seqCountPE;
		printf STAT "\t%s\t%d\t(%d / %d)\n", "Total number of trimmed reads containing homopolymer", $trimCountPE[2], $trimCountPE[0], $trimCountPE[1];
		printf STAT "\t%s\t%d\t(%d / %d)\n", "Total number of trashed reads (<$lowestValidLen bp in length after trimming)", $lt100PE[2], $lt100PE[0], $lt100PE[1];
		printf STAT "\t%s\t%d\t(%d / %d)\n", "Total number of low quality reads (excluding <$lowestValidLen reads)", $lQCountPE[2], $lQCountPE[0], $lQCountPE[1];
		printf STAT "\t%s\t%d\t(%d / %d)\n", "Total number of HQ reads", $hQCountPE[2], $hQCountPE[0], $hQCountPE[1];
		$tmpPer = ($seqCountPE)?(sprintf "%0.2f", $hQCountPE[2]/$seqCountPE*100):"0";
		$tmpPer2 = ($seqCountPE)?(sprintf "%0.2f", $hQCountPE[0]/$seqCountPE*100):"0";
		$tmpPer3 = ($seqCountPE)?(sprintf "%0.2f", $hQCountPE[1]/$seqCountPE*100):"0";
		printf STAT "\t%s\t%s\t(%s / %s)\n", "Percentage of HQ reads", $tmpPer."%", $tmpPer2."%", $tmpPer3."%";
		printf STAT "\t%s\t%d\t(%d / %d)\n", "Total number of bases", $totalBasesPE[2], $totalBasesPE[0], $totalBasesPE[1];
		printf STAT "\t%s\t%d\t(%d / %d)\n", "Total number of bases in HQ reads", $totalBasesAfterHQPE[2], $totalBasesAfterHQPE[0], $totalBasesAfterHQPE[1];
		printf STAT "\t%s\t%d\t(%d / %d)\n", "Total number of HQ bases in HQ reads", $totalHQBasesAfterHQPE[2], $totalHQBasesAfterHQPE[0], $totalHQBasesAfterHQPE[1];
		$tmpPer = ($totalBasesAfterHQPE[2])?(sprintf "%0.2f", $totalHQBasesAfterHQPE[2]/$totalBasesAfterHQPE[2]*100):"0";
		$tmpPer2 = ($totalBasesAfterHQPE[2])?(sprintf "%0.2f", $totalHQBasesAfterHQPE[0]/$totalBasesAfterHQPE[0]*100):"0";
		$tmpPer3 = ($totalBasesAfterHQPE[2])?(sprintf "%0.2f", $totalHQBasesAfterHQPE[1]/$totalBasesAfterHQPE[1]*100):"0";
		printf STAT "\t%s\t%s\t(%s / %s)\n", "Percentage of HQ bases in HQ reads", $tmpPer."%", $tmpPer2."%", $tmpPer3."%";
		if(defined($priAdaLib)) {
			printf STAT "\t%s\t%d\t(%d / %d)\n", "Number of Primer/Adaptor trimmed reads", $totalValidReadsWithPriAdaPE[2], $totalValidReadsWithPriAdaPE[0], $totalValidReadsWithPriAdaPE[1];
		}
		else {
			printf STAT "\t%s\t%s\t(%s / %s)\n", "Number of Primer/Adaptor trimmed reads", "NA", "NA", "NA";
		}
		printf STAT "\t%s\t%d\t(%d / %d)\n", "Total number of HQ filtered reads", $totalReadsFinalPE[2], $totalReadsFinalPE[0], $totalReadsFinalPE[1];
		$tmpPer = ($seqCountPE)?(sprintf "%0.2f", $totalReadsFinalPE[2]/$seqCountPE*100):"0";
		$tmpPer2 = ($seqCountPE)?(sprintf "%0.2f", $totalReadsFinalPE[0]/$seqCountPE*100):"0";
		$tmpPer3 = ($seqCountPE)?(sprintf "%0.2f", $totalReadsFinalPE[1]/$seqCountPE*100):"0";
		printf STAT "\t%s\t%s\t(%s / %s)\n", "Percentage of HQ filtered reads", $tmpPer."%", $tmpPer2."%", $tmpPer3."%";
		printf STAT "\t%s\t%d\n", "Total number of HQ filtered reads (Unpaired)", $totalReadsFinalUP[0];
		$tmpPer = ($seqCountPE)?(sprintf "%0.2f", $totalReadsFinalUP[0]/$seqCountPE*100):"0";
		printf STAT "\t%s\t%s\n", "Percentage of HQ filtered reads (Unpaired)", $tmpPer."%";

		print STAT "\n";
		printf STAT "\t%s\t%s\n", "QC analysis of Unpaired (UPOri*) reads:", "Total";
		printf STAT "\t%s\t%d\n", "Total number of Unpaired reads", $seqCount;
		printf STAT "\t%s\t%d\n", "Total number of trimmed reads containing homopolymer", $trimCount;
		printf STAT "\t%s\t%d\n", "Total number of trashed reads (<$lowestValidLen bp in length after trimming)", $lt100;
		printf STAT "\t%s\t%d\n", "Total number of low quality reads (excluding <$lowestValidLen reads)", $lQCount;
		printf STAT "\t%s\t%d\n", "Total number of HQ reads", $hQCount;
		$tmpPer = sprintf "%0.2f", $hQCount/$seqCount*100;
		printf STAT "\t%s\t%s\n", "Percentage of HQ reads", $tmpPer."%";
		printf STAT "\t%s\t%d\n", "Total number of bases", $totalBasesUPOri;
		printf STAT "\t%s\t%d\n", "Total number of bases in HQ reads", $totalBasesAfterHQ;
		printf STAT "\t%s\t%d\n", "Total number of HQ bases in HQ reads", $totalHQBasesAfterHQ;
		$tmpPer = sprintf "%0.2f", $totalHQBasesAfterHQ/$totalBasesAfterHQ*100;
		printf STAT "\t%s\t%s\n", "Percentage of HQ bases in HQ reads", $tmpPer."%";
		if(defined($priAdaLib)) {
			printf STAT "\t%s\t%d\n", "Number of Primer/Adaptor trimmed reads", $totalValidReadsWithPriAda;
		}
		else {
			printf STAT "\t%s\t%s\n", "Number of Primer/Adaptor trimmed reads", "NA";
		}
		printf STAT "\t%s\t%d\n", "Total number of HQ filtered reads (Unpaired)", $totalReadsFinal;
		$tmpPer = sprintf "%0.2f", $totalReadsFinal/$seqCount*100;
		printf STAT "\t%s\t%s\n", "Percentage of HQ filtered reads (Unpaired)", $tmpPer."%";
		print STAT ("-"x100)."\n";
		$tmpPer = sprintf "%0.2f", $totalReadsFinalPE[2]/($seqCount+$seqCountPE)*100;
		$tmpPer2 = sprintf "%0.2f", $totalReadsFinalUP[0]/($seqCount+$seqCountPE)*100;
		$tmpPer3 = sprintf "%0.2f", $totalReadsFinal/($seqCount+$seqCountPE)*100;
		printf STAT "\t%s\t%d\t(%s)\n", "Number of HQ filtered reads (Paired)", $totalReadsFinalPE[2], $tmpPer."%";
		printf STAT "\t%s\t%d\t(%s)\n", "Number of HQ filtered reads (UPPair*)", $totalReadsFinalUP[0], $tmpPer2."%";
		printf STAT "\t%s\t%d\t(%s)\n", "Number of HQ filtered reads (UPOri*)", $totalReadsFinal, $tmpPer3."%";
		print STAT ("-"x100)."\n";
		$tmpPer = sprintf "%0.2f", ($totalReadsFinalPE[2]+$totalReadsFinalUP[0]+$totalReadsFinal)/($seqCount+$seqCountPE)*100;
		printf STAT "\t%s\t%d\t(%s)\n", "Total number of HQ filtered reads", ($totalReadsFinalPE[2]+$totalReadsFinalUP[0]+$totalReadsFinal), $tmpPer."%";
		print STAT "\n";
		print STAT "* UPOri: Unpaired reads (i.e. reads without linker sequence) found in the input file\n";
		print STAT "  UPPair: One of the paired reads which passed QC\n";

		print STAT "\n\n";
	
######### Adding the statistics of PE, UPPair and UPOri reads
		$seqCount = $seqCount + $seqCountPE;
		$totalReadsFinal = $totalReadsFinalPE[2]+$totalReadsFinalUP[0]+$totalReadsFinal;
######### Done adding
		print STAT "Detailed QC statistics\n";
		my @arr = (
			["File name", $seqFileName, (fileparse($outSeqFile))[0]],
			["Total number of reads", $seqCount, $totalReadsFinal],
			["Minimum read length", $minRawLen, $minHQLen],
			["Maximum read length", $maxRawLen, $maxHQLen],
			["Average read length", (sprintf "%0.2f", $totalBases/$seqCount), (sprintf "%0.2f", $totalBasesFinal/$totalReadsFinal)],
			["Median read length", calcMedian(@rawLen), calcMedian(@hQLen)],
			["N25 length", calcN50(\@rawLen, 25), calcN50(\@hQLen, 25)],
			["N50 length", calcN50(\@rawLen, 50), calcN50(\@hQLen, 50)],
			["N75 length", calcN50(\@rawLen, 75), calcN50(\@hQLen, 75)],
			["N90 length", calcN50(\@rawLen, 90), calcN50(\@hQLen, 90)],
			["N95 length", calcN50(\@rawLen, 95), calcN50(\@hQLen, 95)],
			["Total number of bases", $totalBases, $totalBasesFinal],
			["Total number of HQ bases", $totalHQBases, $totalHQBasesFinal],
			["Percentage of HQ bases", (sprintf "%0.2f", $totalHQBases/$totalBases*100)."%", (sprintf "%0.2f", $totalHQBasesFinal/$totalBasesFinal*100)."%"],
			["Average quality score (Overall)", (sprintf "%0.2f", $avgQual/$seqCount), (sprintf "%0.2f", $avgQualFinal/$totalReadsFinal)],
		);
		for(my $i=0; $i<@arr; $i++) {
			if(!defined($isOnlyStat)) {
				printf STAT "\t%s\t%s\t%s\n", $arr[$i][0], $arr[$i][1], $arr[$i][2];
			}
			else {
				printf STAT "\t%s\t%s\n", $arr[$i][0], $arr[$i][1];				
			}
		}
	
		print STAT "\n\n";		
######### Subtracting the values back
		$seqCount = $seqCount - $seqCountPE;
		$totalReadsFinal = $totalReadsFinal-$totalReadsFinalPE[2]-$totalReadsFinalUP[0];
######### Done
	}

	my $lenDistF1 = getFileName($seqFile)."_lenDistribution.png";
	my $qualDistF1 = getFileName($seqFile)."_qualDistribution.png";
	my $sumPieFPE = getFileName($seqFile). "_PE_summary.png";
	my $sumPieF = getFileName($seqFile). "_UP_summary.png";
	my $gcDistF1 = getFileName($seqFile)."_gcDistribution.png";
	my $baseCntF1 = getFileName($seqFile)."_baseCompostion.png";

	my $c = 0;
	my @lenLabels = ();
	foreach my $arrRef (@lenDistrib) {
		my $str = "";
		foreach my $val (@{$arrRef}) {
			if($c == 0) {
				$str = "0-$lenInterval";
			}
			else {
				$str = $lenInterval*$c . "-" . $lenInterval*($c+1);
			}
			$c++;
			push(@lenLabels, $str);
		}
		last;
	}
	
	unshift(@lenDistrib, \@lenLabels);

	if($isGDMod) {
		drawLenDist(\@lenDistrib, $outFolder.$lenDistF1, getFileName($seqFile), 550, 350);
	}

	$c = 0;
	my @qualLabels = ();
	foreach my $arrRef (@qualDistrib) {
		my $str = "";
		foreach my $val (@{$arrRef}) {
			if($c == 0) {
				$str = "0";
				$str .= "-$qualInterval" if($qualInterval>1);
			}
			else {
				$str = $qualInterval*$c;
				$str .=  "-" . $qualInterval*($c) if($qualInterval>1);
			}
			push(@qualLabels, $str);
			$c++;
		}
		last;
	}
	
	unshift(@qualDistrib, \@qualLabels);

	if($isGDMod) {
		drawQualDist(\@qualDistrib, $outFolder.$qualDistF1, getFileName($seqFile), 650, 350);
	}

	my $trashedReads = $lt100PE[2];
	my $trimmedHP = $trimCountPE[2];
	my $trimmedPA = $totalValidReadsWithPriAdaPE[2];
	my $hQreadsExcptHP_PATrimmed = $totalReadsFinalPE[2] - $trimmedHP - $trimmedPA;
	my $UPPairCount = $totalReadsFinalUP[0];
	my $lQreadsGT100 = $seqCountPE - $totalReadsFinalPE[2] - $trashedReads - $UPPairCount;
	my @summaryData = (["", "", "", "", "", ""], [$trashedReads, $trimmedHP, $trimmedPA, $hQreadsExcptHP_PATrimmed, $lQreadsGT100, $UPPairCount]);
	
	if($isGDMod) {
		drawSummaryPiePE(\@summaryData, $outFolder.$sumPieFPE, 520, 350);
	}

	$trashedReads = $lt100;
	$trimmedHP = $trimCount;
	$trimmedPA = $totalValidReadsWithPriAda;
	$hQreadsExcptHP_PATrimmed = $totalReadsFinal - $trimmedHP - $trimmedPA;
	$lQreadsGT100 = $seqCount - $totalReadsFinal - $trashedReads;
	@summaryData = (["", "", "", "", ""], [$trashedReads, $trimmedHP, $trimmedPA, $hQreadsExcptHP_PATrimmed, $lQreadsGT100]);

	if($isGDMod) {
		drawSummaryPie(\@summaryData, $outFolder.$sumPieF, 520, 350);
	}

	$c=0;
	my @gcLabel;
	foreach my $ref (@gcDistrib) {
		foreach my $val (@{$ref}) {
			my $str = "";
			if($c == 0) {
				$str = "0-$gcInterval";
			}
			else {
				$str = $gcInterval*$c . "-" . $gcInterval*($c+1);
			}
			$c++;
			push(@gcLabel, $str);
		}
		last;
	}

	unshift(@gcDistrib, \@gcLabel);
	if($isGDMod) {
		drawGCDist(\@gcDistrib, $outFolder.$gcDistF1, getFileName($seqFile), 550, 350);
	}


	my @file1 = (["A", "T", "G", "C", "Non-ATGC"], $charCount[0]);
	@file1 = (["A", "T", "G", "C", "Non-ATGC"], $charCount[0], $charCount[1]) if(!$isOnlyStat);
	if($isGDMod) {
		drawBaseComp(\@file1, $outFolder.$baseCntF1, getFileName($seqFile), 500, 300);
	}
	

	close(I);
	close(Q);
	close(OI);
	close(OQ);
	close(STAT);

	my $iFol = getFilePath(abs_path($seqFile));
	my $oFol = abs_path($outFolder) . "/";
	my $inpFs = getFileName($seqFile);
	$inpFs .= ":::::" . getFileName($qualFile);
	my $htF = $oFol . "output_" . getFileName($seqFile);
	$htF .= ".html";
	my @fileNames4HTML;
	@fileNames4HTML = ($outSeqFile, $outQualFile, $lenDistF1, $baseCntF1, $gcDistF1, $qualDistF1, $sumPieFPE, $sumPieF);
	htmlPrint(getFilePath(abs_path($0)), getFileName($0), $htF, $iFol, $isOnlyStat, $inpFs, $statFile, $oFol, \@fileNames4HTML);

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
	my ($file, $rOrw, $ref) = @_;
	if($file =~ /\.gz$/i) {
		$$ref = new IO::Zlib;
		$$ref->open("$file", "rb") or die "Can not open file $file" if($rOrw eq "r");
		$$ref->open("$file", "wb") or die "Can not create file $file" if($rOrw eq "w");
	}
	else {
		open($$ref, "<$file") or die "Can not open file $file" if($rOrw eq "r");
		open($$ref, ">$file") or die "Can not create file $file" if($rOrw eq "w");
	}
}

sub splitPEReads {
	my ($seq, $qual, $id) = @_;
	my $len = length $seq;
	$linker = $linkerF;
	if($linkerR) {
		$linker = $linkerR if($seq =~ /$linkerR/);
	}
	my @seqs = split(/$linker/, $seq);
	my ($fRead, $rRead, $fQual, $rQual, $lQual);
	my $readStatus = 0;			# 0: Paired 1: Forward only 2: Reverse only
	if(@seqs == 1) {
		$fRead = $seqs[0];
		$readStatus = 1;
	}
	elsif(@seqs == 2) {
		if(!$seqs[1]) {
			$fRead = $seqs[0];
			$readStatus = 1;
		}
		elsif(!$seqs[0]) {
			$fRead = $seqs[1];
			$readStatus = 2;
		}
		else {
			$fRead = $seqs[0];
			$rRead = $seqs[1];
			$readStatus = 0;
		}
	}
	elsif(@seqs == 3) {
		if(!$seqs[0] && !$seqs[1]) {
			$fRead = $seqs[2];
			$readStatus = 2;
		}
		else {
			$fRead = $seqs[0];
			$rRead = $seqs[2];
			$readStatus = 0;
		}
	}
	else {
		print "There is some problem with the linker sequence found in the read: $id\n";
		return;
	}
	if($readStatus == 0) {
		my $fLen = length $fRead;
		my $rLen = length $rRead;
		my $lLen = $len - $fLen - $rLen;
		my ($t1, $t2);
		($fQual, $t1, $lQual, $t2, $rQual) = $qual =~ /^((\d{1,2}\s+){$fLen})((\d{1,2}\s+){$lLen})((\d{1,2}\s+){$rLen})$/;
	}
	elsif($readStatus == 1) {
		my $fLen = length $fRead;
		my $lLen = $len - $fLen;
		if($lLen == 0) {
			$fQual = $qual;
		}
		else {
			($fQual, $lQual) = $qual =~ /^((\d{1,2}\s+){$fLen})((\d{1,2}\s+){$lLen})$/;
		}
	}
	elsif($readStatus == 2) {
		my $fLen = length $fRead;
		my $lLen = $len - $fLen;
		($lQual, $fQual) = $qual =~ /^((\d{1,2}\s+){$lLen})((\d{1,2}\s+){$fLen})$/;
	}
	$fQual =~ s/\s+$// if($fQual);
	$rQual =~ s/\s+$// if($rQual);
	$lQual =~ s/\s+$// if($lQual);
	my $linkerCount = $seq =~ /$linker/gi;
	return ($readStatus, $fRead, $fQual, $linkerCount, $lQual, $rRead, $rQual);
}

sub doQC {
	my ($read, $qual) = @_;
	my $len = length $read;
	my ($basesAfterHQ);
	my $qcStatus = 1;   # 0: Failed QC (Read trashed), 1: Passed QC
	my $lenFltrStatus = 1;	# 1: Passed length filter
	my $homoStatus = 0;	# 0: absent homopolymer, 1: trimmed homopolymer
	my $primStatus = 0;	# 0: absent contam, 1: trimmed contam
	my $hqStatus = -1;	# 0: low quality read, 1: high quality read, -1: default (that means QC could not reach upto HQ filtering due to previous trashing in length filter)
	my $validBases = 0;
	if($len < $lowestValidLen && $isLenFilterOn) {
		$lenFltrStatus = 0;
		$qcStatus = 0;
	}
	else {
		if($homoPolyLen != 0) {
			if(hasPolyChar(\$read)) {
				$homoStatus = 1;
				$len = length $read;
				if($len >= $lowestValidLen || !$isLenFilterOn) {
					$qual = trimQualSeq($qual, $len, -1);
				}
			}
		}
		if($len < $lowestValidLen && $isLenFilterOn) {
			$lenFltrStatus = 0;
			$qcStatus = 0;
		}
		else {
			$validBases = isReadOfHQ($qual);
			if($validBases) {
				$hqStatus = 1;
				$basesAfterHQ = $len;
				if(defined $priAdaLib) {
					my $t=isWOPriAda(\$read);
					$len = length $read;
					if($t > -1) {
						$qual = trimQualSeq($qual, $len, $t);
						$primStatus = 1;
					}
					if(length $fastaSeq < $lowestValidLen && $isLenFilterOn) {
						$lenFltrStatus = 0;
						$qcStatus = 0;
					}
					else {
					}
				}
				else {
				}
			}
			else {
				$hqStatus = 0;
				$qcStatus = 0;
			}
		}
	}
	return ($read, $qual, $qcStatus, $lenFltrStatus, $primStatus, $homoStatus, $hqStatus, $basesAfterHQ, $validBases);
}

sub processSeq {
	$fastaSeq =~ s/\s//g;
	my ($readStatus, $fRead, $fQual, $linkerCount, $lQual, $rRead, $rQual) = splitPEReads($fastaSeq, $qualSeq, $prevFastaSeqId);
	my $len = length $fastaSeq;
##########Calculating statistics for Input Data
	$maxRawLen = max($maxRawLen, $len);
	$minRawLen = min($minRawLen, $len);
	push(@rawLen, $len);
	$qualSeq =~ s/\s+$//;				# To remove the last space added in 'else' part;
	my @tmpArr = getQualBases($qualSeq);
	$totalBases += $tmpArr[0];
	$totalHQBases += $tmpArr[1];
	$avgQual += $tmpArr[2];
	$lenDistrib[0][getIndex($len,$lenInterval)]++;
	$qualDistrib[0][getIndex($tmpArr[2],$qualInterval)]++;
	my $As = $fastaSeq =~ s/A/A/gi;
	my $Ts = $fastaSeq =~ s/T/T/gi;
	my $Gs = $fastaSeq =~ s/G/G/gi;
	my $Cs = $fastaSeq =~ s/C/C/gi;
	my $Ns = $len - $As - $Ts - $Gs - $Cs;
	my $gcPercent = ($Gs + $Cs)/$len*100;
	$gcDistrib[0][getIndex($gcPercent,$gcInterval)]++;
	$charCount[0][0] += $As;
	$charCount[0][1] += $Ts;
	$charCount[0][2] += $Gs;
	$charCount[0][3] += $Cs;
	$charCount[0][4] += $Ns;
#########Done calculating stat
	my ($readOut, $qualOut);
	my $outReadType = 0; # 0: Paired, 1: Unpaired
	if($readStatus == 0) {
		my @fReadQC = doQC($fRead, $fQual);			#($read, $qual, $qcStatus, $lenFltrStatus, $homoStatus, $hqStatus)
		my @rReadQC = doQC($rRead, $rQual);			#($read, $qual, $qcStatus, $lenFltrStatus, $homoStatus, $hqStatus)
		my $fLen = length $fRead;
		my $rLen = length $rRead;
		$seqCountPE++;
		$lt100PE[0]++ if($fReadQC[3] == 0);
		$lt100PE[1]++ if($rReadQC[3] == 0);
		$lt100PE[2]++ if($fReadQC[3] == 0 && $rReadQC[3] == 0);
		$totalValidReadsWithPriAdaPE[0]++ if($fReadQC[4] == 1);
		$totalValidReadsWithPriAdaPE[1]++ if($rReadQC[4] == 1);
		$totalValidReadsWithPriAdaPE[2]++ if($fReadQC[4] == 1 && $rReadQC[4] == 1);
		$trimCountPE[0]++ if($fReadQC[5] == 1);
		$trimCountPE[1]++ if($rReadQC[5] == 1);
		$trimCountPE[2]++ if($fReadQC[5] == 1 && $rReadQC[5] == 1);
		$lQCountPE[0]++ if($fReadQC[6] == 0);
		$lQCountPE[1]++ if($rReadQC[6] == 0);
		$lQCountPE[2]++ if($fReadQC[6] == 0 && $rReadQC[6] == 0);
		$hQCountPE[0]++ if($fReadQC[6] == 1);
		$hQCountPE[1]++ if($rReadQC[6] == 1);
		$hQCountPE[2]++ if($fReadQC[6] == 1 && $rReadQC[6] == 1);
		$totalBasesPE[0] += $fLen;
		$totalBasesPE[1] += $rLen;
		$totalBasesPE[2] += $fLen + $rLen;
		$totalBasesAfterHQPE[0] += $fReadQC[7] if($fReadQC[6] == 1);
		$totalBasesAfterHQPE[1] += $rReadQC[7] if($rReadQC[6] == 1);
		$totalBasesAfterHQPE[2] += $fReadQC[7]+$rReadQC[7] if($fReadQC[6] == 1 && $rReadQC[6] == 1);
		$totalHQBasesAfterHQPE[0] += $fReadQC[8] if($fReadQC[6] == 1);
		$totalHQBasesAfterHQPE[1] += $rReadQC[8] if($rReadQC[6] == 1);
		$totalHQBasesAfterHQPE[2] += $fReadQC[8]+$rReadQC[8] if($fReadQC[6] == 1 && $rReadQC[6] == 1);
		$totalReadsFinalPE[0]++ if($fReadQC[2] == 1);
		$totalReadsFinalPE[1]++ if($rReadQC[2] == 1);
		$totalReadsFinalPE[2]++ if($fReadQC[2] == 1 && $rReadQC[2] == 1);
		if($fReadQC[2] == 1 && $rReadQC[2] == 0) {
			$totalReadsFinalUP[0]++;
			$readOut = $fReadQC[0];
			$qualOut = $fReadQC[1];
			$outReadType = 1;
		}
		elsif($fReadQC[2] == 0 && $rReadQC[2] == 1) {
			$totalReadsFinalUP[0]++;
			$readOut = $rReadQC[0];
			$qualOut = $rReadQC[1];
			$outReadType = 1;
		}
		elsif($fReadQC[2] && $rReadQC[2]) {
			$readOut = $fReadQC[0] . ($linker x $linkerCount) . $rReadQC[0];
			$qualOut = $fReadQC[1] . " " . $lQual . " " . $rReadQC[1];
		}
	}
	else {
		my @fReadQC = doQC($fRead, $fQual);			#($read, $qual, $qcStatus, $lenFltrStatus, $homoStatus, $hqStatus)
		my $fLen = length $fRead;
		$seqCount++;
		$lt100++ if($fReadQC[3] == 0);
		$totalValidReadsWithPriAda++ if($fReadQC[4] == 1);
		$trimCount++ if($fReadQC[5] == 1);
		$lQCount++ if($fReadQC[6] == 0);
		$hQCount++ if($fReadQC[6] == 1);
		$totalBasesUPOri += $fLen;
		$totalBasesAfterHQ += $fReadQC[7] if($fReadQC[6] == 1);
		$totalHQBasesAfterHQ += $fReadQC[8] if($fReadQC[6] == 1);
		$totalReadsFinal++ if($fReadQC[2] == 1);
		if($fReadQC[2]) {
			$readOut = $fReadQC[0];
			$qualOut = $fReadQC[1];
			$outReadType = 1;
		}
	}
	if(defined($readOut)) {
##########Calculating statistics for Output Data
		$len = length $readOut;
		$maxHQLen = max($maxHQLen, $len);
		$minHQLen = min($minHQLen, $len);
		push(@hQLen, $len);
		$qualOut =~ s/\s+$//;				# To remove the last space, if present;
		my @tmpArr = getQualBases($qualOut);
		$totalBasesFinal += $tmpArr[0];
		$totalHQBasesFinal += $tmpArr[1];
		$avgQualFinal += $tmpArr[2];
		if(!defined($isOnlyStat)) {
			$lenDistrib[1][getIndex($len,$lenInterval)]++;
			$qualDistrib[1][getIndex($tmpArr[2],$qualInterval)]++;
			my $As = $readOut =~ s/A/A/gi;
			my $Ts = $readOut =~ s/T/T/gi;
			my $Gs = $readOut =~ s/G/G/gi;
			my $Cs = $readOut =~ s/C/C/gi;
			my $Ns = $len - $As - $Ts - $Gs - $Cs;
			my $gcPercent = ($len)?(($Gs + $Cs)/$len*100):0;
			$gcDistrib[1][getIndex($gcPercent,$gcInterval)]++;
			$charCount[1][0] += $As;
			$charCount[1][1] += $Ts;
			$charCount[1][2] += $Gs;
			$charCount[1][3] += $Cs;
			$charCount[1][4] += $Ns;
			print OI "$prevFastaSeqId\n";
			print OI formatSeq($readOut), "\n";
			print OQ "$prevFastaSeqId\n";
			print OQ formatQualSeq($qualOut), "\n";
		}
#########Done calculating stat
	}
	if(($seqCount+$seqCountPE) % (10000) == 0) {
		my $tmpP = sprintf "%0.0f", (($seqCount+$seqCountPE)/$ttlSeqCount*100);
		print "$indOfAnalysis: Number of reads processed: " . ($seqCount+$seqCountPE) . "/$ttlSeqCount ($tmpP\%)...\n";
	}
	return;
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

sub calcN50 {
	my @x = @{$_[0]};
	my $n = $_[1];
	@x=sort{$b<=>$a} @x;
	my $total = sum(@x);
	my ($count, $n50)=(0,0);
	for (my $j=0; $j<@x; $j++){
        $count+=$x[$j];
        if(($count>=$total*$n/100)){
            $n50=$x[$j];
            last;
        }
	}
	return $n50;
}

sub calcMedian {
	my @arr = @_;
	my @sArr = sort{$a<=>$b} @arr;
	my $arrLen = @arr;
	my $median;
	if($arrLen % 2 == 0) {
		$median = ($sArr[$arrLen/2-1] + $sArr[$arrLen/2])/2;
	}
	else {
		$median = $sArr[$arrLen/2];
	}
	return $median;
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
	my $priAdaStart = $_[2];
	my $trimmedQualSeq;
	if($priAdaStart != -1) {
		if($priAdaStart < 50) {
			my $t = $seqLen-1;
			$qualSeq =~ /((\d{1,2}\s+){$t}\d{1,2})$/;
			$trimmedQualSeq = $1;
		}
		else {
			$qualSeq =~ /^((\d{1,2}\s+){$seqLen})/;
			$trimmedQualSeq = $1;
		}
	}
	else {
		$qualSeq =~ /^((\d{1,2}\s+){$seqLen})/;
		$trimmedQualSeq = $1;
	}
	$trimmedQualSeq =~ s/\s+$//;
	return $trimmedQualSeq;
}

sub isReadOfHQ {	# Criteria for HQ is greater than or equal to 70% of bases have phred score > 20
	my $read = $_[0];
	my $validBaseCount = 0;
	my @ASCII = split(/\s+/, $read);
	my $readLen = scalar @ASCII;
	my $cutOffLen = sprintf("%0.0f", $readLen * $cutOffReadLen4HQ / 100);	# 70% length of read length is calculated.
	foreach my $val (@ASCII) {
		if($val >= $cutOffPhScore) {
			$validBaseCount++;
		}
	}
	if($validBaseCount >= $cutOffLen) {
		return $validBaseCount;				# Return true.
	}
	else {
		return 0;				# Return false.
	}
}

sub getQualBases {			# This will return an array. 1) Total bases 2) HQ bases 3) Average quality
	my $read = $_[0];
	my $qualSum = 0;
	my @retArr = ();
	my $validBaseCount = 0;
	my @ASCII = split(/\s+/, $read);
	my $readLen = scalar @ASCII;
	foreach my $val (@ASCII) {
		$qualSum += $val;
		if($val >= $cutOffPhScore) {
			$validBaseCount++;
		}
	}
	$retArr[0] = $readLen;
	$retArr[1] = $validBaseCount;
	$retArr[2] = ($readLen)?(sprintf "%0.2f", $qualSum/$readLen):0;
	return @retArr;
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

sub prtHelp {
	print "\n$0 options:\n\n";
	print "### Input reads (FASTA format; .fna and .qual files) (Required)\n";
	print "  -i <Read file> <Quality file> <Primer/Adaptor library>\n";
	print "    Read and quality file in FASTA format with primer/adaptor library\n";
	print "    User may choose from the provided primer/adaptor library or can give a file containing primer/adaptor sequences, one per line\n";
	print "    Multiple libraries can be given using multiple '-i' options\n";
	print "      For eg.: -i read1.fna read1.qual 3 -i read2.fna read2.qual 2\n\n";
	print "    Primer/Adaptor libraries:\n";
	my $c = 1;
	foreach my $lib (@priAdaLibNames) {
		print "      $c = $lib\n";
		$c++;
	}
	print "      N = Do not filter for Primer/Adaptor\n";
	print "      <File> = File for user defined primer/adaptor sequences, one per line\n";
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
	print "  -n | -homoPolyLen <Integer>\n";
	print "    Minimum length of the homopolymer to be trimmed (0: to skip the homopolymer trimming)\n";
	print "      For eg.: -n 8, will trim the right end of read from the homopolymer of at least 8 bases long\n";
	print "    default: 0 (homopolymer trimming is off)\n";
	print "  -m | -minLen <Integer>\n";
	print "    Filter sequences shorter than the given minimum length\n";
	print "    default: 100\n";
	print "  -f | -lenFilter <Y/N>\n";
	print "    Are sequences to be filtered on the basis of length: (Y)es or (N)o\n";
	print "    default: Y\n";
	print "  -linker <Linker Sequence>\n";
	print "    Linker sequence used while preparing the paired-end library for sequencing using Roche 454\n";
	print "    default: GTTGGAACCGAAAGGGTTTGAATTCAAACCCTTTCGGTTCCAAC\n";
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
	print "    default: By default, output folder (454QC_Filtered_files) will be generated where the input files are\n";
	print "  -z | -outputDataCompression <Character>\n";
	print "    Output format for HQ filtered data\n";
	print "    Formats:\n";
	print "      t = text FASTA files\n";
	print "      g = gzip compressed files\n";
	print "    default: t\n";
	print "\n";
}

sub prtUsage {
	print "\nUsage: perl $0 <options>\n";
	prtHelp();
}

sub getFileName {	# This sub takes a path of a file and returns just its name after separating the path from it.
	my $path = $_[0];
	my $name = "";
	$path =~ /([^\/]+)$/;
	$name = $1;
	return $name;	
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





	open(IMG, ">$fileNameWPath") or print STDERR "Error:\n\tCan not create image file: $fileNameWPath\n";
	binmode IMG;
	print IMG $myImage->png;
	close(IMG);
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
		x_labels_vertical	=> 1,
		legend_placement	=> 'BR',
		x_labels_vertical	=> 1,
		x_label_position	=> 1/2,
		long_ticks			=> 1,
		fgclr               => '#dddddd',
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
	open(IMG, ">$fileNameWPath") or print STDERR "Error:\n\tCan not create image file: $fileNameWPath\n";
	binmode IMG;
	print IMG $myImage->png;
	close(IMG);
}

sub drawSummaryPiePE {
	my $dataRef = $_[0];
	my $fileNameWPath = $_[1];
	my $width = $_[2];
	my $height = $_[3];
	my $mygraph = new GD::Graph::pie($width, $height+15);
	
	$mygraph->set( 
		title => "Summary of quality check and filtering of Paired reads",
		axislabelclr => 'black',
		pie_height => 40,
	
		l_margin => 15,
		r_margin => 15,
		b_margin => 70,
		start_angle => -10,
		dclrs => [ qw(lred cyan lyellow lgreen purple dblue) ],
		transparent => 0,
	) or warn $mygraph->error;
	
    $mygraph->set_label_font($f, 8);
    $mygraph->set_value_font(['verdana', 'arial'],14);
    $mygraph->set_title_font($f, 11);

	my $myImage = $mygraph->plot($dataRef);
	
	my $black = $myImage->colorAllocate(0,0,0);			# To set the color for the next time printing on the image.
	my $lred = $myImage->colorAllocate(255,0,0);	
	my $lyellow = $myImage->colorAllocate(255,255,0);	
	my $lgreen = $myImage->colorAllocate(0,255,0);
	my $cyan = $myImage->colorAllocate(0,255,255);
	my $purple = $myImage->colorAllocate(191,0,191);
	my $dblue = $myImage->colorAllocate(0,0,127);
	
	my $sum = sum(@{$$dataRef[1]});
	
	my $wrapbox = GD::Text::Wrap->new( $myImage,
	        line_space  => 4,
	        text        => (sprintf "Trashed reads (shorter than $lowestValidLen bp) (%0.2f", @{$$dataRef[1]}[0]/$sum*100) . "\%)",
	        color		=> $black,
	);
	
	$wrapbox->set(align => 'left', width => 300);
	$wrapbox->set_font($f, 8);
	$wrapbox->draw(20,$height-45);

	$wrapbox = GD::Text::Wrap->new( $myImage,
	        line_space  => 4,
	        text        => (sprintf "Homopolymer trimmed reads (%0.2f", @{$$dataRef[1]}[1]/$sum*100) . "\%)",
	        color		=> $black,
	);
	
	$wrapbox->set(align => 'left', width => 300);
	$wrapbox->set_font($f, 8);
	$wrapbox->draw($width/2+40,$height-45);

	$wrapbox = GD::Text::Wrap->new( $myImage,
	        line_space  => 4,
	        text        => (sprintf "Trashed reads (low quality reads) (%0.2f", @{$$dataRef[1]}[4]/$sum*100) . "\%)",
	        color		=> $black,
	);
	
	$wrapbox->set(align => 'left', width => 300);
	$wrapbox->set_font($f, 8);
	$wrapbox->draw(20,$height-30);

	$wrapbox = GD::Text::Wrap->new( $myImage,
	        line_space  => 4,
	        text        => (sprintf "Primer/Adaptor trimmed reads (%0.2f", @{$$dataRef[1]}[2]/$sum*100) . "\%)",
	        color		=> $black,
	);
	
	$wrapbox->set(align => 'left', width => 300);
	$wrapbox->set_font($f, 8);
	$wrapbox->draw($width/2+40,$height-30);

	$wrapbox = GD::Text::Wrap->new( $myImage,
	        line_space  => 4,
	        text        => (sprintf "High quality reads other than homopolymer and primer/adaptor trimmed (%0.2f", @{$$dataRef[1]}[3]/$sum*100) . "\%)",
	        color		=> $black,
	);
	
	$wrapbox->set(align => 'left', width => 500);
	$wrapbox->set_font($f, 8);
	$wrapbox->draw(20,$height-15);

	$wrapbox = GD::Text::Wrap->new( $myImage,
	        line_space  => 4,
	        text        => (sprintf "Unpaired reads (one of the paired reads which passed QC) (%0.2f", @{$$dataRef[1]}[5]/$sum*100) . "\%)",
	        color		=> $black,
	);
	
	$wrapbox->set(align => 'left', width => 500);
	$wrapbox->set_font($f, 8);
	$wrapbox->draw(20,$height);

	my $startRectX1 = 10;
	my $startRectX2 = $width/2+30;
	my $startRectY = $height-45;
	$myImage->filledRectangle($startRectX1,$startRectY,$startRectX1+8,$startRectY+8,$lred);
	$myImage->rectangle($startRectX1,$startRectY,$startRectX1+8,$startRectY+8,$black);

	$myImage->filledRectangle($startRectX2,$startRectY,$startRectX2+8,$startRectY+8,$cyan);
	$myImage->rectangle($startRectX2,$startRectY,$startRectX2+8,$startRectY+8,$black);

	$startRectY += 15;
	$myImage->filledRectangle($startRectX1,$startRectY,$startRectX1+8,$startRectY+8,$purple);
	$myImage->rectangle($startRectX1,$startRectY,$startRectX1+8,$startRectY+8,$black);

	$myImage->filledRectangle($startRectX2,$startRectY,$startRectX2+8,$startRectY+8,$lyellow);
	$myImage->rectangle($startRectX2,$startRectY,$startRectX2+8,$startRectY+8,$black);

	$startRectY += 15;
	$myImage->filledRectangle($startRectX1,$startRectY,$startRectX1+8,$startRectY+8,$lgreen);
	$myImage->rectangle($startRectX1,$startRectY,$startRectX1+8,$startRectY+8,$black);

	$startRectY += 15;
	$myImage->filledRectangle($startRectX1,$startRectY,$startRectX1+8,$startRectY+8,$dblue);
	$myImage->rectangle($startRectX1,$startRectY,$startRectX1+8,$startRectY+8,$black);

	open(IMG, ">$fileNameWPath") or print STDERR "Error:\n\tCan not create image file: $fileNameWPath\n";
	binmode IMG;
	print IMG $myImage->png;
	close(IMG);
	
}

sub drawSummaryPie {
	my $dataRef = $_[0];
	my $fileNameWPath = $_[1];
	my $width = $_[2];
	my $height = $_[3];
	my $mygraph = new GD::Graph::pie($width, $height);
	
	$mygraph->set( 
		title => "Summary of quality check and filtering of Unpaired (UPOri) reads",
		axislabelclr => 'black',
		pie_height => 40,
	
		l_margin => 15,
		r_margin => 15,
		b_margin => 50,
		start_angle => -10,
		dclrs => [ qw(lred cyan lyellow lgreen purple) ],
		transparent => 0,
	) or warn $mygraph->error;
	
    $mygraph->set_label_font($f, 8);
    $mygraph->set_value_font(['verdana', 'arial'],14);
    $mygraph->set_title_font($f, 11);

	my $myImage = $mygraph->plot($dataRef);
	
	my $black = $myImage->colorAllocate(0,0,0);			# To set the color for the next time printing on the image.
	my $lred = $myImage->colorAllocate(255,0,0);	
	my $lyellow = $myImage->colorAllocate(255,255,0);	
	my $lgreen = $myImage->colorAllocate(0,255,0);
	my $cyan = $myImage->colorAllocate(0,255,255);
	my $purple = $myImage->colorAllocate(191,0,191);
	
	my $sum = sum(@{$$dataRef[1]});
	
	my $wrapbox = GD::Text::Wrap->new( $myImage,
	        line_space  => 4,
	        text        => (sprintf "Trashed reads (shorter than $lowestValidLen bp) (%0.2f", @{$$dataRef[1]}[0]/$sum*100) . "\%)",
	        color		=> $black,
	);
	
	$wrapbox->set(align => 'left', width => 300);
	$wrapbox->set_font($f, 8);
	$wrapbox->draw(20,$height-45);

	$wrapbox = GD::Text::Wrap->new( $myImage,
	        line_space  => 4,
	        text        => (sprintf "Homopolymer trimmed reads (%0.2f", @{$$dataRef[1]}[1]/$sum*100) . "\%)",
	        color		=> $black,
	);
	
	$wrapbox->set(align => 'left', width => 300);
	$wrapbox->set_font($f, 8);
	$wrapbox->draw($width/2+40,$height-45);

	$wrapbox = GD::Text::Wrap->new( $myImage,
	        line_space  => 4,
	        text        => (sprintf "Trashed reads (low quality reads) (%0.2f", @{$$dataRef[1]}[4]/$sum*100) . "\%)",
	        color		=> $black,
	);
	
	$wrapbox->set(align => 'left', width => 300);
	$wrapbox->set_font($f, 8);
	$wrapbox->draw(20,$height-30);

	$wrapbox = GD::Text::Wrap->new( $myImage,
	        line_space  => 4,
	        text        => (sprintf "Primer/Adaptor trimmed reads (%0.2f", @{$$dataRef[1]}[2]/$sum*100) . "\%)",
	        color		=> $black,
	);
	
	$wrapbox->set(align => 'left', width => 300);
	$wrapbox->set_font($f, 8);
	$wrapbox->draw($width/2+40,$height-30);

	$wrapbox = GD::Text::Wrap->new( $myImage,
	        line_space  => 4,
	        text        => (sprintf "High quality reads other than homopolymer and primer/adaptor trimmed (%0.2f", @{$$dataRef[1]}[3]/$sum*100) . "\%)",
	        color		=> $black,
	);
	
	$wrapbox->set(align => 'left', width => 500);
	$wrapbox->set_font($f, 8);
	$wrapbox->draw(20,$height-15);

	my $startRectX1 = 10;
	my $startRectX2 = $width/2+30;
	my $startRectY = $height-45;
	$myImage->filledRectangle($startRectX1,$startRectY,$startRectX1+8,$startRectY+8,$lred);
	$myImage->rectangle($startRectX1,$startRectY,$startRectX1+8,$startRectY+8,$black);

	$myImage->filledRectangle($startRectX2,$startRectY,$startRectX2+8,$startRectY+8,$cyan);
	$myImage->rectangle($startRectX2,$startRectY,$startRectX2+8,$startRectY+8,$black);

	$startRectY += 15;
	$myImage->filledRectangle($startRectX1,$startRectY,$startRectX1+8,$startRectY+8,$purple);
	$myImage->rectangle($startRectX1,$startRectY,$startRectX1+8,$startRectY+8,$black);

	$myImage->filledRectangle($startRectX2,$startRectY,$startRectX2+8,$startRectY+8,$lyellow);
	$myImage->rectangle($startRectX2,$startRectY,$startRectX2+8,$startRectY+8,$black);

	$startRectY += 15;
	$myImage->filledRectangle($startRectX1,$startRectY,$startRectX1+8,$startRectY+8,$lgreen);
	$myImage->rectangle($startRectX1,$startRectY,$startRectX1+8,$startRectY+8,$black);

	open(IMG, ">$fileNameWPath") or print STDERR "Error:\n\tCan not create image file: $fileNameWPath\n";
	binmode IMG;
	print IMG $myImage->png;
	close(IMG);
	
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
	open(IMG, ">$fileNameWPath") or print STDERR "Error:\n\tCan not create image file: $fileNameWPath\n";
	binmode IMG;
	print IMG $myImage->png;
	close(IMG);
}

sub drawLenDist {
	my $dataRef = $_[0];
	my $fileNameWPath = $_[1];
	my $fileName = $_[2];
	my $width = $_[3];
	my $height = $_[4];
	
	my $mygraph = GD::Graph::bars->new($width, $height);
	
	$mygraph->set( 
		x_label => 'Read length (bp)',
		y_label => 'Number of reads',
		title => "Length distribution for $fileName",
		y_min_value => 0,
		box_axis => 0,
		line_width => 3,
		transparent => 0,
		dclrs => [ qw(lred dgreen) ],
		legend_placement	=> 'BR',
		x_labels_vertical	=> 1,
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
    $mygraph->set_x_axis_font($font_spec, 9);
    $mygraph->set_title_font($f, 11);
    $mygraph->set_legend_font($f, 8);

	my $myImage = $mygraph->plot($dataRef);
	open(IMG, ">$fileNameWPath") or print STDERR "Error:\n\tCan not create image file: $fileNameWPath\n";
	binmode IMG;
	print IMG $myImage->png;
	close(IMG);
}


sub isWOPriAda {
	my $seq = $_[0];
	chomp($$seq);

	my @rapid = (
		"CCATCTCATCCCTGCGTGTC",
		"CCATCTCATCCCTGCGTGTCTCCGACTCAG",
		"CTGAGTCGGAGA",
		"CCTATCCCCTGTGTGCCTTG",
		"CCTATCCCCTGTGTGCCTTGGCAGTCTCAG",
		"CTGAGACTGCCA",
	);

	my @arrPE = (
		"GCCTCCCTCGCGCCATCAG",
		"CTGATGGCGCGAGGG",
		"GCCTTGCCAGCCCGCTCAG",
		"CTGAGCGGGCTGGCA",
		"GCCTCCCTCGCGCCA",
		"GCCTTGCCAGCCCGC",
		"CCATCTCATCCCTGCGTGTC",
		"CCTATCCCCTGTGTGCCTTG",
	);

	my @arrAmplicon = (
		"CGTATCGCCTCCCTCGCGCCATCAG",
		"CGTATCGCCTCCCTCGCGCCATCAG",
		"CCATCTCATCCCTGCGTGTC",
		"CCTATCCCCTGTGTGCCTTG",
	);
	
	my @arrsmRna = (
		"GCCTCCCTCGCGCCATCAGTATCGTAGGCACCTGAGA",
		"GCCTTGCCAGCCCGCTCAGTATTGATGGTGCCTACAG",
		"CCATCTCATCCCTGCGTGTC",
		"CCTATCCCCTGTGTGCCTTG",
	);
	
	my @priAdas = (\@rapid, \@arrPE, \@arrAmplicon, \@arrsmRna);
	my %checkedPriStr = ();	# The 20 bp from start and end are stored in this hash as key. So that next time when another pri/ada seq

	my @priAdaSeqs = ();
	if($priAdaLib eq "u") {
		@priAdaSeqs = @usrDefinedPriAda;
	}
	else {
		@priAdaSeqs = @{$priAdas[$priAdaLib]};
	}

	my $priInd = 0;
	my $priAdaStart = 1;
	
	my $isMatched = 0;
	foreach my $priAda (@priAdaSeqs) {
		$priAdaStart = findSeq($priAda, $$seq, \%checkedPriStr);
		if($priAdaStart) {
			if($priAdaStart < 50) {
				$$seq = substr($$seq, $priAdaStart+$substrlen, length($$seq)-($priAdaStart+$substrlen));	
			}
			else {
				$$seq = substr($$seq, 0, $priAdaStart);
			}
			$isMatched = 1;
			last;
		}
	}
	
	if($isMatched) {
		return $priAdaStart;
	}
	else {
		return -1;
	}
}

sub findSeq {
	my $pri = $_[0];
	my $seq = $_[1];
	my $hashRef = $_[2];
	my $subsl = $substrlen;
	$subsl = length $pri if(length($pri) < $substrlen);
	my $spri = substr($pri, 0, $subsl);
	my $epri = substr($pri, (length $pri) - $subsl, $subsl);
	my $sseq = substr($seq, 0, 50);
	my $tmpInd = (length $seq) - 50;
	$tmpInd = 0 if($tmpInd < 0);
	my $eseq = substr($seq, $tmpInd, 50);
	my $ans;
	if(!defined($$hashRef{$spri})) {
		my @catches = String::Approx::amatch($spri, ['I0 D0 S1'], $sseq);
		if(@catches != 0) {
			return findStart($sseq, $spri);
		}
		@catches = String::Approx::amatch($spri, ['I0 D0 S1'], $eseq);
		if(@catches != 0) {
			return findStart($eseq, $spri) + length($seq) - 50;
		}
		$$hashRef{$spri} = 1;
	}
	if(!defined($$hashRef{$epri})) {
		my @catches = String::Approx::amatch($epri, ['I0 D0 S1'], $sseq);
		if(@catches != 0) {
			return findStart($sseq, $epri);
		}
		@catches = String::Approx::amatch($epri, ['I0 D0 S1'], $eseq);
		if(@catches != 0) {
			return findStart($eseq, $epri) + length($seq) - 50;
		}
		$$hashRef{$epri} = 1;
	}
	return 0;
}

use re qw(eval);
use vars qw($matchStart);

sub findStart {
   my $pattern;
   local $_;
   ($_, $pattern) = @_;
   $pattern = fuzzy_pattern($pattern, $mismLim);
   my @results;
   local $matchStart;
   my $instrumentedPattern = qr/(?{ $matchStart = pos() })$pattern/;
   while (/$instrumentedPattern/g) {
      my $nextStart = pos();
      return $matchStart;
      push @results, "[$matchStart..$nextStart)";
      pos() = $matchStart+1;
   }
}

sub fuzzy_pattern {
   my ($original_pattern, $mismatches_allowed) = @_;
   $mismatches_allowed >= 0
      or die "Number of mismatches must be greater than or equal to zero\n";
   my $new_pattern = make_approximate($original_pattern, $mismatches_allowed);
   return qr/$new_pattern/;
}

sub make_approximate {
   my ($pattern, $mismatches_allowed) = @_;
   if ($mismatches_allowed == 0) { return $pattern }
   elsif (length($pattern) <= $mismatches_allowed)
      { $pattern =~ tr/ACTG/./; return $pattern }
   else {
      my ($first, $rest) = $pattern =~ /^(.)(.*)/;
      my $after_match = make_approximate($rest, $mismatches_allowed);
      if ($first =~ /[ACGT]/) {
         my $after_miss = make_approximate($rest, $mismatches_allowed-1);
         return "(?:$first$after_match|.$after_miss)";
      }
      else { return "$first$after_match" }
   }
}
