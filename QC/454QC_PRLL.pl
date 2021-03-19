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
require "454html.pl";
use threads('yield');
use File::Path;
use Thread::Queue;
my $DataQueue;
my $ProcessingQueue;
my $thr;


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
		print "Error:\n\tCan not find module 'String::Approx'\n";
		print "\tInstall it and try again\n\n";
		exit;
	}
}


# Setting parameters
my $lowestValidLen = 100;
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

GetOptions(
			"i=s{$noOfInp}" => \@files,
			"h|help" => \$helpAsked,
			"l|cutOffReadLen4HQ=f" => \$cutOffReadLen4HQ,
			"n|homoPolyLen=i" => \$homoPolyLen,
			"o|outputFolder=s" => \$outFolder,
			"z|outputDataCompression=s" => \$outputDataFmt,
			"t|statOutFmt=i" => \$statOutFmt,
			"onlyStat" => \$isOnlyStat,
			"c|cpus=i" => \$noOfProcesses,
			"s|cutOffQualScore=i" => \$cutOffPhScore,
			"m|minLen=i" => \$lowestValidLen,
			"f|lenFilter=s" => \$isLenFilterOn,
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

#my $pm = new Parallel::ForkManager($noOfProcesses);


my $seqCount = 0;
my $substrlen = 20;						# For removePriAda
my $mismLim = 1;						# For removePriAda

my $trimCount = 0;
my $lt100 = 0;
my $hQCount = 0;
my $lQCount = 0;
my $maxRawLen = 0;
my $minRawLen = 1000000000000;
#my $avgRawLen = 0;
my $maxHQLen = 0;
my $minHQLen = 1000000000000;
my $avgHQLen = 0;
my @rawLen = ();
my @hQLen = ();
my $totalBases = 0;
my $totalHQBases = 0;
my $totalBasesAfterHQ = 0;
my $totalHQBasesAfterHQ = 0;
my $totalBasesFinal = 0;
my $totalHQBasesFinal = 0;
my $totalReadsFinal = 0;
my $avgQual = 0;
my $avgQualFinal = 0;
my $totalValidReadsWithPriAda = 0;
my $totalValidReadsNoPriAda = 0;
my @lenDistrib = ();
my $lenInterval = 40;
my @qualDistrib = ();
my $qualInterval = 1;
my @gcDistrib = ();
my $gcInterval = 5;
my @charCount = ();


my $cmaxRawLen = 0;
my $cminRawLen = 1000000000000;
my @crawLen = ();
my $ctotalBases = 0;
my $ctotalHQBases = 0;
my $cavgQual = 0;
my $clt100 = 0;
my $ctrimCount = 0;
my $chQCount = 0;
my $ctotalBasesAfterHQ = 0;
my $cmaxHQLen = 0;
my $cminHQLen = 1000000000000;
my $cavgHQLen = 0;
my @chQLen = ();
my $ctotalReadsFinal = 0;
my $ctotalBasesFinal = 0;
my $ctotalHQBasesFinal = 0;
my $cavgQualFinal = 0;
my $clQCount = 0;
my $ctotalHQBasesAfterHQ = 0;
my $ctotalValidReadsWithPriAda = 0;
my $ctotalValidReadsNoPriAda = 0;
my @clenDistrib = ();
my @cqualDistrib = ();
my @cgcDistrib = ();
my @ccharCount = ();


my $fastaSeqId = "";
my $fastaSeq = "";
my $qualSeqId = "";
my $qualSeq = "";
my $prevFastaSeqId = "";
my $indOfAnalysis = 0;
my $uniqFolder = "";
my $isInpGzip = 0;

my @idArr = ();
my @seqArr = ();
my @qualArr = ();

my $font_spec = getFilePath($0) . "lib/Fonts/Dustismo_Sans.ttf";
my $f = getFilePath($0) . "lib/Fonts/LucidaSansDemiBold.ttf";


#Temp
my $c=0;

foreach my $inpData (@files) {
	$indOfAnalysis++;
#my $pid = $pm->start and next;

	$fastaSeqId = "";
	$fastaSeq = "";
	$qualSeqId = "";
	$qualSeq = "";
	$seqCount = 0;

	$trimCount = 0;
	$lt100 = 0;
	$hQCount = 0;
	$lQCount = 0;
	$maxRawLen = 0;
	$minRawLen = 1000000000000;
	$maxHQLen = 0;
	$minHQLen = 1000000000000;
	$avgHQLen = 0;
	@rawLen = ();
	@hQLen = ();
	$totalBases = 0;
	$totalHQBases = 0;
	$totalBasesAfterHQ = 0;
	$totalHQBasesAfterHQ = 0;
	$totalBasesFinal = 0;
	$totalHQBasesFinal = 0;
	$totalReadsFinal = 0;
	$avgQual = 0;
	$avgQualFinal = 0;
	$totalValidReadsWithPriAda = 0;
	$totalValidReadsNoPriAda = 0;
	@lenDistrib = ();
	@qualDistrib = ();
	@gcDistrib = ();
	@charCount = ();

	@idArr = ();
	@seqArr = ();
	@qualArr = ();

	$inpData =~ s/\\([A-Za-z_\.])/\/$1/g;		# To remove '\' from the path of windows file
	my @iData = split(" ", $inpData);
	my $seqFile = $iData[0];
	my $qualFile = $iData[1];
	if($seqFile =~ /\.gz$/i || $qualFile =~ /\.gz$/i) {
		$isInpGzip = 1;
	}
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

    $DataQueue = Thread::Queue->new();
    unlink($outSeqFile) if(-e $outSeqFile);
    unlink($outQualFile) if(-e $outQualFile);
    $thr = threads->create(sub {
        while (my $DataElement = $DataQueue->dequeue()) {
        	$DataElement =~ s/([sq]$)//;
        	my $readType = $1;
			my $outH;
			openFileGetHandle($outSeqFile, "a", \$outH) if($readType eq "s");
			openFileGetHandle($outQualFile, "a", \$outH) if($readType eq "q");
			*OOO = $outH;
			print OOO "$DataElement";
			close(OOO);
        }
    });
    $ProcessingQueue = Thread::Queue->new();

	my $iH;
	openFileGetHandle($seqFile, "r", \$iH);
	*I = $iH;

		do {
			$uniqFolder = "";
			for(my $i=0; $i<5; $i++) {
				$uniqFolder .= int(rand(10));
			}
			$uniqFolder = $outFolder . $uniqFolder;
		}
		while(-e $uniqFolder);
		mkdir($uniqFolder) or die "Can not create folder for temporary files\n";
	open(STAT, ">$statFile") or die "Can not open file: $statFile\n";
	while(my $line = <I>) {
		$seqCount++ if($line =~ /^>/);
	}
	close(I);
	
	
	if($seqFile =~ /\.gz$/i || $qualFile =~ /\.gz$/i) {
		my @fileNames = ($seqFile, $qualFile);
		my $thRef = threads->create('readDivideGzip', @fileNames);
		threading4Processing();
		$thRef->join();
		print "$indOfAnalysis: Number of reads processed: " . $seqCount . "/$seqCount (100\%)...\n";
	}
	else {
		print "$indOfAnalysis: Number of reads processed: " . "0/$seqCount (0\%)...\n";
		undef $iH;
		openFileGetHandle($seqFile, "r", \$iH);
		*I = $iH;
		my $qH;
		openFileGetHandle($qualFile, "r", \$qH);
		*Q = $qH;
		
		my @thArr = ();
		my $noOfSeqPerThread = int(($seqCount-1)/$noOfProcesses); #20000;
		my $roughSeqCounter = 0;
		my $sCounter = 0;
		my $jobCounter = 0;
		my $ttlJobCounter = 0;
		my $fileEOF = 0;
	
		while(1) {
			$jobCounter = 0;
			OUTERLOOP:
			for(my $i=0; $i<$noOfProcesses; $i++) {
				$jobCounter++;
				$ttlJobCounter++;
				for(my $j=0; $j<$noOfSeqPerThread;) {
					my $line = <I>;
					chomp $line;
					my $qualLine = <Q>;
					chomp($qualLine);
					if($line =~ /^>/) {
						$sCounter++;
						$prevFastaSeqId = $fastaSeqId;
						$fastaSeqId = $line;
						$qualSeqId = $qualLine;
						if($fastaSeqId ne $qualSeqId) {
							print "Error: Read Id doesn't match in sequence and quality file for read number $seqCount in sequence file.\n";
							exit(-1);
						}
						if($fastaSeq ne "") {
							push(@idArr, $prevFastaSeqId);
							push(@seqArr, $fastaSeq);
							push(@qualArr, $qualSeq);
							$j++;
							if($j == $noOfSeqPerThread) {
								my $id = sprintf "%05s", $ttlJobCounter;
								my @refArr = (\@idArr, \@seqArr, \@qualArr, $id);
								$thArr[$i] = threads->create('passSeq', @refArr);
								@idArr = ();
								@seqArr = ();
								@qualArr = ();
							}
						}
						$fastaSeq = "";
						$qualSeq = "";
						if($sCounter == $seqCount) {
							$jobCounter-- if((scalar @idArr) != 0);
							$ttlJobCounter-- if((scalar @idArr) != 0);		#This is not commented because this is used to create part files. 
							$fileEOF = 1;
							last OUTERLOOP;
						}
					}
					else {
						$fastaSeq .= $line;
						$qualSeq .= $qualLine . " ";
					}
				}
			}
			my @tmpArr = ();
			for(my $i=0; $i<$jobCounter; $i++) {
		        my $refArr = $thArr[$i]->join;
		        &updateData(@{$refArr});
		        $roughSeqCounter += $noOfSeqPerThread;
		        if($roughSeqCounter%$noOfSeqPerThread == 0) {
					my $tmpP = sprintf "%0.0f", ($roughSeqCounter/$seqCount*100);
					print "$indOfAnalysis: Number of reads processed: " . $roughSeqCounter . "/$seqCount ($tmpP\%)...\n";
		        }
			}
			last if($fileEOF);
		}
		while(my $line = <I>) {
			my $qualLine = <Q>;
			chomp $line;
			chomp($qualLine);
			$fastaSeq .= $line;
			$qualSeq .= $qualLine . " ";
		}
			$ttlJobCounter++;
			my $id = sprintf "%05s", $ttlJobCounter;
			$prevFastaSeqId = $fastaSeqId;
			push(@idArr, $prevFastaSeqId);
			push(@seqArr, $fastaSeq);
			push(@qualArr, $qualSeq);
			my @refArr = (\@idArr, \@seqArr, \@qualArr, $id);
			my $thId = threads->create('passSeq', @refArr);
			my $refRefArr = $thId->join; #passSeq(@refArr);
			&updateData(@{$refRefArr});
		close(I);
		close(Q);
		print "$indOfAnalysis: Number of reads processed: " . $seqCount . "/$seqCount (100\%)...\n";
	}
	if(!defined($isOnlyStat)) {
		my ($sHndl, $qHndl);
			openFileGetHandle($outSeqFile, "w", \$sHndl);
			openFileGetHandle($outQualFile, "w", \$qHndl);
			*OOS = $sHndl;
			*OOQ = $qHndl;
		print "$indOfAnalysis: Printing filtered data...\n";
		opendir(DIR, $uniqFolder);
		my @partFiles = readdir(DIR);
		@partFiles = sort @partFiles;
		foreach my $pFile (@partFiles) {
			next if($pFile =~ /\./);
			my $npFile = "$uniqFolder/$pFile";
			open(P, "<$npFile") or die "Can not open part file\n";
			while(<P>) {
				print OOS if($pFile =~ /seq[^\n]+out/);
				print OOQ if($pFile =~ /qual[^\n]+out/);
			}
			close(P);
		}
		closedir(DIR);
		close(OOS);
		close(OOQ);
	}	
	print "$indOfAnalysis: Analysis completed\n";
	print "$indOfAnalysis: Printing Statistics...\n";

	if($statOutFmt == 1) {
		my $inde = " " x 1;
		my $tmpPer = 0;
		printf STAT "Parameters\n";
		printf STAT "$inde %-40s %s  %s\n", "Input files ", $seqFile, $qualFile;
		printf STAT "$inde %-40s %s\n", "Primer/Adaptor library", defined($priAdaLib)?(($priAdaLib ne "u")?$priAdaLibNames[$priAdaLib]:"User defined ($priAdaFile)"):"NA";
		printf STAT "$inde %-40s %s\n", "Homopolymer trimming", "Off" if($homoPolyLen == 0);
		printf STAT "$inde %-40s %s\n", "Homopolymer trimming", "On" if($homoPolyLen != 0);
		printf STAT "$inde %-40s %s\n", "Length of the homopolymer to be removed", $homoPolyLen if($homoPolyLen != 0);
		printf STAT "$inde %-40s %s\n", "Length filter", ($isLenFilterOn)?"On":"Off";
		printf STAT "$inde %-40s %s\n", "Cut-off for minimum read length", $lowestValidLen if($isLenFilterOn);
		printf STAT "$inde %-40s %s\n", "Cut-off read length for HQ", $cutOffReadLen4HQ."%";
		printf STAT "$inde %-40s %s\n", "Cut-off quality score", $cutOffPhScore;
		printf STAT "$inde %-40s %s\n", "Only statistics", defined($isOnlyStat)?"On":"Off";
		printf STAT "$inde %-40s %s\n", "Number of CPUs", $noOfProcesses;
	
		print STAT "\n\n";
	
		print STAT "QC statistics\n";
		printf STAT "$inde %-70s %s\n", "File name", $seqFileName;
		printf STAT "$inde %-70s %d\n", "Total number of reads", $seqCount;
		printf STAT "$inde %-70s %d\n", "Total number of trimmed reads containing homopolymer", $trimCount;
		printf STAT "$inde %-70s %d\n", "Total number of trashed reads (<$lowestValidLen bp in length after trimming)", $lt100;
		printf STAT "$inde %-70s %d\n", "Total number of low quality reads (excluding <$lowestValidLen reads)", $lQCount;
		printf STAT "$inde %-70s %d\n", "Total number of HQ reads", $hQCount;
		$tmpPer = sprintf "%0.2f", $hQCount/$seqCount*100;
		printf STAT "$inde %-70s %s\n", "Percentage of HQ reads", $tmpPer."%";
		printf STAT "$inde %-70s %d\n", "Total number of bases", $totalBases;
		printf STAT "$inde %-70s %d\n", "Total number of bases in HQ reads", $totalBasesAfterHQ;
		printf STAT "$inde %-70s %d\n", "Total number of HQ bases in HQ reads", $totalHQBasesAfterHQ;
		$tmpPer = sprintf "%0.2f", $totalHQBasesAfterHQ/$totalBasesAfterHQ*100;
		printf STAT "$inde %-70s %s\n", "Percentage of HQ bases in HQ reads", $tmpPer."%";
		if(defined($priAdaLib)) {
			printf STAT "$inde %-70s %d\n", "Number of Primer/Adaptor trimmed reads", $totalValidReadsWithPriAda;
		}
		else {
			printf STAT "$inde %-70s %s\n", "Number of Primer/Adaptor trimmed reads", "NA", "NA";
		}
		printf STAT "$inde %-70s %d\n", "Total number of HQ filtered reads", $totalReadsFinal;
		$tmpPer = sprintf "%0.2f", $totalReadsFinal/$seqCount*100;
		printf STAT "$inde %-70s %s\n", "Percentage of HQ filtered reads", $tmpPer."%";
	
		print STAT "\n\n";
	
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
	}
	elsif($statOutFmt == 2) {
		my $inde = " " x 1;
		my $tmpPer = 0;
		printf STAT "Parameters\n";
		printf STAT "\t%s\t%s\t%s\n", "Input files ", $seqFile, $qualFile;
		printf STAT "\t%s\t%s\n", "Primer/Adaptor library", defined($priAdaLib)?(($priAdaLib ne "u")?$priAdaLibNames[$priAdaLib]:"User defined ($priAdaFile)"):"NA";
		printf STAT "\t%s\t%s\n", "Homopolymer trimming", "Off" if($homoPolyLen == 0);
		printf STAT "\t%s\t%s\n", "Homopolymer trimming", "On" if($homoPolyLen != 0);
		printf STAT "\t%s\t%s\n", "Length of the homopolymer to be removed", $homoPolyLen if($homoPolyLen != 0);
		printf STAT "\t%s\t%s\n", "Length of the homopolymer to be removed", $homoPolyLen;
		printf STAT "\t%s\t%s\n", "Length filter", ($isLenFilterOn)?"On":"Off";
		printf STAT "\t%s\t%s\n", "Cut-off for minimum read length", $lowestValidLen if($isLenFilterOn);
		printf STAT "\t%s\t%s\n", "Cut-off read length for HQ", $cutOffReadLen4HQ."%";
		printf STAT "\t%s\t%s\n", "Cut-off quality score", $cutOffPhScore;
		printf STAT "\t%s\t%s\n", "Only statistics", defined($isOnlyStat)?"On":"Off";
		printf STAT "\t%s\t%s\n", "Number of CPUs", $noOfProcesses;
	
		print STAT "\n\n";
	
		print STAT "QC statistics\n";
		printf STAT "\t%s\t%s\n", "File name", $seqFileName;
		printf STAT "\t%s\t%d\n", "Total number of reads", $seqCount;
		printf STAT "\t%s\t%d\n", "Total number of trimmed reads containing homopolymer", $trimCount;
		printf STAT "\t%s\t%d\n", "Total number of trashed reads (<$lowestValidLen bp in length after trimming)", $lt100;
		printf STAT "\t%s\t%d\n", "Total number of low quality reads (excluding <$lowestValidLen reads)", $lQCount;
		printf STAT "\t%s\t%d\n", "Total number of HQ reads", $hQCount;
		$tmpPer = sprintf "%0.2f", $hQCount/$seqCount*100;
		printf STAT "\t%s\t%s\n", "Percentage of HQ reads", $tmpPer."%";
		printf STAT "\t%s\t%d\n", "Total number of bases", $totalBases;
		printf STAT "\t%s\t%d\n", "Total number of bases in HQ reads", $totalBasesAfterHQ;
		printf STAT "\t%s\t%d\n", "Total number of HQ bases in HQ reads", $totalHQBasesAfterHQ;
		$tmpPer = sprintf "%0.2f", $totalHQBasesAfterHQ/$totalBasesAfterHQ*100;
		printf STAT "\t%s\t%s\n", "Percentage of HQ bases in HQ reads", $tmpPer."%";
		if(defined($priAdaLib)) {
			printf STAT "\t%s\t%d\n", "Number of Primer/Adaptor trimmed reads", $totalValidReadsWithPriAda;
		}
		else {
			printf STAT "\t%s\t%s\n", "Number of Primer/Adaptor trimmed reads", "NA", "NA";
		}
		printf STAT "\t%s\t%d\n", "Total number of HQ filtered reads", $totalReadsFinal;
		$tmpPer = sprintf "%0.2f", $totalReadsFinal/$seqCount*100;
		printf STAT "\t%s\t%s\n", "Percentage of HQ filtered reads", $tmpPer."%";
	
		print STAT "\n\n";
	
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
	}

	my $lenDistF1 = getFileName($seqFile)."_lenDistribution.png";
	my $qualDistF1 = getFileName($seqFile)."_qualDistribution.png";
	my $sumPieF = getFileName($seqFile). "_summary.png";
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

	my $trashedReads = $lt100;
	my $trimmedHP = $trimCount;
	my $trimmedPA = $totalValidReadsWithPriAda;
	my $hQreadsExcptHP_PATrimmed = $totalReadsFinal - $trimmedHP - $trimmedPA;
	my $lQreadsGT100 = $seqCount - $totalReadsFinal - $trashedReads;
	my @summaryData = (["", "", "", "", ""], [$trashedReads, $trimmedHP, $trimmedPA, $hQreadsExcptHP_PATrimmed, $lQreadsGT100]);
	
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
		drawBaseComp(\@file1, $outFolder.getFileName($seqFile)."_baseCompostion.png", getFileName($seqFile), 500, 300);
	}


	close(STAT);

	my $iFol = getFilePath(abs_path($seqFile));
	my $oFol = abs_path($outFolder) . "/";
	my $inpFs = getFileName($seqFile);
	$inpFs .= ":::::" . getFileName($qualFile);
	my $htF = $oFol . "output_" . getFileName($seqFile);
	$htF .= ".html";
	my @fileNames4HTML;
	@fileNames4HTML = ($outSeqFile, $outQualFile, $lenDistF1, $baseCntF1, $gcDistF1, $qualDistF1, $sumPieF);
	htmlPrint(getFilePath(abs_path($0)), getFileName($0), $htF, $iFol, $isOnlyStat, $inpFs, $statFile, $oFol, \@fileNames4HTML);

    $DataQueue->enqueue(undef);
    $thr->join();
	rmtree($uniqFolder, 0, 0);
#$pm->finish;
}
#$pm->wait_all_children;

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
		$$ref->open("$file", "ab") or die "Can not open-append file $file" if($rOrw eq "a");
	}
	else {
		open($$ref, "<$file") or die "Can not open file $file" if($rOrw eq "r");
		open($$ref, ">$file") or die "Can not create file $file" if($rOrw eq "w");
		open($$ref, ">>$file") or die "Can not open-append file $file" if($rOrw eq "a");
	}
}


sub updateData() {
	my @arr = @_;
	$maxRawLen = max($maxRawLen, $arr[0]);
	$minRawLen = min($minRawLen, $arr[1]);
	push(@rawLen, @{$arr[2]});
	$totalBases += $arr[3];
	$totalHQBases += $arr[4];
	$avgQual += $arr[5];
	$lt100 += $arr[6];
	$trimCount += $arr[7];
	$hQCount += $arr[8];
	$totalBasesAfterHQ += $arr[9];
	$maxHQLen = max($maxHQLen, $arr[10]);
	$minHQLen = min($minHQLen, $arr[11]);
	$avgHQLen += $arr[12];
	push(@hQLen, @{$arr[13]});
	$totalReadsFinal += $arr[14];
	$totalBasesFinal += $arr[15];
	$totalHQBasesFinal += $arr[16];
	$avgQualFinal += $arr[17];
	$lQCount += $arr[18];
	$totalHQBasesAfterHQ += $arr[19];
	$totalValidReadsWithPriAda += $arr[20];
	$totalValidReadsNoPriAda += $arr[21];
	addTwoArrays($arr[22], \@lenDistrib);
	addTwoArrays($arr[23], \@qualDistrib);
	addTwoArrays($arr[24], \@gcDistrib);
	addTwoArrays($arr[25], \@charCount);
}

sub resetVariables() {
	$cmaxRawLen = 0;
	$cminRawLen = 1000000000000;
	@crawLen = ();
	$ctotalBases = 0;
	$ctotalHQBases = 0;
	$cavgQual = 0;
	$clt100 = 0;
	$ctrimCount = 0;
	$chQCount = 0;
	$ctotalBasesAfterHQ = 0;
	$cmaxHQLen = 0;
	$cminHQLen = 1000000000000;
	$cavgHQLen = 0;
	@chQLen = ();
	$ctotalReadsFinal = 0;
	$ctotalBasesFinal = 0;
	$ctotalHQBasesFinal = 0;
	$cavgQualFinal = 0;
	$clQCount = 0;
	$ctotalHQBasesAfterHQ = 0;
	$ctotalValidReadsWithPriAda = 0;
	$ctotalValidReadsNoPriAda = 0;	
	@clenDistrib = ();
	@cqualDistrib = ();
	@cgcDistrib = ();
	@ccharCount = ();
}

sub readDivideGzip {
	my ($seqFile, $qualFile) = @_;
	my $chunkCounter = 0;
	my $sCounter = 0;
	my $fileEOF = 0;
	my ($seqH, $qualH, $id);
	my $isFileOpen = 0;
	my $iH;
	openFileGetHandle($seqFile, "r", \$iH);
	my $qH;
	openFileGetHandle($qualFile, "r", \$qH);
	my $noOfSeqPerThread = int(($seqCount-1)/$noOfProcesses); #20000;
	while(1) {
		OUTERLOOP:
		for(my $i=0; $i<$noOfProcesses; $i++) {
			$chunkCounter++;
			$id = sprintf "%05s", $chunkCounter;
			undef $seqH;
			openFileGetHandle("$uniqFolder/part_seq_$id", "w", \$seqH);
			undef $qualH;
			openFileGetHandle("$uniqFolder/part_qual_$id", "w", \$qualH);
			$isFileOpen = 1;
			for(my $j=0; $j<$noOfSeqPerThread;) {
				my $line = <$iH>;
				chomp $line;
				my $qualLine = <$qH>;
				chomp($qualLine);
				if($line =~ /^>/) {
					$sCounter++;
					$prevFastaSeqId = $fastaSeqId;
					$fastaSeqId = $line;
					$qualSeqId = $qualLine;
					if($fastaSeqId ne $qualSeqId) {
						print "Error: Read Id doesn't match in sequence and quality file for read number $seqCount in sequence file.\n";
						exit(-1);
					}
					if($fastaSeq ne "") {
						print $seqH "$prevFastaSeqId\n";
						print $seqH "$fastaSeq\n";
						print $qualH "$prevFastaSeqId\n";
						print $qualH "$qualSeq\n";
						$j++;
						if($j == $noOfSeqPerThread) {
							close($seqH);
							close($qualH);
							$ProcessingQueue->enqueue("$uniqFolder/part_seq_$id"."\t"."$uniqFolder/part_qual_$id");
							$isFileOpen = 0;
						}
					}
					$fastaSeq = "";
					$qualSeq = "";
					if($sCounter == $seqCount) {
						$chunkCounter-- if((scalar @idArr) != 0);		#This is not commented because this is used to create part files. 
						$fileEOF = 1;
						last OUTERLOOP;
					}
				}
				else {
					$fastaSeq .= $line;
					$qualSeq .= $qualLine . " ";
				}
			}
		}
		last if($fileEOF);
	}
	while(my $line = <$iH>) {
		my $qualLine = <$qH>;
		chomp $line;
		chomp($qualLine);
		$fastaSeq .= $line;
		$qualSeq .= $qualLine . " ";
	}
		if(! $isFileOpen) {
			$chunkCounter++;
			$id = sprintf "%05s", $chunkCounter;
			undef $seqH;
			openFileGetHandle("$uniqFolder/part_seq_$id", "w", \$seqH);
			undef $qualH;
			openFileGetHandle("$uniqFolder/part_qual_$id", "w", \$qualH);
		}
		$prevFastaSeqId = $fastaSeqId;
		print $seqH "$prevFastaSeqId\n";
		print $seqH "$fastaSeq\n";
		print $qualH "$prevFastaSeqId\n";
		print $qualH "$qualSeq\n";
		close($seqH);
		close($qualH);
		$ProcessingQueue->enqueue("$uniqFolder/part_seq_$id"."\t"."$uniqFolder/part_qual_$id");
		
	close($iH);
	close($qH);
	$ProcessingQueue->enqueue(undef);
}

sub fireMyJob {
	my $lineCount = $_[0];
	my $fileName = $ProcessingQueue->dequeue();
	return undef if(!defined($fileName));
	my ($file1, $file2) = split(/\t/, $fileName);
	my @idArr = ();
	my @seqArr = ();
	my @qualArr = ();
	open(CHKS, "<$file1") or die "Can't open chunk file containing input reads for processing: $file1\n";
	open(CHKQ, "<$file2") or die "Can't open chunk file containing input quality for processing: $file2\n";
	while(my $id = <CHKS>) {
		<CHKQ>;
		chomp $id;
		my $seq = <CHKS>;
		my $qual = <CHKQ>;
		chomp $seq;
		chomp $qual;
		push(@idArr, $id);
		push(@seqArr, $seq);
		push(@qualArr, $qual);
		$$lineCount++;
	}
	close(CHKS);
	close(CHKQ);
	my ($id) = $file1=~/(\d+)$/;			
	my @reads = (\@idArr, \@seqArr, \@qualArr, $id);
	my $thRef = threads->create('passSeq', @reads);
	return $thRef;
}


sub threading4Processing {
	my @thArr = ();
	my $done = 0;
	my $processedSeqCount = 0;
	my $roughSeqCounter = 0;
	my $sCounter = 0;
	my $jobCounter = 0;
	my $ttlJobCounter = 0;
	my $fileEOF = 0;

	while(1) {
		if($processedSeqCount % 10000 == 0) {
			my $tmpP = sprintf "%0.0f", ($processedSeqCount/$seqCount*100);
			print "$indOfAnalysis: Number of reads processed: " . $processedSeqCount . "/$seqCount ($tmpP\%)...\n";
		}
		my $i;
		for($i=0; $i<$noOfProcesses; $i++) {
			my $thRef = fireMyJob(\$processedSeqCount);
			if(!defined($thRef)) {
				$done = 1;
				last;
			}
			$thArr[$i] = $thRef;
		}
		for(my $j=0; $j<$i; $j++) {
	        my $refArr = $thArr[$j]->join;
	        &updateData(@{$refArr});
		}
		last if($done);
	}
}

sub passSeq {
	yield;
	my ($idArrRef, $seqArrRef, $qualArrRef, $id) = @_;
	resetVariables();
	open(PSEQ, ">$uniqFolder/part_seq_$id"."_out") or die "Can not open part_seq_$id"."_out file\n";
	open(PQUAL, ">$uniqFolder/part_qual_$id"."_out") or die "Can not open part_qual_$id"."_out file\n";
	$c++;
	my @idArr = ();
	my @seqArr = ();
	my @qualArr = ();
	for(my $i=0; $i<@{$idArrRef}; $i++) {
		processSeq($$idArrRef[$i], $$seqArrRef[$i], $$qualArrRef[$i]);
	}
	close(PSEQ);
	close(PQUAL);
	my @retArrRef = ($cmaxRawLen, $cminRawLen, \@crawLen, $ctotalBases, $ctotalHQBases, $cavgQual, $clt100, $ctrimCount, $chQCount, $ctotalBasesAfterHQ, $cmaxHQLen, $cminHQLen, $cavgHQLen, \@chQLen, $ctotalReadsFinal, $ctotalBasesFinal, $ctotalHQBasesFinal, $cavgQualFinal, $clQCount, $ctotalHQBasesAfterHQ, $ctotalValidReadsWithPriAda, $ctotalValidReadsNoPriAda, \@clenDistrib, \@cqualDistrib, \@cgcDistrib, \@ccharCount);
	return \@retArrRef;
}

sub processSeq {
	my ($prevFastaSeqId, $fastaSeq, $qualSeq) = @_;
	$fastaSeq =~ s/\s//g;
	my $len = length $fastaSeq;
	$cmaxRawLen = max($cmaxRawLen, $len);
	$cminRawLen = min($cminRawLen, $len);
	push(@crawLen, $len);
	$qualSeq =~ s/\s+$//;				# To remove the last space added in 'else' part;
	my @tmpArr = getQualBases($qualSeq);
	$ctotalBases += $tmpArr[0];
	$ctotalHQBases += $tmpArr[1];
	$cavgQual += $tmpArr[2];
	$clenDistrib[0][getIndex($len,$lenInterval)]++;
	$cqualDistrib[0][getIndex($tmpArr[2],$qualInterval)]++;
	my $As = $fastaSeq =~ s/A/A/gi;
	my $Ts = $fastaSeq =~ s/T/T/gi;
	my $Gs = $fastaSeq =~ s/G/G/gi;
	my $Cs = $fastaSeq =~ s/C/C/gi;
	my $Ns = $len - $As - $Ts - $Gs - $Cs;
	my $gcPercent = ($Gs + $Cs)/$len*100;
	$cgcDistrib[0][getIndex($gcPercent,$gcInterval)]++;
	$ccharCount[0][0] += $As;
	$ccharCount[0][1] += $Ts;
	$ccharCount[0][2] += $Gs;
	$ccharCount[0][3] += $Cs;
	$ccharCount[0][4] += $Ns;
	if(length $fastaSeq < $lowestValidLen && $isLenFilterOn) {
		$clt100++;
	}
	else {
		if($homoPolyLen != 0) {
			if(hasPolyChar(\$fastaSeq)) {
				$ctrimCount++;
				if(length $fastaSeq >= $lowestValidLen || !$isLenFilterOn) {
					$qualSeq = trimQualSeq($qualSeq, length $fastaSeq, -1);
				}
			}
		}
		if(length $fastaSeq < $lowestValidLen && $isLenFilterOn) {
			$clt100++;
		}
		else {
			if(isReadOfHQ($qualSeq)) {
				$chQCount++;
				$ctotalBasesAfterHQ += length $fastaSeq;
				if(defined $priAdaLib) {
					my $t=isWOPriAda(\$fastaSeq);
					if($t > -1) {
						$qualSeq = trimQualSeq($qualSeq, length $fastaSeq, $t);
					}
					if(length $fastaSeq < $lowestValidLen && $isLenFilterOn) {
						$clt100++;
					}
					else {
						my $len = length $fastaSeq;
						$cmaxHQLen = max($cmaxHQLen, $len);
						$cminHQLen = min($cminHQLen, $len);
						$cavgHQLen += $len;
						push(@chQLen, $len);
						$ctotalReadsFinal++;
						@tmpArr = getQualBases($qualSeq);
						$ctotalBasesFinal += $tmpArr[0];
						$ctotalHQBasesFinal += $tmpArr[1];
						$cavgQualFinal += $tmpArr[2];
						if(!defined($isOnlyStat)) {
							$clenDistrib[1][getIndex($len,$lenInterval)]++;
							$cqualDistrib[1][getIndex($tmpArr[2],$qualInterval)]++;
							my $As = $fastaSeq =~ s/A/A/gi;
							my $Ts = $fastaSeq =~ s/T/T/gi;
							my $Gs = $fastaSeq =~ s/G/G/gi;
							my $Cs = $fastaSeq =~ s/C/C/gi;
							my $Ns = $len - $As - $Ts - $Gs - $Cs;
							my $gcPercent = ($len)?(($Gs + $Cs)/$len*100):0;
							$cgcDistrib[1][getIndex($gcPercent,$gcInterval)]++;
							$ccharCount[1][0] += $As;
							$ccharCount[1][1] += $Ts;
							$ccharCount[1][2] += $Gs;
							$ccharCount[1][3] += $Cs;
							$ccharCount[1][4] += $Ns;
							print PSEQ "$prevFastaSeqId\n";
							print PSEQ formatSeq($fastaSeq), "\n";
							print PQUAL "$prevFastaSeqId\n";
							print PQUAL formatQualSeq($qualSeq), "\n";
						}
					}
				}
				else {
					my $len = length $fastaSeq;
					$cmaxHQLen = max($cmaxHQLen, $len);
					$cminHQLen = min($cminHQLen, $len);
					$cavgHQLen += $len;
					push(@chQLen, $len);
					$ctotalReadsFinal++;
					@tmpArr = getQualBases($qualSeq);
					$ctotalBasesFinal += $tmpArr[0];
					$ctotalHQBasesFinal += $tmpArr[1];
					$cavgQualFinal += $tmpArr[2];
					if(!defined($isOnlyStat)) {
						$clenDistrib[1][getIndex($len,$lenInterval)]++;
						$cqualDistrib[1][getIndex($tmpArr[2],$qualInterval)]++;
						my $As = $fastaSeq =~ s/A/A/gi;
						my $Ts = $fastaSeq =~ s/T/T/gi;
						my $Gs = $fastaSeq =~ s/G/G/gi;
						my $Cs = $fastaSeq =~ s/C/C/gi;
						my $Ns = $len - $As - $Ts - $Gs - $Cs;
						my $gcPercent = ($Gs + $Cs)/$len*100;
						$cgcDistrib[1][getIndex($gcPercent,$gcInterval)]++;
						$ccharCount[1][0] += $As;
						$ccharCount[1][1] += $Ts;
						$ccharCount[1][2] += $Gs;
						$ccharCount[1][3] += $Cs;
						$ccharCount[1][4] += $Ns;
						print PSEQ "$prevFastaSeqId\n";
						print PSEQ formatSeq($fastaSeq), "\n";
						print PQUAL "$prevFastaSeqId\n";
						print PQUAL formatQualSeq($qualSeq), "\n";
					}
				}
			}
			else {
				$clQCount++;
			}
		}
	}
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

sub addTwoArrays {
	my $arr1Ref = $_[0];
	my $arr2Ref = $_[1];
    my $c=0;
    my $i=0;
    foreach my $arrRef (@{$arr1Ref}) {
    	$c=0;
	    foreach my $val (@{$arrRef}) {
	    	@{$$arr2Ref[$i]}[$c] += $val if($val);
	    	@{$$arr2Ref[$i]}[$c] = 0 if(!defined(@{$$arr2Ref[$i]}[$c]));
	    	$c++;
	    }
    	$i++
    }
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
		$ctotalHQBasesAfterHQ += $validBaseCount;
		return 1;				# Return true.
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
	print "----------------------------- Processing Options -----------------------------\n";
	print "  -c | -cpus <Integer>\n";
	print "    Number of CPUs to be used\n";
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
	my @stat = ();
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
		$ctotalValidReadsWithPriAda++;
		return $priAdaStart;
	}
	else {
		$ctotalValidReadsNoPriAda++;
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
