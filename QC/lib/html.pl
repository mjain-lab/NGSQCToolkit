sub htmlPrint{
	my ($progPath, $prog, $htF, $iFol, $isPairedEnd, $isOnlyStat, $inpFs, $seqFormatName, $statFile, $oFol, $fileNames4HTML) = @_;
	my $imgPath = $progPath . "lib/imgs";
	my $cssPath = $progPath . "lib";
	my $analMsg1 = ($isOnlyStat)?"":"and filtering";
	my $analMsg2 = ($isPairedEnd)?"paired end data":"single end data";
	my ($file1, $file2) = split(":::::", $inpFs);
	open(I, "<$statFile") or print "Can not open statistics file: $statFile\n";
	my @statFData = <I>;
	close(I);
	my $statFileOnlyName = getFileName($statFile);
	my ($t, $priAdaLib, $cutLen, $cutQual, $nCPUs, $onlySOnOff);
	($t, $t, $priAdaLib) = split(/ {2,}|\t/, $statFData[3]);
	($t, $t, $cutLen) = split(/ {2,}|\t/, $statFData[4]);
	($t, $t, $cutQual) = split(/ {2,}|\t/, $statFData[5]);
	($t, $t, $nCPUs) = split(/ {2,}|\t/, $statFData[7]);
	$onlySOnOff = ($isOnlyStat)?"On":"Off";
	my ($outFile1, $outFile2, $unPaired, $avgQF1, $avgQF2, $avgQRangeF1, $avgQRangeF2, $baseCntF1, $baseCntF2, $gcDistF1, $gcDistF2, $qualDistF1, $qualDistF2, $sumPieF);
	($outFile1, $avgQF1, $baseCntF1, $gcDistF1, $qualDistF1, $sumPieF, $avgQRangeF1) = @$fileNames4HTML;
	($outFile1, $outFile2, $unPaired, $avgQF1, $avgQF2, $baseCntF1, $baseCntF2, $gcDistF1, $gcDistF2, $qualDistF1, $qualDistF2, $sumPieF, $avgQRangeF1, $avgQRangeF2) = @$fileNames4HTML if($isPairedEnd);
	$outFile1 = getFileName($outFile1);
	$outFile2 = getFileName($outFile2) if($isPairedEnd);
	$unPaired = getFileName($unPaired) if($isPairedEnd);
	my ($avgQRangeRF1, $avgQRangeFF1) = split(":::", $avgQRangeF1);
	my ($avgQRangeRF2, $avgQRangeFF2) = split(":::", $avgQRangeF2);
	my $inpFilesMsg;
	$inpFilesMsg = "input file,  $file1(A),";
	$inpFilesMsg = "both input files,  $file1(A) and $file2(B)," if($isPairedEnd);
	my $b4A4Msg;
	$b4A4Msg = "before and after QC";
	$b4A4Msg = "before QC" if($isOnlyStat);
	#### Getting current time
	my @months = qw(Jan Feb Mar Apr May Jun Jul Aug Sep Oct Nov Dec);
	my @weekDays = qw(Sun Mon Tue Wed Thu Fri Sat Sun);
	my ($second, $minute, $hour, $dayOfMonth, $month, $yearOffset, $dayOfWeek, $dayOfYear, $daylightSavings) = localtime();
	my $year = 1900 + $yearOffset;
	#my $theTime = "$hour:$minute:$second, $weekDays[$dayOfWeek] $months[$month] $dayOfMonth, $year";
	my $theTime = "$weekDays[$dayOfWeek] $months[$month] $dayOfMonth, $year";
	open(O,">$htF") or die "Can not create HTML file: $htF\n";
	open(V, $progPath."lib/version");
	my $version = <V>;
	close(V);

print O <<EOF;
	<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN">
	<html>
	<head>
	<title>NGS QC Toolkit</title>
	<meta http-equiv="Content-Type" content="text/html; charset=iso-8859-1">
	
	
	<style>
		BODY {
			margin-top: 0;
			background-repeat: repeat;
			font-family: Arial, Helvetica, sans-serif;
			font-size: 5px;
			margin-bottom: 0;
		}
		.cnt {
			font-family: Verdana, Arial, Helvetica, sans-serif;
			font-size: 12px;
			line-height: 20px;
			padding: 5px 20px 0px;
		}
		Table .exp {
			font-family: Verdana, Arial, Helvetica, sans-serif;
			font-size: 13px;
			line-height: 20px;
		}
		TD .padding {
			padding: 0px 20px;
		}
		.head1 {
			font-size: 18px;
			font-weight: bold;
			padding: 5px 0px;
		}
		.head2 {
			font-size: 14px;
			font-weight: bold;
			padding: 5px 0px;
		}
		.head3 {
			font-size: 12px;
			font-weight: bold;
		}
		A {
			text-decoration: none;
			color: #0000FF;
		}
		.tblBg TABLE TD {
			background-color: #EEEEEE;
		}
		.tblBg2 TABLE TD {
			background-color: #E1E1E1;
		}
	</style>

	
	
	</head>
	
	<body bgcolor="#bbbbbb">
	<table align="center" width="900" cellspacing="0" cellpadding="0" bgcolor="#ffffff" border="0">
	
		<tr>
			<td width="17" rowspan="6">&nbsp;</td>
			<td height="150" bgcolor="#E1E1E1"><center><a href="http://www.nipgr.res.in/ngsqctoolkit.html"><b><font  style="font-size: 50px;">NGS QC T</font><font style="font-size: 40px;">OOLKIT</font></b></center></a>
			</td>
			<td width="17" rowspan="6"></td>
		</tr>
		<tr>
			<td valign="top" class="tblBg">
			<div class="cnt">
				<table class="cnt" width="100%" border="0">
					<tr>
						<td class="head1">Results of quality control (QC) using $prog v$version <font style="font-size:10px">($theTime)</font></td>
					</tr>
					<tr>
						<td class="head2">Input files and parameters:</td>
					</tr>
					<tr><td class="tblBg2">
						<table width="100%" border="0" class="cnt">
						<tr>
							<td>Analysis type</td><td>Quality check $analMsg1 of $analMsg2</td>
						</tr>
						<tr>
							<td>Input file directory</td><td><a href="file://$iFol" target="_blank">$iFol</a></td>
						</tr>
						<tr>
							<td>Input file 1</td><td>$file1</td>
						</tr>
EOF
if($isPairedEnd) {
print O <<EOF;
						<tr>
							<td>Input file 2</td><td>$file2</td>
						</tr>
EOF
}
print O <<EOF;
						<tr>
							<td>Input file format</td><td>FASTQ ($seqFormatName variant)</td>
						</tr>
						<tr>
							<td>Primer/Adaptor library</td><td>$priAdaLib</td>
						</tr>
						<tr>
							<td>Cut-off read length for HQ</td><td>$cutLen</td>
						</tr>
						<tr>
							<td>Cut-off quality score</td><td>$cutQual</td>
						</tr>
						<tr>
							<td>Only statistics</td><td>$onlySOnOff</td>
						</tr>
						<tr>
EOF
print O "
							<td>Number of ". (($prog=~/_PRLL/)?"CPUs":"processes") ."</td><td>$nCPUs</td>
\n";
print O <<EOF;
						</tr>
						</table>
					</td>
					</tr>
					<tr>
						<td class="head2">Output files:</td>
					</tr>
					<tr><td class="tblBg2">
						<table width="100%" border="0" class="cnt">
						<tr>
							<td>Output folder</td><td><a href="file://$oFol" target="_blank">$oFol</a></td>
						</tr>
						<tr>
							<td>QC statistics</td><td><a href="$statFileOnlyName" target="_blank">$statFileOnlyName</a></td>
						</tr>
EOF
if(!$isOnlyStat) {
print O <<EOF;
						<tr>
							<td>High quality filtered file 1</td><td>$outFile1</td>
						</tr>
EOF
if($isPairedEnd) {
print O <<EOF;
						<tr>
							<td>High quality filtered file 2</td><td>$outFile2</td>
						</tr>
						<tr>
							<td>High quality filtered un-paired data</td><td>$unPaired</td>
						</tr>
EOF
}
}
print O <<EOF;
						<tr>
							<td>Per base average quality score for file 1</td><td><a href="$avgQF1" target="_blank">$avgQF1</a></td>
						</tr>
EOF
if($isPairedEnd) {
print O <<EOF;
						<tr>
							<td>Per base average quality score for file 2</td><td><a href="$avgQF2" target="_blank">$avgQF2</a></td>
						</tr>
EOF
}
print O <<EOF;
						<tr>
							<td>Read count (%) per base for different quality ranges for file 1 (Raw)</td><td><a href="$avgQRangeRF1" target="_blank">$avgQRangeRF1</a></td>
						</tr>
EOF
if(!$isOnlyStat) {
print O <<EOF;
						<tr>
							<td>Read count (%) per base for different quality ranges for file 1 (Filtered)</td><td><a href="$avgQRangeFF1" target="_blank">$avgQRangeFF1</a></td>
						</tr>
EOF
}
if($isPairedEnd) {
print O <<EOF;
						<tr>
							<td>Read count (%) per base for different quality ranges for file 2 (Raw)</td><td><a href="$avgQRangeRF2" target="_blank">$avgQRangeRF2</a></td>
						</tr>
EOF
if(!$isOnlyStat) {
print O <<EOF;
						<tr>
							<td>Read count (%) per base for different quality ranges for file 2 (Filtered)</td><td><a href="$avgQRangeFF2" target="_blank">$avgQRangeFF2</a></td>
						</tr>
EOF
}
}
print O <<EOF;
						<tr>
							<td>Base composition for file 1</td><td><a href="$baseCntF1" target="_blank">$baseCntF1</a></td>
						</tr>
EOF
if($isPairedEnd) {
print O <<EOF;
						<tr>
							<td>Base composition for file 2</td><td><a href="$baseCntF2" target="_blank">$baseCntF2</a></td>
						</tr>
EOF
}
print O <<EOF;
						<tr>
							<td>GC content distribution for file 1</td><td><a href="$gcDistF1" target="_blank">$gcDistF1</a></td>
						</tr>
EOF
if($isPairedEnd) {
print O <<EOF;
						<tr>
							<td>GC content distribution for file 2</td><td><a href="$gcDistF2" target="_blank">$gcDistF2</a></td>
						</tr>
EOF
}
print O <<EOF;
						<tr>
							<td>Quality distribution for file 1</td><td><a href="$qualDistF1" target="_blank">$qualDistF1</a></td>
						</tr>
EOF
if($isPairedEnd) {
print O <<EOF;
						<tr>
							<td>Quality distribution for file 2</td><td><a href="$qualDistF2" target="_blank">$qualDistF2</a></td>
						</tr>
EOF
}
print O <<EOF;
						<tr>
							<td>Summary of QC</td><td><a href="$sumPieF" target="_blank">$sumPieF</a></td>
						</tr>
						</table>
					</td>
					</tr>
					<tr>
						<td class="head2" style="background-color: #ffffff;">&nbsp;</td>
					</tr>
					<tr>
						<td class="head2">Results of QC</td>
					</tr>
EOF
my $flag = 0;
for(my $i=0; $i<@statFData; $i++) {
	my $line = $statFData[$i];
	if($line =~ /^QC statistics/) {
		$flag = 1;
		print O "<tr><td class=\"head3\">QC statistics</td></tr>\n<tr><td class=\"tblBg2\"><table width=\"100%\" border=\"0\" class=\"cnt\">\n";
		next;
	}
	if($line =~ /^Detailed QC statistics/) {
		$flag = 2;
		print O "</table></td></tr>\n<tr><td class=\"head3\">Detailed QC statistics</td></tr>\n<tr><td class=\"tblBg2\"><table width=\"100%\" border=\"0\" class=\"cnt\">\n";
		next;
	}
	if($line =~ /^Average quality score/) {
		$flag = 0;
		next;
	}	
	chomp($line);
	if($flag != 0 && $line ne "") {
		my @clms = split(/ {2,}|\t/, $line);
		shift(@clms);
		print O "<tr>";
		foreach my $f (@clms) {
			print O "<td>$f</td>";
		}
		print O "</tr>\n";		
	}
}
print O "</table></td></tr>\n";
print O <<EOF;							
					<tr>
						<td class="head3">Summary of QC</td>
					</tr>
					<tr>
						<td>
							<table width="100%" border="0" class="cnt">
								<tr>
									<td><div align="justify">Following pie chart shows the summary of QC depicting percentage of high quality, low quality and contaminated reads.</div></td>
								</tr>
								<tr>
									<td align="center"><img src="$sumPieF" border="1"><br><b>(A)</b></td>
								</tr>
							</table>
						</td>
					</tr>
					<tr>
						<td class="head3">Per base average quality scores</td>
					</tr>
					<tr>
						<td>
							<table width="100%" border="0" class="cnt">
								<tr>
									<td><div align="justify">Following graph(s) show per base average PHRED quality scores for $inpFilesMsg $b4A4Msg.</div></td>
								</tr>
								<tr>
									<td align="center"><img src="$avgQF1" border="1"><br><b>(A)</b>
EOF
if($isPairedEnd) {
print O <<EOF;
									<br><br><img src="$avgQF2" border="1"><br><b>(B)</b></td>
EOF
}
print O <<EOF;
								</tr>
							</table>
						</td>
					</tr>
					<tr>
						<td class="head3">Read count (%) per base for different quality score ranges</td>
					</tr>
					<tr>
						<td>
							<table width="100%" border="0" class="cnt">
								<tr>
									<td><div align="justify">Following graph(s) show per base read count (%) for different quality score ranges for $inpFilesMsg $b4A4Msg.</div></td>
								</tr>
								<tr>
									<td align="center"><img src="$avgQRangeRF1" border="1">
EOF
if(!defined($isOnlyStat)) {
print O <<EOF;
									<br><br><img src="$avgQRangeFF1" border="1">
EOF
}
print O <<EOF;
									<br><b>(A)</b>
EOF
if($isPairedEnd) {
print O <<EOF;
									<br><br><img src="$avgQRangeRF2" border="1">
EOF
if(!defined($isOnlyStat)) {
print O <<EOF;
									<br><br><img src="$avgQRangeFF2" border="1">
EOF
}
print O <<EOF;
									<br><b>(B)</b></td>
EOF
}
print O <<EOF;
								</tr>
							</table>
						</td>
					</tr>
					<tr>
						<td class="head3">GC content distribution</td>
					</tr>
					<tr>
						<td>
							<table width="100%" border="0" class="cnt">
								<tr>
									<td><div align="justify">Following graph(s) show number of reads for distinct average GC content (%) ranges for $inpFilesMsg $b4A4Msg.</div></td>
								</tr>
								<tr>
									<td align="center"><img src="$gcDistF1" border="1"><br><b>(A)</b>
EOF
if($isPairedEnd) {
print O <<EOF;
									<br><br><img src="$gcDistF2" border="1"><br><b>(B)</b></td>
EOF
}
print O <<EOF;
								</tr>
							</table>
						</td>
					</tr>
					<tr>
						<td class="head3">Quality distribution</td>
					</tr>
					<tr>
						<td>
							<table width="100%" border="0" class="cnt">
								<tr>
									<td><div align="justify">Following graph(s) show number of reads for different average PHRED quality scores for $inpFilesMsg $b4A4Msg.</div></td>
								</tr>
								<tr>
									<td align="center"><img src="$qualDistF1" border="1"><br><b>(A)</b>
EOF
if($isPairedEnd) {
print O <<EOF;
									<br><br><img src="$qualDistF2" border="1"><br><b>(B)</b></td>
EOF
}
print O <<EOF;
								</tr>
							</table>
						</td>
					</tr>
					<tr>
						<td class="head3">Base composition</td>
					</tr>
					<tr>
						<td>
							<table width="100%" border="0" class="cnt">
								<tr>
									<td><div align="justify">Following graph(s) show base composition (count) for $inpFilesMsg $b4A4Msg with percentage of bases at the bottom.</div></td>
								</tr>
								<tr>
									<td align="center"><img src="$baseCntF1" border="1"><br><b>(A)</b>
EOF
if($isPairedEnd) {
print O <<EOF;
									<br><br><img src="$baseCntF2" border="1"><br><b>(B)</b></td>
EOF
}
print O <<EOF;
								</tr>
							</table>
						</td>
					</tr>
				</table>
			</div>
			</td>
		</tr>
		<tr><td style="font-size: 11px;" align="center" valign="middle" height="25" background="$imgPath/btmLine.png" bgcolor="#EEEEEE">For Questions and Suggestions, contact <a href="mailto:mjain\@nipgr.res.in">Mukesh Jain (mjain\@nipgr.res.in)</a>; <a href="mailto:ravipatel\@nipgr.res.in">Ravi Patel (ravipatel\@nipgr.res.in)</a></td></tr>
	
		<!--- <tr><td colspan="3"><img src="$imgPath/down.png"></td></tr> --->
	</table>
	</body>
	</html>
EOF

close(O);
}
1;
