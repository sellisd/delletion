#!/usr/bin/perl
use warnings;
use strict;
#calculate average normalized values for homozygous/heterozygous files within the same group
#using the quantile normalization and excluding values under the detection threshold
#imput data from normalized files and groupings from common.dat (output of match.pl)
#output are .dat files for each condition with 3 columns, ORFname, homozygote Mean, heterozygote Mean
#also prints in stdout the condition and the overall mean (across all ORFs) for each condition

open IN, 'common.dat' or die $!;
my $outPath = "../analysis/fitComp/";

# read files and make list of common strains
open HOM, '../fitness_scripts/files/straintag_Tag3hom.tab' or die $!;
open HET, '../fitness_scripts/files/straintag_Tag3het.tab' or die $!;
my %common;
my %homORFs;
my $lineCounter = 0;
while (my $line = <HOM>){
    if($lineCounter > 0){
	chomp $line;
	my @ar = split "\t", $line;
	$homORFs{$ar[0]} = $ar[3];
    }
    $lineCounter++;
}
$lineCounter = 0;
while (my $line = <HET>){
    if($lineCounter >0){
	chomp $line;
	my @ar = split "\t", $line;
	if(defined($homORFs{$ar[0]})){
	    $common{$ar[0]} = 1;
	}
    }
    $lineCounter++;
}
close HOM;
close HET;

#read common.dat line by line (one group at a time)
while (my $line = <IN>){
    chomp $line;
    $line =~ /\d+\.\s+(.*) het: ([\d_\s]+), hom:([\d_\s]+)/;
    my $condition = $1;
    my $het = $2;
    my $hom = $3;
    my @het = split " ", $het;
    my @hom = split " ", $hom;
    my %hom;
    my %het;
#parse conditions and make 4 output files for each one 
    $condition =~s/[\s:+_\-\(\),]//g;
    my $outputFile = $condition; 
    mkdir $outPath;
    open OUTUA,'>'.$outPath.$outputFile.'.ua.dat' or die $!;
    open OUTUS,'>'.$outPath.$outputFile.'.us.dat' or die $!;
    open OUTDA,'>'.$outPath.$outputFile.'.da.dat' or die $!;
    open OUTDS,'>'.$outPath.$outputFile.'.ds.dat' or die $!;
#if multiple files calculate averages
#exclude if under threshold build hash for both homozygote and heterozygotes
#print out and calculate means only of shared keys
    foreach my $o (@hom){
	open my $oF, '../fitdb/sgtc_hom/'.$o.'.dnorm' or die $!;
	my $c=0;
	while(my $l = <$oF>){
	    chomp $l;
	    my @ar = split " ", $l;
	    if ($c>0){
		#only include if common ORFs
		$ar[0] =~ /^(.*?::.*?):.*$/;
		next unless (defined($common{$1}));
		#exclude missing values (in diff these are under threshold)
		next if $ar[1] eq 'NA';
		if(defined($hom{$ar[0]})){
		    push @{$hom{$ar[0]}},$ar[1];
		}else{
		    $hom{$ar[0]}=[$ar[1]];
		}
	    }
	    $c++;
	}
	close $oF;
    }
    
    foreach my $e (@het){
	open my $eF, '../fitdb/sgtc_het/'.$e.'.dnorm' or die $!;
	my $c = 0;
	while(my $l = <$eF>){
	    chomp $l;
	    my @ar = split " ", $l;
	    if($c>0){
		#only include if common ORFs
		$ar[0] =~ /^(.*?::.*?):.*$/;
		next unless (defined($common{$1}));
		next if $ar[1] eq 'NA';
		if(defined($het{$ar[0]})){
		    push @{$het{$ar[0]}},$ar[1];
		}else{
		    $het{$ar[0]}=[$ar[1]];
		}
	    }
	    $c++;
	}
    }
#filter by common genes
    my @homList;
    my @hetList;
    foreach my $orf (sort keys %hom){
	if (defined $het{$orf}){#common ORFs
	    my $homMean = &mean($hom{$orf});
	    my $hetMean = &mean($het{$orf});
	    unless($homMean eq 'NA' or $hetMean eq 'NA'){ #only print if both not NA
		if($orf=~/.*:uptag:antisense/){
		    print OUTUA $orf,"\t",$homMean,"\t",$hetMean,"\n";
		}elsif($orf=~/.*:uptag:sense/){
		    print OUTUS $orf,"\t",$homMean,"\t",$hetMean,"\n";
		}elsif($orf=~/.*:downtag:antisense/){
		    print OUTDA $orf,"\t",$homMean,"\t",$hetMean,"\n";
		}elsif($orf=~/.*:downtag:sense/){
		    print OUTDS $orf,"\t",$homMean,"\t",$hetMean,"\n";
		}else{
		    die;
		}
	    }
	    push @homList,$homMean;
	    push @hetList,$hetMean;
	}
    }
    print $condition,"\t",&mean(\@homList),"\t",&mean(\@hetList),"\n";
    close OUTUA;
    close OUTUS;
    close OUTDA;
    close OUTDS;
}
close IN;


sub mean{
    my $arref= shift @_;
    my $sum=0;
    my $c=0;
    foreach (@$arref){
	$sum+=$_;
	$c++;
    }
    return $sum/$c;
}
