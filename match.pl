#!/usr/bin/perl
use warnings;
use strict;
#Hillenmeyer et al,. 2008
#Read data raw data from fitdb and group by condition, generations and scanner used. Then normalize data by a system call to normalize_data.pl script

#match conditions to their control set and hom/het. make two new hashes with filename and control set for the het and hom
#do the same for hetC and homC
#when printing out the final data find the control set for each file name and the corresponging control file

#read het files and
my $hetfile ='../fitdb/het.txt';
my $homfile ='../fitdb/hom.txt';
my $homControlfile = '../fitdb/hom_controls.txt';
my $hetControlfile = '../fitdb/het_controls.txt';

my $outputFile = 'common.dat';
open HET, $hetfile or die $!;
open HOM, $homfile or die $!;
open OUT, '>',$outputFile or die $!;

my %het;
my $lineNumber = 0;
while (my $line = <HET>){
    chomp $line;
    $lineNumber++;
    next if $lineNumber == 1;
    (my $filename, my $condition, my $control_set) = split "\t", $line;
    $condition = uc($condition);
    (my $name, my $scanner, my $generations, my $tag) = split '::', $control_set;
    if(defined($het{$condition.'::'.$generations.'::'.$scanner})){
	push @{$het{$condition.'::'.$generations.'::'.$scanner}},$filename;
    }else{
	$het{$condition.'::'.$generations.'::'.$scanner}=[$filename];
    }
}

my %hom;
$lineNumber = 0; 
while (my $line = <HOM>){
    chomp $line;
    $lineNumber++;
    next if $lineNumber == 1;
    (my $filename, my $condition, my $control_set) = split "\t", $line;
    $condition = uc($condition);
    (my $name, my $scanner, my $generations, my $tag) = split '::', $control_set;
    if(defined($hom{$condition.'::'.$generations.'::'.$scanner})){ #common conditions etc
	push @{$hom{$condition.'::'.$generations.'::'.$scanner}},$filename;
    }else{
	$hom{$condition.'::'.$generations.'::'.$scanner}=[$filename];
    }
}
close HET;
close HOM;

#normalize in groups that share conditions, generations and instrument
my $groupCounter = 0;
foreach my $entries (sort keys %het){
    if(defined($hom{$entries})){ #only print shared conditions
	print OUT "$groupCounter. $entries ";
	print OUT " het: @{$het{$entries}},";
	print OUT " hom: @{$hom{$entries}}";
	print OUT "\n";
	my @het = @{$het{$entries}};
	my @hom = @{$hom{$entries}};
	my $commandString = '../fitness_scripts/normalize_data.pl ';
	foreach (@het){
	    $commandString .= "../fitdb/sgtc_het/$_.raw ";
	}
	foreach (@hom){
	    $commandString .= "../fitdb/sgtc_hom/$_.raw ";
	}
	$commandString .= " --savedistr ../comp/$groupCounter.tab --imgdir ../comp";
	$groupCounter++;
	print $entries,"\t", $#het+1,"\t",$#hom+1,"\n";
###############################
	#Normalize data
#	print $commandString,"\n";
#	system $commandString;
###############################
    }
}
close OUT;
#foreach group if multiple heterozygote or homozygotes present replace with average values
