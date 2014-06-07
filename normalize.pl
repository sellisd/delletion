#!/usr/bin/perl
use warnings;
use strict;

#for each het/hom the corresponding control group and normalize them. Then calculate the difference from the het/hom average of the control group
my $hetfile ='../fitdb/het.txt';
my $homfile ='../fitdb/hom.txt';
my $homControlfile = '../fitdb/hom_controls.txt';
my $hetControlfile = '../fitdb/het_controls.txt';

open HET, $hetfile or die $!;
open HOM, $homfile or die $!;
open HETC, $hetControlfile or die $!;
open HOMC, $homControlfile or die $!;

# reat control files and build hashes for control sets and filenames
my %hetC;
my $lineNumber = 0;
while(my $line = <HETC>){
    chomp $line;
    $lineNumber++;
    next if $lineNumber == 1;
    (my $control_set, my $filename) = split "\t", $line;
    if(defined($hetC{$control_set})){
	push @{$hetC{$control_set}},$filename;
    }else{
	$hetC{$control_set} = [$filename];	
    }
}

my %homC;
$lineNumber = 0;
while(my $line = <HOMC>){
    chomp $line;
    $lineNumber++;
    next if $lineNumber == 1;
    (my $control_set, my $filename) = split "\t", $line;
    if(defined($homC{$control_set})){
	push @{$homC{$control_set}},$filename;
    }else{
	$homC{$control_set} = [$filename];	
    }
}
close HETC;
close HOMC;

# read homozygote files normalize and calculate averages
$lineNumber = 0;
if(0){my $line;
#while (my $line = <HOM>){
    chomp $line;
    $lineNumber++;
    next if $lineNumber == 1;
    (my $filename, my $condition, my $control_set) = split "\t", $line;
    my @hom = @{$homC{$control_set}};
    my $commandString = '../fitness_scripts/normalize_data.pl ';
    $commandString .= "../fitdb/sgtc_hom/$filename.raw ";
    my @av;
    foreach (@hom){
	$commandString .= "../fitdb/sgtc_hom_controls/$_.raw ";
	push @av, "../fitdb/sgtc_hom_controls/$_.norm";
    }
    print " $commandString\n";
#    system "$commandString";
    #average
    (my $homCAvRef, my $strainTagsRef) = &averageC(\@av);
    #calculate difference (meanControl - hom)
    #read normlized hom
    my @normHom;
    my @normHomT; #threshold
    my @strainTags;
    open HOMN, "../fitdb/sgtc_hom/$filename.norm" or die $!;
    my $lineCounter = 0;
    while (my $l = <HOMN>){
	chomp $l;
	$lineCounter++;
	next if $lineCounter==1;
	(my $strain_tag, my $norm_quantile, my $above_thresh) = (split "\t", $l)[0,1,3];
	push @normHom, $norm_quantile;
	push @normHomT, $above_thresh;
	push @strainTags, $strain_tag;
    }
    close HOMN;
    die unless (@{$strainTagsRef} ~~ @strainTags); #perl >= 5.10.1
    #calculate difference
    my @diff;
    for(my $i = 0; $i <= $#normHomT; $i++){
	my $above = $normHomT[$i];
	if ($above == 1){
	    if(${$homCAvRef}[$i] ne "NA"){
		push @diff, ${$homCAvRef}[$i] - $normHom[$i];
	    }else{
		push @diff, "NA";
	    }
	}elsif($above == 0){
	    push @diff, "NA";
	}else{
	    die;
	}
    }
    #save output in new file $file.dnorm
    open HOMOUT, ">../fitdb/sgtc_hom/$filename.dnorm";
    print HOMOUT "strain_tag\tdiff\n";
    for(my $i = 0; $i <= $#strainTags; $i++){
	print HOMOUT $strainTags[$i],"\t",$diff[$i],"\n";
    }
    close HOMOUT;
}
close HOM;

#repeat for the heterozygotes
$lineNumber = 0;
while (my $line = <HET>){
    chomp $line;
    $lineNumber++;
    next if $lineNumber == 1;
    (my $filename, my $condition, my $control_set) = split "\t", $line;
    my @het = @{$hetC{$control_set}};
    my $commandString = '../fitness_scripts/normalize_data.pl ';
    $commandString .= "../fitdb/sgtc_het/$filename.raw ";
    my @av;
    foreach (@het){
	$commandString .= "../fitdb/sgtc_het_controls/$_.raw ";
	push @av, "../fitdb/sgtc_het_controls/$_.norm";
    }
    print " $commandString\n";
    system "$commandString";
    #average
    (my $hetCAvRef, my $strainTagsRef) = &averageC(\@av);
    #calculate difference (meanControl - het)
    #read normlized het
    my @normHet;
    my @normHetT; #threshold
    my @strainTags;
    open HETN, "../fitdb/sgtc_het/$filename.norm" or die $!;
    my $lineCounter = 0;
    while (my $l = <HETN>){
	chomp $l;
	$lineCounter++;
	next if $lineCounter==1;
	(my $strain_tag, my $norm_quantile, my $above_thresh) = (split "\t", $l)[0,1,3];
	push @normHet, $norm_quantile;
	push @normHetT, $above_thresh;
	push @strainTags, $strain_tag;
    }
    close HETN;
    die unless (@{$strainTagsRef} ~~ @strainTags); #perl >= 5.10.1
    #calculate difference
    my @diff;
    for(my $i = 0; $i <= $#normHetT; $i++){
	my $above = $normHetT[$i];
	if ($above == 1){
	    if(${$hetCAvRef}[$i] ne "NA"){
		push @diff, ${$hetCAvRef}[$i] - $normHet[$i];
	    }else{
		push @diff, "NA";
	    }
	}elsif($above == 0){
	    push @diff, "NA";
	}else{
	    die;
	}
    }
    #save output in new file $file.dnorm
    open HETOUT, ">../fitdb/sgtc_het/$filename.dnorm";
    print HETOUT "strain_tag\tdiff\n";
    for(my $i = 0; $i <= $#strainTags; $i++){
	print HETOUT $strainTags[$i],"\t",$diff[$i],"\n";
    }
    close HETOUT;
}
close HOM;

close HET;


sub averageC{
#average quantile normalized values from control files excluding values that are under the threshold
    my $fileRef = shift @_;
    my $i = 0;
    my @m; #matrix of quantile values to be averaged
    my @rowMeans;
    my @strain_tag;
    #build matrix
    foreach my $fileName (@$fileRef){
	my $j = 0;
	open IN, $fileName or die $!;
	my $lineCounter = 0;
	while (my $line = <IN>){
	    $lineCounter++;
	    next if $lineCounter == 1;
	    chomp $line;
	    (my $strain_tag, my $norm_quantile, my $above_threshold) =  (split "\t", $line)[0,1,3];
	    if($fileName eq ${$fileRef}[0]){ #if it is the first file
		push @strain_tag, $strain_tag; #build array
	    }else{
		if($strain_tag ne $strain_tag[$lineCounter-2]){
		    die;
		}
	    }
	    if($above_threshold == 1){
		$m[$j][$i] = $norm_quantile;
	    }elsif ($above_threshold == 0){
		#keep these undefined
	    }else{
		die;
	    }
	    $j++;
	}
	close IN;
	$i++;
    }
    #average by row
    for(my $x = 0; $x <= $#m; $x++){
	my $counter = 0;
	my $sum = 0;
	for(my $y = 0; $y <= $#{$m[0]}; $y++){
	    if(defined($m[$x][$y])){
		$sum+= $m[$x][$y];
		$counter++;
	    }
	}
	if ($counter !=0){
	    push @rowMeans, $sum/$counter;
	}else{
	    push @rowMeans, "NA";
	}
    }
    return (\@rowMeans,\@strain_tag);
}
