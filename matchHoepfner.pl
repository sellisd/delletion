#!/usr/bin/perl
use warnings;
use strict;
#match HIP and HOP experiments by compound id and concentration
#first pass #read both files and save first line
#  make hash (cmbid.conc => [zscore.index, zscore.index])
#for each cmbid.conc shared across HIP and HOP print
# hip indexes hipZ indexes hop indexes and hopZ indexes

open HIP, "../Hoepfner/HIP_scores.txt" or die $!;
open HOP, "../Hoepfner/HOP_scores.txt" or die $!;
my %hop;
my %hip;
my $HIPline = readline(HIP);
chomp $HIPline;
my @HIPar = split "\t", $HIPline;
shift @HIPar;

my $HOPline = readline(HOP);
chomp $HOPline;
my @HOPar = split "\t", $HOPline;
shift @HOPar;

my %groups;
my $index = 2;
foreach my $entry (@HIPar){
    die unless $entry =~ /"(Ad.|MADL) scores for Exp. (\d+)_([\d.]+)_(HIP|HOP)_([\dAB]+)\s*(z\-score|)"/;
    my $adj = $1;
    my $cmbid = $2; #compound id
    my $conc = $3; #concentration
    my $type = $4; #HIP|HOP
    my $study = $5; #study number
    my $suffix = $6; #z-score or MADL
    my $key = $cmbid.'_'.$conc;
    my $zscore =($suffix eq 'z-score'?1:0);
    my $value = $zscore.'_'.$index;
    if(defined($hip{$key})){
	push @{$hip{$key}},$value;
    }else{
	$hip{$key} = [$value];
    }
    $index++;
}

$index = 2;
foreach my $entry (@HOPar){
    die unless $entry =~ /"(Ad.|MADL) scores for Exp. (\d+)_([\d.]+)_(HIP|HOP)_([\dAB]+)\s*(z\-score|)"/;
    my $adj = $1;
    my $cmbid = $2; #compound id
    my $conc = $3; #concentration
    my $type = $4; #HIP|HOP
    my $study = $5; #study number
    my $suffix = $6; #z-score or MADL
   my $key = $cmbid.'_'.$conc;
    my $zscore =($suffix eq 'z-score'?1:0);
    my $value = $zscore.'_'.$index;
    if(defined($hop{$key})){
	push @{$hop{$key}},$value;
    }else{
	$hop{$key} = [$value];
    }
    $index++;
}

foreach my $hopExp (keys %hop){
    $hopExp =~ /^(.*?)_(.*?)$/;
    my $cmbid = $1;
    my $conc = $2;
#    print $cmbid,' ', $conc,"\n";
#    die $hopExp;
    foreach my $hipExp (keys %hip){
	if ($hopExp eq $hipExp){
	    foreach my $hopI (@{$hop{$hopExp}}){
		$hopI =~ /(.*?)_(.*)/;
		my $Izscore = $1;
		my $Iindex = $2;
		print $cmbid,"\t";
		print $conc,"\t";
		print $Izscore,"\t";
		print "HOP\t";
		print $Iindex,"\n";
	    }
	    foreach my $hipI (@{$hip{$hipExp}}){
		$hipI =~ /(.*?)_(.*)/;
		my $Izscore = $1;
		my $Iindex = $2;
		print $cmbid,"\t";
		print $conc,"\t";
		print $Izscore,"\t";
		print "HIP\t";
		print $Iindex,"\n";
	    }
	}
    }
}
