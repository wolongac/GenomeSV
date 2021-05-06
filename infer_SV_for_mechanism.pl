#########################################################################
#      File Name: infer_SV_for_mechanism.pl
#    > Author: hwlu
#    > Mail: hongweilu@genetics.ac.cn 
#    Created Time: Wed 15 Jul 2020 01:56:45 AM CST
#########################################################################

#!/usr/bin/perl -w
use strict;

my $sv=shift; #02428.syri.out.SV.bed.mergeDel.CPV2InDel
my $infer_result=shift;	    #Merge_Syri_bed.result.0.9.0.9.rename.tag.infer_1.infer_2.Add.vntr
my $spe = (split /\//,$sv)[-1];
$spe =~ s/\..*//g;

print "$spe\n";
my %hash;
open(IN,"$infer_result");
while(<IN>){
    chomp;
    my @line = split /\s+/,$_;
    next if ($line[4] ne $spe);
    $hash{$line[1]}{$line[2]}{$line[3]}{$line[0]}=$line[5];
    #print "$line[1]\t$line[2]\t$line[3]\t$line[0]\t$line[5]\n";
}
close IN;

open(IN,"$sv");
open(OUT,">$sv.infered");
while(<IN>){
    chomp;
    my $line = $_;
    my @line = split /\s+/,$line;
    if($line[6] eq "Deletion"){
	$line[6]="DEL";
    }
    if($line[6] eq "Insertion"){
	$line[6] = "INS";
    }
    my $type="small:$line[6]";
    my $sv_length;
    if($line[1]>$line[2]){
	($line[1],$line[2]) = ($line[2],$line[1]);
    }
    if($line[6] eq "DEL" or $line[6] eq "INS"){
	$sv_length=abs(abs($line[2]-$line[1])-abs($line[5]-$line[4]));
    }else{
	$sv_length=abs($line[2]-$line[1]);
    }
    if(exists $hash{$line[0]}{$line[1]}{$line[2]}{$line[6]} and $sv_length>=50){
	$type=$hash{$line[0]}{$line[1]}{$line[2]}{$line[6]};
    }
    print OUT "$line[0]\t$line[1]\t$line[2]\t$line[3]\t$line[4]\t$line[5]\t$type\n";
}
close IN;
close OUT;
