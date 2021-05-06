#########################################################################
#      File Name: split_Merged_SV_For_VNTR.pl
#    > Author: hwlu
#    > Mail: hongweilu@genetics.ac.cn 
#    Created Time: Thu 17 Dec 2020 05:29:51 PM CST
#########################################################################

#!/usr/bin/perl -w
use strict;

my $file = "Merge_Syri_bed.result.0.9.0.9.rename.tag.infer_1.infer_2.Add";

open(IN,"$file");
open(OUT,">$file.vntr");

while(<IN>){
    chomp;
    my $line = $_;
    my @line = split /\s+/,$line;
    my @spe = split /\|/,$line[10];
    $line[0]=~s/Undefined_//g;
    foreach my $spe (@spe){
	my @info = split /_/,$spe;
	print OUT "$info[0]\t$info[1]\t$info[2]\t$info[3]\t$info[4]\t$line[0]\t$line[1]\n";
    }

}
close IN;
close OUT;
