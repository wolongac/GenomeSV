#########################################################################
#      File Name: Infer_1.pl
#    > Author: hwlu
#    > Mail: hongweilu@genetics.ac.cn 
#    Created Time: Tue 15 Dec 2020 10:50:33 AM CST
#########################################################################

#!/usr/bin/perl -w
use strict;

my $file = "Merge_Syri_bed.result.0.9.0.9.tag";

open(IN,"$file");
open(OUT,">$file.infer1");
while(<IN>){
    chomp;
    my $line = $_;
    my @line = split /\s+/,$line;
    if($line =~ /CJ14/ and $line[-1] ne "CJ14"){
	    if($line[0] eq "INS"){
		$line[0]="Infered_DEL";
	    }elsif($line[0] eq "DEL"){
		$line[0] = "Infered_INS";
	    }else{
		$line[0] = "Infered_$line[0]";
	    }
    }elsif($line[-1] eq "CJ14"){
	$line[0]="Undefined_$line[0]";
    }
    $line = join "\t",@line;
    print OUT  "$line\n";
}
close IN;
close OUT;
