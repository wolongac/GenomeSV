#########################################################################
#      File Name: Infer_2.pl
#    > Author: hwlu
#    > Mail: hongweilu@genetics.ac.cn 
#    Created Time: Tue 15 Dec 2020 12:05:07 PM CST
#########################################################################

#!/usr/bin/perl -w
use strict;

my $file = "Merge_Syri_bed.result.0.9.0.9.tag.infer1";
my $log = "Infer_2.log";
my %hash;

open(IN,"$log");
while(<IN>){
    chomp;
    my @line = split /\s+/,$_;
    next if ($line[1] eq "Keep");
    $hash{$line[0]}=$line[1];
    $hash{$line[2]}=$line[1];
}
close IN;

open(IN,"$file");
open(OUT,">$file.infer_2");
while(<IN>){
    chomp;
    my $line = $_;
    my @line = split /\s+/,$line;
    if($line[0] eq "Translocation"){
	    if(exists $hash{$line[1]}){
		$line[0] = "Undefined_Translocation";
	    }
	$line = join "\t",@line;
	print OUT "$line\n";
	next;
    }
    if($line[0] =~ /Infered/ or !exists $hash{$line[1]}){
	print OUT "$line\n";
    }else{
	if($hash{$line[1]} eq "Abandon"){
	    if($line[0] !~ /_/){
		$line[0]="Undefined_$line[0]";
	    }
	}else{
	    $line[0]=(split /_/,$line[0])[-1];
	    if($line[0] eq "INS"){
		$line[0]="Infered_DEL";
	    }else{
		$line[0] = "Infered_$line[0]";
	    }
	}
	$line = join "\t",@line;
	print OUT "$line\n";
    }
}
close OUT;
close IN;
