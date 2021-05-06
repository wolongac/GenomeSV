#########################################################################
#      File Name: Filter_overlaped_SV.pl
#    > Author: hwlu
#    > Mail: hongweilu@genetics.ac.cn 
#    Created Time: Sun 09 Aug 2020 12:27:02 AM CST
#########################################################################

#!/usr/bin/perl -w
use strict;

my $file = shift;
open(IN,"$file");
open(OUT,">$file.uniq");
my %hash_count;
my %hash_info;
while(<IN>){
    chomp;
    my $line = $_;
    my @line = split /\s+/,$_;
    my $pos="$line[0]\t$line[1]\t$line[2]";
    $hash_info{$pos}=$line;
    $hash_count{$pos}++;
}
close IN;

foreach my $pos (keys %hash_info){
    next if ($hash_count{$pos}>1);
    print OUT "$hash_info{$pos}\n";
}
close OUT;
