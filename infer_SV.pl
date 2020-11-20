#########################################################################
#      File Name: infer_SV.pl
#    > Author: hwlu
#    > Mail: hongweilu@genetics.ac.cn 
#    Created Time: Wed 15 Jul 2020 12:12:26 AM CST
#########################################################################

#!/usr/bin/perl -w
use strict;

my $sv_merged = shift;

open(IN,"$sv_merged");
while(<IN>){
    chomp;
    my @line = split /\s+/,$_;
    my @sv = split /\|/,$line[-1];
    foreach my $sv (@sv){
	my @tmp = split /\|/,$sv;
	my $new_type = "$line[0]";
	if($line[-1] =~ /CJ14/){
	    $new_type ="Infered:$line[0]";
	    if($line[0] eq "DEL"){
		$new_type = "Infered:INS";
	    }elsif($line[0] eq "INS"){
		$new_type = "Infered:DEL";
	    }
	}
	$sv =~ s/_/\t/g;
	print "$sv\t$new_type\n";
    }
}
close IN;
