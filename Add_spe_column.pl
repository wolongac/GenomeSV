#########################################################################
#      File Name: tmp.pl
#    > Author: hwlu
#    > Mail: hongweilu@genetics.ac.cn 
#    Created Time: Tue 15 Dec 2020 09:45:08 AM CST
#########################################################################

#!/usr/bin/perl -w
use strict;

my $file = shift;
open(IN,"$file");
open(OUT,">$file.tag");

while(<IN>){
	chomp;
    my $line = $_;
    my @line = split /\s+/,$line;
    my @spe = split /\|/,$line[-1];
    my @spe_sim;
    foreach my $spe (@spe){
	my $spe_new = (split /_/,$spe)[-1];
	push(@spe_sim,$spe_new);
    }
    my $spe = join "|",@spe_sim;
    print OUT "$line\t$spe\n";
}
