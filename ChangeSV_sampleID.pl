#########################################################################
#      File Name: ChangeSV_sampleID.pl
#    > Author: hwlu
#    > Mail: hongweilu@genetics.ac.cn 
#    Created Time: Sat 12 Dec 2020 06:53:35 PM CST
#########################################################################

#!/usr/bin/perl -w
use strict;

my %hash_genome;
open(IN,"genome.list");

while(<IN>){
    chomp;
    $hash_genome{$_}=1;
}
close IN;

open(IN,"Merge_Syri_bed.result.0.9.0.9.tag.infer_1.infer_2");
open(OUT,">Merge_Syri_bed.result.0.9.0.9.tag.infer_1.infer_2.Add");
while(<IN>){
    chomp;
    my $line = $_;
    my @line = split /\s+/,$line;
    my @spe_tmp = split /\|/,$line[10];
    my @spe_tmp2;
    foreach my $spe_tmp (@spe_tmp){
	$spe_tmp=(split /_/,$spe_tmp)[-1];
	push(@spe_tmp2,$spe_tmp);

    }
    my $spe_tmp = join "|",@spe_tmp2;
    if($line[0] !~/_/ ){
	push(@line,$spe_tmp);
	my $line_new1=join "\t",@line;
	print OUT "$line_new1\n";
	next;
    }
    if($line[0] =~ /Undefined/){
	push(@line,"-");
	my $line_new1=join "\t",@line;
	print OUT "$line_new1\n";
	next;
    }
    my @spe=split /\|/,$spe_tmp;
    my %hash_tmp=();
    foreach my $spe (@spe){
	$hash_tmp{$spe}=1;
    }
    my @spe_new;
    foreach my $spe (keys %hash_genome){
	next if (exists $hash_tmp{$spe});
	push(@spe_new,$spe);
    }
    my $spe = join "|",@spe_new;
    push(@line,$spe);
    my $line_new=join "\t",@line;
    print OUT "$line_new\n";
    
}
close IN;
close OUT;
