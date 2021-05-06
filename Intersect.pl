#########################################################################
#      File Name: Intersect.pl
#    > Author: hwlu
#    > Mail: hongweilu@genetics.ac.cn 
#    Created Time: Tue 15 Dec 2020 09:58:49 AM CST
#########################################################################

#!/usr/bin/perl -w
use strict;

my $file = "Merge_Syri_bed.result.0.9.0.9.tag.infer1";

my %hash_c;
my %hash_o;

open(IN,"$file");
while(<IN>){
    chomp;
    my @line = split /\s+/,$_;
    if($line[-1] eq "CJ14"){
	$line[0] = (split /_/,$line[0])[-1];
	$hash_c{$line[0]}{$line[2]}{$line[3]}{$line[4]}=$line[1];
    }else{
	next if ($line[0] =~ /Infered_/);
	$hash_o{$line[0]}{$line[2]}{$line[3]}{$line[4]}=$line[1];
    }
}
close IN;

open(OUT,">Infer_2.log");
foreach my $type (keys %hash_o){
    foreach my $chr (keys %{$hash_o{$type}}){
	foreach my $s (keys %{$hash_o{$type}{$chr}}){
	    foreach my $e (keys %{$hash_o{$type}{$chr}{$s}}){
		my $SV1=$hash_o{$type}{$chr}{$s}{$e};
		my $tag="Keep"; #no change
		if($type eq "DEL"){
		    foreach my $s1 (sort {$a <=> $b} keys %{$hash_c{$type}{$chr}}){
			last if ($s1 > $s);
			foreach my $e1 (keys %{$hash_c{$type}{$chr}{$s1}}){
			    next if ($e1 < $e);
			    $tag="Abandon";   #contain by CG14_only; abandon
			    my $SV2 = $hash_c{$type}{$chr}{$s1}{$e1};
			    print OUT "$SV1\t$tag\t$SV2\n";
			}
		    }
		    if($tag eq "Keep"){
			print OUT "$SV1\t$tag\t-\n";
		    }
		}elsif($type eq "INS"){
		    if(exists $hash_c{"INS"}{$chr}{$s}{$e}){
			$tag="Infer";  #covery by CH14_only ; infer;
			my $SV2=$hash_c{$type}{$chr}{$s}{$e};
			print OUT "$SV1\t$tag\t$SV2\n";
		    }else{
			print OUT "$SV1\t$tag\t-\n";
		    }
		}elsif($type eq "Inversion"){
		    foreach my $s1 (sort {$a <=> $b} keys %{$hash_c{"DEL"}{$chr}}){
			last if ($s1>$s);
			foreach my $e1 (keys %{$hash_c{'DEL'}{$chr}{$s1}}){
			    next if ($e1<$e);
			    $tag="Abandon"; #contained by DEL, cannot infer the type 
			    my $SV2=$hash_c{'DEL'}{$chr}{$s1}{$e1};
			    print OUT "$SV1\t$tag\t$SV2\n";
			}
		    }
		    last if ($tag eq "Abandon");
		    foreach my $s1 (sort {$a <=> $b} keys %{$hash_c{$type}{$chr}}){
			last if ($s1 > $e);
			foreach my $e1 (sort keys %{$hash_c{$type}{$chr}{$s1}}){
			    next if ($e1 < $s);
			    $tag="Infer"; #covery by CG14_only, infer;
			    my $SV2=$hash_c{$type}{$chr}{$s1}{$e1};
			    print OUT "$SV1\t$tag\t$SV2\n";
			}
		    }
		    if($tag eq "Keep"){
			print OUT "$SV1\t$tag\t-\n";
		    }
		}else{     #Translocation
		    foreach my $s1 (sort {$a <=> $b} keys %{$hash_c{"DEL"}{$chr}}){
			last if ($s1>$s);
			foreach my $e1 (keys %{$hash_c{'DEL'}{$chr}{$s1}}){
			    next if ($e1<$e);
			    $tag="Abandon"; #contained by DEL of CG14_only, cannot infer the type 
			    my $SV2=$hash_c{'DEL'}{$chr}{$s1}{$e1};
			    print OUT "$SV1\t$tag\t$SV2\tTranslocation\n";
			}
		    }
		}
	    }
	}
    }
}
close OUT;
