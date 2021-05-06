#########################################################################
#	File Name: formatTRF.pl
#	> Author: CaoYinghao
#	> Mail: yhcao@genetics.ac.cn 
#	Created Time: Tue May 26 16:06:47 2015
#########################################################################

#!/usr/bin/perl -w
use strict;
use Getopt::Long;

###trf genome.scount.fasta 2 7 7 80 10 50 500 -f -d -m
my @files = glob("*.txt.html");
my %gene;
my $c = 1;
foreach my $f(@files){
	#../573Mb-PBJelly.fasta.s1.2.7.7.80.10.50.500.1.txt.html
	open IN, $f;
	my ($name);
	while(<IN>){
		chomp;
		if(/^Sequence:\s(.*)/){
			$name = $1;
		}
		if(/Found at/){
			my ($match, $indel, $start,$end,$score,$psize,$copynum, $seq);
			while(<IN>){
				chomp;
				if(/^Pmatch=(.*),Pindel=(.*)/){
					$match = $1;
					$indel = $2;
				}elsif(/^\s+Indices:\s(\d+)--(\d+)\s+Score:\s(\d+)/){
					$start = $1;
					$end = $2;
					$score = $3;
				}elsif(/^\s+Period\ssize:\s(\d+).*Copynumber:\s(\S+)\s/){
					$psize = $1;
					$copynum = $2;
				}elsif(/^Consensus\spattern/){
					$seq = <IN>;
					$seq =~ s/\s+$//;
				}elsif(/^Right\sflanking/){
					push @{$gene{$name}{"str"}},$start."\t".$end."\t".$score."\t".$psize."\t".$copynum."\t".$seq;
					my $tname = $name;
					$tname =~ s/\D+//g;
					$gene{$name}{"name"} = $tname;
					last;
				}
			}
		}
	}
	close IN;
	#last;
}
open OUT, ">trf.gff";
print OUT "##gff-version 3\n";
foreach my $na(sort {$gene{$a}{"name"} <=> $gene{$b}{"name"}}keys %gene){
	my @array = @{$gene{$na}{"str"}};
	my %hash;
	foreach my $slice(@array){
		my @part = split "\t",$slice;
		$hash{$part[0]."\t".$part[1]}{"s"} = $part[0];
		$hash{$part[0]."\t".$part[1]}{"str"} = $slice;
	}
	foreach my $loc(sort {$hash{$a}{"s"} <=> $hash{$b}{"s"} || $a cmp $b} keys %hash){
		my @part = split "\t",$hash{$loc}{"str"};
		printf OUT "%s%06d%s", $na."\tTRF\tTandemRepeat\t".$part[0]."\t".$part[1]."\t".$part[2]."\t+\t.\tID=TR",$c,";PeriodSize=$part[3];CopyNumber=$part[4];Consensus=$part[5];\n";
		$c++;
	}
}

close OUT;

#my @files = glob("gff/*.gff");
#open OUT, ">trf.gff";
#my $count = 1;
#foreach my $p (@files){
#	my $name = $p;
#	#TRAP_file1_573Mb-PBJelly.fasta.s100023_Contig100022.gff
#	$name =~ s/.*fasta.*_(.*).gff/$1/;
#	#print $name."\n";
#	#Scaffold4811	TRF	TandemRepeat	1	481	761	+	.	ID=TR000008;PeriodSize=164;CopyNumber=2.9;PercentMatches=91;PercentIndels=1;Consensus=AAAACGGACCACTACTGATTTAGCCTGTATTCCAGCTTAAATCATATTTCTTCAAAAACGCAAACGCTTTCTCGACTACATCTAAAAATACGACACAGCTCATGGCAAAGCTATGACTATTGATCCTTGACCCGTTTCTGACGAAATTAGCCTACCGATTTGTCT;
#	#scaffold739     TRF     TandemRepeat    77220   77767   1096    +       .       ID=TR000024;PeriodSize=2;CopyNumber=274.0;PercentMatches=100;PercentIndels=0;Consensus=GA;
##.	TRF	satellite	44	86	100	+	.	note "satellite sequence" "TRF parameters   2 7 7 80 10 50 500" "repeat unit size = 11" "copy number = 3.9" "predicted by Tandem Repeats Finder 4.07b" ; label "satellite" ; rpt_type "tandem" ; rpt_unit "GGCCTTGTCAC" ; color 9
#	open IN,$p;
#	while(<IN>){
#		chomp;
#		my @d = split "\t";
#		my @part = split ";", $d[-1];
#		foreach my $slice (@part){
#			if($slice =~ /^note/){
#
#			}
#		}
#		print OUT $name."\tTRF\tTandemRepeat\t".$d[3]."\t".$d[4]."\t".$d[5]."\t".$d[6]."\t".$d[7]."\t";
#		$count++;
#	}
#}
