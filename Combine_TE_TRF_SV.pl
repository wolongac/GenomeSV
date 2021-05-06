#########################################################################
#      File Name: Combine_TE_TRF_SV.pl
#    > Author: hwlu
#    > Mail: hongweilu@genetics.ac.cn 
#    Created Time: Wed 15 Jul 2020 03:02:35 PM CST
#########################################################################

#!/usr/bin/perl -w
use strict;

my $sv=shift;
my @trf=`ls *.trf`;
#my @TE_gff=`ls /public-dss/Project/Project_hwlu/hwlu/Core_set/genomes/TE_repeat/*.genome.out.gff`;
my @TE_gff=`ls *.genome.out.gff`;



print "$sv\n";
my $spe = $sv; 	
$spe =~ s/\..*//g;


my %hash_trf;
my %hash_TE;

foreach my $trf (@trf){
	chomp($trf);
	&read_trf($trf);
}

foreach my $te (@TE_gff){
	chomp($te);
	&read_TE($te);
}



open(OUT,">$sv.TRF_TE");
open(IN,$sv);
while(<IN>){
	chomp;
	my @line = split /\s+/,$_;
	if($line[1]>$line[2]){
		($line[2],$line[1]) = ($line[1],$line[2]);
	}
	
	if($line[4] > $line[5]){
		($line[4],$line[5]) = ($line[5],$line[4]);
	}
	my $trf_ref  = check_trf("MSU",$line[0],$line[1],$line[2]);
	my $TE_ref   = check_TE("MSU",$line[0],$line[1],$line[2]);
	my $trf_query= check_trf($spe,$line[3],$line[4],$line[5]);
	my $TE_query = check_TE($spe,$line[3],$line[4],$line[5]);
	print OUT "$line[0]\t$line[1]\t$line[2]\t$trf_ref\t$TE_ref\t$line[3]\t$line[4]\t$line[5]\t$trf_query\t$TE_query\t$line[6]\n";
}
close OUT;

sub check_trf{
	my $spe=shift;
	my $chr=shift;
	my $start=shift;
	my $end = shift;
	my $overlap_last=0;
	foreach my $start_tmp (sort {$a <=> $b} keys %{$hash_trf{$spe}{$chr}}){
		my $end_tmp = $hash_trf{$spe}{$chr}{$start_tmp};
		my $overlap=0;
		next if ($end_tmp < $start);
		last if ($start_tmp > $end);
		if($start_tmp >= $start){
			if($end_tmp > $end){
				$overlap = $end - $start_tmp + 1;
			}else{
				$overlap = $end_tmp - $start_tmp + 1;
			}
		}else{
			if($end_tmp > $end){
				$overlap = $end - $start + 1;
			}else{
				$overlap = $end_tmp - $start + 1;
			}
		}
		if($overlap > $overlap_last){
			$overlap_last = $overlap;
		}
	}
	my $ratio = $overlap_last/($end-$start+1);
	$ratio = sprintf("%.2f", $ratio);
	return "VNTR:$ratio";
}

sub check_TE{
	my $spe=shift;
	my $chr=shift;
	my $start=shift;
	my $end = shift;
	my $overlap_last=0;
	foreach my $start_tmp (sort {$a <=> $b} keys %{$hash_TE{$spe}{$chr}}){
		my $end_tmp = $hash_TE{$spe}{$chr}{$start_tmp};
		my $overlap=0;
		next if ($end_tmp < $start);
		last if ($start_tmp > $end);
		if($start_tmp >= $start){
			if($end_tmp > $end){
				$overlap = $end - $start_tmp + 1;
			}else{
				$overlap = $end_tmp - $start_tmp + 1;
			}
		}else{
			if($end_tmp > $end){
				$overlap = $end - $start + 1;
			}else{
				$overlap = $end_tmp - $start + 1;
			}
		}
		if($overlap > $overlap_last){
			$overlap_last = $overlap;
			#print "$spe\t$chr\t$start\t$end\t$start_tmp\t$end_tmp\n";
		}
	}
	my $ratio = $overlap_last/($end-$start+1);
	$ratio = sprintf("%.2f", $ratio);
	return "TE:$ratio";
}


sub read_trf{
    my $file = shift;
	my $spe_tmp = (split /\//,$file)[-1];
	$spe_tmp =~ s/\.trf//g;
	$spe_tmp =~ s/v2//g;
	next if ($spe_tmp ne $spe and $spe_tmp ne "MSU");
	print "$spe\t$file\n";
	open(IN,"$file");
	while(<IN>){
	    chomp;
	    next if ($_ =~ /^#/);
	    my $line = $_;
	    my @line = split /\s+/,$line;
	    next if (@line != 9);
	    my ($copy_num,$seq)=$line[-1]=~/CopyNumber=(.*?);Consensus=(.*?);/;
	    my $unit_length=length($seq);
	    if($line[3]>$line[4]){
			($line[4],$line[3]) = ($line[3],$line[4]);
	    }
	    my $total_length=$line[4]-$line[3]+1;
	    next if ($copy_num<10 or $copy_num > 1500);
	    next if ($unit_length < 10 or $unit_length > 100);
	    next if ($total_length < 500 or $total_length > 15000);
			$hash_trf{$spe_tmp}{$line[0]}{$line[3]}=$line[4];
	    }
	    close IN;

}

sub read_TE{
	my $file = shift;
	my $spe_tmp = (split /\//,$file)[-1];
	print "$file\t$spe\n";
	$spe_tmp =~ s/\..*//g;
	$spe_tmp =~ s/v2//g;
	next if ($spe_tmp ne $spe and $spe_tmp ne "MSU");
	print "$spe\t$file\n";
	open(IN,"$file");
	while(<IN>){
		chomp;
		my $line = $_;
		next if ($line =~/^#/);
		my @line = split /\s+/,$line;
		next if (@line < 9);
		$hash_TE{$spe_tmp}{$line[0]}{$line[3]}=$line[4];
	}
	close IN;
}

