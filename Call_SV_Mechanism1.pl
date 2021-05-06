#########################################################################
#      File Name: Call_SV_Mechanism1.pl
#    > Author: hwlu
#    > Mail: hongweilu@genetics.ac.cn 
#    Created Time: Tue 14 Jul 2020 11:31:32 PM CST
#########################################################################

#!/usr/bin/perl -w
use strict;

my $sv=shift;
my $spe_input=shift;
my $distance=1;
my $num=2;
my @trf=`ls /public-dss/Project/Project_hwlu/hwlu/Core_set/TRF/result/*.trf`;
my @TE_gff=`ls /public-dss/Project/Project_hwlu/hwlu/Core_set/genomes/TE_repeat/*.genome.out.gff`;

my %hash_trf;
my %hash_TE;   #$hash_TE{$spe}{$chr}{$start}=$end;
my %hash_sv_ref;	#$hash_sv{$}
my %hash_sv_query;	#$hash_sv{$}

print "Reading TRF...\n";
foreach my $trf (@trf){
	chomp($trf);
	#&read_trf($trf);
}
print "Done!!!\n";

print "Reading TE result...\n";
foreach my $TE_gff (@TE_gff){
	chomp($TE_gff);
	#$read_TE($TE_gff);
}
print "Done!!!\n";

print "Reading SV...\n";
open(IN,"$sv");
while(<IN>){
	chomp;
	my $line = $_;
	my @line = split /\s+/,$line;
	$hash_sv_ref{$line[0]}{$line[1]}{$line[2]}{'type'}=$line[-1];
	$hash_sv_ref{$line[0]}{$line[1]}{$line[2]}{'pos'}="$line[3]\t$line[4]\t$line[5]";
	$hash_sv_ref{$line[0]}{$line[1]}{$line[2]}{'info'}=$line;
	next if ($line[3] eq "-");
	$hash_sv_query{$line[3]}{$line[4]}{$line[5]}{'type'}=$line[-1];
	$hash_sv_query{$line[3]}{$line[4]}{$line[5]}{'pos'}="$line[0]\t$line[1]\t$line[2]";
	$hash_sv_query{$line[3]}{$line[4]}{$line[5]}{'info'}=$line;
}
close IN;

# Get SV Cluster: by REF;
my %hash_sv_ref_array;
foreach my $chr (keys %hash_sv_ref){
	foreach my $start (sort {$a <=> $b} keys %{$hash_sv_ref{$chr}}){
		foreach my $end (sort {$a <=> $b} keys %{$hash_sv_ref{$chr}{$start}}){
			push(@{$hash_sv_ref_array{$chr}},[$chr,$start,$end]);
		}
	}
}
my $cluster_ID=1;
my %hash_sv_ref_Cluster_sv;
my %hash_sv_ref_Cluster;
foreach my $chr (keys %hash_sv_ref_array){
	my @svs=@{$hash_sv_ref_array{$chr}};
	my $sv_tmp_id=join "\t",@{$svs[0]};
	$hash_sv_ref_Cluster_sv{$sv_tmp_id}=$cluster_ID;
	push(@{$hash_sv_ref_Cluster{$cluster_ID}},$sv_tmp_id);
	for(my $i=1;$i<@svs;$i++){
		if($svs[$i][1] - $svs[$i-1][2] < $distance ){
			$sv_tmp_id=join "\t",@{$svs[$i]};
			#print "$sv_tmp_id\n";
			$hash_sv_ref_Cluster_sv{$sv_tmp_id}=$cluster_ID;
			push(@{$hash_sv_ref_Cluster{$cluster_ID}},$sv_tmp_id);
		}else{
			$cluster_ID++;
			$sv_tmp_id=join "\t",@{$svs[$i]};
			#print "$sv_tmp_id\n";
			$hash_sv_ref_Cluster_sv{$sv_tmp_id}=$cluster_ID;
			push(@{$hash_sv_ref_Cluster{$cluster_ID}},$sv_tmp_id);
		}
	}
}

open(OUT1,">$sv.ref.Cluster_Type");
foreach my $Cluster (sort {$a <=> $b} keys %hash_sv_ref_Cluster){
	my @svs = @{$hash_sv_ref_Cluster{$Cluster}};
	my $sv_num=@svs;
	my $Cluster_type="solo";
	if(@svs>$num){
		$Cluster_type="Cluster";
	}elsif(@svs==1){
		$Cluster_type="solo";
	}else{
		$Cluster_type="coupling";
	}
	print OUT1 "$Cluster\t$Cluster_type\t$sv_num";
	foreach my $sv (@svs){
		my @sv_info = split /\t/,$sv;
		my $sv_type = $hash_sv_ref{$sv_info[0]}{$sv_info[1]}{$sv_info[2]}{'type'};
		print OUT1 "\t$sv_info[0]:$sv_info[1]:$sv_info[2]:$sv_type";
	}
	print OUT1 "\n";
}
close OUT1;

# Get SV Cluster: by Query;
my %hash_sv_query_array;
foreach my $chr (keys %hash_sv_query){
	foreach my $start (sort {$a <=> $b} keys %{$hash_sv_query{$chr}}){
		foreach my $end (sort {$a <=> $b} keys %{$hash_sv_query{$chr}{$start}}){
			push(@{$hash_sv_query_array{$chr}},[$chr,$start,$end]);
		}
	}
}
my $cluster_ID=1;
my %hash_sv_query_Cluster_sv;
my %hash_sv_query_Cluster;
foreach my $chr (keys %hash_sv_query_array){
	my @svs=@{$hash_sv_query_array{$chr}};
	my $sv_tmp_id=join "\t",@{$svs[0]};
	$hash_sv_query_Cluster_sv{$sv_tmp_id}=$cluster_ID;
	push(@{$hash_sv_query_Cluster{$cluster_ID}},$sv_tmp_id);
	for(my $i=1;$i<@svs;$i++){
		if($svs[$i][1] - $svs[$i-1][2] < $distance ){
			$sv_tmp_id=join "\t",@{$svs[$i]};
			#print "$sv_tmp_id\n";
			$hash_sv_query_Cluster_sv{$sv_tmp_id}=$cluster_ID;
			push(@{$hash_sv_query_Cluster{$cluster_ID}},$sv_tmp_id);
		}else{
			$cluster_ID++;
			$sv_tmp_id=join "\t",@{$svs[$i]};
			#print "$sv_tmp_id\n";
			$hash_sv_query_Cluster_sv{$sv_tmp_id}=$cluster_ID;
			push(@{$hash_sv_query_Cluster{$cluster_ID}},$sv_tmp_id);
		}
	}
}

open(OUT1,">$sv.query.Cluster_Type");
foreach my $Cluster (sort {$a <=> $b} keys %hash_sv_query_Cluster){
	my @svs = @{$hash_sv_query_Cluster{$Cluster}};
	my $sv_num=@svs;
	my $Cluster_type="solo";
	if(@svs>$num){
		$Cluster_type="Cluster";
	}elsif(@svs==1){
		$Cluster_type="solo";
	}else{
		$Cluster_type="coupling";
	}
	print OUT1 "$Cluster\t$Cluster_type\t$sv_num";
	foreach my $sv (@svs){
		my @sv_info = split /\t/,$sv;
		my $sv_type = $hash_sv_query{$sv_info[0]}{$sv_info[1]}{$sv_info[2]}{'type'};
		print OUT1 "\t$sv_info[0]:$sv_info[1]:$sv_info[2]:$sv_type";
	}
	print OUT1 "\n";
}
close OUT1;




###Type:
###DEL
###Infered:DEL
###Infered:INS
###Infered:Inversion
###Infered:Translocation
###INS
###Inversion
###small:DEL
###small:INS
###Translocation

###Inversion;
#foreach my $chr (keys %hash_sv){
#	foreach my $start (sort {$a <=> $b} keys %{$hash_sv{$chr}}){
#		foreach my $end (sort ($a <=> $b) keys %{$hash_sv}{$chr}{$start}){
#			my $type = $hash_sv{$chr}{$start}{$end}{'type'};
#			my $pos = $hash_sv{$chr}{$start}{$end}{'pos'};
#			next if ($type !~ "Inversion");
#			
#			
#		}
#	}
#}





sub read_TE{
	my $file = shift;
	my $spe = (split /\//,$file)[-1];
	my $spe =~s/\..*//g;
	next if ($spe ne $spe_input and $spe ne "MSU");
	print "$spe\t$file\n";
	open(IN,"$file");
	while(<IN>){
		chomp;
		my $line = $_;
		next if ($line =~/^#/);
		my @line = split /\s+/,$line;
		next if (@line < 9);
		$hash_TE{$spe}{$line[0]}{$line[3]}={$line[4]};
	}
	close IN;
}




sub read_trf{
    my $file = shift;
	my $spe = (split /\//,$file)[-1];
	$spe =~ s/\.trf//g;
	next if ($spe ne $spe_input and $spe ne "MSU");
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
			$hash_trf{$spe}{$line[0]}{$line[3]}=$line[4];
	    }
	    close IN;

}
