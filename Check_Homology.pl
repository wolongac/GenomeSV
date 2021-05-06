#########################################################################
#      File Name: Check_Homology.pl
#    > Author: hwlu
#    > Mail: hongweilu@genetics.ac.cn 
#    Created Time: Wed 15 Jul 2020 02:27:01 PM CST
#########################################################################

#!/usr/bin/perl -w
use strict;
use Bio::Tools::Run::Alignment::Blat;
use Bio::Seq;

my $sv=shift;
my @genomes=`ls *.genome`;
my $spe = $sv;
$spe =~ s/\..*//g;

my %hash_genome;
foreach my $genome (@genomes){
    chomp($genome);
    my $spe_genome = (split /\//,$genome)[-1];
    $spe_genome =~ s/\.genome//g;
    $spe_genome =~ s/v2//g;
    next if ($spe_genome ne $spe and $spe_genome ne "MSU");
    &read_genome($genome,$spe_genome);
}

open(IN,"$sv");
open(OUT,">$sv.homo");
while(<IN>){
    chomp;
    my $line = $_;
    my @line = split /\s+/,$line;
    my $homo1;
    my $homo2;
    my $ref_seq1;
    my $ref_seq2;
    my $query_seq1;
    my $query_seq2;
    if ($line[-1] =~ /small/){
	$homo1="-";
	$homo2="-";
    }else{
	if($line[2]-$line[1]>=500){
		$ref_seq1=substr($hash_genome{'MSU'}{$line[0]},$line[1]-250,500);
		$ref_seq2=substr($hash_genome{'MSU'}{$line[0]},$line[2]-250,500);
	}else{
		my $length_tmp = $line[2]-$line[1]+1;
		my $cut = int($length_tmp/2);
		my $cut_length=250+$cut;
		$ref_seq1=substr($hash_genome{'MSU'}{$line[0]},$line[1]-250,$cut_length);
		$ref_seq2=substr($hash_genome{'MSU'}{$line[0]},$line[2]-$cut,$cut_length);
	}
	if($line[7]-$line[6]>500){
		$query_seq1=substr($hash_genome{$spe}{$line[5]},$line[6]-250,500);
		$query_seq2=substr($hash_genome{$spe}{$line[5]},$line[7]-250,500);		
	}else{
		my $length_tmp = $line[7]-$line[6]+1;
		my $cut = int($length_tmp/2);
		my $cut_length=250+$cut;
		$query_seq1=substr($hash_genome{$spe}{$line[5]},$line[6]-250,$cut_length);
		$query_seq2=substr($hash_genome{$spe}{$line[5]},$line[7]-$cut,$cut_length);
	}
	if($line[2]-$line[1]>=50){
	   $homo1=&Blat($ref_seq1,$ref_seq2);
	}else{
	    $homo1=0;
	}
	if($line[7]-$line[6]>=50){
	   $homo2=&Blat($query_seq1,$query_seq2);
	}else{
	    $homo2=0;
	}
	if($line[5] eq "-"){
	    $homo2=0;
	}
    }
    print OUT "$line[0]\t$line[1]\t$line[2]\t$line[3]\t$line[4]\t$homo1\t$line[5]\t$line[6]\t$line[7]\t$line[8]\t$line[9]\t$homo2\t$line[10]\n";
}
close OUT;


sub Blat{
    my $fasta1=shift;
    my $fasta2=shift;
    my $seq1=Bio::Seq->new(-seq => $fasta1,-alphabet => 'dna',-id=>"seq1");
    my $seq2=Bio::Seq->new(-seq => $fasta2,-alphabet => 'dna',-id=>"seq2");
    my $out=Bio::SeqIO->new(-file => ">$spe.db.blat",-format => 'Fasta');
    $out->write_seq($seq1);
    my $blat_db=Bio::Tools::Run::Alignment::Blat->new(-db => "$spe.db.blat", -t => 'dna', -q => 'dna');
    $blat_db->outfile_name("$spe.out.blat");
    my $blat_result=$blat_db->run($seq2);
    #my $match=`cat $spe.out.blat |grep -P "^\\d"|awk '{print \$1}'|sort -k 1,1nr|head -n 1`;
	my $match=0;
    open(IN1,"$spe.out.blat");
    while(<IN1>){
		chomp;
		next if ($_ !~ /^\d/);
		my @line = split /\s+/,$_;
		#next if ($line[8] eq "+");
		next if ($line[11]>251 or  $line[12]<249);
		#next if ($line[11]>250 or  $line[12]<250);
		next if ($line[15]>251 or  $line[16]<249);
		#next if ($line[15]>250 or  $line[16]>250);
		if($match < $line[0]){
			$match = $line[0];
		}
    }
	close IN1;
    if($match==0){
	$match=0;
    }
    return $match;
}


sub read_genome{
    my $file = shift;
    my $spe_genome=shift;
    open(IN,"$file");
    my $name;
    while(<IN>){
	chomp;
	if($_ =~ />/){
	    $name = $_;
	    $name =~s/>//;
	    $name =~ s/\s.*//g;
	}else{
	    $hash_genome{$spe_genome}{$name}.=$_;
	}
    }
    close IN;
}
