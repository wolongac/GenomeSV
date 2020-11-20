#########################################################################
#      File Name: GetInsertionSeq.pl
#    > Author: hwlu
#    > Mail: hongweilu@genetics.ac.cn 
#    Created Time: Mon 23 Mar 2020 07:50:33 PM CST
#########################################################################

#!/usr/bin/perl -w
use strict;

#my @bed = `ls *.syri.out.SV.bed`;
my @bed = `ls *.syri.out.SV.bed.mergeDel.CPV2InDel.50.bp`;

my $ins_distance=0;
my $overlap=0.9;
my $coverage=0.9;

my @spe;
my %hash_SV;
my %length_SV;
print "Reading SV bed files...\n";
foreach my $bed (@bed){
    chomp($bed);
    my $spe = $bed;
    #$spe =~ s/\.syri.out.SV.bed//g;
    $spe =~ s/\.syri\.out\.SV\.bed\.mergeDel\.CPV2InDel\.50\.bp//g;
    print "\t".$spe."\n";
    push(@spe,$spe);
    
    open(IN,"$bed");
    while(<IN>){
	chomp;
	my $line=$_;
	my @line = split /\s+/,$line;
	if($line[4]>$line[5]){
	    ($line[4],$line[5])=($line[5],$line[4]);
	}
	if($line[1]>$line[2]){
	    ($line[1],$line[2]) = ($line[2],$line[1]);
	}
	if($line[6] eq "Inversion"){
	    $hash_SV{$line[6]}{$line[0]}{$line[1]}{$line[2]}{$spe}="$spe\_$line[3]\_$line[4]\_$line[5]";
	    $length_SV{"$line[6]\_$line[0]\_$line[1]\_$line[2]\_$spe"}=abs($line[2]-$line[1]+1);
	}elsif($line[6] eq "DEL" or $line[6] eq "Deletion"){
	    $hash_SV{"DEL"}{$line[0]}{$line[1]}{$line[2]}{$spe}=1;
	    $length_SV{"DEL\_$line[0]\_$line[1]\_$line[2]\_$spe"}=abs($line[2]-$line[1]+1);
	}elsif($line[6] eq "INS" or $line[6] eq "Insertion" ){
	    $hash_SV{"INS"}{$line[0]}{$line[1]}{$line[2]}{$spe}="$spe\_$line[3]\_$line[4]\_$line[5]";
	    $length_SV{"INS\_$line[0]\_$line[1]\_$line[2]\_$spe"}=abs($line[5]-$line[4]+1);
	}elsif($line[6] =~/Translocation/ or $line[6] =~ /TR/){
	    $hash_SV{$line[6]}{$line[0]}{$line[1]}{$line[2]}{$spe}="$spe\_$line[3]\_$line[4]\_$line[5]";
	    $length_SV{"$line[6]\_$line[0]\_$line[1]\_$line[2]\_$spe"}=abs($line[2]-$line[1]);
	}else{
		next;
	}
    }
    close IN;
}
print "Done...\n";

print "Format data ...\n";
my %hash_SV_array;
foreach my $type (keys %hash_SV){
    foreach my $chr (keys %{$hash_SV{$type}}){
	foreach my $start (sort {$a <=> $b} keys %{$hash_SV{$type}{$chr}}){
	    foreach my $end (sort {$a <=> $b} keys %{$hash_SV{$type}{$chr}{$start}}){
		    foreach my $spe (keys %{$hash_SV{$type}{$chr}{$start}{$end}}){
			push(@{$hash_SV_array{$type}{$chr}},"$start\t$end\t$spe");
		    }
	    }
	}
    }
}
print "Done...\n";

print "Merge SV ... \n";
my %hash_merged_result;
foreach my $type (keys %hash_SV_array){
    print "\t$type...\n";
    my $SV_ID=1;
    foreach my $chr (keys %{$hash_SV_array{$type}}){
	my @tmp_chr = @{$hash_SV_array{$type}{$chr}};
	my $num_sv=@tmp_chr;
	print "\t\t$chr contains $num_sv  SVs ...\n";
	for(my $i=0;$i<@tmp_chr;$i++){
	    my ($start_1,$end_1,$spe_1)=split /\t/,$tmp_chr[$i];
	    next if (exists $hash_merged_result{$type}{$chr}{$start_1}{$end_1}{$spe_1});
	    $SV_ID++;
	    $hash_merged_result{$type}{$chr}{$start_1}{$end_1}{$spe_1}=$SV_ID;
	    for(my $j=$i+1;$j<@tmp_chr;$j++){
		my ($start_2,$end_2,$spe_2)=split /\t/,$tmp_chr[$j];
		next if (exists $hash_merged_result{$type}{$chr}{$start_2}{$end_2}{$spe_2});
		last if ($start_2>$end_1);
		if($type eq "INS" or $type eq "Insertion"){     #ref overlap and query RO >0.5
			my $ins_1=$hash_SV{$type}{$chr}{$start_1}{$end_1}{$spe_1};
			my $ins_2=$hash_SV{$type}{$chr}{$start_2}{$end_2}{$spe_2};
			my @ins_1=split /_/,$ins_1;
			my @ins_2=split /_/,$ins_2;
			my $length1=$ins_1[3]-$ins_1[2]+1;
			my $length2=$ins_2[3]-$ins_2[2]+1;
			    if($length1>$length2){
				($length1,$length2)=($length2,$length1);
			    }
			    next if($length1/$length2<$coverage);

		}elsif($type eq "DEL" or $type eq "INV" or $type eq "TRANS" or $type eq "INVTR" or $type eq "Inversion" or $type eq "Translocation"){  #ref RO >0.5
			my $RO=0;
			if($end_2 > $end_1) {
			    $RO=($end_1-$start_2+1)/($end_2-$start_1+1);
			}else{
			    $RO=($end_2-$start_2+1)/($end_1-$start_1+1);
			}
			next if ($RO<$overlap);
		}else{    #HDR,CPG,CPL,TDM : ref RO >0.5 and query RO > 0.5
			my $RO=0;
			if($end_2>$end_1){
			    $RO=($end_1-$start_2+1)/($end_2-$start_1+1);
			}else{
			    $RO=($end_2-$start_2+1)/($end_1-$start_1+1);
			}
			next if ($RO<$overlap);
			my $ins_1=$hash_SV{$type}{$chr}{$start_1}{$end_1}{$spe_1};
			my $ins_2=$hash_SV{$type}{$chr}{$start_2}{$end_2}{$spe_2};
			my @ins_1=split /_/,$ins_1;
			my @ins_2=split /_/,$ins_2;
			my $length1=$ins_1[3]-$ins_1[2]+1;
			my $length2=$ins_2[3]-$ins_2[2]+1;
			if($length1>$length2){
			    ($length1,$length2)=($length2,$length1);
			}
			next if($length1/$length2<$coverage);
		}
		$hash_merged_result{$type}{$chr}{$start_2}{$end_2}{$spe_2}=$SV_ID;
	    }
	}
	    print "\t\t\tNow $type\tSV_ID is :\t$SV_ID\n";
    }
}
print "Done...\n";
my %result;
foreach my $type (keys %hash_merged_result){
    foreach my $chr (sort keys %{$hash_merged_result{$type}}){
	foreach my $start (sort {$a <=> $b} keys %{$hash_merged_result{$type}{$chr}}){
	    foreach my $end (sort {$b <=> $a} keys %{$hash_merged_result{$type}{$chr}{$start}}){
		foreach my $spe (keys %{$hash_merged_result{$type}{$chr}{$start}{$end}}){
		    my $SV_ID=$hash_merged_result{$type}{$chr}{$start}{$end}{$spe};
		    push(@{$result{$type}{$SV_ID}},"$type\_$chr\_$start\_$end\_$spe");
		}
	    }
	}
    }
}

print "OUTPUT...\n";
open(OUT,">Merge_Syri_bed.result.$overlap.$coverage.v2");
foreach my $type (keys %result){
    foreach my $ID (sort {$a <=> $b} keys %{$result{$type}}){
	my $present_sv;
	foreach my $sv (sort {$length_SV{$b} <=> $length_SV{$a}} @{$result{$type}{$ID}}){
		$present_sv=$sv;
		last;
	}
	foreach my $sv (sort {$length_SV{$b} <=> $length_SV{$a}} @{$result{$type}{$ID}}){
		next if ($sv !~ /R527/);
		$present_sv=$sv;
		last;
	}
	my @present_sv=split /_/,$present_sv;
	if($type eq "DEL"){
	    print OUT "$type\t$ID\t$present_sv[1]\t$present_sv[2]\t$present_sv[3]\t$present_sv[4]\t-\t-\t-\t$length_SV{$present_sv}\t";
	    my @spe=();
	    foreach my $sv (@{$result{$type}{$ID}}){
		my @tmp = split /_/,$sv;
		push(@spe,$tmp[-1]);
	    }
	    my $spe_out= join "|",@spe;
	    print OUT "$spe_out\n";
	}else{
	    my $target=$hash_SV{$type}{$present_sv[1]}{$present_sv[2]}{$present_sv[3]}{$present_sv[4]};
	    my @target=split /\_/,$target;
	    my @spe=();
	    foreach my $sv (@{$result{$type}{$ID}}){
		my @tmp = split /_/,$sv;
		push(@spe,$tmp[-1]);
	    }
	    my $spe_out=join "|",@spe;
	    print OUT "$type\t$ID\t$present_sv[1]\t$present_sv[2]\t$present_sv[3]\t$present_sv[4]\t$target[1]\t$target[2]\t$target[3]\t$length_SV{$present_sv}\t$spe_out\n";
	}
    }
}
close OUT;
