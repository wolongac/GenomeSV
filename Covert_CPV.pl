#########################################################################
#      File Name: Covert_CPV.pl
#    > Author: hwlu
#    > Mail: hongweilu@genetics.ac.cn 
#    Created Time: Thu 23 Apr 2020 10:10:24 AM CST
#########################################################################

#!/usr/bin/perl -w
use strict;
#CPG,CPL,TDM infer the real position variated using the alignment information.
#Merge joining Deletion
#the input is the output of script : SyRI_Parse.pl

my $bed = shift;

my %hash_CPV;
my %hash;
my %hash_ins;
open(IN,"$bed");
open(OUT,">$bed.mergeDel.CPV2InDel");
while(<IN>){
    chomp;
    my $line=$_;
    my @line = split /\s+/,$line;
    if($line[6] eq "DEL" or $line[6] eq "Deletion"){
		$line[6]="Deletion";
		if($line[7]=~/copyloss/){
			#$line[6]="copyloss";
		}
		$line=join "\t",@line;
		$line[2]+=1;
		$line[1]-=1;
		$hash{'DEL'}{$line[0]}{$line[1]}{$line[2]}=$line;
    }elsif($line[7] eq "CPG" or $line[7] eq "CPL" or $line[7] eq "TDM"){
		$hash_CPV{$line[0]}{$line[1]}{$line[2]}{'target'}="$line[3]\t$line[4]\t$line[5]";
		$hash_CPV{$line[0]}{$line[1]}{$line[2]}{'parent'}="$line[9]";
		$hash_CPV{$line[0]}{$line[1]}{$line[2]}{'type'}="$line[7]\t$line[8]";
    }else{
		next if ($line[7] eq "HDR" or $line[7] eq "SYN");
		if($line[7] eq "DUP_copygain"){
		    $line[6] = "copygain";
		}elsif($line[7] eq "INV"){
		    $line[6] = "Inversion";
		}elsif($line[6] eq "INVTR" or $line[6] eq "TRANS"){
		    $line[6] ="Translocation";
		}elsif($line[7] eq "SNP"){
		    next;
		}elsif($line[6] eq "INS"){
		    $line[6]="Insertion";
		}
		$line = join "\t",@line;
		if($line[6] eq "Insertion"){
		    $hash_ins{$line[3]}{$line[4]}{$line[5]}=$line;
		}else{
			print OUT  "$line\n";
	}
    }
}

close IN;

#Merge Deletion:
foreach my $chr (keys %{$hash{"DEL"}}){
	my $start_now=0;
	my $end_now=0;
	 foreach my $start (sort {$a <=> $b}  keys %{$hash{"DEL"}{$chr}}){
	    foreach my $end (sort {$a <=> $b} keys %{$hash{"DEL"}{$chr}{$start}}){
		    if($start > $end_now){
			#End and output:
			if($end_now==0){
			    $start_now = $start;
			    $end_now = $end;
			    next;
			}else{
			    if(exists $hash{"DEL"}{$chr}{$start_now}{$end_now}){
				print OUT  "$hash{'DEL'}{$chr}{$start_now}{$end_now}\n";
			    }else{
				my $start_tmp=$start_now+1;
				my $end_tmp = $end_now-1;
				print OUT  "$chr\t$start_tmp\t$end_tmp\t-\t-\t-\tDeletion\tCombine\tCombine\t-\n";
				}
			    $start_now=$start;
			    $end_now = $end;
			    next;
			}
		    }else{
			next if ($end < $end_now);
			$end_now = $end;
		    }
	    }
	 }
}

my $syri_out=$bed;
open(ERR,">$syri_out.CPVerr");
$syri_out  =~ s/\.SV\.bed//g;
print "$syri_out\n";
my %hash_SYNAL;
my %hash_ins_raw;
open(IN,"$syri_out");
while(<IN>){
    chomp;
    my @line = split /\s+/,$_;
    if ($line[8] =~ /SYNAL/){
    $hash_SYNAL{$line[9]}{$line[0]}{$line[1]}{$line[2]}="$line[5]\t$line[6]\t$line[7]";
    }elsif($line[-1] eq "copygain"){
	$hash_ins_raw{$line[5]}{$line[6]}{$line[7]}=1;
    }
}
foreach my $chr (keys %hash_ins){
    foreach my $start (sort {$a <=> $b} keys %{$hash_ins{$chr}}){
	foreach my $end (keys %{$hash_ins{$chr}{$start}}){
	    if(exists $hash_ins_raw{$chr}{$start}{$end}){
		my @tmp =split /\t/,$hash_ins{$chr}{$start}{$end};
		#$tmp[6]="copygain";
		$hash_ins{$chr}{$start}{$end} = join "\t",@tmp;
	    }
	    print OUT "$hash_ins{$chr}{$start}{$end}\n";
	}
    }
}
#extract CPV:
foreach my $chr (keys %hash_CPV) {
    foreach my $start (keys %{$hash_CPV{$chr}}) {
	foreach my $end (keys %{$hash_CPV{$chr}{$start}}){
	    my $target = $hash_CPV{$chr}{$start}{$end}{'target'};
	    my @target=split /\t/,$target;
	    my $parent=$hash_CPV{$chr}{$start}{$end}{'parent'};
	    #print "test\t$parent\t$chr\t$start\t$end\n";
	    my @synal_1=();
	    my @synal_2=();
	    my $stat=0;
	    foreach my $start_syn (sort {$a <=> $b} keys %{$hash_SYNAL{$parent}{$chr}}){
		    last if ($start_syn > $end );
		foreach my $end_syn (keys %{$hash_SYNAL{$parent}{$chr}{$start_syn}}){
		    next if ($end_syn < $start);
		    if($start_syn<$start){
			if (@synal_1>1){
				#print "exists synal_1 $chr $start $end @synal_1 $start_syn $end_syn ";
			    $stat=1;
			}
			@synal_1=($start_syn,$end_syn);
		    }else{
			if(@synal_2>1){
				#print "exists synal_2 $chr $start $end @synal_2 $start_syn $end_syn ";
			    $stat=1;
			}
			@synal_2=($start_syn,$end_syn);
		    }
		}
	    }
	    if (@synal_1 ==0 or @synal_2==0 or $stat==1){
		    print ERR  "Can not parse CPV : $chr\t$start\t$end\n";
		next;
	    }
	    if($target[2]-$target[1]>$end-$start){    #more 
		my $target_syn_1=$hash_SYNAL{$parent}{$chr}{$synal_1[0]}{$synal_1[1]};
		my $target_syn_2=$hash_SYNAL{$parent}{$chr}{$synal_2[0]}{$synal_2[1]};
		my @target_1=split /\t/,$target_syn_1;
		my @target_2=split /\t/,$target_syn_2;
		#if($target_2[1]-$target_1[2]>0){
		  #  $hash_pav{'query'}{$chr}{$target_1[2]}{$target_2[1]}="gain";
		  	print OUT "$chr\t$synal_2[0]\t$synal_2[0]\t$chr\t$target_1[2]\t$target[2]\tInsertion\t$hash_CPV{$chr}{$start}{$end}{'type'}\t$parent\n";
			
		  #}
		  #print OUT "$chr\t$synal_1[1]\t$synal_1[1]\t$chr\t$target_2[1]\t$target[2]\tCPG\t$hash_CPV{$chr}{$start}{$end}{'type'}\t$parent\n";
		#$hash_cnv{'query'}{$chr}{$target_2[1]}{$target[2]}="gain";
		#print "query\t$chr\t$target_2[1]\t$target[2]\n";    
	}else{    #less
		#if($synal_1[1]-$synal_2[0]>0){
		    #$hash_pav{'ref'}{$chr}{$synal_2[0]}{$synal_1[2]}="loss";
			print OUT "$chr\t$synal_1[1]\t$end\t-\t-\t-\tDeletion\t$hash_CPV{$chr}{$start}{$end}{'type'}\t$parent\n";
			#}
		#print OUT "$chr\t$synal_2[0]\t$end\t-\t-\t-\tCPL\t$hash_CPV{$chr}{$start}{$end}{'type'}\t$parent\n";
		#$hash_cnv{'ref'}{$chr}{$synal_1[1]}{$end}="loss";
		#print "ref\t$chr\t$synal_1[1]\t$end\n";    
	    }
	}
    }
}
close OUT;
close IN;
close ERR;
