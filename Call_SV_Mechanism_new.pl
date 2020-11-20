#########################################################################
#      File Name: Call_SV_Mechanism_version2.pl
#    > Author: hwlu
#    > Mail: hongweilu@genetics.ac.cn 
#    Created Time: Thu 16 Jul 2020 02:52:32 PM CST
#########################################################################

#!/usr/bin/perl -w
use strict;

#Call By SV Cluster;

my $sv=shift;
my $cluster_ref=shift;
my $cluster_query=shift;
my $prefix=$sv;
$prefix=~s/\..*//g;


my %hash_Cluster;

open(IN,"$cluster_ref");
while(<IN>){
    chomp;
    my @line = split /\s+/,$_;
    for(my $i=3;$i<@line;$i++){
		$line[$i] =~ s/small:/small_/g;
		$line[$i] =~ s/Infered:/Infered_/g;
		my @tmp = split /:/,$line[$i];
		if($tmp[1]>$tmp[2]){
			($tmp[1],$tmp[2]) = ($tmp[2],$tmp[1]);
			$line[$i] = join ":",@tmp;
		}
		if($line[1] eq "coupling"){
		    $line[1]="Coupling";
		}
		$hash_Cluster{'REF'}{'type'}{$line[$i]}=$line[1];
		$hash_Cluster{'REF'}{'Cluster_ID'}{$line[$i]}=$line[0];
		$hash_Cluster{'REF'}{'Cluster'}{$line[0]}{$line[$i]}=1;
    }
}
close IN;

open(IN,"$cluster_query");
while(<IN>){
    chomp;
    my @line = split /\s+/,$_;
    for(my $i=3;$i<@line;$i++){
		$line[$i] =~ s/small:/small_/g;
		$line[$i] =~ s/Infered:/Infered_/g;
		my @tmp = split /:/,$line[$i];
		if($tmp[1]>$tmp[2]){
			($tmp[1],$tmp[2]) = ($tmp[2],$tmp[1]);
			$line[$i] = join ":",@tmp;
		}
		if($line[1] eq "coupling"){
		    $line[1]="Coupling";
		}
		$hash_Cluster{'Query'}{'type'}{$line[$i]}=$line[1];
		$hash_Cluster{'Query'}{'Cluster_ID'}{$line[$i]}=$line[0];
		$hash_Cluster{'Query'}{'Cluster'}{$line[0]}{$line[$i]}=1;
    }
}
close IN;

my %hash_SV_info;
open(IN,"$sv");
my $total_sv_num=0;
while(<IN>){
    chomp;
	my $line = $_;
    my @line = split /\t/,$line;

	my $sv_length_ref=0;
	my $sv_length_query=0;
	$line[-1] =~ s/:/_/g;
	if($line[-1] =~ /Inversion/ or $line[-1] eq "Translocation"){
		$sv_length_ref=$line[2]-$line[1]+1;
		$sv_length_query=$line[8]-$line[7]+1;
	}elsif($line[-1] eq "DEL" or $line[-1] eq "Infered_INS" or $line[-1] eq "small_DEL") {
		$sv_length_ref=$line[2]-$line[1]+1;
		$sv_length_query=$sv_length_ref;
	}else{
		$sv_length_ref=$line[8]-$line[7]+1;
		$sv_length_query=$sv_length_ref;
	}
	if($line[-1] !~/small/){
		$total_sv_num++;
	}
    $hash_SV_info{'REF'}{"$line[0]:$line[1]:$line[2]:$line[-1]"}{"VNTR"}=$line[3];
    $hash_SV_info{'REF'}{"$line[0]:$line[1]:$line[2]:$line[-1]"}{"TE"}=$line[4];
    $hash_SV_info{'REF'}{"$line[0]:$line[1]:$line[2]:$line[-1]"}{"Homo"}=$line[5];
    $hash_SV_info{'REF'}{"$line[0]:$line[1]:$line[2]:$line[-1]"}{"SV_Type"}=$line[-1];
    $hash_SV_info{'REF'}{"$line[0]:$line[1]:$line[2]:$line[-1]"}{"SV_Length"}=$sv_length_ref;
    $hash_SV_info{'REF'}{"$line[0]:$line[1]:$line[2]:$line[-1]"}{'Target'}="$line[6]:$line[7]:$line[8]:$line[-1]";
#	next if ($line[6] eq "-");
    $hash_SV_info{'Query'}{"$line[6]:$line[7]:$line[8]:$line[-1]"}{'VNTR'}=$line[9];
    $hash_SV_info{'Query'}{"$line[6]:$line[7]:$line[8]:$line[-1]"}{'TE'}=$line[10];
    $hash_SV_info{'Query'}{"$line[6]:$line[7]:$line[8]:$line[-1]"}{'Homo'}=$line[11];
    $hash_SV_info{'Query'}{"$line[6]:$line[7]:$line[8]:$line[-1]"}{'SV_Type'}=$line[-1];
    $hash_SV_info{'Query'}{"$line[6]:$line[7]:$line[8]:$line[-1]"}{'SV_Length'}=$sv_length_query;
    $hash_SV_info{'Query'}{"$line[6]:$line[7]:$line[8]:$line[-1]"}{'Target'}="$line[0]:$line[1]:$line[2]:$line[-1]";

}
close IN;

my %hash_used;
open(OUT,">$prefix.SV_Mechanism.result");
my $Mechanism_ID=1;
my $total_sv_num_2=0;
foreach my $tmp (keys %{$hash_SV_info{'REF'}}){
	next if ($hash_SV_info{'REF'}{$tmp}{'SV_Type'} =~ /small/);
	$total_sv_num_2++;
}
print "$total_sv_num_2\n";

my $out=0;

#Do Inversion:
foreach my $SV (keys %{$hash_SV_info{'REF'}}){
	next if (exists $hash_used{$SV});
	next if ($hash_SV_info{'REF'}{$SV}{'SV_Type'} =~ /small/);
	my $Query_sv=$hash_SV_info{'REF'}{$SV}{'Target'};
	my $SV_Type=$hash_SV_info{'REF'}{$SV}{"SV_Type"};
	my $genome="Query";
	if($SV_Type eq "DEL" or ($SV_Type =~/Infered/ and $SV_Type ne "Infered_DEL")){
	    $genome = "REF";
	}
	my $test_SV=$SV;
	if($genome eq "Query"){
		$test_SV = $Query_sv;
	}

	my $Cluster_ID=$hash_Cluster{$genome}{"Cluster_ID"}{$test_SV};
	my $Cluster_Type=$hash_Cluster{$genome}{'type'}{$test_SV};
	next if ($SV_Type !~  /Inversion/);
	&Call_SV_Mechanism($SV);
}


#Do Translocation:
foreach my $SV (keys %{$hash_SV_info{'REF'}}){
	next if (exists $hash_used{$SV});
	next if ($hash_SV_info{'REF'}{$SV}{'SV_Type'} =~ /small/);
	my $Query_sv=$hash_SV_info{'REF'}{$SV}{'Target'};
	my $SV_Type=$hash_SV_info{'REF'}{$SV}{"SV_Type"};
	my $genome="Query";
	if($SV_Type eq "DEL" or ($SV_Type =~/Infered/ and $SV_Type ne "Infered_DEL")){
	    $genome = "REF";
	}
	my $test_SV=$SV;
	if($genome eq "Query"){
		$test_SV = $Query_sv;
	}

	my $Cluster_ID=$hash_Cluster{$genome}{"Cluster_ID"}{$test_SV};
	my $Cluster_Type=$hash_Cluster{$genome}{'type'}{$test_SV};
	next if ($SV_Type !~  /Translocation/);
	&Call_SV_Mechanism($SV);
}

#Do Cluster:
foreach my $SV (keys %{$hash_SV_info{'REF'}}){
	next if (exists $hash_used{$SV});
	next if ($hash_SV_info{'REF'}{$SV}{'SV_Type'} =~ /small/);
	my $Query_sv=$hash_SV_info{'REF'}{$SV}{'Target'};
	my $SV_Type=$hash_SV_info{'REF'}{$SV}{"SV_Type"};
	my $genome="Query";
	if($SV_Type eq "DEL" or ($SV_Type =~/Infered/ and $SV_Type ne "Infered_DEL")){
	    $genome = "REF";
	}
	my $test_SV=$SV;
	if($genome eq "Query"){
		$test_SV = $Query_sv;
	}

	my $Cluster_ID=$hash_Cluster{$genome}{"Cluster_ID"}{$test_SV};
	my $Cluster_Type=$hash_Cluster{$genome}{'type'}{$test_SV};
	next if ($Cluster_Type ne "Cluster");
	&Call_SV_Mechanism($SV);
}






#Do DEL and coupling:
foreach my $SV (keys %{$hash_SV_info{'REF'}}){
	if($SV eq "Chr10:11747031:11747031:Infered_DEL"){
		#print "Now1..$SV\n";
	}
	next if (exists $hash_used{$SV});
	next if ($hash_SV_info{'REF'}{$SV}{'SV_Type'} =~ /small/);
	my $Query_sv=$hash_SV_info{'REF'}{$SV}{'Target'};
	my $SV_Type=$hash_SV_info{'REF'}{$SV}{"SV_Type"};
	my $genome="Query";
	if($SV_Type eq "DEL" or ($SV_Type =~/Infered/ and $SV_Type ne "Infered_DEL")){
	    $genome = "REF";
	}
	my $test_SV=$SV;
	if($genome eq "Query"){
		$test_SV = $Query_sv;
	}

	my $Cluster_ID=$hash_Cluster{$genome}{"Cluster_ID"}{$test_SV};
	my $Cluster_Type=$hash_Cluster{$genome}{'type'}{$test_SV};

	next if ($SV_Type !~  /DEL/ or $Cluster_Type ne "Coupling");

	&Call_SV_Mechanism($SV);
}

#Do DEL and solo:
foreach my $SV (keys %{$hash_SV_info{'REF'}}){
	next if (exists $hash_used{$SV});
	next if ($hash_SV_info{'REF'}{$SV}{'SV_Type'} =~ /small/);
	my $Query_sv=$hash_SV_info{'REF'}{$SV}{'Target'};
	my $SV_Type=$hash_SV_info{'REF'}{$SV}{"SV_Type"};
	my $genome="Query";
	if($SV_Type eq "DEL" or ($SV_Type =~/Infered/ and $SV_Type ne "Infered_DEL")){
	    $genome = "REF";
	}
	my $test_SV=$SV;
	if($genome eq "Query"){
		$test_SV = $Query_sv;
	}

	my $Cluster_ID=$hash_Cluster{$genome}{"Cluster_ID"}{$test_SV};
	my $Cluster_Type=$hash_Cluster{$genome}{'type'}{$test_SV};
	next if ($SV_Type !~  /DEL/);

	&Call_SV_Mechanism($SV);
}
#Do INS and coupling:
foreach my $SV (keys %{$hash_SV_info{'REF'}}){
	next if (exists $hash_used{$SV});
	next if ($hash_SV_info{'REF'}{$SV}{'SV_Type'} =~ /small/);
	my $Query_sv=$hash_SV_info{'REF'}{$SV}{'Target'};
	my $SV_Type=$hash_SV_info{'REF'}{$SV}{"SV_Type"};
	my $genome="Query";
	if($SV_Type eq "DEL" or ($SV_Type =~/Infered/ and $SV_Type ne "Infered_DEL")){
	    $genome = "REF";
	}
	my $test_SV=$SV;
	if($genome eq "Query"){
		$test_SV = $Query_sv;
	}

	my $Cluster_ID=$hash_Cluster{$genome}{"Cluster_ID"}{$test_SV};
	my $Cluster_Type=$hash_Cluster{$genome}{'type'}{$test_SV};
	next if ($SV_Type !~  /INS/ or $Cluster_Type ne "Coupling");

	&Call_SV_Mechanism($SV);
}



#Do INS:
foreach my $SV (keys %{$hash_SV_info{'REF'}}){
	next if (exists $hash_used{$SV});
	next if ($hash_SV_info{'REF'}{$SV}{'SV_Type'} =~ /small/);
	my $Query_sv=$hash_SV_info{'REF'}{$SV}{'Target'};
	my $SV_Type=$hash_SV_info{'REF'}{$SV}{"SV_Type"};
	my $genome="Query";
	if($SV_Type eq "DEL" or ($SV_Type =~/Infered/ and $SV_Type ne "Infered_DEL")){
	    $genome = "REF";
	}
	my $test_SV=$SV;
	if($genome eq "Query"){
		$test_SV = $Query_sv;
	}

	my $Cluster_ID=$hash_Cluster{$genome}{"Cluster_ID"}{$test_SV};
	my $Cluster_Type=$hash_Cluster{$genome}{'type'}{$test_SV};
	next if ($SV_Type !~  /INS/);
	&Call_SV_Mechanism($SV);
}


my $used = keys %hash_used;
print "$total_sv_num\t$total_sv_num_2\t$used\t$out\n";

foreach my $sv (keys %hash_used){
	if(!exists $hash_SV_info{'REF'}{$sv}){
		print "$sv\n";
	}
}


sub Call_SV_Mechanism{
	my $SV=shift;
	my $Query_sv=$hash_SV_info{'REF'}{$SV}{'Target'};
	my $SV_Type=$hash_SV_info{'REF'}{$SV}{"SV_Type"};
	my $genome="Query";
	if($SV_Type eq "DEL" or ($SV_Type =~/Infered/ and $SV_Type ne "Infered_DEL")){
	    $genome = "REF";
	}
	my $test_SV=$SV;
	if($genome eq "Query"){
		$test_SV = $Query_sv;
	}

	my $Cluster_ID=$hash_Cluster{$genome}{"Cluster_ID"}{$test_SV};
	my $Cluster_Type=$hash_Cluster{$genome}{'type'}{$test_SV};
	$Cluster_Type=~s/coupling/Coupling/;
	
	if($SV_Type =~ /Inversion/){
		&Call_Inversion($SV);
		next;
	}
	
	if($SV_Type =~ /Translocation/){
		&Call_Translocation($SV);
		next;		
	}

	my @svs=keys %{$hash_Cluster{$genome}{'Cluster'}{$Cluster_ID}};
	if ($hash_Cluster{$genome}{'type'}{$test_SV} eq "Cluster"){
		foreach my $sv_tmp (@svs){
			next if ($sv_tmp =~ /small/);
			if($genome eq "Query"){
				$sv_tmp = $hash_SV_info{$genome}{$sv_tmp}{'Target'};
			}
			next if (exists $hash_used{$sv_tmp});
			$hash_used{$sv_tmp}=1;
			$out++;
			print OUT "Mechanism_$Mechanism_ID\t$sv_tmp\t$hash_SV_info{'REF'}{$sv_tmp}{'Target'}\t$hash_Cluster{'REF'}{'Cluster_ID'}{$sv_tmp}\t$hash_SV_info{'REF'}{$sv_tmp}{'SV_Length'}\t$hash_SV_info{'REF'}{$sv_tmp}{'Homo'}\t$hash_SV_info{'Query'}{$hash_SV_info{'REF'}{$sv_tmp}{'Target'}}{'Homo'}\t$Cluster_Type\tNA\t$total_sv_num\n";
		}
		$Mechanism_ID++;
		next;
	}
	
	my $index=0;
	my $flag=0 ; #0: this is highest level ; more than 0 next;

		
	for(my $i=0;$i<@svs;$i++){
		if($genome eq "REF"){
			#$hash_used{$svs[$i]}=1;
		}else{
			my $tmp_target=$hash_SV_info{$genome}{$svs[$i]}{'Target'};
			#$hash_used{$tmp_target}=1;
		}
		if($svs[$i] eq $test_SV){
			$index=$i;
		}
	}
		
	my $attaching_SV="-";    #共断点SV
	if($Cluster_Type eq "Coupling"){
		my $left=$index-0;
		my $right=@svs-$index-1;
		if($left==0){
			$attaching_SV=$svs[$index+1];
		}elsif($right==0){
			$attaching_SV=$svs[$index-1];
		}else{
			$attaching_SV=$svs[$index-1];
			if($hash_SV_info{$genome}{$attaching_SV}{'SV_Length'}<$hash_SV_info{$genome}{$svs[$index+1]}{"SV_Length"}){
				$attaching_SV=$svs[$index-1];
			}
		}
	}

	my $Mechanism;
	my $homo=$hash_SV_info{$genome}{$test_SV}{'Homo'};
	my $TE = $hash_SV_info{$genome}{$test_SV}{'TE'};
	$TE=~s/TE://g;
	my $VNTR = $hash_SV_info{$genome}{$test_SV}{'VNTR'};
	$VNTR=~s/VNTR://g;
	my $attaching_sv_size=0;
	if($Cluster_Type eq "solo"){
		$attaching_sv_size=0;
	}else{
		$attaching_sv_size = $hash_SV_info{$genome}{$attaching_SV}{"SV_Length"};
	}
	if($SV_Type =~  /DEL/ ){
		if(($SV_Type =~ /Infered/ and $attaching_SV =~ /:DEL/) or ($SV_Type !~/Infered/ and $attaching_SV =~/Infered_DEL/)){
			$attaching_SV=~s/DEL/INS/g;
		}
		if($SV_Type !~/Infered/ and $attaching_SV =~/Infered_INS/){
			$attaching_SV=~s/INS/DEL/g;
		}
		$Mechanism=&Call_Deletion($genome,$TE,$VNTR,$attaching_SV,$Cluster_Type,$homo);
	}else{
		if($attaching_SV =~/INS/){
			$attaching_sv_size=0;
		}
		$Mechanism=&Call_Insertion($TE,$VNTR,$attaching_sv_size);
	}

	for(my $i=0;$i<@svs;$i++){
		my $sv_type_tmp = (split /:/,$svs[$i])[-1];
		next if ($sv_type_tmp =~ /small/);
		if($genome eq "REF"){
			next if ($hash_used{$svs[$i]}==1);
			$hash_used{$svs[$i]}=1;
			$out++;
			print OUT "Mechanism_$Mechanism_ID\t$svs[$i]\t$hash_SV_info{'REF'}{$svs[$i]}{'Target'}\t$hash_Cluster{'REF'}{'Cluster_ID'}{$svs[$i]}\t$hash_SV_info{'REF'}{$svs[$i]}{'SV_Length'}\t$hash_SV_info{'REF'}{$svs[$i]}{'Homo'}\t$hash_SV_info{'Query'}{$hash_SV_info{'REF'}{$svs[$i]}{'Target'}}{'Homo'}\t$Cluster_Type\t$Mechanism\t$total_sv_num\n";
		}else{
			if(!exists $hash_SV_info{$genome}{$svs[$i]}{'Target'}){
				print "no target: $svs[$i]\n";
			}
			my $target_tmp = $hash_SV_info{$genome}{$svs[$i]}{'Target'};
			next if ($hash_used{$target_tmp}==1);
			$hash_used{$target_tmp}=1;
			$out++;
			print OUT "Mechanism_$Mechanism_ID\t$target_tmp\t$hash_SV_info{'REF'}{$target_tmp}{'Target'}\t$hash_Cluster{'REF'}{'Cluster_ID'}{$target_tmp}\t$hash_SV_info{'REF'}{$target_tmp}{'SV_Length'}\t$hash_SV_info{'REF'}{$target_tmp}{'Homo'}\t$hash_SV_info{'Query'}{$hash_SV_info{'REF'}{$target_tmp}{'Target'}}{'Homo'}\t$Cluster_Type\t$Mechanism\t$total_sv_num\n";

		}	
	}
	$Mechanism_ID++;
}
close OUT;


sub Call_Inversion{
	my $sv = shift;
	my $query_sv = $hash_SV_info{'REF'}{$sv}{'Target'};
	my $SV_Type = $hash_SV_info{'REF'}{$sv}{'SV_Type'};
	my $genome = "Query";
	if($SV_Type =~ /Infered/){
		$genome = "REF";
	}
	my $Cluster_ID_ref = $hash_Cluster{'REF'}{'Cluster_ID'}{$sv};
	my $Cluster_ID_query = $hash_Cluster{'Query'}{'Cluster_ID'}{$query_sv};
	
	my @svs_ref = keys %{$hash_Cluster{'REF'}{'Cluster'}{$Cluster_ID_ref}};
	my @svs_query = keys %{$hash_Cluster{'Query'}{'Cluster'}{$Cluster_ID_query}};
	
	my $Cluster_Type_ref = $hash_Cluster{'REF'}{'type'}{$sv};
	my $Cluster_Type_query = $hash_Cluster{'Query'}{'type'}{$query_sv};
	my $mechanism;
	my $Cluster_Type;
	if(@svs_ref>2 or @svs_query>2){
		$Cluster_Type="Cluster";
		my @test_svs=@svs_query;
		if($SV_Type =~ /Infered/){
			@test_svs  =  @svs_ref;
		}
		my $max=0;
		foreach my $tmp_svs (@test_svs){
			if($max < $hash_SV_info{$genome}{$tmp_svs}{'SV_Length'}){
				$max = $hash_SV_info{$genome}{$tmp_svs}{'SV_Length'};
			}
		}
		$mechanism = "NA";
		if($max <10){
			$mechanism = "NHEJ";
		}
	}else{
		if(@svs_ref  ==1){
			if(@svs_query ==1 ){   #REF 
				$Cluster_Type = "solo";
				if($hash_SV_info{'REF'}{$sv}{'Homo'}>200){
					$mechanism = "NAHR";
				}else{
					#$mechanism = "alt_EJ/NHEJ";
					if($hash_SV_info{'REF'}{$sv}{'Homo'}>2){
						$mechanism = "alt_EJ";
					}else{
						$mechanism = "NHEJ";
					}
				}
			}else{      #Query  coupling 
				$Cluster_Type = "Coupling";
				my $sv_coupling = $svs_query[0];
				if($sv_coupling eq $query_sv){
					$sv_coupling = $svs_query[1];
				}
				if($SV_Type !~ /Infered/){
					if($hash_SV_info{'Query'}{$sv_coupling}{'SV_Type'} =~ /small_INS/ or $hash_SV_info{'Query'}{$sv_coupling}{'SV_Type'} eq "INS" or $hash_SV_info{'Query'}{$sv_coupling}{'SV_Type'} eq "Translocation" or $hash_SV_info{'Query'}{$sv_coupling}{'SV_Type'} eq "Infered_DEL"){
						if($hash_SV_info{'Query'}{$sv_coupling}{'SV_Length'}>10){
							$mechanism = "FoSTeS/MMBIR";
						}else{
							$mechanism = "NHEJ";
						}
					}else{
						$mechanism = "NA";
					}
				}else{
					if($hash_SV_info{'Query'}{$sv_coupling}{'SV_Type'} eq "DEL" or $hash_SV_info{'Query'}{$sv_coupling}{'SV_Type'} eq "small_DEL" or $hash_SV_info{'Query'}{$sv_coupling}{'SV_Type'} eq "Infered_INS"){
						if($hash_SV_info{'Query'}{$sv_coupling}{'SV_Length'}>10){
							$mechanism = "FoSTeS/MMBIR";
						}else{
							$mechanism = "NHEJ";
						}
					}else{
						$mechanism = "NA";
					}
				}
			}
		}else{   #@svs_ref >1 
			if(@svs_query	>2 or @svs_ref >2){   #REF 
				$Cluster_Type ="Cluster";
				$mechanism = "NA";
			}elsif(@svs_query==2){   #REF coupling 
				my $attach_ref=$svs_ref[0];
				if($attach_ref eq $sv){
					$attach_ref=$svs_ref[1];
				}
				my $attach_query=$svs_query[0];
				if($attach_query eq $query_sv){
					$attach_query=$svs_query[1];
				}
				if($hash_SV_info{'REF'}{$attach_ref}{'Target'} ne $attach_query){
					$Cluster_Type ="Cluster";
					$mechanism = "NA";
				}else{
					$Cluster_Type ="Coupling";
					my $sv_coupling = $svs_ref[0];
					if($sv_coupling eq $sv){
						$sv_coupling = $svs_ref[1];
					}
					if($SV_Type !~ /Infered/){ 
						if($hash_SV_info{'REF'}{$sv_coupling}{'SV_Type'} =~ /small_INS/ or $hash_SV_info{'REF'}{$sv_coupling}{'SV_Type'} eq "INS" or $hash_SV_info{'REF'}{$sv_coupling}{'SV_Type'} eq "Translocation"){
							if($hash_SV_info{'REF'}{$sv_coupling}{'SV_Length'}>10){
								$mechanism = "FoSTeS/MMBIR";
							}else{
								$mechanism = "NHEJ";
							}
						}else{
							$mechanism = "NA";
						}	
					}else{
						if($hash_SV_info{'REF'}{$sv_coupling}{'SV_Type'} eq "small_INS" or $hash_SV_info{'REF'}{$sv_coupling}{'SV_Type'} eq "INS" or $hash_SV_info{'REF'}{$sv_coupling}{'SV_Type'} eq "Infered_DEL"){
							if($hash_SV_info{'REF'}{$sv_coupling}{'SV_Length'}>10){
							$mechanism = "FoSTeS/MMBIR";
							}else{
								$mechanism = "NHEJ";
							}
						}else{
							$mechanism = "NA";
						}
					}
				}
			}else{      #REF  coupling 
				$Cluster_Type ="Coupling";
				my $sv_coupling = $svs_ref[0];
				if($sv_coupling eq $sv){
					$sv_coupling = $svs_ref[1];
				}
				if($SV_Type !~ /Infered/){ 
					if($hash_SV_info{'REF'}{$sv_coupling}{'SV_Type'} =~ /small_INS/ or $hash_SV_info{'REF'}{$sv_coupling}{'SV_Type'} eq "INS" or $hash_SV_info{'REF'}{$sv_coupling}{'SV_Type'} eq "Translocation"){
						if($hash_SV_info{'REF'}{$sv_coupling}{'SV_Length'}>10){
							$mechanism = "FoSTeS/MMBIR";
						}else{
							$mechanism = "NHEJ";
						}
					}else{
						$mechanism = "NA";
					}	
				}else{
					if($hash_SV_info{'REF'}{$sv_coupling}{'SV_Type'} eq "small_INS" or $hash_SV_info{'REF'}{$sv_coupling}{'SV_Type'} eq "INS" or $hash_SV_info{'REF'}{$sv_coupling}{'SV_Type'} eq "Infered_DEL"){
						if($hash_SV_info{'REF'}{$sv_coupling}{'SV_Length'}>10){
							$mechanism = "FoSTeS/MMBIR";
						}else{
							$mechanism = "NHEJ";
						}
					}else{
						$mechanism = "NA";
					}
				}
			}	
		}
	}
	
	##NAHR：
	
	if($hash_SV_info{'REF'}{$sv}{'Homo'} >= 200){
		$mechanism = "NAHR";
	}
	foreach my $sv_tmp (@svs_ref){
		next if (exists $hash_used{$sv_tmp});
		next if ($sv_tmp =~ /small/);
		$hash_used{$sv_tmp}=1;
		$out++;
		if($mechanism eq "NAHR" and $sv_tmp ne $sv){
			print OUT "Mechanism_$Mechanism_ID\t$sv_tmp\t$hash_SV_info{'REF'}{$sv_tmp}{'Target'}\t$hash_Cluster{'REF'}{'Cluster_ID'}{$sv_tmp}\t$hash_SV_info{'REF'}{$sv_tmp}{'SV_Length'}\t$hash_SV_info{'REF'}{$sv_tmp}{'Homo'}\t$hash_SV_info{'Query'}{$hash_SV_info{'REF'}{$sv_tmp}{'Target'}}{'Homo'}\t$Cluster_Type\tNA\t$total_sv_num\n";
		}else{
		print OUT "Mechanism_$Mechanism_ID\t$sv_tmp\t$hash_SV_info{'REF'}{$sv_tmp}{'Target'}\t$hash_Cluster{'REF'}{'Cluster_ID'}{$sv_tmp}\t$hash_SV_info{'REF'}{$sv_tmp}{'SV_Length'}\t$hash_SV_info{'REF'}{$sv_tmp}{'Homo'}\t$hash_SV_info{'Query'}{$hash_SV_info{'REF'}{$sv_tmp}{'Target'}}{'Homo'}\t$Cluster_Type\t$mechanism\t$total_sv_num\n";
	}
}
	foreach my $sv_tmp (@svs_query){
		my $target_tmp = $hash_SV_info{"Query"}{$sv_tmp}{'Target'};
		next if ($target_tmp =~ /small/);
		next if (exists $hash_used{$target_tmp});
		$hash_used{$target_tmp}=1;
		$out++;
		if($mechanism eq "NAHR" and $target_tmp ne $sv){
			print OUT "Mechanism_$Mechanism_ID\t$target_tmp\t$hash_SV_info{'REF'}{$target_tmp}{'Target'}\t$hash_Cluster{'REF'}{'Cluster_ID'}{$target_tmp}\t$hash_SV_info{'REF'}{$target_tmp}{'SV_Length'}\t$hash_SV_info{'REF'}{$target_tmp}{'Homo'}\t$hash_SV_info{'Query'}{$hash_SV_info{'REF'}{$target_tmp}{'Target'}}{'Homo'}\t$Cluster_Type\tNA\t$total_sv_num\n";
		}else{
		print OUT "Mechanism_$Mechanism_ID\t$target_tmp\t$hash_SV_info{'REF'}{$target_tmp}{'Target'}\t$hash_Cluster{'REF'}{'Cluster_ID'}{$target_tmp}\t$hash_SV_info{'REF'}{$target_tmp}{'SV_Length'}\t$hash_SV_info{'REF'}{$target_tmp}{'Homo'}\t$hash_SV_info{'Query'}{$hash_SV_info{'REF'}{$target_tmp}{'Target'}}{'Homo'}\t$Cluster_Type\t$mechanism\t$total_sv_num\n";
		}
	}
	$Mechanism_ID++;
}



sub Call_Translocation{
	my $sv = shift;
	my $sv_query = $hash_SV_info{'REF'}{$sv}{'Target'};
	
	my $Cluster_ID_ref = $hash_Cluster{'REF'}{"Cluster_ID"}{$sv};
	my $Cluster_ID_query = $hash_Cluster{'Query'}{"Cluster_ID"}{$sv_query};
	
	my @svs_ref = keys %{$hash_Cluster{'REF'}{"Cluster"}{$Cluster_ID_ref}};
	my @svs_query = keys %{$hash_Cluster{'Query'}{"Cluster"}{$Cluster_ID_query}};
	
	my $TE;
	my $Homo;
	my $genome="Query";
	if($sv=~/Infered/){
		$TE=$hash_SV_info{'REF'}{$sv}{'TE'};
		$Homo = $hash_SV_info{'REF'}{$sv}{'Homo'};
		$genome = "REF";
	}else{
		$TE=$hash_SV_info{'Query'}{$sv_query}{'TE'};
		$TE=$hash_SV_info{'Query'}{$sv_query}{'Homo'};
	}
	my $Cluster_Type="Coupling";
	if(@svs_ref>2 or  @svs_query > 2){
		$Cluster_Type = "Cluster";
	}elsif(@svs_ref==1 and @svs_query ==1){
		$Cluster_Type = "solo";
	}elsif(@svs_ref==2 and @svs_query==2){
		my $attach_ref=$svs_ref[0];
		if($attach_ref eq $sv){
			$attach_ref=$svs_ref[1];
		}
		my $attach_query=$svs_query[0];
		if($attach_query eq $sv_query){
			$attach_query=$svs_query[1];
		}
		if($hash_SV_info{'REF'}{$attach_ref}{'Target'} ne $attach_query){
			$Cluster_Type ="Cluster";
		}else{
			$Cluster_Type="Coupling";
		}
	}else{
		$Cluster_Type="Coupling";
	}
	my $mechanism;
	
	if($TE>0.9){
	    $mechanism =  "TEI";
	}else{
		if(@svs_ref==1 and @svs_query ==1){
			if($Homo >200){
				$mechanism = "NAHR";
			}else{
				#$mechanism = "alt_EJ/NHEJ";
				if($Homo > 2){
					$mechanism = "alt_EJ";
				}else{
					$mechanism = "NHEJ";
				}

			}
		}else{
			if($Cluster_Type eq  "Cluster"){
				my @test_svs = @svs_ref;
				if($sv !~ /Infered/){
					@test_svs = @svs_query;
				}
				my $max=0;
				foreach my $test_svs (@test_svs){
					if($max < $hash_SV_info{$genome}{$test_svs}{"SV_Length"}){
						$max = $hash_SV_info{$genome}{$test_svs}{"SV_Length"};
					}
				}
				$mechanism="NA";
				if($max < 10){
					$mechanism = "NHEJ";
				}
			}else{
			if(@svs_ref ==1 ) {
				my $attaching_sv = $svs_query[0];
				if($attaching_sv eq $sv_query){
					$attaching_sv = $svs_query[1];
				}
				if($sv=~/Infered/){
					if($hash_SV_info{'Query'}{$attaching_sv}{"SV_Type"} eq "DEL" or $hash_SV_info{'Query'}{$attaching_sv}{"SV_Type"} eq "small_DEL" or $hash_SV_info{'Query'}{$attaching_sv}{"SV_Type"} eq "Infered_INS"){
						if($hash_SV_info{'Query'}{$attaching_sv}{"SV_Length"}>10){
							$mechanism = "FoSTeS/MMBIR";
						}else{
							$mechanism = "NHEJ";
						}
					}else{
						$mechanism = "NA";
					}
				}else{
					if($hash_SV_info{'Query'}{$attaching_sv}{"SV_Type"} eq "INS" or $hash_SV_info{'Query'}{$attaching_sv}{"SV_Type"} eq "small_INS" or $hash_SV_info{'Query'}{$attaching_sv}{"SV_Type"} eq "Infered_DEL"){
						if($hash_SV_info{'Query'}{$attaching_sv}{"SV_Length"}>10){
							$mechanism = "FoSTeS/MMBIR";
						}else{
							$mechanism = "NHEJ";
						}
					}else{
						$mechanism = "NA";
					}
				}
			}elsif (@svs_query ==1 ){
				my $attaching_sv = $svs_ref[0];
				if($attaching_sv eq $sv_query){
					$attaching_sv = $svs_ref[1];
				}
				if($sv=~/Infered/){
					if($hash_SV_info{'REF'}{$attaching_sv}{"SV_Type"} eq "INS" or $hash_SV_info{'REF'}{$attaching_sv}{"SV_Type"} eq "small_INS" or $hash_SV_info{'REF'}{$attaching_sv}{"SV_Type"} eq "Infered_DEL"){
						if($hash_SV_info{'REF'}{$attaching_sv}{"SV_Length"}>10){
							$mechanism = "FoSTeS/MMBIR";
						}else{
							$mechanism = "NHEJ";
						}
					}else{
						$mechanism = "NA";
					}
				}else{
					if($hash_SV_info{'REF'}{$attaching_sv}{"SV_Type"} eq "DEL" or $hash_SV_info{'REF'}{$attaching_sv}{"SV_Type"} eq "small_DEL" or $hash_SV_info{'REF'}{$attaching_sv}{"SV_Type"} eq "Infered_INS"){
						if($hash_SV_info{'REF'}{$attaching_sv}{"SV_Length"}>10){
							$mechanism = "FoSTeS/MMBIR";
						}else{
							$mechanism = "NHEJ";
						}
					}else{
						$mechanism = "NA";
					}
				}
			}else{
				$mechanism = "NA";
				#$Cluster_Type = "Cluster";
			}
			}
		
		}
	}
	foreach my $sv_tmp (@svs_ref){
		next if (exists $hash_used{$sv_tmp});
		next if ($sv_tmp =~ /small/);
		$hash_used{$sv_tmp}=1;
		$out++;
		print OUT "Mechanism_$Mechanism_ID\t$sv_tmp\t$hash_SV_info{'REF'}{$sv_tmp}{'Target'}\t$hash_Cluster{'REF'}{'Cluster_ID'}{$sv_tmp}\t$hash_SV_info{'REF'}{$sv_tmp}{'SV_Length'}\t$hash_SV_info{'REF'}{$sv_tmp}{'Homo'}\t$hash_SV_info{'Query'}{$hash_SV_info{'REF'}{$sv_tmp}{'Target'}}{'Homo'}\t$Cluster_Type\t$mechanism\t$total_sv_num\n";

	}
	foreach my $sv_tmp (@svs_query){

		my $target_tmp = $hash_SV_info{"Query"}{$sv_tmp}{'Target'};
		next if ($target_tmp =~ /small/);
		next if (exists $hash_used{$target_tmp});
		$hash_used{$target_tmp}=1;
		$out++;
		print OUT "Mechanism_$Mechanism_ID\t$target_tmp\t$hash_SV_info{'REF'}{$target_tmp}{'Target'}\t$hash_Cluster{'REF'}{'Cluster_ID'}{$target_tmp}\t$hash_SV_info{'REF'}{$target_tmp}{'SV_Length'}\t$hash_SV_info{'REF'}{$target_tmp}{'Homo'}\t$hash_SV_info{'Query'}{$hash_SV_info{'REF'}{$target_tmp}{'Target'}}{'Homo'}\t$Cluster_Type\t$mechanism\t$total_sv_num\n";
	}
	$Mechanism_ID++;
}

sub Call_Deletion{
    my $genome=shift;  #REF or Query 
    my $TE=shift;
    my $VNTR=shift;
    my $attaching_SV=shift;  #sv_ID
    my $type=shift;
    my $homo=shift;
    if($TE>0.9){
		return "TEI";
    }
    if($VNTR>0.9){
		return "VNTR";
    }
    
    if($type eq "solo"){
		if($homo > 200){
			return "NAHR";
		}else{
			if($homo >1){
				if($homo>50){
					return "alt_EJ";
				}else{
					return "SSA";
				}
			}else{
				return "NHEJ";
			}
		}
    }else{
		my $attaching_SV_type = (split /:/,$attaching_SV)[-1];

		if($attaching_SV_type !~ /INS/){
			return "NA";
		}else{
			my $attaching_SV_TE=$hash_SV_info{$genome}{$attaching_SV}{'TE'};
			my $attaching_SV_Length=$hash_SV_info{$genome}{$attaching_SV}{"SV_Length"};
			if($attaching_SV_TE>0.9){
				return "TE_mediated DEL";
			}else{
				if($attaching_SV_Length >10){
					return "FoSTeS/MMBIR";
				}else{
					return "NHEJ";
				}
			}
		}	
    }
}

sub Call_Insertion{
	my $TE=shift;
	my $VNTR=shift;
	my $attaching_sv_size=shift;
	#if($attaching_sv_size>0){
	#	if($attaching_sv_size > 10){
	#		return "FoSTes/MMBIR";
	#	}else{
	#		return "NHEJ";
	#	}
	#}
	if($TE>0.9){
		return "TEI";
	}else{
		if($VNTR > 0.9){
			return "VNTR";
		}else{
			return "NA";
		}
	}
}
