#Syri out Parse 

#In Details :
#1. Insertion(Trans),Insertion(copygain) and Insertion(Notal) infer the corresponding insertion position  on REF genome using the alignment information.
#3. For the contiguous insertion-type SVs,for example: Insertion, Insertion(copygain), Insertion(Notal), and Insertion(Trans), merged into a larger SVs. For TRANS/INVTR, the corresponding Deletions part on ref are converted to Deletions.

#####################Input file format###########################
#1	chromosome ID in genome ref	string
#2	genome ref start position (1-based, includes start position)	int
#3	genome ref end position (1-based, includes end position)	int
#4	Sequence in genome ref (Only for SNPs and indels)	string
#5	Sequence in genome query (Only for SNPs and indels)	string
#6	chromosome ID in genome query	string
#7	genome query start position (1-based, includes start position)	int
#8	genome query end position (1-based, includes end position)	int
#9	Unique ID (annotation type + number)	string
#10	Parent ID (annotation type + number)	string  (for structual variation, parent ID  is "-")
#11	Annotation type	string
#12	Copy status (for duplications)	string
#
###################################################################

####################output file format###########################
#1	chromosome ID in genome ref	string
#2	genome ref start position (1-based, includes start position)	int
#3	genome ref end position (1-based, includes end position)	int
#4	chromosome ID in genome query	string
#5	genome query start position (1-based, includes start position)	int
#6	genome query end position (1-based, includes end position)	int
#7	New SV type
#8	Raw SV type from SyRI
#9	Uniq ID
#10	Parent ID
###################################################################

my $syri_out=shift;

open(IN,"$syri_out");
open(OUT,">$syri_out.SV.bed");

my %hash_ref;				#hash_ref{$chr}{$start}{$end}="SV_Type_info"
##########################################################################
#The SV_Type_info is :
#For DUP/INVDP or NOTAL : "DUP_Deletion","NOTAL_Deletion"
#For INV and SYN : SV_Type_info is "Chr_query\tStart_query\tEnd_query"
#For TRANS/INVTR : "Struc_Deletion|SV_ID"
###########################################################################

my %hash_query;				#$hash_query{$chr}{$start}{$end}="SV_Type_info"
###########################################################################
#The SV_Type_info is : 
#For DUP/INVDP and NOTAL : "DUP_Insertion|SV_ID","NOTAL_Insertion|SV_ID"
#For INV and SYN : SV_Type_info is "Chr_ref\tStart_ref\tEnd_ref"
#For TRANS/INVTR : "Struc_Insertion|SV_ID"
############################################################################

my %hash_syn_var;			
#$hash_syn_var{$sv_name}{'ref'}=$position_on_query
#$hash_syn_var{$sv_name}{'query'}=$position_on_ref
#$hash_syn_var{$sv_name}{'info'}="$line[0]\t$line[1]\t$line[2]\t$line[5]\t$line[6]\t$line[7]\t$line[10]\t$line[10]\t$line[8]\t$line[9]";
my %hash_replacement;  #replacement region : $hash_replacement{'ref'}{$chr}{$start}{$end}  or  $hash_replacement{'query'}{$chr}{$start}{$end}
my %hash_trans;
while(<IN>){
	chomp;
	my $line = $_;
	my @line = split /\s+/,$line;
	#####################################
	#SYN INV DUP INVDP TRANS INVTR NOTAL  are Structual Difference(StrucD)
	#SNP INS DEL CPG CPL HDR TDM are Sequencial Difference(SequentialD)
	#INVAL INVDPAL INVTRAL SYNAL TRANSAL DUPAL are Alignment(AL)
	#####################################
	
	#skip SV on the non_synteny regions like : DUP,TRANs or INV et al
	#Only keep syntenic information for the following analysis
	next if ($line[9] !~ /SYN/ and $line[9] ne "-");
	
	if (($line[10] eq "DUP" or $line[10] eq "INVDP") and $line[11] eq "copyloss"){
		#DUP is StrucD,Convert DUP copyloss to Deletion 
		print OUT "$line[0]\t$line[1]\t$line[2]\t-\t-\t-\tDeletion\tDUP_copyloss\t$line[8]\t$line[9]\n";
		$hash_ref{$line[0]}{$line[1]}{$line[2]}="DUP_Deletion";
	}elsif(($line[10] eq "DUP" or $line[10] eq "INVDP") and $line[11] eq "copygain"){
		#DUP is StrucD,Convert DUP copyloss to Insertion
		$hash_query{$line[5]}{$line[6]}{$line[7]}="DUP_Insertion|$line[8]";
	}elsif($line[10] eq "NOTAL" and $line[0] ne "-"){
		#NOTAL is StrucD,Convert NOTAL of Ref to  Deletion
		print OUT "$line[0]\t$line[1]\t$line[2]\t-\t-\t-\tDeletion\tNOTAL_Deletion\t$line[8]\t$line[9]\n";
		$hash_ref{$line[0]}{$line[1]}{$line[2]}="NOTAL_Deletion";
	}elsif($line[10] eq "NOTAL" and $line[1] eq "-"){
		#NOTAL is StrucD,Convert NOTAL of Query to Insertion
		$hash_query{$line[5]}{$line[6]}{$line[7]}="NOTAL_Insertion|$line[8]";
	}elsif($line[10] eq "CPG" or $line[10] eq "CPL" or $line[10] eq "TDM" or $line[10] eq "HDR"){
		#CGP,CGL,TDM,HDR of SYN are SequentialD, Convert then to Region_replacement 
		print OUT "$line[0]\t$line[1]\t$line[2]\t$line[5]\t$line[6]\t$line[7]\tRegion_replacement\t$line[10]\t$line[8]\t$line[9]\n";
		$hash_replacement{'ref'}{$line[0]}{$line[1]}=$line[2];
		$hash_replacement{'query'}{$line[5]}{$line[6]}=$line[7];
	}elsif($line[10] eq "INV"){
		#INV is StrucD,keep the INVs information directly
		print OUT "$line[0]\t$line[1]\t$line[2]\t$line[5]\t$line[6]\t$line[7]\tRegion_replacement\t$line[10]\t$line[8]\t$line[9]\n";
		$hash_ref{$line[0]}{$line[1]}{$line[2]}="$line[5]\t$line[6]\t$line[7]";
		$hash_query{$line[5]}{$line[6]}{$line[7]}="$line[0]\t$line[1]\t$line[2]";
	}elsif($line[10] eq "TRANS" or $line[10] eq "INVTR"){
		#Translocation and Inverted Translocation are StrucD, Convert TRANS/INVTR to Deletion and Insertion
		$hash_ref{$line[0]}{$line[1]}{$line[2]}="Struc_Deletion|$line[8]";
		$hash_query{$line[5]}{$line[6]}{$line[7]}="Struc_Insertion|$line[8]";
		$hash_trans{$line[8]}{"Combine"}="$line[0]\t$line[1]\t$line[2]\t-\t-\t-\tDeletion\t$line[10]\t$line[8]\t$line[9]";
		$hash_trans{$line[8]}{"No_Combine"}="$line[0]\t$line[1]\t$line[2]\t$line[5]\t$line[6]\t$line[7]\t$line[10]\t$line[10]\t$line[8]\t$line[9]";
	}elsif($line[10] eq "DEL" or $line[10] eq "INS" or $line[10] eq "SNP"){
		#DEL,INS,SNP of SYN are SequentialD, keep them directly
		#print OUT "$line[0]\t$line[1]\t$line[2]\t$line[5]\t$line[6]\t$line[7]\t$line[10]\t$line[10]\t$line[8]\n";
		$hash_syn_var{$line[8]}{'ref'}="$line[0]\t$line[1]\t$line[2]";
		$hash_syn_var{$line[8]}{'query'}="$line[5]\t$line[6]\t$line[7]";
		$hash_syn_var{$line[8]}{'info'}="$line[0]\t$line[1]\t$line[2]\t$line[5]\t$line[6]\t$line[7]\t$line[10]\t$line[10]\t$line[8]\t$line[9]";
	}elsif($line[10] eq "SYN"){
		$hash_ref{$line[0]}{$line[1]}{$line[2]}="$line[5]\t$line[6]\t$line[7]";
		$hash_query{$line[5]}{$line[6]}{$line[7]}="$line[0]\t$line[1]\t$line[2]";
		print OUT "$line[0]\t$line[1]\t$line[2]\t$line[5]\t$line[6]\t$line[7]\t$line[10]\t$line[10]\t$line[8]\t$line[9]\n";
	}
}

#明确起始位置的变异：CPG，CPL，copyloss，TDM，HDR，Translocation_Deletion,INV,NOTAL_Deletion
#需要推算其实位置的变异：Translocation_Insertion,copygain,NOTAL_Insertion
#推算都需要根据query的前后位置的比对情况

#output Structural Defference 
my $chr_now_ref="-";
my $chr_now_query="-";
my $start_now_ref=0;
my $start_now_query=0;
my $end_now_ref=0;
my $end_now_query=0;
my $end_flag="NO";
my @ID;
my %hash_trans_used;
foreach my $chr (keys %hash_query){
	if($chr ne $chr_now){
		$chr_now_query=$chr;
		$start_now_query=0;
		$end_now_query=0;
		$end_flag="NO";
	}
	foreach my $start  (sort {$a <=> $b} keys %{$hash_query{$chr}}){
		foreach my $end (sort {$b <=> $a} keys %{$hash_query{$chr}{$start}}){
			if($hash_query{$chr}{$start}{$end} =~ /Insertion/){
				#extend or start a new insertion
				my $tmp = (split /\|/,$hash_query{$chr}{$start}{$end})[-1];
				if($end_flag eq "INS"){
					#extend
					$end_now_query=$end;
					push(@ID,$tmp);
				}else{
					#start a new insertion
					$start_now_query=$start;
					$end_now_query=$end;
					$end_flag="INS";
					@ID=();
					push(@ID,$tmp);
				}
				next;
			}else{
				my $target = $hash_query{$chr}{$start}{$end};
				my ($chr_tmp,$start_tmp,$end_tmp)=split /\t/,$target;
				#end a insertion or extend syntenic region
				if($end_flag eq "INS"){
					#end a insertion
					#if this insertion is a first insert, attach it to followed systenic region  or  to the previous region
					if ($start_now_query==1){
						#attach to the followd region
						my $ID=join "|",@ID;
						if(@ID==1 and $ID =~/TR/){    #only one INS_type SV and it's a TRANS, then keep the raw SV Type : TRANS
							#next;
						}else{
						print OUT "$chr_tmp\t$start_tmp\t$start_tmp\t$chr\t$start_now_query\t$end_now_query\tInsertion\tCombine\t$ID\t-\n";   #or print OUT the merged SV and mark the related SV to be used to merge
						}
						$end_now_ref=$end_tmp;
						$chr_now_ref=$chr_tmp;
						if(@ID>1){
						    foreach my $tmp_ID (@ID){
							$hash_trans_used{$tmp_ID}=1;
						    }
						}
					}else{
						#attach to the previous region
						my $ID=join "|",@ID;
						if(@ID==1 and $ID =~ /TR/){
							#next;
						}else{
						print OUT "$chr_now_ref\t$end_now_ref\t$end_now_ref\t$chr\t$start_now_query\t$end_now_query\tInsertion\tCombine\t$ID\t-\n";
						}
						$end_now_ref=$end_tmp;
						$chr_now_ref=$chr_tmp;
						if(@ID>1){
						    foreach my $tmp_ID (@ID){
							$hash_trans_used{$tmp_ID}=1;
						    }
						}
					}
					$end_flag="SYN";
				}else{
					$end_now_ref=$end_tmp;
					$chr_now_ref=$chr_tmp;
					next;
				}
			}
		}
	}
}

foreach my $trans (keys %hash_trans){
	if($hash_trans_used{$trans}==1){
		print OUT "$hash_trans{$trans}{'Combine'}\n";
	}else{
		print OUT "$hash_trans{$trans}{'No_Combine'}\n";
	}
}

#output sequencial differences:
foreach my $sv (keys %hash_syn_var){
	my @ref=split /\t/,$hash_syn_var{$sv}{'ref'};
	my @query=split /\t/,$hash_syn_var{$sv}{'query'};
	$flag=0;
	foreach my $start (sort {$a <=> $b} keys %{$hash_replacement{$ref[0]}}){
		next if ($start < $ref[1]);
		my $end = $hash_replacement{$ref[0]}{$start};
		last if ($end > $ref[2]);
		$flag=1;
	}
	if($flag==0){
		foreach my $start (sort {$a <=> $b} keys %{$hash_replacement{$query[0]}}){
			next if ($start < $query[1]);
			my $end = $hash_replacement{$query[0]}{$start};
			last if ($end > $query[2]);
			$flag=1;	
		}
	}
	if($flag==0){
		print OUT "$hash_syn_var{$sv}{'info'}\n";
	}
}
close OUT;



