# Rice SV Analysis


#wd = "/public-dss/Project/Project_hwlu/hwlu/Core_set/SV/SV_pipline"


## 1. data required:

### 1.1. genome fasta files :

   ```Bash
   # named as $spe.genome
   ```


### 1.2. trf annotation result in gff3 format produced by command like:

   ```Bash
   trf genome.fasta 2 7 7 80 10 50 500 -f -d -m
   perl trfRepeatMask.pl
   
   # named as $spe.trf 
   ```


 

### 1.3. TE annotation result in gff3 format produced by command like : 

   ```Bash
   RepeatMasker -pa 36 -q -no_is -norna -nolow -div 40 -lib $TE_library -cutoff 225 -gff genome.fasta
   
   # named as $spe.genome.out.gff 
   ```



### 1.4. Syri result produced by commands like :

 ```Bash
 nucmer --mum -l 50 -c 100 ref.fa query.fa
 delta-filter -m -i 90 $delta > $delta.filtered
 show-coords -THrd $delta.filtered > $delta.filtered.coords 
 syri  -c $delta.filtered.coords -r $ref -q $query -d $delta --all --no-chrmatch
 
 # named by $spe.syri.out 
 ```



## 2. processing syri out for individuals ：

 ```Bash
 input = $spe.syri.out
 perl SyRI_Parse.pl $input
 perl Covert_CPV.pl $input.SV.bed
 ```


## 3. merge results for all individuals:

 ```Bash
 for i in $(*.syri.out.SV.bed.mergeDel.CPV2InDel);do awk '{len1=0;if($7=="Insertion"){len1=$6-$5+1}else{len1=$3-$2+1};if(len1>=50 || len1<=-50){print}}' $i >$i.50bp;done
 
 perl Merge_SyRI_SV_spe.pl 
 ```



## 4. Infer SV Type (Derived SV) :

### 4.1 rename SV ID:

   ```纯文本
   awk 'BEGIN{OFS="\t"}{$2="SV"NR;print $0}'   Merge_Syri_bed.result.0.9.0.9 > Merge_Syri_bed.result.0.9.0.9.rename
   mv  Merge_Syri_bed.result.0.9.0.9.rename  Merge_Syri_bed.result.0.9.0.9
   
   ```



### 4.2 add spe column

   ```Bash
   # out file is : Merge_Syri_bed.result.0.9.0.9.tag
   perl Add_spe_column.pl Merge_Syri_bed.result.0.9.0.9 
   ```



### 4.3 Infer SV derived type 

   ```Bash
   perl Infer_1.pl
   perl Intersect.pl
   perl Infer_2.pl 
   ```



## 5. infer SV mechanism:

   ```Bash
   perl ChangeSV_sampleID.pl
   awk '{if($1!~/Undefined/){print}}'  Merge_Syri_bed.result.0.9.0.9.rename.tag.infer_1.infer_2.Add >Merge_Syri_bed.result.0.9.0.9.rename.tag.infer_1.infer_2.Add.Defined
    
   #split SV for VNTR :
   perl split_Merged_SV_For_VNTR.pl 
   
   #infer SV for spe: 
   for i in $(ls *.syri.out.SV.bed.mergeDel.CPV2InDel);do perl infer_SV_for_mechanism.pl $i Merge_Syri_bed.result.0.9.0.9.tag.infer_1.infer_2.Add.vntr ;done 
   
   # filter ambiguous SVs
   for i in $(ls *.syri.out.SV.bed.mergeDel.CPV2InDel.infered);do perl Filter_overlaped_SV.pl  $i;done
   
   # Check TE and tandem repeat : 
   for i in $(ls *.syri.out.SV.bed.mergeDel.CPV2InDel.infered.uniq);do perl Combine_TE_TRF_SV.pl  $i;done
   
   # call SV cluster type :
   for i in $(ls *.syri.out.SV.bed.mergeDel.CPV2InDel.infered.uniq);do a=${i%%\.*};perl Call_SV_Mechanism1.pl  $i $a 1 ; done
   
   # Check SV breakpoint homology: 
   for i in $(ls *.syri.out.SV.bed.mergeDel.CPV2InDel.infered.uniq.TRF_TE);do perl Check_Homology.pl $i;done
   
   # call SV mechanism :
   for i in $(ls *.syri.out.SV.bed.mergeDel.CPV2InDel.infered.uniq.TRF_TE.homo);do a=${i%.syri.out.SV.bed.mergeDel.CPV2InDel.infered.uniq.TRF_TE.homo};echo $a;perl Call_SV_Mechanism_new.pl   $i $a.syri.out.SV.bed.mergeDel.CPV2InDel.infered.uniq.ref.Cluster_Type $a.syri.out.SV.bed.mergeDel.CPV2InDel.infered.uniq.query.Cluster_Type;done
   ```




