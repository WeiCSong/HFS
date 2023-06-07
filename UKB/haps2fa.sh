export PATH=$PATH:/dssg/home/acct-bmelgn/bmelgn-3/FLAT/code/
stat=$1

state="${stat//"hapchunk/"}"
fasta=""
mkdir $state.int
for NUM in {1..30000}; do
  line=`sed "${NUM}q;d" $stat `
  str=(${line// / })
  j=${str[0]}
  i=${str[1]}
  if [[ "$j" == 1 ]]; then
   cat $i/ref.fa | tail -n +2 | awk -v header=$i.int 'BEGIN{print ">"header}1' > $state.int/$i.int.fa 
   else
    colnum=$((j+8))
    zcat $i/haplotype.vcf.gz | awk -v col="$colnum" -v OFS="\t" '{ if ($col == 1) print $1,$2,$3,$4,$5,$6,$7,$8 }'  | awk 'BEGIN{print "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"}1' | awk 'BEGIN{print "##fileformat=VCFv4.2"}1'  | bgzip > $i/$j.vcf.gz
    bcftools index $i/$j.vcf.gz 
    bcftools consensus -f $i/ref.fa $i/$j.vcf.gz | tail -n +2 | awk -v header=$i.$j 'BEGIN{print ">"header}1' > $state.int/$i.$j.fa
    rm $i/$j.vcf*
  fi 
  

  echo $line
done

cat $state.int/*fa | gzip > $state.fa.gz
rm -r $state.int/

#echo $fasta >$state.fa
#sed -i 's/ /\n/g' $state.fa
#gzip $state.fa

#| awk -v header=">$i.$j" 'BEGIN{print header}1'
  #  mv $i/$i.$j.fa $state.int/$i.$j.fa

  
  #  if [[ "$j" == 1 ]]; then
  #  new=`cat $i/ref.fa | tail -n +2 | awk -v header=$i.int 'BEGIN{print ">"header}1' `
  #else
  #  new=`bcftools consensus -f $i/ref.fa -s $j -H 1 $i/haplotype.vcf.gz | tail -n +2 | awk -v header=$i.$j 'BEGIN{print ">"header}1' `
  #fi
  #fasta+=" $new"
