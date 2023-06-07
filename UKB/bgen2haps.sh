export PATH=$PATH:/dssg/home/acct-bmelgn/bmelgn-3/FLAT/code/
state=$1
str=(${state// / })
i=${str[3]}

#if [ ! -f "$file" ]; then
  echo $state
  range=${str[0]}:${str[1]}-${str[2]}
  chr=${str[0]}
  CHR="${chr//"chr"}"
  mkdir -p $i
  cd $i

	samtools faidx ~/FLAT/data/$chr.fa $range > int.fa
	samtools faidx int.fa
  #picard CreateSequenceDictionary -R int.fa
  bgenix -g /dssg/home/acct-bmelgn/share/ukb22828_c${CHR}_b0_v3.bgen -incl-rsids ../../int/loci$CHR/$i > int.bgen  
	plink2 --bgen int.bgen ref-first --sample /dssg/home/acct-bmelgn/share/ukb22828_c${CHR}_b0_v3.sample --hwe 1e-6 --mac 10 --max-alleles 2 --geno 0.1 --mind 0.1 --recode vcf --out int --output-chr chrM
  CrossMap.py  vcf  ~/FLAT/data/hg19ToHg38.over.chain  int.vcf  ~/FLAT/data/$chr.fa  lift.vcf
  #picard LiftoverVcf -I int.vcf -O lift.vcf -C ~/FLAT/data/hg19ToHg38.over.chain -R ~/FLAT/data/$chr.fa --REJECT rej.vcf 
  bgzip lift.vcf
  bcftools index lift.vcf.gz
  if shapeit4 --input lift.vcf.gz --map /dssg/home/acct-bmelgn/bmelgn-3/FLAT/data/$chr.b38.gmap.gz --output phased.vcf --region $range ;then
    plink2 --vcf phased.vcf --export haps ref-first --out geno
  else
    bgzip -d lift.vcf.gz
    sed -i 's/\//|/g' lift.vcf
    sed -i 's/\.|./0|0/g' lift.vcf
    plink2 --vcf lift.vcf --export haps ref-first --out geno
  fi  
	awk 'BEGIN{OFS="\t"}{print "chr"$1,$3,$2,$4,$5,".",".",".","GT"}' geno.haps | awk 'BEGIN{print "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT"}1'  > variant.txt
	cut -d " " -f 6- geno.haps | datamash transpose -W > int.txt
  cut -s -f 2 -d " " geno.sample | tail -n +3 | awk -F"_" '{for(i=0;i<2;i++)print$1}' | paste -d"," - int.txt > int1.txt
  sort int1.txt -k 2 -t "," > int.txt
	awk -F"," '!_[$2]++' int.txt | awk -F"," '{$0=$2","NR} 1' >haps.txt
	join -1 2 -2 1 -t "," -o1.1 -o2.2 int.txt haps.txt | sort -k 1 -t "," | datamash -t"," -g1 collapse 2 > haps.geno
	sed -i 's/,/\t/g' haps.txt
  nr=$(wc -l haps.txt)
  nr=(${nr// / })
  nr=${nr[0]}  
  yes $i | head -n $nr | nl > hap.index
     
	awk 'BEGIN{OFS="\t";} NF{NF--};1' haps.txt | nl | datamash -W transpose | paste variant.txt - | awk 'BEGIN{print "##fileformat=VCFv4.2"}1' | bgzip > haplotype.vcf.gz
  bcftools index haplotype.vcf.gz

  cp int.fa ref.fa
  rm *txt
  rm geno*
  rm int*
  rm lift*
  rm phased.vcf  
  #awk '/^>/{print ">" substr(FILENAME,1,length(FILENAME)-3); next} 1' *.fa
  #rm 1_*
#fi  
#fi

