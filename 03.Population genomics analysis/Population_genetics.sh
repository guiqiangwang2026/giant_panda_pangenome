vg giraffe -p -t 64 -Z  ref.d2.gbz -d ref.d2.dist -m ref.d2.min \
-f /fq/s1_1.fq.gz -f /fq/s1_2.fq.gz -N s1 > s1.gam
vg pack -x ref.d2.gbz -g s1.sorted.gam -Q 5 -s 5 -o s1.pack
vg call s1.pack -a -t 16 > s1.vcf

bcftools view -i 'F_MISSING <= 0.2 & MAF > 0.05 & N_ALT = 1 & TYPE="snp" & F_PASS(FORMAT/DP >=5 & GT!="mis") >0.8' --threads 30 Panda.vcf.gz | bcftools annotate -x INFO,FORMAT -Oz -o Panda.snp.vcf.gz
bcftools index Panda.snp.vcf.gz
bcftools view -i 'F_MISSING <= 0.2 & MAF > 0.05 & N_ALT = 1 & TYPE="indel" & F_PASS(FORMAT/DP >=5 & GT!="mis") >0.8' --threads 30 Panda.vcf.gz | bcftools annotate -x INFO,FORMAT -Oz -o Panda.indel.vcf.gz
bcftools index Panda.indel.vcf.gz

#pi
vcftools --gzvcf vcf.gz --window-pi 100000 --out result

#Het
vcftools --gzvcf vcf.gz --het --out result

#Fst
vcftools --gzvcf vcf.gz --weir-fst-pop pop1.list --weir-fst-pop pop2.list --fst-window-size 50000 --fst-window-step 10000 --out --fst-window-size 50000 --fst-window-step 10000 --out 


#PCA 
vcftools  --gzvcf vcf.gz  --plink --out Panda
plink --noweb --file Panda --make-bed --out Panda.bfile
plink --bfile Panda.bfile --pca 10 --out pca


#tree
iqtree -s Panda.min4.phy -nt 36 -bb 1000 -alrt 1000  -m MFP

#admixture
admixture -j10 --cv Panda.bfile.bed n

#IBD
java -jar -Xmx50g -Xms32g -Xss5m  refined-ibd.jar gt=vcf.gz out=Panda map=map.gz

#ROH
plink --allow-extra-chr --bfile result --homozyg --homozyg-window-snp 20 --homozyg-kb 100 --homozyg-density 50 --out ROH

#mutation load
snpEff ann -c snpEff.config Panda.v1 vcf.gz > Panda.snpeff

