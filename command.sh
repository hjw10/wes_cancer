
#建文件夹
mkdir wes_cancer
cd wes_cancer
#建三个文件夹，软件用老师下好的
mkdir biosoft project data
cd project
mkdir -p 0.sra 1.raw_fq 2.clean_fq 3.qc/{raw_qc,clean_qc} 4.align/{qualimap,flagstats,stats} 5.gatk/gvcf 6.mutect 7.annotation/{vep,annovar,funcotator,snpeff} 8.cnv/{gatk,cnvkit,gistic,facet} 9.pyclone 10.signature

#有些分析软件服务器没有，我们自己搭个conda环境安装一下，用时补下
#所需软件：sra-tools fastqc trim-galore multiqc bwa samtools gnuplot qualimap subread vcftools bedtools cnvkit hcc GATK
#conda create -n wes python=3
#conda activate wes 
#conda install -y sra-tools fastqc trim-galore multiqc bwa samtools gnuplot qualimap subread vcftools bedtools cnvkit
#conda install -y -c hcc aspera-cli=3.7.7
#conda install bioconda::qualimap
#conda install bioconda/label/cf201901::qualimap

#下载GATK
cd ../biosoft
wget https://github.com/broadinstitute/gatk/releases/download/4.5.0.0/gatk-4.5.0.0.zip
unzip gatk-4.5.0.0.zip
#下载gatk资源文件
#google cloud下载，提前安装gdown软件（下载文件大小没有限制）
gdown -O ~/reference/gatk4/ https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta
gdown -O ~/reference/gatk4/ https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf
gdown -O ~/reference/gatk4/ https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/1000G_phase1.snps.high_confidence.hg38.vcf.gz
gdown -O ~/reference/gatk4/ https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/1000G_phase1.snps.high_confidence.hg38.vcf.gz.tbi
gdown -O ~/reference/gatk4/ https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
gdown -O ~/reference/gatk4/ https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi
gdown -O ~/reference/gatk4/ https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta.fai
gdown -O ~/reference/gatk4/ https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dict

#下载sra SRP070662下list下载过来
prefetch --option-file SRR_Acc_List.txt -O ./0.sra

#文件重命名
## 首先提取 sra 文件的 name
cat SraRunTable.txt | grep "WXS" | cut -d ',' -f 1 > sra
## 然后提取 sra 文件对应的样本 name
cat SraRunTable.txt | grep "WXS" | cut -d ',' -f 37 | sed s/_1// > config
## 合并，结果如下：两列，第一列为 sra 文件的 name，第二列为重命名后的样本 name
paste sra config > sra2case.txt
cat sra2case.txt
# SRR3182418 case5_germline
# SRR3182419 case2_germline
# SRR3182420 case4_germline
# SRR3182421 case3_germline
# SRR3182422 case6_germline
# SRR3182423 case1_germline
# SRR3182425 case4_biorep_B_techrep
# SRR3182426 case4_biorep_C
# SRR3182427 case3_biorep_A
# SRR3182428 case3_biorep_B
# ......
## 用下面代码重命名
cat sra2case.txt | while read name
do
 arr=(${name})
 sample=${arr[0]}
 case=${arr[1]}
 mv ./0.sra/${sample}/${sample}.sra ./0.sra/${case}.sra
done
## 最后删掉没用的文件夹


## srr转fastq
cd ./0.sra
cat ../config | while read id
do
time fasterq-dump -3 -e 16 ${id}.sra -O ../1.raw_fq --outfile ${id}.fastq >> ../1.raw_fq/${id}_sra2fq.log 2>&1
time pigz -p 10 -f ../1.raw_fq/${id}_1.fastq >../1.raw_fq/${id}_1.fastq.gz
time pigz -p 10 -f ../1.raw_fq/${id}_2.fastq >../1.raw_fq/${id}_2.fastq.gz
done
cd 1.raw_fq
time fastq-dump --gzip --split-files ../0.sra/*.sra

##fastqc
for name in `cat config`
do
fastqc --outdir ./3.qc/raw_qc/ --threads 16 \
./1.raw_fq/${name}*.fastq.gz >> ./3.qc/raw_qc/${name}_fastqc.log 2>&1
done

multiqc ./3.qc/raw_qc/*zip -o ./3.qc/raw_qc/multiqc

#质控、去接头
in_dir=1.raw_fq
out_dir=2.clean_fq
for name in `cat config`
do
    echo $name
	echo `date`
	echo "** Starting to clean $name **"
	time fastp -i $in_dir/${name}'_1.fastq.gz' -o $out_dir/${name}'_1.fastq.gz' \
        -I $in_dir/${name}'_2.fastq.gz' -O $out_dir/${name}'_2.fastq.gz' \
		-w 10
done

## trim_galore.sh
cat config | while read id
do
fq1=./1.raw_fq/${id}_1.fastq.gz
fq2=./1.raw_fq/${id}_2.fastq.gz
~/tools/TrimGalore-0.6.10/trim_galore \
--paired -q 28 --phred33 --length 30 --stringency 3 --gzip --cores 8 \
-o ./2.clean_fq $fq1 $fq2 >> ./2.clean_fq/${id}_trim.log 2>&1
done



##fastqc
for name in `cat config`
do
fastqc --outdir ./3.qc/clean_qc/ --threads 16 \
./2.clean_fq/${name}*.fq.gz >> ./3.qc/clean_qc/${name}_fastqc.log 2>&1
done

multiqc ./3.qc/clean_qc/*zip -o ./3.qc/clean_qc/multiqc

##序列比对
#构建参考基因组
time bwa index -a bwtsw -p gatk_hg38 ~/wes_cancer/data/Homo_sapiens_assembly38.fasta
#real    71m52.085s
#user    69m12.803s
#sys     1m55.061s

INDEX=~/wes_cancer/data/gatk_hg38
cat config | while read id
#id=case2_techrep_2
do
echo "start bwa for ${id}" `date`
fq1=./2.clean_fq/${id}_1_val_1.fq.gz
fq2=./2.clean_fq/${id}_2_val_2.fq.gz
bwa mem -M -t 16 -R "@RG\tID:${id}\tSM:${id}\tLB:WXS\tPL:Illumina" ${INDEX} ${fq1} ${fq2} -o ./4.align/${id}.sam 
samtools sort -@ 10 -m 1G ./4.align/${id}.sam -o ./4.align/${id}.bam 
echo "end bwa for ${id}" `date`
done

##比对后质控
for name in `sort config | cat`
do
    echo $name
	echo `date`
	echo "** Starting to stat $name **"
	samtools stats -@ 16 --reference ~/wes_cancer/data/Homo_sapiens_assembly38.fasta \
    ./4.align/${name}.bam \
	> ./4.align/stats/${name}.stat
	plot-bamstats -p ./4.align/stats/${name} ./4.align/stats/${name}.stat
done

##用qualimap查看基因组覆盖深度等信息
#下载参考文件
cd ~/reference
wget https://ftp.ncbi.nlm.nih.gov/pub/CCDS/current_human/CCDS.current.txt
cat CCDS.current.txt | grep "Public" | perl -alne '{/\[(.*?)\]/;next unless $1;$gene=$F[2];$exons=$1;$exons=~s/\s//g;$exons=~s/-/\t/g;print "$F[0]\t$_\t$gene" foreach split/,/,$exons;}' | sort -u | bedtools sort -i | awk '{if($3>$2) print "chr"$0"\t0\t+"}' > hg38.exon.bed
#会报错，原因是染色体名字不对，chr1 -> 1
#cat CCDS.current.txt | grep  "Public"| perl -alne '{/\[(.*?)\]/;next unless $1;$gene=$F[2];$exons=$1;$exons=~s/\s//g;$exons=~s/-/\t/g;print "$F[0]\t$_\t$gene" foreach split/,/,$exons;}' | sort -u | bedtools sort -i | awk '{print $0"\t0\t+"}'  > hg38.exon2.bed
for name in `sort config | cat | tail -25`
do
echo $name
qualimap bamqc --java-mem-size=10G -gff ~/reference/hg38.exon.bed -nr 100000 -nw 500 -nt 16 -bam ./4.align/${name}.bam -outdir ./4.align/qualimap/${name}
done

#GATK去重复
GATK=~/wes_cancer/biosoft/gatk-4.5.0.0/gatk
cat config | while read id
do
BAM=./4.align/${id}.bam
if [ ! -f ./5.gatk/ok.${id}_marked.status ]
then
echo "start MarkDuplicates for ${id}" `date`
gatk --java-options "-Xmx20G -Djava.io.tmpdir=./" MarkDuplicates \
-I ${BAM} \
--REMOVE_DUPLICATES true \
-O ./5.gatk/${id}_marked.bam \
-M ./5.gatk/${id}.metrics \
1>./5.gatk/${id}_log.mark 2>&1 
if [ $? -eq 0 ]
then
touch ./5.gatk/ok.${id}_marked.status
fi
echo "end MarkDuplicates for ${id}" `date`
samtools index -@ 16 -m 4G -b ./5.gatk/${id}_marked.bam ./5.gatk/${id}_mar
fi
done


##碱基质量重校正BQSR
GATK=~/wes_cancer/biosoft/gatk-4.5.0.0/gatk
snp=~/wes_cancer/data/1000G_phase1.snps.high_confidence.hg38.vcf.gz
indel=~/wes_cancer/data/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
ref=~/wes_cancer/data/Homo_sapiens_assembly38.fasta
cat config | while read id
do
if [ ! -f ./5.gatk/${id}_bqsr.bam ]
then
echo "start BQSR for ${id}" `date`
gatk --java-options "-Xmx20G -Djava.io.tmpdir=./" BaseRecalibrator \
-R $ref \
-I ./5.gatk/${id}_marked.bam \
--known-sites ${snp} \
--known-sites ${indel} \
-O ./5.gatk/${id}_recal.table \
1>./5.gatk/${id}_log.recal 2>&1 
gatk --java-options "-Xmx20G -Djava.io.tmpdir=./" ApplyBQSR \
-R $ref \
-I ./5.gatk/${id}_marked.bam \
-bqsr ./5.gatk/${id}_recal.table \
-O ./5.gatk/${id}_bqsr.bam \
1>./5.gatk/${id}_log.ApplyBQSR 2>&1 
echo "end BQSR for ${id}" `date`
fi
done

#报错后处理
#for name in `head -23 config`
#do 
#rm 5.gatk/${name}_bqsr.bam
#done

#单样本找 Germline SNPS + INDELs
GATK=~/wes_cancer/biosoft/gatk-4.5.0.0/gatk
snp=~/wes_cancer/data/1000G_phase1.snps.high_confidence.hg38.vcf.gz
indel=~/wes_cancer/data/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
ref=~/wes_cancer/data/Homo_sapiens_assembly38.fasta
bed=~/wes_cancer/data/hg38.exon.bed
cat config | while read id
do
echo "start HC for ${id}" `date`
gatk --java-options "-Xmx20G -Djava.io.tmpdir=./" HaplotypeCaller -ERC GVCF \
-R ${ref} \
-I ./5.gatk/${id}_bqsr.bam \
--dbsnp ${snp} \
-L ${bed} \
-O ./5.gatk/${id}_raw.vcf \
1>./5.gatk/${id}_log.HC 2>&1
echo "end HC for ${id}" `date`
done

cd ./5.gatk/gvcf
for chr in chr{1..22} chrX chrY chrM
do
time gatk --java-options "-Xmx20G -Djava.io.tmpdir=./" GenomicsDBImport \
-R ${ref} \
$(ls ../*raw.vcf | awk '{print "-V "$0" "}') \
-L ${chr} \
--genomicsdb-workspace-path gvcfs_${chr}.db
time gatk --java-options "-Xmx20G -Djava.io.tmpdir=./" GenotypeGVCFs \
-R ${ref} \
-V gendb://gvcfs_${chr}.db \
-O gvcfs_${chr}.vcf
done

gatk --java-options "-Xmx20G -Djava.io.tmpdir=./" GatherVcfs \
$(for i in {1..22} X Y M;do echo "-I gvcfs_chr${i}.vcf" ;done) \
-O merge.vcf

#mutect 找 Somatic mutations
GATK=~/wes_cancer/biosoft/gatk-4.5.0.0/gatk
ref=~/wes_cancer/data/Homo_sapiens_assembly38.fasta
bed=~/wes_cancer/data/hg38.exon.bed
cat config2 | while read id
do
arr=(${id})
sample=${arr[1]}
T=./5.gatk/${arr[1]}_bqsr.bam
N=./5.gatk/${arr[0]}_bqsr.bam
echo "start Mutect2 for ${id}" `date`
gatk --java-options "-Xmx20G -Djava.io.tmpdir=./" Mutect2 -R ${ref} \
-I ${T} -tumor $(basename "$T" _bqsr.bam) \
-I ${N} -normal $(basename "$N" _bqsr.bam) \
-L ${bed} \
-O ./6.mutect/${sample}_mutect2.vcf
gatk FilterMutectCalls \
 -R ${ref} \
-V ./6.mutect/${sample}_mutect2.vcf \
-O ./6.mutect/${sample}_somatic.vcf
echo "end Mutect2 for ${id}" `date`
cat ./6.mutect/${sample}_somatic.vcf | perl -alne '{if(/^#/){print}else{next unless $F[6] eq "PASS";next if $F[0] =~/_/;print } }' > ./6.mutect/${sample}_filter.vcf
done

#annovar
cat config2 | cut -d ' ' -f 2 > config3
cd annovar
./annotate_variation.pl -downdb -webfrom annovar gnomad_genome --buildver hg38 humandb/
./annotate_variation.pl -downdb -webfrom annovar refGene -buildver hg38 humandb/
./annotate_variation.pl -downdb -webfrom annovar knownGene -buildver hg38 humandb/
./annotate_variation.pl -downdb -webfrom annovar exac03 -buildver hg38 humandb/
./annotate_variation.pl -downdb -webfrom annovar cosmic70 -buildver hg38 humandb/
./annotate_variation.pl -downdb -webfrom annovar nci60 -buildver hg38 humandb/
./annotate_variation.pl -downdb -webfrom annovar avsnp150 -buildver hg38 humandb/
./annotate_variation.pl -downdb -webfrom annovar 1000g2015aug -buildver hg38 humandb/
./annotate_variation.pl -downdb -webfrom annovar clinvar_20221231 -buildver hg38 humandb/

./annotate_variation.pl -downdb -buildver hg38 -webfrom annovar avdblist hg38_list/
wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/funcotator//funcotator_dataSources.v1.8.hg38.20230908s.tar.gz

cat config3 | while read id
do
echo "start ANNOVAR for ${id} " `date`
~/tools/annovar/table_annovar.pl ./6.mutect/${id}_filter.vcf ~/tools/annovar/humandb/ \
-buildver hg38 \
-out ./7.annotation/annovar/${id} \
-remove \
-protocol refGene,knownGene,clinvar_20221231 \
-operation g,g,f \
-nastring . \
-vcfinput
echo "end ANNOVAR for ${id} " `date`
done


#gatk funcatator
wget -c ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/funcotator//funcotator_dataSources.v1.8.hg38.20230908s.tar.gz
tar -zxvf 

bed=~/wes_cancer/data/hg38.exon.bed
dict=~/wes_cancer/data/Homo_sapiens_assembly38.dict
gatk BedToIntervalList -I ${bed} -O hg38.exon.interval_list -SD ${dict}

GATK=~/wes_cancer/biosoft/gatk-4.5.0.0/gatk
ref=~/wes_cancer/data/Homo_sapiens_assembly38.fasta
interval=~/wes_cancer/data//hg38.exon.interval_list
source=~/wes_cancer/data/funcotator_dataSources.v1.8.hg38.20230908s
cat config3 | while read id
do
echo "start Funcotator for ${id} " `date`
gatk IndexFeatureFile --input ./6.mutect/${id}_filter.vcf
gatk Funcotator -R $ref \
-V ./6.mutect/${id}_filter.vcf \
-O ./7.annotation/funcotator/${id}_funcotator.tmp.maf \
--data-sources-path ${source} \
--intervals ${interval} \
--output-file-format MAF \
--ref-version hg38 
echo "end Funcotator for ${id} " `date`
done


#vep
vep_install -a cf -s homo_sapiens -y GRCh38 -c ~/wes_cancer/data/
curl -O https://ftp.ensembl.org/pub/release-105/variation/indexed_vep_cache/homo_sapiens_vep_105_GRCh38.tar.gz

cat ./config3 | while read id
do
echo "start vep_annotation for ${id} " `date`
vep --cache --format vcf --vcf --force_overwrite \
--everything \
--species homo_sapiens \
--dir ~/wes_cancer/data/vep \
--af  --af_gnomad --af_exac \
--fasta ../data/Homo_sapiens_assembly38.fasta \
-i ./6.mutect/${id}_filter.vcf \
-o ./7.annotation/vep/${id}_vep.vcf 
echo "end vep_annotation for ${id} " `date`
done

filter_vep -i 7.annotation/vep/case1_biorep_A_techrep_vep.vcf \
--format vcf --filter "AF < 0.005 or not AF" \
-o 7.annotation/vep/case1_biorep_A_techrep_vep_filtered.vcf


#snpeff
wget -c https://snpeff.blob.core.windows.net/versions/snpEff_latest_core.zip
unzip 
java -jar snpEff.jar download -v GRCh38.99

dir=~/wes_cancer/project/7.annotation/
config=~/wes_cancer/project/config3
cat ${config} | while read id
do
echo "start snpeff_annotation for ${id} " `date`
java -Xmx4g -jar ../biosoft/snpEff/snpEff.jar -v \
-cancer -cancerSamples ./7.annotation/snpeff/samples_cancer.txt \
GRCh38.86 \
./6.mutect/${id}_filter.vcf > ./7.annotation/snpeff/${id}_eff.vcf 
echo "end snpeff_annotation for ${id} " `date`
done

#不同 vcf 注释结果转 maf
#ANNOVAR 注释后拿到的结果，还需要做一个处理，就是在每一个 ${id}.hg38_multianno.txt 文件中增加两列，为 Tumor_Sample_Barcode 和 Matched_Norm_Sample_Barcode ，而最
#后的一列 Otherinfo 我们暂时丢弃掉（如果你需要就不要丢弃），添加上这两列信息之后再把所有的 vcf 文件合并成一个 annovar_merge.vcf 。处理的方法如下
cat config3 | while read id
do
grep -v '^Chr' ./7.annotation/annovar/${id}.hg38_multianno.txt | cut -f 1-20 | awk -v T=${id} -v N=${id:0:5}_germline '{print $0"\t"T"\t"N}' >./7.annotation/annovar/${id}.annovar.vcf
done
head -1 ./7.annotation/annovar/case1_biorep_A_techrep.hg38_multianno.txt | sed 's/Otherinfo/Tumor_Sample_Barcode\tMatched_Norm_Sample_Barcode/' >./7.annotation/annovar/head
cat ./7.annotation/annovar/head ./7.annotation/annovar/*annovar.vcf >./7.annotation/annovar/annovar_merge.vcf
#R

#GATK
## 添加Barcode
cat config3 | while read id
do 
	grep -v '^#' ./7.annotation/funcotator/${id}_funcotator.tmp.maf| grep -v '^Hugo_Symbol'| awk -v T=${id} -v N=${id:0:5}_germline 'BEGIN{FS="\t";OFS="\t"}{$16=T;$17=N;print $0}' >./7.annotation/funcotator/${id}_funcotator.maf 
done
## 取出一个列名
grep 'Hugo_Symbol' ./7.annotation/funcotator/case1_biorep_A_techrep_funcotator.tmp.maf >./7.annotation/funcotator/header
## 删除掉旧文件
rm ./7.annotation/funcotator/*tmp.maf
## 合并所有的样本的maf
cat ./7.annotation/funcotator/header ./7.annotation/funcotator/*maf >./7.annotation/funcotator/funcotator_merge.maf

#vep
export VCF2MAF_URL=`curl -sL https://api.github.com/repos/mskcc/vcf2maf/releases | grep -m1 tarball_url | cut -d\" -f4`
curl -L -o mskcc-vcf2maf.tar.gz $VCF2MAF_URL
tar -zxf mskcc-vcf2maf.tar.gz
cd mskcc-vcf2maf-*

cat config3 | while read id
do
	sample=./7.annotation/vep/${id}_vep.vcf
	normal=${id%%_*}_germline
	echo ${sample}
	echo ${normal}
	echo ${id}
	perl ~/wes_cancer/biosoft/mskcc-vcf2maf-754d68a/vcf2maf.pl \
	--input-vcf  ${sample}   \
	--output-maf ./7.annotation/vep/${id}_vep.maf  \
	--ref-fasta ~/wes_cancer/data/Homo_sapiens_assembly38.fasta \
	--tumor-id  ${id}  \
	--normal-id ${normal}  \
	--vep-path ~/anaconda3/envs/vep/bin \
	--vep-data ~/wes_cancer/data/vep \
	--ncbi-build GRCh38
done

cat ./7.annotation/vep/*maf | grep -v '^#'| grep -v '^Hugo_Symbol' >./7.annotation/vep/tmp 
grep 'Hugo_Symbol' ./7.annotation/vep/case1_biorep_A_techrep_vep.maf >./7.annotation/vep/header
cat ./7.annotation/vep/header ./7.annotation/vep/tmp >./7.annotation/vep/vep_merge.maf


#snpeff

#外显子坐标的interval文件
ref=~/wes_cancer/data/Homo_sapiens_assembly38.fasta
bed=~/wes_cancer/data/hg38.exon.bed
dict=~/wes_cancer/data/Homo_sapiens_assembly38.dict

## bed to intervals_list
gatk BedToIntervalList -I ${bed} -O ~/wes_cancer/data/hg38.exon.interval_list -SD ${dict}

## Preprocess Intervals
gatk  PreprocessIntervals \
-L ~/wes_cancer/data/hg38.exon.interval_list \
--sequence-dictionary ${dict} \
--reference ${ref}  \
--padding 250 \
--bin-length 0 \
--interval-merging-rule OVERLAPPING_ONLY \
--output ~/wes_cancer/data/targets.preprocessed.interval_list


#获取样本的read counts
GENOME=~/wes_cancer/data/Homo_sapiens_assembly38.fasta
dict=~/wes_cancer/data/Homo_sapiens_assembly38.dict
INDEX=~/wes_cancer/data/bwa_index/gatk_hg38
DBSNP=~/wes_cancer/data/Homo_sapiens_assembly38.dbsnp138.vcf.gz
kgSNP=~/wes_cancer/data/1000G_phase1.snps.high_confidence.hg38.vcf.gz
kgINDEL=~/wes_cancer/data/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
interval=~/wes_cancer/data/targets.preprocessed.interval_list
cat ~/wes_cancer/project/config | while read id
do
	i=~/wes_cancer/project/5.gatk/${id}_bqsr.bam
	echo ${i}
	## step1 : CollectReadCounts
	time gatk  --java-options "-Xmx20G -Djava.io.tmpdir=./" 	CollectReadCounts \
	-I ${i} \
	-L ${interval} \
	-R ${GENOME} \
	--format HDF5  \
	--interval-merging-rule OVERLAPPING_ONLY \
	--output ${id}.clean_counts.hdf5
	## step2 : CollectAllelicCounts
	time gatk  --java-options "-Xmx20G -Djava.io.tmpdir=./"  CollectAllelicCounts \
	-I ${i} \
	-L ${interval} \
	-R ${GENOME} \
	-O ${id}.allelicCounts.tsv
done


mkdir allelicCounts
mv *.allelicCounts.tsv ./allelicCounts
mkdir counts
mv *.clean_counts.hdf5  ./counts
##################################################
# 接着合并所有的normal样本的数据创建 cnvponM.pon.hdf5 #
##################################################

gatk  --java-options "-Xmx20g" CreateReadCountPanelOfNormals \
--minimum-interval-median-percentile 5.0 \
--output cnvponM.pon.hdf5 \
--input counts/case1_germline.clean_counts.hdf5 \
--input counts/case2_germline.clean_counts.hdf5 \
--input counts/case3_germline.clean_counts.hdf5 \
--input counts/case4_germline.clean_counts.hdf5 \
--input counts/case5_germline.clean_counts.hdf5 \
--input counts/case6_germline.clean_counts.hdf5
 

############################################
############# 最后走真正的CNV流程 #############
############################################

cat config | while read id
do
	i=./counts/${id}.clean_counts.hdf5
	gatk  --java-options "-Xmx20g" DenoiseReadCounts \
	-I $i \
	--count-panel-of-normals cnvponM.pon.hdf5 \
	--standardized-copy-ratios ${id}.clean.standardizedCR.tsv \
	--denoised-copy-ratios ${id}.clean.denoisedCR.tsv
done

mkdir denoisedCR standardizedCR segments cnv_plots
mv *denoisedCR.tsv ./denoisedCR
mv *standardizedCR.tsv ./standardizedCR

cat config | while read id
do
	i=./denoisedCR/${id}.clean.denoisedCR.tsv
	## ModelSegments的时候有两个策略，是否利用CollectAllelicCounts的结果
	gatk   --java-options "-Xmx20g" ModelSegments \
	--denoised-copy-ratios $i \
	--output segments \
	--output-prefix ${id}
done
	## 如果要利用CollectAllelicCounts的结果就需要增加两个参数，这里就不讲解了。
cat config | while read id
do
    i=./denoisedCR/${id}.clean.denoisedCR.tsv
	gatk   --java-options "-Xmx20g" CallCopyRatioSegments \
	-I segments/${id}.cr.seg \
	-O segments/${id}.clean.called.seg
done

##非常重要！！！！把R版本搞到4.0！！！
cat config | while read id
do
	## 这里面有两个绘图函数，PlotDenoisedCopyRatios 和 PlotModeledSegments ，可以选择性运行。
    i=./denoisedCR/${id}.clean.denoisedCR.tsv
	gatk   --java-options "-Xmx20g" PlotDenoisedCopyRatios \
	--standardized-copy-ratios 	./standardizedCR/${id}.clean.standardizedCR.tsv \
	--denoised-copy-ratios $i \
	--sequence-dictionary ${dict} \
	--output-prefix ${id} \
	--output cnv_plots 
done

cat config | while read id
do
    i=./denoisedCR/${id}.clean.denoisedCR.tsv
	gatk   --java-options "-Xmx20g" PlotModeledSegments \
	--denoised-copy-ratios $i \
	--segments segments/${id}.modelFinal.seg \
	--sequence-dictionary ${dict} \
	-O cnv_plots \
	--output-prefix ${id}.clean
done





#17
#pyclone
cat config3 | while read id
do 
cat ./7.annotation/vep/${id}_vep.maf | sed '1,2d' | awk -F '\t' '{print $5"\t"$6"\t"$7"\t"$5":"$6":"$1"\t"$41"\t"$42"\t"2"\t"0"\t"2"\t"}' >./9.pyclone/${id}.tmp.tsv
done

head ./9.pyclone/case1_biorep_A_techrep.tmp.tsv

cat config3 | while read id
do 
cat ./8.cnv/gatk/segments/${id}.cr.igv.seg | sed '1d' | awk 'BEGIN{OFS="\t"}{print $0"\t"int((2^$6)*2+0.5)}'| awk 'BEGIN{OFS="\t"}{if ($7!=0)print $0}' | cut -f 2-7  >./9.pyclone/${id}.bed
done
## 生成的文件如下
head ./9.pyclone/case1_biorep_A_techrep.bed

cat config3 | while read id;
do 
bedtools window -a ./9.pyclone/${id}.tmp.tsv  -b ./9.pyclone/${id}.bed | cut -f 4-8,15 | awk 'BEGIN{OFS="\t";print "mutation_id\tref_counts\tvar_counts\tnormal_cn\tminor_cn\tmajor_cn"}{print $0}' >./9.pyclone/${id}.tsv
done


#pyclone
for i in case{1..6}
do
PyClone run_analysis_pipeline --prior total_copy_number --in_files ./9.pyclone/${i}*tsv --working_dir ./9.pyclone/${i}_pyclone_analysis 1>./9.pyclone/${i}.log 2>&1
done

tree ./9.pyclone/case1_pyclone_analysis/
head ./9.pyclone/case1_pyclone_analysis/tables/loci.tsv
cat  ./9.pyclone/case1_pyclone_analysis/tables/loci.tsv | sed '1d' | cut -f 1 | sort | uniq -c
cat  ./9.pyclone/case1_pyclone_analysis/tables/loci.tsv | sed '/mutation_id/d' | cut -f 1|sort |uniq -c |wc -l



#18 肿瘤进化分析
conda config --add channels http://conda.anaconda.org/dranew
conda config --add channels https://conda.anaconda.org/IBMDecisionOptimization/linux-64

conda install citup

for case in case{1..6}
do
cat ./9.pyclone/${case}_pyclone_analysis/tables/loci.tsv | cut -f 6 | sed '1d' | paste - - - -  >./9.pyclone/${case}_pyclone_analysis/freq.txt
cat ./9.pyclone/${case}_pyclone_analysis/tables/loci.tsv | cut -f 3 | sed '1d' | paste - - - - |cut -f 1 >./9.pyclone/${case}_pyclone_analysis/cluster.txt
## 获取 sample_id 信息，后面画图要用到
cat ./9.pyclone/${case}_pyclone_analysis/tables/loci.tsv | cut -f 2 | sed '1d' | head -4 >./9.pyclone/${case}_pyclone_analysis/sample_id
done

for case in case{1..6}
do
run_citup_qip.py --submit local \
./9.pyclone/${case}_pyclone_analysis/freq.txt \
./9.pyclone/${case}_pyclone_analysis/cluster.txt \
./9.pyclone/${case}_pyclone_analysis/results.h5
done

for case in case{1..6}
do
python citup.py ./9.pyclone/${case}_pyclone_analysis/results.h5 | sed 's/^ \[//;s/\[//g;s/\]//g' | tr ' ' '\t'| grep '\.' >./9.pyclone/${case}_pyclone_analysis/cellfreq.txt
python citup.py ./9.pyclone/${case}_pyclone_analysis/results.h5 | sed 's/^ \[//;s/\[//g;s/\]//g' | tr ' ' '\t'| grep -v '\.' >./9.pyclone/${case}_pyclone_analysis/tree.txt
done


