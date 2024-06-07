##annovar
rm(list = ls())
require(maftools)
options(stringsAsFactors = F)
## annovar
annovar.laml <- annovarToMaf(annovar = "./7.annotation/annovar/annovar_merge.vcf", 
                     refBuild = 'hg38',
                     tsbCol = 'Tumor_Sample_Barcode', 
                     table = 'refGene',
                     MAFobj = T)

## 统计一下case1_biorep_A_techrep样本的突变类型
table(annovar.laml@data[annovar.laml@data$Tumor_Sample_Barcode=='case1_biorep_A_techrep',]$Func.knownGene)
## exonic intergenic 
##   172          1

table(annovar.laml@data[annovar.laml@data$Tumor_Sample_Barcode=='case1_biorep_A_techrep',]$Variant_Type)
## DEL INS SNP 
##  5   0 168

table(annovar.laml@data[annovar.laml@data$Tumor_Sample_Barcode=='case1_biorep_A_techrep',]$Variant_Classification)
##   Frame_Shift_Del 
##                 4 
##   Frame_Shift_Ins 
##                 0 
##      In_Frame_Del 
##                 1 
##      In_Frame_Ins 
##                 1 
## Missense_Mutation 
##               153 
## Nonsense_Mutation 
##                14 
##  Nonstop_Mutation 
##                 0

rm(list = ls())
require(maftools)
options(stringsAsFactors = F)
## annovar
annovar.laml <- annovarToMaf(annovar = "./7.annotation/annovar/annovar_merge.vcf",
 refBuild = 'hg38',
 tsbCol = 'Tumor_Sample_Barcode', 
 table = 'refGene',
 MAFobj = T)
## gatk
# gatk.laml = read.maf(maf = 'gatk/gatk4.1.4.0_merge.maf')
library(data.table)
tmp=fread('./7.annotation/funcotator/funcotator_merge.maf')
gatk.laml = read.maf(maf = tmp)
## vep
vep.laml = read.maf(maf = './7.annotation/vep/vep_merge.maf')
## for vep.laml
library(stringr)
vep.laml@data$Protein_Change = paste0("p.",
 str_sub(vep.laml@data$Amino_acids,1,1),
 vep.laml@data$Protein_position,
 str_sub(vep.laml@data$Amino_acids,3,3))
## save Rdata
save(annovar.laml, gatk.laml, vep.laml, file = 'laml.Rdata')
## Summary
laml=annovar.laml
unique(laml@data$Tumor_Sample_Barcode)
getSampleSummary(laml)
getGeneSummary(laml)
getFields(laml)
laml=gatk.laml
unique(laml@data$Tumor_Sample_Barcode)
getSampleSummary(laml)
getGeneSummary(laml)
getFields(laml)
laml=vep.laml
unique(laml@data$Tumor_Sample_Barcode)
getSampleSummary(laml)
getGeneSummary(laml)
getFields(laml)


rm(list = ls())
require(maftools)
options(stringsAsFactors = F)
load(file = 'laml.Rdata')
laml=c(annovar.laml,gatk.laml,vep.laml)
## mafsummary
anno=c('annovar','gatk','vep')
for (i in 1:3) {
 #i = 1
 png(paste0('plotmafSummary_', anno[i], '.png'),res = 150,width = 1080,height = 1080)
 plotmafSummary(maf = laml[[i]],
 rmOutlier = TRUE,
 showBarcodes = T,
 textSize = 0.4,
 addStat = 'median',
 dashboard = TRUE,
 titvRaw = FALSE)
 dev.off()
}

## oncoplot_top30
for (i in 1:3) {
 #i = 1
 png(paste0('oncoplot_top30_', anno[i], '.png'),res = 150,width = 1080,height = 1080)
 oncoplot(maf = laml[[i]],
 top = 30,
 fontSize = 0.5,
 sampleOrder = laml[[i]]@clinical.data$Tumor_Sample_Barcode,
 showTumorSampleBarcodes = T)
 dev.off()
}

## lollipopPlot for SP140
gene='SP140'
protein=c("AAChange.refGene","Protein_Change","Protein_Change")
for (i in 1:3) {
 #i=3
 png(paste0(gene,'_', anno[i], '.png'),res = 150,width = 1080,height = 1080)
 maftools::lollipopPlot(maf = laml[[i]],
 gene = gene,
 AACol = protein[i],
 labelPos = 'all')
 dev.off()
}

## lollipopPlot for TP53
gene='TP53'
protein=c("AAChange.refGene","Protein_Change","Protein_Change")
for (i in 1:3) {
 #i=3
 png(paste0(gene,'_', anno[i], '.png'),res = 150,width = 1080,height = 1080)
 maftools::lollipopPlot(maf = laml[[i]],
 gene = gene,
 AACol = protein[i],
 labelPos = 'all')
 dev.off()
}



#13
#突变 Somatic Signature 图谱
rm(list=ls())
options(stringsAsFactors=FALSE)
## 切换镜像
options(BioC_mirror="https://mirrors.ustc.edu.cn/bioc/")
options("repos" = c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
## 安装R包
install.packages('deconstructSigs')
BiocManager::install('BSgenome')
BiocManager::install('BSgenome.Hsapiens.UCSC.hg38')
library(deconstructSigs)
## https://github.com/raerose01/deconstructSigs

## 读入数据
maf=read.table('./7.annotation/vep/vep_merge.maf',header = T,sep = '\t',quote = "")
maf[1:5,1:5]

## 构建肿瘤突变数据框，需要5列信息: sample.ID,chr,pos,ref,alt 
sample.mut.ref <- data.frame(Sample=maf[,16], 
							chr = maf[,5],
							pos = maf[,6],
							ref = maf[,11],
							alt = maf[,13])

sample.mut.ref[1:5,1:5]
				  Sample  chr      pos ref alt
1 case1_biorep_A_techrep chr1  6146376   G   T
2 case1_biorep_A_techrep chr1  6461445   G   T
3 case1_biorep_A_techrep chr1 31756671   C   A
4 case1_biorep_A_techrep chr1 32672798   A   T
5 case1_biorep_A_techrep chr1 39441098   G   T

sigs.input <- mut.to.sigs.input(mut.ref = sample.mut.ref,
								sample.id = "Sample", 
								chr = "chr", 
								pos = "pos", 
								ref = "ref", 
								alt = "alt",
								bsg = BSgenome.Hsapiens.UCSC.hg38)

> sigs.input[1:5,1:5]
					   A[C>A]A A[C>A]C A[C>A]G A[C>A]T C[C>A]A
case1_biorep_A_techrep       8       6      22       7       9
case1_biorep_B              10       4      12       8       7
case1_biorep_C              14       4      17       4       7
case1_techrep_2              7       5      10       1       7
case2_biorep_A              15       5       5       6       6


sigs.output <- whichSignatures(tumor.ref = sigs.input,
								signatures.ref = signatures.cosmic, 
								sample.id = 'case1_biorep_A_techrep',
								contexts.needed = TRUE)


# Plot output
plotSignatures(sigs.output)

# pheatmap
df = data.frame()
for (i in rownames(sigs.input)) {
  sigs.output <- whichSignatures(tumor.ref = sigs.input,
                                 signatures.ref = signatures.cosmic, 
                                 sample.id = i,
                                 contexts.needed = TRUE)
  df = rbind(df,sigs.output$weights)
}
df = df[ , apply(df, 2, function(x){sum(x>0)})>5]
pheatmap::pheatmap(df,cluster_cols = F,cluster_rows = F,fontsize = 20)

#14 初探肿瘤异质性
rm(list = ls())
options(stringsAsFactors = F)
library(dplyr)
library(stringr)
# 读入数据
laml = read.maf('./7.annotation/vep/vep_merge.maf')
laml@data = laml@data[!grepl('^MT-', laml@data$Hugo_Symbol),] 
# 增加一列t_vaf，即肿瘤样本中突变位点的覆盖深度t_alt_count占测序覆盖深度t_depth的比值
laml@data$t_vaf = (laml@data$t_alt_count/laml@data$t_depth)
unique(laml@data$Tumor_Sample_Barcode)
getSampleSummary(laml) 
getGeneSummary(laml) 
getFields(laml)


mut = laml@data[laml@data$t_alt_count >= 5 &
				  laml@data$t_vaf >= 0.05, c("Hugo_Symbol",
											 "Chromosome",
											 "Start_Position",
											 "Tumor_Sample_Barcode",
											 "t_vaf")]

mut$patient = substr(mut$Tumor_Sample_Barcode, 1, 5)

pid = unique(mut$patient)

lapply(pid , function(p){
	# p = 'case1'
	print(p)
	mat = unique(mut[mut$patient == p,c("Tumor_Sample_Barcode",'Hugo_Symbol')]) 
	mat$tmp = 1
	# 长变扁
	mat = spread(mat,Tumor_Sample_Barcode,tmp,fill = 0)
	class(mat)
	mat = as.data.frame(mat)
	rownames(mat) = mat$Hugo_Symbol
	mat = mat[,-1]
	dat = mat[order(mat[,1],mat[,2],mat[,3],mat[,4]),]
	png(paste0('overlap_', p, '.png'),width = 300, height = 1500)
	pheatmap::pheatmap(mat = dat, cluster_cols = F, cluster_rows = F, show_rownames = T, legend = F,fontsize_row = 16,fontsize_col = 30)
	dev.off()
})

rm(list = ls())
options(stringsAsFactors = F)
library(dplyr)
library(stringr)
# 读入数据
laml = read.maf('./7.annotation/vep/vep_merge.maf')
laml@data = laml@data[!grepl('^MT-',laml@data$Hugo_Symbol),] 
# 增加一列 t_vaf，即肿瘤样本中突变位点的覆盖深度t_alt_count占测序覆盖深度t_depth的比值
laml@data$t_vaf = (laml@data$t_alt_count/laml@data$t_depth)
heter = inferHeterogeneity(maf = laml, top = 24,vafCol = 't_vaf')

plotClusters(clusters = heter,
			 showCNvars = T,
			 tsb = 'case1_biorep_A_techrep')

for (i in unique(laml@data$Tumor_Sample_Barcode)) {
  # i = 'case1_biorep_A_techrep'
  png(paste0('vaf_', i, '.png'),width = 500, height = 330)
  plotClusters(clusters = heter,
			   tsb = i,
			   showCNvars = T)
  dev.off()
}


##15 对突变位点所在的基因进行 KEGG 注释
rm(list=ls())
options(stringsAsFactors = F)
library(dplyr)
library(stringr)
# 读入数据
laml = read.maf('./7.annotation/vep/vep_merge.maf')
laml@data=laml@data[!grepl('^MT-',laml@data$Hugo_Symbol),] 
# 增加一列t_vaf，即肿瘤样本中突变位点的覆盖深度t_alt_count占测序覆盖深度t_depth的比值
laml@data$t_vaf = (laml@data$t_alt_count/laml@data$t_depth)
unique(laml@data$Tumor_Sample_Barcode)
getSampleSummary(laml) 
getGeneSummary(laml) 
getFields(laml)

mut = laml@data[laml@data$t_alt_count >= 5 &
				  laml@data$t_vaf >= 0.05, c("Hugo_Symbol",
											 "Chromosome",
											 "Start_Position",
											 "Tumor_Sample_Barcode",
											 "t_vaf")]

mut$patient = substr(mut$Tumor_Sample_Barcode, 1, 5)
pid = unique(mut$patient)

all_snv = lapply(pid , function(p){
	# p='case1'
	print(p)
	mat=unique(mut[mut$patient %in% p,c("Tumor_Sample_Barcode",'Hugo_Symbol')]) 
	mat$tmp = 1
	# 长变扁
	mat = spread(mat,Tumor_Sample_Barcode,tmp,fill = 0)
	class(mat)
	mat = as.data.frame(mat)
	rownames(mat) = mat$Hugo_Symbol
	mat=mat[,-1]
	dat = mat[order(mat[,1],mat[,2],mat[,3],mat[,4]),]
	return(dat)
})

trunk_gene = unlist(sapply(all_snv, function(x) rownames(x[rowSums(x) == 4,])))
branch_gene = unlist(sapply(all_snv, function(x) rownames(x[rowSums(x) == 3|2,])))
private_gene = unlist(sapply(all_snv, function(x) rownames(x[rowSums(x) == 1,])))

BiocManager::install("org.Hs.eg.db")
library(org.Hs.eg.db)
library(clusterProfiler)
kegg_SYMBOL_hsa <- function(genes){ 
  gene.df <- bitr(genes, fromType = "SYMBOL",
				  toType = c("SYMBOL", "ENTREZID"),
				  OrgDb = org.Hs.eg.db)
  head(gene.df) 
  diff.kk <- enrichKEGG(gene         = gene.df$ENTREZID,
						organism     = 'hsa',
						pvalueCutoff = 0.99,
						qvalueCutoff = 0.99
  )
  return( setReadable(diff.kk, OrgDb = org.Hs.eg.db,keyType = 'ENTREZID'))
}

trunk_kk=kegg_SYMBOL_hsa(trunk_gene)
trunk_df=trunk_kk@result
write.csv(trunk_df,file = 'trunk_kegg.csv')
png(paste0('trunk_kegg', '.png'),width = 1080,height = 540)
barplot(trunk_kk,font.size = 20)
dev.off()

branch_kk=kegg_SYMBOL_hsa(branch_gene)
branch_df=branch_kk@result
write.csv(branch_df,file = 'branch_kegg.csv')
png(paste0('branch_kegg', '.png'),width = 1080,height = 540)
barplot(branch_kk,font.size = 20)
dev.off()

private_kk=kegg_SYMBOL_hsa(private_gene)
private_df=private_kk@result
write.csv(private_df,file = 'private_kegg.csv')
png(paste0('private_kegg', '.png'),width = 1080,height = 540)
barplot(private_kk,font.size = 20)
dev.off()



#pyclone
rm(list = ls())
options(stringsAsFactors = F)
case1_loci = read.table("./9.pyclone/case1_pyclone_analysis/tables/loci.tsv",
						header = T)
# 获取clusters 的分组信息
clusters_list = unique(case1_loci[, c(1, 3)])
rownames(clusters_list) = clusters_list[, 1]
cluster_id = data.frame(cluster_id = as.character(clusters_list$cluster_id))
rownames(cluster_id) = clusters_list[, 1]
# 获取同一个突变位点在不同样本中的cellular_prevalence，然后画热图可视化
library(tidyr)
cellular_prevalence = spread(case1_loci[, c(1, 2, 4)], key = sample_id, value = cellular_prevalence)
rownames(cellular_prevalence) = cellular_prevalence[, 1]
cellular_prevalence = cellular_prevalence[,-1]
sampe_id = colnames(cellular_prevalence)
cellular_prevalence = as.data.frame(t(apply(cellular_prevalence, 1, as.numeric)))
colnames(cellular_prevalence) = sampe_id
pheatmap::pheatmap(
  cellular_prevalence,
  annotation_row = cluster_id,
  show_rownames = F,
  clustering_method = 'median',
  angle_col = 0
)
# 获取同一个突变位点在不同样本中的variant_allele_frequency，也就是vaf，同样可视化，为了聚类，采用了不同的聚类方法

library(tidyr)
variant_allele_frequency = spread(case1_loci[, c(1, 2, 6)], key = sample_id, value = variant_allele_frequency)
rownames(variant_allele_frequency) = variant_allele_frequency[, 1]
variant_allele_frequency = variant_allele_frequency[,-1]
sampe_id = colnames(variant_allele_frequency)
variant_allele_frequency = as.data.frame(t(apply(variant_allele_frequency, 1, as.numeric)))
colnames(variant_allele_frequency) = sampe_id

pheatmap::pheatmap(
  cellular_prevalence,
  annotation_row = cluster_id,
  show_rownames = F,
  clustering_method = 'average',
  angle_col = 0
)


#Timescape 可视化

biocmanager::install("timescape")
library(timescape)
options(stringsAsFactors = F)
options(browser="firefox")
example("timescape")
browseVignettes("timescape") 
library(plotly)
library(htmlwidgets)
library(webshot)
for (i in 1:6) {
  # i = 1
  # tree 
  tree_edges = read.table(paste0("9.pyclone/case", i, "_pyclone_analysis/tree.txt"))
  colnames(tree_edges) = c("source","target")
  # clonal prevalences
  cellfreq = read.table(paste0("9.pyclone/case", i, "_pyclone_analysis/cellfreq.txt"))
  colnames(cellfreq) = 0:(length(cellfreq)-1)
  sample_id = read.table(paste0("9.pyclone/case", i, "_pyclone_analysis/sample_id"))
  cellfreq$timepoint = sample_id[ , 1]
  library(tidyr)
  clonal_prev = gather(cellfreq, key="clone_id", value = "clonal_prev", -timepoint)
  clonal_prev = clonal_prev[order(clonal_prev$timepoint),]
  # targeted mutations
  # mutations <- read.csv(system.file("extdata", "AML_mutations.csv", package = "timescape"))
  p = timescape(clonal_prev = clonal_prev, tree_edges = tree_edges,height=260)
  saveWidget(p, paste0("case", i,"_timescape", ".html"))
  }


