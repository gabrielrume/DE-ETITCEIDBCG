

#MAP reads - using markdup

# fastq_15182Byr_N403_L003

for i in `ls *R1_001.fastq.gz`;
do
~/bin/STAR-2.5.3a/bin/Linux_x86_64/STAR --runThreadN 8 --genomeDir /home/thales/Bioinfo/genomes/C_c/ --readFilesIn $i `basename -s R1_001.fastq.gz $i`R2_001.fastq.gz --outFileNamePrefix `basename -s R1_001.fastq.gz $i` --readFilesCommand zcat --outSAMstrandField intronMotif --outSAMtype BAM Unsorted;

done

#Mark dup with picard


#The program can take either coordinate-sorted or query-sorted inputs, however the behavior is slightly different. When the input is coordinate-sorted, unmapped mates of mapped records and supplementary/secondary alignments are not marked as duplicates. However, when the input is query-sorted (actually query-grouped), then unmapped mates and secondary/supplementary reads are not excluded from the duplication test and can be marked as duplicate reads.


for i in `ls *Aligned.out.bam`;
do
java -XX:ParallelGCThreads=8 -jar ~/bin/picard.jar SortSam I=$i O=`basename -s .bam $i`.sorted.pikard.query.bam SORT_ORDER=queryname &
done 


for i in `ls *sorted.pikard.query.bam`;
do
java -XX:ParallelGCThreads=8 -jar ~/bin/picard.jar MarkDuplicates I=$i O=`basename -s .bam $i`.markdup_via_pikard.bam REMOVE_DUPLICATES=true M=`basename -s .bam $i`.markdup_via_pikard.txt
done


for i in `ls *markdup_via_pikard.bam`
do
htseq-count -a 10 -t exon -i Parent -f bam --stranded=no $i /home/thales/Bioinfo/genomes/C_c/coffea_canephora.gff3   > `basename -s .bam $i`_rmdup_counts_picard.txt &
#htseq-count -a 10 -t exon -i Parent -f bam --stranded=no $i /home/thales/Bioinfo/genomes/C_c/coffea_canephora.gff3   > `basename -s .bam $i`_rmdup_counts_picard_no_tranded.txt &
#htseq-count -a 10 -t exon -i Parent -f bam --stranded=yes $i /home/thales/Bioinfo/genomes/C_c/coffea_canephora.gff3   > `basename -s .bam $i`_rmdup_counts_picard_yes.txt &
#htseq-count -a 10 -t exon -i Parent -f bam --stranded=reverse $i /home/thales/Bioinfo/genomes/C_c/coffea_canephora.gff3   > `basename -s .bam $i`_rmdup_counts_picard_reverse.txt &
done

######################################################
################ edgeR rmdup #########################
######################################################


#create a file Targets.txt:
files	group	description
./counts/C1P1-2_S1_L003_Aligned.out.sorted.pikard.query.markdup_via_pikard_rmdup_counts_picard.txt	CC	CC
./counts/C1P2-1_S2_L003_Aligned.out.sorted.pikard.query.markdup_via_pikard_rmdup_counts_picard.txt	CC	CC
./counts/C1P3-2_S3_L003_Aligned.out.sorted.pikard.query.markdup_via_pikard_rmdup_counts_picard.txt	CC	CC
./counts/C1P4-2_S4_L003_Aligned.out.sorted.pikard.query.markdup_via_pikard_rmdup_counts_picard.txt	CC	CC
./counts/C1P5-2_S5_L003_Aligned.out.sorted.pikard.query.markdup_via_pikard_rmdup_counts_picard.txt	CC	CC
./counts/C1P6-1_S6_L003_Aligned.out.sorted.pikard.query.markdup_via_pikard_rmdup_counts_picard.txt	SC	SC
./counts/C1P7-2_S7_L003_Aligned.out.sorted.pikard.query.markdup_via_pikard_rmdup_counts_picard.txt	SC	SC
./counts/C1P8-2_S8_L003_Aligned.out.sorted.pikard.query.markdup_via_pikard_rmdup_counts_picard.txt	SC	SC
./counts/C1P9-1_S9_L003_Aligned.out.sorted.pikard.query.markdup_via_pikard_rmdup_counts_picard.txt	SC	SC
./counts/C1P10-2_S10_L003_Aligned.out.sorted.pikard.query.markdup_via_pikard_rmdup_counts_picard.txt	SC	SC
#./counts/C2P1-2_S11_L003_Aligned.out.sorted.pikard.query.markdup_via_pikard_rmdup_counts_picard.txt	CA	CA
./counts/C2P2-2_S12_L003_Aligned.out.sorted.pikard.query.markdup_via_pikard_rmdup_counts_picard.txt	CA	CA
./counts/C2P3-2_S13_L003_Aligned.out.sorted.pikard.query.markdup_via_pikard_rmdup_counts_picard.txt	CA	CA
./counts/C2P4-1_S14_L003_Aligned.out.sorted.pikard.query.markdup_via_pikard_rmdup_counts_picard.txt	CA	CA
./counts/C2P5-1_S15_L003_Aligned.out.sorted.pikard.query.markdup_via_pikard_rmdup_counts_picard.txt	CA	CA
./counts/C2P6-2_S16_L003_Aligned.out.sorted.pikard.query.markdup_via_pikard_rmdup_counts_picard.txt	SA	SA
./counts/C2P7-1_S17_L003_Aligned.out.sorted.pikard.query.markdup_via_pikard_rmdup_counts_picard.txt	SA	SA
./counts/C2P8-1_S18_L003_Aligned.out.sorted.pikard.query.markdup_via_pikard_rmdup_counts_picard.txt	SA	SA
./counts/C2P9-2_S19_L003_Aligned.out.sorted.pikard.query.markdup_via_pikard_rmdup_counts_picard.txt	SA	SA
./counts/C2P10-1_S20_L003_Aligned.out.sorted.pikard.query.markdup_via_pikard_rmdup_counts_picard.txt	SA	SA


# reads the tables of counts, calculates the sizes of the count libraries and produces a DGEList object for use by subsequent functions.
#GLM_edgeR approach - Generalised Linear Model
 library("edgeR")
#get the path to the counts files
targets <- readTargets()
names <- targets$description

#create a DGE matrix
matrix_input <- readDGE(targets, comment.char = "!")
#remove meta Tags
MetaTags <- grep("^__", rownames(matrix_input))
matrix_input <- matrix_input[-MetaTags, ]

reads_before <- sum(matrix_input$counts)

#remove low expressed genes 
rnaseqmatrix <- matrix_input$counts
rnaseqmatrix <- rnaseqmatrix[rowMeans(rnaseqmatrix) >=10,]
conditions = matrix_input$samples[,2]
analysis_matrix <- DGEList(counts = rnaseqmatrix,group = conditions)
colnames(analysis_matrix$counts) <- names
#create dea design (or model) matrix.
#Here, the 0+ in the model formula is an instruction not to include an intercept column and instead to include a column for each group. Beter in situation of "pair-wise" comparisons

design <- model.matrix(~0+group, data=analysis_matrix$samples)
colnames(design) <- levels(analysis_matrix$samples$group)

#NORMALIZATIONS
analysis_matrix <- calcNormFactors(analysis_matrix)

#To estimate common dispersion: 
analysis_matrix <- estimateGLMCommonDisp(analysis_matrix, design)
#To estimate trended dispersions:
analysis_matrix <- estimateGLMTrendedDisp(analysis_matrix, design)
#To estimate tagwise dispersions:
analysis_matrix <- estimateGLMTagwiseDisp(analysis_matrix, design)

#Fit a negative binomial generalized log-linear model to the read counts for each gene.
fit <- glmFit(analysis_matrix,design)

############edgeR Plots##################
#take care with very high number of genes >20k
#norm.vals <- cpm(analysis_matrix)
#library(heatmap3)
#pdf(file = "edgeR_normalization_heatmap_rmdup_picard.pdf")
#heatmap3(norm.vals)
#dev.off()

#visualize the libraries sizes
pdf("libraries_sizes.pdf", wi=12,he=8)
barplot(fit$samples$lib.size, names.arg=names,col=c(rep("firebrick2",5),rep("orange1",5),rep("chartreuse",4),rep("mediumblue",5)))
abline(h=mean(fit$samples$lib.size)-2*sd(fit$samples$lib.size),col="black",lwd = 3)
abline(h=mean(fit$samples$lib.size),lty = 2,lwd = 3)
dev.off()
#The square root of dispersion is the coe cient of biological variation (BCV)

pdf(file = "edgeR_BCV.pdf")
plotBCV(analysis_matrix)
dev.off()

#An MDS plots shows distances, in terms of biological coeficient of variation (BCV) - An MDS plot shows the relative similarities of the samples.
pdf(file = "edgeR_MDS.pdf")
plotMDS(analysis_matrix)
dev.off()

#samples_trees
cpm.matrix <- cpm(analysis_matrix,normalized.lib.sizes=TRUE)
colnames(cpm.matrix) <- names
t.cpm.matrix <- t(cpm.matrix)
sampleTree <- hclust(dist(t.cpm.matrix), method = "average");


# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
# The user should change the dimensions if the window is too large or too small.
#sizeGrWindow(12,9)
pdf(file = "sampleClustering.pdf", width = 12, height = 9);
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
cex.axis = 1.5, cex.main = 2)
# Plot a line to show the cut
#abline(h = 15, col = "red")
dev.off()


#analyse densisties before normalization - input
pdf("Densities_input.pdf")
plotDensities( log(matrix_input$counts), legend = "topright")
dev.off()

#analyse densisties low expression filter
pdf("Densities_low_expression_fiter.pdf")
plotDensities( log(rnaseqmatrix), legend = "topright")
dev.off()

#analyse densisties after normalization
pdf("Densities_normalization.pdf")
plotDensities( log(cpm.matrix), legend = "topright")
dev.off()

pdf("Densities_log_cpm_fitted_norm.pdf")
plotDensities(log(cpm(fit$fitted.values, normalized.lib.sizes=TRUE)), legend = "topright")
dev.off()

#histogram of densities log10
pdf("Log10_histogram_normilized.pdf", h=10,w=10)
hist(log(cpm.matrix+1,10), col=gray.colors(19, start = 0.9, end = 0.3))
dev.off()

#histogram of densities Log2
pdf("Log2_histogram_normilized.pdf", h=10,w=10)
hist(log(cpm.matrix+1,2), col=gray.colors(19, start = 0.9, end = 0.3))
dev.off()

#histogram of densities no log
pdf("histogram_normilized.pdf", h=10,w=10)
hist(cpm.matrix, col=gray.colors(19, start = 0.9, end = 0.3))
dev.off()

# Load library for pheatmap
cpm.matrix.corr <- cor(cpm.matrix, method="spearman",use="pairwise.complete.obs")

library("pheatmap")
library(gplots)
pdf(file = "sampleClusteringHeatmap.pdf", width = 25, height = 20);
pheatmap(cpm.matrix.corr, fontsize=20,cellwidth=70, cellheight=60, treeheight_col= 200, treeheight_row=200,angle_col=0, legend=T)
#heatmap.2(cpm.matrix.corr*10000,trace = "none",margins = c(5, 11),keysize =1 , key.title="",col ="bluered",density.info="none")
dev.off()

#print(paste0("There were ",reads_before ," reads aligned in all dataset"))
#print(paste0("There are ", sum(analysis_matrix$counts)," after the removal of low expressed elements (", (sum(analysis_matrix$counts)/reads_before)*100, "%)"))
#print(paste("There are ", nrow(as.matrix(dedpuplicated_keeped_names)), " non-duplicatede miRNAs DE"))


#gene expression profile
normalized.cpm.matrix <- cpm(analysis_matrix,normalized.lib.sizes=T)
colnames(normalized.cpm.matrix) <- names

system("mkdir expression_profiles")
for (gene in rownames(normalized.cpm.matrix)){
pdf(paste0("./expression_profiles/expression_profile_",gene,".pdf"),he=8,wi=12)
sdgene <- sd(normalized.cpm.matrix[gene,])

plot(type="b", x=c(1:19), y=normalized.cpm.matrix[gene,],ylab="",xlab = "", ylim = c(0-(sdgene),max(normalized.cpm.matrix[gene,])+sdgene))
axis(side=1,at=seq(0,19,1))

arrows(x0=c(1:19),y0=normalized.cpm.matrix[gene,], y1= normalized.cpm.matrix[gene,] + sdgene, lwd=1.3,angle =90, length=.05, code =2 )
arrows(x0=c(1:19),y0=normalized.cpm.matrix[gene,], y1= normalized.cpm.matrix[gene,] - sdgene, lwd=1.3,angle =90, length=.05, code =2 )

abline(h=0,col="red",lwd = 2 )
abline(v= 5.5, col="blue", lty = 2, lwd=1.5)
abline(v= 10.5, col="blue", lty = 2, lwd=1.5)
abline(v= 14.5, col="blue", lty = 2, lwd=1.5)

mtext(cex=.95,side = 3,c("CC1      CC2      CC3      CC4      CC5      SC1      SC2      SC3      SC4      SC5      CA1      CA2      CA3      CA4      SA1      SA2      SA3      SA4      SA5"))
dev.off()
}


save.image("edgeR_rmdup_via_picard.RData")

############edgeR DE######################


#Control Acauã (CA) x Stress Acauã (SA)

lrt_ca_x_sa <- glmLRT(fit, contrast = c(-1,0,1,0))
tTags_lrt_ca_x_sa <- topTags(lrt_ca_x_sa, n= NULL)
keep_sig <- matrix(tTags_lrt_ca_x_sa$table$FDR <= 0.05 & abs(tTags_lrt_ca_x_sa$table$logFC) >= 1)
sig_kept_ca_x_sa <- tTags_lrt_ca_x_sa$table[keep_sig[,1],]
write.table(sig_kept_ca_x_sa, file = "ca_x_sa.tab", row.names = T, quote = F, sep = "\t")

#Control Catuaí (CA) x Control acauã (CA)

lrt_cc_x_ca <- glmLRT(fit, contrast = c(1,-1,0,0))
tTags_lrt_cc_x_ca <- topTags(lrt_cc_x_ca, n= NULL)
keep_sig <- matrix(tTags_lrt_cc_x_ca$table$FDR <= 0.05 & abs(tTags_lrt_cc_x_ca$table$logFC) >= 1)
sig_kept_cc_x_ca <- tTags_lrt_cc_x_ca$table[keep_sig[,1],]
write.table(sig_kept_cc_x_ca, file = "cc_x_ca.tab", row.names = T, quote = F, sep = "\t")

#Control Catuai (CC) x Stress catuai (SC)

lrt_cc_x_sc <- glmLRT(fit, contrast = c(0,-1,0,1))
tTags_lrt_cc_x_sc <- topTags(lrt_cc_x_sc, n= NULL)
keep_sig <- matrix(tTags_lrt_cc_x_sc$table$FDR <= 0.05 & abs(tTags_lrt_cc_x_sc$table$logFC) >= 1)
sig_kept_cc_x_sc <- tTags_lrt_cc_x_sc$table[keep_sig[,1],]
write.table(sig_kept_cc_x_sc, file = "cc_x_sc.tab", row.names = T, quote = F, sep = "\t")

#Stress Catuai (SC) x Stress acaua (SA)

lrt_sc_x_sa <- glmLRT(fit, contrast = c(0,0,1,-1))
tTags_lrt_sc_x_sa <- topTags(lrt_sc_x_sa, n= NULL)
keep_sig <- matrix(tTags_lrt_sc_x_sa$table$FDR <= 0.05 & abs(tTags_lrt_sc_x_sa$table$logFC) >= 1)
sig_kept_sc_x_sa <- tTags_lrt_sc_x_sa$table[keep_sig[,1],]
write.table(sig_kept_sc_x_sa, file = "sc_x_sa.tab", row.names = T, quote = F, sep = "\t")

#find the dimensions of all significant DE matrix
dim(sig_kept_cc_x_ca); dim(sig_kept_ca_x_sa); dim(sig_kept_cc_x_sc); dim(sig_kept_sc_x_sa);

 
 #find all mRNAs DE expressed
dedpuplicated_keeped_names <-unique(c(rownames(sig_kept_cc_x_ca),rownames(sig_kept_ca_x_sa),rownames(sig_kept_cc_x_sc),rownames(sig_kept_sc_x_sa)))
write.table(as.matrix(dedpuplicated_keeped_names),row.names=F, file = "dedup_sequence_names.tab")


nrow(as.matrix(dedpuplicated_keeped_names))


#create a matrix with the values of expression of the DE mRNAs
FC_matrix <- matrix(c(lrt_cc_x_ca$table[intersect(dedpuplicated_keeped_names,rownames(lrt_cc_x_ca)),1],
lrt_ca_x_sa$table[intersect(dedpuplicated_keeped_names, rownames(lrt_ca_x_sa$table)),1],
lrt_cc_x_sc$table[intersect(dedpuplicated_keeped_names, rownames(lrt_cc_x_sc$table)),1],
lrt_sc_x_sa$table[intersect(dedpuplicated_keeped_names, rownames(lrt_sc_x_sa$table)),1]),ncol=4)

colnames(FC_matrix) <- c("cc_x_ca","ca_x_sa","cc_x_sc","sc_x_sa")
rownames(FC_matrix) <- dedpuplicated_keeped_names

library(gplots)

plotCol <- colorRampPalette(c("#000dff","white","#ff0004"))(16)


pdf("FC_heatmap_only_DE.pdf", title="created by THALES HENRIQUE CHERUBINO RIBEIRO (THC.RIBEIRO) - LFMP - UFLA",width =7.09, height=7.09)
heatmap.2(FC_matrix, trace="none", col =plotCol,srtCol=75, margins = c(6, 1),cexCol=1, key.title="",key.ylab="",  key.xlab="Log2 FC")
dev.off()

system("open FC_heatmap_only_DE.pdf")

#veen

library(VennDiagram)
pdf("Veen_Genotype_differences.pdf",width=7.09,height=7.09)
draw.pairwise.venn(area1 = 182, area2 = 25, cross.area = 21, fill = c("#7291ff","#ff7581"), scaled = TRUE, cex = rep (3,3), cat.cex=rep(2,2),ext.text=F)
dev.off()
system("open Veen_Genotype_differences.pdf")

library(VennDiagram)
pdf("Veen_Stress_differences.pdf",width=7.09,height=7.09)
draw.pairwise.venn(area2 = 22, area1 = 36, cross.area = 6, fill = c("#ffb472","#80ff72"), scaled = TRUE, cex = rep (3,3), cat.cex=rep(2,2),ext.text=F,inverted=T)
dev.off()
system("open Veen_Stress_differences.pdf")


#cat.pos=c(1,1)
#category=c("Down_G5_Cat","Down_G5_Sir")

#find number of expressed DE mRNAs

expression_barplot_DE <- as.matrix(c(nrow(sig_kept_cc_x_ca),nrow(sig_kept_sc_x_sa), nrow(sig_kept_ca_x_sa),nrow(sig_kept_cc_x_sc)))


rownames(expression_barplot_DE) <- c("cc_x_ca","sc_x_sa","ca_x_sa","cc_x_sc")

#plot number of expressed DE miRNAs

pdf("DE_all_conditions_barplot.pdf", width = 9, height = 7,, title="created by THALES HENRIQUE CHERUBINO RIBEIRO (THC.RIBEIRO) - LFMP - UFLA")
barplot(expression_barplot_DE, beside = T, space = c(.1,1), main="number of DE mRNAs per comparison", col = c("firebrick2", "orange1","chartreuse","mediumblue"),legend.text=rownames(expression_barplot_DE), ylim=c(0,500))
abline(h=0)
dev.off()

#barplot with number of UP or Down expresed mRNAs


pdf(file="number_of_DE_mRNAs_all_comparisons.pdf", height = 6,width= 7.09, title="created by THALES HENRIQUE CHERUBINO RIBEIRO (THC.RIBEIRO) - LFMP - UFLA")

barplot(nrow(sig_kept_cc_x_ca[sig_kept_cc_x_ca$logFC>=1,]), col = "#3344ff",
beside = F, ylim=c(-200,100), space = 1,xlim = c(0,8),
main ="Number of DE mRNAs up and down", ylab="",xaxt='n',yaxt='n')
barplot(-nrow(sig_kept_cc_x_ca[sig_kept_cc_x_ca$logFC<=1,]), col = "#3344ff",
beside = F, ylim=c(-200,400), space = 1,
add=T, ylab="",xaxt='n',yaxt='n')

barplot(nrow(sig_kept_sc_x_sa[sig_kept_sc_x_sa$logFC>=1,]), col = "#ff333a",
beside = F, ylim=c(-200,400), space = 3,
add=T, ylab="",xaxt='n',yaxt='n')
barplot(-nrow(sig_kept_sc_x_sa[sig_kept_sc_x_sa$logFC<=1,]), col = "#ff333a",
beside = F, ylim=c(-200,400), space = 3,
 add=T, ylab="",xaxt='n',yaxt='n')

barplot(nrow(sig_kept_ca_x_sa[sig_kept_ca_x_sa$logFC>=1,]), col = "#ff9633",
beside = F, ylim=c(-200,400), space = 7,
add=T, ylab="",xaxt='n',yaxt='n')
barplot(-nrow(sig_kept_ca_x_sa[sig_kept_ca_x_sa$logFC<=1,]), col = "#ff9633",
beside = F, ylim=c(-200,400), space = 7,
 add=T, ylab="",xaxt='n',yaxt='n')
 
 barplot(nrow(sig_kept_cc_x_sc[sig_kept_cc_x_sc$logFC>=1,]), col = "#3aff33",
beside = F, ylim=c(-200,400), space = 5,
add=T, ylab="",xaxt='n',yaxt='n')
barplot(-nrow(sig_kept_cc_x_sc[sig_kept_cc_x_sc$logFC<=1,]), col = "#3aff33",
beside = F, ylim=c(-200,400), space = 5,
 add=T, ylab="",xaxt='n',yaxt='n')
 
text(seq(1.5,7.5,by=2),y= -180,srt = 60, adj= 1, xpd = TRUE,labels = c("cc_x_ca","sc_x_sa","cc_x_sc","ca_x_sa"),cex=1.5)
    
 abline(h=0)

axis(side=2, at=(c(-200,-150,-100,-50,0,50,100)),labels=F)
mtext(side=2,at=c(-200,-150,-100,-50,0,50,100),text=c(200,150,100,50,0,50,100),line = 1, las=1,cex = 1)

    
 dev.off()

system("open number_of_DE_mRNAs_all_comparisons.pdf")

#per gene expression graph

#create a matrix with expression values in counts per millon

#cpm.matrix <- cpm(analysis_matrix,normalized.lib.sizes=TRUE)
#load the namesEquivalences produced from the gff file
namesEquivalences <- read.table("/Users/johnjoyce/Dropbox/Bioinformatica/LFMP_coffe/coffee_c_new_annot/names_equivalence.tab" ,sep = "\t", header = F)
rownames(namesEquivalences) <- namesEquivalences$V2

system("mkdir absolute_expression_OpT_x_WaT")
for (gene in rownames(normalized.cpm.matrix)){
    pdf(paste0("./absolute_expression_OpT_x_WaT/absoluteExpression_",namesEquivalences[gene,4],".pdf"),he=7.09,wi=7.09)
    #sdgene <- sd(normalized.cpm.matrix[gene,])

    #Catuaí
    mean_stress_CA <- mean(normalized.cpm.matrix[gene,c(6,7,8,9,10)])
    mean_control_CA <- mean(normalized.cpm.matrix[gene,c(1,2,3,4,5)])
    se_stress_CA <- sd(normalized.cpm.matrix[gene,c(6,7,8,9,10)])/sqrt(5)
    se_control_CA  <- sd(normalized.cpm.matrix[gene,c(1,2,3,4,5)])/sqrt(5)

    #Acauã
    mean_stress_AC <- mean(normalized.cpm.matrix[gene,c(15,16,17,18,19)])
    mean_control_AC <- mean(normalized.cpm.matrix[gene,c(11,12,13,14)])
    se_stress_AC <- sd(normalized.cpm.matrix[gene,c(15,16,17,18,19)])/sqrt(5)
    se_control_AC  <- sd(normalized.cpm.matrix[gene,c(11,12,13,14)])/sqrt(4)


    barplot <- barplot(c(mean_stress_AC, mean_control_AC, mean_stress_CA, mean_control_CA),beside=T,ylim=c(0,max(c(mean_stress_AC + se_stress_AC, mean_control_AC + se_control_AC,mean_stress_CA + se_stress_CA, mean_control_CA + se_control_CA )*1.20)),
    ylab="Expression values in CPM",names.arg = c("Acauã", "Acauã", "Catuaí", "Catuaí"),legend.text=c("WaT","OpT"),
    args.legend=list(x = "topright", bty="n", cex =2),main= paste0("Expression of gene ",  namesEquivalences[gene,4]),col=c("#ff7581","#7291ff"), cex.names = 2)

#col=c("grey26","grey66")

    points(barplot,c(mean_stress_AC, mean_control_AC, mean_stress_CA, mean_control_CA),pch=19)

    arrows(x0 = barplot, y0 = c(mean_stress_AC-se_stress_AC, mean_control_AC-se_control_AC,mean_stress_CA -se_stress_CA,mean_control_CA - se_control_CA),x1 = barplot, y1=c(mean_stress_AC+se_stress_AC, mean_control_AC+se_control_AC,mean_stress_CA +se_stress_CA,mean_control_CA + se_control_CA),code=3,angle=90,length=0.05)
#mtext(x = 1, y = -2, "Catuaí")
#mtext(x = 3, y = -2, "Acauã")
    abline(h=0)
    dev.off()
}
#log2FC of a specifc gene
#log(mean(cpm.matrix["GSCOCT00006570001",c(6,7,8,9,10)])/mean(cpm.matrix["GSCOCT00006570001",c(1,2,3,4,5)]),2)
save.image("new.edgeR.RData")

#################################################
################ Grep anotations ################
#################################################

#cc_x_sc
for transcript in `awk '{print $1}' RNA_IDs_cc_x_sc.tab` ;
    do
    echo grepping transcript $transcript;
    grep $transcript /Users/johnjoyce/Dropbox/Bioinformatica/LFMP_coffe/coffee_c_new_annot/names_equivalence.tab >> all_ids_cc_x_sc.tab;
done
for gene in `awk '{print $4}' all_ids_cc_x_sc.tab`;
    do
    echo grepping gene $gene;
    grep $gene /Users/johnjoyce/Dropbox/Bioinformatica/LFMP_coffe/coffee_c_new_annot/nr_plants/coffe_prots_nr.annot >> cc_x_sc_all.annot.tab
done

#sc_x_sa
for transcript in `awk '{print $1}' RNA_IDs_sc_x_sa.tab` ;
    do
    echo grepping transcript $transcript;
    grep $transcript /Users/johnjoyce/Dropbox/Bioinformatica/LFMP_coffe/coffee_c_new_annot/names_equivalence.tab >> all_ids_sc_x_sa.tab;
done

for gene in `awk '{print $4}' all_ids_sc_x_sa.tab`;
    do
    echo grepping gene $gene;
    grep $gene /Users/johnjoyce/Dropbox/Bioinformatica/LFMP_coffe/coffee_c_new_annot/nr_plants/coffe_prots_nr.annot >> sc_x_sa_all.annot.tab
done


#cc_x_ca
for transcript in `awk '{print $1}' RNA_IDs_cc_x_ca.tab` ;
    do
    echo grepping transcript $transcript;
    grep $transcript /Users/johnjoyce/Dropbox/Bioinformatica/LFMP_coffe/coffee_c_new_annot/names_equivalence.tab >> all_ids_cc_x_ca.tab;
done
for gene in `awk '{print $4}' all_ids_cc_x_ca.tab`;
    do
    echo grepping gene $gene;
    grep $gene /Users/johnjoyce/Dropbox/Bioinformatica/LFMP_coffe/coffee_c_new_annot/nr_plants/coffe_prots_nr.annot >> cc_x_ca_all.annot.tab
done

#ca_x_sa

for transcript in `awk '{print $1}' RNA_IDs_ca_x_sa.tab` ;
do
echo grepping transcript $transcript;
grep $transcript /Users/johnjoyce/Dropbox/Bioinformatica/LFMP_coffe/coffee_c_new_annot/names_equivalence.tab >> all_ids_ca_x_sa.tab;
done
for gene in `awk '{print $4}' all_ids_ca_x_sa.tab`;
do
echo grepping gene $gene;
grep $gene /Users/johnjoyce/Dropbox/Bioinformatica/LFMP_coffe/coffee_c_new_annot/nr_plants/coffe_prots_nr.annot >> ca_x_sa_all.annot.tab
done


#get GO terms

#sc_x_sa_all.annot.tab

for gene in `awk '{print $1}' sc_x_sa_all.annot.tab`
do
    echo $gene
    grep $gene  ~/Dropbox/Bioinformatica/LFMP_coffe/coffee_c_new_annot/nr_plants/custon_ontology_background.tab >> GO_sc_x_sa_all.tab
done

#cc_x_ca_all.annot.tab

for gene in `awk '{print $1}' cc_x_ca_all.annot.tab`
do
echo $gene
grep $gene  ~/Dropbox/Bioinformatica/LFMP_coffe/coffee_c_new_annot/nr_plants/custon_ontology_background.tab >> GO_cc_x_ca_all.tab
done

#ca_x_sa_all.annot.tab

for gene in `awk '{print $1}' ca_x_sa_all.annot.tab`
do
echo $gene
grep $gene  ~/Dropbox/Bioinformatica/LFMP_coffe/coffee_c_new_annot/nr_plants/custon_ontology_background.tab >> GO_ca_x_sa_all.tab
done

#cc_x_sc_all.annot.tab

for gene in `awk '{print $1}' cc_x_sc_all.annot.tab`
do
echo $gene
grep $gene  ~/Dropbox/Bioinformatica/LFMP_coffe/coffee_c_new_annot/nr_plants/custon_ontology_background.tab >> GO_cc_x_sc_all.tab
done

#get sequences
#sc_x_sa
for sequence in `awk '{print $4}' all_ids_sc_x_sa.tab`
do
    grep -A 1 $sequence ~/genomes/C_c_genome/coffea_pep.single_lines.fna >> sc_x_sa_all.sequences.fasta
done

#cc_x_ca
for sequence in `awk '{print $4}' all_ids_cc_x_ca.tab`
do
    grep -A 1 $sequence ~/genomes/C_c_genome/coffea_pep.single_lines.fna >> cc_x_ca_all.sequences.fasta
done

#ca_x_sa
for sequence in `awk '{print $4}' all_ids_ca_x_sa.tab`
do
    grep -A 1 $sequence ~/genomes/C_c_genome/coffea_pep.single_lines.fna >> ca_x_sa_all.sequences.fasta
done

#cc_x_sc

for sequence in `awk '{print $4}' all_ids_cc_x_sc.tab`
do
    grep -A 1 $sequence ~/genomes/C_c_genome/coffea_pep.single_lines.fna >> cc_x_sc_all.sequences.fasta
done
