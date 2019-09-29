getwd()
setwd("D:/Introduction to R")
setwd("P:/Introduction-to-R-master/Introduction-to-R-master")



library(RDocumentation)
library(dint)
library(Biobase)
library(genefilter)
library(affy)
#library("rSFFreader")
library(oligo)
library(GEOquery)
library(normtest)
library(limma)
library(reshape2)
library(gplots)
library(ggplot2)
library(seqinr)
library(ape) 
library(stringr)
library(Peptides) 
library(Biostrings)
library(timeSeries)
library(MASS)
library(rgl)
library(CrossValidate)
library(random)
library(edgeR)
library(limma)
library(Glimma)
#library(org.Mm.eg.db)
library(RColorBrewer)
#library(Mus.musculus)
library(lavaan)
library(semPlot)
library(corrplot)
library(multcomp)
install.packages("ROCR")
library(ROCR)

#cell2label
cell2label<- function(cell){
label<- rep(0, length(cell))
lb <- unique(cell)
N<- length(lb)
#sapply(seq(1, length(cell)), label[i]<-{for (j in lb){if(cell[i] == j) return(j)}})
for (i in seq(1, length(cell))){for (j in seq(1,N)) {if(cell[i] ==lb[j]) label[i]<-j}}
}

#read in the DC MC
data<- read.csv('D:/TUT/Medical/biophysics/experiment/WiML/raw_data_DGinMCout.tsv',sep ='	', header= TRUE)
#data<- read.csv('P:/Introduction-to-R-master/Introduction-to-R-master/raw_data_DGinMCout.tsv',sep ='	', header= TRUE)
#gene <- cell2label(data[,2])
#cell <- cell2label(data[,3])
cell<- data[,2]
label<- rep(0, length(cell))
lb <- unique(cell)
N<- length(lb)
for (i in seq(1, length(cell))){for (j in seq(1,N)) {if(cell[i] ==lb[j]) label[i]<-j}}
gene<- label

cell<- data[,3]
label<- rep(0, length(cell))
lb <- unique(cell)
N<- length(lb)
for (i in seq(1, length(cell))){for (j in seq(1,N)) {if(cell[i] ==lb[j]) label[i]<-j}}
cell<- label

rawData = data[,4]
normData = data[,5]
#general look into the data before any process
#boxplot(data[,4]~data[,3], main = "rawData Gene Expressions", xlab = 'cell', ylab ='fpkm')
df = data.frame('cell' = cell, 'expr' = rawData, 'gene' = gene)
#heatmap(as.matrix(df), xlab = 'cell', ylab = 'gene', main = 'heatmap of genes of gc mc')

GCMC<- data.frame('cell1' = df[1:(length(cell)/3),2],'cell2' = df[(length(cell)/3+1):(2*length(cell)/3),2],'cell3' = df[(2*length(cell)/3+1):length(cell),2])
#layout(matrix(c(1,2), 1,2, byrow = TRUE))
boxplot(GCMC, main = "rawData Gene Expressions", xlab = 'cell', ylab ='fpkm')
hm<- heatmap.2(as.matrix(GCMC), xlab = 'cell', ylab = 'gene', main = 'heatmap of genes of gc mc', sideColors = c('darkgreen','yellowgreen'), col = terrain.colors(12))
#dendrogram
hc<- as.hclust(hm$rowDendrogram)
ct<- cutree(hc, h = 6)
data.frame(ct)
colors = c('red','blue','green','yellow','purple','black')
plot(as.phylo(hc), type= 'fan', tip.color = colors[ct], label.offset = 1)

#QC test
d<- density(as.matrix(GCMC))  
hist(as.matrix(GCMC), main = "histogram of rawData Gene Expressions", freq = F)
lines(d)

#check MAplot
logData <- log(as.vector(GCMC), 2)
logData[is.na(logData)] <- max(logData)/1000
M<- logData
MAData<- logData
A<- logData

logData<- cbind(logData,logData[,1])
M<- sapply(seq(1,3), function(i) M[,i]<- (logData[,i+1]-logData[,i]))
A<- sapply(seq(1,3), function(i) A[,i]<- 0.5*(logData[,i+1]+logData[,i]))

stat(logData)
outliers <- boxplot(M, plot=FALSE)$out
M[M[,3] %in% outliers[outliers>0],3]<- max(M[,1])
M[M[,2] %in% outliers[outliers>0],2]<- max(M[,1])
M[M[,3] %in% outliers[outliers<0],3]<- min(M[,1])
M[M[,2] %in% outliers[outliers<0],2]<- min(M[,1])

outliers <- boxplot(A, plot=FALSE)$out
A[A[,2] %in% outliers[outliers<0],2]<- min(A[,1])
A[A[,3] %in% outliers[outliers<0],3]<- min(A[,1])
A[A[,2] %in% outliers[outliers>0],2]<- max(A[,1])
A[A[,3] %in% outliers[outliers>0],3]<- max(A[,1])

#MA Plot 
smoothScatter(A, M, col = 1, main="MA plot ", xlab="A", ylab="M")
abline(h=c(-1,1), col="red")
# according to the result, there is severe bias of M thus need to do the LOESS normalization

#LOESSnorm
Mc<- A[,1]
#MAData<- sapply(seq(1,3), function(i) M[,i]<- (M[,i] - loess(M[,i]~A[,i])))
MAData<- sapply(seq(1,3), function(i) {l<- loess(M[,i]~A[,i])
Mc<- predict(l, A[,i])
Mc[is.na(Mc)]<- 0
print(Mc)
return(MAData[,i]<- (M[,i] - Mc))})
logData<- log2(GCMC)

#QC test
d<- density(MAData)  
hist(MAData, , main = "histogram of LOESSNormData Gene Expressions", freq = F)
lines(d) 
kurtosis.norm.test(MAData)

##################################################################

        Kurtosis test for normality

data:  MAData
T = 5.3885, p-value = 0.0035
######################################################################
boxplot(MAData, main = "MAData Gene Expressions", xlab = 'cell', ylab ='fpkm')
hm<- heatmap.2(as.matrix(MAData), xlab = 'cell', ylab = 'gene', main = 'heatmap of genes of gc mc', sideColors = c('darkgreen','yellowgreen'), col = terrain.colors(12))
#dendrogram
hc<- as.hclust(hm$rowDendrogram)
ct<- cutree(hc, h = 6)
data.frame(ct)
colors = c('red','blue','green','yellow','purple','black')
plot(as.phylo(hc), type= 'fan', tip.color = colors[ct], label.offset = 0.2)


#generate more synthetic data
test_data<- t(MAData)
data0<- test_data
data = data0
test_data<- cbind(data0, c(1,2,3))
for (i in seq(1,3)){
	mu<- mean(data0[i,])
	sigma<- sqrt(var(data0[i,]))
#	h<-hist(data0[i,])
	set.seed(123)
	n<- 0
	while (n < 24){
#		data<- dnorm(h$density, mu, sigma)
		data<- rnorm(25, mu, sigma)
		data<- cbind(t(data), i)
		test_data<- rbind(test_data, as.numeric(unlist(data)))
		print(n)
		n<- n+1
	}		
}
data1 <- rbind(test_data[1,],test_data[4:27,])
data2 <- rbind(test_data[2,],test_data[28:51,])
data3 <- rbind(test_data[3,],test_data[52:75,])
dataset <- rbind(data1, data2)
dataset<- rbind(dataset, data3)
MAData2<- dataset[,-26]

d2<- density(MAData2)  
hist(MAData2, , main = "histogram of LOESSNormData Gene Expressions", freq = F)
lines(d2) 
kurtosis.norm.test(MAData2)

boxplot(MAData2, main = "MAData_Extented Gene Expressions", xlab = 'cell', ylab ='fpkm')
hm<- heatmap.2(as.matrix(MAData2), xlab = 'cell', ylab = 'gene', main = 'heatmap of genes(synthetic) of gc mc', sideColors = c('darkgreen','yellowgreen'), col = terrain.colors(12))
#dendrogram
hc<- as.hclust(hm$rowDendrogram)
ct<- cutree(hc, h = 6)
data.frame(ct)
colors = c('red','blue','green','yellow','purple','black')
plot(as.phylo(hc), type= 'fan', tip.color = colors[ct], label.offset = 0.2)

#
dg_d_mean2 <- colMeans(data1)
dg_v_mean2 <- colMeans(data2)
mc_mean2 <- colMeans(data3)
#
dg_d_mean <- mean(MAData2[, 1])
dg_v_mean <- mean(MAData2[, 2])
mc_mean <- mean(MAData2[, 3])

dg_mean_dtov<- dg_d_mean - dg_v_mean
cell_mean_dgdtomc<- dg_d_mean - mc_mean
cell_mean_dgvtomc<- dg_v_mean - mc_mean

dg_mean_dtov2<- dg_d_mean2 - dg_v_mean2
cell_mean_dgdtomc2<- dg_d_mean2 - mc_mean2
cell_mean_dgvtomc2<- dg_v_mean2 - mc_mean2

all_data <- cbind(dg_d_mean, dg_v_mean, mc_mean, dg_mean_dtov, cell_mean_dgdtomc, cell_mean_dgvtomc)
all_data2 <- cbind(dg_d_mean2, dg_v_mean2, mc_mean2, dg_mean_dtov2, cell_mean_dgdtomc2, cell_mean_dgvtomc2)
write.table(all_data, file = "MA_mean_ALL.txt", quote=F, sep="\t")
write.table(all_data2, file = "MA_mean_ALL_extented.txt", quote=F, sep="\t")

#see if there are any differentially expressed genes
#T-test
t_test_dg <- t.test(MAData[, 1], MAData[, 2], "two.sided")
t_test_dgd_dgv <- t.test(dg_d_mean2,dg_v_mean2, "two.sided")
####################################################
  Welch Two Sample t-test

data:  MAData[, 1] and MAData[, 2]
t = -0.034779, df = 45.019, p-value = 0.9724
alternative hypothesis: true difference in means is not equal to 0
95 percent confidence interval:
 -0.6319860  0.6105304
sample estimates:
  mean of x   mean of y 
-0.02709765 -0.01636982 



        Welch Two Sample t-test

data:  dg_d_mean2 and dg_v_mean2
t = -0.56043, df = 39.924, p-value = 0.5783
alternative hypothesis: true difference in means is not equal to 0
95 percent confidence interval:
 -0.2508925  0.1419621
sample estimates:
 mean of x  mean of y 
0.03129729 0.08576246 
####################################################

t_test_dgd_mc <- t.test(MAData[, 1], MAData[, 3], "two.sided")
t_test_dgd_mc <- t.test(dg_d_mean2,  mc_mean2, "two.sided")
####################################################
  Welch Two Sample t-test

data:  MAData[, 1] and MAData[, 3]
t = -0.24874, df = 47.999, p-value = 0.8046
alternative hypothesis: true difference in means is not equal to 0
95 percent confidence interval:
 -0.5994478  0.4674583
sample estimates:
  mean of x   mean of y 
-0.02709765  0.03889710 



        Welch Two Sample t-test

data:  dg_d_mean2 and mc_mean2
t = -1.1141, df = 33.441, p-value = 0.2732
alternative hypothesis: true difference in means is not equal to 0
95 percent confidence interval:
 -0.3964088  0.1157955
sample estimates:
 mean of x  mean of y 
0.03129729 0.17160398 
####################################################

t_test_dgv_mc <- t.test(MAData[, 2], MAData[, 3], "two.sided")
t_test_dgv_mc <- t.test(dg_v_mean2,  mc_mean2, "two.sided")

####################################################
        Welch Two Sample t-test

data:  MAData[, 2] and MAData[, 3]
t = -0.17943, df = 44.941, p-value = 0.8584
alternative hypothesis: true difference in means is not equal to 0
95 percent confidence interval:
 -0.6756662  0.5651324
sample estimates:
  mean of x   mean of y 
-0.01636982  0.03889710 

        Welch Two Sample t-test

data:  dg_v_mean2 and mc_mean2
t = -0.598, df = 45.58, p-value = 0.5528
alternative hypothesis: true difference in means is not equal to 0
95 percent confidence interval:
 -0.3748602  0.2031771
sample estimates:
 mean of x  mean of y 
0.08576246 0.17160398 
####################################################

#AP <- apply(MAData, 1, paste, collapse="")
#genesP = names(AP[AP != "AAAAAAAA"])
#length(genesP)
#MAData <- MAData[genesP,]
 
#adjust for the FDR, or false discovery rate(FDR correction is necessary when looking at many thousands of tests as a way of minimizing the
#false positives that inevitably are found when so many samples are analyzed, e.g. 5% of the time a
#95% confidence interval is wrong)

t_test_dgd_dgv_p <- t_test_dgd_dgv$p.value
t_test_dgd_mc_p <- t_test_dgd_mc$p.value
t_test_dgv_mc_p <- t_test_dgv_mc$p.value

#The p values are automatically adjusted to compensate for the false discovery rate 
t_test_dg_fdr_p_valsP <- p.adjust(t_test_dgd_dgv_p, method="fdr")
t_test_dgd_mc_fdr_p_valsP <- p.adjust(t_test_dgd_mc_p, method="fdr")
t_test_dgv_mc_fdr_p_valsP <- p.adjust(t_test_dgv_mc_p, method="fdr")

t_test_dg_fdr_p_valsP_sorted <- t_test_dg_fdr_p_valsP[order(t_test_dg_fdr_p_valsP)]
t_test_dgd_mc_fdr_p_valsP_sorted <- t_test_dgd_mc_fdr_p_valsP[order(t_test_dgd_mc_fdr_p_valsP)]
t_test_dgv_mc_fdr_p_valsP_sorted <- t_test_dgv_mc_fdr_p_valsP[order(t_test_dgv_mc_fdr_p_valsP)]
t_test_dg_DE_probesets <- names(t_test_dgd_dgv_p[t_test_dgd_dgv_p < 0.1])
t_test_dgd_mc_DE_probesets <- names(t_test_dgd_mc_p[t_test_dgd_mc_p < 0.1])
t_test_dgv_mc_DE_probesets <- names(t_test_dgv_mc_p[t_test_dgv_mc_p < 0.1])

t_test_dg_DE_log2ratios <- all_data2[t_test_dg_DE_probesets, c(1, 2)]
t_test_dgd_mc_DE_log2ratios <- all_data2[t_test_dgd_mc_DE_probesets, c(1, 3)]
t_test_dgv_mc_DE_log2ratios <- all_data2[t_test_dgv_mc_DE_probesets, c(2, 3)]

dg_d<- data1
dg_v <- data2
mc <- data3

x1data<- MAData[, 1]
x2data<- MAData[, 2]
x3data<- MAData[, 3]
layout(matrix(c(1,2,3)),3,1)
plot(x1data, x2data, main = "Log2 expression in dg_d vs dg_v", xlab="dg_d", ylab="dg_v", col="blue", type = 'p', cex=2)
abline(0,1)
plot(x2data, x3data, main = "Log2 expression in dg_v vs mc", xlab="dg_v", ylab="mc", col="red", type = 'p', cex=2)
abline(0,1)
plot(x1data, x3data, main = "Log2 expression in dg_d vs mc", xlab="dg_v", ylab="mc", col="black", type = 'p', cex=2)
abline(0,1)

layout(matrix(c(1,2,3)),3,1)
plot(dg_d, dg_v, main = "Log2 expression in dg_d vs dg_v", xlab="dg_d", ylab="dg_v", col="blue", type = 'p', cex=2)
abline(0,1)
plot(dg_v, mc, main = "Log2 expression in dg_v vs mc", xlab="dg_v", ylab="mc", col="red", type = 'p', cex=2)
abline(0,1)
plot(dg_d, mc, main = "Log2 expression in dg_d vs mc", xlab="dg_v", ylab="mc", col="black", type = 'p', cex=2)
abline(0,1)

M1<- x2data - x1data
A1<- 0.5*(x2data + x1data)
M2<- x3data - x2data
A2<- 0.5*(x2data + x3data)
M2<- x3data - x1data
A3<- 0.5*(x3data + x1data)
layout(matrix(c(1,2,3)),3,1)
plot(A1, M1, main = "MA in dg_d vs dg_v", xlab="dg_d_A", ylab="dg_v_M", col="blue", pch = 19, cex=2)
abline(h =c( -1,1))
plot(A2, M2, main = "MA in dg_v vs mc", xlab="dg_v_A", ylab="mc_M", col="red", pch = 19, cex=2)
abline(h =c( -1,1))
plot(A1, M1, main = "MA in dg_d vs mc", xlab="dg_v_A", ylab="mc_M", col="black", pch = 19, cex=2)
abline(h =c( -1,1))

M1<- dg_d- dg_v
A1<- 0.5*(dg_d + dg_d)
M2<- mc - dg_v
A2<- 0.5*(mc + dg_v)
M3<- mc - dg_d
A3<- 0.5*(mc + dg_d)

layout(matrix(c(1,2,3)),3,1)
plot(A1, M1, main = "MA in dg_d vs dg_v", xlab="dg_d_A", ylab="dg_v_M", col="blue", pch = 19, cex=2)
abline(h =c( -1,1))
plot(A2, M2, main = "MA in dg_v vs mc", xlab="dg_v_A", ylab="mc_M", col="red", pch = 19, cex=2)
abline(h =c( -1,1))
plot(A3, M3, main = "MA in dg_d vs mc", xlab="dg_d_A", ylab="mc_M", col="black", pch = 19, cex=2)
abline(h =c( -1,1))


#volcano plot
expr_pvals<- cbind(MAData2, t_test_dgd_dgv_p, t_test_dgd_mc_p, t_test_dgv_mc_p, t_test_dg_fdr_p_valsP, t_test_dgd_mc_fdr_p_valsP, t_test_dgv_mc_fdr_p_valsP)
log2_ratios1<- expr_pvals[, 3] - expr_pvals[, 2]
log2_ratios2<- expr_pvals[, 4] - expr_pvals[, 3]
log2_ratios3<- expr_pvals[, 4] - expr_pvals[, 2]
log2_ratios1_fdr<- expr_pvals[, 5] - expr_pvals[, 4]
log2_ratios2_fdr<- expr_pvals[, 6] - expr_pvals[, 5]
log2_ratios3_fdr<- expr_pvals[, 6] - expr_pvals[, 4]
p_values1<- expr_pvals[, 2]
p_values1_fdr<- expr_pvals[, 4]
p_values2<- expr_pvals[, 3]
p_values2_fdr<- expr_pvals[, 5]
p_values3<- expr_pvals[, 4]
p_values3_fdr<- expr_pvals[, 6]

layout(matrix(c(1,2,3,4,5,6,7,7),4,2, byrow = TRUE))
plot(log2_ratios1, p_values1)
plot(log2_ratios1, -log(p_values1, 10))
plot(log2_ratios2, p_values2)
plot(log2_ratios2, -log(p_values2, 10))
plot(log2_ratios3, p_values3)
plot(log2_ratios3, -log(p_values3, 10))
plot(c(log2_ratios1_fdr,log2_ratios1_fdr,log2_ratios2_fdr, log2_ratios2_fdr, log2_ratios3_fdr, log2_ratios3_fdr),c(p_values1_fdr, -log(p_values1_fdr, 10), p_values2_fdr, -log(p_values2_fdr, 10), p_values3_fdr, -log(p_values3_fdr, 10)), col = c(1,2,3,4,5,6))

#plot(log2_ratios1_fdr, p_values1_fdr, col = 1)
#points(log2_ratios1_fdr, -log(p_values1_fdr, 10), col = 2)
#points(log2_ratios2_fdr, p_values2_fdr, col = 3)
#points(log2_ratios2_fdr, -log(p_values2_fdr, 10), col = 4)
#points(log2_ratios3_fdr, p_values1_fdr, col = 5)
#points(log2_ratios3_fdr, -log(p_values1_fdr, 10), col = 6)

#MDS plot 

dataorig<- read.csv('D:/TUT/Medical/biophysics/experiment/WiML/raw_data_DGinMCout.tsv',sep ='	', header= TRUE)
gene <- cell2label(data[,2])
cell <- cell2label(data[,3])

#MDSplot – this stands for Multidimensional Scaling, and is
#essentially a visual representation of a principle components analysis.  This means that it can be used
#to see which individual elements contribute most strongly to any variation seen within the data.
par(mfrow=c(2,1))
col_cells <- c("blue","red","green")
plotMDS(MAData2, cex = 1.8, col = col_cells) 
legend("topleft",fill=c("blue","red","green"),legend=c(dg_d, dg_v, mc),cex=0.6)
title("Cell type")
col_genes <- seq(1,length(MAData2[,1]))
plotMDS(t(MAData2),cex= 1, col=col_genes)
legend("bottomright",fill=col_genes, legend=levels(seq(1,dataorig[,2]),cex=0.6)
title("Gene type")

#dimensionality
par(mfrow=c(2,1))
char_celltype <- c("X", "O", "*")
plotMDS(t(MAData2),col=col_genes,pch=char_celltype,cex=2)
legend("topright",legend=levels(dataorig[,3]),pch=char_celltype)
legend("bottomleft",legend=levels(dataorig[,2]),col=col_genes, pch= 19, cex= 0.4)
title("MDS Plot for first two dimension")
plotMDS(t(MAData2),dim=c(1,3), col=col_genes,pch=char_celltype,cex=2)
legend("topright",legend=levels(dataorig[,3]),pch=char_celltype)
legend("bottomleft",legend=levels(dataorig[,2]),col=col_genes, pch= 19, cex= 0.4)
title("MDS Plot for dimensions 1 and 3") 

var_genes <- apply(MAData2, 1, var)
select_var <-(sort(var_genes, decreasing=FALSE, index.return = TRUE)$ix)[1:20]
highly_variable_lcpm <- MAData2[select_var,]
dim(highly_variable_lcpm)
#mypalette <- brewer_pal(11,"RdYlBu")
#morecols <- colorRampPalette
#generate a series of interpolated colors that will give us a nice full spectrum of colors that fall
between the 11 we chose for my palette – this is crucial for displaying data that is "analog", where
there are essentially an infinite number of values
#heatmap.2(highly_variable_lcpm,col=rev(morecols(50)),trace="none", main="Top 20 mostvariable genes across samples", ColSideColors = col_cells,scale="row")
#hm<- heatmap.2(as.matrix(highly_variable_lcpm), xlab = 'cell', ylab = 'gene', main="Top 20 mostvariable genes across samples", ColSideColors = c('darkgreen','yellowgreen','green'), col = terrain.colors(12))
hm<- heatmap.2(as.matrix(highly_variable_lcpm), xlab = 'cell', ylab = 'gene', main="Top 20 mostvariable genes across samples",  col = terrain.colors(12))

#dendrogram
hc<- as.hclust(hm$rowDendrogram)
ct<- cutree(hc, h = 2)
data.frame(ct)
colors = c('red','blue')
plot(as.phylo(hc), type= 'fan', tip.color = colors[ct], label.offset = 0.1, cex= 1)
png(file="High_var_genes.heatmap.png")
dev.off()



test_data<- t(MAData)
data0<- test_data
data = data0
test_data<- cbind(data0, c(1,2,3))
for (i in seq(1,3)){
	mu<- mean(data0[i,])
	sigma<- sqrt(var(data0[i,]))
#	h<-hist(data0[i,])
	set.seed(123)
	n<- 0
	while (n < 24){
#		data<- dnorm(h$density, mu, sigma)
		data<- rnorm(25, mu, sigma)
		data<- cbind(t(data), i)
		test_data<- rbind(test_data, as.numeric(unlist(data)))
		print(n)
		n<- n+1
	}		
}
data1 <- rbind(test_data[1,],test_data[4:27,])
data2 <- rbind(test_data[2,],test_data[28:51,])
data3 <- rbind(test_data[3,],test_data[52:75,])
dataset <- rbind(data1, data2)
dataset <- rbind(dataset, data3)

#LDA(binary)


# classify between dg_d and dg_v; dg_d and mc

LDA_pred<- function(test_data, test_id, splitrate) {
smp_size <- floor(splitrate * nrow(test_data))
train_ind <- sample(nrow(test_data), smp_size)
tr_df <- as.data.frame(test_data[train_ind, ])
if (splitrate == 1)
	{te_df<- tr_df}
else
	{te_df <- as.data.frame(test_data[-train_ind, ])}
print(paste("ID = ", train_ind, test_id))

f <- paste(names(tr_df)[26], "~", paste(names(tr_df)[-26], collapse=" + "))
Mlda <- lda(as.formula(paste(f)), data = tr_df)

pred <- predict(Mlda, te_df)

### CONSTRUCTING ROC AUC PLOT:
# Get the posteriors as a dataframe.
y <- as.data.frame(pred$posterior)
# Evaluate the model
pred <- prediction(pred$x, te_df[26])
roc = performance(pred, measure = "tpr", x.measure = "fpr")
#acc = performamce(pred, measure = "acc")
#acc = (y[,1]< y[,2]) == (te_df[26] == 2)
#auc_train <- performance(pred, measure = "auc")
#auc_train <- auc_train@y.values
AUC<- as.numeric(unlist(roc@alpha.values[[1]]))
AUC[is.na(AUC)]<-0

y1<- as.numeric(unlist(pred@predictions))
y2<- as.numeric(unlist(pred@labels))
y1[y1 > 0]<- 2
y1[y1 < 0]<- 1 

# Plot
layout(matrix(c(1,2),2,1, byrow = TRUE))
plot(y1, col = c(1,2), type = 'p')
#plot(y1[y1 == 2], col = 1, type = 'p')
#points(y1[y1 == 1], col = 2)
plot(roc@y.values[[1]], type = 'b')
abline(a=0, b= 1)
text(x = .25, y = .65 ,paste("AUC = ", round(AUC[2],3), sep = ""))
#result <- [roc, y1, train_ind]
#return(test_id[-which(test_id == train_ind)],y1)
return(train_ind)
}



test_data <- rbind(data1, data2)
id1<- c(1, seq(4,27))
#Metrice1<- LDA_pred(test_data)
Tr1_id<- LDA_pred(test_data, id1, 0.8)
#Te1_id<- id1[-Tr1_id]
Tr11_id<- LDA_pred(test_data[seq(4,13)], id1, 1)

test_data_dgdmc <- rbind(data1, data3)
#id2<- c(2, seq(28,37))
Tr2_id<- LDA_pred(test_data_dgdmc, id2, 0.8)
#Tr2_id<- LDA_pred(test_data_dgdmc)
#Metrice2<- id2[-Tr2_id]
Tr12_id<- LDA_pred(test_data[seq(28,51)], id2, 1)

test_data_dgvmc <- rbind(data2, data3)
#id3<- c(3, seq(52,75))
Tr3_id<- LDA_pred(test_data_dgdmc, id3, 0.8)
Metrice3<- LDA_pred(test_data_dgvmc)
#Te3_id<- id3[-Tr3_id]
Tr13_id<- LDA_pred(test_data[seq(56,65)], id3, 1)



class<-rbind(test_data, data3)[,26] 
y<- c(1,4,5,11,15,16,20,22,41,31,28,29,32,33,36,37,39,40,43,47,56,2,12,58,60,67,68,70,74,75)
x<- c(seq(4,13),seq(28,37), seq(56,65))
Class<- class[y]
Class2<- c(1,2,1,2,2,2,2,1,1,1,3,3,3,3,3,3,3,3,3,3,3,2,2,2,2,2,2,2,2,2)
plotdata<-data.frame(Class, x, y)
plotdata2<-data.frame(Class2, x, y)
centroids <- aggregate(cbind(x,y)~Class,plotdata,mean)
centroids2 <- aggregate(cbind(x,y)~Class2,plotdata2,mean)
CentroidDistances <- dist(centroids, method = "euclidean", diag = TRUE, upper = FALSE, p = 2)
CentroidDistances2 <- dist(centroids2, method = "euclidean", diag = TRUE, upper = FALSE, p = 2)
attr(CentroidDistances, "Labels") <- centroids$Class
attr(CentroidDistances2, "Labels") <- centroids2$Class2


centroids<- centroids[,1]
layout(matrix(c(1, 2), 2, 1, byrow = TRUE))
ggplot(plotdata,aes(x,y,color=factor(Class))) + geom_point(size=3)+ geom_point(data=centroids,size=7) + geom_text(data=centroids, size=7, label= c(1,2,3), colour="black")
ggplot(plotdata2,aes(x,y,color=factor(Class2))) + geom_point(size=3)+ geom_point(data=centroids2,size=7) + geom_text(data=centroids2, size=7, label = c(1, 2, 3), colour="black")
#multiplot(p1, p2)

