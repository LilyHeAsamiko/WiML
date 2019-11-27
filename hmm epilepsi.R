getwd()
ls()
setwd("D:/Introduction to R/HMM")
setwd("P:/Introduction-to-R-master/Introduction-to-R-master")

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("STAN")

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("mgsa")

bioManager.install("BSgenome")
library(BSgenome)
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
#install.packages("ROCR")
library(ROCR)
library(STAN)
install.packages("permute")
library(permute)
install.packages("vegan")
library(vegan)

data<- read.table("D:/Introduction to R/HMM/Dock7.txt",header = TRUE, sep=",",stringsAsFactors=FALSE, quote="")
colnames(data)
data1<- read.fasta(file = "D:/Introduction to R/HMM/DOCK7seq.fa",  as.string = TRUE, seqtype = "AA")[[1]]
#  stopifnot( data1== "SEQINRSEQINRSEQINRSEQINR*") # ret = 'sequence',
#data1<- system.file('D:/Introduction to R/HMM/DOCK7.fa', package = 'seqinr')
seq1 <- data1[1]

data2<- read.fasta(file = "D:/Introduction to R/HMM/STXBP1seq.fa",  as.string = TRUE, seqtype = "AA")[[1]]
seq2 <-data2[1]
vnames <- c('seqname','source','feature','start','end','score','strand','frame','hid','hstart','hend','genscan','gene_id','transcript_id','exon_id','gene_type','variation_name','probe_name')
data3<- read.table("D:/Introduction to R/HMM/epilepsy STXBP1 phenotype 1.txt", sep=",", na.strings=c("", "NA"))



snp1_dis <- data3
colnames(snp1_dis) <- vnames
Snp1_dis<- snp1_dis

Snp1_dis$feature<- seq(1,dim(Snp1_dis)[1])
Snp1_dis<- Snp1_dis[,colSums(is.na(snp1_dis))== 0] 

#dependent, bidirectional (ex: one motif: 5'-binding: TGCAAG, 5'-splice: CAAGTG and AAGTGC)
#
#snp1<- aregexec("TGCAAG", seq1, max.distance = 2)

snp1<- gregexpr("TGCA*", seq1, perl = TRUE)
str_sub(seq1, 188, 193)

snpstxbp1.motif<- snp1[[1]]
snpstxbp1.len<-attr(snp1[[1]],"match.length")

#write.fasta(as.list(snp1_dis), names = 'STXBP1pheno', file.out="D:/Introduction to R/HMM/epilepsy STXBP1 phenotype.fasta")
write.table(seq1, file  = 'STXBP1pheno.txt', sep = ",")

snp1_00 <- Snp1_dis[Snp1_dis$frame == -1,]
snp1_10 <- Snp1_dis[Snp1_dis$frame == 0,]
snp1_11 <- Snp1_dis[Snp1_dis$frame == 1,]
snp1_12 <- Snp1_dis[Snp1_dis$frame == 2,]
N100<- dim(snp1_00)[1] 
N110<- dim(snp1_10)[1] 
N111<- dim(snp1_11)[1] 
N112<- dim(snp1_12)[1] 
N1<- N100+N110+N111+N112 

steps<- 50

#snp1_00[1:50]

data4<- read.table(file = "D:/Introduction to R/HMM/epilepsy Dock7 phenotype 0.txt", sep=",", na.strings=c("", "NA"), header = FALSE)

snp2_dis <- data4[2:dim(data4)[1],]
colnames(snp2_dis) <- vnames
Snp2_dis<- snp2_dis
Snp2_dis$feature<- seq(1, dim(Snp2_dis)[1])
Snp2_dis<- Snp2_dis[,colSums(is.na(snp2_dis))== 0] 
Snp2_dis$frame<- as.numeric(Snp2_dis$frame)-2

snp2_00 <- Snp2_dis[Snp2_dis$frame == -1,]
snp2_10 <- Snp2_dis[Snp2_dis$frame == 0,]
snp2_11 <- Snp2_dis[Snp2_dis$frame == 1,]
snp2_12 <- Snp2_dis[Snp2_dis$frame == 2,]
N200<- dim(snp2_00)[1] 
N210<- dim(snp2_10)[1]
N211<- dim(snp2_11)[1] 
N212<- dim(snp2_12)[1]
N2<- N200+N210+N211+N212 

# k = 0,1, i,j = -1,0,1,2 T = 6, iter = 50
#independent, contingency table
eps1<- 0.01
w <- ones(9, 1)./9;
lmd <- w;

steps<- 50
N1_1A<- str_count(str_sub(seq1,snp1_11$feature[1],max(snp1_11$feature[snp1_11$feature<= min(nchar(seq1))])),'A')
N1_1C<- str_count(str_sub(seq1,snp1_11$feature[1],max(snp1_11$feature[snp1_11$feature<= min(nchar(seq1))])),'C')
N1_1T<- str_count(str_sub(seq1,snp1_11$feature[1],max(snp1_11$feature[snp1_11$feature<= min(nchar(seq1))])),'T')
N1_1G<- str_count(str_sub(seq1,snp1_11$feature[1],max(snp1_11$feature[snp1_11$feature<= min(nchar(seq1))])),'G')

N1_2A<- str_count(str_sub(seq1,snp1_12$feature[1],max(snp1_12$feature[snp1_12$feature<= min(nchar(seq1))])),'A')
N1_2C<- str_count(str_sub(seq1,snp1_12$feature[1],max(snp1_12$feature[snp1_12$feature<= min(nchar(seq1))])),'C')
N1_2T<- str_count(str_sub(seq1,snp1_12$feature[1],max(snp1_12$feature[snp1_12$feature<= min(nchar(seq1))])),'T')
N1_2G<- str_count(str_sub(seq1,snp1_12$feature[1],max(snp1_12$feature[snp1_12$feature<= min(nchar(seq1))])),'G')

kse<- c(N110/N1*N210/N2, (N111+N100)/N1*N210/N2, N112/N1*N210/N2, N110/N1*N211/N2, (N100+N111)/N1*N211/N2, N112/N1*N211/N2, (N100+N111)/N1*N212/N2, (N100+N111)/N1*N212/N2, N112/N1*N212/N2) 
#1.additive
alphaa0sq<- kse[7:9]/kse[1:3]
alphaa0<- mean(sqrt(alphaa0sq))
betaa0sq<- kse[c(3,6,9)]/kse[c(1,4,7)]
betaa0<- mean(sqrt(betaa0sq))

#2.multiplicative
betab0<- 0.5*sqrt(kse[3]/kse[1])+0.5*kse[2]/kse[1]
alphab0<- 0.5*sqrt(kse[4]/kse[1])+0.5*kse[7]/kse[1]
deltab0<-0.5*sqrt(kse[8]/kse[7]/(kse[5]/kse[4]))+0.5*(kse[9]/kse[8]/(kse[6]/kse[5]))

N<- N1+N2
ksec<- matrix(0, 9, floor(min(N, snpstxbp1.motif[length(snpstxbp1.motif)-4])/6))
ksea<- ksec
ksem<- ksec
c<- ksec
a<- ksec
m<- ksec
rc<-ksec
ra<- ksec
rm<- ksec

lambda1<- matrix(0, 1, floor(min(N, snpstxbp1.motif[length(snpstxbp1.motif)-4])/6))
lambda2<- lambda1
P00<- lambda1
P01<- lambda1
P10<- lambda1
P11<- lambda1
mutation<-lambda1
chi2c<- lambda1
chi2a<- lambda1
chi2m<- lambda1

i<-1

#i, j = -1,0,1,2 k = 0,1 
N<- min(N1, N2)
while (i < min(N, snpstxbp1.motif[length(snpstxbp1.motif)-4])) {
#?	temp1<- Snp1_dis$frame[i:i+5] 
	temp1<- Snp1_dis$frame[seq(i, i+5)]
	temp2<- Snp2_dis$frame[seq(i, i+5)]
	if (sum(snpstxbp1.motif >= i & snpstxbp1.motif <= i+5)>0){
		mutation[floor(i/6)+1]<- 1}
	#control
	temp_000<- sum(temp1 == 0 & temp2 == 0 & mutation[floor(i/6)+1]==1) 
	temp_001<- sum(temp1 == 0 & temp2 == 0 & mutation[floor(i/6)+1]==0)
	tp_00<- sum(temp1 == 0 & temp2 == 0) 
	temp_010<- sum(temp1 == 0 & abs(temp2) == 1 & mutation[floor(i/6)+1]==1)
	temp_011<- sum(temp1 == 0 & abs(temp2) == 1 & mutation[floor(i/6)+1]==0)
	tp_01<- sum(temp1 == 0 & temp2 == 1) 
	temp_020<- sum(temp1 == 0 & temp2 == 2 & mutation[floor(i/6)+1]==1) 
	temp_021<- sum(temp1 == 0 & temp2 == 2 & mutation[floor(i/6)+1]==0)
	tp_02<- sum(temp1 == 0 & temp2 == 2) 
	temp_100<- sum(abs(temp1) == 1 & temp2 == 0 & mutation[floor(i/6)+1]==1) 
	temp_101<- sum(abs(temp1) == 1 & temp2 == 0 & mutation[floor(i/6)+1]==0)
	tp_10<- sum(temp1 == 1 & temp2 == 0) 
	temp_110<- sum(abs(temp1) == 1 & abs(temp2) == 1 & mutation[floor(i/6)+1]==1)
	temp_111<- sum(abs(temp1) == 1 & abs(temp2) == 1 & mutation[floor(i/6)+1]==0)
	tp_11<- sum(temp1 == 1 & temp2 == 1)
	temp_120<- sum(abs(temp1) == 1 & temp2 == 2 & mutation[floor(i/6)+1]==1) 
	temp_121<- sum(abs(temp1) == 1 & temp2 == 2 & mutation[floor(i/6)+1]==0)
	tp_12<- sum(temp1 == 1 & temp2 == 2)
	temp_200<- sum(temp1 == 2 & temp2 == 0 & mutation[floor(i/6)+1]==1) 
	temp_201<- sum(temp1 == 2 & temp2 == 0 & mutation[floor(i/6)+1]==0)
	tp_20<- sum(temp1 == 2 & temp2 == 0)
	temp_210<- sum(temp1 == 2 & abs(temp2) == 1 & mutation[floor(i/6)+1]==1)
	tp_21<- sum(temp1 == 2 & temp2 == 1)
	temp_220<- sum(temp1 == 2 & temp2 == 2 & mutation[floor(i/6)+1]==1) 
	temp_221<- sum(temp1 == 2 & temp2 == 2 & mutation[floor(i/6)+1]==0)
	tp_22<- sum(temp1 == 2 & temp2 == 2)

	#odds
	print(floor(i/6)+1)
	print(min(N, snpstxbp1.motif[length(snpstxbp1.motif)-4]))
	ksec[,floor(i/6)+1] <- c((temp_001+1)/(temp_000+1), (temp_011+1)/(temp_010+1), (temp_021+1)/(temp_020+1), (temp_101+1)/(temp_100+1), (temp_111+1)/(temp_110+1), (temp_121+1)/(temp_120+1), (temp_201+1)/(temp_200+1), (temp_211+1)/(temp_210+1), (temp_221+1)/(temp_220+1))
	ksea[,floor(i/6)+1] <- c((temp_001+1)/(temp_000+1), (betaa0*temp_011+0.0001)/(temp_010+0.0001), (betaa0^2*temp_021+0.0001)/(temp_020+0.0001), (alphaa0*temp_101+0.0001)/(temp_100+0.0001), (alphaa0*betaa0*temp_111+0.0001)/(temp_110+0.0001), (alphaa0*betaa0^2*temp_121+1)/(temp_120+1), (alphaa0^2*temp_201+1)/(temp_200+1), (alphaa0^2*betaa0*temp_211+1)/(temp_210+1), (alphaa0^2*betaa0^2*temp_221+1)/(temp_220+1))
	ksem[,floor(i/6)+1] <- c((temp_001+0.0001)/(temp_000+0.0001), (betab0*temp_011+0.0001)/(temp_010+0.0001), (betab0^2*temp_021+0.0001)/(temp_020+0.0001), (alphab0*temp_101+0.0001)/(temp_100+0.0001), (alphab0*betab0*deltab0*temp_111+0.0001)/(temp_110+0.0001), (alphab0*betaa0^2*deltab0^2*temp_121+0.0001)/(temp_120+0.0001), (alphab0^2*temp_201+0.0001)/(temp_200+0.0001), (alphab0^2*betab0*deltab0^2*temp_211+0.0001)/(temp_210+0.0001), (alphab0^2*betab0^2*deltab0^4*temp_221+0.0001)/(temp_220+0.0001))		
	
	c[,floor(i/6)+1] <- c((tp_00+0.0001)/(6-tp_00+0.0001), (tp_01+0.0001)/(6-tp_01+0.0001), (tp_02+0.0001)/(6-tp_02+0.0001), (tp_10+0.0001)/(6-tp_10+0.0001), (tp_11+0.0001)/(6-tp_11+0.0001), (6-tp_12+0.0001)/(6-tp_12+0.0001), (tp_20+0.0001)/(6-tp_20+0.0001), (tp_21+0.0001)/(6-tp_21+0.0001), (tp_22+0.0001)/(6-tp_22+0.0001))
	a[,floor(i/6)+1] <- c((tp_00+0.0001)/(6-tp_00+0.0001), (betaa0*tp_01+0.0001)/(6-tp_01+0.0001), (betaa0^2*tp_02+0.0001)/(6-tp_02+0.0001), (alphaa0*tp_10+0.0001)/(6-tp_10+0.0001), (alphaa0*betaa0*tp_11+0.0001)/(6-tp_11+0.0001), (alphaa0*betaa0^2*tp_12+0.0001)/(6-tp_12+0.0001), (alphaa0^2*tp_20+0.0001)/(6-tp_20+0.0001), (alphaa0^2*betaa0*tp_21+0.0001)/(6-tp_21+0.0001), (alphaa0^2*betaa0^2*tp_22+0.0001)/(6-tp_22+0.0001))
	m[,floor(i/6)+1] <- c((tp_00+0.0001)/(6-tp_00+0.0001), (betab0*tp_01+0.0001)/(6-tp_01+0.0001), (betab0^2*tp_02+0.0001)/(6-tp_02+0.0001), (alphab0*tp_10+0.0001)/(6-tp_10+0.0001), (alphab0*betab0*deltab0*tp_11+0.0001)/(6-tp_11+0.0001), (alphab0*betaa0^2*deltab0^2*tp_12+0.0001)/(6-tp_12+0.0001), (alphab0^2*tp_20+0.0001)/(6-tp_20+0.0001), (alphab0^2*betab0*deltab0^2*tp_21+0.0001)/(6-tp_21+0.0001), (alphab0^2*betab0^2*deltab0^4*tp_22+0.0001)/(6-tp_22+0.0001))		

	rc[,floor(i/6)+1] <- ksec[,floor(i/6)+1]/(ksec[,floor(i/6)+1] + c[,floor(i/6)+1]+0.0001)
	ra[,floor(i/6)+1] <- ksea[,floor(i/6)+1]/(ksea[,floor(i/6)+1] + a[,floor(i/6)+1]+0.0001)
	rm[,floor(i/6)+1] <- ksem[,floor(i/6)+1]/(ksem[,floor(i/6)+1] + m[,floor(i/6)+1]+0.0001)
	
	#prevalence	00 01 10 11 
	P00[floor(i/6)+1]<- tp_00/(temp_101 + temp_201 + temp_001 + temp_111 + temp_211 + temp_011 + temp_121 + temp_221 + temp_021 +0.00001)	
	P01[floor(i/6)+1]<- tp_01/(temp_101 + temp_201 + temp_001 + temp_111 + temp_211 + temp_011 + temp_121 + temp_221 + temp_021 +0.00001)
	P10[floor(i/6)+1]<- tp_10/(temp_101 + temp_201 + temp_001 + temp_111 + temp_211 + temp_011 + temp_121 + temp_221 + temp_021 +0.00001)
	P11[floor(i/6)+1]<- tp_11/(temp_101 + temp_201 + temp_001 + temp_111 + temp_211 + temp_011 + temp_121 + temp_221 + temp_021 +0.00001)

	#p-value of chi test
	chi2c[floor(i/6)+1] <- chisq.test(rbind(abs(ksec[,floor(i/6)+1]), abs(log(c[,floor(i/6)+1]))))
	chi2a[floor(i/6)+1] <- chisq.test(rbind(abs(ksea[,floor(i/6)+1]), abs(log(a[,floor(i/6)+1]))))
	chi2m[floor(i/6)+1] <- chisq.test(rbind(abs(ksem[,floor(i/6)+1]), abs(log(m[,floor(i/6)+1]))))

if(floor(i/6)+1< min(N, snpstxbp1.motif[length(snpstxbp1.motif)-10])/6){i<- i+6}else{
break}
}

P00[is.na(P00) or is.infinite(P00)]<-1
P01[is.na(P01) or is.infinite(P01)]<-1
P10[is.na(P10) or is.infinite(P10)]<-1
P11[is.na(P11) or is.infinite(P11)]<-1

y<-mutation
eps0<- 0.01
wc <- rep(1/9, 9)
lmdc <-  wc
wa <- rep(1/9, 9)
lmda <-  wa
wm <- rep(1/9, 9)
lmdm <-  wm
stepsc <- 50
stepsa <- 50
stepsm <- 50

Lc <- matrix(0, dim(rc)[1], dim(rc)[2])
La <- Lc
Lm <- Lc
lc <- rep(0, 9)
la <- rep(0, 9)
lm <- rep(0, 9)
ypredc <- lc
ypreda <- la 
ypredm <- lm

for (kc in seq(1, dim(rc)[2]) & ka in seq(1, dim(ra)[2]) & km in seq(1, dim(rc)[2])) {	
	while (stepsc >=0 | stepsa >= 0 | stepsm >= 0) {
		flag<-1 
		trainc<- ksec[,kc]
		trainc<-(trainc - mean(trainc))/sqrt(var(trainc)+0.00001)		
		traina<- ksea[,ka]
		traina<-(traina - mean(traina))/sqrt(var(traina)+0.00001)
		trainm<- ksem[,km]
		trainm<-(trainm - mean(trainm))/sqrt(var(trainm)+0.00001)
 		scale<- 2
		shape<- 0.2
		Yc<-y[kc]
		Ya<-y[ka]
		Ym<-y[km]

		bc <- rgamma(dgamma(trainc, 0.2, 2), 0.2, 2)
		ba <- rgamma(dgamma(traina, 0.2, 2), 0.2, 2)
		bm <- rgamma(dgamma(trainm, 0.2, 2), 0.2, 2)
		wc <- wc/(sum(wc)+0.000001)
		wa <- wa/(sum(wa)+0.000001)
		wm <- wm/(sum(wm)+0.000001)
		Kc<- abs(trainc)
		Ka<- abs(traina)
		Km<- abs(trainm)
		condc <- log(1 + exp(-Yc* wc[1:length(Kc)]* trainc)) + lmdc[1:length(Kc)]*wc[1:length(Kc)]**2/2 + 2*Kc*bc* trainc/(eps0*seq(1,length(Kc)))> lc
            conda <- log(1 + exp(-Ya* wa[1:length(Ka)]* traina)) + lmda[1:length(Ka)]*wa[1:length(Ka)]**2/2 + 2*Ka*ba* traina/(eps0*seq(1,length(Ka)))> la
		condm <- log(1 + exp(-Ym* wm[1:length(Km)]* trainm)) + lmdm[1:length(Km)]*wm[1:length(Km)]**2/2 + 2*Km*bm* trainm/(eps0*seq(1,length(Km)))> lm
 		if (sum(is.nan(condc))!= length(condc) & sum(condc==0)!= length(condc)){
            	wc[1:length(Kc)] <- -Yc* trainc/(1 + exp(-Yc* wc[1:length(Kc)]* trainc))
		      lc<- log(1 + exp(-Yc* wc[1:length(Kc)]* trainc)) + lmdc[1:length(Kc)]*wc[1:length(Kc)]**2/2 + 2*Kc*bc* trainc /(eps0*seq(1,length(Kc)))}
		if (sum(is.nan(conda))!= length(conda) & sum(conda==0)!= length(conda)){
            	wa[1:length(Ka)] <- -Ya* traina/(1 + exp(-Ya* wa[1:length(Ka)]* traina))
            	la<- log(1 + exp(-Ya* wa[1:length(Ka)]* traina)) + lmda[1:length(Ka)]*wa[1:length(Ka)]**2/2 + 2*Ka*ba* traina /(eps0*seq(1,length(Ka)))}
		if (sum(is.nan(condm))!= length(condm) & sum(condm==0)!= length(condm)){
            	wm[1:length(Km)] <- -Ym* trainm/(1 + exp(-Ym* wm[1:length(Km)]* trainm))
			lm<- log(1 + exp(-Ym* wm[1:length(Km)]* trainm)) + lmdm[1:length(Km)]*wm[1:length(Km)]**2/2 + 2*Km*bm* trainm /(eps0*seq(1,length(Km)))
		}
		if (sum(is.nan(condc))== length(Kc) | sum(condc==0)==  length(Kc)){
			kc <- kc+1
			stepsc <- 50
			flag <- 0}
	      if (sum(is.nan(conda))== length(Ka) | sum(conda==0)==  length(Ka)){
   		      ka <- ka+1
            	stepsa <- 50
			flag <- 0}
	      if (sum(is.nan(condm))== length(Km)| sum(condm==0)==  length(Km)){
 		  	km <- km+1
                  stepsm <- 50
			flag <- 0}
            if (flag== 0) break

            if (stepsc == 0){
			if (sum(is.nan(lc))==length(lc)){
				lc <- 0.5*Lc[,max(kc-1,1)]+0.5*Lc[,min(kc+1,length(Kc))]}
			if (sum(is.nan(lc))>0){
				lc[is.nan(lc)]<- 0.5*lc[max(is.nan(lc)-1,1)]+0.5*lc[min(is.nan(lc)+1,length(Kc))]}
			flag<- 0
			stepcs <- 50
			kc<- kc +1
			if (sum(lc != 0)> length(lc)/2)  ypredc[kc-1] <- 1
			Lc[,kc-1] <- lc}	else{stepsc <- stepsc -1}
            if (stepsa == 0){
			if (sum(is.nan(la))==length(la)){
				la <- 0.5*La[,max(ka-1,1)]+0.5*La[,min(ka+1,length(Ka))]}
			if (sum(is.nan(la))>0 ){
				la[is.nan(la)]<- 0.5*la[max(is.nan(la)-1,1)]+0.5*la[min(is.nan(la)+1,length(Ka))]}
			flag<- 0
			stepsa <- 50
			ka<- ka +1
			if (sum(la != 0)> length(la)/2)  ypreda[ka-1] <- 1
			La[,ka-1] <- la}     
			else{stepsa <- stepsa -1}
            if (stepsm == 0){                	
			if (sum(is.nan(lm))==length(lm)){
				lm <- 0.5*Lm[,max(km-1,1)]+0.5*Lm[,min(km+1,length(Km))]}
                	if (sum(is.nan(lm))>0) {
				lm[is.nan(lm)]<- 0.5*lm[max(is.nan(lm)-1,1)]+0.5*lm[min(is.nan(lm)+1,length(Km))]}
      		flag<- 0  
			stepsm <- 50    
			km<- km +1
			if (sum(lm != 0)> length(lm)/2)  ypredm[km-1] <- 1
			Lm[,km-1] <- lm}     else{stepsm <- stepsm -1}
		if (flag== 0){break}           
	}
}

Pr<-rbind(P00, 0.5*P01, 0.5*P01, 0.5*P10, 0.25*P11, 0.25*P11, 0.5*P10, 0.25*P11,0.25*P11)
ypredm<- sapply(seq(1, 101), function(i) sum((Pr*Lm)[,i]!=0)> dim(Lm)[1]/2)
ypreda<- sapply(seq(1, 101), function(i) sum((Pr*La)[,i]!=0)> dim(La)[1]/2)
ypredc<- sapply(seq(1, 101), function(i) sum((Pr*Lc)[,i]!=0)> dim(Lc)[1]/2)

accc<- sum(ypredc==y)/101
acca<- sum(ypreda==y)/101
accm<- sum(ypredm==y)/101

mr<- sapply(seq(1, 101), function(i) sum(y[1:i]==ypredm[1:i])/i)
ar<- sapply(seq(1, 101), function(i) sum(y[1:i]==ypreda[1:i])/i)
cr<- sapply(seq(1, 101), function(i) sum(y[1:i]==ypredc[1:i])/i)

YMIN <- rbind(cr, ar, mr)-rbind(stdc, stda, stdm)
YMAX <- rbind(cr, ar, mr)+rbind(stdc, stda, stdm)

stdm<- sapply(seq(1, 101), function(i) sqrt(var(mr[1:i])))
stda<- sapply(seq(1, 101), function(i) sqrt(var(ar[1:i])))
stdc<- sapply(seq(1, 101), function(i) sqrt(var(cr[1:i])))

stdc[is.na(stdc)]<- 0
stda[is.na(stda)]<- 0
stdm[is.na(stdm)]<- 0

#sensitivity(tpr = tp/(tp+fn))
TPRc<- sum(y == 1 & ypredc == 1)/(sum(y == 1 & ypredc == 1)+sum(y == 1 & ypredc == 0))
TPRa<- sum(y == 1 & ypreda == 1)/(sum(y == 1 & ypreda == 1)+sum(y == 1 & ypreda == 0))
TPRm<- sum(y == 1 & ypredm == 1)/(sum(y == 1 & ypredm == 1)+sum(y == 1 & ypredm == 0))

#specitficiy(tnr = tn/(tn+fp))
TNRc<- sum(y == 0 & ypredc == 0)/(sum(y == 0 & ypredc == 0)+sum(y == 0 & ypredc == 1))
TNRa<- sum(y == 0 & ypreda == 0)/(sum(y == 0 & ypreda == 0)+sum(y == 0 & ypreda == 1))
TNRm<- sum(y == 0 & ypredm == 0)/(sum(y == 0 & ypredm == 0)+sum(y == 0 & ypredm == 1))

#precision(ppv = tp/(tp+fp))
PPVc<- sum(y == 1 & ypredc == 1)/(sum(y == 1 & ypredc == 1)+sum(y == 0 & ypredc == 1))
PPVa<- sum(y == 1 & ypreda == 1)/(sum(y == 1 & ypreda == 1)+sum(y == 0 & ypreda == 1))
PPVm<- sum(y == 1 & ypredm == 1)/(sum(y == 1 & ypredm == 1)+sum(y == 0 & ypredm == 1))

#miss rate(fnr = fn/(fn+tp))
FNRc<- sum(y == 1 & ypredc == 0)/(sum(y == 1 & ypredc == 0)+sum(y == 0 & ypredc == 0))
FNRa<- sum(y == 1 & ypreda == 0)/(sum(y == 1 & ypreda == 0)+sum(y == 0 & ypreda == 0))
FNRm<- sum(y == 1 & ypredm == 0)/(sum(y == 1 & ypredm == 0)+sum(y == 0 & ypredm == 0))

#false discovety rate(fdr = fp/(tp+fp))
FDRc<- sum(y == 0 & ypredc == 1)/(sum(y == 0 & ypredc == 1)+sum(y == 1 & ypredc == 1))
FDRa<- sum(y == 0 & ypreda == 1)/(sum(y == 0 & ypreda == 1)+sum(y == 1 & ypreda == 1))
FDRm<- sum(y == 0 & ypredm == 1)/(sum(y == 0 & ypredm == 1)+sum(y == 1 & ypredm == 1))



	#parameter
#	nominator1<-(temp_101 + temp_111 + temp_121 + temp_201 + temp_211 + temp_221)/(temp_101 + temp_111 + temp_121 + temp_201 + temp_211 + temp_221 + temp_100 + temp_110 + temp_120 + temp_200 + temp_210 + temp_220)*(temp_000 + temp_010 + temp_020)/(temp_001 + temp_011 + temp_021 + temp_000 + temp_010 + temp_020)
#	denominator1<- (temp_100 + temp_110 + temp_120 + temp_200 + temp_210 + temp_220)/(temp_100 + temp_110 + temp_120 + temp_200 + temp_210 + temp_220 + temp_100 + temp_110 + temp_120 + temp_200 + temp_210 + temp_220)*(temp_001 + temp_011 + temp_021)/(temp_001 + temp_011 + temp_021 + temp_101 + temp_111 + temp_121 + temp_201 + temp_211 + temp_221)
#	nominator2<- (temp_011 + temp_111 + temp_211 + temp_021 + temp_121 + temp_221)/(temp_011 + temp_111 + temp_211 + temp_021 + temp_121 + temp_221 + temp_010 + temp_110 + temp_210 + temp_020 + temp_120 + temp_220)*(temp_000 + temp_010 + temp_200)/(temp_001 + temp_101 + temp_201 + temp_000 + temp_100 + temp_200)
#	denominator2<- (temp_010 + temp_110 + temp_210 + temp_020 + temp_120 + temp_220)/(temp_010 + temp_110 + temp_210 + temp_020 + temp_120 + temp_220 + temp_010 + temp_110 + temp_210 + temp_020 + temp_120 + temp_220)*(temp_001 + temp_101 + temp_201)/(temp_001 + temp_101 + temp_201 + temp_011 + temp_111 + temp_211 + temp_021 + temp_121 + temp_221)
#	lambda1[floor(i/6)+1]<- nominator1/denominator1 -1
#	lambda2[floor(i/6)+1]<- nominator2/denominator2 -1


par(mfrow = c(4,2))
boxplot(t(Pr),main = 'prevalence',names=c('00', '01', '02', '10','11','12', '20', '21', '22'))

#p<- ggplot(rbind(cr, ar, mr), aes(x=motifs, y=accuracy, group=c(control, additive, multiplicative), color=c(control, additive, multiplicative)) + 
#  geom_line() +
#  geom_point()+
#  geom_errorbar(aes(ymin= t(YMIN), ymax = t(YMAX)), width=0.2)	
#print(p)
#  position=position_dodge(0.05)

plot(1:101, cr,
    ylim=range(c(cr-stdc, cr+stdc)),
    pch=19, xlab="Motif", ylab="Mean +/- SD",
    main="Scatter plot with std.dev error bars for control model"
)
# hack: we draw arrows but with very special "arrowheads"
arrows(1:101, cr-stdc, 1:101, cr+stdc, length=0.05, angle=90, code=3)

barplot(t(rbind(accc, acca, accm)), main = 'model accuracies')
barplot(t(rbind(TPRc, TPRa, TPRm)), main = 'sensitivity')
barplot(t(rbind(TNRc, TNRa, TNRm)), main = 'specitficiy')
barplot(t(rbind(PPRc, PPTa, PPTm)), main = 'precision')
barplot(t(rbind(FDRc, FDRa, FDRm)), main = 'fixed modle') 
barplot(t(rbind(FNRc, FNRa, FNRm)), main = 'false discovery rate') 

for(i in seq(1,min(N, length(snpstxbp1.motif)))){
	eps0<- 0.01
	w1 <- rep(1/9, 6)
	lmd1 <-  w1
#	w2 <- rep(1/9, 6)
#	lmd2 <-  w2
#	steps <- 50
#	L1 <- rep(-100,6)
#	L2 <- L1 
#      l <- rep(0, 6*(snpstxbp1.motif[i+1]-snpstxbp1.motif[i]-6))
	#two stage logistic
#	ypred <- rep(0, length= 6, times = snpstxbp1.motif[i+1]-snpstxbp1.motif[i])
#	ypred <- rep(0, 6*(snpstxbp1.motif[i+1]-snpstxbp1.motif[i]-6))
	#regression
#	ypred_2<- ypred
#	W1 <- rep(0, 6*(snpstxbp1.motif[i+1]-snpstxbp1.motif[i])-6)
#	W2 <- W1	
	#chr
#	Y<- Snp1_dis$frame[snpstxbp1.motif[i]:snpstxbp1.motif[i+1]]

#	for (k in seq(snpstxbp1.motif[i], snpstxbp1.motif[i + 1])){
#		id<- seq(k,k+5)
#		y<- Snp1_dis$frame[id]
#		while steps>0{ 
#			train01<- y[y!=2]
#			train01id<- which(y %in% c(-1,0)) 
#			train02id<- which(y %in% c(1,2))
#			train02<- y[train02id]			
#			train01p<-(train01 - mean(train01))/sqrt(var(train01)+0.00001)
#			train02p<-(train02 - mean(train02))/sqrt(var(train02)+0.00001)
#			scale<- 2
#			shape<- 0.2
			
#			b1 <- rgamma(length(train01p), sqrt(var(train01p)),mean(train01p))
#			b1 <- rgamma(length(train01p), shape, scale)
#			b2 <- rgamma(length(train02p), shape, scale)
#			K1<- abs(train01p)
#			K2<- abs(train02p)
#            	l_temp1<- log(1 + exp(-train01* w1[1:length(K1)]* train01p)) + lmd1[1:length(K1)]*w1[1:length(K1)]**2/2 + 2*K1*b1* train01p /(eps0*length(K1))
#            	l_temp2<- log(1 + exp(-train02* w2[1:length(K2)]* train02p)) + lmd2[1:length(K2)]*w2[1:length(K2)]**2/2 + 2*K2*b1* train02p /(eps0*length(K2))
#			cond1 <- log(1 + exp(-train01* w1[1:length(K1)]* train01p)) + lmd1[1:length(K1)]*w1[1:length(K1)]**2/2 + 2*K1*b1* train01/(eps0*seq(1,length(K1)))> L1
#            	cond2 <- log(1 + exp(-train02* w2[1:length(K2)]* train02p)) + lmd2[1:length(K2)]*w2[1:length(K2)]**2/2 + 2*K2*b2* train02/(eps0*seq(1,length(K2)))> L2
#			if (sum(is.nan(cond1))!= length(cond1) & sum(cond1==0)!= length(cond1)){
#                		w1 <- (-train01)* train01p/(1 + exp(-train01 * w1[1:length(K1)]* train01p))}
#			if (sum(is.nan(cond2))!= length(cond2) & sum(cond2==0)!= length(cond2)){
#                		w2 <- (-train02)* train02p/(1 + exp(-train02 * w2[1:length(K2)]* train02p))}
#            	if (sum(is.nan(cond1))== length(cond1) | sum(cond1==0)== length(cond1) | (sum(is.nan(cond2))== length(cond2) | sum(cond2==0)== length(cond2))){
#   	            	k <- k+1
#                		steps <- 50
#                		break
#           		}
#			W1[snpstxbp1.motif[i]+train01id] <- w1 
#			W2[snpstxbp1.motif[i]+train02id] <- w2
#			L1 <- log(1 + exp(train01 * w1[1:length(K1)]* train01p)) + lmd1[1:length(K1)]*w1[1:length(K1)]**2/2 + 2*K1*b1[1:length(K1)]* train01p/(eps0*seq(1,length(K1)))
#            	L2 <- log(1 + exp(train02 * w2[1:length(K2)]* train02p)) + lmd2[1:length(K2)]*w2[1:length(K2)]**2/2 + 2*K2*b2[1:length(K2)]* train02p/(eps0*seq(1,length(K2)))
#			steps <- steps -1
#                  if (steps == 0){
#                		l_temp1[is.nan(l_temp1[(1+length(K1)*(k-1)):(length(K1)*k)])] <- 0.5*l_temp1[max(is.nan(l_temp1[1+length(K1)*(k-1):length(K1)*k]-1),1)]+0.5*1_temp1[is.nan(1_temp1[min(1_temp1+length(K1)*(k-1):length(K1)*k+1,n)])]
#                		l_temp1[sum(is.nan(l[1+length(K)*(k-1):length(K)*k])==length(l),k] <- 0.5*l[sum(is.nan(l[1+length(K)*(k-1):length(K)*k])]==length(l),max(k-1,1))+0.5*l[sum(is.nan(l[1+length(K)*(k-1):length(K)*k]))]==length[l],min(k+1,n))
#                		l_temp2[is.nan(l_temp2[(1+length(K2)*(k-1)):(length(K2)*k)])] <- 0.5*l_temp2[max(is.nan(l_temp2[1+length(K2)*(k-1):length(K2)*k]-1),1)]+0.5*l_temp2[is.nan(l_temp2[min(1_temp2+length(K2)*(k-1):length(K2)*k+1,n)])]
#	            	ypred[l_temp1>0.5] <- 0
#				ypred[l_temp1<=0.5] <- -1
#				ypred[l_temp2>0.5] <- 2
#				ypred[l_temp2<=0.5] <- 1
#				rep(0, 6*(snpstxbp1.motif[i+1]-snpstxbp1.motif[i]-6))
#				k <- k+1
#                	      steps <- 50
#                        break
#				}             
#           if k > length(ctx)
#                 break
#	      	if k > snpstxbp1.motif[i+1]   break
#	}
		
#}


#dependent, numeric
for 



#TGCAAG
for (i in seq(1,min(N100, length(snpstxbp1.motif)))){
	# chr	
	train0<- str_sub(seq1, snpstxbp1.motif[i], snpstxbp1.motif[i] + snpstxbp1.len[i]-1)
	initHMM<-()
	str_sub(seq1, snpstxbp1.motif[i], snpstxbp1.motif[i] + snpstxbp1.len[i]-1)
	c<- sum(kse-(kse-mean(kse))/sqrt(var(kse)))
	w<- rep(1/length(train0), length(train0))
	alpha<- sum(w*kse)
	trainT<- str_sub(seq1, snpstxbp1.motif[[1]][i+1], snpstxbp1.motif[[1]][i+1] + snpstxbp1.len[i+1]-1)
	beta<- sum(w*kse) 
	Kse<- alpha*beta/
}

#tutorial: 
#Genomic state annotation of Roadmap Epigenomics
#Sequencing data

## Loading library and data


data(trainRegions)
str(trainRegions)

dim(pilot.hg19)
data(pilot.hg19)

celltypes = list("E123"=grep("E123", names(trainRegions)),
"E116"=grep("E116", names(trainRegions)))

sizeFactors <- getSizeFactors(trainRegions, celltypes)

nStates = 10
hmm_poilog <- initHMM(trainRegions, nStates, "PoissonLogNormal", sizeFactors)
hmm_fitted_poilog <- fitHMM(trainRegions, hmm_poilog, sizeFactors=sizeFactors, maxIters=10)
viterbi_poilog <- getViterbi(hmm_fitted_poilog, trainRegions, sizeFactors)
viterbi_poilog_gm12878 <- viterbi2GRanges(viterbi_poilog[1:3], regions=pilot.hg19, binSize=200)


hmm_nb <- initHMM(trainRegions, nStates, "NegativeBinomial", sizeFactors)
hmm_fitted_nb <- fitHMM(trainRegions, hmm_nb, sizeFactors=sizeFactors, maxIters=10)
viterbi_nb <- getViterbi(hmm_fitted_nb, trainRegions, sizeFactors=sizeFactors)
viterbi_nb_gm12878 <- viterbi2GRanges(viterbi_nb[1:3], pilot.hg19, 200)


avg_cov_nb <- getAvgSignal(viterbi_nb, trainRegions)
avg_cov_poilog <- getAvgSignal(viterbi_poilog, trainRegions)

library(gplots)
heat <- c("dark blue", "dodgerblue4", "darkred", "red", "orange", "gold", "yellow")
colfct <- colorRampPalette(heat)
colpal_statemeans <- colfct(200)

## define state order and colors
ord_nb <- order(apply(avg_cov_nb,1,max), decreasing=TRUE)
statecols_nb <- rainbow(nStates)
names(statecols_nb) <- ord_nb
heatmap.2(log(avg_cov_nb+1)[as.character(ord_nb),], margins=c(8,7), srtCol=45,
RowSideColors<- statecols_nb[as.character(ord_nb)], dendrogram="none",
	Rowv=FALSE, Colv=FALSE, col=colpal_statemeans, trace="none",
	cellnote=round(avg_cov_nb,1)[as.character(ord_nb),], notecol="black")
## define state order and colors
ord_poilog <- order(apply(avg_cov_poilog,1,max), decreasing=TRUE)
statecols_poilog <- rainbow(nStates)
names(statecols_poilog) <- ord_poilog
heatmap.2(log(avg_cov_poilog+1)[as.character(ord_poilog),], margins=c(8,7), srtCol=45,
	RowSideColors=statecols_poilog[as.character(ord_poilog)], dendrogram="none",
	Rowv=FALSE, Colv=FALSE, col=colpal_statemeans, trace="none",
	cellnote=round(avg_cov_poilog,1)[as.character(ord_poilog),], notecol="black")

library(Gviz)
from <- start(pilot.hg19)[3]
to <- from+300000
gvizViterbi_nb <- viterbi2Gviz(viterbi_nb_gm12878, "chr11", "hg19", from, to, statecols_nb)
gvizViterbi_poilog <- viterbi2Gviz(viterbi_poilog_gm12878, "chr11", "hg19", from, to,
statecols_poilog)
gvizData <- data2Gviz(trainRegions[[3]], pilot.hg19[3], binSize = 200, gen = "hg19", col="black", chrom = "chr11")
#Then, we use the plotTracks function to plot everything (see Figure 2).
gaxis <- GenomeAxisTrack()
data(ucscGenes)
mySize <- c(1,rep(1.2,9), 0.5,0.5,3)
plotTracks(c(list(gaxis), gvizData,gvizViterbi_nb,gvizViterbi_poilog,ucscGenes["chr11"]),
	from=from, to=to, showFeatureId=FALSE, featureAnnotation="id", fontcolor.feature="black",
	cex.feature=0.7, background.title="darkgrey", lwd=2, sizes=mySize)

## Modeling Sequencing data using other emission functions

hmm_pois <- initHMM(trainRegions, nStates, "Poisson")
hmm_fitted_pois <- fitHMM(trainRegions, hmm_pois, maxIters=10)
viterbi_pois <- getViterbi(hmm_fitted_pois, trainRegions)

simData_nmn <- lapply(trainRegions, function(x) cbind(apply(x,1,sum), x))
hmm_nmn <- initHMM(simData_nmn, nStates, "NegativeMultinomial")
hmm_fitted_nmn <- fitHMM(simData_nmn, hmm_nmn, maxIters=10)
viterbi_nmn <- getViterbi(hmm_fitted_nmn, simData_nmn)

trainRegions_smooth <- lapply(trainRegions, function(x)
apply(log(x+sqrt(x^2+1)), 2, runningMean, 2))
hmm_gauss <- initHMM(trainRegions_smooth, nStates, "IndependentGaussian", sharedCov=TRUE)
hmm_fitted_gauss <- fitHMM(trainRegions_smooth, hmm_gauss, maxIters=10)
viterbi_gauss <- getViterbi(hmm_fitted_gauss, trainRegions_smooth)

trainRegions_binary = binarizeData(trainRegions)
hmm_ber = initHMM(trainRegions_binary, nStates, "Bernoulli")
hmm_fitted_ber = fitHMM(trainRegions_binary, hmm_ber, maxIters=10)
viterbi_ber = getViterbi(hmm_fitted_ber, trainRegions_binary)

#We calculate the mean read coverage for each method and segmentation:
avg_cov_gauss = getAvgSignal(viterbi_gauss, trainRegions)
avg_cov_nmn = getAvgSignal(viterbi_nmn, trainRegions)
avg_cov_ber = getAvgSignal(viterbi_ber, trainRegions)
avg_cov_pois = getAvgSignal(viterbi_pois, trainRegions)

#These are again plotted using the heatmap.2 function (see Figure 3).
heatmap.2(log(avg_cov_gauss+1), margins=c(8,7),srtCol=45, dendrogram="row", Rowv=TRUE,
	Colv=FALSE, col=colpal_statemeans, trace="none", notecex=0.7, cexRow=0.75, cexCol=1,
	cellnote=round(avg_cov_gauss,1), notecol="black")
heatmap.2(log(avg_cov_nmn+1), margins=c(8,7),srtCol=45, dendrogram="row", Rowv=TRUE,
	Colv=FALSE, col=colpal_statemeans, trace="none", notecex=0.7, cexRow=0.75, cexCol=1,
	cellnote=round(avg_cov_nmn,1), notecol="black")
heatmap.2(log(avg_cov_ber+1), margins=c(8,7),srtCol=45, dendrogram="row", Rowv=TRUE,
	Colv=FALSE, col=colpal_statemeans, trace="none", notecex=0.7, cexRow=0.75, cexCol=1,
	cellnote=round(avg_cov_ber,1), notecol="black")
heatmap.2(log(avg_cov_pois+1), margins=c(8,7),srtCol=45, dendrogram="row", Rowv=TRUE,
	Colv=FALSE, col=colpal_statemeans, trace="none", notecex=0.7, cexRow=0.75, cexCol=1,
	cellnote=round(avg_cov_pois,1), notecol="black")



bdhmm_gauss <- initBdHMM(yeastTF_databychrom_ex, dStates = dStates, method = "Gaussian", directedObs=dirobs)
bdhmm_fitted_gauss <- fitHMM(yeastTF_databychrom_ex, bdhmm_gauss)
viterbi_bdhmm_gauss <- getViterbi(bdhmm_fitted_gauss, yeastTF_databychrom_ex)


