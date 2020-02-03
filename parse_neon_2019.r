#Version 1 data
ALL_mono_data <- read.csv("/home/liu/EDGE/monoallele_2/ALL_monoalleles_data_f5.csv",header=F,stringsAsFactors=F)
ALL_mono_data_A11 <- read.csv("/home/liu/EDGE/A1101_MPTAC_rmdup_blast_flanking_neonTPM.csv",header=T,stringsAsFactors=F)
ALL_mono_data_A2  <- read.csv("/home/liu/EDGE/A0201_MPTAC_rmdup_blast_flanking_neonTPM.csv",header=T,stringsAsFactors=F)

#includ MAPTAC monoallele
MPTAC_A11 <- read.csv("/home/liu/EDGE/A1101_MPTAC.csv",header=T,stringsAsFactors=F)
MPTAC_A2  <- read.csv("/home/liu/EDGE/A0201_MPTAC.csv",header=T,stringsAsFactors=F)
colnames(ALL_mono_data_A11) <- colnames(ALL_mono_data)
colnames(ALL_mono_data_A2) <- colnames(ALL_mono_data)
ALL_mono_data <- rbind(ALL_mono_data,ALL_mono_data_A11,ALL_mono_data_A2)
#load deconvlution peptidome
ALL_multi_data <- read.csv("/home/liu/EDGE/monoallele_2/deconv_CCLETPM_f5.csv",stringsAsFactors=F,header=T)
colnames(ALL_multi_data) <- colnames(ALL_mono_data)
ALL_multi_data <- ALL_multi_data[-which(ALL_multi_data[,5] %in% c("HLA-C*07:02","HLA-C*01:02","HLA-C*03:04","HLA-C*03:03")),]
tmp.negcon <- read.csv("/home/liu/EDGE/monoallele_9/negcon.csv",header=F,stringsAsFactors=F)
ALL_mono_data <- ALL_mono_data[!(ALL_mono_data[,1] %in% tmp.negcon[,1]),]
ALL_multi_data <- ALL_multi_data[!(ALL_multi_data[,1] %in% tmp.negcon[,1]),]
#moving deconv overlaping with MAPTAC
index.mv2mono <- which((ALL_multi_data[,1] %in% MPTAC_A11[,1]) & (ALL_multi_data[,5]=="HLA-A*11:01"))
ALL_mono_data <- rbind(ALL_mono_data,ALL_multi_data[index.mv2mono,])
ALL_multi_data <- ALL_multi_data[-index.mv2mono,]
index.mv2mono <- which((ALL_multi_data[,1] %in% MPTAC_A2[,1]) & (ALL_multi_data[,5]=="HLA-A*02:01"))
ALL_mono_data <- rbind(ALL_mono_data,ALL_multi_data[index.mv2mono,])
ALL_multi_data <- ALL_multi_data[-index.mv2mono,]

mono2019_unfiltered <- read.csv("/home/liu/EDGE/neon2019/mono_unfiltered.csv",header=T,stringsAsFactors=F)
mono2019_unfiltered <- mono2019_unfiltered[,c(1,5,7,8)]
colnames(mono2019_unfiltered) <- c("Allele","Peptide","geneSymbol","accession_numbers")
#blasting
mono2019_data <- read.csv("/home/liu/EDGE/mono_filtered_2019.csv",header=T,stringsAsFactors=F)
outname="/home/liu/EDGE/neon2019/neon2019mono_31.csv"
TMP.name="/tmp/neon2019mono_31.fsa"
RANGE=1:186464
print(RANGE)
print(TMP.name)
print(outname)
ncbi_blastp_short_peptide <- function(query,tmp.name,MATRIX="PAM30",gap=9,ext=1, DB="/home/liu/ncbi/public/refseq/uniport"){
    write.table(rbind(">",query),tmp.name,col.names=FALSE,row.names=FALSE,quote=FALSE)
    tmp <- system(paste("blastp -db ",DB," -query ",tmp.name," -evalue 200000 -max_target_seqs 5000 -word_size 2 -matrix ",MATRIX," -gapopen ",gap," -gapextend ",ext," -outfmt \"10 qseqid sseqid sseq pident qlen length mismatch gapope evalue bitscore sstart send\" "),intern = TRUE)
    tmp2 <- t(as.data.frame(lapply(tmp,strsplit,",")))
    PM.index <- which(tmp2[,3]==query)    
	if(length(PM.index)==0){
	    return(c(NA,NA,NA))
	}else{
	    return(tmp2[PM.index[1],c(2,10,11)])
	}
}
res.rmdup <- mono2019_data[RANGE,] #186464
blast_result <- matrix(ncol=3,nrow=nrow(res.rmdup))
for(i in 1:nrow(res.rmdup)){
blast_result[i,] <- ncbi_blastp_short_peptide(res.rmdup[i,3],tmp.name=TMP.name)
print(i)
}
res.rmdup_blast_res<-data.frame(res.rmdup,"info"=blast_result[,1],"start_pos"=as.numeric(blast_result[,2]),"end_pos"=as.numeric(blast_result[,3]),stringsAsFactors=FALSE)
write.table(res.rmdup_blast_res,outname,col.names=T,row.names=F,quote=F,sep=",")

for(i in 1:38){
tmp.name <- paste(c("/home/liu/EDGE/neon2019/neon2019_",i,".csv"),collapse="")
tmp.res <- read.csv(tmp.name,header=TRUE,stringsAsFactors=FALSE)
if(i==1){res.rmdup_blast_res <- tmp.res}else{res.rmdup_blast_res <- rbind(res.rmdup_blast_res,tmp.res)}
}
res.rmdup_blast_res <- res.rmdup_blast_res[!duplicated(res.rmdup_blast_res),]
res.rmdup_blast_res.tmp <- merge(res.rmdup_blast_res,mono2019_unfiltered,all.x=TRUE)
res.rmdup_blast_res.unIDable <- res.rmdup_blast_res.tmp[is.na(res.rmdup_blast_res.tmp$info),]
write.table(res.rmdup_blast_res.unIDable,"/home/liu/EDGE/neon2019/res.rmdup_blast_res.unIDable.csv",row.names=F,quote=F,sep=",")
res.rmdup_blast_res.IDable <- res.rmdup_blast_res.tmp[!is.na(res.rmdup_blast_res.tmp$info),]
write.table(res.rmdup_blast_res.IDable,"/home/liu/EDGE/neon2019/res.rmdup_blast_res.IDable.csv",row.names=F,quote=F,sep=",")
res.rmdup_blast_res2 <- res.rmdup_blast_res.IDable[which(nchar(res.rmdup_blast_res.IDable$Peptide)== res.rmdup_blast_res.IDable$end_pos - res.rmdup_blast_res.IDable$start_pos + 1),]
write.table(res.rmdup_blast_res2,"/home/liu/EDGE/neon2019/new_mono_rmdup_blast.csv",col.names=T,row.names=F,quote=F,sep=",")

#get flanking
library(seqinr)
UNIPORT_FASTA <- read.fasta(file = "/home/liu/uniport/uniprot_reviewed_human_UBremoved.fasta", 
  seqtype = c("AA"), as.string = FALSE, forceDNAtolower = FALSE,
  set.attributes = FALSE, legacy.mode = FALSE, seqonly = FALSE, strip.desc = FALSE,
  bfa = FALSE, sizeof.longlong = .Machine$sizeof.longlong,
  endian = .Platform$endian, apply.mask = TRUE)
index.uniport.fasta <- names(UNIPORT_FASTA)

Neon_2019_blast_res <- read.csv("/home/liu/EDGE/neon2019/new_mono_rmdup_blast.csv",header=T,stringsAsFactors=F)
for(i in 1:nrow(Neon_2019_blast_res)){
    print(i)
    rm("protein_seq")
    index.uni <- which(index.uniport.fasta==Neon_2019_blast_res$info[i])
	if(length(index.uni)>0){
    protein_seq <- UNIPORT_FASTA[index.uni][[1]]
    U_start <- (Neon_2019_blast_res$start_pos[i]-5)
    U_end <- (Neon_2019_blast_res$start_pos[i]-1)
    D_start <- (Neon_2019_blast_res$end_pos[i]+1)
    D_end <- (Neon_2019_blast_res$end_pos[i]+5)
    if(U_start<1){U_start<-0}
    if(U_end<1){U_end<-0}
    if(D_start>length(protein_seq)){D_start<-0}
    if(D_end>length(protein_seq)){D_end<-length(protein_seq)}
	U_seq <- protein_seq[U_start:U_end]
	D_seq <- protein_seq[D_start:D_end]
	Neon_2019_blast_res$U5_seq[i] <- paste(U_seq,collapse="")
	Neon_2019_blast_res$D5_seq[i] <- paste(D_seq,collapse="")
	}
}
#compensate margin blank
for(i in 1:nrow(Neon_2019_blast_res)){
print(i)
if(nchar(Neon_2019_blast_res$U5_seq[i])<5){
    short <- nchar(Neon_2019_blast_res$U5_seq[i])
    Neon_2019_blast_res$U5_seq[i] <- paste(c(rep("-",5-short),Neon_2019_blast_res$U5_seq[i]),collapse="")
}
if(nchar(Neon_2019_blast_res$D5_seq[i])<5){
    short <- nchar(Neon_2019_blast_res$D5_seq[i])
    Neon_2019_blast_res$D5_seq[i] <- paste(c(Neon_2019_blast_res$D5_seq[i],rep("-",5-short)),collapse="")
}
}

f5 <- gsub(" ","",paste(Neon_2019_blast_res$U5_seq,Neon_2019_blast_res$D5_seq))
Neon_2019_blast_res$f5 <- f5
Neon_2019_blast_res2 <- Neon_2019_blast_res[nchar(f5)==10,]
write.table(Neon_2019_blast_res2,"/home/liu/EDGE/neon2019/new_mono_rmdup_blast_flanking.csv",col.names=F,row.names=F,quote=F,sep=",")

#assign TPM
neon_expression <- read.csv("/home/liu/EDGE/neon2019/neon_expression.csv",header=T,stringsAsFactors=F)
tmp.multi <- read.csv("/home/liu/EDGE/neon2019/new_mono_rmdup_blast_flanking.csv",header=F,stringsAsFactors=F)
Neon_2019_blast_res2 <- tmp.multi[,c(1,2,3,4,7,8,11)]
#Neon_2019_blast_res2[,4] <- unlist(lapply(Neon_2019_blast_res2[,4],function(x){unlist(strsplit(x,"\\|"))[2]}))
for(i in 1:nrow(Neon_2019_blast_res2)){
    geneid <- Neon_2019_blast_res2[i,6]
    geneids <- unlist(strsplit(geneid,"\\|"))
    geneids.index <- which(neon_expression[,1] %in% geneids)
	if(length(geneids.index)!=0){
	    Neon_2019_blast_res2$exprs[i] <- max(neon_expression[geneids.index,4])
	} else {
	    Neon_2019_blast_res2$exprs[i] <- 0.00001
	}
}
write.table(Neon_2019_blast_res2,"/home/liu/EDGE/neon2019/new_mono_rmdup_blast_flanking_TPM.csv",col.names=F,row.names=F,quote=F,sep=",")

Neon_2019_blast_res2 <- read.csv("/home/liu/EDGE/neon2019/new_mono_rmdup_blast_flanking_TPM.csv",header=FALSE, stringsAsFactors=FALSE)
Neon_2019_blast_res2 <- Neon_2019_blast_res2[Neon_2019_blast_res2[,3] %in% 8:11,]

for(i in 1:length(hla)){
    print(hla[i])
    tmp.hla <- gsub(":","",(gsub("\\*","",gsub("HLA-","",hla[i]))))
    tmp.pep.2019 <- Neon_2019_blast_res2[Neon_2019_blast_res2[,1]==tmp.hla,]
    tmp.pep.2019.tmp <- tmp.pep.2019[,c(2,4,7,8,1)]
    tmp.pep.2019.tmp[,4] <- log10(tmp.pep.2019.tmp[,4])
    tmp.pep.2019.tmp[,6] <- "mono"
    tmp.pep.2019.tmp[,5] <- hla[i]
	colnames(tmp.pep.2019.tmp) <- colnames(ALL_mono_data)
    tmp.pep.acer.mono <- ALL_mono_data[ALL_mono_data[,5]==hla[i],]
    tmp.pep.acer.mono <- tmp.pep.acer.mono[!duplicated(tmp.pep.acer.mono),]
    tmp.pep.acer.mult <- ALL_multi_data[ALL_multi_data[,5]==hla[i],]
    tmp.pep.acer.mult <- tmp.pep.acer.mult[!duplicated(tmp.pep.acer.mult),]
    
    #if multi matched, move 2019 to mono
    index.multi <- which(tmp.pep.acer.mult[,1] %in% tmp.pep.2019.tmp[,1])
    index.2019 <- which(tmp.pep.2019.tmp[,1] %in% tmp.pep.acer.mult[,1])
    if(length(index.multi)!=0){
        tmp.pep.acer.mult <- tmp.pep.acer.mult[-index.multi,] #very important, do not omit it !
        tmp.pep.acer.mono_2019 <- rbind(tmp.pep.acer.mono,tmp.pep.2019.tmp[index.2019,])
        tmp.pep.acer.mono_2019 <- tmp.pep.acer.mono_2019[!duplicated(tmp.pep.acer.mono_2019),]
        #sum(tmp.pep.acer.mult[,1] %in% tmp.pep.acer.mono_2019[,1])
    } else {
	    tmp.pep.acer.mono_2019 <- tmp.pep.acer.mono
	}
    #if mono_2017 matched, take maximum expression record
	
    total.mono <- rbind(tmp.pep.acer.mono_2019,tmp.pep.2019.tmp)
    unique.pep <- unique(total.mono[,1])
	res.mono <- data.frame(matrix(ncol=6,nrow=length(unique.pep)))
    for(j in 1:length(unique.pep)){
	    index.pep <- which(total.mono[,1]==unique.pep[j])
		if(length(index.pep)==1){
		    res.mono[j,] <- total.mono[index.pep,]
		} else {
		    highest.expression <- max(total.mono[index.pep,4])
			tmp.res <- total.mono[index.pep[total.mono[index.pep,4] == highest.expression],]
			res.mono[j,] <- tmp.res[1,]
	    }
	}
	if(nrow(tmp.pep.acer.mult)>0){
	    colnames(tmp.pep.acer.mult) <- colnames(res.mono) 
	    res.allele <- rbind(res.mono,tmp.pep.acer.mult)
	} else {
	    res.allele <- res.mono
	}
	if(i==1){res.all <- res.allele} else {res.all <- rbind(res.all,res.allele)}
}

#put all mono allele into it
all.hla <- unique(Neon_2019_blast_res2[,1])
tmp.hla <- unlist(lapply(hla,function(x)gsub(":","",(gsub("\\*","",gsub("HLA-","",x))))))
for(i in 1:length(all.hla)){
    if(all.hla[i] %in% tmp.hla){
	    print(all.hla[i])
	    print("done")
	} else {
	    print(all.hla[i])
        tmp.pep.2019 <- Neon_2019_blast_res2[Neon_2019_blast_res2[,1]==all.hla[i],]
        tmp.pep.2019.tmp <- tmp.pep.2019[,c(2,4,7,8,1)]
        tmp.pep.2019.tmp[,4] <- log10(tmp.pep.2019.tmp[,4])
        tmp.pep.2019.tmp[,6] <- "mono"
		tmp.hla.sep <- unlist(strsplit(all.hla[i],""))
        tmp.pep.2019.tmp[,5] <- paste(c("HLA-",tmp.hla.sep[1],"*",tmp.hla.sep[2:3],":",tmp.hla.sep[4:5]),collapse="")
		colnames(tmp.pep.2019.tmp) <- colnames(res.all)
		res.all <- rbind(res.all,tmp.pep.2019.tmp)
	}
}
write.table(res.all,"/home/liu/EDGE/neon2019/positive_pep.csv",row.names=F,quote=F,sep=",")
res.all <- read.csv("/home/liu/EDGE/neon2019/positive_pep.csv",header=T,stringsAsFactors=F)
