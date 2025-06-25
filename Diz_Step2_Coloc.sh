library(data.table)
library(dplyr)
library(arrow)
library(coloc)
library(data.table)
library(stringr)

arg1_tissue <- "Adipose_Visceral_Omentum"
arg4_diz_sub <- "central"
gene_name <- "ASCC2"


gff_file <- "gencode.v47.annotation.gff3"

# 주석줄 아닌 부분만 추출
lines <- readLines(gff_file)
lines <- lines[!grepl("^#", lines)]  # 주석라인 제거

# 임시파일로 저장
temp_file <- tempfile(fileext = ".gff3")
writeLines(lines, temp_file)

# fread로 읽기
gff <- fread(temp_file, sep = "\t", header = FALSE, fill = TRUE, data.table = FALSE)
colnames(gff)[1:9] <- c("chr", "source", "type", "start", "end", "score", "strand", "phase", "attribute")

# 이제부터 기존 코드와 동일
target_row <- gff[
  gff$type == "gene" & 
  grepl("gene_type=protein_coding", gff$attribute) & 
  grepl(paste0("gene_name=", gene_name), gff$attribute),
]

if (nrow(target_row) == 0) stop("Gene not found!")

gene_id <- str_match(target_row$attribute, "gene_id=([^;]+)")[,2]
chr <- gsub("chr", "", target_row$chr)
start_bp <- as.numeric(target_row$start)
end_bp <- as.numeric(target_row$end)

cat(paste("chr:", chr, "start:", start_bp, "end:", end_bp, "gene_id:", gene_id, "\n"))

# 변수 자동 할당
arg2_chr <- chr
arg3_gene <- gene_name       
arg5_gene_ENS <- gene_id
arg6_start_BP <- start_bp
arg7_end_BP <- end_bp



var11 <- paste(arg3_gene,"_clump.clumped",sep="")
d <- fread(var11)
d <- data.frame(d$SNP)
names(d) <- c("variant_id")
var12 <- paste("/BiO/hae/000037_Predixcan/",arg4_diz_sub,".step1",sep="")
hg38 <- fread(var12)


m <- left_join(d,hg38,by="variant_id")
m2 <- na.omit(m)
var13 <- paste(arg3_gene,"_clump.clumped.hg38",sep="")
fwrite(m2,var13,quote=FALSE,row.names=FALSE,sep="\t")


var10 <- paste("/BiO/hae/000019_EUR_eQTL/GTEx_Analysis_v10_eQTL_all_associations/",arg1_tissue,".v10.allpairs.chr",arg2_chr,".parquet",sep="")
eqtl_data <- read_parquet(var10)
var14 <- paste("/BiO/hae/000037_Predixcan/",arg4_diz_sub,".step1",sep="")
gwas <- fread(var14)
eqtl_data2 <- eqtl_data[grepl(paste0("^", arg5_gene_ENS, "(\\.|$)"), eqtl_data$gene_id), ]
eqtl_data2 <- data.frame(eqtl_data2)
start <- arg6_start_BP
end <- arg7_end_BP
cis_start <- start - 500000
cis_end <- end + 500000
chr <- arg2_chr

names(eqtl_data2)[2] <- c("SNP")
names(m2)[2] <- c("SNP")

m3 <- left_join(m2,eqtl_data2,by="SNP")
m4 <- na.omit(m3)

m5 <- data.frame(m4$SNP)
names(m5) <- c("SNP")
haha <- dim(m5)[1]

coloc_summary <- data.frame(
  SNP = character(),
  nsnps = numeric(),
  H0 = numeric(),
  H1 = numeric(),
  H2 = numeric(),
  H3 = numeric(),
  H4 = numeric(),
  stringsAsFactors = FALSE
)


for( i in 1:haha){
zz <- m5[i,]
zz2 <- m4 [ m4$SNP == zz,]
aa2 <- zz2$position
real_start <- aa2 - 500000
real_end <- aa2 + 500000

m_tmp <- m4 [ m4$position > real_start,]
m6 <- m_tmp [ m_tmp$position < real_end,]

dataset1 <- list(snp=m6$SNP,beta=m6$slope,varbeta=(m6$slope_se)^2,pval=m6$pval_nominal,MAF = pmin(m6$af, 1-m6$af),N=15201,type="quant")
dataset2 <- list(snp=m6$SNP,beta=m6$effect_size,varbeta=(m6$standard_error)^2,pval=m6$pvalue,MAF=pmin(m6$af,1-m6$af),N=889004,type="quant")


my.res <- coloc.abf(dataset1=dataset1,dataset2=dataset2, p1=1e-4, p2=1e-4, p12=1e-5)
s <- my.res$summary

  coloc_summary <- rbind(coloc_summary, data.frame(
    SNP = as.character(selected_SNP$SNP),
    nsnps = s["nsnps"],
    H0 = s["PP.H0.abf"],
    H1 = s["PP.H1.abf"],
    H2 = s["PP.H2.abf"],
    H3 = s["PP.H3.abf"],
    H4 = s["PP.H4.abf"]
  ))
}
var20 <- paste("Result.",arg3_gene,".",arg1_tissue,sep="")
names(gwas)[2] <- c("SNP")
tmp2 <- left_join(coloc_summary,gwas,by="SNP")
fwrite(tmp2,var20,quote=FALSE,row.names=FALSE,sep="\t")
