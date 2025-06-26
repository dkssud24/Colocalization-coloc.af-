#!/bin/bash

# 사용법: ./run_plink_clump.sh 유전자이름
GENE=$1

# 1. GENCODE GFF3에서 유전자 정보 추출 (chr, start, end, strand, gene_name, gene_id)
INFO=$(awk '$3=="gene" && $0 ~ /gene_type "protein_coding"/ && $0 ~ /gene_name "BPTF"/ {print $1, $4, $5, $7, $9}' gencode.v19.annotation.gtf)
CHR=$(echo $INFO | awk '{print $1}' | sed 's/chr//')
START=$(echo $INFO | awk '{print $2}')
END=$(echo $INFO | awk '{print $3}')
STRAND=$(echo $INFO | awk '{print $4}')
GENE_ID=$(echo $INFO | awk -F 'gene_id "' '{print $2}' | awk -F '"' '{print $1}')
GENE_NAME=$(echo $INFO | awk -F 'gene_name "' '{print $2}' | awk -F '"' '{print $1}')


# chr start end strand gene_name gene_id
CHR=$(echo $INFO | awk '{print $1}' | sed 's/chr//')
START=$(echo $INFO | awk '{print $2}')
END=$(echo $INFO | awk '{print $3}')

# 2. 500kb window 계산
FROM_BP=$((START - 500000))
TO_BP=$((END + 500000))

# 3. PLINK extract 및 clump 수행
# (필요한 경우 FROM_BP, TO_BP 값이 0 미만이거나 음수가 되는 유전자 예외 처리 가능)
OUT=${GENE}

plink --bfile /BiO/hae/000006_ref_1000G/ref \
      --chr $CHR \
      --from-bp $FROM_BP \
      --to-bp $TO_BP \
      --make-bed --out ${OUT}

plink --bfile ${OUT} \
      --clump /BiO/hae/000029_GenomicSEM/central_MAF005/Result.GWAS.input.fuma \
      --clump-p1 1 --clump-p2 1 --clump-kb 1000 --clump-r2 0.01 \
      --out ${OUT}_clump

# (필요시 intermediate/최종 output 경로나 정리 추가 가능)
