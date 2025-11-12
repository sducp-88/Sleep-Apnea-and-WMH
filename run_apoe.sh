#!/usr/bin/env bash
set -euo pipefail

# ==================================================================
# APOE genotype extractor (v4, POS-based header matching)
# - UKB imputed v3 (BGEN 1.2, GRCh37) support
# - Exports PLINK2 .raw, then finds rs429358/rs7412 columns by POS
#   regardless of header ID decoration (:_ or _A1 or _<ALT> etc.)
# - NA-safe; bi-allelic; outputs APOE_calls.csv (headered only)
# ==================================================================

BGEN_IN="/home/dnanexus/in/in/0/ukb22828_c19_b0_v3.bgen"
BGI_IN="/home/dnanexus/in/in/1/ukb22828_c19_b0_v3.bgen.bgi"
SAMPLE_IN="/home/dnanexus/in/in/2/ukb22828_c19_b0_v3.sample"

cp -f "$BGEN_IN"   ./chr19.bgen
cp -f "$BGI_IN"    ./chr19.bgen.bgi
cp -f "$SAMPLE_IN" ./chr19.sample

# Windows
HG19_RANGE="19:45411000-45413000"
HG38_RANGE="19:44908000-44911000"

bgenix -g chr19.bgen -incl-range ${HG19_RANGE} > apoe_hg19.bgen || true
bgenix -g chr19.bgen -incl-range ${HG38_RANGE} > apoe_hg38.bgen || true

plink2 --bgen apoe_hg19.bgen ref-first --sample chr19.sample --max-alleles 2 \
  --make-pgen --set-all-var-ids @:#:\$r:\$a --new-id-max-allele-len 100 \
  --out apoe_hg19 >/dev/null 2>&1 || true

plink2 --bgen apoe_hg38.bgen ref-first --sample chr19.sample --max-alleles 2 \
  --make-pgen --set-all-var-ids @:#:\$r:\$a --new-id-max-allele-len 100 \
  --out apoe_hg38 >/dev/null 2>&1 || true

CNT19=0; CNT38=0
[[ -f apoe_hg19.pvar ]] && CNT19=$(grep -vc '^#' apoe_hg19.pvar || true)
[[ -f apoe_hg38.pvar ]] && CNT38=$(grep -vc '^#' apoe_hg38.pvar || true)

if [[ ${CNT19} -gt 0 ]]; then
  CHUNK="apoe_hg19"
elif [[ ${CNT38} -gt 0 ]]; then
  CHUNK="apoe_hg38"
else
  echo "ERROR: No variants in APOE windows." >&2
  exit 2
fi

cp -f ${CHUNK}.pgen apoe.pgen
cp -f ${CHUNK}.pvar apoe.pvar
cp -f ${CHUNK}.psam apoe.psam

plink2 --pfile apoe --export A --out apoe >/dev/null

# Build CSV using POS-based matching from .raw header only.
# We accept header tokens like:
#   19:45411941:T:C
#   19:45411941:T:C_A1
#   19:45411941:T:C_C
#   19_45411941_T_C
#   19_45411941_T_C_C
# Parsing rule: replace ":" and "_" with spaces, split ->
#   [chr pos ref alt (opt a1)]
awk -v OFS=',' '
function tokenize(tok, arr,    s){ s=tok; gsub("[:_]", " ", s); n=split(s, arr, /[ ]+/); return n }
function find_by_pos(h, pos_target,    i, tok, parts, n, pos){
  for(i=1;i<=NF;i++){
    tok=$i
    n=tokenize(tok, parts)
    if(n<4) continue
    pos=parts[2]+0
    if(pos==pos_target) return i
  }
  return 0
}
function a1_from_token(tok,    parts, n){
  n=tokenize(tok, parts)
  if(n>=5) return parts[5]
  # If no explicit A1 suffix, assume A1 = ALT
  if(n>=4) return parts[4]
  return ""
}
function allele_at(tok, idx,  parts, n){
  n=tokenize(tok, parts)
  if(idx<=n) return parts[idx]
  return ""
}

BEGIN{
  pos429_19=45411941; pos7412_19=45412079;
  pos429_38=44909813; pos7412_38=44908822;
  header_dumped=0;
}
# First file is .raw: header on FNR==1
FNR==1 && NR==1 {
  # store header line in H[1..NF]
  for(i=1;i<=NF;i++){ H[i]=$i }
  # try hg19 first
  c429=find_by_pos(H, pos429_19)
  c7412=find_by_pos(H, pos7412_19)
  build="hg19"
  if(!(c429>0 && c7412>0)){
    # try hg38
    c429=find_by_pos(H, pos429_38)
    c7412=find_by_pos(H, pos7412_38)
    build="hg38"
  }
  if(!(c429>0 && c7412>0)){
    # dump the first 40 tokens for debugging
    msg="ERROR: Cannot find header columns for rs429358/rs7412 by POS. Header preview: "
    for(i=1;i<=40 && i<=NF;i++) msg=msg H[i] "|"
    print msg > "/dev/stderr"
    exit 4
  }
  # Derive A1 for each locus from the token itself
  a1_429=a1_from_token(H[c429])
  a1_7412=a1_from_token(H[c7412])
  # Also record normalized IDs (chr:pos:ref:alt) for output
  # Use ref/alt from token
  ref429=allele_at(H[c429],3); alt429=allele_at(H[c429],4)
  ref7412=allele_at(H[c7412],3); alt7412=allele_at(H[c7412],4)
  id429 = allele_at(H[c429],1) ":" allele_at(H[c429],2) ":" ref429 ":" alt429
  id7412= allele_at(H[c7412],1) ":" allele_at(H[c7412],2) ":" ref7412 ":" alt7412

  print "FID","IID","rs429358_id","rs7412_id","e4_count","e2_count","APOE_genotype","APOE_e4_carrier","APOE_e2_carrier"
  next
}

# Data rows (skip first 6 fixed columns if present; but indexing by c429/c7412 is safe)
NR>1 {
  fid=$1; iid=$2;
  d429s=$(c429); d7412s=$(c7412);
  e4=""; e2=""; geno="NA"; e4car=""; e2car="";
  if(d429s!="NA" && d7412s!="NA"){
    d429=d429s+0; d7412=d7412s+0;
    # dosage in .raw is for A1; count e4 if A1 is C at rs429358; else 2-dosage
    e4 = (a1_429=="C") ? d429 : (2 - d429);
    e2 = (a1_7412=="T") ? d7412 : (2 - d7412);
    if(e4<0.0005 && e4>-0.0005) e4=0; if(e4>1.9995 && e4<2.0005) e4=2;
    if(e2<0.0005 && e2>-0.0005) e2=0; if(e2>1.9995 && e2<2.0005) e2=2;
    if(e4==2) geno="ε4/ε4";
    else if(e2==2) geno="ε2/ε2";
    else if(e4==1 && e2==0) geno="ε3/ε4";
    else if(e4==0 && e2==1) geno="ε2/ε3";
    else if(e4==1 && e2==1) geno="ε2/ε4";
    else if(e4==0 && e2==0) geno="ε3/ε3";
    e4car=(e4>0)?1:0; e2car=(e2>0)?1:0;
  }
  print fid,iid,id429,id7412,e4,e2,geno,e4car,e2car
}
' apoe.raw > APOE_calls.csv

# Clean up
rm -f apoe.* apoe_hg19.* apoe_hg38.* chr19.bgen chr19.bgen.bgi chr19.sample

echo "APOE_calls.csv written."
