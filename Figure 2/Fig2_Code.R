setwd("~/Box Sync/Metagenomics_PRISM_BC")
#final script
#libraries
library(readstata13)
library(dplyr)
library(ggplot2)
library(tibble)
library(knitr)
library(tidyr)
library(tidyverse)
library(fs)
library(lubridate)
library(readr)
library(readxl)
library(patchwork)
library(pals)
library(plyr)
library(stringr)

#attaining metrics for fig 2A and 2B
#Final Run - Batch 1_Plasma
batch_1p<- fs::dir_ls("prism-3-_nextseq_1722/plasma/", regexp="\\.csv$")
#reading the files just by calling a function purr
batch_1p %>%
  map_dfr(read_csv, .id = "source") -> batch_1pk
#renaming the samples
batch_1pk$source <-gsub("prism-3-_nextseq_1722/plasma/","",as.character((batch_1pk$source)))
batch_1pk$source <-gsub(".csv","",as.character((batch_1pk$source)))

## Threshold filters
CUTOFF.NT_L <- 50
CUTOFF.NT_rPM <- 10
CUTOFF.NR_rPM <- 5
CUTOFF.NT_Z <- 1

#filtering only at the species level
batch_1pk %>%
  filter(tax_level==1) -> batch_1pk

#uncategorized reads - species level
batch_1pk%>%
  filter(is.na(category)) -> na_1pk
na_1pk$category=c('uncategorized')
#summing all the na by sample type
na_top_1pk <- aggregate(nt_count~source+category, data = na_1pk, FUN = sum)

#total reads
total_reads_per_sample_batch_1_plasma <- aggregate(nt_count ~ source, data = batch_1pk, FUN = sum)

#pre filtering
merged_pre_filter_1pk <- aggregate(nt_count ~ source + category, data = batch_1pk, FUN =sum)
#uncategorized reads + categorized reads no filtering step
merged_pre_filter_1pk_s <- rbind(na_top_1pk,merged_pre_filter_1pk)
#final merge with reads per sample
merged_pre_filter_1pk_reads <- merge(merged_pre_filter_1pk_s,total_reads_per_sample_batch_1_plasma, by ='source')
#reads per microbes
batch_1_plasma_microbes_unfiltered <- aggregate(nt_count ~ category, data = merged_pre_filter_1pk_s , FUN = sum)

#filtering according to threshold filters
batch_1pk %>%
  filter(category=="archaea") %>%
  filter(tax_level==1) %>%
  filter(is_phage=='FALSE') %>%
  filter(nt_alignment_length >= CUTOFF.NT_L &
           nt_rpm >= CUTOFF.NT_rPM &
           nr_rpm >= CUTOFF.NR_rPM &
           nt_z_score >= CUTOFF.NT_Z) -> arc1p
batch_1pk %>%
  filter(category=="bacteria") %>%
  filter(source %in% c("H1NKMS",
                       "H1QQVV",
                       "H1ZBSP",
                       "H1ZB32",
                       "H1ZMG8",
                       "H1Q3XE",
                       "H1ZJZH",
                       "H1Q7PY",
                       "H1RB6S",
                       "H1U8MB",
                       "H1PQS7",
                       "H1P7E6",
                       "H1RKSU",
                       "H1V9XG",
                       "H1TKPC",
                       "H1P33H",
                       "H1QU58",
                       "H1QYAZ",
                       "H1NQQ8",
                       "H1PRBT",
                       "H1PHE3",
                       "H1RFAG",
                       "H1NWL9",
                       "H1RKVJ",
                       "H1UYUH",
                       "H1Q778",
                       "H1NM3P",
                       "H1RS2K",
                       "H1ZHDX",
                       "H1P4X4",
                       "H1PRY2",
                       "H1PW7M",
                       "H1NZVG",
                       "H1RQA4",
                       "H1NLRG",
                       "H1NYPW",
                       "H1RKQN",
                       "H1RTDA",
                       "H1Q5BZ",
                       "H1QW7M",
                       "H1RHAN",
                       "H1QJF2",
                       "H1PYJN",
                       "H1PSWB",
                       "H1NUVL",
                       "H1QS8W",
                       "H1H7VV",
                       "H1NZ9J",
                       "H1QS6F",
                       "H1P9RY",
                       "H1QNSE",
                       "H1QZPN",
                       "H1V4WP",
                       "H1TL86",
                       "H1PLCT",
                       "H1PRRD",
                       "H1NMJL",
                       "H1QS2J",
                       "H1QQUY",
                       "H1377X",
                       "H1RQ8B",
                       "H1PEFS",
                       "H1RCZC",
                       "H1UT9Q",
                       "H1UD26",
                       "H1R4EP",
                       "H1QCEQ",
                       "H1P2VJ",
                       "H1P7XJ",
                       "H1R8MY",
                       "H1NSNV",
                       "H1ZLSW",
                       "H1PX9Z",
                       "H1RFLX",
                       "H1ZLHW",
                       "H1RB6Z",
                       "H1ZTS8",
                       "H1SZGP",
                       "H1QKHA",
                       "H1NPSL",
                       "H1ZJ4X",
                       "H1VAKB",
                       "H1NXDL",
                       "H1QCBJ",
                       "H1R5DZ",
                       "H1QCPT",
                       "H1PB96",
                       "H1RHT3",
                       "H1ZTNU",
                       "H1RQVU",
                       "H1RBCT",
                       "H1ZX7S",
                       "H1AELC",
                       "H1V492",
                       "H1ZH4F",
                       "H1NPHK",
                       "H1PURH",
                       "H1RP3H",
                       "H1RK3L",
                       "H1PQNW",
                       "H1ZWPT",
                       "H1HNUQ",
                       "H1ZBDV",
                       "H1PBK4",
                       "H1R7L7",
                       "H1UJ3U",
                       "H1STVV",
                       "H1UTXL",
                       "H1T2TN",
                       "H1PUQB",
                       "H1QNPA",
                       "H1ZVPX",
                       "H1Y48F",
                       "H1P4LP",
                       "H1U7YU",
                       "H1NP7W",
                       "H1U2WP",
                       "H1SMS4",
                       "H1QSNP",
                       "H1R2QD",
                       "H1DGZ4",
                       "H1T7JT",
                       "H1QPZ8",
                       "H1TTC9",
                       "H1Q6SH",
                       "H1VLDB",
                       "H1RACD",
                       "H1JJ9R",
                       "H1U9BT",
                       "H1QLG8",
                       "H1QMFH",
                       "H1NLMW",
                       "H1Q7PK",
                       "H1QPMU",
                       "H1TF98",
                       "H1ULAZ",
                       "H1NPT4",
                       "H1QNEQ",
                       "H1VMQ6",
                       "H1ZULV",
                       "H1QZLQ",
                       "H1PG4P",
                       "H1QJEG",
                       "H1NHEL",
                       "H1QJZN",
                       "H1NSNA",
                       "H1P5VK",
                       "H1KTNZ",
                       "H1RSQN",
                       "H1TCTD",
                       "H1ZKAT",
                       "H1FXVK",
                       "H1RSYZ",
                       "H1QRR7",
                       "H1R2UP",
                       "H1PXAQ",
                       "H1P8U7",
                       "H1NJAT",
                       "H1UH1P",
                       "H1RYNQ",
                       "H1VL2Q",
                       "H1T5QJ",
                       "H1Q896",
                       "H13YDQ",
                       "H1V9E6",
                       "H1P9BX",
                       "H1RTWP",
                       "H1EN34",
                       "H1UBSB",
                       "H1NS8T",
                       "H1VM26",
                       "H1ZVAN",
                       "H1V6DY",
                       "H1NRET",
                       "H1RZCP",
                       "H1Q4F6",
                       "H1STAL",
                       "H1R8X8",
                       "H1VK8W",
                       "H1PS8D",
                       "H1ZN52",
                       "H1RD6K",
                       "H1SBGN",
                       "H1P3HA",
                       "H1PSU4",
                       "H1PG2U",
                       "H1FYCC",
                       "H1RA4N",
                       "H1PQU8",
                       "H1NLQ3",
                       "H1QDYR",
                       "H19URT",
                       "H1RC3Y",
                       "H1287P",
                       "H1PKUW",
                       "H17U5Y",
                       "H1RLCR",
                       "H1V5VY",
                       "H1PUU9",
                       "H1QR4B")==FALSE)%>%
  filter(tax_level==1) %>%
  filter(is_phage=='FALSE') %>%
  filter(nt_alignment_length >= CUTOFF.NT_L &
           nt_rpm >= CUTOFF.NT_rPM &
           nr_rpm >= CUTOFF.NR_rPM &
           nt_z_score >= CUTOFF.NT_Z) -> bak1p
batch_1pk %>%
  filter(category=="eukaryota") %>%
  filter(tax_level==1) %>%
  filter(is_phage=='FALSE') %>%
  filter(nt_alignment_length >= CUTOFF.NT_L &
           nt_rpm >= CUTOFF.NT_rPM &
           nr_rpm >= CUTOFF.NR_rPM &
           nt_z_score >= CUTOFF.NT_Z) -> euk1p 
batch_1pk %>%
  filter(category=="viruses") %>%
  filter(tax_level==1) %>%
  filter(is_phage=='FALSE') %>%
  filter(nt_alignment_length >= CUTOFF.NT_L &
           nt_rpm >= CUTOFF.NT_rPM &
           nr_rpm >= CUTOFF.NR_rPM &
           nt_z_score >= CUTOFF.NT_Z) -> vik1p

#merging metrics - NT rPM 10
merged_1p <- rbind(arc1p,bak1p,euk1p,vik1p)%>%
  group_by(source,category,nt_count)%>%
  tally()%>%
  print(n=Inf) ->bk1p

#aggregating to get total number of reads per microbe type- NT rPM 10
p1s <- aggregate(nt_count ~ source + category, data = bk1p, FUN = sum)
p1sf <- merge(total_reads_per_sample_batch_1_plasma,p1s, by='source')
#merging filtered and uncategorized reads. 
p1k_uncat_filtered_reads <- rbind(na_top_1pk,p1s)
#merging the data with total reads
p1k_filtered_final_reads <- merge(p1k_uncat_filtered_reads,total_reads_per_sample_batch_1_plasma, by ='source')
#Figure 2 - filtered
batch_1_plasma_microbes_filtered <- aggregate(nt_count ~ category, data = p1k_uncat_filtered_reads, FUN = sum)

#attaining metrics for fig 2A and 2B
##### plasma batch 2
batch_2p<- fs::dir_ls("prism_3_batch_2_novaseq_2373/plasma/", regexp="\\.csv$")
#reading the files just by calling a function purr
batch_2p %>%
  map_dfr(read_csv, .id = "source") -> batch_2pk

batch_2pk$source <-gsub("prism_3_batch_2_novaseq_2373/plasma/","",as.character((batch_2pk$source)))
batch_2pk$source <-gsub(".csv","",as.character((batch_2pk$source)))

#filtering only at the species level
batch_2pk %>%
  filter(tax_level==1) -> batch_2pk

#uncategorized reads - species level
batch_2pk%>%
  filter(is.na(category)) -> na_2pk
na_2pk$category=c('uncategorized')
#summing all the na by sample type
na_top_2pk <- aggregate(nt_count~source+category, data = na_2pk, FUN = sum)
#total reads
total_reads_per_sample_batch_2_plasma <- aggregate(nt_count ~ source, data = batch_2pk, FUN = sum)

#pre filtering
merged_pre_filter_2pk <- aggregate(nt_count ~ source + category, data = batch_2pk, FUN =sum)
#uncategorized reads + categorized reads no filtering step
merged_pre_filter_2pk_s <- rbind(na_top_2pk,merged_pre_filter_2pk)
#final merge with reads per sample
merged_pre_filter_2pk_reads <- merge(merged_pre_filter_2pk_s,total_reads_per_sample_batch_2_plasma, by ='source')
#reads per microbes
batch_2_plasma_microbes_unfiltered <- aggregate(nt_count ~ category, data = merged_pre_filter_2pk_s , FUN = sum)

#archea
batch_2pk %>%
  filter(category=="archaea") %>%
  filter(tax_level==1) %>%
  filter(is_phage=='FALSE') %>%
  filter(nt_alignment_length >= CUTOFF.NT_L &
           nt_rpm >= CUTOFF.NT_rPM &
           nr_rpm >= CUTOFF.NR_rPM &
           nt_z_score >= CUTOFF.NT_Z) -> arc2p
batch_2pk %>%
  filter(category=="bacteria") %>%
  filter(source %in% c("H1NKMS",
                       "H1QQVV",
                       "H1ZBSP",
                       "H1ZB32",
                       "H1ZMG8",
                       "H1Q3XE",
                       "H1ZJZH",
                       "H1Q7PY",
                       "H1RB6S",
                       "H1U8MB",
                       "H1PQS7",
                       "H1P7E6",
                       "H1RKSU",
                       "H1V9XG",
                       "H1TKPC",
                       "H1P33H",
                       "H1QU58",
                       "H1QYAZ",
                       "H1NQQ8",
                       "H1PRBT",
                       "H1PHE3",
                       "H1RFAG",
                       "H1NWL9",
                       "H1RKVJ",
                       "H1UYUH",
                       "H1Q778",
                       "H1NM3P",
                       "H1RS2K",
                       "H1ZHDX",
                       "H1P4X4",
                       "H1PRY2",
                       "H1PW7M",
                       "H1NZVG",
                       "H1RQA4",
                       "H1NLRG",
                       "H1NYPW",
                       "H1RKQN",
                       "H1RTDA",
                       "H1Q5BZ",
                       "H1QW7M",
                       "H1RHAN",
                       "H1QJF2",
                       "H1PYJN",
                       "H1PSWB",
                       "H1NUVL",
                       "H1QS8W",
                       "H1H7VV",
                       "H1NZ9J",
                       "H1QS6F",
                       "H1P9RY",
                       "H1QNSE",
                       "H1QZPN",
                       "H1V4WP",
                       "H1TL86",
                       "H1PLCT",
                       "H1PRRD",
                       "H1NMJL",
                       "H1QS2J",
                       "H1QQUY",
                       "H1377X",
                       "H1RQ8B",
                       "H1PEFS",
                       "H1RCZC",
                       "H1UT9Q",
                       "H1UD26",
                       "H1R4EP",
                       "H1QCEQ",
                       "H1P2VJ",
                       "H1P7XJ",
                       "H1R8MY",
                       "H1NSNV",
                       "H1ZLSW",
                       "H1PX9Z",
                       "H1RFLX",
                       "H1ZLHW",
                       "H1RB6Z",
                       "H1ZTS8",
                       "H1SZGP",
                       "H1QKHA",
                       "H1NPSL",
                       "H1ZJ4X",
                       "H1VAKB",
                       "H1NXDL",
                       "H1QCBJ",
                       "H1R5DZ",
                       "H1QCPT",
                       "H1PB96",
                       "H1RHT3",
                       "H1ZTNU",
                       "H1RQVU",
                       "H1RBCT",
                       "H1ZX7S",
                       "H1AELC",
                       "H1V492",
                       "H1ZH4F",
                       "H1NPHK",
                       "H1PURH",
                       "H1RP3H",
                       "H1RK3L",
                       "H1PQNW",
                       "H1ZWPT",
                       "H1HNUQ",
                       "H1ZBDV",
                       "H1PBK4",
                       "H1R7L7",
                       "H1UJ3U",
                       "H1STVV",
                       "H1UTXL",
                       "H1T2TN",
                       "H1PUQB",
                       "H1QNPA",
                       "H1ZVPX",
                       "H1Y48F",
                       "H1P4LP",
                       "H1U7YU",
                       "H1NP7W",
                       "H1U2WP",
                       "H1SMS4",
                       "H1QSNP",
                       "H1R2QD",
                       "H1DGZ4",
                       "H1T7JT",
                       "H1QPZ8",
                       "H1TTC9",
                       "H1Q6SH",
                       "H1VLDB",
                       "H1RACD",
                       "H1JJ9R",
                       "H1U9BT",
                       "H1QLG8",
                       "H1QMFH",
                       "H1NLMW",
                       "H1Q7PK",
                       "H1QPMU",
                       "H1TF98",
                       "H1ULAZ",
                       "H1NPT4",
                       "H1QNEQ",
                       "H1VMQ6",
                       "H1ZULV",
                       "H1QZLQ",
                       "H1PG4P",
                       "H1QJEG",
                       "H1NHEL",
                       "H1QJZN",
                       "H1NSNA",
                       "H1P5VK",
                       "H1KTNZ",
                       "H1RSQN",
                       "H1TCTD",
                       "H1ZKAT",
                       "H1FXVK",
                       "H1RSYZ",
                       "H1QRR7",
                       "H1R2UP",
                       "H1PXAQ",
                       "H1P8U7",
                       "H1NJAT",
                       "H1UH1P",
                       "H1RYNQ",
                       "H1VL2Q",
                       "H1T5QJ",
                       "H1Q896",
                       "H13YDQ",
                       "H1V9E6",
                       "H1P9BX",
                       "H1RTWP",
                       "H1EN34",
                       "H1UBSB",
                       "H1NS8T",
                       "H1VM26",
                       "H1ZVAN",
                       "H1V6DY",
                       "H1NRET",
                       "H1RZCP",
                       "H1Q4F6",
                       "H1STAL",
                       "H1R8X8",
                       "H1VK8W",
                       "H1PS8D",
                       "H1ZN52",
                       "H1RD6K",
                       "H1SBGN",
                       "H1P3HA",
                       "H1PSU4",
                       "H1PG2U",
                       "H1FYCC",
                       "H1RA4N",
                       "H1PQU8",
                       "H1NLQ3",
                       "H1QDYR",
                       "H19URT",
                       "H1RC3Y",
                       "H1287P",
                       "H1PKUW",
                       "H17U5Y",
                       "H1RLCR",
                       "H1V5VY",
                       "H1PUU9",
                       "H1QR4B")==FALSE)%>%
  filter(tax_level==1) %>%
  filter(is_phage=='FALSE') %>%
  filter(nt_alignment_length >= CUTOFF.NT_L &
           nt_rpm >= CUTOFF.NT_rPM &
           nr_rpm >= CUTOFF.NR_rPM &
           nt_z_score >= CUTOFF.NT_Z) -> bak2p
batch_2pk %>%
  filter(category=="eukaryota") %>%
  filter(tax_level==1) %>%
  filter(is_phage=='FALSE') %>%
  filter(nt_alignment_length >= CUTOFF.NT_L &
           nt_rpm >= CUTOFF.NT_rPM &
           nr_rpm >= CUTOFF.NR_rPM &
           nt_z_score >= CUTOFF.NT_Z) -> euk2p 
batch_2pk %>%
  filter(category=="viruses") %>%
  filter(tax_level==1) %>%
  filter(is_phage=='FALSE') %>%
  filter(nt_alignment_length >= CUTOFF.NT_L &
           nt_rpm >= CUTOFF.NT_rPM &
           nr_rpm >= CUTOFF.NR_rPM &
           nt_z_score >= CUTOFF.NT_Z) -> vik2p
#merging metrics - nt_rpm_10
merged_2p <- rbind(arc2p,bak2p,euk2p,vik2p)%>%
  group_by(source,category,nt_count)%>%
  tally()%>%
  print(n=Inf) ->bk2p

merged_2p_ <- rbind(arc2p,bak2p,euk2p,vik2p) 
merged_1p_ <- rbind(arc1p,bak1p,euk1p,vik1p)                 
#aggregating to get total number of reads per microbe type - nt rpm 10
p2s <- aggregate(nt_count ~ source + category, data = bk2p, FUN = sum)
p2sf <- merge(total_reads_per_sample_batch_2_plasma,p2s, by='source')

#merging filtered and uncategorized reads. 
p2k_uncat_filtered_reads <- rbind(na_top_2pk,p2s)
#merging the data with total reads
p2k_filtered_final_reads <- merge(p2k_uncat_filtered_reads,total_reads_per_sample_batch_2_plasma, by ='source')
#Figure 2 - filtered
batch_2_plasma_microbes_filtered <- aggregate(nt_count ~ category, data = p2k_uncat_filtered_reads, FUN = sum)

#attaining metrics for fig 2C and 2D
#Swabs
#Final Run - Batch 1_swabs
batch_1sk<- fs::dir_ls("prism-3-_nextseq_1722/swab/", regexp="\\.csv$")
#reading the files just by calling a function purr
batch_1sk %>%
  map_dfr(read_csv, .id = "source") -> batch_1sk
#renaming the samples
batch_1sk$source <-gsub("prism-3-_nextseq_1722/swab/","",as.character((batch_1sk$source)))
batch_1sk$source <-gsub(".csv","",as.character((batch_1sk$source)))

#filtering only at the species level
batch_1sk %>%
  filter(tax_level==1) -> batch_1sk

#uncategorized reads - species level
batch_1sk%>%
  filter(is.na(category)) -> na_1sk
na_1sk$category=c('uncategorized')
#summing all the na by sample
na_top_1sk <- aggregate(nt_count~source+category, data = na_1sk, FUN = sum)
#total reads
total_reads_per_sample_batch_1_swab <- aggregate(nt_count ~ source, data = batch_1sk, FUN = sum)

#pre filtering
merged_pre_filter_1sk <- aggregate(nt_count ~ source + category, data = batch_1sk, FUN =sum)
#uncategorized reads + categorized reads no filtering step
merged_pre_filter_1sk_s <- rbind(na_top_1sk,merged_pre_filter_1sk)
#final merge with reads per sample
merged_pre_filter_1sk_reads <- merge(merged_pre_filter_1sk_s,total_reads_per_sample_batch_1_swab, by ='source')
#reads per microbes
batch_1_swab_microbes_unfiltered <- aggregate(nt_count ~ category, data = merged_pre_filter_1sk_s , FUN = sum)

#filtering step
batch_1sk %>%
  filter(category=="archaea") %>%
  filter(tax_level==1) %>%
  filter(is_phage=='FALSE') %>%
  filter(nt_alignment_length >= CUTOFF.NT_L &
           nt_rpm >= CUTOFF.NT_rPM &
           nr_rpm >= CUTOFF.NR_rPM &
           nt_z_score >= CUTOFF.NT_Z) -> arc1s
batch_1sk %>%
  filter(category=="bacteria") %>%
  filter(tax_level==1) %>%
  filter(is_phage=='FALSE') %>%
  filter(nt_alignment_length >= CUTOFF.NT_L &
           nt_rpm >= CUTOFF.NT_rPM &
           nr_rpm >= CUTOFF.NR_rPM &
           nt_z_score >= CUTOFF.NT_Z) -> bak1s
batch_1sk %>%
  filter(category=="eukaryota") %>%
  filter(tax_level==1) %>%
  filter(is_phage=='FALSE') %>%
  filter(nt_alignment_length >= CUTOFF.NT_L &
           nt_rpm >= CUTOFF.NT_rPM &
           nr_rpm >= CUTOFF.NR_rPM &
           nt_z_score >= CUTOFF.NT_Z) -> euk1s 
batch_1sk %>%
  filter(category=="viruses") %>%
  filter(tax_level==1) %>%
  filter(is_phage=='FALSE') %>%
  filter(nt_alignment_length >= CUTOFF.NT_L &
           nt_rpm >= CUTOFF.NT_rPM &
           nr_rpm >= CUTOFF.NR_rPM &
           nt_z_score >= CUTOFF.NT_Z) -> vik1s

#merging metrics - NT rPM 10
merged_1s <- rbind(arc1s,bak1s,euk1s,vik1s)%>%
  group_by(source,category,nt_count)%>%
  tally()%>%
  print(n=Inf) ->bk1s

#aggregating to get total number of reads per microbe type- NT rPM 10
s1s <- aggregate(nt_count ~ source + category, data = bk1s, FUN = sum)
s1sf <- merge(total_reads_per_sample_batch_1_swab,s1s, by='source')
#merging filtered and uncategorized reads. 
s1k_uncat_filtered_reads <- rbind(na_top_1sk,s1s)
#merging the data with total reads
s1k_filtered_final_reads <- merge(s1k_uncat_filtered_reads,total_reads_per_sample_batch_1_swab, by ='source')
#Figure 2 - filtered
batch_1_swab_microbes_filtered <- aggregate(nt_count ~ category, data = s1k_uncat_filtered_reads, FUN = sum)

#Final Run - Batch 2 
batch_2sk<- fs::dir_ls("prism_3_batch_2_novaseq_2373/swab/", regexp="\\.csv$")
#reading the files just by calling a function purr
batch_2sk %>%
  map_dfr(read_csv, .id = "source") -> batch_2sk
#renaming the samples
batch_2sk$source <-gsub("prism_3_batch_2_novaseq_2373/swab/","",as.character((batch_2sk$source)))
batch_2sk$source <-gsub(".csv","",as.character((batch_2sk$source)))

#filtering only at the species level
batch_2sk %>%
  filter(tax_level==1) -> batch_2sk

#uncategorized reads - species level
batch_2sk%>%
  filter(is.na(category)) -> na_2sk
na_2sk$category=c('uncategorized')
#summing all the na by sample type
na_top_2sk <- aggregate(nt_count~source+category, data = na_2sk, FUN = sum)

#total reads
total_reads_per_sample_batch_2_swab <- aggregate(nt_count ~ source, data = batch_2sk, FUN = sum)
#pre filtering
merged_pre_filter_2sk <- aggregate(nt_count ~ source + category, data = batch_2sk, FUN =sum)
#uncategorized reads + categorized reads no filtering step
merged_pre_filter_2sk_s <- rbind(na_top_2sk,merged_pre_filter_2sk)
#final merge with reads per sample
merged_pre_filter_2sk_reads <- merge(merged_pre_filter_2sk_s,total_reads_per_sample_batch_2_swab, by ='source')
#reads per microbes
batch_2_swab_microbes_unfiltered <- aggregate(nt_count ~ category, data = merged_pre_filter_2sk_s , FUN = sum)

#filtering
batch_2sk %>%
  filter(category=="archaea") %>%
  filter(tax_level==1) %>%
  filter(is_phage=='FALSE') %>%
  filter(nt_alignment_length >= CUTOFF.NT_L &
           nt_rpm >= CUTOFF.NT_rPM &
           nr_rpm >= CUTOFF.NR_rPM &
           nt_z_score >= CUTOFF.NT_Z) -> arc2s
batch_2sk %>%
  filter(category=="bacteria") %>%
  filter(tax_level==1) %>%
  filter(is_phage=='FALSE') %>%
  filter(nt_alignment_length >= CUTOFF.NT_L &
           nt_rpm >= CUTOFF.NT_rPM &
           nr_rpm >= CUTOFF.NR_rPM &
           nt_z_score >= CUTOFF.NT_Z) -> bak2s
batch_2sk %>%
  filter(category=="eukaryota") %>%
  filter(tax_level==1) %>%
  filter(is_phage=='FALSE') %>%
  filter(nt_alignment_length >= CUTOFF.NT_L &
           nt_rpm >= CUTOFF.NT_rPM &
           nr_rpm >= CUTOFF.NR_rPM &
           nt_z_score >= CUTOFF.NT_Z) -> euk2s 
batch_2sk %>%
  filter(category=="viruses") %>%
  filter(tax_level==1) %>%
  filter(is_phage=='FALSE') %>%
  filter(nt_alignment_length >= CUTOFF.NT_L &
           nt_rpm >= CUTOFF.NT_rPM &
           nr_rpm >= CUTOFF.NR_rPM &
           nt_z_score >= CUTOFF.NT_Z) -> vik2s

#merging metrics - NT rPM 10
merged_2s <- rbind(arc2s,bak2s,euk2s,vik2s)%>%
  group_by(source,category,nt_count)%>%
  tally()%>%
  print(n=Inf) ->bk2s

#aggregating to get total number of reads per microbe type- NT rPM 10
s2s <- aggregate(nt_count ~ source + category, data = bk2s, FUN = sum)
s2sf <- merge(total_reads_per_sample_batch_2_swab,s2s, by='source')
#merging filtered and uncategorized reads. 
s2k_uncat_filtered_reads <- rbind(na_top_2sk,s2s)
#merging the data with total reads
s2k_filtered_final_reads <- merge(s2k_uncat_filtered_reads,total_reads_per_sample_batch_2_swab, by ='source')
#Figure 2 - filtered
batch_2_swab_microbes_filtered <- aggregate(nt_count ~ category, data = s2k_uncat_filtered_reads, FUN = sum)


#Figure 2E
#importing metadata
###################################
## Read in mNGS sample meta-data_Batch 1 - qPCR updated results ##
###################################

dat_meta_mNGS_1_pf <- read_excel("Batch_1_2021-05/PRISM_mNGS_manifest_2021-05_final_with_Grant_variables_qPCR_update.xlsx") %>%
  mutate(reqdate = ymd(reqdate)) %>%
  mutate(dob = ymd(dob)) %>%
  mutate(enrolldate = ymd(enrolldate)) %>%
  mutate(agecat = factor(agecat, levels=c("< 5 years", "5-15 years", "16 years or older"))) %>%
  mutate(agecat_enrollment = factor(agecat_enrollment, levels=c("< 5 years", "5-15 years", "16 years or older"))) %>%
  mutate(hhid = factor(hhid))

###################################
## Read in mNGS sample meta-data_Batch 2 - qPCR updated results ##
###################################

dat_meta_mNGS_2_pf <- read_excel("Batch_2_2021-09/PRISM_mNGS_manifest_2021-09_final_with_Grant_variables_qPCR_update_2022-04-25.xlsx") %>%
  mutate(reqdate = ymd(reqdate)) %>%
  mutate(dob = ymd(dob)) %>%
  mutate(enrolldate = ymd(enrolldate)) %>%
  mutate(agecat = factor(agecat, levels=c("< 5 years", "5-15 years", "16 years or older"))) %>%
  mutate(agecat_enrollment = factor(agecat_enrollment, levels=c("< 5 years", "5-15 years", "16 years or older"))) %>%
  mutate(hhid = factor(hhid))

#combining metadata from batch 1 and batch 2 - double check this
final_metadata <- bind_rows(dat_meta_mNGS_1_pf,dat_meta_mNGS_2_pf)

#virus graph, merging viral microbes detect with metadata. 
#final dataframe for graph
viral_list <- rbind(vik1p,vik2p,vik1s,vik2s)
names(viral_list)[1] <- "barcode"
final_viral_list <- left_join(viral_list,final_metadata, by="barcode")

final_viral_list%>%
  filter(tax_level=="1") %>%
  group_by(name, barcode, diagnosis1, dx1code_name, dx2code_name, dx3code_name, nt_rpm, sampletype, cohortid, reqdate) %>%
  tally() %>%
  print(n=Inf) -> fvl_1

fvl_1$name[which(fvl_1$name=="Betacoronavirus 1")] <- "Human coronavirus OC43"
fvl_1$name[which(fvl_1$name=="Severe acute respiratory syndrome-related coronavirus")] <- "SARS-CoV-2"

#reading coverage
coverage <- read.csv("coverage_viral.csv")
## Rename "Betacoronavirus 1" to "Human coronavirus OC43"
coverage$name[which(coverage$name=="Betacoronavirus 1")] <- "Human coronavirus OC43"
coverage$name[which(coverage$name=="Severe acute respiratory syndrome-related coronavirus")] <- "SARS-CoV-2"

ff <- merge(fvl_1,coverage, by = c("barcode","name"))
#getting n_count
ff$n_detected <- table(ff$name)[ff$name]
ff$name_detect <- paste0(ff$name," (",ff$n_detected,")")
ff %>%
  group_by(name,cohortid,reqdate)%>%
  tally() -> unique_visits

unique_visits%>%
  group_by(name)%>%
  tally()-> unique_visits_f

#merging the unique visits
ff_f <- merge(ff,unique_visits_f,by='name')
#new column
ff_f$name_detect_visit <-paste0(ff_f$name_detect," (",ff_f$n.y,")")

#change the ordering of the viruses to relatedness on the y-axis.on for the name_detect column
ff_f%>%
  arrange(nt_rpm)%>%
  mutate(name_detect=factor(name_detect_visit, levels = c("Cardiovirus B (1) (1)","Human polyomavirus 3 (1) (1)","Mamastrovirus 1 (1) (1)","Norwalk virus (1) (1)","Human coronavirus OC43 (4) (4)","Human coronavirus NL63 (1) (1)","Human coronavirus HKU1 (1) (1)","Severe acute respiratory syndrome-related coronavirus (11) (11)","Pegivirus A (14) (14)","Pegivirus C (19) (19)","Influenza A virus (14) (14)","Human respirovirus 1 (12) (12)","Human respirovirus 3 (7) (6)","Enterovirus A (7) (5)","Enterovirus B (3) (2)","Rhinovirus A (8) (8)","Rhinovirus B (7) (7)", "Rhinovirus C (28) (25)","Human metapneumovirus (9) (9)","Human orthopneumovirus (12) (12)","Human orthorubulavirus 2 (1) (1)","Rotavirus A (7) (4)","Human betaherpesvirus 5 (2) (2)","Hepatitis GB virus B (1) (1)","Human bocavirus (1) (1)","Human mastadenovirus C (5) (5)"))) ->ttf

###### - FINAL Graph
ttt<-ggplot(data=ttf, mapping = aes(x = nt_rpm, y = name_detect_visit, col = breadth, shape = sampletype)) + geom_point(size=4,colour="black") + geom_point(size=3) +  
  scale_x_log10() + ylab ("Viral microbe") + xlab ("Nucleotide Reads Per Million (NT rPM)") + scale_colour_gradient2(midpoint = 20) + theme(axis.text = element_text(size=12)) 
ggsave("figure_viral_pathogens_species_nc_outline_group_sample_type_count_unique_visits_v2.png", ttt, width = 40 , height = 20, units ="cm", device = 'png', dpi=700)

