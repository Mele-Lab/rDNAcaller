#!/usr/bin/env Rscript
# @Author: Jose Miguel Ramirez
# @E-mail: jose.ramirez1@bsc.es
# @Description: Code to prepare the supplementary table
# @software version: R=4.2.2

rm(list=ls())

library(dplyr)
library(fuzzyjoin)

#The pipeline output:
variants <- read.delim("calling_1KG/merged_with_DP.txt")


#Compute IGF
samples <- grep(".DP", colnames(variants), fixed = T, value=T)
samples <- gsub(".final.DP", "", samples)
# Compute AF columns
for (sample in samples) {
  alt <- variants[[paste0(sample, ".final.ALT_D")]]
  dp  <- variants[[paste0(sample, ".final.DP")]]
  variants[[paste0(sample, ".AF")]] <- ifelse(dp > 0, alt / dp, NA)  # avoid div by zero
}

# Keep only POS, REF, ALT, and AF columns
table <- variants[, c("POS", "REF", "ALT", paste0(samples, ".AF"))] #2470 samples, since we do not include the trios

#Sort table
table <- table[order(table$POS, table$REF, table$ALT), ]

#Add position starting from TSS:
table$POS_TSS <- table$POS - 9338 #We use chrR, which starts from a point in the IGS


#Add the region and the position with respect to the start of the region:
region_starts <- c(
  "IGS" = 1,
  "5'ETS" = 9339,
  "18S" = 12996,
  "ITS1" = 14865,
  "5.8S" = 15935,
  "ITS2" = 16092,
  "28S" = 17259,
  "3'ETS" = 22310
)

assign_region <- function(position){
  if(position <= 9338){"IGS"
  }else if(position <= 12995){
    region <- "5'ETS"
  }else if(position <= 14864){
    region <- "18S"
  }else if(position <= 15934){
    region <- "ITS1"
  }else if(position <= 16091){
    region <- "5.8S"
  }else if(position <= 17258){
    region <- "ITS2"
  }else if(position <= 22309){
    region <- "28S"
  }else if(position <= 22670){
    region <- "3'ETS"
  }else{
    region <- "IGS"
  }
  return(list(region = region, relative_position = position - region_starts[region] + 1))
}

table <- table %>%
  rowwise() %>%
  mutate(
    tmp = list(assign_region(POS)),
    Region = tmp$region,
    Region_position = tmp$relative_position
  ) %>%
  select(-tmp)


#Add the expansion segments and the position with respect to the start of the expansion segments:
expansion_segments <- read.delim("../04_Pipeline/Plots/expansion_segments.txt", sep = "")[,c(1,5,6)]
colnames(expansion_segments)[1] <- "Expansion_segment"

#The exapansion segments have been obtained from https://www.rcsb.org/structure/4V6X

# Example: table has POS, REF, ALT; expansion_segments has ES, Start_chrR, End_chrR
table <- fuzzy_left_join(
  table, 
  expansion_segments,
  by = c("POS" = "Start_chrR", "POS" = "End_chrR"),
  match_fun = list(`>=`, `<=`)  # POS >= Start_chrR & POS <= End_chrR
) %>%
  mutate(
    Expansion_segment_position = POS - Start_chrR + 1
  ) %>%
  select(-c(Start_chrR, End_chrR))


#Reorder the columns of the tables so the information of the variants is before the IGF
table <- table[,c(1:3,2471:2475, 4:2470)]
colnames(table)[1] <- "Position"
colnames(table)[4] <- "Position_from_TSS"


write.csv(table, "Supplementary_Table_1.csv", row.names = F)
