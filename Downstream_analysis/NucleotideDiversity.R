#Code to compute nucleotide diversity, as most methods do not work with polyploidy
#!/usr/bin/env Rscript
# @Author: Jose Miguel Ramirez; Adapted by Winona Oliveros
# @E-mail: jose.ramirez1@bsc.es
# @Description: Code to prepare the necessary metadata per tissue
# @software version: R=4.2.2
rm(list=ls())
#Load libraries
library(tidyr)
library(dplyr)
library(purrr)
library(ggplot2)
library(patchwork)

# setwd(system("pwd", intern = T)) #If in linux
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) #If using Rstudio

#Reading 1KG variants
# data <- read.delim("../../Winona/variantCalling/1000G/merged.txt")
data <- read.delim("calling_1KG/merged_with_DP.txt")
data <- data[order(data$POS), ]

#Computing allele frequency per donor:
computing_allele_frequency <- function(data) {
  # Identify GT columns
  gt_cols <- grep(".GT", colnames(data), fixed = TRUE)
  new_colnames <- gsub(".GT", ".AF", colnames(data)[gt_cols], fixed = TRUE)
  
  # Initialize output
  output <- data.frame(POS = data$POS,
                       REF = data$REF,
                       ALT = data$ALT,
                       matrix(NA_character_, nrow = nrow(data), ncol = length(new_colnames)))
  colnames(output)[4:ncol(output)] <- new_colnames
  
  # Process donors (step by 4: REF, ALT, GT, DP — adjust if different)
  for (sample in seq(6, ncol(data), by = 4)) {
    donor <- gsub(".ALT_D", ".AF", colnames(data)[sample], fixed = TRUE)
    
    # Get allele counts
    alt_counts <- data[[sample]]
    ref_counts <- data[[sample - 1]]
    dp <- data[[sample + 1]]
    
    # Replace NAs with 0
    alt_counts[is.na(alt_counts)] <- 0
    ref_counts[is.na(ref_counts)] <- 0
    dp[is.na(dp)] <- NA
    
    # Compute frequencies
    ref_freq <- ref_counts / dp
    alt_freq <- alt_counts / dp
    
    # Store in donor column as "ref,alt"
    output[[donor]] <- paste0(ref_freq, ",", alt_freq)
  }
  
  return(output)
}

af <- computing_allele_frequency(data)


# Computing nucleotide diversity score:
# Number of comparisons. The value to normalize for if we used the original formular Number of pairwise differences / Number of combinations
# n_samples <- ncol(af) -3
# N <- (n_samples) * (n_samples-1) / 2

pis <- data.frame()
for(pos in unique(af$POS)){ #We get one value per variant and then we average over certain regions
  print(pos)
  subset <- af[af$POS==pos,]
  frequencies <- subset[,-c(1:3)]
  
  # Transpose so each sample is a row
  freq_long <- as.data.frame(t(frequencies))
  freq_long$sample <- rownames(freq_long)
  
  # Split each cell "REF,ALT" into two columns
  split_freq <- purrr::reduce(colnames(freq_long[!colnames(freq_long) %in% "sample"]), function(df, col) { #Split a column of type REF,ALT into two columns: REF and ALT
    separate(df, col, into = paste0(col, c("_REF", "_ALT")), sep = ",", remove = TRUE)
  }, .init = freq_long)
  
  #Get as the reference allele, the one that is not NA, if all are NA, then 1:
  if(nrow(frequencies)==1){
    freq <- split_freq[,c(1,grep("ALT",colnames(split_freq)))] #keep only one reference allele
    freq[,1][freq[,1]=="NA"] <- 1
    freq[,2][freq[,2]=="NA"] <- 0
  }else{
    ref <- split_freq[,grep("REF",colnames(split_freq))]
    ref$first_non_na <- apply(ref, 1, function(x) x[which(x!="NA")[1]])
    ref$first_non_na[is.na(ref$first_non_na)] <- 1
    
    freq <- cbind(ref$first_non_na, split_freq[,c(grep("ALT",colnames(split_freq)))]) #keep only one reference allele
  }
  freq <- freq %>%
    mutate(across(everything(), as.numeric)) #transform into numeric, NA characters are transformed into real NAs
  freq[is.na(freq)] <- 0
  mean_af <- colMeans(freq) #Computing mean allele frequency per allele
  
  #Traditionally, for nucleotide diversity, the number of pairwise difference is computed, but to get the equivalent using frequencies we compute:
  #2*(all pairwise combinatios of allele products)
  #pi <- 2*sum(combn(mean_af, 2, prod))
  #which is equivalent to:
  pi <- 1 - sum(mean_af^2)
  # pi <- pairwise_differences/N #The original formula would be number of pairwise differences divided by number of pairwise combinations
  pis <- rbind(pis, c(pos, pi))
  
}

saveRDS(pis, "pis_mean.rds")
pis <- readRDS("pis_mean.rds")

# Compare average nucleotide diversity in regions, as it is more reliable than per variant
assign_region <- function(position){
  # <= only in IGS at start but all other should be > only but then variants are 0 bed based
  #so variant at the end of region be annotated to be in wrong region so we should use >=
  if(position <= 9338){
    "IGS"
  }else if(position <= 12995){
    "5'ETS"
  }else if(position <= 14864){
    "18S"
  }else if(position <= 15934){
    "ITS1"
  }else if(position <= 16091){
    "5.8S"
  }else if(position <= 17258){
    "ITS2"
  }else if(position <= 22309){
    "28S"
  }else if(position <= 22670){
    "3'ETS"
  }else{
    "IGS"
  }
}

colnames(pis) <- c("Pos", "Pi")

data <- pis %>%
  rowwise() %>%
  mutate(Region = assign_region(Pos))

### Using all the positions

data <- data %>%
  rowwise() %>%
  mutate(Region = assign_region(Pos))

data$Region <- factor(data$Region, levels = unique(data$Region))
colors_regions <- c("5'ETS" = "#5A0002","18S" = "#A40606","ITS1" = "#437C90","5.8S" = "#D98324","ITS2" = "#255957","28S" = "#60D394", "3'ETS" = "#DAD6D6")

average <- data %>%
  group_by(Region) %>%
  summarize(Mean = mean(Pi, na.rm=TRUE))


plot <- ggplot(average[average$Region!='IGS',]) + geom_col(aes(Region, Mean, fill=Region)) + ylab("Mean nucleotide diversity") + 
  scale_fill_manual(values=colors_regions) + theme_bw() +
  theme(
    axis.text.x = element_text(size = 12),  # Increase x-axis text size
    axis.text.y = element_text(size = 12),  # Increase y-axis text size
    axis.title.x = element_text(size = 13), # Increase x-axis title size
    axis.title.y = element_text(size = 13), # Increase y-axis title size
    plot.title = element_text(size = 16),   # Increase plot title size
    legend.position = "none"               # Place legend on the right
  ) 
pdf("barplot_nucleotide_diversity_per_region.pdf", width = 4, height = 3)

print(plot)

dev.off()

#### nucleotide diversity ES ####
#Plot number of positions overlapping expansion segments
expansion_segments <- read.delim("../04_Pipeline/Plots/expansion_segments.txt", sep = "") #Obtained from code
#Obtained from: 13.Overlap_with_expansion_segments.R
data$ES <- "No"

for(pos in unique(data$Pos)){
  es <- expansion_segments$ES[pos > expansion_segments$Start_chrR & pos <= expansion_segments$End_chrR]
  if(length(es)!=0){
    data$ES[data$Pos==pos] <- es
  }
}
data$ES_any <- data$ES!="No"

average <- data %>%
  group_by(Region,ES_any) %>%
  summarize(Mean = mean(Pi, na.rm=TRUE))


plot <- ggplot(average[average$Region %in% c('5.8S','28S','18S'),]) + geom_col(aes(Region, Mean, fill=ES_any),position = "dodge") + ylab("Mean nucleotide diversity") + 
  scale_fill_manual(values = c("TRUE" = "dodgerblue4", "FALSE" = "firebrick"),
                    labels = c("TRUE" = "Yes", "FALSE" = "No")) +
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 12),  # Increase x-axis text size
    axis.text.y = element_text(size = 12),  # Increase y-axis text size
    axis.title.x = element_text(size = 13), # Increase x-axis title size
    axis.title.y = element_text(size = 13), # Increase y-axis title size
    plot.title = element_text(size = 16)#,   # Increase plot title size
    #legend.position = "none"               # Place legend on the right
  ) 
pdf("barplot_nucleotide_diversity_per_region_ES_segments.pdf", width = 4, height = 3)

print(plot)

dev.off()

