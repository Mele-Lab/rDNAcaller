#Code to merge individual vcf calls into a merged table

setwd(system("pwd", intern = T)) #If in linux
# setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) #If using Rstudio

args = commandArgs(trailingOnly=TRUE)

list <- args[1] #This is a .txt file with the name of all samples that we want to merge
output <- args[2] #This is the output directory

list <- read.table(list)


suppressWarnings(library(data.table))     # For efficient data handling

process_sample <- function(sample){
  message("Processing sample: ", grep(sample, samples))
  #Reading vcf file
  path <- paste0(output, "/", sample, ".vcf.gz")
  if(file.size(path)<4122){ #If the size file is 0, we return NULL, if less than 4122 in Haplotypecaller it only has headers
    return(NULL)
  }
  vcf <- fread(path,select = c(1, 2, 4, 5, 10),
               col.names = c("CHR", "POS", "REF", "ALT", sample))
  #Procesing vcf file
  # Split genotype fields using fast data.table::tstrsplit
  vcf[, c("GT", "AD", "DP", "GQ") := tstrsplit(get(sample), ":", fixed = TRUE)] #If HaplotypeCaller
  # vcf[, c("GT", "AD", "DP") := tstrsplit(get(sample), ":", fixed = TRUE)] #If mutect
  vcf[, c("REF_D", "ALT_D") := tstrsplit(AD, ",", fixed = TRUE)]
  
  # Keep minimal columns
  vcf <- vcf[, .(POS, REF, ALT,
                 GT, REF_D, ALT_D, DP)]
  colnames(vcf)[4:7] <- paste0(sample, ".", colnames(vcf)[4:7])
  return(vcf)
}
samples <- list$V1
print("about to start")
result_list <- lapply(samples, process_sample)
# Remove NULLs in case any files didn't exist
result_list <- result_list[!sapply(result_list, is.null)]

print("about to reduce")
# Combine all data frames, including the first sample. 
# Collect all variant keys first
dt_list <- lapply(result_list, as.data.table)
positions <- unique(rbindlist(lapply(dt_list, function(x) x[, .(POS, REF, ALT)])))

# Initialize final table
final_vcf <- copy(positions)

# Loop through samples
for (i in seq_along(dt_list)) {
  message("Merging sample: ", i)
  dt <- dt_list[[i]]
  setkey(dt, POS, REF, ALT)
  
  # Merge this sample’s 4 columns
  final_vcf <- dt[final_vcf]
}


#Save final vcf (even though it is no longer a vcf format, so I will save it as txt)
print(paste0(output, "/merged_with_DP.txt"))
write.table(final_vcf, paste0(output, "/merged_with_DP.txt"), row.names = F, quote = F, sep = "\t")
