suppressPackageStartupMessages(library(tidyverse))
library(stringr)

ind = "S051"
# trainIND <- c("S001", "S002", "S004", "S006", "S007", "S010", "S011", "S012", "S015",
#               "S016", "S017", "S025", "S028", "S030", "S031", "S032", "S033", "S038",
#               "S040", "S042", "S044", "S045", "S046", "S047", "S049", "S051")
# test <- c("S003", "S005", "S008", "S009", "S013", "S014", "S018", "S019", "S020", "S021",
#           "S022", "S023", "S024", "S026", "S027", "S029", "S034", "S035", "S036", "S037",
#           "S039", "S041", "S043", "S048", "S050")

# Import list of nucleotide positions with their Fst estimates and count
list <- read_table2(paste0("~/src/Fstiminator/FstCalculations/UpdatedLASSO/TrainingFiles/Fst/Train_FstList_HeldOut", ind, ".fst"))



list <- list %>% mutate(avgfst = SumFST/Count)
filterlist <- list %>% filter(avgfst >= 0.1)
maxlist <- max(list$Count)
filterlist <- filterlist %>% filter(Count >= maxlist*0.95)

fileout <- paste0("~/src/Fstiminator/FstCalculations/UpdatedLASSO/TrainingFiles/Train_AllelePull/AlleleFreq_", ind, "/",
                             "HeldOut", ind, "_top95_12192021.csv", sep = "")
write.csv(filterlist, fileout, row.names = FALSE)

# Select the specific columns need to format marker names
listofsnps<- filterlist %>% select("Marker", "SNP")

# Split column into multiple columns contained sections of the file name
markerposition <- separate(
  data = listofsnps, col = "Marker", 
  into = c("gi", "nu", "ref", "mr", "loc"), sep = "\\|", remove = FALSE)
# Bring part of the sections back together for marker name
markerposition$com <- pmap(markerposition[5:7], paste, sep = '-')

# Make a list of just the marker names
myvar <- list(markerposition)
for (i in seq(nrow(markerposition)) ) { myvar[[i]]  <-  unclass( markerposition[i,] )}
# Add system date
d = Sys.Date()
# Make lines for searching for the markers in all Bam files
capture.output(cat(map_chr(
  myvar, ~ (paste(
    'zgrep -n "', .x$Marker, " ", .x$SNP,
    ' " *27.gz > ', '"', 
    "~/src/Fstiminator/FstCalculations/UpdatedLASSO/TrainingFiles/Train_AllelePull/AlleleFreq_", ind, "/",
    .x$com, '_12142021.csv"', sep = ""))), sep = "\n"), 
    file = paste0("~/src/Fstiminator/FstCalculations/UpdatedLASSO/TrainingFiles/Train_AllelePull/AlleleFreq_", 
    ind, "/", "AlleleFreq_HeldOut", ind, ".sh")) 

#may have do reset ulimit to 'ulimit -s 1000000' it is orginially 'ulimit -s 8192'