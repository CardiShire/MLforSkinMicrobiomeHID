library(tidyverse)

#sample(lt = c(1:51), 26) was used to choose random samples
#Training individuals
trainIND <- c("S001", "S002", "S004", "S006", "S007", "S010", "S011", "S012", "S015",
              "S016", "S017", "S025", "S028", "S030", "S031", "S032", "S033", "S038",
              "S040", "S042", "S044", "S045", "S046", "S047", "S049", "S051")

#Test individuals
test <- c("S003", "S005", "S008", "S009", "S013", "S014", "S018", "S019", "S020", "S021",
          "S022", "S023", "S024", "S026", "S027", "S029", "S034", "S035", "S036", "S037",
          "S039", "S041", "S043", "S048", "S050")

comblist <- combn(trainIND, 2)
comblist <- as.data.frame(t(comblist))

combout <- comblist %>% mutate(reg = paste0(V1, ".+", V2))

createfilenames <- function(x){

  allfiles <- list.files(path = "~/src/Fstiminator/FstCalculations/AllBamFiles/", pattern= x)
  return(allfiles)
}

fn <- map(combout$reg, createfilenames)
filenames <- do.call(c, fn)



