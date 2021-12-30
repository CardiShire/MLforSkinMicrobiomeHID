library(tidyverse)
library(glmnet)
library(doParallel)
registerDoParallel(9)

args = commandArgs(trailingOnly=TRUE)

ind = args[1]


# TestIND <- c("S001", "S002", "S004", "S006", "S007", "S010", "S011", "S012", "S015",
#               "S016", "S017", "S025", "S028", "S030", "S031", "S032", "S033", "S038",
#               "S040", "S042", "S044", "S045", "S046", "S047", "S049", "S051")
# test <- c("S003", "S005", "S008", "S009", "S013", "S014", "S018", "S019", "S020", "S021",
#           "S022", "S023", "S024", "S026", "S027", "S029", "S034", "S035", "S036", "S037",
#           "S039", "S041", "S043", "S048", "S050")

# import data frame of all potential SNPs and there nucleotide frequencies for Test data set
SNPmatrixML <- read.csv(paste0("~/src/Fstiminator/FstCalculations/UpdatedLASSO/TrainingFiles/Train_AllelePull/AlleleFreq_", ind,"/CombinedAlleleFreq_HeldOut", ind, "12212021.csv"))

# SNPmatrixML <- SNPmatrixML[-226:-234,]

SNPmatrixML<- arrange(SNPmatrixML, desc(Sample))
# Users can explicitly control the fold that each observation is assigned to the via the foldid argument.
SNPmatrixML <- SNPmatrixML %>% mutate(foldid = rep(1:9, times=26, each=1), .after = Sample)
# separate sample information so individual number has a column
matrixML <- SNPmatrixML %>% separate(Sample, c("Ind", "Bodysite", "Time", "Rep"), "-")
# make any NAs in the df 0
matrixML[is.na(matrixML)] <- 0


samplelist <- matrixML$Ind
samplelist <- unique(samplelist)


ID <- dplyr::mutate(matrixML, ID = row_number())
ID <- ID %>% select(Ind, ID)


reducefeatures <- function(holdout = samplelist, matrixML) {
  # convert to matrix for analysis 
  # remove first column with sample names
  x = as.matrix(matrixML[-1:-5])
  
  # Remove one sample from matrix for initial Test
  trainmatrix <- matrixML %>% filter(Ind != holdout)
  
  trainy <- as.vector(trainmatrix$Ind)
  trainx <- as.matrix(trainmatrix[-1:-6])
  
  

  # Making custom list of 100 lambda's to test for optimization
  grid=1.1^seq(1,-200,length=100)
  
  set.seed(1)
  fit.lo <- cv.glmnet(x = trainx, y = trainy, family="multinomial", 
                      alpha=1, standardize = TRUE, keep = TRUE,
                      lambda = grid, type.multinomial = "grouped", 
                      maxit = 1000000, foldid = trainmatrix$foldid, parallel = TRUE)
  slam <- as.numeric(which.max(fit.lo$nzero)) + 1
  
  myCoefs <- coef(fit.lo, s = "lambda.min")

  # dimnames is the same for every fold of the test, but the coefficient numbers change
  myCoefs1 <- myCoefs[[1]]
  coefsmin <- myCoefs1[which(myCoefs1 != 0)]
  coefslabel <- myCoefs1@Dimnames[[1]][which(myCoefs1 != 0)]
  reducedsnplist <- as.data.frame(cbind(coefslabel, coefsmin))
  # making data frame for specific SNPs to add all allele frequencies
  outlist <- as.data.frame(coefslabel)
  outlist <- suppressWarnings(outlist %>% separate(coefslabel, into = c("Marker", "secn", "allele"), sep = "_"))
  outlist <- outlist %>% select("Marker", "secn")
  outmarker <- outlist %>% unite('file', Marker:secn, remove = TRUE)
  outmarker <- unique(outmarker)
  outmarker <- outmarker[-1,]
  alleles <- as.list(c("A", "T", "C", "G"))
  outsnp <- expand.grid(outmarker, alleles)
  #outsnp <- outsnp %>% unite('snp', Var1:Var2, remove = TRUE)
  outsnp$ind <- holdout
  outsnp <- as.data.frame(outsnp)
  
  return(outsnp)

}

secondcv_holdout <- function(outsnp, holdout = ht, matrixML) {
  outsnp1 <- as.data.frame(outsnp)

  outsnp1 <- outsnp1 %>% mutate(adt = "12212021.tsv")
  outsnpcount <- paste(outsnp1$Var1, outsnp1$adt, outsnp1$Var2, sep = "_")
  
  # preparing list of markers to be used with all initial markers 
  matchingList <- as.list(outsnpcount)
  matchingList <- unlist(matchingList)

  # Remove one sample from matrix for initial Test
  trainmatrix <- matrixML %>% filter(Ind != holdout)
  # subset list of markers for variable analysis
  newDF <- trainmatrix[ ,which((colnames(trainmatrix) %in% matchingList)==TRUE)]
  newDF <- as.matrix(newDF)
  
  trainy <- unlist(trainmatrix$Ind)

  # Making custom list of 100 lambda's to test for optimization
  grid=1.1^seq(1,-200,length=100)

  
  # now to test data classification with first removed individual
  # the line below should produce a dgCMatrix
  # alpha set to 0 to use all selected markers
  fit.pre <- cv.glmnet(x = newDF, y = trainy, family="multinomial", alpha=0, standardize = TRUE, keep = TRUE,
                      lambda = grid, type.multinomial = "grouped", maxit = 1000000, foldid = trainmatrix$foldid, parallel = TRUE)

  lam <- fit.pre$lambda.min

  return(lam)
}

predictsample <- function(holdout, matrixML) {

  ID <- dplyr::mutate(matrixML, ID = row_number())
  ID <- ID %>% select(Ind, ID)
  rowID <- ID %>% filter(Ind == holdout)
  rowID <- unlist(rowID$ID)
  
  return(rowID)
}

sampleID <- function(outsnp, rowID, lam, matrixML)  {

  samplelist <- matrixML$Ind
  samplelist <- unique(samplelist)
  
  outsnp1 <- as.data.frame(outsnp)
  outsnp1 <- outsnp1 %>% mutate(adt = "12212021.tsv")
  outsnpcount <- paste(outsnp1$Var1, outsnp1$adt, outsnp1$Var2, sep = "_")
  
  # preparing list of markers to be used with all initial markers 
  matchingList <- as.list(outsnpcount)
  matchingList <- unlist(matchingList)
  
  newDF <- matrixML[ ,which((colnames(matrixML) %in% matchingList)==TRUE)]
  newDF <- as.matrix(newDF)
  testy <- unlist(matrixML$Ind)
  testy <- as.data.frame(testy)
  
  
  set.seed(1)
  fpre <- glmnet(x = newDF[-rowID,], y = testy[-rowID,], family="multinomial", 
                 alpha=0, lambda = lam, type.multinomial = "grouped", maxit = 1000000)
  
  predict <-predict(fpre, newx = newDF, s = lam, type = "class")

  predictions <- cbind(testy, predict)
  singlepredict <- predictions[rowID,]
  return(singlepredict)
}


outsnp <- map(.x = ind, ~ reducefeatures(.x, matrixML))

lam <- map2(outsnp, ind, ~ secondcv_holdout(.x, .y, matrixML))

rowID <- map(ind, ~ predictsample(.x, matrixML))



providenumber <- function(i){
  foo <- map_df(.x = rowID[[i]],  ~ sampleID(outsnp[[i]], .x, lam[[i]], matrixML))
  return(foo)
}

i = (1:length(ind))

class_predict <- map(i, providenumber)
##### works with single ind = "SOXX" to this point. out puts class_predict as list of lists


snplist <- mapply(cbind, outsnp[i], "HeldSample"= ind, SIMPLIFY=F)
outsnplist <- bind_rows(snplist, .id = "column_label")
fullpanel <- do.call(rbind, outsnp)
fullpanel <- fullpanel %>% group_by(Var2) %>% mutate(count = n())
fullpanel2 <- apply(fullpanel,2,as.character)
# fullpanel <- unique(fullpanel)
write.csv(fullpanel2, paste0("~/src/Fstiminator/FstCalculations/UpdatedLASSO/TrainingFiles/Train_MLOuputPanels/SNP_Panel_HeldOut", ind, "_12262021.csv"), row.names = FALSE)

outsnp1 <- fullpanel %>% mutate(adt = "12212021.tsv")
outsnpcount <- paste(outsnp1$Var1, outsnp1$adt, outsnp1$Var2, sep = "_")
newDF <- SNPmatrixML[ ,which((colnames(SNPmatrixML) %in% outsnpcount)==TRUE)]
newDF <- as.matrix(newDF)
st <- SNPmatrixML$Sample
fulldf <- cbind(st, newDF)
fulldf <- as.data.frame(fulldf)


forNA <- fulldf
# determine number of NA in each row
forNA$na_count <- apply(fulldf, 1, function(x) sum(is.na(x)))
forNA <- forNA %>% select(st, na_count)
# total number of ATCG missing for each sample that has greater than 0

# total number of markers missing for each sample
forNA <- forNA %>% mutate(na_marker = na_count/4)
forNA <- forNA %>% rename(Sample = st)
forNA$HeldOut <- ind
write.csv(forNA, paste0(
  "~/src/Fstiminator/FstCalculations/UpdatedLASSO/TrainingFiles/Train_MLOuputPanels/NA_HeldOut", ind, "_12262021.csv"),
  row.names = FALSE)



class_predict_all <- do.call(rbind, class_predict)

# class_predict_all <- cbind(st, class_predict_all)
class_predict_all <- class_predict_all %>% rename(Class = testy, Classification = `1`)
# coblam <- left_join(class_predict_all, forNA)
write.csv(class_predict_all, paste0(
          "~/src/Fstiminator/FstCalculations/UpdatedLASSO/TrainingFiles/Train_MLOuputPanels/Prediction_HeldOut", ind, "_12262021.csv"),
          row.names = FALSE)

lam2 <- unlist(lam)
lamb <- as.data.frame(cbind(samplelist, lam2))
write.csv(lamb,paste0(
  "~/src/Fstiminator/FstCalculations/UpdatedLASSO/TrainingFiles/Train_MLOuputPanels/Lam_HeldOut", ind, "_12262021.csv"),
  row.names = FALSE)
