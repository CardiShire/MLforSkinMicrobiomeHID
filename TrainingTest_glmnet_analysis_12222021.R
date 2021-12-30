library(tidyverse)
library(glmnet)
library(doParallel)
registerDoParallel(30)
# holdout = "S001"
# rowheld = 1

# import data frame of all potential SNPs and there nucleotide frequencies for training data set
SNPmatrixML <- read.csv("CombinedAlleleFreq_All_12222021.csv")



SNPmatrixML<- arrange(SNPmatrixML, desc(Sample))
# Users can explicitly control the fold that each observation is assigned to the via the foldid argument.
SNPmatrixML <- SNPmatrixML %>% mutate(foldid = rep(1:9, times=51, each=1), .after = Sample)
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
  print("Inside reducefeatures")
  print(holdout)
  # Remove one sample from matrix for initial training
  trainmatrix <- matrixML %>% filter(Ind != holdout)
  trainy <- unlist(trainmatrix$Ind)
  trainx <- as.matrix(trainmatrix[-1:-6])
  
  # Making custom list of 100 lambda's to test for optimization
  grid=1.1^seq(1,-200,length=100)
  
  set.seed(1)
  fit.lo <- cv.glmnet(x = trainx, y = trainy, family="multinomial", alpha=1, standardize = TRUE, keep = TRUE,
                      lambda = grid, type.multinomial = "grouped", maxit = 1000000, foldid = trainmatrix$foldid, parallel = TRUE)
  
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
  outsnp <- outsnp %>% unite('snp', Var1:Var2, remove = TRUE)
  
  return(outsnp)
  
}

secondcv_holdout <- function(outsnp, holdout = ht, matrixML) {
  outsnp1 <- as.data.frame(outsnp)
  print("Inside secondcv")
  outsnp1 <- outsnp1 %>% mutate(adt = "12222021.tsv")
  outsnp1 <- outsnp1 %>% separate(snp, c("NC", "snp", "nuc"), sep = "_")
 
  outsnpcount <- paste(outsnp1$NC, outsnp1$snp, outsnp1$adt, outsnp1$nuc, sep = "_")
  # preparing list of markers to be used with all initial markers 
  matchingList <- as.list(outsnpcount)
  matchingList <- unlist(matchingList)
  
  
  # Remove one sample from matrix for initial test
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
  print(holdout)
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
  print("sampleID")
  samplelist <- matrixML$Ind
  samplelist <- unique(samplelist)
  
  outsnp1 <- as.data.frame(outsnp)
  
  outsnp1 <- outsnp1 %>% mutate(adt = "12222021.tsv")
  outsnp1 <- outsnp1 %>% separate(snp, c("NC", "snp", "nuc"), sep = "_")
  
  outsnpcount <- paste(outsnp1$NC, outsnp1$snp, outsnp1$adt, outsnp1$nuc, sep = "_")
  # preparing list of markers to be used with all initial markers 
  matchingList <- as.list(outsnpcount)
  matchingList <- unlist(matchingList)
  
  
  newDF <- matrixML[ ,which((colnames(matrixML) %in% matchingList)==TRUE)]
  newDF <- as.matrix(newDF)
  testy <- unlist(matrixML$Ind)
  testy <- as.data.frame(testy)
  
  grid=1.1^seq(1,-200,length=100)
  set.seed(1)
  fpre <- glmnet(x = newDF[-rowID,], y = testy[-rowID,], family="multinomial", 
                 alpha=0, type.multinomial = "grouped", lambda = grid)
  
  predict <-predict(fpre, newx = newDF, s = lam, type = "class")
  
  predictions <- cbind(testy, predict)
  singlepredict <- predictions[rowID,]
  print(rowID)
  return(singlepredict)
}

print("Hello")
outsnp <- map(.x = samplelist, ~ reducefeatures(.x, matrixML))
print(outsnp)


print("After outsnp")
lam <- map2(outsnp, samplelist, ~ secondcv_holdout(.x, .y, matrixML))
print("Lam")
rowID <- map(samplelist, ~ predictsample(.x, matrixML))
print("rowID")
providenumber <- function(i){
  foo <- map_df(.x = rowID[[i]],  ~ sampleID(outsnp[[i]], .x, lam[[i]], matrixML))
  return(foo)
}
print("providenumber")
i = (1:length(samplelist))
print("i")
class_predict <- map(i, providenumber)
print("class_predict")


snplist <- mapply(cbind, outsnp[i], "HeldSample"=samplelist[i], SIMPLIFY=F)
outsnplist <- bind_rows(snplist, .id = "column_label")
fullpanel <- do.call(rbind, outsnp)
fullpanel <- fullpanel %>% group_by(snp) %>% mutate(count = n())
# fullpanel <- unique(fullpanel)
write.csv(fullpanel, 
          "~/src/Fstiminator/FstCalculations/UpdatedLASSO/AllPotentialSNPs_Both.csv",
          row.names = FALSE)


fullpanel <- unlist(fullpanel[1])
newDF <- SNPmatrixML[ ,which((colnames(SNPmatrixML) %in% fullpanel)==TRUE)]
newDF <- as.matrix(newDF)
st <- SNPmatrixML$Sample
fulldf <- cbind(st, newDF)
fulldf <- as.data.frame(fulldf)


forNA <- fulldf
# determine number of NA in each row
forNA$na_count <- apply(fulldf, 1, function(x) sum(is.na(x)))
forNA <- forNA %>% select(st, na_count)
# total number of ATCG missing for each sample that has greater than 0
forNA <- forNA %>% filter(na_count > 0)
# total number of markers missing for each sample
forNA <- forNA %>% mutate(na_marker = na_count/4)
forNA <- forNA %>% rename(Sample = st)



class_predict_all <- do.call(rbind, class_predict)

class_predict_all <- cbind(st, class_predict_all)
class_predict_all <- class_predict_all %>% rename(Sample = st, Class = testy, Classification = `1`)
coblam <- left_join(class_predict_all, forNA)
write.csv(coblam, 
          "~/src/Fstiminator/FstCalculations/UpdatedLASS)/Both_Classification_Results.csv", 
          row.names = FALSE)

lam2 <- unlist(lam)
lamb <- as.data.frame(cbind(samplelist, lam2))
write.csv(lamb, 
          "~/src/Fstiminator/FstCalculations/Bam/Both_Lambda_Results.csv", 
          row.names = FALSE)


# Calculating missing SNPs for intial list of 757 SNPs
# A <- SNPmatrixML[ , grepl( "_A$" , names( SNPmatrixML ) ) ]
# T <- SNPmatrixML[ , grepl( "_T$" , names( SNPmatrixML ) ) ]
# G <- SNPmatrixML[ , grepl( "_G$" , names( SNPmatrixML ) ) ]
# C <- SNPmatrixML[ , grepl( "_C$" , names( SNPmatrixML ) ) ]
# Samples <- SNPmatrixML$Sample
# Samples <- as.data.frame(Samples)
# Ac <- cbind(Samples, A)
# Tc <- cbind(Samples, T)
# Gc <- cbind(Samples, G)
# Cc <- cbind(Samples, C)
# Ac$na_count <- apply(Ac, 1, function(x) sum(is.na(x)))
# Ac <- Ac %>% select(Samples, na_count)
# Tc$na_count <- apply(Tc, 1, function(x) sum(is.na(x)))
# Tc <- Tc %>% select(Samples, na_count)
# Gc$na_count <- apply(Gc, 1, function(x) sum(is.na(x)))
# Gc <- Gc %>% select(Samples, na_count)
# Cc$na_count <- apply(Cc, 1, function(x) sum(is.na(x)))
# Cc <- Cc %>% select(Samples, na_count)
# Acred <- Ac %>% filter(na_count > 0)
# Tcred <- Tc %>% filter(na_count > 0)
# Gcred <- Gc %>% filter(na_count > 0)
# Ccred <- Cc %>% filter(na_count > 0)
# summary(Acred)
# sd(Acred$na_count)
