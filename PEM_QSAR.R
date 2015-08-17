### obtaining data from live Google Sheet (LIVE Data Streaming)
library(googlesheets)
(my_sheets <- gs_ls())
my_sheets %>% glimpse()
data <- gs_title("PEM_data.xlsx")
data_parsed <- data %>% gs_read(ws = "Sheet1")
smiles <- data_parsed$smile
Activity <- data.frame(data_parsed$type)
colnames(Activity) <- c("Activity")
des_quantum <- data_parsed[,c("MeanAbs", "Energy", "Dipole", "HOMO", 
                             "LUMO", "GAP", "EA", "IP", "Electronegativity", 
                             "Hardness", "Electrophilicity", "Softness", "EI")]
library(rcdk)
# Prased SMILES 
mols <- parse.smiles(smiles)
library(Rcpi)
# Molecular Descriptors
des <- suppressWarnings(cbind(
  extractDrugALOGP(mols),
  extractDrugApol(mols),
  extractDrugECI(mols),
  extractDrugTPSA(mols),
  extractDrugWeight(mols),
  extractDrugZagrebIndex(mols),
  extractDrugWienerNumbers(mols),
  extractDrugXLogP(mols),
  extractDrugWienerNumbers(mols),
  extractDrugWHIM(mols),
  extractDrugVAdjMa(mols),
  extractDrugVABC(mols),
  extractDrugRuleOfFive(mols),
  extractDrugRotatableBondsCount(mols),
  extractDrugPetitjeanShapeIndex(mols),
  extractDrugPetitjeanNumber(mols),
  extractDrugMomentOfInertia(mols),
  extractDrugMDE(mols),
  extractDrugMannholdLogP(mols),
  extractDrugLongestAliphaticChain(mols),
  extractDrugLengthOverBreadth(mols),
  extractDrugLargestPiSystem(mols),
  extractDrugLargestChain(mols),
  extractDrugKierHallSmarts(mols),
  extractDrugKappaShapeIndices(mols),
  extractDrugIPMolecularLearning(mols),
  extractDrugHybridizationRatio(mols),
  extractDrugHBondDonorCount(mols),
  extractDrugHBondAcceptorCount(mols),
  extractDrugGravitationalIndex(mols),
  extractDrugFragmentComplexity(mols),
  extractDrugFMF(mols),
  extractDrugCPSA(mols),
  extractDrugChiPathCluster(mols),
  extractDrugChiPath(mols),
  extractDrugChiCluster(mols),
  extractDrugChiChain(mols), 
  extractDrugCarbonTypes(mols), 
  extractDrugBPol(mols), 
  extractDrugBondCount(mols), 
  extractDrugBCUT(mols),
  extractDrugAutocorrelationPolarizability(mols), 
  extractDrugAutocorrelationMass(mols), 
  extractDrugAutocorrelationCharge(mols), 
  extractDrugAtomCount(mols), 
  extractDrugAromaticBondsCount(mols), 
  extractDrugAromaticAtomsCount(mols), 
  extractDrugAminoAcidCount(mols)
))

des_mol <- des[, colSums(is.na(des)) != nrow(des)]
des_mol = des_mol[, -nearZeroVar(des_mol)]
removed_cor <- function(x) {
  library(caret)
  corelation <- cor(x)
  des_correlation <- corelation[1:dim(corelation)[1], 1:dim(corelation)[1]]
  high_correlation <- findCorrelation(des_correlation, cutoff = .07)
  filtered_des <- x[, -high_correlation]
  return(filtered_des)
}
des_mol_filtered <- removed_cor(des_mol)
des_quantum_filtered <- removed_cor(des_quantum)
data_ready <- cbind(Activity, des_mol, des_quantum_filtered)
### Model Training

J48_training <- function(x) {
  library(caret)
  results <- list(100)
  for (i in 1:100) {
    set.seed(i)
    trainIndex <- createDataPartition(x$Activity, p = .8, list = FALSE, times = 1)
    train <- x[ trainIndex,]
    test <- x[ -trainIndex, ]
    model_train <- J48(Activity~., data = train)
    summary <- summary(model_train)
    confusionmatrix <- summary$confusionMatrix
    results[[i]] <- as.numeric(confusionmatrix)
  }
  return(results)
}

mean_and_sd <- function(x) {
  c(round(mean(x, na.rm = TRUE), digits = 4),
    round(sd(x, na.rm = TRUE), digits = 4))
}

J48_train <- function(x) {
  ok <- J48_training(x)
  results <- data.frame(ok)
  data <- data.frame(results)
  m = ncol(data)
  ACC  <- matrix(nrow = m, ncol = 1)
  SENS  <- matrix(nrow = m, ncol = 1)
  SPEC  <-matrix(nrow = m, ncol = 1)
  MCC <- matrix(nrow = m, ncol = 1)
  
  for(i in 1:m){ 
    ACC[i,1]  = (data[1,i]+data[4,i])/(data[1,i]+data[2,i]+data[3,i]+data[4,i])*100
    SENS[i,1]  =  (data[4,i])/(data[3,i]+data[4,i])*100
    SPEC[i,1]  = (data[1,i]/(data[1,i]+data[2,i]))*100
    MCC1      = (data[1,i]*data[4,i]) - (data[2,i]*data[3,i])
    MCC2      =  (data[4,i]+data[2,i])*(data[4,i]+data[3,i])
    MCC3      =  (data[1,i]+data[2,i])*(data[1,i]+data[3,i])
    MCC4  =  sqrt(MCC2)*sqrt(MCC3)
    
    
    MCC[i,1]  = MCC1/MCC4
  }
  results_ACC <- mean_and_sd(ACC)
  results_SENS <- mean_and_sd(SENS)
  results_SPEC <- mean_and_sd(SPEC)
  results_MCC <- mean_and_sd(MCC)
  return(data.frame(c(results_ACC, results_SENS, results_SPEC, results_MCC)))
}
### 10 fold cross validation
J48_10fold <- function(x) {
  results <- list(100)
  for (i in 1:100) {
    trainIndex <- createDataPartition(x$Activity, p = .8, list = FALSE, times = 1)
    train <- x[ trainIndex,]
    test <- x[ -trainIndex, ]
    model_train <- J48(Activity~., data = train)
    eval_j48 <- evaluate_Weka_classifier(model_train, numFolds = 10, complexity = FALSE, seed = 1, class = TRUE)
    confusionmatrix <- eval_j48$confusionMatrix
    results[[i]] <- as.numeric(confusionmatrix)
  }
  return(results)
}

J48_cross_validation <- function(x) {
  ok <- J48_10fold(x)
  results <- data.frame(ok)
  data <- data.frame(results)
  m = ncol(data)
  ACC  <- matrix(nrow = m, ncol = 1)
  SENS  <- matrix(nrow = m, ncol = 1)
  SPEC  <-matrix(nrow = m, ncol = 1)
  MCC <- matrix(nrow = m, ncol = 1)
  
  for(i in 1:m){ 
    ACC[i,1]  = (data[1,i]+data[4,i])/(data[1,i]+data[2,i]+data[3,i]+data[4,i])*100
    SENS[i,1]  =  (data[4,i])/(data[3,i]+data[4,i])*100
    SPEC[i,1]  = (data[1,i]/(data[1,i]+data[2,i]))*100
    MCC1      = (data[1,i]*data[4,i]) - (data[2,i]*data[3,i])
    MCC2      =  (data[4,i]+data[2,i])*(data[4,i]+data[3,i])
    MCC3      =  (data[1,i]+data[2,i])*(data[1,i]+data[3,i])
    MCC4  =  sqrt(MCC2)*sqrt(MCC3)
    
    
    MCC[i,1]  = MCC1/MCC4
  }
  results_ACC <- mean_and_sd(ACC)
  results_SENS <- mean_and_sd(SENS)
  results_SPEC <- mean_and_sd(SPEC)
  results_MCC <- mean_and_sd(MCC)
  return(data.frame(c(results_ACC, results_SENS, results_SPEC, results_MCC)))
}
### Testing Results 

J48_testing <- function(x) {
  results <- list(100)
  for (i in 1:100) {
    trainIndex <- createDataPartition(x$Activity, p = .8, list = FALSE, times = 1)
    train <- x[ trainIndex,]
    test <- x[ -trainIndex, ]
    model_train <- J48(Activity~., data = train)
    eval_external <- evaluate_Weka_classifier(model_train, newdata = test, numFolds = 0, complexity = FALSE, seed = 1, class = TRUE)
    confusionmatrix <- eval_external$confusionMatrix
    results[[i]] <- as.numeric(confusionmatrix)
  }
  return(results)
}

J48_external <- function(x) {
  ok <- J48_testing(x)
  results <- data.frame(ok)
  data <- data.frame(results)
  m = ncol(data)
  ACC  <- matrix(nrow = m, ncol = 1)
  SENS  <- matrix(nrow = m, ncol = 1)
  SPEC  <-matrix(nrow = m, ncol = 1)
  MCC <- matrix(nrow = m, ncol = 1)
  
  for(i in 1:m){ 
    ACC[i,1]  = (data[1,i]+data[4,i])/(data[1,i]+data[2,i]+data[3,i]+data[4,i])*100
    SENS[i,1]  =  (data[4,i])/(data[3,i]+data[4,i])*100
    SPEC[i,1]  = (data[1,i]/(data[1,i]+data[2,i]))*100
    MCC1      = (data[1,i]*data[4,i]) - (data[2,i]*data[3,i])
    MCC2      =  (data[4,i]+data[2,i])*(data[4,i]+data[3,i])
    MCC3      =  (data[1,i]+data[2,i])*(data[1,i]+data[3,i])
    MCC4  =  sqrt(MCC2)*sqrt(MCC3)
    
    
    MCC[i,1]  = MCC1/MCC4
  }
  results_ACC <- mean_and_sd(ACC)
  results_SENS <- mean_and_sd(SENS)
  results_SPEC <- mean_and_sd(SPEC)
  results_MCC <- mean_and_sd(MCC)
  return(data.frame(c(results_ACC, results_SENS, results_SPEC, results_MCC)))
}




J48_train_results <- J48_train(data_ready)
J48_10_fold_results <- J48_cross_validation(data_ready)
J48_test_results <- J48_external(data_ready)

### PCA Analysis 
data_pca <- data_ready[-1]

df.sc <- scale(data_pca)
df.svd <- svd(df.sc)
df.scores <- df.svd$u %*% diag(df.svd$d)
df.loadings <- df.svd$v
df.vars <- df.svd$d^2 / (nrow(df) -1)
df.totalvar <- sum(df.vars)
df.relvars <- df.vars / df.totalvar
variances <- 100 * round(df.relvars, digits = 3)
variances[1:5]
par(mfrow = c(2, 2))
barplot(df.vars[1:10], main = "Variances",
        names.arg = paste("PC", 1:10))
barplot(log(variances[1:10]), main = "log(Variances)",
        names.arg = paste("PC", 1:10))
barplot(100*df.relvars[1:10], main = "Relative variances (%)",
        names.arg = paste("PC", 1:10))
barplot(cumsum(100*df.relvars[1:10]),
        main = "Cumulative variances (%)",
        names.arg = paste("PC", 1:10), ylim = c(0, 100))
### Feature Importance
# C5.0 Algorithm is used to obtain feature importance of each AAC, DPC and PCP
set.seed(34440)
Importance <- lapply(data_ready, function(x) {
  library(C50)
  Model <- C5.0(Activity~., data = data_ready, rules=TRUE)
  Importance <- C5imp(Model)
  return(Importance)
})
set.seed(333)
features <- data.frame(Importance)
importance <- cbind( Descriptors = rownames(features), features[1])
set.seed(1)
top3Importance <- head(Importance,3)




