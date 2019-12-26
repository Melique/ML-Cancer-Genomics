### Machine Learning for Cancer Genomics 
##  K-NN & Random Forest Methods with Leukemia Cancer 
## Bioconductor Installation & Data Installation 

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# Acute Lymphoblastic Leukemia Data installation 
BiocManager::install("ALL")
library(ALL)
data(ALL)

# Segregating the B-Cell Cancer Type 
B_Cell <- grep("^B",as.character(ALL$BT))

# Philadelphia Chromsome Fusion 
gene_Phila <- which(as.character(ALL$mol.biol) %in% c("BCR/ABL","NEG")) #comapre cancer and non-cancer
ALL_Phila <- ALL[,intersect(gene_Phila,B_Cell)]

# Following Step ONE to filter gene with little variation
BiocManager::install("genefilter")
library(genefilter)


ALL_Phila <- nsFilter(ALL_Phila,var.cutoff = 0.99)$eset #extracting gene expression profiles

# STEP 2 , Features Standardization 
# Define the function to measure the IQR for each raw 

rowIQR <- function(eSet) {
  numSamp <- ncol(eSet)
  lowQ <- rowQ(eSet,floor(0.25*numSamp))
  upQ <- rowQ(eSet,ceiling(0.75*numSamp))
  upQ-lowQ
}

standardize <- function(x) (x-rowMedians(x))/rowIQR(x)
ALL_Phila_exprs <- standardize(exprs(ALL_Phila)) #extracts all the gene expression profiles

#STEP 3 , Selecting a Distance & Similarity Metrics 

# Euclidean Distance 
allD <- dist(exprs(ALL_Phila))
allM <- as.matrix(allD)

eucD<- dist(exprs(ALL_Phila),method = "euclidean")
eucMat <- as.matrix(eucD)

manD<- dist(exprs(ALL_Phila),method = "manhattan")
manMat <- as.matrix(manD)

# Color Paletter for the distance visualization 
install.packages("RColorBrewer")
library(RColorBrewer)

install.packages("gplots")
library("gplots")


hmcol <- colorRampPalette(brewer.pal(10,"RdBu"))(256)
heatmap.2(eucMat,main="Euclidean",tracecol=NA,
          symm = TRUE,col=hmcol,
          distfun = as.dist)
heatmap.2(manMat,main="Manhattan",tracecol=NA,
          symm = TRUE,col=hmcol,
          distfun = as.dist)


# Localizing the genes of interets in the proximity 
BiocManager::install("bioDist")
library(bioDist)

### Supervised Machine Learning 

NEGS <- which(ALL_Phila$mol.biol=="NEG")
BCR <- which(ALL_Phila$mol.biol=="BCR/ABL")
sample_1 <- sample(NEGS,22,replace = FALSE)
sample_2 <- sample(BCR,22, replace = FALSE)
Training_set <- c(sample_1,sample_2)
Test_Set <- setdiff(1:80,Training_set)

# Machine Learning using the "MLInterfaces" 
BiocManager::install("MLInterfaces")
library(MLInterfaces)

# KNN Algorithm 
#mol.biol is the gene reagrermanet 
knn1 <- MLearn(formula =mol.biol~.,
               data=ALL_Phila,
               .method = knnI(k=1,l=0),
               trainInd = Training_set)
confuMat(knn1)
error_rate_knn <- round(100-sum(diag(confuMat(knn1)))/sum(confuMat(knn1))*100,2)


# Diagnol Linear Descriminant Analysis Algorithm 
dlda1 <- MLearn(formula=mol.biol~.,
                data = ALL_Phila,
                .method = dldaI,
                trainInd = Training_set)
confuMat(dlda1)
error_rate_dlda1<- round(100-sum(diag(confuMat(dlda1)))/sum(confuMat(dlda1))*100,2)


### Performance Improvements & Tuning 
##  Manual / Statistical Methods 
#   T-Statistics Method on the training set 

training_ttets <- rowttests(ALL_Phila[,Training_set],"mol.biol")
order_training_ttest <- order(abs(training_ttets), decreasing = TRUE)
training_selected <- featureNames(ALL_Phila)[order_training_ttest[1:20]]

# KNN with the reduced genes 

knn2 <- MLearn(formula = mol.biol~.,
               data = ALL_Phila[training_selected,],
               .method = knnI(k=1,l=0),
               trainInd = Training_set)
confuMat(knn2)
error_rate_knn2 <- round(100-sum(diag(confuMat(knn2)))/sum(confuMat(knn2))*100,2)

# Diagnal Linear Desriminant Analysis 

dlda2 <- MLearn(formula=mol.biol~.,
                data = ALL_Phila[training_selected,],
                .method = dldaI,
                trainInd = Training_set)
confuMat(dlda2)
error_rate_dlda2 <- round(100-sum(diag(confuMat(dlda2)))/sum(confuMat(dlda2))*100,2)

### Cross Validation 
##  Cross Validation applying xvalSpec
#   KNN Algorithm 

knn3 <- MLearn(formula =mol.biol~.,
               data=ALL_Phila,
               .method = knnI(k=1,l=0),
               xvalSpec("LOO"))

confuMat(knn3)
error_rate_knn3 <- round(100-sum(diag(confuMat(knn3)))/sum(confuMat(knn3))*100,2)

# KNN with the feature selection using cross validation 

knn4 <- MLearn(formula =mol.biol~.,
               data=ALL_Phila,
               .method = knnI(k=1,l=0),
               xvalSpec("LOO",fsFun =fs.absT(20))) #most revlevant genes that cause cancer
confuMat(knn4)
error_rate_knn4 <- round(100-sum(diag(confuMat(knn4)))/sum(confuMat(knn4))*100,2)

table(unlist(fsHistory(knn4)))
write.csv(table(unlist(fsHistory(knn4))), "top_20.csv")

# KNN with the feature selection with lower selected features 
knn5 <- MLearn(formula =mol.biol~.,
               data=ALL_Phila,
               .method = knnI(k=1,l=0),
               xvalSpec("LOO",fsFun =fs.absT(5)))
confuMat(knn5)
error_rate_knn5 <- round(100-sum(diag(confuMat(knn5)))/sum(confuMat(knn4))*100,2)
write.csv(table(unlist(fsHistory(knn5))), "top_5.csv")

### Ensemble Methods & Random Forest Algorithms 
#   Install Package the random Forest 
install.packages("randomForest")
library(randomForest)

Rf1 <- MLearn(mol.biol~.,
              ALL_Phila,
              randomForestI,
              Training_set,
              ntree=1200,
              mtry=50,
              importance=TRUE) #relevent genes
confuMat(Rf1,"test")
error_rate_knn6 <- round(100-sum(diag(confuMat(Rf1,"test")))/sum(confuMat(Rf1,"test"))*100,2)
 
Rf2 <- MLearn(mol.biol~.,
              ALL_Phila,
              randomForestI,
              Training_set,
              ntree=1200,
              mtry=5,
              importance=TRUE)

confuMat(Rf2,"test")
error_rate_knn7 <- round(100-sum(diag(confuMat(Rf2,"test")))/sum(confuMat(Rf2,"test"))*100,2)

# Feature Selection using Random Forest 
BiocManager::install("hgu95av2.db")
library(hgu95av2.db)

BiocManager::install("org.Hs.eg.db")
library(org.Hs.eg.db)

getVarImp(Rf1) #returns relative importance of different variables in the models using
par(mar=c(5,6,4,1)+.1)
plot(getVarImp(Rf1),main="Rf1 Predictors",n=15,plat="hgu95av2",toktype="SYMBOL",las=1)

getVarImp(Rf2)
par(mar=c(5,6,4,1)+.1)
plot(getVarImp(Rf2),main="Rf2 Predictors",n=15,plat="hgu95av2",toktype="SYMBOL", las=1)
MeanDecreaseGini()
