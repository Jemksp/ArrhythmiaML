#Set a working directory
out_dir = paste0("Directory Containing Ca2+ processed excel file")
dir.create(out_dir, showWarnings = F)
setwd(out_dir)

##################### Preloading all libraries required ####
library(caret)
library(caretEnsemble)
library(doSNOW)
library(tictoc)

library(openxlsx)
library(ggplot2)
library(tidyr)
library(dplyr)
library(cowplot)

library(Rmisc)
library(caroline)
library(tibble)
library(magrittr)

library(grid)
library(gridExtra)
library(Matrix)
library(scatterplot3d)
library(rgl)
library(plotly)
library(gridExtra)
library(gtable)
library(lattice)
library(stringr)
library(ggpubr)

library(e1071)
library(MLeval)

############################# Preprocessing ####
#Loading in data
data = read.csv(file.choose())

#Setting Title of Project
title = "Full analysis"

#Dropping indexes
# data <- data[(data$Index == 2),] #Keep only this index
data <- data[(data$Index == 1 | data$Index == 3),] #Keep these indexes

#Relabelling columns
colnames(data)[1] = "Sample"
colnames(data)[2] = "Class"
colnames(data)[4] = "BeatToBeat_M"
colnames(data)[5] = "BeatToBeat_V"

#Dropping unimportant parameters
data <- data[c(1:38)]

#Dropping fluorescent parameters
# data <- data[c(1:31)]

#Dropping incomplete samples
data <- data[complete.cases(data),]

#Setting Arrhythmic Classifications as Factors (0 = Not Arrhythmic, 1 = Arrhythmic)
for (factors in 1:length(data$Class)) {
  if (data$Class[factors] == 0) {
    data$Class[factors] = "No"
  } else {
    data$Class[factors] = "Yes"
  }
}

data$Class <- as.factor(data$Class)
data$Class <- factor(data$Class, levels = c("Yes","No"))

#data partition
set.seed(1234)
trainIndex <- createDataPartition(data$Class, p = 0.8, #Change this number to alter training dataset size
                                  list = FALSE,
                                  times = 1)

#Create Train Test sets and drop labels
dataTrain <- data[-3][-1][trainIndex,]
dataTest <- data[-3][-1][-trainIndex,]

################Training ML classifiers and outputing ROC and VarImp plots######
######Training Base learners
#Set-up
set.seed(147)
seeds <- vector(mode = "list", length = 31)
for(i in 1:30) {
  seeds[[i]] <- sample.int(5000, 5000)
}
set.seed(2)
seeds[[31]] <- sample.int(1000,1)

train.control <- trainControl(method = "repeatedcv",
                              number = 10,
                              repeats = 3,
                              savePredictions = "final",
                              index = createFolds(dataTrain$Class, 10),
                              classProbs = TRUE,
                              allowParallel = T,
                              verboseIter = TRUE,
                              seeds = seeds,
                              summaryFunction = twoClassSummary)

algorithmList = c('nnet','rf','xgbTree','svmRadial')
tic()

#Number of Cores to use
cl <- makeCluster(7, type = "SOCK") #Change this number according to the number of cores to be used for training

#Training the learners
registerDoSNOW(cl)
models <- caretList(x = dataTrain[-1],y = dataTrain$Class,
                    trControl = train.control,
                    preProcess = c("center","scale"),
                    metric = "ROC",
                    tuneLength = 8,
                    methodList = algorithmList)
stopCluster(cl)

###### Ensemble Stacking
fullpred <- function(ensemble, test = dataTest) {
  a = list()
  namevector = list()
  results <- resamples(ensemble)
  for (i in 1:(length(ensemble))) {
    test.pred <- predict(ensemble[[i]], test[-1])
    test.cm <- confusionMatrix(test.pred, test$Class)
    a[[i]] <- c(test.cm$overall[1:2],test.cm$byClass[1:6])
    namevector[[i]] <- ensemble[[i]]$method
  }
  names(a) = namevector
  df <- data.frame(a)
  out_list <- list(df, summary(results), dotplot(results), modelCor(results))
  return(out_list)
}
ensemble <- function(fullpred, models, trControl = train.control, newname = "Ensemble", seed = 1234, test = dataTest) {
  set.seed(seed)
  ensemble <- caretEnsemble(models, metric = "ROC", trControl = trControl)
  add.pred <- predict(ensemble,test[-1])
  add <- confusionMatrix(add.pred, test$Class)
  df <- fullpred[[1]]
  df$ens <- c(add$overall[1:2],add$byClass[1:6])
  names(df)[names(df) == "ens"] <- newname
  a = list(df, ensemble, fullpred[[2]], fullpred[[3]], fullpred[[4]])
  return(a)
}

train.control <- trainControl(method = "repeatedcv",
                              number = 10,
                              repeats = 3,
                              savePredictions = TRUE,
                              classProbs = TRUE,
                              verboseIter = TRUE,
                              summaryFunction = twoClassSummary)



cl <- makeCluster(1, type = "SOCK")
registerDoSNOW(cl)

models.backup <- models
models <- models.backup
set.seed(159)
c <- ensemble(fullpred(models), models, newname = "FullEnsemble")
saveRDS(c, file = "model.RDS") #Saving the trained model

#Computing variable importances
vImp <- varImp(c[[2]])

#barplot of variable importance
a <- data.frame(Features = factor(rownames(vImp), levels = rownames(vImp)), Importance = vImp$overall)
out <- ggplot(data = tail(a,10), aes(x = Features)) +
  geom_bar(aes(y = Importance), stat = "identity", colour = "black", size = 1) +
  theme_classic()+
  theme(axis.line = element_line(size = 1, colour = "black"),
        axis.text.y = element_text(size = 18, colour = "black"),
        axis.text.x = element_text(size = 18, colour = "black"),
        axis.title = element_text(size = 18),
        plot.title = element_text(size = 24),
        plot.margin = margin(r = 20),
        text = element_text(face = "bold"),
        axis.ticks.length = unit(.15, "cm"),
        axis.ticks = element_line(size = 1, colour = "black")) +
  labs(title = title) +
  coord_flip()

#ROC Plots
rocscreen <- function(ensemble, title) {
  cols = c("red","purple","springgreen4","darkorchid4","blue")
  screen.nnet <- ensemble[[2]]$models$nnet
  screen.rf <- ensemble[[2]]$models$rf
  screen.xgboost <- ensemble[[2]]$models$xgbTree
  screen.svm <- ensemble[[2]]$models$svmRadial
  screen.ens <- ensemble[[2]]$ens_model
  df.nnet <- data.frame(No = screen.nnet$pred$No, Yes = screen.nnet$pred$Yes, obs = screen.nnet$pred$obs, Group = "nnet")
  df.rf <- data.frame(No = screen.rf$pred$No, Yes = screen.rf$pred$Yes, obs = screen.rf$pred$obs, Group = "rf")
  df.xgboost <- data.frame(No = screen.xgboost$pred$No, Yes = screen.xgboost$pred$Yes, obs = screen.xgboost$pred$obs, Group = "xgbTree")
  df.svm <- data.frame(No = screen.svm$pred$No, Yes = screen.svm$pred$Yes, obs = screen.svm$pred$obs, Group = "svmRadial")
  df.ens <- data.frame(No = screen.ens$pred$No, Yes = screen.ens$pred$Yes, obs = screen.ens$pred$obs, Group = "Stacked")
  df.compiled <- rbind(df.nnet,df.rf,df.xgboost,df.svm,df.ens)
  out <- evalm(df.compiled, fsize = 20, rlinethick = 2, cols = cols, title = title)
  return(out)
}

out2 <- rocscreen(ensemble = c, title = title)

#Output
out_list <- marrangeGrob(list(out, out2$roc), 
             layout_matrix = matrix(1:1, nrow = 1, ncol = 1, byrow = TRUE),
             top = textGrob("", gp = gpar(fontsize = 10)))
ggsave("Summary Combined 3.pdf", out_list,
       width = 10, height = 7, limitsize = F)
toc()



#################Calculate ROC, Sens, Spec CI from training data####################
rocci <- function(mod) {
  a = data.frame()
  a <- rbind(CI(mod[[2]]$models$nnet$resample$ROC),
             CI(mod[[2]]$models$nnet$resample$Sens),
             CI(mod[[2]]$models$nnet$resample$Spec),
             CI(mod[[2]]$models$rf$resample$ROC),
             CI(mod[[2]]$models$rf$resample$Sens),
             CI(mod[[2]]$models$rf$resample$Spec),
             CI(mod[[2]]$models$xgbTree$resample$ROC),
             CI(mod[[2]]$models$xgbTree$resample$Sens),
             CI(mod[[2]]$models$xgbTree$resample$Spec),
             CI(mod[[2]]$models$svmRadial$resample$ROC),
             CI(mod[[2]]$models$svmRadial$resample$Sens),
             CI(mod[[2]]$models$svmRadial$resample$Spec),
             CI(mod[[2]]$ens_model$resample$ROC),
             CI(mod[[2]]$ens_model$resample$Sens),
             CI(mod[[2]]$ens_model$resample$Spec))
  return(a)
}

b <- rocci(c)

b

################Computing Accuracy, Sensitivity and Specificity of test data#########

c[[1]]

############################## Testing a new dataset #############
#Read in the model
c <- readRDS(file.choose())

#Loading in data
data = read.csv(file.choose())

#Dropping indexes
data <- data[(data$Index == 3),] #Keep only this index
# data <- data[(data$Index == 2 | data$Index == 3),] #Keep these indexes

#Relabelling columns
colnames(data)[1] = "Sample"
colnames(data)[2] = "Class"
colnames(data)[4] = "BeatToBeat_M"
colnames(data)[5] = "BeatToBeat_V"

#Dropping unimportant parameters
# data <- data[c(1:38)]

#Dropping fluorescent parameters
data <- data[c(1:31)]

#Dropping incomplete samples
data <- data[complete.cases(data),]

#Setting Arrhythmic Classifications as Factors (0 = Not Arrhythmic, 1 = Arrhythmic)
for (factors in 1:length(data$Class)) {
  if (data$Class[factors] == 0) {
    data$Class[factors] = "No"
  } else {
    data$Class[factors] = "Yes"
  }
}

data$Class <- as.factor(data$Class)
data$Class <- factor(data$Class, levels = c("Yes","No"))

# #Predicting new data
dataTest <- data[-3][-1]

#Test a new dataset
test.pred <- predict(c[[2]], dataTest[-1])
test.cm <- confusionMatrix(test.pred, dataTest$Class)
test.cm
