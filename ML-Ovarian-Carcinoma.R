if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("curatedOvarianData")
library(curatedOvarianData)
library(tidyverse)
data(package="curatedOvarianData")
data(TCGA_eset)
mat <- exprs(TCGA_eset)
matT <- as.data.frame(t(mat))
[1]   578 13104
grade <- TCGA_eset$grade
stage <- TCGA_eset$tumorstage
substage <- TCGA_eset$substage
recurrence <- TCGA_eset$recurrence_status
normal_cells <- TCGA_eset$percent_normal_cells
tumor_cells <- TCGA_eset$percent_tumor_cells

# Predict grade/stage/etc in patients using NN only with selected 14 sig genes
Final <- matT %>% select(KRAS, BRAF, PIK3CA, ERBB2, CTNNB1, ARID1A, PPP2R1A, TP53, BRCA1, BRCA2, RAD51, PALB2, CHEK2, BARD1)
Final$grade <- grade
Final$stage <- stage
Final$substage <- substage
Final$recurrence <- recurrence
Final$normal_cells <- normal_cells
Final$tumor_cells <- tumor_cells
data <- na.omit(Final) 
[1] 457  20

write.csv(Final,"/Users/yalegenomecenter/Desktop/Expression_Matrix.csv", row.names = FALSE)


####################Neural Net###################################
samplesize = 0.60 * nrow(data)
set.seed(80)
index = sample( seq_len ( nrow ( data ) ), size = samplesize )

datatrain = data[ index, ]
datatest = data[ -index, ]

install.packages("neuralnet")
library(neuralnet)

data_norm<- function(x) {((x-min(x))/(max(x)-min(x)))}
scaled <- as.data.frame(lapply(data[1:16], data_norm)) 
summary(variables_norm)

trainNN = scaled[index , ]
testNN = scaled[-index , ]

set.seed(2)
NN = neuralnet(stage ~ KRAS + BRAF + PIK3CA + ERBB2 + CTNNB1 + ARID1A + PPP2R1A + TP53 + BRCA1 + BRCA2 + RAD51 + PALB2 + CHEK2 + BARD1, trainNN, hidden = 3 , linear.output = T )
plot(NN)

predict_testNN = compute(NN, testNN)
predict_testNN = (predict_testNN$net.result * (max(data$stage) - min(data$stage))) + min(data$stage)

plot(datatest$stage, predict_testNN, col='blue', pch=16, ylab = "Predicted Stage", xlab = "Actual Stage")

abline(0,1)

# Calculate Root Mean Square Error (RMSE)
RMSE.NN = (sum((datatest$stage - predict_testNN)^2) / nrow(datatest)) ^ 0.5

# Load libraries
library(boot)
library(plyr)

# Initialize variables
set.seed(50)
k = 100
RMSE.NN = NULL

List = list( )

# Fit neural network model within nested for loop
for(j in 1:20){
    for (i in 1:k) {
        index = sample(1:nrow(data),j )

        trainNN = scaled[index,]
        testNN = scaled[-index,]
        datatest = data[-index,]

        NN = neuralnet(stage ~ KRAS + BRAF + PIK3CA + ERBB2 + CTNNB1 + ARID1A + PPP2R1A + TP53 + BRCA1 + BRCA2 + RAD51 + PALB2 + CHEK2 + BARD1, trainNN, hidden = 3 , linear.output = T )
        predict_testNN = compute(NN,testNN)
        predict_testNN = (predict_testNN$net.result*(max(data$stage)-min(data$stage)))+min(data$stage)

        RMSE.NN [i]<- (sum((datatest$stage - predict_testNN)^2)/nrow(datatest))^0.5
    }
    List[[j]] = RMSE.NN
}

Matrix.RMSE = do.call(cbind, List)
library(matrixStats)

med = colMedians(Matrix.RMSE)

X = seq(1,20)

plot (med~X, type = "l", xlab = "Length Of Training Set", ylab = "Median RMSE", main = "Variation of RMSE with length of training set")

##########################SVM#############################
library(e1071) 
classifier = svm(formula = grade ~ ., 
                 data = datatrain, 
                 type = 'C-classification', 
                 kernel = 'linear')
y_pred = predict(classifier, newdata = datatest) 
cm = table(datatest$grade, y_pred)


