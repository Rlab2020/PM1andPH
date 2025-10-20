#######################################################################
###--------------------------------------------------------------------
### LOADING PACKAGES ###
setwd("/Users/guuuuu/Desktop")

packages=c('reshape2','dplyr','openair','zoo','car','foreign','stats','data.table',
           'stringr','scales','readxl','writexl','compareGroups','survival','survminer',
           'lme4','lmerTest','car','glmnet','broom','caret','rms','xgboost','impute','mediation')
# lapply(packages, install.packages, character.only=T)
lapply(packages, require, character.only=T)
# Appoint Func
select <- dplyr::select
rename <- dplyr::rename


#######################################################################
###--------------------------------------------------------------------
### LOADING DATA ###
load("/Users/guuuuu/Desktop/data_pm01_m.rda")


#######################################################################
###--------------------------------------------------------------------
### LINEAR REGRESSION WITH FDR

FUN4 <- function(x, y, dataset){
  eval(parse(text=paste("model1 <- glm(",y,"~",x,"+age+sex+ses+ethnic+bmi+whr+nosmoking+drinkingM+sleeph+regular+dieth+occup+fsp0610+ndvi0610+lden+fastt+sampleage+season, family='gaussian', ",dataset,")",sep='')))
  re1a=summary(model1)$coefficients;re1b=confint(model1)
  result1 = re1a[2,1]; ci1 = re1b[2,1]; ci2 = re1b[2,2]; p1=re1a[2,4]
  result = cbind(result1,ci1,ci2,p1); result}

variables <- c("mean0610b1000")

table4a <- matrix(NA, nrow = length(nmrlist), ncol = 4*length(variables))
colnames(table4a) <- c("coef1","ci1","ci2","p1")
rownames(table4a) <- nmrlist

for (i in seq_along(variables)) {
  variable_i <- variables[i]
  for (j in seq_along(nmrlist)) {
    nmrlist_j <- nmrlist[j]
    table4a[j, (4 * i - 3):(4 * i)] <- FUN4(variable_i,nmrlist_j,"data0")
    print(paste0("Done ",variables[i]," ",nmrlist[j]))}}

table4a <- as.data.frame(table4a, stringsAsFactors = FALSE)
table4a$var <- row.names(table4a)  
table4a <- table4a %>% select(5,1:4)

names(table4a) <- c("name","coef1","ci1","ci2","p1")
table4a$adjusted_p <- p.adjust(table4a$p1, method = "BH")
table4a <- table4a %>% filter(adjusted_p<=0.05)

library(writexl)
write_xlsx(table4a, path = "/Users/guuuuu/Desktop/TABLE - PM01-NMR BY LINEAR.xlsx")


#######################################################################
###--------------------------------------------------------------------
### ELASTIC NET REGRESSION
table4a <- read_excel("/Users/guuuuu/Desktop/TABLE - PM01-NMR BY LINEAR.xlsx")

set.seed(123)
index <- sample(1:nrow(data0), 0.7 * nrow(data0))
table(data0$death[index])
table(data0$death[-index])

data0_training <- data0[index, ]
data0_testing <- data0[-index, ]

data0$subset <- "train"
data0[-index, ]$subset <- "test"
data0$subset <- factor(data0$subset, levels = c("train", "test"))

covariates <- c("age","sex","ses","bmi","whr","ethnic","nosmoking","drinkingM","sleeph","regular","dieth","occup","fsp0610","ndvi0610","lden","fastt","sampleage","season")
nmrlist1 <- table4a$name

X <- data0_training %>% select(all_of(nmrlist1),all_of(covariates))
X_matrix <- makeX(X)
Y <- data0_training[["mean0610b1000"]]

set.seed(123)  
cv_model <- cv.glmnet(X_matrix, Y, alpha = 0.5, nfolds = 10, 
                      penalty.factor = ifelse(colnames(X_matrix) %in% nmrlist1, 1, 0))

lambda_min <- cv_model$lambda.min  
lambda_1se <- cv_model$lambda.1se 

elastic_net_model <- glmnet(X_matrix, Y, alpha = 0.5, lambda = lambda_min, 
                            penalty.factor = ifelse(colnames(X_matrix) %in% nmrlist1, 1, 0))

table4b <- as.data.frame(as.matrix(coef(elastic_net_model)))
table4b$name <- rownames(table4b)
table4b <- table4b[table4b$s0 != 0, ] 

rownames(table4b) <- NULL 
table4b <- table4b %>%
  filter(name != "(Intercept)") %>%
  filter(name %in% nmrlist1) %>%
  rename(coef = s0) %>%
  select(name, coef)

library(writexl)
write_xlsx(table4b, path = "/Users/guuuuu/Desktop/TABLE - PM01-NMR BY ENR.xlsx")


#######################################################################
###--------------------------------------------------------------------
### LASSO REGRESSION
table4a <- read_excel("/Users/guuuuu/Desktop/TABLE - PM01-NMR BY LINEAR.xlsx")

set.seed(123)
index <- sample(1:nrow(data0), 0.7 * nrow(data0))
table(data0$death[index])
table(data0$death[-index])

data0_training <- data0[index, ]
data0_testing <- data0[-index, ]

data0$subset <- "train"
data0[-index, ]$subset <- "test"
data0$subset <- factor(data0$subset, levels = c("train", "test"))

covariates <- c("age","sex","ses","bmi","whr","ethnic","nosmoking","drinkingM","sleeph","regular","dieth","occup","fsp0610","ndvi0610","lden","fastt","sampleage","season")
nmrlist1 <- table4a$name

X <- data0_training %>% select(all_of(nmrlist1),all_of(covariates))
X_matrix <- makeX(X)
Y <- data0_training[["mean0610b1000"]]

set.seed(123)  
cv_model <- cv.glmnet(X_matrix, Y, alpha = 1, nfolds = 10, 
                      penalty.factor = ifelse(colnames(X_matrix) %in% nmrlist1, 1, 0))

lambda_min <- cv_model$lambda.min  
lambda_1se <- cv_model$lambda.1se 

lasso_model <- glmnet(X_matrix, Y, alpha = 1, lambda = lambda_min, 
                            penalty.factor = ifelse(colnames(X_matrix) %in% nmrlist1, 1, 0))

table4b <- as.data.frame(as.matrix(coef(lasso_model)))
table4b$name <- rownames(table4b)
table4b <- table4b[table4b$s0 != 0, ] 

rownames(table4b) <- NULL 
table4b <- table4b %>%
  filter(name != "(Intercept)") %>%
  filter(name %in% nmrlist1) %>%
  rename(coef = s0) %>%
  select(name, coef)

library(writexl)
write_xlsx(table4b, path = "/Users/guuuuu/Desktop/TABLE - PM01-NMR BY LR.xlsx")
