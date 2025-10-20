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
### LOADING FINAL DATA ###
load("/Users/guuuuu/Desktop/pm01/data/data_pm01_all.rda")


#######################################################################
###--------------------------------------------------------------------
### VARIABLE TRANSFORMATION ###
### AGE: <=59 and >=60
data0 <- data0 %>% mutate(age0 = case_when(age<=59~"level1",T~"level2"))
data0$age0 <- factor(data0$age0, levels=c("level1","level2"))
### SES
data0$ses <- factor(data0$ses, levels = c("1","2","3"), labels = c("high", "moderate", "low"))
### WHR: poor - male+0.9 and female+0.85
data0 <- data0 %>% mutate(whr0 = case_when(is.na(whr)~NA_character_,whr>=0.9&sex=="male"~"poor",whr>=0.85&sex=="female"~"poor",T~"ideal")) 
data0$whr0 <- factor(data0$whr0, levels=c("ideal","poor"))
### BMI: <18.5, 18.5-24.9, 25.0-29.9, >=30.0
data0 <- data0 %>% mutate(bmi0 = case_when(is.na(bmi)~NA_character_,bmi<18.5~"level1", bmi>=18.5&bmi<=24.9~"level2", bmi>=25.0&bmi<=29.9~"level3", T~"level4")) 
data0$bmi0 <- factor(data0$bmi0, levels=c("level1","level2","level3","level4"))
### Ethnic: white, mixed, black and asian
data0$ethnic <- ifelse(data0$ethnic == "white", "white", "nonwhite")
data0$ethnic <- factor(data0$ethnic, levels=c("white","nonwhite"))


#######################################################################
###--------------------------------------------------------------------
### DATA PREPARATION 
data0 <- data0 %>% filter(!is.na(mean0610b1000)) %>% filter(city=="urban") %>% filter(unfollow==0) 
P995 <- quantile(data0$mean0610b1000, 0.995, na.rm = TRUE)
data0 <- data0[data0$mean0610b1000 <= P995, ] ### 425165

data0 <- data0 %>% mutate(
  transport0610 = transmean0610b1000 + othertransmean0610b1000,
  industry0610 = energymean0610b1000 + industmean0610b1000 + processmean0610b1000 + solventsmean0610b1000,
  resident0610 = nonindustmean0610b1000,
  others0610 = wastemean0610b1000 + naturemean0610b1000)

data_london <- data0  %>%  filter(center=="Barts"|center=="Croydon"|center=="Hounslow")
data_mid <- data0  %>%  filter(center=="Birmingham"|center=="Stoke")
data_nw <- data0  %>%  filter(center=="Manchester"|center=="Liverpool"|center=="Bury"|center=="Stockport(pilot)")
data_york <- data0  %>%  filter(center=="Leeds"|center=="Sheffield")
data_ne <- data0  %>%  filter(center=="Newcastle"|center=="Middlesbrough")


#######################################################################
###--------------------------------------------------------------------
### CONTINUOUS VARIABLE ###
data0$mean0610b1000i <- data0$mean0610b1000/IQR(data0$mean0610b1000)


#######################################################################
###--------------------------------------------------------------------
### CATEGORICAL VARIABLE TERTILE ###
data0$mean0610b1000c3 <- cut(data0$mean0610b1000,breaks = quantile(data0$mean0610b1000, probs = c(0, 1/3, 2/3, 1), na.rm = TRUE),labels = c("L1", "L2", "L3"),include.lowest = TRUE)
mean0610b1000m <- quantile(data0$mean0610b1000, c(0, 1/3, 2/3, 1), na.rm = TRUE)
data0$mean0610b1000c30 <- cut(data0$mean0610b1000, breaks = mean0610b1000m, labels = c(rollapplyr(as.numeric(mean0610b1000m), list(-1:0), median, fill = NA)[2:4]), na.rm = TRUE)
data0$mean0610b1000c30 <- as.numeric(as.character(data0$mean0610b1000c30)) ### mean0610b1000c30 is used to calculate P for trend


#######################################################################
###--------------------------------------------------------------------
### TABLE 1
library(compareGroups)
data_table1 <- data0 %>% select(mean0610b1000c3, age, sex, ethnic, ses, whr, bmi, nosmoking, drinkingM, sleeph, regular, dieth, occup, fsp0610, ndvi0610, lden)
table1 <- compareGroups(mean0610b1000c3 ~ ., data = data_table1)
table1 <- createTable(table1, show.all = TRUE, hide.no = "no", show.p.overall = TRUE)
export2xls(table1, file = "table1.xlsx")

#######################################################################
###--------------------------------------------------------------------
### SUB DATASET PREPARATION ###
event_vars <- c("cancer", "dm2", "deme", "mdd", "gad","ihd", "ar", "af", "hf", "cevd", "copd", "cld", "op", "ckd")

for (var in event_vars) {data0[[var]] <- as.numeric(data0[[var]] == 1)
time_col <- paste0(var, "_time")
data0[[time_col]] <- as.numeric(data0[[time_col]])}

subset_list <- lapply(event_vars, function(var){
  pre_col  <- paste0(var, "p")
  time_col <- paste0(var, "_time")
  subset(data0, data0[[pre_col]] == 0 & data0[[time_col]] > 0)})

names(subset_list) <- event_vars

data_cancer <- subset_list[["cancer"]]
data_dm2    <- subset_list[["dm2"]]
data_deme   <- subset_list[["deme"]]
data_mdd    <- subset_list[["mdd"]]
data_gad    <- subset_list[["gad"]]
data_ihd    <- subset_list[["ihd"]]
data_ar     <- subset_list[["ar"]]
data_af     <- subset_list[["af"]]
data_hf     <- subset_list[["hf"]]
data_cevd   <- subset_list[["cevd"]]
data_copd   <- subset_list[["copd"]]
data_cld    <- subset_list[["cld"]]
data_op     <- subset_list[["op"]]
data_ckd    <- subset_list[["ckd"]]

#######################################################################
###--------------------------------------------------------------------
### COX FOR DISEASE - CONTINUOUS

FUNa <- function(x,y,z,dataset){
  eval(parse(text=paste("model1<-coxph(Surv(",y,",",x,")~",z,"+age+sex+ethnic+ses+bmi+whr+nosmoking+drinkingM+sleeph+regular+dieth+fsp0610+ndvi0610+lden+occup,",dataset,")",sep='')))
  re1=summary(model1)$conf.int
  result1 = paste0(round(re1[1,1],2)," (",round(re1[1,3],2),", ",round(re1[1,4],2),")")
  result=cbind(result1);result}

FUNa1 <- function(x,y,z,dataset){
  eval(parse(text=paste("model1<-coxph(Surv(",y,",",x,")~",z,"*relevel(age0, ref='level1')+sex+ethnic+ses+bmi+whr+nosmoking+drinkingM+sleeph+regular+dieth+fsp0610+ndvi0610+lden+occup,",dataset,")",sep='')))
  eval(parse(text=paste("model2<-coxph(Surv(",y,",",x,")~",z,"*relevel(age0, ref='level2')+sex+ethnic+ses+bmi+whr+nosmoking+drinkingM+sleeph+regular+dieth+fsp0610+ndvi0610+lden+occup,",dataset,")",sep='')))
  re1=summary(model1)$conf.int; re2=summary(model2)$conf.int; re3=summary(model1)$coefficient
  result1=paste0(round(re1[1,1],2)," (",round(re1[1,3],2),", ",round(re1[1,4],2),")")
  result2=paste0(round(re2[1,1],2)," (",round(re2[1,3],2),", ",round(re2[1,4],2),")")
  p1=round(re3[nrow(re3),5],3)
  result=c(result1,result2,p1); result}

FUNa2 <- function(x,y,z,dataset){
  eval(parse(text=paste("model1<-coxph(Surv(",y,",",x,")~",z,"*relevel(sex, ref='male')+age+ethnic+ses+bmi+whr+nosmoking+drinkingM+sleeph+regular+dieth+fsp0610+ndvi0610+lden+occup,",dataset,")",sep='')))
  eval(parse(text=paste("model2<-coxph(Surv(",y,",",x,")~",z,"*relevel(sex, ref='female')+age+ethnic+ses+bmi+whr+nosmoking+drinkingM+sleeph+regular+dieth+fsp0610+ndvi0610+lden+occup,",dataset,")",sep='')))
  re1=summary(model1)$conf.int; re2=summary(model2)$conf.int; re3=summary(model1)$coefficient
  result1=paste0(round(re1[1,1],2)," (",round(re1[1,3],2),", ",round(re1[1,4],2),")")
  result2=paste0(round(re2[1,1],2)," (",round(re2[1,3],2),", ",round(re2[1,4],2),")")
  p1=round(re3[nrow(re3),5],3)
  result=c(result1,result2,p1); result}

table2a <- matrix(NA, nrow = length(event_vars), ncol = 7)
colnames(table2a) <- c("ALL","AGE1","AGE2","P_AGE","MALE","FEMALE","P_SEX")
rownames(table2a) <- event_vars

expo_var <- "mean0610b1000i"

for (i in seq_along(event_vars)) {ev <- event_vars[i]
time_col <- paste0(ev, "_time")
data_name <- paste0("data_", ev)
table2a[i, ] <- c(FUNa(ev, time_col, expo_var, data_name),FUNa1(ev, time_col, expo_var, data_name),FUNa2(ev, time_col, expo_var, data_name))}

table2a <- as.data.frame(table2a, stringsAsFactors = FALSE)
table2a$var <- rownames(table2a)
table2a <- table2a %>% select(var, everything())
write.csv(table2a, file = "/Users/guuuuu/Desktop/table2a.csv", row.names = FALSE)


#######################################################################
###--------------------------------------------------------------------
### COX FOR DISEASE - TERTILE

FUNb <- function(x,y,dataset){
  eval(parse(text=paste("model1<-coxph(Surv(",y,",",x,")~mean0610b1000c3+age+sex+ethnic+ses+bmi+whr+nosmoking+drinkingM+sleeph+regular+dieth+fsp0610,",dataset,")",sep='')))
  eval(parse(text=paste("model2<-coxph(Surv(",y,",",x,")~mean0610b1000c30+age+sex+ethnic+ses+bmi+whr+nosmoking+drinkingM+sleeph+regular+dieth+fsp0610,",dataset,")",sep='')))
  re1=summary(model1)$conf.int;re1p=round(summary(model2)$coef[1,5],3)
  result1a = paste0(round(re1[1,1],2)," (",round(re1[1,3],2),", ",round(re1[1,4],2),")")
  result1b = paste0(round(re1[2,1],2)," (",round(re1[2,3],2),", ",round(re1[2,4],2),")")
  result = matrix(c(result1a,result1b,re1p),ncol=3,nrow=1,byrow=T)
  result}

table2b <- t(sapply(event_vars, function(ev){
  time_col <- paste0(ev, "_time")
  data_name <- paste0("data_", ev)
  FUNb(ev, time_col, data_name)}))

table2b <- as.data.frame(table2b, stringsAsFactors = FALSE)
colnames(table2b) <- c("total_middle", "total_high", "p_for_trend")
table2b$var <- rownames(table2b)
table2b <- table2b %>% select(var, everything())
write.csv(table2b, file = "/Users/guuuuu/Desktop/table2b.csv", row.names = FALSE)

