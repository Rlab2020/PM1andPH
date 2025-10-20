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
### COX FOR DEATH - CONTINUOUS
death_labels <- c(
  "death"    = "all-cause",
  "deathAB"  = "infectious",
  "deathC"   = "neoplasms",
  "deathD"   = "blood",
  "deathE"   = "endocrine",
  "deathF"   = "mental",
  "deathG"   = "neurological",
  "deathI"   = "circulatory",
  "deathJ"   = "respiratory",
  "deathK"   = "digestive",
  "deathM"   = "musculoskeletal",
  "deathN"   = "genitourinary")

death_vars <- c("death","deathAB","deathC","deathD","deathE","deathF","deathG","deathI","deathJ","deathK","deathM","deathN")

table2c <- matrix(NA, nrow = 12, ncol = 7)
colnames(table2c) <- c("ALL","AGE1","AGE2","P_AGE","MALE","FEMALE","P_SEX")
rownames(table2c) <- unname(death_labels)

for (i in seq_along(death_vars)) {
  var <- death_vars[i]
  table2c[i, ] <- c(FUNa(var, "death_time", "mean0610b1000i", "data0"),FUNa1(var, "death_time", "mean0610b1000i", "data0"),FUNa2(var, "death_time", "mean0610b1000i", "data0"))}

table2c <- as.data.frame(table2c, stringsAsFactors = FALSE)
table2c$var <- row.names(table2c)
table2c <- table2c %>% select(var, everything())
write.csv(table2c, file = "/Users/guuuuu/Desktop/table2c.csv", row.names = FALSE)


#######################################################################
###--------------------------------------------------------------------
### COX FOR DEATH - TERTILE
table2d <- matrix(NA, nrow = 12, ncol = 3)
colnames(table2d) <- c("total_middle","total_high","p for trend")
rownames(table2d) <- unname(death_labels)

for (i in seq_along(death_vars)) {
  var <- death_vars[i]
  table2d[i, ] <- FUNb(var, "death_time", "data0")
}

table2d <- as.data.frame(table2d, stringsAsFactors = FALSE)
table2d$var <- row.names(table2d)
table2d <- table2d %>% select(var, everything())
write.csv(table2d, file = "/Users/guuuuu/Desktop/table2d.csv", row.names = FALSE)

