install.packages("ggplot2")
install.packages("psyphy")
library (ggplot2)
library (psyphy)
setwd('/Volumes/emmt/data/beh/')

sub.name <- c("MG", "AW", "KJ", "SP", "LR", "KC", "SW", "AK", "NR", "CH", "KH", "MM", "BL", "SP", "SS",
              "CS", "MR", "JW", "RW", "MT", "SC", "SS", "GH", "HG", "SA", "AU", "EL", "NM", "SL", 
              "CP", "PK", "GS", "MH", "HB", "SR", "JL", "LS", "ME", "LR", "ZA", "PS", "AB", "AB",
              "AR", "BC", "JL", "YB", "SF", "JB", "MH", "JG", "MR")
sub.num <- c("007","008","009","012", "015", "016", "020", "025", "030", "040", "045", "050", "055", "059",
             "060", "061", "062", "063", "058", "064", "065", "067", "070", "071", "072", "073",
             "074", "075", "076", "077", "078", "079", "081", "082", "083", "084", "085", "086", "087", 
             "088", "089", "090", "091", "092", "093", "094", "095", "097", "098", "100", "101", "102")
sub.num.wm <- c("034","039","042","054","056","057","058","060","061","062","063","064","065","066","067",
             "069","071","072","073","074","075","076","077","078","079","081","082","084","086","087","088",
             "090","091","092","093","094","095","096","097","098","100","101","102","103","104","105","106",
             "107","108","109","110","111","112","113","114","116","117","118","119")
total.study.trials <- 6
total.test.trials <- 3

##Importing participant study data into one dataframe (study.sub.data)
study.sub.data <- c("")
for (isub in 1:length(sub.num)){
  for (trial.num in 1:total.study.trials){
    curr.sub <- read.csv(file = paste("s", sub.num[isub], "/", "eMMT",sub.num[isub], "_study", toString(trial.num), "_", sub.name[isub],
                                       "_out.txt", sep = ""), header = TRUE, sep = "\t", quote = "\"",dec = ".")
    study.sub.data <- merge(study.sub.data, curr.sub, all=T)
  }
}

##Importing participant test data into one dataframe (test.sub.data)
test.sub.data <- c("")
for (isub in 1:length(sub.num)){
  for (trial.num in 1:total.test.trials){
    curr.sub <- read.csv(file = paste("s", sub.num[isub], "/", "eMMT",sub.num[isub], "_test", toString(trial.num), ".", sub.name[isub],
                                      ".out.txt", sep = ""), header = TRUE, sep = "\t", quote = "\"",dec = ".")
    test.sub.data <- merge(test.sub.data, curr.sub, all=T)
  }
}
##Adding mem to test data 
test.sub.data$mem <- ""
for(isub in 1:nrow(test.sub.data)){
  if((test.sub.data$OldNew[isub]=="New" && test.sub.data$resp[isub]=="j") || 
      (test.sub.data$OldNew[isub]=="Old" && test.sub.data$resp[isub]=="f")){
    test.sub.data$mem[isub] <- "hi_miss"
  }
  else if((test.sub.data$OldNew[isub]=="New" && test.sub.data$resp[isub]=="k") || 
     (test.sub.data$OldNew[isub]=="Old" && test.sub.data$resp[isub]=="d")){
    test.sub.data$mem[isub] <- "low_miss"
  }
  else if((test.sub.data$OldNew[isub]=="New" && test.sub.data$resp[isub]=="f") || 
          (test.sub.data$OldNew[isub]=="Old" && test.sub.data$resp[isub]=="j")){
    test.sub.data$mem[isub] <- "hi_hit"
  }
  else{
    test.sub.data$mem[isub] <- "low_hit"
  }
}

##Importing participant WM data 
wm.sub.data <- c("")
for (isub in 1:length(sub.num)){
  curr.sub <- read.csv(file = paste("s", sub.num.wm[isub], "/", "eMMT",sub.num.wm[isub], "_WM",
                                    ".tsv", sep = ""), header = TRUE, sep = "\t", quote = "\"",dec = ".")
  curr.sub$subnum <- sub.num.wm[isub]
  wm.sub.data <- merge(wm.sub.data, curr.sub, all=T)
}

##   Defining summarySE function
##   Gives count, mean, standard deviation, standard error of the mean, and confidence interval (default 95%).
##   data: a data frame.
##   measurevar: the name of a column that contains the variable to be summariezed
##   groupvars: a vector containing names of columns that contain grouping variables
##   na.rm: a boolean that indicates whether to ignore NA's
##   conf.interval: the percent range of the confidence interval (default is 95%)
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}
##Calculating individual d' 

##Plotting study RT by distCond
summary.stud.distCond <- summarySE(data=study.sub.data, measurevar = "RT", groupvars = c("distCond", "resp"), 
                          na.rm = FALSE, conf.interval = .95, .drop = TRUE)
##cleaning dataframe
for (row in 1:nrow(summary.stud.distCond)){
  if(summary.stud.distCond$resp[row] != "j" && summary.stud.distCond$resp[row] != "k"){
    summary.stud.distCond <- summary.stud.distCond[-c(row), ]
  }
}
ggplot(data=summary.stud.distCond, aes(x=distCond, y=RT, group=resp, colour=resp)) + geom_line() + geom_point() + 
  geom_errorbar(aes(ymin=RT-se, ymax=RT+se), width=.1)

##Plotting study RT by distCat
summary.stud.distCat <- summarySE(data=study.sub.data, measurevar = "RT", groupvars = c("distCond", "resp", "distCat"), 
                          na.rm = FALSE, conf.interval = .95, .drop = TRUE)
##cleaning dataframe
for (row in 1:nrow(summary.stud.distCat)){
  if(summary.stud.distCat$resp[row] != "j" && summary.stud.distCat$resp[row] != "k"){
    summary.stud.distCat <- summary.stud.distCat[-c(row), ]
  }
}
ggplot(data=summary.stud.distCat, aes(x=distCat, y=RT, group=resp, colour=resp)) + geom_line() + geom_point() + 
  geom_errorbar(aes(ymin=RT-se, ymax=RT+se), width=.1)

##Plotting jk responses during study
ggplot(data=study.sub.data, aes(x=resp)) + geom_bar()

##Plotting study RT by mem
study.sub.data$mem <- ""
for(curr in 1:nrow(study.sub.data)){
    study.sub.data$mem[curr] <- test.sub.data$mem[test.sub.data$subjNum==study.sub.data$subjNum[curr] 
                                                  & test.sub.data$word==study.sub.data$word[curr]]
}
summary.study.mem <- summarySE(data=study.sub.data, measurevar = "RT", groupvars = c("mem"), 
                              na.rm = FALSE, conf.interval = .95, .drop = TRUE)
ggplot(data=summary.study.mem, aes(x=mem, y=RT, group=1)) + geom_line() + geom_point() + 
  geom_errorbar(aes(ymin=RT-se, ymax=RT+se), width=.1)

##Plotting study RT by mem (jk)
summary.study.mem.jk <- summarySE(data=study.sub.data, measurevar = "RT", groupvars = c("mem","resp"), 
                               na.rm = FALSE, conf.interval = .95, .drop = TRUE)
for (row in 1:nrow(summary.study.mem.jk)){
  if(summary.study.mem.jk$resp[row] != "j" & summary.study.mem.jk$resp[row] != "k"){
    summary.study.mem.jk <- summary.study.mem.jk[-c(row), ]
  }
}
ggplot(data=summary.study.mem.jk, aes(x=mem, y=RT, group=resp, colour=resp)) + geom_line() + geom_point() + 
  geom_errorbar(aes(ymin=RT-se, ymax=RT+se), width=.1)

##Plotting test RT by distCond
summary.test.distCond <- summarySE(data=test.sub.data, measurevar = "RT", groupvars = c("distCond"), 
                                   na.rm = FALSE, conf.interval = .95, .drop = TRUE)
summary.test.distCond <- summary.test.distCond[-c(3), ]
ggplot(data=summary.test.distCond, aes(x=distCond, y=RT, group=1)) + geom_line() + geom_point() + 
  geom_errorbar(aes(ymin=RT-se, ymax=RT+se), width=.1)

##Plotting test RT by distCat
summary.test.distCat <- summarySE(data=test.sub.data, measurevar = "RT", groupvars = c("distCond", "distCat"), 
                                  na.rm = FALSE, conf.interval = .95, .drop = TRUE)
summary.test.distCat <- summary.test.distCat[-c(4), ]
ggplot(data=summary.test.distCat, aes(x=distCat, y=RT, group=1)) + geom_line() + geom_point() + 
  geom_errorbar(aes(ymin=RT-se, ymax=RT+se), width=.1)

##Plotting test RT by old/new
summary.test.oldnew <- summarySE(data=test.sub.data, measurevar = "RT", groupvars = c("OldNew"), 
                                  na.rm = FALSE, conf.interval = .95, .drop = TRUE)
ggplot(data=summary.test.oldnew, aes(x=OldNew, y=RT, group=1)) + geom_line() + geom_point() + 
  geom_errorbar(aes(ymin=RT-se, ymax=RT+se), width=.1)

##Plotting test RT by mem
summary.test.mem <- summarySE(data=test.sub.data, measurevar = "RT", groupvars = c("mem"), 
                              na.rm = FALSE, conf.interval = .95, .drop = TRUE)
ggplot(data=summary.test.mem, aes(x=mem, y=RT, group=1)) + geom_line() + geom_point() + 
  geom_errorbar(aes(ymin=RT-se, ymax=RT+se), width=.1)







