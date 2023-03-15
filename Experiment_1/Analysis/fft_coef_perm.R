#--------------------------------------------------------------------
# Multivariate Permutation Paired T-tests of
# Fourier-Transformed Data Coefficients (for R 2.5.1)
#
# Written By: Jennifer Urbano Blackford, PhD
#         Ronald Salomon, MD
#         Niels Waller, PhD
#
# User instructions:
#   User changes may be made in the first section of the program
#   entitled “User Specifications”. The program expects 2
#   separate raw data files, one for pre and one for post.
#   Data should be comma delimited with subjects as rows and time
#   in the columns. For 10,000 permutations on 20 subjects and
#   144 time points processing time is approximately 5 min
#   on a Dell Precision (XP version 2002, 3 GHz,
#   3.25 GB RAM, R version 2.5.1)
#
#   1. Raw Data Files:
#      Put raw pre- and post-treatment data into 2 separate comma
#      delimited files names “predata.csv” and “postdata.csv”.
#      Either put the files in the C:/directory or change the drive
#      below in “User Specifications”.
#   2. Output File:
#      Default output file is C:/output.txt.
#   3. Specify the number of permutations (Number.of.loops) below.
#      Default is 10,000.
#-------------------------------------------------------------------- 
#----------------User Specifications-----------------------
#specify pre-treatment data
data.pre <- as.matrix(read.csv("c: / predata.csv"))
#specify post-treatment data
data.post <- as.matrix(read.csv("c: / postdata.csv"))
#specify output filename
file.name <- "c: /  / Output.txt"
sink(file.name)
#specify number of permutations, should be 5000–10000
Number.of.loops <- 10000
#----------------------------------------------------------
#load necessary libraries
library(Hmisc)
library(MASS)
#-------------------Create Variables----------------------
#create data object
data <- list(data.pre = data.pre, data.post = data.post)
#create variables
subjects <- nrow(data.pre)
time <- ncol(data.pre)
#----------------------------------------------------------
#----------------------------Functions---------------------
#hanning window
hanning <- function(n) {
  t(.5 * (1-cos(2 * pi * t(1:n) / (n + 1))))
}
#calculate Fourier coefficients for pre and post data
Fourier <- function(subjects, time, data.pre, data.post) {
  f.pre <- matrix(0, subjects, time / 2)
  f.post <- matrix(0, subjects, time / 2)
  for (i in 1:subjects) {
    a1 <- fft(han * data.pre[i, ])[1:(time / 2)]
    a2 <- fft(han * data.post[i, ])[1:(time / 2)]
    f.pre[i, ] <- Re(a1 * Conj(a1))
    f.post[i, ] <- Re(a2 * Conj(a2))
  }
  Fourier <- list(f.pre = f.pre, f.post = f.post)
  return(Fourier)
}
#----------------------------------------------------------------------
#-----------------------------Analysis---------------------------------
#variables
nvar <- time / 2
han <- hanning(time)
nah <- 1 / han
output <- 0
continue <- 1
#Fourier Transform data
f.data <- Fourier(
  subjects = subjects,
  time = time,
  data.pre = data$data.pre,
  data.post = data$data.post
)
#create difference score for paired t-test
f.diff <- f.data$f.post - f.data$f.pre
#name columns
colnames(f.diff) <- colnames(f.diff, do.NULL = FALSE, prefix = "Var")
colnames <- dimnames(f.diff)[[2]]
#Multivariate Permutation Test
#create matrices to store values of t(1), t significance(2)
# t.test<-rep(0,nvar)
obs.t <- matrix(0, 2, nvar)
t.test <- apply(f.diff, 2, t.test)
for (i in 1:nvar) {
  obs.t[1, i] <- as.matrix(t.test[[i]]$statistic)
  obs.t[2, i] <- as.matrix(t.test[[i]]$p.value)
}
writeLines ("Observed t and p values")
print(t(round(obs.t, 2)))
writeLines("")
writeLines ("Significant Variables from MPT for difference t - test")
writeLines("")
#create simulated distribution of diff values
t.sim.all <- matrix(0, Number.of.loops, nvar)
for (loop in 1:Number.of.loops) {
  sign <- rep(1, subjects)
  sign[runif(subjects, 0, 1) <= .5] <- -1
  newdata <- f.diff * sign
  #save the t values for each variable and each loop into a matrix called
  t.simulated.all
  new.t.test <- apply(newdata, 2, t.test)
  for (i in 1:nvar) {
    t.sim.all[loop, i] <- as.matrix(new.t.test[[i]]$statistic)
  }
} #end of this do loop
#create matrices of maximum ts and use in each permutation
#loop for each variable to test significance--step down procedure
while (continue == 1) {
  #select largest absolute value of simulated ts from each loop--save column
  #position
  max.col <- matrix(0, Number.of.loops)
  max.sim.t <- matrix(0, Number.of.loops)
  for (loop in 1:Number.of.loops) {
    #select the largest absolute t
    max.col[loop] <- (order(abs(t.sim.all[loop, ]))[nvar])
    #save the largest value--maintain original sign
    max.sim.t[loop] <- t.sim.all[loop, max.col[loop, ]]
  }
  #save t value and traditional p value information for the largest t value
  max.col <- (order(abs(obs.t[1, ]))[nvar])
  largest.obs.t <- obs.t[1, max.col]
  smallest.obs.p <- obs.t[2, max.col]
  variable <- colnames[max.col]
  #compare largest t to the distribution to determine significance by
  #creating an array of 0s and 1s representing whether the simulated m is larger 
  #than the observed t
  p.dist <- rep(0, Number.of.loops)
  p.dist[(abs(max.sim.t)) >= (abs(largest.obs.t))] <- 1
  #get p value as number of permutations above/below obs t divided by number of 
  #permutations
  p.value <- 0
  if (sum(p.dist) > 0)
    (p.value <- (sum(p.dist) / Number.of.loops))
  print ("Variable, t value, p value, MPT p value")
  print (variable)
  print (largest.obs.t)
  print (smallest.obs.p)
  print (p.value)
  print ("--------------------------------")
  #if the p value was significant then delete the column and continue
  #otherwise stop the process
  if (p.value > .05)
    (continue <- 0) else {
      #delete the largest column and continue
      #delete the column associated with the maximum observed t 
      t.sim.all<-as.matrix(t.sim.all[,-max.col])
      obs.t<-as.matrix(obs.t[,-max.col])
      colnames<-as.matrix(colnames[-max.col])
      #reset number of variables
      nvar<-nvar-1
      if (nvar==0) (continue<-0)}
} #end of program