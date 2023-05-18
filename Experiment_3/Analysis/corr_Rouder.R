# {r,loadData}
link="https://raw.githubusercontent.com/amthapar/data-east/main/antiSaccade/as12/as12All.csv"
dat=read.csv(url(link))
dat$sub=as.integer(as.factor(dat$sub))



# {r,cleanFixed}
goodParticipants= (dat$sub != 26)
goodBlock=dat$blk>3
goodTrial.Fix=dat$trl>15
clean=dat[goodParticipants &goodBlock &goodTrial.Fix,]


# {r preprocess}
# sub: part_id
#blk: block no, three first was training 
# blkType: prosaccade = 0, antisaccade = 0
# trl: no of trial 15 first was training IN EACH BLOCK
# dur: duration of what? presentation time?
# isLet 0 or 1, IDK
# oriT oriC: positions of a que around center, unrelevant
# resp: 0 and 1, somehow related to reaction letters X and M
# rt: reaction time
# acc: correctness

meanDur=tapply(clean$dur,list(clean$sub,clean$blkType),mean)*1000/75
# mean duration on pro and antisaccade task, per participant
meanBlkDur=tapply(clean$dur,list(clean$sub,clean$blk),mean)*1000/75
# mean duration per participant per block (4 of them, 4, 5, 6, 7)
blkType=tapply(clean$blkType,list(clean$sub,clean$blk),mean)
# Just a type of a block (anti or pro) per participant 
#   4 5 6 7
#1  0 1 1 0
#2  0 1 1 0
S=dim(meanBlkDur)[1] # no of participants
pro=anti=matrix(nrow=S,ncol=2)
for (s in 1:S){ # for each participants
  pro[s,]=meanBlkDur[s,blkType[s,]==0]
  anti[s,]=meanBlkDur[s,blkType[s,]==1]}
# po 2 bloki są pro a anty, powyższa pętla rozbija na bloki tylko pro i 
# tylko anty, po 2 
y=round(c(as.vector(t(pro)),as.vector(t(anti))),2) 
# all pro mean duration and all mean anty dur after that
sub=rep(rep(1:S,each=2),2)
task=rep(1:2,each=S*2)
datByBlk=data.frame(sub,task,y)
# subject pro/anty and mean dur
write.table(file='datForModel.dat',quote=F,row.names=F,datByBlk)

blk=rep(1:2,S*2)
m=tapply(datByBlk$y,list(datByBlk$sub,blk,datByBlk$task),mean)
# two df's each for pro/anti task
# row = subject and two columns corr to 2 blocks for participant





#{R correlations}
allTimes=cbind(meanDur,meanDur[,2]-meanDur[,1])
colnames(allTimes)=c("pro","anti","diff")
# for each subject pro,anti and diff mean dur 
cors=round(cor(allTimes),2)
ci=round(cor.test(allTimes[,1],allTimes[,2])$conf.int,2)
