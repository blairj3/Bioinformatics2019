library(ggplot2)
df=data.frame(status=waterstats$status,
              dayssurv=waterstats$days.survived)
ggplot(myData,aes(x=myData$status, y=myData$mean, fill= factor(status)))+stat_summary( geom="bar")+
  ylab("Days Survived at 0% RH")+xlab("Reproductive Status")+
  geom_errorbar(limits, position = dodge, width=0.25)

myData <- aggregate(waterstats$days.survived,
                    by = list(status= waterstats$status),
                    FUN = function(x) c(mean = mean(x), sd = sd(x),
                                        n = length(x)))
myData<-do.call(data.frame, myData)
myData$se<-myData$x.sd/sqrt(myData$x.n) #calc std error
colnames(myData)<-c("status","mean","sd","n","se")
dodge<-position_dodge(width=0.9)
limits<-aes(ymax=myData$mean+myData$se,
            ymin=myData$mean-myData$se)
