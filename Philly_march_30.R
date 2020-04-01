d<-read.csv('/Users/mzlevy/Dropbox/coronavirus/Philly/Mapping_data.csv')
d10<-d[d$All.Tests>10,]
attach(d)

Tot<-d10$All.Tests
Pos<-d10$Positives

counts <- table(Tot,Pos)
barplot(counts, main="Preavlence Covid by Zip",
  xlab="Number Positive our of Total Tests", col=c("darkblue","red"),
  legend = rownames(counts))

plot(Pos,Tot, xlab="Number Positives", ylab="Number Tested", ylim=c(0,400), main="Philly Covid Testing March 30, 2020")