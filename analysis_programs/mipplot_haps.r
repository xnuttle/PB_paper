#Xander Nuttle
#mipplot_haps.r
#Call: source(mipplot_haps.r)
#
#General purpose program for plotting haplotype/paralog-specific MIP read count frequencies for
#scenarios involving only two haplotypes/paralogs. Called by pdf_pb_mips.r (used in analysis of PB MIP pool).
#Requires the following exisiting variables in R: experiment_name, guide_locations, individual, gene

#load ggplot2 library
library(ggplot2)

#input mip count data
data<-read.table(experiment_name,header=T)

#input locations of guide RNAs
if(file.exists(guide_locations))
{
	gdata<-read.table(guide_locations,header=F)
}

#extract rows having data for the individual being analyzed
rows<-which(data[,1]==individual)
data2<-data[1:length(rows),]
for(i in 1:length(rows))
{
	data2[i,]<-data[rows[i],]
}

#extract rows having data for the gene of interest
rows<-which(data2[,2]==gene)
data3<-data2[1:length(rows),]
for(i in 1:length(rows))
{
	data3[i,]<-data2[rows[i],]
}

#create column having total read counts per MIP and calculate median number of reads per MIP
for(i in 1:length(data3[,1]))
{
  data3[i,7]<-sum(data3[i,5:6],na.rm=TRUE)
}
median_count<-median(data3[,7])

#calculate paralog-specific count frequencies for each MIP and create a data frame to store this info
#also calculate and store information on MIP performance in terms of reads per MIP relative to the median
x<-data3[,3]
f1<-x
f2<-x
alphavec<-x
for(i in 1:length(x))
{
  f1[i]<-(data3[i,5]/(sum(data3[i,5:6],na.rm=TRUE)))
  f2[i]<-(data3[i,6]/(sum(data3[i,5:6],na.rm=TRUE)))
  alphaval<-(data3[i,7]/median_count)
  if(is.nan(alphaval))
  {
    alphaval<-0
  }
  if(alphaval>=1)
  {
    alphavec[i]<-1
  }
  else
  {
    alphavec[i]<-alphaval
  }
}
freqs<-data.frame(x,f1,f2,alphavec)

#create the plot
p<-ggplot(freqs,aes(x,f1,alpha=alphavec))
p<-p+geom_point(color="#FF0000") #haplotype 1/paralog 1 = red
p<-p+geom_point(data=freqs,aes(x,f2,alpha=alphavec),color="#0000FF") #haplotype 2/paralog 2 = blue
p<-p+theme(panel.background = element_rect(fill='darkgray'))
p<-p+coord_cartesian(ylim=c(-0.25,1.1))
p<-p+xlab("Alignment Coordinate")
p<-p+ylab("Haplotype/Paralog-Specific Count Frequency")
p<-p+theme(legend.position="none")
p<-p+ggtitle(individual)
p<-p+scale_y_continuous(minor_breaks=seq(0.1,0.9,0.2),breaks=seq(0,1,0.2))
if(file.exists(guide_locations))
{
	p<-p+geom_vline(xintercept=gdata[,1],linetype='dashed',color="#000000")
}

