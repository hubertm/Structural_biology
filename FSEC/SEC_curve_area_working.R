## HM 150330 v 0.01
## In order to be able to pass arguments from the terminal
args <- commandArgs(TRUE)

library("pracma")
library("KernSmooth")
## Function I found online which allows me to return a vector of values
speciallist <- structure(NA,class="result")
"[<-.result" <- function(x,...,value) {
      args <- as.list(match.call())
      args <- args[-c(1:2,length(args))]
      length(value) <- length(args)
      for(i in seq(along=args)) {
            a <- args[[i]]
            if(!missing(a)) eval.parent(substitute(a <- v,list(a=a,v=value[[i]])))
      }
      x
}

SEC<-function (time, intensity,  color = "blue",
               ylab = "intensity", ylim = c(0,1000), las = 0, plotname="SECrun", peakestimate = 12.2, usefulrangeforpeaks=c(5,18), producepdf=TRUE,...)
{
      producepdf=as.logical(producepdf)
      if (isTRUE(producepdf)){
            pdf(file=paste(plotname,".pdf",sep = "")) 
      } 
      intensitybaseline<-intensity-mean(intensity[time<5])
      ylim<-c(min(intensitybaseline),max(intensitybaseline)*1.1)
      ## lowess is one way of fitting the polynomial seems to be a bit better with close peaks
      ## if peaks are not detected (signal overload) drop this value to f=0.01
      ## SECsmoothed<-lowess(time, intensitybaseline, f=0.01)
      SECsmoothed<-locpoly(time, intensitybaseline, bandwidth=0.025)
      ## Get the 2nd derivative of the smoothed curve, smooth more to get less peaks
      SECpeaks<-findpeaks(SECsmoothed$y, threshold = 2)
      ## Generates the 2nd derivative of the curve and marks the points
      #smoothed2ndderv<-locpoly(SECsmoothed$x,SECsmoothed$y, bandwidth=0.1, drv=2,range.x = usefulrangeforpeaks)
      #seconddervpeaks<-findpeaks(smoothed2ndderv$y, threshold=10)
      #turningpoints<-vector(length=length(seconddervpeaks[,2]))
      #for (h in 1:length(seconddervpeaks[,2])){
      #      turningpoints[h]<- intensitybaseline[which(abs(time-smoothed2ndderv$x[seconddervpeaks[,2][h]])==min(abs(time-smoothed2ndderv$x[seconddervpeaks[,2][h]])))]
      #}
      
      ##get the index of the closest peak to be used in the next step but also later
      closestpeakestimateindex<-which.min(abs(SECsmoothed$x[SECpeaks[,2]]-peakestimate))
      ## Write the peak retention time to peakestimate
      calpeakestimate<-as.numeric(SECsmoothed$x[SECpeaks[closestpeakestimateindex,2]])
      ### ------------- Critical Parameter ------------- ###
      ## CHECK if needs to be changed as this assumes a peak between 9 and 14!
      if (calpeakestimate<11 | calpeakestimate>16){
            print("The estimated peak is outside the useful limits. Stopping")
            print(calpeakestimate)
            plot(time, intensitybaseline, type = "l", ylim = ylim, ylab = ylab, las = las, main = plotname)
            abline(v=SECsmoothed$x[SECpeaks[,2]])  ## Shows the positions of the detected peaks
            lines(SECsmoothed$x, SECsmoothed$y,type='l', col="red")
            print("------------------------------------------------------------")
            return(0)
      }
      ## Calculate the left border of the peak
      ## Checks if there is a peak before, if yes finds the minimal intensity between the 2 peaks
      if (closestpeakestimateindex-1>0){
            leftpeak<-as.numeric(SECsmoothed$x[SECpeaks[(closestpeakestimateindex-1),2]])
            rangestart<-time[time<calpeakestimate & time>leftpeak][which.min(intensity[time<calpeakestimate & time>leftpeak])]
      }
      ## Otherwise sets the peak between -2 and +2 ml of the calculated peak
      ### ------------- Critical Parameter ------------- ###
      ## If the peak should be wider or smaller this parameter has to be changed (the 2 at the end of the line)
      else{
            rangestart<-time[time<calpeakestimate & time>(calpeakestimate-2)][which.min(intensity[time<calpeakestimate & time>(calpeakestimate-2)])]
      }
      rangestop<-(calpeakestimate+(calpeakestimate-rangestart))
      
      retentionTime <- vector(mode = "numeric", length = length(rangestart))
      peakArea <- vector(mode = "numeric", length = length(rangestart))
      maxIntensity <- vector(mode = "numeric", length = length(rangestart))
      peakTime <- time[time >= rangestart & time <= rangestop]
      #range to search for the maximum Intensity of the selected peak
      MaxIntensitySelPeakRange<- intensitybaseline[time >= (calpeakestimate-0.1) & time <= (calpeakestimate+0.1)]
      MaxIntensitySelPeakTimeRange<- time[time >= (calpeakestimate-0.1) & time <= (calpeakestimate+0.1)]
      peakIntensity <- intensitybaseline[time >= rangestart & time <= rangestop]
      n <- length(peakTime)-1
      x <- vector(mode = "numeric", length = n)
      for (j in 1:n) {
            x[j] <- (((peakIntensity[j]+peakIntensity[j+1])/2) * (peakTime[j+1]-peakTime[j]))
      }
      peakArea <- round(abs(sum(x)), digits=2)
      retentionTime <- round(MaxIntensitySelPeakTimeRange[which.max(MaxIntensitySelPeakRange)], digits=2)
      maxIntensity <- round(max(peakIntensity), digits=2)
      estsize_kDa<-round(exp((calpeakestimate-38.245)/-4.8975), digits = 1) ##estimate the size based on the S200 calibration 31/1/2015
            
      ## Plotting commands
      plot(time, intensitybaseline, type = "l", ylim = ylim, xlab = paste("Retention volume =", retentionTime, "kDa= ",estsize_kDa),
           ylab = ylab, las = las, main = plotname)
      abline(v=SECsmoothed$x[SECpeaks[,2]])  ## Shows the positions of the detected peaks
      #points(smoothed2ndderv$x[seconddervpeaks[,2]],turningpoints , col="blue", pch=20)
      polygonx=c(rangestart,peakTime,rangestop)
      polygony=c(0,peakIntensity,0)
      peakColor <- color
      polygon(polygonx, polygony, col = peakColor)
      
      if (isTRUE(producepdf)){
            dev.off() 
      }   
      
      print(c(retentionTime=paste(retentionTime," min"), peakArea=peakArea, maxIntensity=paste(maxIntensity," mAU"), Apparent_size=paste(estsize_kDa," kDa")
              ,plotname=plotname))
      print("------------------------------------------------------------")
      return(c(retentionTime, peakArea, maxIntensity, estsize_kDa,plotname))
}

fluo=read.csv("Ready.csv")
## Initializing the values for the read-out
RretentionTime <- vector(mode = "numeric")
RpeakArea <- vector(mode = "numeric")
RmaxIntensity <- vector(mode = "numeric")
Restsize_kDa <- vector(mode = "numeric")
Rplotname <- vector(mode = "numeric")
for (i in 1:(length(fluo)/2)){     
      time<-fluo[,2*i-1]
      intensity<-fluo[,2*i]
      ## The R at the beginning of the variable was added to indicate the (R)esult and to avoid confusion
      speciallist[RretentionTime[i], RpeakArea[i], RmaxIntensity[i], Restsize_kDa[i],Rplotname[i]]<-SEC(time,intensity, plotname = paste("Run ",i), producepdf=args)
}
## In order to remove the NA values of the curves that didn't give a result
ind<-which(is.na(RmaxIntensity))
RmaxIntensity[ind]<-0
RmaxIntensity<-as.numeric(RmaxIntensity)

ind<-which(is.na(RretentionTime))
RretentionTime[ind]<-0
RretentionTime<-as.numeric(RretentionTime)

ind<-which(is.na(RpeakArea))
RpeakArea[ind]<-0
RpeakArea<-as.numeric(RpeakArea)
##Take a sequence of vectors and combine them as the rows of a matrix
result<-rbind(RpeakArea,RmaxIntensity)
print('Peakarea= blue, MaxIntensity=red')
pdf('barplotAREAmaxINTENSITY.pdf')
barplot(result, names.arg = RretentionTime,las=2,horiz = TRUE,col=c("darkblue","red"),beside=TRUE)
dev.off()