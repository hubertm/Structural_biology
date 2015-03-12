library("pracma")

SEC<-function (time, intensity, range = list(start, stop), color = "blue", 
               xlab = "retention time", ylab = "intensity", ylim = c(0,
               max(intensity) * 1.1), las = 0, plotname="SECrun",...) 
{
  if (length(range$start) != length(range$stop)) 
    stop("The number of start and stop points must be equal")
  if (length(color) != 1 & length(color) != length(range$start)) 
    stop("Number of different colors must match number of peaks")
  intensitybaseline<-intensity-mean(intensity[time<5])
  #SECsmoothed<-lowess(time, intensitybaseline, f=0.001)
  SECsmoothed<-locpoly(time, intensitybaseline, bandwidth=0.25)
  #if peaks are not detected (signal overload) drop this value to f=0.01
  SECsmoothed_x<-unlist(SECsmoothed[1])
  SECsmoothed_y<-unlist(SECsmoothed[2])
  SECpeaks<-findpeaks(SECsmoothed_y, threshold = 5)
  
  plot(time, intensitybaseline, type = "l", ylim = ylim, xlab = xlab, 
       ylab = ylab, las = las, main = plotname)
  points(time[SECpeaks[,2]], SECpeaks[,1], col="red", pch=20)
  peakposition<-vector(length=length(SECpeaks[,2]))
  for (k in 1:length(SECpeaks[,2])){
    peakposition[k]<- intensitybaseline[which(abs(time-SECsmoothed$x[SECpeaks[,2][k]])==min(abs(time-SECsmoothed$x[SECpeaks[,2][k]])))]
  }
  points(SECsmoothed$x[SECpeaks[,2]],peakposition , col="blue", pch=20)
    
  usefulrangeforpeaks<-c(5,18)   #to avoid stuff outside of this range
  smoothed2ndderv<-locpoly(SECsmoothed_x,SECsmoothed_y, bandwidth=0.25, drv=2,range.x = usefulrangeforpeaks)
  seconddervpeaks<-findpeaks(smoothed2ndderv$y, threshold=10)
  #plot(smoothed2ndderv, type="l")
  #points(smoothed2ndderv$x[seconddervpeaks[,2]], seconddervpeaks[,1], col="red", pch=20)
  
  turningpoints<-vector(length=length(seconddervpeaks[,2]))
  for (h in 1:length(seconddervpeaks[,2])){
    turningpoints[h]<- intensitybaseline[which(abs(time-smoothed2ndderv$x[seconddervpeaks[,2][h]])==min(abs(time-smoothed2ndderv$x[seconddervpeaks[,2][h]])))]
  }
  points(smoothed2ndderv$x[seconddervpeaks[,2]],turningpoints , col="blue", pch=20)
  retentionTime <- vector(mode = "numeric", length = length(range$start))
  peakArea <- vector(mode = "numeric", length = length(range$start))
  maxIntensity <- vector(mode = "numeric", length = length(range$start))
  for (i in 1:length(range$start)) {
    peakTime <- time[time >= range$start[i] & time <= range$stop[i]]
    peakIntensity <- intensitybaseline[time >= range$start[i] & time <= range$stop[i]]
    #peakIntensitysmoothed<-SECsmoothed_y[time >= range$start[i] & time <= range$stop[i]]
    polygonx=c(range$start[i],peakTime,range$stop[i])
    polygony=c(0,peakIntensity,0)
    peakTime=peakTime[4:(length(peakTime)-3)]
    if (length(color) == 1) 
      peakColor <- color
    else peakColor <- color[i]
    polygon(polygonx, polygony, col = peakColor)
    n <- length(peakTime)-1
    x <- vector(mode = "numeric", length = n)
    for (j in 1:n) {
      x[j] <- (((peakIntensity[j]+peakIntensity[j+1])/2) * (peakTime[j+1]-peakTime[j]))
    }
    peakArea[i] <- abs(sum(x))
    retentionTime[i] <- peakTime[which.max(peakIntensity)]
    maxIntensity[i] <- max(peakIntensity)
  }
  print(data.frame(retentionTime, peakArea, maxIntensity))
}

fluo=read.csv("processed_t.csv")
#time<-fluo[,1]
#intensity<-fluo[,2]
#SEC(time,intensity, range=list(start=10.3, stop=12))
for (i in 1:(length(fluo)/2)){
#for (i in 1:1){
  time<-fluo[,2*i-1]
  intensity<-fluo[,2*i]
  SEC(time,intensity, range=list(start=10.3, stop=12), plotname = names(fluo[i*2-1]))
}
