SEC<-function (time, intensity, range = list(start, stop), color = "blue", 
               xlab = "retention time", ylab = "intensity", ylim = c(0,
               max(intensity) * 1.1), las = 1, ...) 
{
  if (length(range$start) != length(range$stop)) 
    stop("start and stop vectors in the range list must be equal")
  if (length(color) != 1 & length(color) != length(range$start)) 
    stop("color vector length must be 1 or equal to the vectors in the range list")
  intensitybaseline<-intensity-mean(intensity[time<5])
  plot(time, intensitybaseline, type = "l", ylim = ylim, xlab = xlab, 
       ylab = ylab, las = las, ...)
  retentionTime <- vector(mode = "numeric", length = length(range$start))
  peakArea <- vector(mode = "numeric", length = length(range$start))
  maxIntensity <- vector(mode = "numeric", length = length(range$start))
  for (i in 1:length(range$start)) {
    peakTime <- time[time >= range$start[i] & time <= range$stop[i]]
    peakIntensity <- intensitybaseline[time >= range$start[i] & time <= range$stop[i]]
    polygonx=c(range$start,peakTime,range$stop)
    polygony=c(0,peakIntensity,0)
    peakTime=peakTime[4:(length(peakTime)-3)]
    if (length(color) == 1) 
      peakColor <- color
    else peakColor <- color[i]
    polygon(polygonx, polygony, col = peakColor)
    n <- length(peakTime)-1
    print(n)
    x <- vector(mode = "numeric", length = n)
    for (j in 1:n) {
      x[j] <- (((peakIntensity[j]+peakIntensity[j+1])/2) * (peakTime[j+1]-peakTime[j]))
    }
    peakArea[i] <- abs(sum(x))
    retentionTime[i] <- peakTime[peakIntensity == max(peakIntensity)]
    maxIntensity[i] <- max(peakIntensity)
  }
  return(data.frame(retentionTime, peakArea, maxIntensity))
}

