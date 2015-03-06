SEC<-function (time, intensity, range = list(start, stop), color = "blue", 
               xlab = "retention time", ylab = "intensity", ylim = c(0, 
                                                                     max(intensity) * 1.1), las = 1, ...) 
{
  if (length(range$start) != length(range$stop)) 
    stop("start and stop vectors in the range list must be equal")
  if (length(color) != 1 & length(color) != length(range$start)) 
    stop("color vector length must be 1 or equal to the vectors in the range list")
  plot(time, intensity, type = "l", ylim = ylim, xlab = xlab, 
       ylab = ylab, las = las, ...)
  retentionTime <- vector(mode = "numeric", length = length(range$start))
  peakArea <- vector(mode = "numeric", length = length(range$start))
  apexIntensity <- vector(mode = "numeric", length = length(range$start))
  for (i in 1:length(range$start)) {
    peakTime <- time[time >= range$start[i] & time <= range$stop[i]]
    peakIntensity <- intensity[time >= range$start[i] & time <= 
                                 range$stop[i]]
    if (length(color) == 1) 
      peakColor <- color
    else peakColor <- color[i]
    polygon(peakTime, peakIntensity, col = peakColor)
    n <- length(peakTime)-1
    x <- vector(mode = "numeric", length = n)
    for (j in 1:n) {
      x[j] <- (((peakIntensity[j]+peakIntensity[j+1])/2) * (peakTime[j+1]-peakTime[j]))
    }
    peakArea[i] <- abs(sum(x))
    retentionTime[i] <- peakTime[peakIntensity == max(peakIntensity)]
    print (peakIntensity == max(peakIntensity))
    apexIntensity[i] <- max(peakIntensity)
  }
  return(data.frame(retentionTime, peakArea, apexIntensity))
}

