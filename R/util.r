#' Inspect effects of thresholds on matches over time
#' 
#' If it can be assumed that matches should only occur within a given time range (e.g., event data should match news items after the event occured)
#' a low effort validation can be obtained by looking at whether the matches only occur within this time range. 
#' This function plots the percentage of matches within a given time range (hourdiff) for different thresholds of the weight column.
#' This can be used to determine a good threshold.
#'
#' @param g       The output of newsflow.compare (either as "igraph" or "edgelist")
#' @param breaks  The number of breaks for the weight threshold
#' @param hourdiff_range The time period (hourdiff range) in which the match 'should' occur.
#' @param min_weight   Optionally, filter out all value below the given weight
#' @param min_hourdiff the lowest possible hourdiff value. This is used to estimate noise. If not specified, will be estimated based on data.
#' @param max_hourdiff the highest possible hourdiff value. 
#'
#' @return Nothing... just plots
#' @export
hourdiff_range_thresholds <- function(g, breaks=20, hourdiff_range=c(0,Inf), min_weight=NA, min_hourdiff=NA, max_hourdiff=NA) {
  weight = if(methods::is(g,'igraph')) E(g)$weight else g$weight
  hourdiff = if(methods::is(g,'igraph')) E(g)$hourdiff else g$hourdiff
  
  if (!is.na(min_weight)) {
    filter = weight >= min_weight
    weight = weight[filter]
    hourdiff = hourdiff[filter]
  }

  if (is.na(min_hourdiff)) min_hourdiff = min(hourdiff)
  if (is.na(max_hourdiff)) max_hourdiff = max(hourdiff)
  
  thresholds = seq(min(weight), max(weight), length.out=breaks)
  res = data.frame(threshold=thresholds, n=NA, after_event=NA, est_noise=NA)
  for (i in 1:length(thresholds)) {
    ti = weight > thresholds[i]
    bz = hourdiff[ti] >= hourdiff_range[1] & hourdiff[ti] < hourdiff_range[2]
    if (sum(bz) == 0) break
    res$n[i] = sum(bz)
    res$after_event[i] = sum(bz) / length(bz)
    
    ## estimate noise
    total_hourdiff = max_hourdiff - min_hourdiff
    min_correct_hourdiff = if (hourdiff_range[1] < min_hourdiff) min_hourdiff else hourdiff_range[1]
    max_correct_hourdiff = if (hourdiff_range[2] > max_hourdiff) max_hourdiff else hourdiff_range[2]
    correct_hourdiff = max_correct_hourdiff - min_correct_hourdiff
    r = (total_hourdiff - correct_hourdiff) / correct_hourdiff
    res$est_noise[i] = (sum(!bz) / r) / sum(bz)
    if (res$est_noise[i] > 1) res$est_noise[i] = 1
  }
  res$n[is.na(res$n)] = 0
  
  graphics::layout(matrix(c(1,3,2,3), nrow=2))
  bars = graphics::hist(weight, main=sprintf('Histogram of weight'), breaks=res$threshold)
  for (i in 2:length(bars$breaks)) {
    graphics::rect(xleft = bars$breaks[i-1], ybottom = 0, xright = bars$breaks[i], ytop = res$n[i-1] - sum(res$n[i]), col='blue')
  }
  
  graphics::hist(hourdiff, main=sprintf('Histogram of hourdiff'), breaks=nrow(res))
  graphics::rect(xleft = hourdiff_range[1], ybottom = 0, xright = hourdiff_range[2], ytop = res$n[1], col='blue')
  graphics::rect(xleft = hourdiff_range[1], ybottom = 0, xright = hourdiff_range[2], ytop = res$n[1]*res$est_noise[1], col='red')
  
  graphics::par(mar = c(5,5,2,5))
  graphics::plot(as.numeric(res$threshold), res$after_event, type='l', col='blue', xlab='weight threshold', ylab='pct in range / est noice', ylim=c(0,1))
  graphics::lines(as.numeric(res$threshold), res$est_noise, type='l', lty=2, col='red', xlab='weight threshold', ylab='pct in range / est noice')
  graphics::par(new = T)
  graphics::plot(as.numeric(res$threshold), res$n, pch=10, axes=F, xlab=NA, ylab=NA, cex=1.2)
  graphics::axis(side = 4)
  graphics::mtext(side = 4, line = 3, 'N in range')
  graphics::legend("topleft",
         legend=c("pct in range", "estimated noise", "N in range"),
         lty=c(1,2,0), pch=c(NA, NA, 16), col=c("blue", "red","black"))
}



