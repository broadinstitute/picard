## script to generate histogram of insert sizes from metrics file
## expecting 3 arguments:
## first is the metrics file with the histogram info
## second is the output file
## third is a name for the plot

args <- commandArgs(trailing=TRUE)
metricsFile <- args[1]
pdfFile <- args[2]
bamName <- args[3]
histoWidth <- ifelse(length(args) < 4, 0, as.numeric(args[4]))

startFinder <- scan(metricsFile, what="character", sep="\n", quiet=TRUE, blank.lines.skip=FALSE)

firstBlankLine=0

for (i in 1:length(startFinder)) {
  if (startFinder[i] == "") {
    if (firstBlankLine==0) {
      firstBlankLine=i+1
    } else {
      secondBlankLine=i+1
      break
    }
  }
}

getCumulative <- function(y, yrange) {
  yNew <- rep(0, nrow(y));
  yLength <- nrow(y)
  ySum <- sum(y[,1])
  for (i in 1:yLength) {
    yNew[i] <- (yrange * sum(as.numeric(y[i:yLength,1])) / ySum)
  }
  return (yNew)
}

histogram <- read.table(metricsFile, header=TRUE, sep="\t", skip=secondBlankLine, comment.char="", quote='', check.names=FALSE)

## The histogram has a fr_count/rf_count/tandem_count for each metric "level"
## This code parses out the distinct levels so we can output one graph per level
headers <- sapply(sub(".fr_count","",names(histogram),fixed=TRUE), "[[" ,1)
headers <- sapply(sub(".rf_count","",headers,fixed=TRUE), "[[" ,1)
headers <- sapply(sub(".tandem_count","",headers,fixed=TRUE), "[[" ,1)

## Duplicate header names could cause this to barf.  But it really shouldn't when we have "All_reads.fr_count" and 
## "All_reads.rf_count" for example.  Not sure why this would fail, but I care.
if (any(duplicated(headers))) {
  levels = unique(headers[2:length(headers)]);
} else {
  levels <- c()
  for (i in 2:length(headers)) {
    if (!(headers[i] %in% levels)) {
      levels[length(levels)+1] <- headers[i]
    }
  }
}

pdf(pdfFile)

for (i in 1:length(levels)) {
  ## Reconstitutes the histogram column headers for this level
  fr <- paste(levels[i], "fr_count", sep=".")
  rf <- paste(levels[i], "rf_count", sep=".")
  tandem <- paste(levels[i], "tandem_count", sep=".")

  frrange = ifelse(fr %in% names(histogram), max(histogram[fr]), 0)
  rfrange = ifelse(rf %in% names(histogram), max(histogram[rf]), 0)
  tandemrange = ifelse(tandem %in% names(histogram), max(histogram[tandem]), 0)

  yrange <- max(frrange, rfrange, tandemrange)
  xrange <- ifelse(histoWidth > 0, histoWidth, max(histogram$insert_size))

  par(mar=c(5,4,4,4));
  plot(x=NULL, y=NULL,
       type="n",
       main=paste("Insert Size Histogram for", levels[i], "\nin file", bamName),
       xlab="Insert Size",
       ylab="Count",
       xlim=range(0, xrange),
       ylim=range(0, yrange))
  axis(side=4, at=seq(from=0, to=1, by=0.1)*yrange, labels=seq(from=0, to=1, by=0.10));
  mtext(side=4, line=2, text="cumulative fraction of reads > insert size");

  colors <- c()
  labels <- c()

  if (fr %in% names(histogram) ) {
    lines(histogram$insert_size, as.matrix(histogram[fr]),  type="h", col="red")
    lines(histogram$insert_size, getCumulative(histogram[fr], frrange), col="darkred", lty=2)
    colors <- c(colors, "red")
    labels <- c(labels, "FR")
  }
  if (rf %in% names(histogram)) {
    lines(histogram$insert_size, as.matrix(histogram[rf]),  type="h", col="blue")
    lines(histogram$insert_size, getCumulative(histogram[rf], rfrange), col="darkblue", lty=2)
    colors <- c(colors, "blue")
    labels <- c(labels, "RF")
   }
  if (tandem %in% names(histogram)) {
    lines(histogram$insert_size, as.matrix(histogram[tandem]),  type="h", col="orange")
    lines(histogram$insert_size, getCumulative(histogram[tandem], tandemrange), col="darkorange", lty=2)
    colors <- c(colors, "orange")
    labels <- c(labels, "TANDEM")
  }

  ## Create the legend
  legend("topright", labels, fill=colors, col=colors, cex=0.7)
}

dev.off()
