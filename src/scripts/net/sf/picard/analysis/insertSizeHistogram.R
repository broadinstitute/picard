# script to generate histogram of insert sizes from metrics file
# expecting 3 arguments:
# first is the metrics file with the histogram info
# second is the output file
# third is a name for the plot

args <- commandArgs(trailing=T)
startFinder <- scan(args[1], what="character", sep="\n", quiet=TRUE, blank.lines.skip=FALSE)

firstBlankLine=0

for (i in 1:length(startFinder))
{
        if (startFinder[i] == "") {
                if (firstBlankLine==0) {
                        firstBlankLine=i+1
                } else {
                        secondBlankLine=i+1
                        break
                }
        }
}

metrics <- read.table(args[1], header=T, nrows=1, skip=firstBlankLine)
histogram <- read.table(args[1], header=T, skip=secondBlankLine)

xrange <- max(histogram$insert_size)
yrange <- max(histogram$fr_count, histogram$rf_count, histogram$tandem_count)

if(length(args) >= 4) {
    xrange <- as.numeric(args[4])
}

pdf(args[2])
plot(x=NULL, y=NULL,
    type="n",
    main=paste(args[3], "Insert Size Histogram"),
    xlab="Insert Size",
    ylab="Count",
    xlim=range(0, xrange),
    ylim=range(0, yrange))



colors <- c()
labels <- c()

if( length(histogram$fr_count) > 0 )
{
    lines(histogram$insert_size, histogram$fr_count,  type="h", col="red")
    colors <- c(colors, "red")
    labels <- c(labels, "FR")
}

if( length(histogram$rf_count) > 0 )
{
    lines(histogram$insert_size, histogram$rf_count,  type="h", col="blue")
    colors <- c(colors, "blue")
    labels <- c(labels, "RF")
}

if( length(histogram$tandem_count) > 0 )
{
    lines(histogram$insert_size, histogram$tandem_count,  type="h", col="orange")
    colors <- c(colors, "orange")
    labels <- c(labels, "TANDEM")
}

paste(labels, "=", colors) # Print legend to console

# Create the legend
legend("topright", labels, fill=colors, col=colors, cex=0.7);

dev.off()

