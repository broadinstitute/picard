# Script to generate a normalized coverage vs. position along transcript plot.
#
# @author Tim Fennell

# Parse the arguments
args <- commandArgs(trailing=T)
metricsFile  <- args[1]
outputFile   <- args[2]
datasetName  <- args[3]

# Figure out where the metrics and the histogram are in the file and parse them out
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

data <- read.table(metricsFile, header=T, sep="\t", skip=secondBlankLine)
pdf(outputFile)

# Some constants that are used below
COLORS = c("royalblue", "#FFAAAA", "palegreen3");

# Do the main plot of the normalized coverage by GC
plot(type="o", x=data$normalized_position, y=data$normalized_coverage,
     xlab="Normalized Distance Along Transcript",
     ylab="Normalized Coverage",
     xlim=c(0,100),
     ylim=c(0, max(data$normalized_coverage)),
     col="royalblue",
     main=paste("RNA-Seq Coverage vs. Transcript Position", "\n", datasetName)
);

# Add a horizontal line at coverage=1
abline(h=1, col="lightgrey");

dev.off();