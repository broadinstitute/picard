# Script to generate a chart to display GC bias based upon read starts observed
# in windows along the genome.
#
# @author Tim Fennell
# @author Kylee Bergin
# edited to make multiple charts for multilevel collection

# Parse the arguments
args <- commandArgs(trailing=T)
metricsFile  <- args[1]
summaryMetricsFile <- args[2]
outputFile   <- args[3]
windowSize   <- args[4]

# Figure out where the metrics and the histogram are in the file and parse them out
startFinder <- scan(metricsFile, what="character", sep="\n", quiet=TRUE, blank.lines.skip=FALSE)
startFinder2 <- scan(summaryMetricsFile, what="character", sep="\n", quiet=TRUE, blank.lines.skip=FALSE)
firstBlankLine=0
firstBlankLine2=0
secondBlankLine2=0

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

for (j in 1:length(startFinder2)) {
        if (startFinder2[j] == "") {
                if (firstBlankLine==0) {
                        firstBlankLine2=j+1
                } else {
                        secondBlankLine2=j+1
                        break
                }
        }
}
pdf(outputFile, onefile=TRUE);

summaryMetrics <- read.table(summaryMetricsFile, header=T, sep="\t", skip=firstBlankLine2);
num.plots <- nrow(summaryMetrics);
my.plots <- vector(num.plots, mode='list');

for (k in 1:(num.plots)){
    if (k==1){
        metrics <- read.table(metricsFile, header=T, sep="\t", skip=(firstBlankLine), nrows=101);
        heads = colnames(metrics);
    }else{
        metrics <-read.table(metricsFile, header=F, sep="\t", skip=(firstBlankLine+((k-1)*101+1)), nrows=101);
        colnames(metrics) <- heads;
    }
    heads = colnames(metrics);
    # Some constants that are used below
    Y_AXIS_LIM = 2;
    MAX_QUALITY_SCORE = 40;
    COLORS = c("royalblue", "#FFAAAA", "palegreen3");

    # Adjust to give more margin on the right hand side
    par(mar = c(5, 4, 4, 4));

    accLevel = summaryMetrics[k,"ACCUMULATION_LEVEL"];
    if(accLevel=="All Reads"){datasetName <- "All Reads";}
    else{
    	if(accLevel == "Sample"){datasetName <- summaryMetrics[k, "SAMPLE"]}
    	if(accLevel == "Library"){datasetName <- summaryMetrics[k, "LIBRARY"]}
    	else(datasetName <- summaryMetrics[k, "READ_GROUP"])}
    subtitle = cat("Total clusters: ",summaryMetrics[k,"TOTAL_CLUSTERS"],", Aligned reads: ",summaryMetrics[k, "ALIGNED_READS"]);
    # Do the main plot of the normalized coverage by GC
    plot(type="p", x=metrics$GC, y=metrics$NORMALIZED_COVERAGE,
        xlab=paste(c("GC% of", windowSize, "base windows"), sep=" ", collapse=" "),
        ylab="Fraction of normalized coverage",
        xlim=c(0,100),
        ylim=c(0, Y_AXIS_LIM),
        col=COLORS[1],
        main=paste(accLevel, "Level:", datasetName, "GC Bias Plot", "\n", subtitle)
        );

    # Add lines at the 50% GC and coverage=1
    abline(h=1, v=50, col="lightgrey");

    # Plot count of windows as a separate series near the bottom
    window_ratio = 0.5 / max(metrics$WINDOWS);
    scaled_windows = metrics$WINDOWS * window_ratio;
    lines(metrics$GC, scaled_windows, type="h", col=COLORS[2], lwd=3);

    # Plot the quality series
    lines(metrics$GC, metrics$MEAN_BASE_QUALITY * Y_AXIS_LIM / MAX_QUALITY_SCORE, type="l", col=COLORS[3]);
    axis(4,
        at=c(0, Y_AXIS_LIM/4, Y_AXIS_LIM/4*2, Y_AXIS_LIM/4*3, Y_AXIS_LIM),
        labels=c(0, MAX_QUALITY_SCORE/4, MAX_QUALITY_SCORE/4*2, MAX_QUALITY_SCORE/4*3, MAX_QUALITY_SCORE)
        );
    mtext("Mean base quality", side=4, line=2.5);

    # And finally add a legend
   	legend("topleft", pch=c(1,15, 45), legend=c("Normalized Coverage", "Windows at GC%", "Base Quality at GC%"), col=COLORS)
   	p=recordPlot()
   	my.plots[[k]] <- p
}

for (my.plot in my.plots){
    replayPlot(my.plot)
}

graphics.off();