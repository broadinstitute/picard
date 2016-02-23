##
# Copyright (c) 2016, Nils Homer
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
# 1. Redistributions of source code must retain the above copyright notice,
# this list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright notice,
# this list of conditions and the following disclaimer in the documentation
# and/or other materials provided with the distribution.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.
##

# Script to generate a chart for the fraction of bases at a given mean fold coverage.
# @author Nils Homer

# Parse the arguments
args <- commandArgs(trailing=T)
metricsFile  <- args[1]
outputFile   <- args[2]
bamFile  <- ifelse(length(args) < 3, NA, args[3])
subtitle <- ifelse(length(args) < 4, "", args[4])

# Figure out where the metrics and the histogram are in the file and parse them out
startFinder <- scan(metricsFile, what="character", sep="\n", quiet=TRUE, blank.lines.skip=FALSE)

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

metrics <- read.table(metricsFile, header=T, nrows=2, sep="\t", skip=firstBlankLine)
histogram <- read.table(metricsFile, header=T, sep="\t", skip=secondBlankLine)

coverages = rbind(histogram$coverage, histogram$coverage)
counts = rbind(histogram$count_WHOLE_GENOME, histogram$count_NON_ZERO_REGIONS)
labels = c("Whole Genome", "Non-Zero Regions")
colors = c("blue", "green")

ymins = c();
ymaxs = c();
percentOfMeans = c()
percentCovereds = c();

for (i in 1:2) {
	coverage = coverages[i,];
	count = counts[i,];

	coverage = coverage[!is.na(count)];
	count = count[!is.na(count)];

	meanCoverage = metrics$MEAN_COVERAGE[i];
	percentOfMean <- coverage / meanCoverage; # x-axis
	percentCovered <- rep(0, length(count)); # y-axis

	# must do a cumulative sume of percentCovered
	totalCount = sum(as.numeric(count));
	for (j in 1:length(percentCovered)) {
		percentCovered[j] = 100.0# sum(as.numeric(count[j:length(percentCovered)])) / totalCount;
	}

	ymin = percentCovered[round(meanCoverage+1)]
	ymax = min(100,max(percentCovered));

	ymins = append(ymins, ymin);
	ymaxs = append(ymaxs, ymax);
	percentOfMeans = append(percentOfMeans, list(percentOfMean));
	percentCovereds = append(percentCovereds, list(percentCovered));
}

ymin = min(ymins);
ymax = max(ymaxs);

# Then plot the histogram as a PDF
pdf(outputFile);

plot(x=c(0, 1.0),
		y=c(ymin, ymax),
		xlim=c(0, 1.0),
		ylim=c(ymin, ymax),
		type="n",
		main=paste("WGS Base Coverage Plot", ifelse(is.na(bamFile),"",paste("\nin file",bamFile))," ",ifelse(subtitle == "","",paste("(",subtitle,")",sep="")),sep=""),
		xlab="Fold Coverage of Mean",
		ylab="% of Bases Covered");

for (i in 1:2) {
	label = labels[i];
	color = colors[i]
	percentOfMean = percentOfMeans[[i]];
	percentCovered = percentCovereds[[i]];

	lines(percentOfMean, percentCovered, col=color, lwd=5);
}

legend(x="topright", legend=labels, lwd=5, col=colors);

dev.off()

