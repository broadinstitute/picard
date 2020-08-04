# Script to generate a chart of quality score distribution in a file
# @author Yossi Farjoun

# Parse the arguments
args <- commandArgs(trailing=T)
metricsFile  <- args[1]
outputFile   <- args[2]
bamFile  <- args[3]
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

metrics <- read.table(metricsFile, header=T, nrows=1, sep="\t", skip=firstBlankLine)
histogram <- read.table(metricsFile, header=T, sep="\t", skip=secondBlankLine)



# get maximal value in histogram 
M=max(histogram$READ_LENGTH)

#introduce rows for the missing ones
histogram=merge(histogram,data.frame(y=seq(0,M)),by.x = "READ_LENGTH",by.y="y",all.y = T)

#melt the histogram
melted=reshape(histogram,idvar="READ_LENGTH", varying = list(2:ncol(histogram)),direction = "long", timevar = "variable", v.names="value",  times=names(histogram)[2:ncol(histogram)],new.row.names = NULL)
rownames(melted) <- c()
#complete the missing values to zero
melted[!complete.cases(melted),"value"]=0

#remove the trailing "_LENGTH_COUNT" from the name of the variable
melted$variable=gsub("_LENGTH_COUNT","",melted$variable)

# Then plot as a PDF
pdf(outputFile)

title=paste("Read Length Distribution\nin file ",bamFile," ",ifelse(subtitle == "","",paste("(",subtitle,")",sep="")),sep="")

histograms=factor(unique(melted$variable))
plot(melted$READ_LENGTH,
     melted$value,
     type="n",
     xlab="Read-Length",
     ylab="Count",
     main=title)

cols=rainbow(length(histograms))
for (v in histograms){
    lines(melted$READ_LENGTH[melted$variable==v],melted$value[melted$variable==v],type='l',col=cols[match(v,histograms)],lwd=3)
}
legend("topleft", legend=histograms, col=cols,lwd = 3)

dev.off()

