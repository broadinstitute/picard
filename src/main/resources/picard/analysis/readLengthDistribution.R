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

melted <- tryCatch (
        {
        histogram <- read.table(metricsFile, header=T, sep="\t", skip=secondBlankLine)
        M=max(histogram$READ_LENGTH)
        histogram=merge(histogram,data.frame(y=seq(0,M)),by.x = "READ_LENGTH",by.y="y",all.y = T)
        melted=reshape(histogram,idvar="READ_LENGTH", varying = list(2:ncol(histogram)),direction = "long", timevar = "variable", v.names="value",  times=names(histogram)[2:ncol(histogram)],new.row.names = NULL)
        rownames(melted) <- c()
        melted[!complete.cases(melted),"value"]=0
        melted$variable=gsub("_LENGTH_COUNT","",melted$variable)
        melted
        },
  error = function(cond) {
          if (grepl("no lines available in input", cond$message)) {
                melted <- data.frame(READ_LENGTH=c(100,100), value=c(0,0), variable=c("PAIRD_TOTAL", "PAIRED_ALIGNED"))
          } else {
                stop(cond)
          }
        }
  )

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

