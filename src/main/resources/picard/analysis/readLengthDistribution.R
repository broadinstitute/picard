# Script to generate a chart of quality score distribution in a file
# @author Yossi Farjoun

# library(ggplot2)
library(reshape2)

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
melted=melt(histogram,id.vars = "READ_LENGTH")

#complete the missing values to zero
melted[!complete.cases(melted),"value"]=0

#remove the trailing "_LENGTH_COUNT" from the name of the variable
melted$variable=gsub("_LENGTH_COUNT","",melted$variable)

# Then plot as a PDF
pdf(outputFile)

title=paste("Read Length Distribution\nin file ",bamFile," ",ifelse(subtitle == "","",paste("(",subtitle,")",sep="")),sep="")

# ggplot(melted)+theme(legend.position="bottom")+geom_line(stat="identity",aes(x=READ_LENGTH,y=value,color=variable))+labs(y="Count",title=title)

dev.off()

