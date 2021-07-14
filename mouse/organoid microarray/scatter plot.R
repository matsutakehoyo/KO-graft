#maek scatter plots of micro array data
library(ggplot2)
library(tidyr)                #gather(), separate()


pstyle = theme(plot.title = element_text(lineheight=.8, face="bold"),
               legend.title = element_blank())

load("data.Rda")
length(unique(data$ProbeName))
length(unique(data$GeneName))

#subset data, remove control & low signal rows
data.relevant <- subset(data, ControlType==FALSE)
length(unique(data.relevant$ProbeName))
length(unique(data.relevant$GeneName))
data.relevant <- subset(data.relevant, gIsWellAboveBG==TRUE)
length(unique(data.relevant$ProbeName))
length(unique(data.relevant$GeneName))

#function to make scatter plots
plotscatter <- function(dataset, xval, name){
  p <- ggplot(data=dataset, aes_string(x=xval, y="Normalized")) +
    geom_point(size=0.1, colour="skyblue")+
    geom_line(data=data.frame(x=c(0.01, 1000), y=c(0.01, 1000)), aes(x=x, y=y), colour="deepskyblue", size=0.3) +
    geom_line(data=data.frame(x=c(0.01, 1000), y=c(0.02, 2000)), aes(x=x, y=y), colour="violet", size=0.3) +
    geom_line(data=data.frame(x=c(0.01, 1000), y=c(0.005, 500)), aes(x=x, y=y), colour="violet", size=0.3) +
    scale_x_log10()+
    scale_y_log10()+
    labs(title = paste("Scatter Plot vs", xval, 
                       "(", name, ":", length(unique(dataset$GeneName)), "genes)"))+
    facet_grid(dif.day~Cell.Line)+
    pstyle
  p
  ggsave(paste("Scatter Plot vs", xval, 
               "(", name, length(unique(dataset$GeneName)), "genes).png"))
}

#Scatter Plot: all
data.all.day <- subset(data.relevant, select=grep("Raw", grep("^gIs", colnames(data), value=TRUE, invert=TRUE), value=TRUE, invert=TRUE)) #drop flags & Raw columns
data.all.day <- data.all.day %>% spread(key=dif.day, value=Normalized)
data.all.day <- data.all.day %>% gather(key=dif.day, value=Normalized, DD16:DD23)
plotscatter(data.all.day, "DD10", "All data")

data.all.cell <- subset(data.relevant, select=grep("Raw", grep("^gIs", colnames(data), value=TRUE, invert=TRUE), value=TRUE, invert=TRUE)) #drop flags & Raw columns
data.all.cell <- data.all.cell %>% spread(key=Cell.Line, value=Normalized)
data.all.cell <- data.all.cell %>% gather(key=Cell.Line, value=Normalized, Bhlhb4:Isl1)
plotscatter(data.all.cell, "wt", "All data")

#function to extract entries conteining keyword in description
scatter.plot.key <- function(keyword){
  key <- deparse(substitute(keyword))
  data.subset <- subset(data.relevant, grepl(key, data.relevant$Description))
  data.subset <- subset(data.subset, select=grep("Raw", grep("^gIs", colnames(data.subset), value=TRUE, invert=TRUE), value=TRUE, invert=TRUE)) #drop flags & Raw columns
  data.subset <- data.subset %>% spread(key=dif.day, value=Normalized)
  data.subset <- data.subset %>% gather(key=dif.day, value=Normalized, DD16:DD23)
  plotscatter(data.subset, "DD10", as.character(key))
  
  data.subset <- subset(data.relevant, grepl(key, data.relevant$Description))
  data.subset <- subset(data.subset, select=grep("Raw", grep("^gIs", colnames(data.subset), value=TRUE, invert=TRUE), value=TRUE, invert=TRUE)) #drop flags & Raw columns
  data.subset <- data.subset %>% spread(key=Cell.Line, value=Normalized)
  data.subset <- data.subset %>% gather(key=Cell.Line, value=Normalized, Bhlhb4:Isl1)
  plotscatter(data.subset, "wt", as.character(key))
}

#subset specific genes
keygenes <- subset(data.relevant, grepl("G protein", data.relevant$Description))
keygenes <- subset(keygenes, select=grep("Raw", grep("^gIs", colnames(keygenes), value=TRUE, invert=TRUE), value=TRUE, invert=TRUE)) #drop flags & Raw columns
keygenes <- keygenes %>% spread(key=dif.day, value=Normalized)
keygenes <- keygenes %>% gather(key=dif.day, value=Normalized, DD16:DD23)
plotscatter(keygenes, "DD10", "G protein")

keygenes <- subset(data.relevant, grepl("G protein", data.relevant$Description))
keygenes <- subset(keygenes, select=grep("Raw", grep("^gIs", colnames(keygenes), value=TRUE, invert=TRUE), value=TRUE, invert=TRUE)) #drop flags & Raw columns
keygenes <- keygenes %>% spread(key=Cell.Line, value=Normalized)
keygenes <- keygenes %>% gather(key=Cell.Line, value=Normalized, Bhlhb4:Isl1)
plotscatter(keygenes, "wt", "G protein")
