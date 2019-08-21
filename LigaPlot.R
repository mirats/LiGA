# Clear plots
if(!is.null(dev.list())) dev.off()
# Clear console
cat("\014") 
# Clean workspace
rm(list=ls())


library(tidyverse)
library(readxl)
library(edgeR)
library(RColorBrewer)
library(scales)
library(cowplot)

campaignName<- "DCSIGNcells_vs_YZ" #Do not put extention. File should be saves in .xlsx format ONLY.
dirFiles <- "~/OneDrive - ualberta.ca/DerdaLab/LiGA/LiGA/" #Folder where all the .txt liga files are stored. No need to put name of files
dirCampaign <- "~/Dropbox/Database/Campaign/" #Folder where the xlsx file is stored with information about experiment. Look at Dropbox/Database/Campaign/CD22_vs_YZ.xlsx for examples. 
dirSave<- dirCampaign #Folder where all the images will be saved. 
dirMaldi<- "~/Dropbox/Database/" #Place where MALDI file is stored. Default is Dropbox/Database/
x_axis<- 3 ##NOT WORKING. #Define the X-axis. Options: Mod, Glytoucan, IUPAC.
#-------------------------------------------####################-----------------------------------------------------------------
###Do not change anything beyond this point--------------------------------------------------------------------------------------
setwd(dirSave)
#load the campaign file.
fileC<-read_excel(paste0(dirCampaign, campaignName, ".xlsx", sep=""))  #reads the Campaign file
#load the test data--------------------------------------------------------------------------------------------
testFiles<-fileC$Test
testFiles<-testFiles[!is.na(testFiles)]
test<-lapply(paste0(dirFiles, testFiles, ".txt", sep=""), function(x) 
  read.table(x, header = TRUE, stringsAsFactors = FALSE,fill = TRUE))
test2<-lapply(test, function(x) x[!(names(x) %in% c("index", "mindex", "Primer","Nuc", "AA"))]) #drops the uncessary
for(i in 1:length(test2)){
  colnames(test2[[i]]) <- c("Mod",rep(paste0(substr(testFiles[i], 1, 20), 
                                             1:(sapply(test2[i], NCOL)-1)[1])))
}
test2<- test2 %>% 
  reduce(full_join)
test2<-aggregate(. ~Mod, test2, sum)
test2<-test2[!test2$Mod == "XX", ]
colnames(test2) <- c("Mod", rep(paste0("test", 1:(NCOL(test2)-1))))
#Load the control data in the environment--------------------------------------------------------------------------------------
controlFiles<-fileC$Control
controlFiles<-controlFiles[!is.na(controlFiles)]
control<-lapply(paste0(dirFiles, fileC$Control, ".txt", sep=""), function(x) 
  read.table(x, header = TRUE, stringsAsFactors = FALSE,fill = TRUE))   #reads the Campaign file
control2<-lapply(control, function(x) x[!(names(x) %in% c("index", "mindex", "Primer","Nuc", "AA"))]) #drops the uncessary stuff
for(i in 1:length(control2)){
  colnames(control2[[i]]) <- c("Mod",rep(paste0(substr(controlFiles[i], 1, 20), 
                                                1:(sapply(control2[i], NCOL)-1)[1])))
}
control2<- control2 %>% 
  reduce(full_join)
control2<-aggregate(. ~Mod, control2, sum)
control2<-control2[!control2$Mod == "XX", ]
colnames(control2) <- c("Mod", rep(paste0("control", 1:(NCOL(control2)-1))))
#load the naive data in the environment--------------------------------------------------------------------------------------
naiveFiles<-fileC$Naïve
naiveFiles<-naiveFiles[!is.na(naiveFiles)]
naive<- lapply(paste0(dirFiles, naiveFiles, ".txt", sep=""), function(x)
  read.table(x, header=T, stringsAsFactors = FALSE,fill = TRUE))
naive2<- lapply(naive, function(x) x[!(names(x) %in% c("index", "mindex", "Primer","Nuc", "AA"))])
for(i in 1:length(naive2)){
  colnames(naive2[[i]]) <- c("Mod",rep(paste0(substr(naiveFiles[i], 1, 20), 
                                              1:(sapply(naive2[i], NCOL)-1)[1])))
}
naive2<- naive2 %>% 
  reduce(full_join)
naive2<-aggregate(. ~Mod, naive2, sum)
naive2<-naive2[!naive2$Mod == "XX", ]
colnames(naive2) <- c("Mod", rep(paste0("naive", 1:(NCOL(naive2)-1))))
listAll<-list(test2, control2, naive2)
mergedData <- Reduce(function(x, y) full_join(x, y, by="Mod"), listAll)
##This is where TMM happens. Need to put a switch to turn it on and off. ------------------------------------------------
pepdat <- sapply(mergedData[,-c(1)], as.numeric)
rownames(pepdat) <- mergedData[,1]
pepdat[is.na(pepdat)] <- 0
barplot(colSums(pepdat), col = c(rep("grey50", 4), rep("grey90", 4)), 
        ylab = "Library sizes", main="")
colnam<-colnames(pepdat)
for (i in 1:NCOL(pepdat)) {
  colnam[i]<-substr((colnames(pepdat)[i]), 1, (nchar(colnames(pepdat)[i])-1))
}
cond<-colnam
cond
design <- model.matrix(~0+cond)
design
DGE <- DGEList(counts=pepdat, group=cond)
DGE
DGE <- calcNormFactors(DGE)
DGE$samples
tmm <- cpm(DGE,normalized.lib.sizes = TRUE)
tmm
c=as.data.frame(tmm)
c$Mod<-row.names(tmm)
rownames(c) <- NULL
mergedDataNorm<-c[c("Mod", setdiff(names(c), "Mod"))]
#End of TMM analysis. --------------------------------------------------------------------------------------------------
mergedDataNorm<-mergedDataNorm[!mergedDataNorm$Mod == "???", ]
testAvg<- apply(mergedDataNorm[,2:NCOL(test2)], 1, mean)
controlAvg<- apply(mergedDataNorm[,(tail(2:NCOL(test2), n=1)+1):(tail(2:NCOL(test2), n=1)+NCOL(control2)-1)], 1, mean)
naiveAvg<- apply(mergedDataNorm[,(tail(2:NCOL(test2), n=1)+NCOL(control2)):(NCOL(mergedDataNorm))], 1, mean)
dataT<- data.frame(glycan=mergedDataNorm[1],testAvg,controlAvg, naiveAvg)
longdataT<- mergedDataNorm %>%
  gather(Sample, Freq, colnames(mergedDataNorm[2:ncol(mergedDataNorm)]))
jitter <- position_jitter(width = 0.2, height = 0.2)
### Load the order of the plotting------------------------------------------------------------------------------------------
liga<-read_excel(paste0(dirMaldi, substr(fileC$Test[1], 12, 13), ".xlsx", sep=""), col_names=T, skip=1)
longdataT$Order <- liga$Order[match(longdataT$Mod, 
                                    liga$Alphanum.)]

ligaFile<-read_excel(paste0(dirMaldi, "MALDI-names.xlsx"), col_names=T, skip=0)
ligaFile = ligaFile[-1,]
head(ligaFile)
longdataT$IUPAC<- ligaFile$IUPAC[match(longdataT$Mod, 
                                       ligaFile$`Glycan Name`)]
longdataT$linker<- ligaFile$Linkage[match(longdataT$Mod, 
                                          ligaFile$`Glycan Name`)]
longdataT$Density<- ligaFile$Density[match(longdataT$Mod, 
                                           ligaFile$`Glycan Name`)]
longdataT$GlycanNum<-(longdataT$Density)*27
longdataT$IUPAC<-paste0(longdataT$IUPAC, longdataT$linker, "-[", longdataT$GlycanNum, "]")
longdataT$Glytoucan<-ligaFile$`GlyTouCan ID`[match(longdataT$Mod, 
                                                   ligaFile$`Glycan Name`)]
longdataT$Glytoucan<-paste0(longdataT$Glytoucan,"-[", longdataT$GlycanNum, "]")
longdataT$Mod2<-paste0(longdataT$Mod,"-[", longdataT$GlycanNum, "]")
longdataT[is.na(longdataT)] <- 0

if (x_axis==1) {
  longdataT$x_label<-longdataT$Mod2
} else if (x_axis==2) {
  longdataT$x_label<-longdataT$Glytoucan
} else {
  longdataT$x_label<-longdataT$IUPAC
}


### Load the plotting parameters------------------------------------------------------------------------------------------
jitter <- position_jitter(width = 0.2, height = 0.2)

ColorP<-c(rep("black", NCOL(control2)-1), rep("#999999", NCOL(naive2)-1), rep("black", NCOL(test2)-1))
names(ColorP) <- levels(factor(longdataT$Sample))
colScale <- scale_colour_manual(name =factor(longdataT$Sample),values = ColorP)

ShapeP<- c(rep(25, NCOL(control2)-1), rep(21, NCOL(naive2)-1), rep(24, NCOL(test2)-1))
names(ShapeP) <- levels(factor(longdataT$Sample))
shapeScale <- scale_shape_manual(name =factor(longdataT$Sample),values = ShapeP)

ShapeFill<-c(rep("white", NCOL(control2)-1), rep("#999999", NCOL(naive2)-1), rep("black", NCOL(test2)-1))
names(ShapeFill) <- levels(factor(longdataT$Sample))
fillScale <- scale_fill_manual(name =factor(longdataT$Sample),values = ShapeFill)
### Load the plotting parameters------------------------------------------------------------------------------------------
### Generate Scatterplot1 Map------------------------------------------------------------------------------------------
scatter1<-ggplot(longdataT)+ 
  theme_bw()+
  geom_point(position = jitter, aes(x=reorder(x_label, +Order), y=Freq, color = Sample, shape=Sample,  fill = factor(Sample)), 
             stroke=0.7, size=3)+
  scale_y_log10(limits=c(1, 1e6), labels = trans_format("log10", math_format(10^.x)))+
  colScale+
  shapeScale+
  fillScale+
  labs(y="PPM", x="Glycan")+
  ggtitle(campaignName)+
  theme(axis.text.x=element_text(family="Arial", color="black", angle=90, size=7,hjust=1,vjust=0.2),
        axis.text.y=element_text(family="Arial", color="black",size=7, face="bold"),
        legend.title=element_text(family="Arial", color="black",size=7),
        legend.text=element_text(family="Arial", color="black", size=7),
        title=element_text(family="Arial", color="black", size=12))+
  coord_equal(ratio =2)
scatter1
ggsave(plot = scatter1, width = 25, height = 6.25, dpi = 300, units="in", 
       filename = paste0(campaignName, "-scatter1.eps", sep=""))
ggsave(plot = scatter1, width = 25, height = 6.25, dpi = 300, units="in", 
       filename = paste0(campaignName, "-scatter1.jpg", sep=""))
### Generate Heat Map------------------------------------------------------------------------------------------
jet.colors <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
hm<-ggplot(longdataT, aes(x=reorder(x_label, +Order), y=fct_rev(factor(Sample)))) + 
  theme_light()+
  geom_tile(aes(fill = log10(Freq)), colour = "black", size=0.3) + 
  scale_fill_gradientn(colours = jet.colors(7))+
  ggtitle(campaignName)+
  labs(y="PPM", x="Glycan")+
  theme(axis.text.x=element_text(family="Arial", color="black", angle=90, size=7,hjust=1,vjust=0.2),
        axis.text.y=element_text(family="Arial", color="black",size=7),
        legend.title=element_text(family="Arial", color="black",size=7),
        legend.text=element_text(family="Arial", color="black", size=7),
        title=element_text(family="Arial", color="black", size=12))+
  coord_equal()
hm
ggsave(plot = hm, width = 25, height = 6.25, dpi = 300, units="in", 
       filename = paste0(campaignName, "-hm.eps", sep=""))
ggsave(plot = hm, width = 25, height = 6.25, dpi = 300, units="in", 
       filename = paste0(campaignName, "-hm.jpg", sep=""))
### Generate Enrichment Data-------------------------------------------------------------------------------------
testEN <- testAvg/naiveAvg
controlEN <- controlAvg/naiveAvg
totalEN<-testEN/controlEN

dataEN<- data.frame(glycan=mergedDataNorm[1],testEN,controlEN, totalEN)
is.na(dataEN) <- sapply(dataEN, is.infinite)
dataEN[is.na(dataEN)] <- 0
dataEN$GlycanID<- longdataT$GlycanID[match(dataEN$Mod, 
                                           longdataT$Mod)]
dataEN$Order<- longdataT$Order[match(dataEN$Mod, 
                                     longdataT$Mod)]
dataEN$Mod2<- longdataT$Mod2[match(dataEN$Mod, 
                                   longdataT$Mod)]
if (x_axis==1) {
  dataEN$x_label<-dataEN$Mod2
} else if (x_axis==2) {
  dataEN$x_label<-dataEN$Glytoucan
} else {
  dataEN$x_label<-dataEN$IUPAC
}
###Plotting parameters for scatter 2------------------------------------------------------------------------------------------
scatter2<-ggplot(data=dataEN)+ 
  theme_bw()+
  geom_point(aes(x=reorder(Mod2, +Order),
                 y=testEN),
             stat='identity', size=6, fill="black", 
             color="black", shape=23)+
  geom_point(aes(x=reorder(Mod2, +Order),
                 y=controlEN),
             stat='identity', size=6, 
             color="black", fill="white", shape=23)+
  geom_segment(aes(x=reorder(Mod2, +Order), 
                   xend=Mod2, 
                   y=testEN,
                   yend=controlEN, size=totalEN))+
  scale_y_log10(labels = trans_format("log10", math_format(10^.x)))+
  scale_size_continuous(range = c(0.1, 1.5))+
  labs(y="Enrichment", x="Glycan")+
  ggtitle(campaignName)+
  scale_fill_manual(values=c("#97CAD8","#DC1452"))+
  theme(axis.text.x=element_text(family="Arial", color="black", angle=90, size=7,hjust=1,vjust=0.2),
        axis.text.y=element_text(family="Arial", color="black",size=7, face="bold"),
        legend.title=element_text(family="Arial", color="black",size=7),
        legend.text=element_text(family="Arial", color="black", size=7),
        title=element_text(family="Arial", color="black", size=12))+
  coord_equal(ratio = 4)
scatter2
ggsave(plot = scatter2, width = 25, height = 6.25, dpi = 300, units="in", 
       filename = paste0(campaignName, "-scatter2.eps", sep=""))
ggsave(plot = scatter2, width = 25, height = 6.25, dpi = 300, units="in", 
       filename = paste0(campaignName, "-scatter2.jpg", sep=""))
#Beyond this point is experimentation
