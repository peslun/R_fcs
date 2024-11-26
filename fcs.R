#install ----
install.packages("ggplot2")
install.packages("BiocManager")
BiocManager::install("flowCore")
BiocManager::install("flowViz")
BiocManager::install("flowAI")
BiocManager::install("ggcyto")

library(ggplot2)
library(flowCore)
library(flowViz)
library(flowAI)

#load libraries ----
library(extrafont)
library(ggplot2)
library(flowCore)
library(flowViz)
library(flowAI)
library(ggcyto)


#one fcs file ----
#file:///home/pt/app/rstudio/13-09-24_reference_3394_A.fcs
#myfile <- "~/app/rstudio/13-09-24_reference_3394_A.fcs"
#fcsfile <- read.FCS(myfile,emptyValue = F)

#fcsfile
#exprs(fcsfile)[1:10,]
#       Event Time  FSC-A   SSC-A  BL1-A  BL3-A  FSC-H  SSC-H  BL1-H  BL3-H FSC-W SSC-W BL1-W BL3-W
# [1,]     1   22 287762  558395 196066  64875 247590 455184 167625  58306    80    80    68    35

#summary(fcsfile)
#str(keyword(fcsfile))
#summary(fcsfile[,5:6])

#plot(fcsfile, c("BL3-A", "BL1-A"))
#plot(fcsfile[,6:5])


#several fcs files/flowset ----
##load files ----

filesR <- list.files(path="~/app/rstudio/10.facstrack/110113-BSSEDTA",pattern="75.*[[:digit:]]{3}$",full.names=T)

filesB <- list.files(path="~/app/rstudio/10.facstrack/110113-BSSEDTA",pattern="72B.*[[:digit:]]{3}$",full.names=T)
                    #,pattern=".fcs$")
filesE <- list.files(path="~/app/rstudio/10.facstrack/110113-BSSEDTA",pattern="72E.*[[:digit:]]{3}$",full.names=T)

##create flow set ----
#fs <- read.flowSet(c(filesB,filesE))
ref <- read.flowSet(filesR)
bss <- read.flowSet(filesB)
edta <- read.flowSet(filesE)

#assign flow set ----
fs <- bss

#fs
#summary(fs[[1]])
#plot(fs[[1]],c("FL3-H","FL1-H"),smooth=F)

#basic graph ----
plot(fs[[1]],c("FL3-H","FL1-H"),smooth=F)
autoplot(fs,x="FL3-H",y="FL1-H",bins=256)


#gating set ----
gs <- GatingSet(fs)
##sperm gate ----
spermC <- matrix(c(120,100,0,0,1000,1000, 0,300,400,1000,1000,0),ncol=2)
colnames(spermC) <- c("FL3-H","FL1-H")
spermG <- polygonGate(filterId="sperm",.gate=spermC)
#gs_pop_remove(gs,node="root/sperm")
gs_pop_add(gs,spermG,parent="root")
gs_get_pop_paths(gs)
recompute(gs)

spermFS <- gs_pop_get_data(gs,"root/sperm")
spermGS <- GatingSet(spermFS)

#everything with the gate
autoplot(gs,x="FL3-H",y="FL1-H",bins=256,"root/sperm")
#only the gate
autoplot(spermFS,x="FL3-H",y="FL1-H",bins=256)

#DFI gate ----
#run sperm gate before this and create spermFS

dfiC <- matrix(c(0,250,1000,1000, 0,1000,1000,0),ncol=2)
colnames(dfiC) <- c("FL3-H","FL1-H")
dfiG <- polygonGate(filterId="dfi",.gate=dfiC)
#gs_pop_remove(spermGS,node="root/dfi")
gs_pop_add(spermGS,dfiG,parent="root")
gs_get_pop_paths(spermGS)
recompute(spermGS)

#everything with the gate
autoplot(spermGS,x="FL3-H",y="FL1-H",bins=256,"root/dfi")

#HGS gate ----
#run sperm gate before this and create spermGS

hgsC <- matrix(c(0,0,1000,1000, 750,1000,1000,750),ncol=2)
colnames(hgsC) <- c("FL3-H","FL1-H")
hgsG <- polygonGate(filterId="hgs",.gate=hgsC)
#gs_pop_remove(spermGS,node="root/dfi")
gs_pop_add(spermGS,hgsG,parent="root")
gs_get_pop_paths(spermGS)
recompute(spermGS)

#everything with the gate
autoplot(spermGS,x="FL3-H",y="FL1-H",bins=256,gate=c("root/dfi","root/hgs"))


#final display ----
autoplot(spermGS,x="FL3-H",y="FL1-H",bins=256,gate=c("root/dfi","root/hgs")) +
  ylab ("AO green") + xlab ("AO red") +
#  theme(panel.background=element_rect(fill="white")) + 
  theme_classic() +
  theme(strip.text=element_text(size=9,face="bold"),axis.text.x=element_blank(),axis.text.y = element_blank(), 
        axis.ticks.x = element_blank(), axis.ticks.y = element_blank(), legend.position="none",
        plot.title=element_blank()
        )
#export image ----
##pipeline ----
files <- list.files(path="~/app/rstudio/10.facstrack/110120-BSSEDTA",
                    pattern="52.*[[:digit:]]{3}$",full.names=T)
gs <- GatingSet(read.flowSet(files))
spermC <- matrix(c(120,100,0,0,1000,1000, 0,300,400,1000,1000,0),ncol=2)
colnames(spermC) <- c("FL3-H","FL1-H")
spermG <- polygonGate(filterId="sperm",.gate=spermC)
gs_pop_add(gs,spermG,parent="root")
gs_get_pop_paths(gs)
recompute(gs)
spermGS <- GatingSet(gs_pop_get_data(gs,"root/sperm"))

dfiC <- matrix(c(0,250,1000,1000, 0,1000,1000,0),ncol=2)
colnames(dfiC) <- c("FL3-H","FL1-H")
dfiG <- polygonGate(filterId="dfi",.gate=dfiC)
gs_pop_add(spermGS,dfiG,parent="root")

hgsC <- matrix(c(0,0,1000,1000, 750,1000,1000,750),ncol=2)
colnames(hgsC) <- c("FL3-H","FL1-H")
hgsG <- polygonGate(filterId="hgs",.gate=hgsC)
#gs_pop_remove(spermGS,node="root/dfi")
gs_pop_add(spermGS,hgsG,parent="root")
gs_get_pop_paths(spermGS)
recompute(spermGS)

###loadfonts(device = "postscript", quiet = TRUE)

setEPS()
postscript("export.eps",onefile = FALSE, horizontal = FALSE, paper= "special")
autoplot(spermGS,x="FL3-H",y="FL1-H",bins=256,gate=c("root/dfi","root/hgs")) +
  ylab ("AO green") + xlab ("AO red") +
  #  theme(panel.background=element_rect(fill="white")) + 
  theme_classic() +
  theme(strip.text=element_text(size=9,face="bold"),axis.text.x=element_blank(),axis.text.y = element_blank(), 
        axis.ticks.x = element_blank(), axis.ticks.y = element_blank(), legend.position="none",
        plot.title=element_blank()
  )
#This actually save the plot in a image
dev.off()

