#load micro array data

source("gen-log.R")
#log results in file
logfile <- gen.start.log("Micro Array Data")
gen.log("Load excel data", logfile, datetime=TRUE)

#load raw data *convert excel file to csv
#data.raw <- read.csv("import test.csv", skip=1, header=TRUE)
data.raw <- read.csv("Individual Sample Data.csv", skip=1, header=TRUE)

#initialize dataframe to store data
data <- data.frame(FeatureNum=double(),
                   ProbeName=character(),
                   GeneName=character(),
                   SystematicName=character(),
                   accessions=character(),
                   chr_coord=character(),
                   Start=double(),
                   Description=character(),
                   ControlType=logical(),
                   Normalized=double(),
                   Raw=double(),
                   gIsSaturated=logical(),
                   gIsFeatNonUnifOL=logical(),
                   gIsBGNonUnifOL=logical(),
                   gIsFeatPopnOL=logical(),
                   gIsBGPopnOL=logical(),
                   gIsPosAndSignif=logical(),
                   gIsWellAboveBG=logical(),
                   Cell.line=character(),
                   dif.day=character())

#load wt dd10
data <- rbind(data, data.frame(
  FeatureNum=data.raw$FeatureNum, 
  ProbeName=data.raw$ProbeName, 
  GeneName=data.raw$GeneName, 
  SystematicName=data.raw$SystematicName, 
  accessions=data.raw$accessions,
  chr_coord=data.raw$chr_coord,
  Start=data.raw$Start,
  Description=data.raw$Description,
  ControlType=as.logical(data.raw$ControlType),
  Normalized=data.raw[,22],
  Raw=data.raw[,23],
  gIsSaturated=as.logical(data.raw[,70]),
  gIsFeatNonUnifOL=as.logical(data.raw[,71]),
  gIsBGNonUnifOL=as.logical(data.raw[,72]),
  gIsFeatPopnOL=as.logical(data.raw[,73]),
  gIsBGPopnOL=as.logical(data.raw[,74]),
  gIsPosAndSignif=as.logical(data.raw[,75]),
  gIsWellAboveBG=as.logical(data.raw[,76]),
  Cell.Line="wt",
  dif.day="DD10"))

#load wt dd16
data <- rbind(data, data.frame(
  FeatureNum=data.raw$FeatureNum, 
  ProbeName=data.raw$ProbeName, 
  GeneName=data.raw$GeneName, 
  SystematicName=data.raw$SystematicName, 
  accessions=data.raw$accessions,
  chr_coord=data.raw$chr_coord,
  Start=data.raw$Start,
  Description=data.raw$Description,
  ControlType=as.logical(data.raw$ControlType),
  Normalized=data.raw[,24],
  Raw=data.raw[,25],
  gIsSaturated=as.logical(data.raw[,77]),
  gIsFeatNonUnifOL=as.logical(data.raw[,78]),
  gIsBGNonUnifOL=as.logical(data.raw[,79]),
  gIsFeatPopnOL=as.logical(data.raw[,80]),
  gIsBGPopnOL=as.logical(data.raw[,81]),
  gIsPosAndSignif=as.logical(data.raw[,82]),
  gIsWellAboveBG=as.logical(data.raw[,83]),
  Cell.Line="wt",
  dif.day="DD16"))

#load wt dd23
data <- rbind(data, data.frame(
  FeatureNum=data.raw$FeatureNum, 
  ProbeName=data.raw$ProbeName, 
  GeneName=data.raw$GeneName, 
  SystematicName=data.raw$SystematicName, 
  accessions=data.raw$accessions,
  chr_coord=data.raw$chr_coord,
  Start=data.raw$Start,
  Description=data.raw$Description,
  ControlType=as.logical(data.raw$ControlType),
  Normalized=data.raw[,26],
  Raw=data.raw[,27],
  gIsSaturated=as.logical(data.raw[,84]),
  gIsFeatNonUnifOL=as.logical(data.raw[,85]),
  gIsBGNonUnifOL=as.logical(data.raw[,86]),
  gIsFeatPopnOL=as.logical(data.raw[,87]),
  gIsBGPopnOL=as.logical(data.raw[,88]),
  gIsPosAndSignif=as.logical(data.raw[,89]),
  gIsWellAboveBG=as.logical(data.raw[,90]),
  Cell.Line="wt",
  dif.day="DD23"))

#load B4KO DD10
data <- rbind(data, data.frame(
              FeatureNum=data.raw$FeatureNum, 
              ProbeName=data.raw$ProbeName, 
              GeneName=data.raw$GeneName, 
              SystematicName=data.raw$SystematicName, 
              accessions=data.raw$accessions,
              chr_coord=data.raw$chr_coord,
              Start=data.raw$Start,
              Description=data.raw$Description,
              ControlType=as.logical(data.raw$ControlType),
              Normalized=data.raw[,10],
              Raw=data.raw[,11],
              gIsSaturated=as.logical(data.raw[,28]),
              gIsFeatNonUnifOL=as.logical(data.raw[,29]),
              gIsBGNonUnifOL=as.logical(data.raw[,30]),
              gIsFeatPopnOL=as.logical(data.raw[,31]),
              gIsBGPopnOL=as.logical(data.raw[,32]),
              gIsPosAndSignif=as.logical(data.raw[,33]),
              gIsWellAboveBG=as.logical(data.raw[,34]),
              Cell.Line="Bhlhb4",
              dif.day="DD10"))

#load B4KO DD16
data <- rbind(data, data.frame(
  FeatureNum=data.raw$FeatureNum, 
  ProbeName=data.raw$ProbeName, 
  GeneName=data.raw$GeneName, 
  SystematicName=data.raw$SystematicName, 
  accessions=data.raw$accessions,
  chr_coord=data.raw$chr_coord,
  Start=data.raw$Start,
  Description=data.raw$Description,
  ControlType=as.logical(data.raw$ControlType),
  Normalized=data.raw[,12],
  Raw=data.raw[,13],
  gIsSaturated=as.logical(data.raw[,35]),
  gIsFeatNonUnifOL=as.logical(data.raw[,36]),
  gIsBGNonUnifOL=as.logical(data.raw[,37]),
  gIsFeatPopnOL=as.logical(data.raw[,38]),
  gIsBGPopnOL=as.logical(data.raw[,39]),
  gIsPosAndSignif=as.logical(data.raw[,40]),
  gIsWellAboveBG=as.logical(data.raw[,41]),
  Cell.Line="Bhlhb4",
  dif.day="DD16"))

#load B4KO DD23
data <- rbind(data, data.frame(
  FeatureNum=data.raw$FeatureNum, 
  ProbeName=data.raw$ProbeName, 
  GeneName=data.raw$GeneName, 
  SystematicName=data.raw$SystematicName, 
  accessions=data.raw$accessions,
  chr_coord=data.raw$chr_coord,
  Start=data.raw$Start,
  Description=data.raw$Description,
  ControlType=as.logical(data.raw$ControlType),
  Normalized=data.raw[,14],
  Raw=data.raw[,15],
  gIsSaturated=as.logical(data.raw[,42]),
  gIsFeatNonUnifOL=as.logical(data.raw[,43]),
  gIsBGNonUnifOL=as.logical(data.raw[,44]),
  gIsFeatPopnOL=as.logical(data.raw[,45]),
  gIsBGPopnOL=as.logical(data.raw[,46]),
  gIsPosAndSignif=as.logical(data.raw[,47]),
  gIsWellAboveBG=as.logical(data.raw[,48]),
  Cell.Line="Bhlhb4",
  dif.day="DD23"))

#load Isl1KO dd10
data <- rbind(data, data.frame(
  FeatureNum=data.raw$FeatureNum, 
  ProbeName=data.raw$ProbeName, 
  GeneName=data.raw$GeneName, 
  SystematicName=data.raw$SystematicName, 
  accessions=data.raw$accessions,
  chr_coord=data.raw$chr_coord,
  Start=data.raw$Start,
  Description=data.raw$Description,
  ControlType=as.logical(data.raw$ControlType),
  Normalized=data.raw[,16],
  Raw=data.raw[,17],
  gIsSaturated=as.logical(data.raw[,49]),
  gIsFeatNonUnifOL=as.logical(data.raw[,50]),
  gIsBGNonUnifOL=as.logical(data.raw[,51]),
  gIsFeatPopnOL=as.logical(data.raw[,52]),
  gIsBGPopnOL=as.logical(data.raw[,53]),
  gIsPosAndSignif=as.logical(data.raw[,54]),
  gIsWellAboveBG=as.logical(data.raw[,55]),
  Cell.Line="Isl1",
  dif.day="DD10"))

#load Isl1KO dd16
data <- rbind(data, data.frame(
  FeatureNum=data.raw$FeatureNum, 
  ProbeName=data.raw$ProbeName, 
  GeneName=data.raw$GeneName, 
  SystematicName=data.raw$SystematicName, 
  accessions=data.raw$accessions,
  chr_coord=data.raw$chr_coord,
  Start=data.raw$Start,
  Description=data.raw$Description,
  ControlType=as.logical(data.raw$ControlType),
  Normalized=data.raw[,18],
  Raw=data.raw[,19],
  gIsSaturated=as.logical(data.raw[,56]),
  gIsFeatNonUnifOL=as.logical(data.raw[,57]),
  gIsBGNonUnifOL=as.logical(data.raw[,58]),
  gIsFeatPopnOL=as.logical(data.raw[,59]),
  gIsBGPopnOL=as.logical(data.raw[,60]),
  gIsPosAndSignif=as.logical(data.raw[,61]),
  gIsWellAboveBG=as.logical(data.raw[,62]),
  Cell.Line="Isl1",
  dif.day="DD16"))

#load Isl1KO dd23
data <- rbind(data, data.frame(
  FeatureNum=data.raw$FeatureNum, 
  ProbeName=data.raw$ProbeName, 
  GeneName=data.raw$GeneName, 
  SystematicName=data.raw$SystematicName, 
  accessions=data.raw$accessions,
  chr_coord=data.raw$chr_coord,
  Start=data.raw$Start,
  Description=data.raw$Description,
  ControlType=as.logical(data.raw$ControlType),
  Normalized=data.raw[,20],
  Raw=data.raw[,21],
  gIsSaturated=as.logical(data.raw[,63]),
  gIsFeatNonUnifOL=as.logical(data.raw[,64]),
  gIsBGNonUnifOL=as.logical(data.raw[,65]),
  gIsFeatPopnOL=as.logical(data.raw[,66]),
  gIsBGPopnOL=as.logical(data.raw[,67]),
  gIsPosAndSignif=as.logical(data.raw[,68]),
  gIsWellAboveBG=as.logical(data.raw[,69]),
  Cell.Line="Isl1",
  dif.day="DD23"))



gen.log("Data laoded", logfile, datetime=TRUE)
gen.log(paste("obs.", dim(data)[1] ), logfile)
gen.log(paste("variables", dim(data)[2] ), logfile)

#load annotation data from eArray 
data.annotations <- read.csv("AllAnnotations_074809_D_AA_20150624.csv", header=TRUE)

#Add GO id 
data$GO <- data.annotations$GO[match(data$ProbeName, data.annotations$ProbeID)]

save(data,file="data.Rda")

