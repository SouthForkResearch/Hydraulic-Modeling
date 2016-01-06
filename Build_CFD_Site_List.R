# Script to find directories with dems and stuff
saved_dir = getwd()
getwd()

# Set the WD to where the CHaMP data is at
#setwd("c://CHaMP/MonitoringData")
setwd("C://Matt-SFR Files/Hydraulic Modeling/champ data from bucket")

# list all the directories with zipped up hydro model data
#dirs = list.files(recursive=T, pattern = "HydroModelInputs")
dirs = list.files(recursive=T, pattern = "WSEDEM.csv")
dirs

# list the dirs where we want to unzip the data into
#hydro.dirs = sub("/HydroModelInputs.zip", "", dirs)
csv_dirs = sub("WSEDEM.csv", "", dirs)
csv_dirs = paste("c://Matt-SFR Files/Hydraulic Modeling/champ data from bucket/", csv_dirs, sep="")
csv_dirs
##unzip the data
#for (i in 1:length(dirs)){
#unzip(dirs[i], exdir=hydro.dirs[i])
#}

## name the dirs where the csv files are
#csv_dirs = paste(hydro.dirs, "/artifacts", sep="")
#dir(csv_dirs[1])

# Do some stuff to extract the VISITID from within the string
VisitIDs = rep("", length(csv_dirs))
i=1
for (i in 1:length(csv_dirs)){
csv_dirs[i]
A=strsplit(csv_dirs[i], "VISIT_")[[1]][2]
VisitIDs[i] = (strsplit(A,"/")[[1]][1])
VisitIDs[i]
}
VisitIDs=as.numeric(VisitIDs)
VisitIDs

# Here's what we've got so far
setwd("C:/Matt-SFR Files/Hydraulic Modeling/R Code to Build Input Files")
dir()

MVI = read.csv("MetricVisitInformation.csv", header=T)
idx = match(VisitIDs, MVI$VisitID)
idx

names(MVI)
CFD_SiteList = data.frame(
"SiteID"=MVI$SiteName[idx],
"Directory" = csv_dirs,
"Model" = rep("Yes", length(idx)),
"Measured.Discharge" = MVI$Q[idx],
"Modeled.Discharge" = MVI$Q[idx],
"D84" = MVI$SubD84[idx],
"Roughness" = MVI$SubD84[idx]*4/1000,
"Year"=MVI$VisitYear[idx],
"WatershedName"=MVI$WatershedName[idx],
"VisitID" = VisitIDs,
"Trim Length" = rep(10, length(idx)),
"HEV" = rep(.01, length(idx)),
"DeltaBC" = rep(0, length(idx)))

#Remove spaces from SiteID
CFD_SiteList$SiteID = gsub(" ","", CFD_SiteList$SiteID)
CFD_SiteList$WatershedName = gsub(" ","", CFD_SiteList$WatershedName)

CFD_SiteList

write.csv(CFD_SiteList,"CFD_Site_List.csv")
CFD_SiteList
setwd(saved_dir)
