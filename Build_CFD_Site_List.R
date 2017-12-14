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

# Screen out any directories that don't contain all the needed input files.
all.files = rep(T, length(csv_dirs))
for (i in 1:length(csv_dirs)){
all.files[i]=("DEM.csv" %in% dir(csv_dirs[i]) & "WSEDEM.csv" %in% dir(csv_dirs[i])  & "Thalweg.csv" %in% dir(csv_dirs[i]))
}
all.files
csv_dirs = csv_dirs[all.files]


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

#MVI = read.csv("MetricVisitInformation.csv", header=T)
MVI = read.csv("C://Matt-SFR Files/Hydraulic Modeling/R Code to Build Input Files/R-Code/CHaMP_and_AEM_Metrics.csv", header=T)
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
"Trim.Length" = rep(2, length(idx)),
"HEV" = rep(.01, length(idx)),
"DeltaBC" = rep(0, length(idx)))

setwd(saved_dir)
# replace trim lengths with those in file "TrimLengths.csv")
TrimLengths = read.csv("TrimLengths.csv", header=T)
idx=match(CFD_SiteList$SiteID,TrimLengths$SiteName)
CFD_SiteList$Trim.Length= TrimLengths$Trim.Length[idx]
CFD_SiteList$Trim.Length[is.na(CFD_SiteList$Trim.Length)] = 2


#Remove spaces from SiteID
CFD_SiteList$SiteID = gsub(" ","", CFD_SiteList$SiteID)
CFD_SiteList$WatershedName = gsub(" ","", CFD_SiteList$WatershedName)


# Set Model to "No" if any NA's or if discharge = 0
complete.cases(CFD_SiteList)
Model = complete.cases(CFD_SiteList)
Model[CFD_SiteList$Measured.Discharge < .0001] = F
CFD_SiteList$Model = as.character(CFD_SiteList$Model)
CFD_SiteList$Model[Model == F] = "No"
##############################

# Set Discharge to 0 if NA, otherwise build_input_files.R craches
CFD_SiteList$Measured.Discharge[is.na(CFD_SiteList$Measured.Discharge)==T]=0
CFD_SiteList$Modeled.Discharge[is.na(CFD_SiteList$Modeled.Discharge)==T]=0


setwd("C:/Matt-SFR Files/Hydraulic Modeling/R Code to Build Input Files")

write.csv(CFD_SiteList,"CFD_Site_List.csv")
CFD_SiteList
setwd(saved_dir)
