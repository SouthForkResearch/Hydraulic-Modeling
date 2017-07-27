# CHaMP Hydraulic Modeling Post-Processing R-Code
# Author: Matt Nahorniak, South Fork Research
# matt@southforkresearch.org
#####################################################33333



Post.Processing.R.Code.Version = 1.1
# Version 1.1 used output from vs.exe rather than Quickplot. Results are exactly the same, except
# for some imperceptable roundoff error in the 4th or 5th or more decimal.

# Set this so scientific notation is not used when writing output files
options(scipen=10)

args = commandArgs(trailingOnly=TRUE)
scriptdir = args[1]
visitdir = args[2]
cfdfile = args[3]
setwd(scriptdir)

library(raster)

############################################################################
savedwd = getwd()
savedwd
# R-Code to Post-Process Delft 3D for CHaMP Hydraulic Modeling
# Take text files created after running Delft 3D flow on CHaMP Sites
# (and running quickplot macro to convert Flow results into text files)
# and transfer those results back onto original DEM Grid, creating
# .csv files at those X-Y points for all key outputs.

######################################################################

# Define "Working Dir" - directory with folders for each site DEM, WSEDEM and Thalweg .csv files
# All results files stored in a single folder.  Assign folder name here.

dir.create("Delft3d Results Files", showWarnings=F)

##########################################################################


# M and N for curvilinear Grid Results to be used later
M=1200
N=51



options(digits=16)
# windows(record=T)
library(RANN)
library(uuid)
#library(R.matlab)

# Assign color palette for plots
colorp= colorRampPalette(c("dark blue","blue","cyan","green","yellow","orange", "red","dark red"))(51) ## (n)


# Set the working directory to read site list info
site.list=read.csv(cfdfile)


#################
# Re-Create directory where Delft3D results files are  stored
# This is needed since we may run multiple discharge variants from the same
# DEM, WSEDEM, and thalweg data files
# Duplicate code from Build_Input_Files.R


site.list$D3D.Input.Folder = rep("", nrow(site.list))
Variant = rep("", nrow(site.list))
k=1
for (k in 1:nrow(site.list)){

# Need a "variant string" 6 characters long, plus one for the leading letter.
Variant.str = as.character(signif(round(site.list$Modeled.Discharge[k],4),4))
lead.zeroes =rep("0", 6-nchar(Variant.str))
lead.zeroes = paste(lead.zeroes, collapse="", sep="")
Variant.str = paste(lead.zeroes, Variant.str,"/", sep="")

Variant.str
Variant.str=gsub(".", "_", Variant.str, fixed=T)
Variant.str=paste("000", Variant.str, sep="")
Variant[k] = Variant.str

if (site.list$Measured.Discharge[k] == site.list$Modeled.Discharge[k]){
site.list$D3D.Input.Folder[k]=paste(site.list$Directory[k],"S", Variant.str, sep="")} else{
site.list$D3D.Input.Folder[k]=paste(site.list$Directory[k],"M", Variant.str, sep="")} 
dir.create(site.list$D3D.Input.Folder[k])

}

site.list$D3D.Input.Folder
names(site.list)

for (k in 1:nrow(site.list)) { 
if (site.list$Measured.Discharge[k] == site.list$Modeled.Discharge[k]){
  Variant[k] = paste("S", Variant[k], sep="")} else {
  Variant[k] = paste("M", Variant[k], sep="")}
}

site.list$Results.Folder = 
  paste(visitdir, "/Hydro/Results/",Variant, sep="")



##################
#################






# match sites in dite_list file to sites in directory
#site.index = match(dir(WorkingDir),site.list$SiteID)
#site.index
# Create a site index for the sites that are to be modeled.
site.index = 1:nrow(site.list)
site.index = (1:nrow(site.list))[site.list$Model == "Yes"]


##################################################################
# Note that this will post-process ALL sites in the specified directory.  To run just
# a subset of sites, comment out the "for (k in site.index) {" line and re-rewite it to inly
# run the index list you want to run.  

#### Loop Through Sites
k = site.index[1]
site.index
for (k in site.index) {


print(paste("k=",k))

# Move into directory for current site to post-process

# Set the working directory
WorkingDir = as.character(site.list$D3D.Input.Folder[k])
setwd(paste(WorkingDir, "/", sep=""))



# Read the meta.data file
Meta.Data = read.csv(paste(WorkingDir,"/Meta.Data.csv", sep=""))
Meta.Data$Post.Processing.R.Code.Version = Post.Processing.R.Code.Version
Meta.Data$Post.Processing.Date.Time = Sys.time()

names(site.list)
# Read the DEM
# Set the working directory
newwd = paste(visitdir,"Hydro/HydroHydroModelInputs", sep="")
setwd(newwd)

data =read.csv("DEM.csv", header=F)


# Re-read if data has headers (i.e. new format from Matt R's C++ code)
if (is.numeric(data[,1])==F){
data= read.csv("DEM.csv", header=T)}

names(data) = c("X", "Y", "grid_code")

min(data$grid_code)
min(data$grid_code[data$grid_code > -9999])

print(paste(site.list$VisitID[k], site.list$SiteID[k], min(data$grid_code[data$grid_code > -9999])))

# Read the WSEDEM
WSEDEM = read.csv("WSEDEM.csv", header=F)

# Re-read if data has headers (i.e. new format from Matt R's C++ code)
if (is.numeric(WSEDEM[,1])==F)
{WSEDEM= read.csv("WSEDEM.csv", header=T)}

WSEDEM= WSEDEM[WSEDEM[,3] > -9999,]


options(digits=12)

# Set the working directory back to where the Delft3D results are stored.
WorkingDir = as.character(site.list$D3D.Input.Folder[k])
setwd(paste(, "/", sep=""))

#################################
#######################################################################
#######################################################################
# Read data files output from vs

Uvel=as.numeric(scan(file= "U_vel.tkl", skip = 4, what="numeric"))
Vvel=as.numeric(scan(file= "V_Vel.tkl", skip = 4, what="numeric"))
Xcoord=as.numeric(scan(file= "X_Coor.tkl", skip = 4, what="numeric"))
Ycoord=as.numeric(scan(file= "Y_Coor.tkl", skip = 4, what="numeric"))
BotDepth=as.numeric(scan(file= "Bot_Depth.tkl", skip = 4, what="numeric"))
WaterLevel=as.numeric(scan(file= "Water_Level.tkl", skip = 4, what="numeric"))
Active=as.numeric(scan(file= "Active.tkl", skip = 4, what="numeric"))

# Get M and N (actually M+1 and N+1)
dim = scan(file = "X_Coor.tkl", n=2, skip=3)
M=dim[1]
N=dim[2]

# Put them into arrays.  Need for interpolating Xvel and Yvel onto depth and LW grid points
Uvel = t(array(Uvel, c(N,M)))
Vvel = t(array(Vvel, c(N,M)))
Xcoord = t(array(Xcoord, c(N,M)))
Ycoord = t(array(Ycoord, c(N,M)))
Active = t(array(Active, c(N,M)))
WaterLevel = t(array(WaterLevel, c(N,M)))
BotDepth = t(array(BotDepth, c(N,M)))


# Now interplate for X and Y velocity components.
# Xvel is defined at a half cell above and below the
# depth point.  Yvel is defined at a half cell left and right of the velocity
# point.  Assume off the edge is zero velocity.

Y_vel = Vvel
Y_vel[1,1:N] = 0.5*Vvel[1,]
for (m in 2:M) {Y_vel[m,] = 0.5*(Vvel[m-1,]+Vvel[m,])}

X_vel = 0*Uvel
X_vel[,1] = 0.5*Uvel[,1]
for (n in 2:N) {X_vel[,n] = 0.5*(Uvel[,n-1]+Uvel[,n])}

# Done interpolating Xvel and Yvel, still in MxN array
#########################


# Move everything to a linear vector, for use in the R-Code.
# This is consistent with how the old version of the code worked when
# using output from Quickplot.
# Initialize vectors
X=0
Y=0
Xvel=0
Yvel=0
WL=0
BD=0
AP=0

# This re-assigns arrays into vectors
idx=1
for (n in 1:N) {
     X[idx:(idx+M-1)] = Xcoord[1:M,n]
     Y[idx:(idx+M-1)] = Ycoord[1:M,n]
     Xvel[idx:(idx+M-1)] = X_vel[1:M,n]
     Yvel[idx:(idx+M-1)] = Y_vel[1:M,n]   
     WL[idx:(idx+M-1)] = WaterLevel[1:M,n]
     BD[idx:(idx+M-1)] = BotDepth[1:M,n]
     AP[idx:(idx+M-1)] = Active[1:M,n]
   idx = idx+M
    }

# Filter everything to only active points
X= X[AP==1]
Y= Y[AP==1]
Xvel=Xvel[AP==1]
Yvel=Yvel[AP==1]
WL=WL[AP==1]
BD=BD[AP==1]

# Plot velocity (check)
#Vmag = sqrt(Xvel^2 + Yvel^2)
#col.idx = 1 + trunc(9*Vmag/max(Vmag))
#color=rainbow(10)[col.idx]
#color[Vmag==0] = "black"
#plot(X, Y, col=color)


# Calculate depth and velocity magnitude
Depth = abs(-WL-BD)

#Vmag = sqrt(Xvel*Xvel + Yvel*Yvel)

# Assign data frames just as used to be generated when reading files generated
# by Quickplot
dataVel = data.frame("X"=X,"Y"=Y, "XVel"=Xvel, "YVel"=Yvel)
dataWaterLevel = data.frame("X"=X[Depth!=0], "Y"=Y[Depth!=0], "Z" = WL[Depth!=0])
dataDepth = data.frame("X"=X, "Y"=Y, "Depth"=Depth)

##########################################################
##########################################################
##########################################################


##################################
if (1==2) {
# Everything below to be replaced with above #############3
# Read the velocity outputs
dataVel = read.table("depth averaged velocity.xyz", header=T)
names(dataVel) = c("X", "Y", "XVel", "YVel")
nrow(dataVel)

# Read the water level outputs
dataWaterLevel = read.table("water level.xyz", header=T)
names(dataWaterLevel) = c("X", "Y", "Z")
nrow(dataWaterLevel)


# Read the depth outputs
dataDepth = read.table("water depth.xyz", header=T)
names(dataDepth) = c("X", "Y","Depth")
nrow(dataDepth)

# Read teh vorticity outputs
#dataVorticity = read.table("vorticity.xyz", header=T)
#names(dataVorticity)  = c("X", "Y", "Z", "Vorticity")
dataVorticity = dataVel[,1:3]
dataVorticity[,3] = rep(0, nrow(dataVorticity))
colnames(dataVorticity) = c("X", "Y", "Vorticity")



# Read the bed shear stress outputs
if ("bed shear stress.xyz" %in% dir()){
dataBedShear = read.table("bed shear stress.xyz", header=T)} else {
dataBedShear = read.table("bead shear stress.xyz", header=T)}
names(dataBedShear) = c("X", "Y", "Bed_Shear_X", "Bed_Shear_Y")
nrow(dataBedShear)
###########################################################
# End replaced section
###########################################################
} # end of "if (1==2)" to skip this section


# Read the offset file created while generating input files.
# this file containing the offsets between true x and y locations
# and x and y locations used in the model.  Because the X and Y points
# have too many digits, I can't pass the entire X and Y locations
# to Delft 3D, or the last digits will be lost, and chaos ensues.  A
# pain in the butt, but that's how it needs to be for now....

offset = read.csv("offset.csv", header=T)

names(dataWaterLevel)
# Correct X and Y for the offset
dataDepth$X = dataDepth$X+offset$X.offset
dataDepth$Y = dataDepth$Y+offset$Y.offset
dataVel$X = dataVel$X+offset$X.offset
dataVel$Y = dataVel$Y+offset$Y.offset
dataWaterLevel$X = dataWaterLevel$X + offset$X.offset
dataWaterLevel$Y = dataWaterLevel$Y + offset$Y.offset
#dataVorticity$X = dataVorticity$X + offset$X.offset
#dataVorticity$Y = dataVorticity$Y + offset$Y.offset
#dataBedShear$X = dataBedShear$X + offset$X.offset
#dataBedShear$Y =dataBedShear$Y + offset$Y.offset

nrow(dataDepth)
nrow(dataWaterLevel)

###################################################
# Now move to the output folder


setwd(site.list$Results.Folder[k])


###################################################################
#
# Check modeled Depth vs WSEDEM depth and modeled WSE vs WSEDEM
# Optional - plotting functions commented out for now.


WSE_check_NN = nn2(WSEDEM[,1:2],dataDepth[,1:2],1)
DEM_check_NN = nn2(data[,1:2],dataDepth[,1:2],1)


Depth_check = WSEDEM[WSE_check_NN$nn.idx,3]-data[DEM_check_NN$nn.idx,3]
Depth_check[WSE_check_NN$nn.dists > .05] =0
Depth_check[Depth_check < 0] = 0

length(Depth_check)
dim(data)

#plot(Depth_check, dataDepth$Depth, pch=19, cex=.1,
#main=paste("Site:", site.list$SiteID[k],"
# Modeled Depths vs WSEDEM Depths"))

#Depth_Error = (Depth_check - dataDepth$Depth)
#color = colorp[round((Depth_Error-(-1*max(abs(Depth_Error)))) / (2*max(abs(Depth_Error)))*48)+1]

### Sould make this a default output plot
#plot(dataDepth$X, dataDepth$Y, col = color, pch=19, 
#main="Depth Error: Measured Depth (WSEDEM) - Modeled Depth (m)",
#ylab = "Northing (m)", xlab="Easting (m)")
#legend.text =as.character(round(seq(-1*max(abs(Depth_Error)),max(abs(Depth_Error)),by=(2*max(abs(Depth_Error)))/10),2))

#legend.col = colorp[seq(1,51, by=5)]
#legend("topright",legend.text, pch=19, col=legend.col, 
#ncol=2, bg="white", cex=.7, title="Depth Error (m)")


#WSE_check2_NN = nn2(WSEDEM[,1:2], dataWaterLevel[,1:2],1)

#WSE_Check = WSEDEM[WSE_check2_NN$nn.idx,3]
#plot(WSE_Check, dataWaterLevel[,3],
#main=paste("Site ",site.list$SiteID[k],": 
#Modeled WSE vs WSEDEM"))


#####################################################################
# get initial bed level from DEM file

# find nearest neighbors between original DEM and computational GRID
data.temp = data[data$grid_code > 0,]
NN = nn2(data.temp[,1:2],dataVel[,1:2],1)

# Assign computed depth and WSE from nearest neighbor search index
dataDepth$BedElevation = data.temp$grid_code[NN$nn.idx]
rm(data.temp)

dataDepth$WSE = dataDepth$Depth+dataDepth$BedElevation
dataDepth$WSE[dataDepth$Depth == 0] = dataDepth$BedElevation[dataDepth$Depth == 0]

# dataDepth = dataDepth[dataDepth$WSE > -9000,]
min(dataDepth$WSE)


#################################################################
# Do some plotting to check stuff out, if desired - for debugging.
# Commented out, since it seems to be working
######
#######
## WSE and BedElev #

maxWSE = max(abs(dataDepth$WSE), na.rm=T)
minWSE = min(abs(dataDepth$WSE), na.rm=T)
dataDepth$WSE[dataDepth$WSE < minWSE] = 0
color = colorp[(round((dataDepth$WSE-minWSE)/(maxWSE-minWSE)*48))+1]
color[is.na(dataDepth$WSE)]="white"
color[dataDepth$WSE == 0]="white"

#color[test$Dep<0.001] = "dark blue"

#plot(dataDepth$X, dataDepth$Y, main="Water Surface Elevation",
#   col=color, pch=19, cex=.1)
#
#legend.text =as.character(round(seq(minWSE, maxWSE,by=(maxWSE-minWSE)/10),2))
#legend.col = colorp[seq(1,51, by=5)]
#legend("topright",legend.text, pch=19, col=legend.col)

########################################################



########### Depth Averaged Velocity ####################################
vmag = sqrt(dataVel$XVel^2 + dataVel$YVel^2)
names(dataVel)

maxV = max(vmag, na.rm=T)


color = colorp[round(vmag/maxV*48)+1]
color[vmag==0]="white"

#plot(dataDepth$X, dataDepth$Y, main="Depth Averaged Velocity (m/s)",
#   col=color, pch=19, cex=.1)

#legend.text =as.character(round(seq(0, maxV,by=maxV/10),2))
#legend.col = colorp[seq(1,51, by=5)]
#legend("topright",legend.text, pch=19, col=legend.col)



########### Depth ####################################
maxd = max(abs(dataDepth$Depth))
color = colorp[round(abs(dataDepth$Depth)/maxd*48)+1]
color[is.na(color)]="white"
color[dataDepth$Depth==0] = "white"

#plot(dataDepth$X, dataDepth$Y, main="Depth (m)",
#   col=color, pch=19, cex=.1)


#legend.text =as.character(round(seq(0, maxd,by=maxd/10),2))
#legend.col = colorp[seq(1,51, by=5)]
#legend("topright",legend.text, pch=19, col=legend.col)


#########
# End of optional plots
########################################################
############################################################


# Make dataframe with results
results= data.frame(
"X"=dataDepth$X, "Y"=dataDepth$Y,
"X.Velocity"= dataVel$XVel, "Y.Velocity"=dataVel$YVel, 
"Velocity.Magnitude" = vmag,
"Depth"= dataDepth$Depth, 
"WSE"=dataDepth$WSE, 
"BedLevel" = dataDepth$BedElevation)
#"BedShear_X" = dataBedShear$Bed_Shear_X, "BedShear_Y" =  dataBedShear$Bed_Shear_Y,
#"Vorticity" = dataVorticity$Vorticity )

results[is.na(results)==T] = -9999

# to some cleanup to free up memory
rm(dataVel)
rm(vmag)
rm(dataDepth)
#rm(dataBedShear)
#rm(dataVorticity)
rm(Depth_check)
rm(dataWaterLevel)
gc()



# Write results file.  This file contains one row for
# each points on the computational Grid.  Can be a BIG file.


#dir.create(paste(savedwd,"/",results.folder,"/",site.list$SiteID[k],"_Delft3D_Results.csv", sep=""),
#showWarnings=F)

# Not going to write the computational grid results for now.
#write.csv(results, paste(results.folder,"/",site.list$SiteID[k],"_",site.list$Year[k],"_Delft3D_Results.csv", sep=""), row.names=F)


#############################################################
# Need to translate results back to DEM grid here ####
# And write new output file on DEM grid

DEM.X.Velocity = rep(0, nrow(data))
DEM.Y.Velocity = rep(0, nrow(data))
DEM.Velocity.Magnitude = rep(0, nrow(data))    
DEM.Depth = rep(0, nrow(data))
DEM.WSE = rep(0, nrow(data))
DEM.BedLevel = rep(0, nrow(data))
#DEM.BedShear_X = rep(0, nrow(data))
#DEM.BedShear_Y = rep(0, nrow(data))
#DEM.Vorticity = rep(0, nrow(data))

# match DEM grid points to computational grid points using nearest neighbor search.
# finding 8 nearest neighbors, and distance to each.
nn.DEM = nn2(results[,1:2], data[,1:2],8)


use = results$Depth >-9999


r=raster(xmn= min(round(results$X,2)-.05), xmx=max(round(results$X,2)+.05), 
             ymn= min(round(results$Y,2)-.05),ymx=max(round(results$Y,2)+.05),
             resolution = 1*Meta.Data$Comp.Grid.Spacing)


# Need to break this into chunks to avoid memory nightmares
chunks = seq(1,nrow(data), by = round(nrow(data)/100))
chunks[length(chunks)] = nrow(data)+1


a=rasterize(data.frame(round(results$X,2), round(results$Y,2)), r, field=results$Depth)

print("Depth")
for (c in 1:(length(chunks)-1)){
print(c)
t.DEM.Depth=extract(a, data.frame(round(data[chunks[c]:(chunks[c+1]-1),]$X,2), round(data[chunks[c]:(chunks[c+1]-1),]$Y,2)), method='bilinear')
if (c==1) {DEM.Depth = t.DEM.Depth} else{DEM.Depth[chunks[c]:(chunks[c+1]-1)] = t.DEM.Depth}
}

rm(t.DEM.Depth)
gc()

idx =  (is.na(DEM.Depth)==F)
plot(a, main=paste("depth", site.list$SiteID[k], site.list$VisitID[k]))


print("X.Velocity")
a=rasterize(data.frame(round(results$X,3), round(results$Y,2)), r, field=results$X.Velocity)

for (c in 1:(length(chunks)-1)){
print(c)
t.DEM.X.Velocity=extract(a, data.frame(round(data[chunks[c]:(chunks[c+1]-1),]$X,2), round(data[chunks[c]:(chunks[c+1]-1),]$Y,2)), method='bilinear')
if (c==1) {DEM.X.Velocity = t.DEM.X.Velocity} else{DEM.X.Velocity[chunks[c]:(chunks[c+1]-1)] = t.DEM.X.Velocity}
}

rm(t.DEM.X.Velocity)
gc()

print("Y.Velocity")
a=rasterize(data.frame(round(results$X,2), round(results$Y,2)), r, field=results$Y.Velocity)
for (c in 1:(length(chunks)-1)){
print(c)
t.DEM.Y.Velocity=extract(a, data.frame(round(data[chunks[c]:(chunks[c+1]-1),]$X,2), round(data[chunks[c]:(chunks[c+1]-1),]$Y,2)), method='bilinear')
if (c==1) {DEM.Y.Velocity = t.DEM.Y.Velocity} else{DEM.Y.Velocity[chunks[c]:(chunks[c+1]-1)] = t.DEM.Y.Velocity}
}

rm(t.DEM.Y.Velocity)
gc()


print("WSE")
a=rasterize(data.frame(round(results$X,2), round(results$Y,2)), r, field=results$WSE)
for (c in 1:(length(chunks)-1)){
print(c)
t.DEM.WSE=extract(a, data.frame(round(data[chunks[c]:(chunks[c+1]-1),]$X,2), round(data[chunks[c]:(chunks[c+1]-1),]$Y,2)), method='bilinear')
if (c==1) {DEM.WSE = t.DEM.WSE} else{DEM.WSE[chunks[c]:(chunks[c+1]-1)] = t.DEM.WSE}
}


rm(t.DEM.WSE)
gc()
plot(a, main=paste("WSE", site.list$SiteID[k], site.list$VisitID[k]))

print("BedLevel")
a=rasterize(data.frame(round(results$X,2), round(results$Y,2)), r, field=results$BedLevel)

for (c in 1:(length(chunks)-1)){
print(c)
t.DEM.BedLevel=extract(a, data.frame(round(data[chunks[c]:(chunks[c+1]-1),]$X,2), round(data[chunks[c]:(chunks[c+1]-1),]$Y,2)), method='bilinear')
if (c==1) {DEM.BedLevel = t.DEM.BedLevel} else{DEM.BedLevel[chunks[c]:(chunks[c+1]-1)] = t.DEM.BedLevel}
}
plot(a, main=paste("BedLevel", site.list$SiteID[k], site.list$VisitID[k]))
rm(t.DEM.BedLevel)
gc()




delta = round(max(abs(results$X[2]-results$X[1]), abs(results$Y[2]-results$Y[1])), 1)+.01

# for any point that is too far from a computational grid point on the DEM, 
# assign value of -9999, so we recognize that as "not calculated"
DEM.X.Velocity[nn.DEM$nn.dists[,1] > delta | is.na(DEM.X.Velocity)] = -9999
DEM.Y.Velocity[nn.DEM$nn.dists[,1] > delta| is.na(DEM.Y.Velocity)] = -9999
DEM.Velocity.Magnitude[nn.DEM$nn.dists[,1] > delta| is.na(DEM.Velocity.Magnitude)] = -9999  
DEM.Depth[nn.DEM$nn.dists[,1] > delta| is.na(DEM.Depth)] = -9999
DEM.BedLevel[nn.DEM$nn.dists[,1] > delta| is.na(DEM.BedLevel)] = -9999

## BedLevel should be taken directly from DEM
#DEM.BedLevel= data$grid_code

#get WSE in same process as get depth, to avoid weird grid alignmment artifacts
#DEM.WSE = DEM.BedLevel + DEM.Depth
DEM.WSE[nn.DEM$nn.dists[,1] > delta| is.na(DEM.X.Velocity)] = -9999


# DEM.WSE[DEM.Depth == 0] = -9999

#DEM.BedShear_X[nn.DEM$nn.dists[,1] > delta| is.na(DEM.BedShear_X)] = -9999
#DEM.BedShear_Y[nn.DEM$nn.dists[,1] > delta| is.na(DEM.BedShear_Y)] = -9999
#DEM.Vorticity[nn.DEM$nn.dists[,1] > delta| is.na(DEM.Vorticity)] = -9999

WSE.idx = nn2(WSEDEM[,1:2], data[,1:2],1)$nn.idx
DEM.WSEDEM = WSEDEM[,3][WSE.idx]
min(DEM.WSEDEM[DEM.WSEDEM > -9999])

DEM.Depth.Error = DEM.WSE - DEM.WSEDEM
DEM.Depth.Error[DEM.WSEDEM ==0]= -9999

DEM.Velocity.Magnitude = sqrt(DEM.X.Velocity^2 + DEM.Y.Velocity^2)

# could add these one at a time, above, to avoid too much garbage in memory.
DEM.Results = data.frame(
"X"= data$X, "Y"=data$Y,
"X.Velocity"=DEM.X.Velocity, "Y.Velocity"=DEM.Y.Velocity,
"Velocity.Magnitude" = DEM.Velocity.Magnitude,
"Depth"= DEM.Depth
, 
"WSE"=DEM.WSE, 
"BedLevel" = DEM.BedLevel,
#"BedShear_X" = DEM.BedShear_X, "BedShear_Y" =  DEM.BedShear_Y,
#"Vorticity" = DEM.Vorticity,
"Depth.Error"= DEM.Depth.Error
 )


rm(data)
rm(DEM.X.Velocity)
rm(DEM.Y.Velocity)
rm(DEM.Velocity.Magnitude)
rm(DEM.Depth)
#rm(DEM.BedShear_X)
#rm(DEM.BedShear_Y)
#rm(DEM.Vorticity)
rm(DEM.Depth.Error)
gc()

#DEM.Results = DEM.Results[(DEM.Results$Depth > 0) ,]

DEM.Results = DEM.Results[DEM.Results$Depth > -9999,]
DEM.Results = DEM.Results[((DEM.Results$Depth > 0) & (abs(DEM.Results$Depth.Error) !=0)),]
DEM.Results = DEM.Results[DEM.Results$BedLevel > -9999,]

# Write DEM Grid Output File!
sub.folder = site.list$Results.Folder[k]
sub.folder


dir.create(sub.folder)

write.csv(DEM.Results, paste(sub.folder,"dem_grid_results.csv",sep=""), row.names=F)


#write.csv(Meta.Data, paste(sub.folder,"Meta_Data_",
#   site.list$SiteID[k],"_",site.list$Year[k],"VisitID_",site.list$VisitID[k],"_Delft3D_Results.csv", sep=""))
write.csv(Meta.Data, paste(sub.folder,"Meta_Data.csv", sep=""))

names(Meta.Data)

# Write the xml output file here!

#print(paste("k=",k,sub.folder))
cat("", file= paste(sub.folder,"summary.xml",sep=""))
Meta.Data

names(site.list)
surv_dish = site.list$Measured.Discharge[k] == site.list$Modeled.Discharge[k]
surv_dish

cat(paste("<?xml version=\"1.0\" encoding=\"utf-8\"?>
<model_results>
       <model>delft3D</model>
       <run_datetime>",Meta.Data$Post.Processing.Date.Time,"</run_datetime>
       <version>",Meta.Data$Build.Input.File.R.Version,"</version>
       <delft3d_version>",Meta.Data$Delft3D.Version,"</delft3d_version>
       <rbt_version>",Meta.Data$rbt.version,"</rbt_version>
       <preprocessing_version>",Meta.Data$Build.Input.File.R.Version,"</preprocessing_version>
       <postprocessing_version>",Meta.Data$Build.Input.File.R.Version,"</postprocessing_version>
       <visit>",Meta.Data$VisitID,"</visit>
       <operator>",Meta.Data$Operator,"</operator>
       <year>",Meta.Data$Year,"</year>
       <watershed>",Meta.Data$WatershedName,"</watershed>
       <run_type>","NA","</run_type>
       <surveyed_discharge>", surv_dish, "</surveyed_discharge>
       <qa_status>pass</qa_status>
       <data>
            <d84>",Meta.Data$D84,"</d84>
            <surface_roughness>",Meta.Data$Roughness.Input,"</surface_roughness>
            <measured_discharge>",Meta.Data$Measured_Discharge,"</measured_discharge>
            <modeled_discharge>",Meta.Data$Discharge,"</modeled_discharge>
            <deltabc>",Meta.Data$DeltaBC,"</deltabc>
            <downstream_boundary>",Meta.Data$Exit_BC,"</downstream_boundary>
            <trim_length>",Meta.Data$TrimLength,"</trim_length>
            <inlet>",Meta.Data$Inlet,"</inlet>
            <outlet>",Meta.Data$Outlet,"</outlet>
            <hev>",Meta.Data$HEV,"</hev>
            <left_reference>",Meta.Data$Left.Reference,"</left_reference>
            <top_reference>",Meta.Data$Top.Reference,"</top_reference>
            <computational_grid_spacing>",Meta.Data$Comp.Grid.Spacing,"</computational_grid_spacing>
       </data> 
</model_results>
",sep=""),file= paste(sub.folder,"summary.xml",sep=""),append=T)
##############################################################


# Also, copy the "diag" file into the results folder.
file.copy(from = paste(WorkingDir,"/tri-diag.test", sep=""),
           to= paste(sub.folder,"tri-diag.test", sep=""), overwrite=T)



################################################
###################################
# Optional:  Plot Results for depth and velocity (and whatever else we want)
# Useful for debugging and checkign results
###################################

names(DEM.Results)
colorp

#dir.create(paste(sub.folder,"plots/", sep=""))
j=12
j=7
j=3
for (j in 3:9) {
print(j)

jpeg(
paste(sub.folder, 
#"plots/", 
#site.list$SiteID[k],"_",site.list$Year[k],
#   "VisitID_",site.list$VisitID[k],"_",
names(DEM.Results)[j],".jpg",sep=""),
    7.5,6.5, units='in', res=600)


title = names(DEM.Results)[j]
names(DEM.Results )

max = max(DEM.Results[,j])
min = min(DEM.Results[,j][DEM.Results[,j] > - 9999])
min
max

# standardize for depth error from -1 to 1
if (title == "Depth.Error"){
 title = "Depth Error (m) "
 min = -.5
 max = .5}

if (title == "Velocity.Magnitude") {title = "Velocity (m/s)"}
if (title =="Depth") {title = "Depth (m)"}
if (title == "WSE"){title = "Surface Elevation (m)"}

color.idx = round((DEM.Results[,j]-min)/(max-min) * 48)+1 
color.idx[color.idx < 1] = 1
color  = colorp[color.idx]

par(mfrow = c(1,2))
layout(matrix(c(1,2),1,2), widths=c(7.8,3.2)) 
par(mar=c(5,5,5,0))


plot(DEM.Results$X, DEM.Results$Y, col = color, pch=20, cex= 0.2, 
main=paste(site.list$SiteID[k], "_VisitID ",site.list$VisitID[k],": ", title, sep=""),
xlab="Northing (m)", ylab= "Easting (m)")

par(mar=c(0,1,5,1))
plot.new()


legend.text =as.character(round(seq(min, max, by=(max-min)/9),2))
legend.col = colorp[seq(1,51, by=5)]
legend("topright",legend.text, pch=19, col=legend.col, title=title, cex=.9, ncol=2)


dev.off()
}


rm(DEM.Results)
gc()

################################################################
# Generate project.rs.xml file starting here!!!!
################################################################

all.data = read.csv(paste(scriptdir,"CHaMP_and_AEM_Metrics.csv", sep=""),header=T)
idx=match(Meta.Data$VisitID,all.data$VisitID)
SiteID = all.data$SiteName[idx]

DEMguid = UUIDgenerate(F)
Realizguid = UUIDgenerate(F)
Metadata_guid = UUIDgenerate(F)
DEMGR_guid = UUIDgenerate(F)

Discharge_Meas = Meta.Data$Measured_Discharge==Meta.Data$Discharge
dateCreated=gsub(" ","T",Meta.Data$Post.Processing.Date.Time)


Flow = substr(sub.folder, nchar(sub.folder)-10, nchar(sub.folder)-1)

#########################################################
# Start wriging project.rs.xml file here
#########################################################

cat("", file=paste(sub.folder,"project.rs.xml",sep=""))

cat(paste("","<?xml version=\"1.0\" encoding=\"utf-8\"?>
<Project xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\"
  xsi:noNamespaceSchemaLocation=\"XSD/V1/Project.xsd\">

  <Name>My Hydraulic Output</Name>
  <!--What's in a name. Do we need to decide on this?-->
  <ProjectType>Hydro</ProjectType>

  <MetaData>
    <!--This first metadata exists only to place this project in the riverscapes project-->
    <Meta name=\"Year\">",Meta.Data$Year,"</Meta>
    <Meta name=\"Watershed\">",Meta.Data$WatershedName,"</Meta>
    <Meta name=\"Visit\">",Meta.Data$VisitID,"</Meta>
    <Meta name=\"Site\">",SiteID,"</Meta>
    <Meta name=\"Flow\">",Flow,"</Meta>
  </MetaData>

  <!-- Inputs must have unique Ids that can be referenced in realizations I'm keeping
    legible here so that it's readable in the realiations the GUID is to keep things unique-->
  <Inputs>
    <CSV id=\"DEM1\" guid=\"",DEMguid,"\">
      <Name>DEM</Name>
      <Path>DEM.csv</Path>
      <Project>../HydroModelInputs/artifacts</Project>
    </CSV>
    <CSV id=\"WSEDEM1\" guid=\"",DEMguid,"\">
      <Name>WSEDEM</Name>
      <Path>WSEDEM.csv</Path>
      <Project>../HydroModelInputs/artifacts</Project>
    </CSV>
    <CSV id=\"THALWEG1\" guid=\"",DEMguid,"\">
      <Name>Thalweg</Name>
      <Path>Thalweg.csv</Path>
      <Project>../HydroModelInputs/artifacts</Project>
    </CSV>
  </Inputs>

  <Realizations>
    <Hydro guid=\"",Realizguid,"\" id=\"hydro1\" dateCreated=\"",dateCreated,"\" productVersion=\"",Meta.Data$Post.Processing.R.Code.Version,"\">
      <Name>Hydraulic Model Run1</Name>

      <MetaData>
        <!--This next metadata relates to this particular realization-->
        <Meta name=\"model\">delft3D</Meta>
        <Meta name=\"run_datetime\">",Meta.Data$Post.Processing.Date.Time,"</Meta>
        <Meta name=\"delft3d_version\">",Meta.Data$Delft3D.Version,"</Meta>
        <Meta name=\"rbt_version\">",Meta.Data$rbt.version,"</Meta>
        <Meta name=\"preprocessing_version\">",Meta.Data$Build.Input.File.R.Version,"</Meta>
        <Meta name=\"postprocessing_version\">",Meta.Data$Post.Processing.R.Code.Version,"</Meta>
        <Meta name=\"operator\">Matt Nahorniak</Meta>
        <Meta name=\"run_type\">NA</Meta>
        <Meta name=\"qa_status\">pass</Meta>
      </MetaData>
      <Parameters>
        <Param name=\"surveyed_discharge\">",as.character(Discharge_Meas),"</Param>
        <!--We need to discuss the other parameters that should go here-->
      </Parameters>

      <Inputs>
        <CSV ref=\"DEM1\"/>
        <CSV ref=\"WSEDEM1\"/>
        <CSV ref=\"THALWEG1\"/>
      </Inputs>

      <Analyses>
        <Analysis>
          <Name>Hydraulic Analysis</Name>

          <Metrics>
            <Metric name=\"d84\">",Meta.Data$D84,"</Metric>
            <Metric name=\"surface_roughness\">",Meta.Data$Roughness.Input,"</Metric>
            <Metric name=\"measured_discharge\">",Meta.Data$Measured_Discharge,"</Metric>
            <Metric name=\"modeled_discharge\">",Meta.Data$Discharge,"</Metric>
            <Metric name=\"deltabc\">",Meta.Data$DeltaBC,"</Metric>
            <Metric name=\"downstream_boundary\">",Meta.Data$Exit_BC,"</Metric>
            <Metric name=\"trim_length\">",Meta.Data$TrimLength,"</Metric>
            <Metric name=\"inlet\">",Meta.Data$Inlet,"</Metric>
            <Metric name=\"outlet\">",Meta.Data$Outlet,"</Metric>
            <Metric name=\"hev\">",Meta.Data$HEV,"</Metric>
            <Metric name=\"left_reference\">",Meta.Data$Left.Reference,"</Metric>
            <Metric name=\"top_reference\">",Meta.Data$Top.Reference,"</Metric>
            <Metric name=\"computational_grid_spacing\">",Meta.Data$Comp.Grid.Spacing,"</Metric>
          </Metrics>

          <Outputs>
            <CSV guid=\"",Metadata_guid,"\" id=\"metadata1\">
              <Name>MetaData</Name>
              <Path>Meta_Data.csv</Path>
            </CSV>
            <CSV guid=\"",DEMGR_guid,"\" id=\"demgridresults1\">
              <Name>DEM Grid Results</Name>
              <Path>dem_grid_results.csv</Path>
            </CSV>
          </Outputs>
        </Analysis>
      </Analyses>
    </Hydro>
  </Realizations>

</Project>
",sep=""), file=paste(sub.folder,"project.rs.xml",sep=""), append=T)





} # end k


setwd(savedwd)



