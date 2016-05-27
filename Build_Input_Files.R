##########################################################
# R-code to build hydraulic model input files.
# Author: Matt Nahorniak, South Fork Research
# matt@southforkresearch.org
# Updated 6_23_2015
##########################################################

Build.Input.File.R.Version = 1.0
Delft3D.Version = "4.01.00.rc.02"
Operator = "Matt Nahorniak"


setwd( "C:/Matt-SFR Files/Hydraulic Modeling/R Code to Build Input Files")

#results.folder = "c:/Matt-SFR Files/Hydraulic Modeling/Delft3D Results Files/"
#QA.folder = "c:/Matt-SFR Files/Hydraulic Modeling/Delft3D QA Files/"


#dir.create(QA.folder)


savedwd = getwd()

# R-Code to Build Input Files for CHaMP Hydraulic Modeling
# General:
# 1. Read "Hydraulic Modeling Site List .csv" file, which contains names of
# Sites to Model, along with D84 and Discharge per site
# 2. Read DEM.csv files, WSEDEM.csv, and thalweg.csv files
# 3. Determine inlet and outlet edges
# 4. Create Computational Grid and tranlsate DEM info onto comp. Grid
# 5. Build Input Files for Delft 3D
# 6. Write Meta-data file

# Define "Working Dir" - directory with folders for each site DEM, WSEDEM and Thalweg .csv files
#WorkingDir = "c:/Matt-SFR Files/Hydraulic Modeling Backup/Hydraulic Modeling/Upper_Grande_Ronde_2012_Annual/"
#WorkingDir ="c:/Matt-SFR Files/Hydraulic Modeling Backup/Hydraulic Modeling/Sites in Progress/"
#WorkingDir = "c:/Matt-SFR Files/Hydraulic Modeling Backup/Hydraulic Modeling/Other_Grande_Ronde_Sites/"
#WorkingDir = "c:/Matt-SFR Files/Hydraulic Modeling Backup/Hydraulic Modeling/Restored or Changed Sites/"

## Set program up to record all plots
#windows(record=T)

library(RANN) # library for nearest neighbor search routine.  Very useful!!!
library(XML) # read input xml file to get rbt_version

# Set Some Hard Coded Numerical Parameters for CFD Solution.  Others are read from .csv file.
  # Hard Coded Time Step used in Delft 3D (to be used later)
  dt = .0025 # time step, in minutes

  # Calculate Grid Spacing: To get ~ NX.NY.Ave number of grid points in each direction, on average
  # Could increase NX.NY.Ave w/ more computer power!
  NX.NY.Ave = 700

  # Horizontal Eddy Viscocity in file HEV
  HEV = .01 # horizontal eddy viscosity

# Set the Working Direction (where the "CFD_Site_List.csv" file is located), abd
# Read the file.  This contains site names, D84, discharges, and some info related to 
# Creating a curvilinear grid (which is done in a separate, post-processing script).
site.list=read.csv("CFD_Site_List.csv")

names(site.list)



######################################
# Create directories where Delft3D input files are to be stored and
# where results are to be stored.
# This is needed since we may run multiple discharge variants from the same
# DEM, WSEDEM, and thalweg data files

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
Variant
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
  paste("c://Matt-SFR Files/Hydraulic Modeling/champ data from bucket/",
         site.list$Year,"/",site.list$WatershedName,"/",
         site.list$SiteID,"/VISIT_",
         site.list$VisitID,"/", "Hydro/Results/",Variant, sep="")

Variant
site.list$Results.Folder

for (k in 1:nrow(site.list)){

# Should I have spaces or not?

  dir.create (paste("c://Matt-SFR Files/Hydraulic Modeling/champ data from bucket/",
         site.list$Year, sep="")[k])

  dir.create (paste("c://Matt-SFR Files/Hydraulic Modeling/champ data from bucket/",
         site.list$Year,"/",site.list$WatershedName, sep="")[k])

  dir.create (paste("c://Matt-SFR Files/Hydraulic Modeling/champ data from bucket/",
         site.list$Year,"/",site.list$WatershedName,
         "/",site.list$SiteID,sep="")[k])

  dir.create (paste("c://Matt-SFR Files/Hydraulic Modeling/champ data from bucket/",
         site.list$Year,"/",site.list$WatershedName,
         "/",site.list$SiteID,"/VISIT_",
         site.list$VisitID, sep="")[k])

  dir.create (paste("c://Matt-SFR Files/Hydraulic Modeling/champ data from bucket/",
         site.list$Year,"/",site.list$WatershedName,
         "/",site.list$SiteID,"/VISIT_",
         site.list$VisitID,"/", "Hydro/", sep="")[k])

  dir.create( paste("c://Matt-SFR Files/Hydraulic Modeling/champ data from bucket/",
         site.list$Year,"/",site.list$WatershedName,"/",site.list$SiteID,
         "/VISIT_",
         site.list$VisitID,"/", "Hydro/Results", sep="")[k])

  dir.create( paste("c://Matt-SFR Files/Hydraulic Modeling/champ data from bucket/",
         site.list$Year,"/",site.list$WatershedName,"/",site.list$SiteID,
         "/VISIT_",
         site.list$VisitID,"/", "Hydro/Results/",Variant, sep="")[k])


}

site.list$Results.Folder
#############

}

names(site.list)

# Create an index to match the directory structure, as read by R, and the site list from 
# the "CFD_Site_List.csv" file, so we can match directories and sites correctly.  Note the
# elements in site.list$SiteID must match directory names exactly.
#site.index = match(dir(WorkingDir),site.list$SiteID)

# Create a site index for the sites that are to be modeled.
site.index = 1:nrow(site.list)
site.index = (1:nrow(site.list))[site.list$Model == "Yes"]
site.index

k=site.index[1]
k

######################################################################################################
# Write matlat macro file.  This will be used after running the Delft3D code to convert the
# Delft 3D output to text files which can be read by the post-processing R scripts.  The quickplot
# program that is called is part of the Delft 3D package.

# Start a new macro file....
cat("", file = "quickplot_macro.m")

# Looping Through All Sites for which we need to build input files, and add to macro for later post-processing use

# Note that this will build input files for ALL sites in the specified directory.  To run just
# a subset of sites, comment out the "for (k in site.index) {" line and re-rewite it to inly
# run the index list you want to run.  
names(site.list)
k  = site.index[1]
for (k in site.index) {
#for (k in c(87)){
WorkingDir = as.character(site.list$D3D.Input.Folder[k])


# Set the working directory
#wd = paste(WorkingDir,"/", as.character(site.list$SiteID[k]),"/",sep="")
wd = paste(WorkingDir,"/", sep="")


# Write the macro steps for this site.  This macro will generate output files for:
# depth averaged velocity, water depth, bed shear strewss, water level, vorticity.  
# Each output file also contains X and Y locations by default

cat("d3d_qp('openfile','",wd,"trim-test.dat') 
d3d_qp('selectfield','depth averaged velocity') 
d3d_qp('component','vector') 
d3d_qp('exporttype','sample file') 
d3d_qp('exportdata','",wd,"depth averaged velocity.xyz')
d3d_qp('selectfield','water depth')
d3d_qp('exporttype','sample file') 
d3d_qp('exportdata','",wd,"water depth.xyz') 
d3d_qp('selectfield','bed shear stress')
d3d_qp('exporttype','sample file') 
d3d_qp('exportdata','",wd,"bed shear stress.xyz')
d3d_qp('selectfield','water level')
d3d_qp('exporttype','sample file') 
d3d_qp('exportdata','",wd,"water level.xyz') 
d3d_qp('selectfield','vorticity')
d3d_qp('exporttype','sample file') 
d3d_qp('exportdata','",wd,"vorticity.xyz') 
",
sep="",
#file = "c://Matt-SFR Files/Hydraulic Modeling Backup/Hydraulic Modeling/quickplot_macro.m",
file = "quickplot_macro.m",
append=T)
}


#### Done looping through sites to build post-processing macro "quickplot_marco.m"
###########################################################################################

###########################################################################################
# Loop through sites to build batch file to run all simulations in batch mode
pathname = "path C:\\Matt-SFR Files\\Delft3D\\Delft3D_Updated\\delft3d_ohmw_4.01.00.rc.02\\delft3d\\win32\\flow2d3d\\bin\n"

cat(pathname, file = "batchprocess.bat")
k
k=1

for (k in site.index) {
#for (k in c(87)){
WorkingDir = as.character(site.list$D3D.Input.Folder[k])
WorkingDir
# Set the working directory
#wd = paste(WorkingDir,"/", as.character(site.list$SiteID[k]),sep="")
wd = WorkingDir
wd=normalizePath(wd)
wd
cat("chdir \"",wd,"\"", sep="", file = "batchprocess.bat",append=T)
cat("
d_hydro test.xml
", file = "batchprocess.bat",append=T)
}
# set batchfile to automatically shutdown to avoid leaving
# on AWS server instance when job finished.
cat("cd c:\\Users\\Administrator
", file = "batchprocess.bat", append=T)
cat("shutdown -f", file = "batchprocess.bat", append=T)





######################################################################################
# Start Looping Through All Sites for Which we need to build input files!!

# Note that this will build input files for ALL sites in the specified directory.  To run just
# a subset of sites, comment out the "for (k in site.index) {" line and re-rewite it to inly
# run the index

k=site.index[1]
k

for (k in site.index) {
#for (k in c(26:27)) {
print(paste("k=",k, "site=", site.list$SiteID[k] ))

# Set Working Dir

WorkingDir = as.character(site.list$Directory[k])
WorkingDir

# Move to the specific directory for site site.list$SiteID[k]:
WorkingDir
#setwd(paste(WorkingDir,"/", as.character(site.list$SiteID[k]), sep=""))
setwd(WorkingDir)


# Need lots of digits!
options(digits=12)

#########################
## read in the xmlfile "hydro_prep_log.xml" and find the RBT version
 
rbt.version = "unknown"
if ("hydro_prep_log.xml" %in% dir()){

xmlfile = xmlParse("hydro_prep_log.xml")

xmltop = xmlRoot(xmlfile) #gives content of root

rbt.version = xmlChildren(xmltop)[c("rbt_version")]$rbt_version[[1]]
rbt.version = as(rbt.version, "character")
}
rbt.version
#
#
#
##########################


# Read the DEM into data`.frame "data"
data =read.csv("DEM.csv", header=F)

# Re-read if data has headers (i.e. new format from Matt R's C++ code)
if (is.numeric(data[,1])==F){
data= read.csv("DEM.csv", header=T)}


names(data) = c("X", "Y", "grid_code")
ref.left = min(data$X)-.05
ref.top = max(data$Y)+.05

max(data[,1])
max(data[,2])
max(data[,3])


#Remove any "-9999" values from DEM
data= data[data$X != (-9999),]
data= data[data$Y != (-9999),]
data= data[data$grid_code != (-9999),]


if ( min(data$grid_code[data$grid_code > -9999]) < 1){
print(paste(site.list$VisitID[k], site.list$SiteID[k], min(data$grid_code[data$grid_code > -9999])))
}

nrow(data)

names(data)

names(data)
# Read teh WSEDEM into data.frame "WSEDEM"
WSEDEM = read.csv("WSEDEM.csv", header=F)

# Re-read if data has headers (i.e. new format from Matt R's C++ code)
if (is.numeric(WSEDEM[,1])==F)
{WSEDEM= read.csv("WSEDEM.csv", header=T)}


names(WSEDEM) = c("ws.X", "ws.Y", "ws.Z")
#max(WSEDEM$ws.Z)
WSEDEM = WSEDEM[WSEDEM$ws.Z > -9999,]
#min(WSEDEM$ws.Z)
#max(WSEDEM2$ws.Z)
#hist(WSEDEM2$ws.Z)

#WSEDEM2 = WSEDEM[WSEDEM[,3] > -9999,]
#nrow(WSEDEM2)

# Discharge is flow rate in cubic m/s
discharge = site.list$Modeled.Discharge[k]

# slop defines how much of the DEM grid is going to be trimmed off to ensure
# our inlet and outlet boundaries fall along one face (East, North, West, or South) of the
# Computation Grid.  This will be improved at some point to allow boundaries more flexibility

slop = site.list$Trim.Length[k]
print(paste("discharge =", discharge))

# Surface Roughness Parameters.  Not defined differently for X and Y directions (u and v on comp. grid)
#Ccofu  =  site.list$Roughness[k] # surface roughness in u direction (m)
#Ccofv  =  site.list$Roughness[k] # surface roughness in v direction (m)
Ccofu  =  site.list$D84[k]*6/1000 # surface roughness in u direction (m)
Ccofv  =  site.list$D84[k]*6/1000 # surface roughness in v direction (m)



# "grid_code" is the elevation (left over name from when I was manually converting GIS data into
# .csv files, and "grid_code" was the default name that was spit out.
data = data[data$grid_code>0,]

# Get max width and height of DEM 
x.length=max(data$X)-min(data$X)
y.length=max(data$Y)-min(data$Y)

# Calculate density of points in x and y directions
X.p.meter = NX.NY.Ave / sqrt(x.length*y.length)
Y.p.meter = NX.NY.Ave / sqrt(x.length*y.length)

# Use only x-density to define computation grid spacing
DX = 1/X.p.meter

# Force Grid Spacing to be evenly divisible by, or into, .1 meter.  This is done so
# that the location of results match, exactly, the location of the original DEM points.
if(DX <= .025) {DX = .025} else
if(DX <= .05) {DX = .05} else
if(DX <= .1) {DX = .1} else
DX = round(DX, 1)

#} # end k - TEMP ######


#Ensure X and Y are to the nearest mm only.  
data$X=round(data$X,2)
data$Y=round(data$Y,2)




# Compute all unique X and Y values in the grid
Xlist = as.numeric(levels(as.factor(data$X)))
Ylist = as.numeric(levels(as.factor(data$Y)))

maxX = max(Xlist)
minX = min(Xlist)
maxY = max(Ylist)
minY = min(Ylist)

### Plot a rectangle that contains all DEM points (and some empty points)
#plot(c(minX, maxX, maxX, minX, minX),
#    c(minY, minY, maxY, maxY, minY), type="l", col="blue",
#    xlab = "Easting", ylab="Northing")
#
#points(data$X, data$Y, col="gray", pch=19, cex=.1)


#################################################################
# Figure out which borders are closest to thalweg inlet and outlet,
# adjust gridX and gridX max's, min's so the don't thalweg crosses
# border of grid.  Will have to cut into grid just a bit ("slop")

# REad the thalweg
thalweg = read.csv("Thalweg.csv", header=F)

# Re-read if data has headers (i.e. new format from Matt R's C++ code)
if (is.numeric(thalweg[,1])==F){
thalweg= read.csv("Thalweg.csv", header=T)}
names(thalweg)
names(thalweg) = c("thalwagX", "thalwagY")

## plot the thalweg
#lines(thalweg)

# match he thalweg points with the water surface points.  With this we'll
# be able to figure out which way the water is flowing
thalwegZ.idx = nn2(as.data.frame(WSEDEM[,1:2]),thalweg,1)$nn.idx
thalwegZ = WSEDEM[thalwegZ.idx,3]

#points(thalweg)

# We that the thalweg points ordered upstream to down stream.  It's now always that way, so,
# if needed, switch thalweg order if currently ordered downstream to upstream

if (mean(thalwegZ[1:10] ) < mean(thalwegZ[length(thalwegZ)-10]:thalwegZ[length(thalwegZ)]))
 {
thalweg = thalweg[nrow(thalweg):1,]
thalwegZ = thalwegZ[length(thalwegZ):1]
}


### Here!!!!##
# Move to the folder where I'm going to write D3D input files
setwd(site.list$D3D.Input.Folder[k])
##############

# write the updated thalweg to a file for use in post-processing.
write.csv(thalweg, "updated_thalweg.csv")
attach(thalweg)

# Now that we know the upstream and downstream sites, we can define an inlet end and an outlet end,
# and fine the nearest grid point to the inlet and outlet
inletX = thalwagX[1] 
inletY = thalwagY[1]
idx = nn2(data.frame(data$X, data$Y),data.frame(inletX, inletY),1)$nn.idx[1]

inletX = data$X[idx]
inletY = data$Y[idx]

##Inlet
#points(inletX, inletY, pch=19, cex=2, col="red")

# find X distance to left and right boundaries, and Y distance to top and bottom boundaries.
# The minimum of these distances will define which boundary is the inlet side of the flow.
dwest = abs(inletX-minX)
deast = abs(inletX-maxX)
dnorth = abs(inletY-maxY)
dsouth = abs(inletY-minY)


if (site.list$SiteID[k]== "LEM00002-00001B") {dsouth = 0}
if (site.list$SiteID[k]== "CBW05583-029535") {dwest = 0}
if (site.list$SiteID[k]== "CBW05583-515058_ModifiedDEM") {dsouth=0}
if (site.list$SiteID[k]== "ENT201301-TyeeSide") {dwest = 0}
if (site.list$SiteID[k]== "ENT00001-1BC11") {dwest = 0}



# Use the minimum distance to define inlet side. Cut some "slop" off the DEM grid to create the
# actual computational grid, to ensure the computation inlet boundary entirely cross the inlet of
# the stream.

# Easier way... could re-write per below
#ord=rank(c(dwest, deast, dnorth, dsouth))
#dirs = c("west","east", "north", "south")
#inlet1 = dirs[ord==1]
#inlet2 = dirs[ord==2]


if (min(dwest, deast, dnorth, dsouth)==dwest){ 
minX = inletX + slop

inlet = "west"}

if (min(dwest, deast, dnorth, dsouth)==deast){ 
maxX = inletX - slop
inlet = "east"}

if (min(dwest, deast, dnorth, dsouth)==dnorth){ 
maxY = inletY - slop
inlet = "north"}

if (min(dwest, deast, dnorth, dsouth)==dsouth){ 
minY = inletY + slop
inlet="south"}

###########
## repeat all of the above for outlet
############




outletX = thalwagX[length(thalwagX)] 
outletY = thalwagY[length(thalwagY)]

idx = nn2(data.frame(data$X, data$Y),
data.frame(outletX, outletY),1)$nn.idx[1]
outletX = data$X[idx]
outletY = data$Y[idx]


#points(outletX, outletY, pch=19, col="green")

dwest = abs(outletX-minX)
deast = abs(outletX-maxX)
dnorth = abs(outletY-maxY)
dsouth = abs(outletY-minY)

###################################################
# A few sites, so far, have broken the code!  Manually for boundaries for
# The following!

if (site.list$SiteID[k]== "CBW05583-086186") {dwest = 0}
if (site.list$SiteID[k]== "CBW05583-142490") {dnorth = 0}
if (site.list$SiteID[k]== "CBW05583-232818") {dnorth = 0}
if (site.list$SiteID[k]== "CBW05583-312265") {dsouth = 0}
if (site.list$SiteID[k]== "LEM00001-Little0Springs-2") {dnorth = 0}
if (site.list$SiteID[k]== "LEM00002-00001B") {dnorth = 0}
#####################################################
if (min(dwest, deast, dnorth, dsouth)==dwest) {
minX = outletX + slop
outlet = "west"}

if (min(dwest, deast, dnorth, dsouth)==deast) {
maxX = outletX - slop
#minX = minX + slop
#maxX = maxX - slop
outlet = "east"}

if (min(dwest, deast, dnorth, dsouth)==dnorth) {
maxY = outletY - slop
#minY = minY + slop
#maxY = maxY - slop
outlet = "north"}

if (min(dwest, deast, dnorth, dsouth)==dsouth) {
minY = outletY + slop
#minY = minY + slop
#maxY = maxY - slop
outlet = "south" }



#lines(c(minX, maxX, maxX, minX, minX),
#    c(minY, minY, maxY, maxY, minY), type="l", col="red")


#plot(data$X, data$Y, col="gray", pch=19, cex=.1,
# ylim=c(min(data$Y), min(data$Y)+1.1*(max(data$Y)-min(data$Y)))
#)
#lines(c(minX, maxX, maxX, minX, minX),
#    c(minY, minY, maxY, maxY, minY), type="l", col="red")


#legend("topleft", legend = c("Reach Extents", "Comp. Grid Extent"), 
#lt=c(1,1), lw=c(4,1), col=c("gray", "red"))

# will use inlet and outlet later when defining inlet and exit
# boundary conditions and writing those files


# plot new extents for grid
#lines(c(minX, maxX, maxX, minX, minX),
#    c(minY, minY, maxY, maxY, minY), col="blue")


xrange = maxX-minX
yrange = maxY-minY


NX = round((maxX - minX)/DX)
NY = round((maxY - minY)/DX)

# Define X and Y values for the computational Grid
Xvals = seq(minX, maxX, by = DX)
Yvals = seq(minY, maxY, by = DX)
NX = length(Xvals)
NY = length(Yvals)

#Define Computation Grid:  A 2-D array for each of X and Y points
GridX = array(NX*NY, c(NX, NY))
GridY = GridX

# X points of computational Grid
length(GridX[,1])
for (j in 1:NY) {
GridX[,j] = Xvals
}

#Y Points of computational Grid
for (i in 1:NX) {
GridY[i,] = Yvals
}



########################################################
# OK, We've got X and Y computational Grids defined.
# Next Step is to define Bed Elevation (Z) at each grid point
########################################################################


#Efficient way to find nearest neighbors and assign elevation..

## find bathymetry of bottom of channel
# Need to build data.frame of X,Y points in correct formaty for function nn 
# (a data.frame of (X,Y) pairs

GXY = data.frame(as.vector(GridX), as.vector(GridY))
DatXY = data.frame(data$X, data$Y)

# find the nearest x-y points for every point in the DEM
nearest = nn2(DatXY, GXY, 1)

# Assign the Z value from the DEM to GZvec - which corresponds to the X,Y points in GXY
GZvec = data$grid_code[nearest$nn.idx]

# If there's no "close" point to match up with, just assign the maximum elevation (some
# points that end up in the computation grid aren't in the DEM, as the DEM isn't rectangular.
# Just give them a "high" elevation so they won't be involved in the fluidics

GZvec[nearest$nn.dists > .5] = max(data$grid_code)

# Now create a GridZ that's in the right format to mach the X and Y matrix formats
GridZ = matrix(GZvec, c(NX, NY))

########################################################


### Also find WS elevation.  This will be used to define exit boundary condition. 
# This is done similarly to above
WSvec = GZvec
WSDatXY = data.frame(WSEDEM$ws.X, WSEDEM$ws.Y)
DatXY
data.frame(WSEDEM$ws.X, WSEDEM$ws.Y)
nearest = nn2(WSDatXY, GXY, 1)
WSvec  = WSEDEM$ws.Z[nearest$nn.idx]

# Prevent max water being higher than max elevation, since this happened once.
# Not sure if I should error proof this or not.
# (3/23/2016)
WSvec[WSvec > max(GZvec)] = max(GZvec)
# if bathymetry is above water surface, set water surface to bathymetry
WSvec[GZvec >= WSvec] = GZvec[GZvec >= WSvec]


#Convert to Depth by subtracting water surface elevation
WSvec = WSvec - max(GZvec)
max(WSvec)
# Convert to matrix in same format as X, Y, and Z matrices
WS.Z = matrix(WSvec, c(NX, NY))



###########################################################
# Do some graphing to check things out.  
# Plot Bathymetry
colorp = topo.colors(50)
colorp = rainbow(50)
GridZ[is.na(GridZ)] = max(na.omit(GridZ))

depth = max(GridZ) -GridZ
max(depth)

color = colorp[round(depth/max(depth)*48)+1]
color[depth == 0] = "white"



sub.folder=site.list$Results.Folder[k]
sub.folder
#dir.create(gsub("Inputs", "Variants",site.list$Directory[k]))
#sub.folder
#dir.create(sub.folder)

QA.folder = sub.folder


# QA plots used to go in separate folder.  Now sticking them here to be consistent
# with output on cm.org.
#paste(sub.folder,"QA Plots/", sep="")
#QA.folder
#dir.create(QA.folder)
#dir(QA.folder)


#jpeg(paste(QA.folder,"Bathymetry_",site.list$SiteID[k],site.list$Year[k],".jpg",sep=""), 6,6, units='in', res=600)
jpeg(paste(QA.folder,"Bathymetry.jpg",sep=""), 6,6, units='in', res=600)
plot(GridX, GridY, col=color, pch=19, cex=.1,
main=paste(site.list$SiteID[k],"Delft 3D Bathymetry Input"),
xlab="easting",
ylab="northing")

legend.text =as.character(0: round(max(depth)))

legend("topright", legend.text, ncol=2, pch=19, 
col= rainbow(length(legend.text)),
title= "Bed Level (m, below max Z")

dev.off()

##########################


#######################################################
# Plot Water Surface Elevation #########################
#

#jpeg(paste(QA.folder,"WS_Elevation_",site.list$SiteID[k],site.list$Year[k],".jpg",sep=""), 6,6, units='in', res=600)
jpeg(paste(QA.folder,"WS_Elevation.jpg",sep=""), 6,6, units='in', res=600)


color = colorp[round((-1*WS.Z)/max(-1*WS.Z)*48)+1]
color[depth == 0] = "white"
plot(GridX, GridY, col=color, pch=19, cex=.1,
main=paste(site.list$SiteID[k],"Delft 3D Water Surface Elevation Input"),
xlab="easting",
ylab="northing")
#lines(thalweg)

#points(TopX, TopY, pch=19, cex=.01,col="blue")
#points(BottomX, BottomY, pch=19, cex=.01, col="blue")


legend.text =as.character(round(seq(min(abs(WS.Z)), max(abs(WS.Z)), by=(abs(max(WS.Z)-min(WS.Z))/9)),2))


legend("topright", legend.text, ncol=2, pch=19, 
col= rainbow(length(legend.text)),
title= "Depth (m) below max Z")
dev.off()
###################################################
###################################################

# Calulate approximate volume of water and estimated simulation
# time required:  Need to "fill" entire volume with discharge,
# and then allow to run until stable.  Use ~ 2X time to fill volume at discharge rate

Wdepth = WS.Z + depth
Wdepth
VOL = sum(Wdepth)*DX*DX
print("VOL")
print(VOL)

# Adjusted this line to ensure sufficient "fill" when we're running at higher discharges,
# which will require higher volume to fill.  It's a swag but should work.
Min.Sim.sec = VOL/discharge * sqrt(discharge/site.list$Measured.Discharge[k]) # time to fill volume, based on discharge, in seconds
if (site.list$Measured.Discharge[k] == 0) {Min.Sim.sec = VOL/discharge+10*60} # To prevent breakage if measured discharge =0

simtime = max(10, round(Min.Sim.sec / 60 * 2))


######################################################
# Plot Water Depth #########################
#

min(Wdepth)
max(Wdepth)
#jpeg(paste(QA.folder,"W_Depth",site.list$SiteID[k],site.list$Year[k],".jpg",sep=""), 6,6, units='in', res=600)
jpeg(paste(QA.folder,"W_Depth.jpg",sep=""), 6,6, units='in', res=600)

color = colorp[round((Wdepth)/max(Wdepth)*48)+1]
color[Wdepth == 0] = "white"
plot(GridX, GridY, col=color, pch=19, cex=.1,
main=paste(site.list$SiteID[k],"WSEDEM & DEM Depth Input"),
xlab="easting",
ylab="northing")
#lines(thalweg)

#points(TopX, TopY, pch=19, cex=.01,col="blue")
#points(BottomX, BottomY, pch=19, cex=.01, col="blue")


legend.text =as.character(round(seq(min(abs(Wdepth)), max(abs(Wdepth)), by=(abs(max(Wdepth)-min(Wdepth))/9)),2))


legend("topright", legend.text, ncol=2, pch=19, 
col= rainbow(length(legend.text)),
title= "Depth (m)")
dev.off()




#######################################################

# Print approximate simulation end point (simtime) required to achieve stable
# simulation results.  This tends to be a very conservative estimate!
print(paste("Site",site.list$SiteID[k],"required simulation time ~", simtime, "minutes"))



##########################################################################
###########################################################################
# Now we have all the information for the computational grid and bathymetry,
# and we can start generating input files for Delft 3D.
#############################################
####################################################################
#define exit boundary condition and write "test.bnd" and "test.src"
# define thalweg.g as only part of thalweg on grid

thalweg.g = thalweg[(thalweg[,1] <= maxX),]
thalweg.g = thalweg.g[thalweg.g[,1] >= minX,]
thalweg.g = thalweg.g[thalweg.g[,2] <= maxY,]
thalweg.g = thalweg.g[thalweg.g[,2] >= minY,]

#thalweg.g
#lines(thalweg.g, lw=2)
outlet

# Define Exit Boundary
if (outlet=="east") {

# find nearest boundary point to outlet
idx= nn2(cbind(GridX[NX,], GridY[NX,]),thalweg.g[nrow(thalweg.g),],1)$nn.idx

 outflow.y.idx.b = idx- match(T, (depth[NX, idx:1] == 0))+4
 outflow.y.idx.t =  idx+ match(T, (depth[NX, idx:NY] == 0))-4

 if(is.na(outflow.y.idx.b)){outflow.y.idx.b =2}
 if(is.na(outflow.y.idx.t)){outflow.y.idx.t =NY-1}
 outflow.y.idx = c(outflow.y.idx.b, outflow.y.idx.t)

   outflow.x.idx = c(NX, NX)
outflow.y.idx
#   outflow.y.idx = c(match(T,depth[NX,] > 0)+1,(NY-match(T,depth[NX,(NY:1)]>0)))
 }
if (outlet=="west") {

idx= nn2(cbind(GridX[1,], GridY[1,]),thalweg.g[nrow(thalweg.g),],1)$nn.idx
   outflow.x.idx = c(1, 1)
 outflow.y.idx.b = idx- match(T, (depth[1, idx:1] == 0))+4
 outflow.y.idx.t =  idx+ match(T, (depth[1, idx:NY] == 0))-4
 if(is.na(outflow.y.idx.b)){outflow.y.idx.b =2}
 if(is.na(outflow.y.idx.t)){outflow.y.idx.t =NY-1}
 outflow.y.idx = c(outflow.y.idx.b, outflow.y.idx.t)
   outflow.y.idx = c(match(T,depth[1,] > 0)+1,(NY-match(T,depth[1,(NY:1)]>0)))
 
 }

if (outlet=="north") {   outflow.y.idx = c(NY, NY)
idx= nn2(cbind(GridX[,1], GridY[,1]),thalweg.g[nrow(thalweg.g),],1)$nn.idx
 outflow.x.idx.l = idx- match(T, (depth[idx:1,NY] == 0))+4
 outflow.x.idx.r =  idx+ match(T, (depth[idx:NX,NY] == 0))-4
 if(is.na(outflow.x.idx.l)){outflow.x.idx.l =2}
 if(is.na(outflow.x.idx.r)){outflow.x.idx.r =NX-1}
 outflow.x.idx = c(outflow.x.idx.l, outflow.x.idx.r)
 outflow.x.idx


#   outflow.x.idx = c(match(T,depth[,NY] > 0)+1,(NX-match(T,depth[NX:1,NY]>0)))
}

if (outlet=="south") {
idx= nn2(cbind(GridX[,1], GridY[,1]),thalweg.g[nrow(thalweg.g),],1)$nn.idx
 outflow.x.idx.l = idx- match(T, (depth[idx:1,1] == 0))+4
 outflow.x.idx.r =  idx+ match(T, (depth[idx:NX,1] == 0))-4
 if(is.na(outflow.x.idx.l)){outflow.x.idx.l =2}
 if(is.na(outflow.x.idx.r)){outflow.x.idx.r =NX-1}
 outflow.x.idx = c(outflow.x.idx.l, outflow.x.idx.r)
 outflow.x.idx
   outflow.y.idx = c(1, 1)
#   outflow.x.idx = c(match(T,depth[,1] > 0)+1,(NX-match(T,depth[NX:1,1]>0)))
 }

## Draw Line along entire exit boundary
#lines(c(minX, maxX, maxX, minX, minX),
#    c(minY, minY, maxY, maxY, minY), type="l", col="blue")

# Draw line alone exit boundary that's within the defined bathymetry
#lines(GridX[cbind(outflow.x.idx, outflow.y.idx)], 
#      GridY[cbind(outflow.x.idx, outflow.y.idx)], col="black", lw=3)

# fix outflow WS level so it's exactly the mean
# WSE along the exit boundary

# Find output points
outflow.idx = nn2(
  WSEDEM[,1:2],data.frame(cbind(GridX[
      outflow.x.idx[1]:outflow.x.idx[2],outflow.y.idx[1]:outflow.y.idx[2]],
        GridY[outflow.x.idx[1]:outflow.x.idx[2], outflow.y.idx[1]:outflow.y.idx[2]])),1)

outflow.ws.level = mean(WSEDEM[outflow.idx$nn.idx[outflow.idx$nn.dists < .11],3])-max(GridZ)

# This is used ONLY is we add a DeltaBC to the cfd.site.list file, which we might do if
# we're running on non-default flows and we've got some info on the change in depth.
outflow.ws.level = outflow.ws.level + site.list$DeltaBC[k]

dim(WSEDEM)
dim(GridZ)

outflow.ws.level

## 7_23_15 Change
## Not adjusting outflow BC ##
## We found that this drastically over-corrected exit boundary condition
## and, likely, caused more error than it reduced.  Better to just leave this
## BC as is, as it likely communicates upstream to a lesser distance than
## adjusting it too much, especially when modeling at flows significantly higher
## than measured.  Still work to do on this front.
if (1==2) {
##########################
### HERE!!!!
# outlet.width is the wetted width at the exit, to be used in conjunction with Manning's 
# equation and a 1-D simplifying assumption of rectangular open channel flow to come up
# with downstream boundary condition at unmeasured discharge values.

X.outlet=(WSEDEM[outflow.idx$nn.idx[outflow.idx$nn.dists < .11],1])
Y.outlet=(WSEDEM[outflow.idx$nn.idx[outflow.idx$nn.dists < .11],2])
Z.idx = nn2(data[,1:2], data.frame(X.outlet, Y.outlet),1)[1]
outlet.dep = outflow.ws.level-(data[Z.idx$nn.idx,3]-max(GridZ))
length(outlet.dep)

W = length((WSEDEM[outflow.idx$nn.idx[outflow.idx$nn.dists < .11],3])-max(GridZ))*DX
D1 = mean(outlet.dep)
# Need
Q1=site.list$Measured.Discharge[k]
Q2=site.list$Modeled.Discharge[k]
# Solve for D2


####################################
# Iterative solution to Manning's equation to get new depth estimate at modeled discharge,
# Assuming 1-D open channel flow in rectangular basin.  Not a great estimate, but we need
# Something, and we've shown we're fairly insensitive to bad BC at exit.

# Here's the function for manning Q1/Q2, assuming |____| shaped
# channel, where W is constant, D1 is known, Q1 and Q2 are known, and
# we want D2
###
A = function(D1, D2, W) {
Ratio = D1/D2 * ((2*D1+W)/(2*D2+W))^(2/3)
return(Ratio)}
###

# Here's the known ratio of measured (Q1) to modeled (Q2)
Actual.Ratio = Q1/Q2

# Starting Guess (D2.Guess set high so we definitely enter while loop
D2.Guess = 9999
D2.Guess.New = D1


# Iterate until our guess for D2 gives us the right ratio within .001
while ((abs(D2.Guess.New-D2.Guess)) > .001) {
D2.Guess = D2.Guess.New
#print(D2.Guess)
Ratio = A(D1, D2.Guess, W)
#print(Ratio)
Actual.Ratio

D2.Guess
D2.Guess.offset =  D2.Guess * Ratio/Actual.Ratio-D2.Guess
D2.Guess.offset
# The 0.1 forces us to creep up slowly, for stability's sake
D2.Guess.New = D2.Guess + .1*D2.Guess.offset
}
###################################
D2 = D2.Guess.New


########################

###############################
# Alternate: Adjust exit boundary condition by scaling exit depth to discharge change
# Assuming no velocity change.  This should be the maximum allowable change; in theory
# the above method could violate this for very non-rectangular channels.
 DisMult = Q2/Q1
outflow.idx
# Find output points
adjust.idx = nn2(
  data[,1:2],data.frame(cbind(GridX[
      outflow.x.idx[1]:outflow.x.idx[2],outflow.y.idx[1]:outflow.y.idx[2]],
        GridY[outflow.x.idx[1]:outflow.x.idx[2], outflow.y.idx[1]:outflow.y.idx[2]])),1)

depth.adjust=WSEDEM[outflow.idx$nn.idx,3]-data[adjust.idx$nn.idx,3]
val = DisMult* sum(depth.adjust[depth.adjust > 0])
cor = seq(-10,10, by=.01)
cor.idx = 1:(length(cor))
val.try = rep(0, length(cor))
c = cor[1]
for (c in 1:length(cor)){
val.try[c] = abs(sum((depth.adjust+cor[c])[((depth.adjust+cor[c])>0)])- val)
}
correction=cor[cor.idx[val.try == min(val.try)]]
correction
D2.max = D1+correction
D2.max
Q1
Q2
D1
D2
D2.max
if (Q2 > Q1) {D2 = min(D2, D2.max)}
if (Q2 <= Q1) {D2 = max(D2, D2.max)}


###################################################
########################
# Here's where we adjust the exit WS level
outflow.ws.level = outflow.ws.level + (D2-D1)
outflow.ws.level

###########################
}# End of "if 1==2" to avoid doing any of the above BC adjustments.  



### Here's my outflow boundary condition ##################


#write file "test.bnd"
cat("outflow             ","Z","T", outflow.x.idx[1], outflow.y.idx[1],
outflow.x.idx[2], outflow.y.idx[2], "0.00",
file = "test.bnd")

outflow.ws.level

#write file "test.bct"

cat("table-name           'Boundary Section : 1'","\n", file="test.bct")
cat("contents             'Uniform             '","\n", file=  "test.bct",append=T)
cat("location             'outflow             '","\n", file=  "test.bct",append=T)
cat("time-function        'non-equidistant'","\n", file=  "test.bct",append=T)
cat("reference-time       20130101","\n", file=  "test.bct",append=T)
cat("time-unit            'minutes'","\n", file=  "test.bct",append=T)
cat("interpolation        'linear'","\n", file=  "test.bct",append=T)
cat("parameter            'time                '                     unit '[min]'","\n", file=  "test.bct",append=T)
cat("parameter            'water elevation (z)  end A'               unit '[m]'","\n", file=  "test.bct",append=T)
cat("parameter            'water elevation (z)  end B'               unit '[m]'","\n", file=  "test.bct",append=T)
cat("records-in-table     2","\n", file=  "test.bct",append=T)
cat(" 0.0000000e+000 ", outflow.ws.level,"  ",outflow.ws.level,"\n", file=  "test.bct",append=T)
cat(simtime, outflow.ws.level,"  ",outflow.ws.level,"\n", file=  "test.bct",append=T)

##############################################################################



###################################################################
# Inlet boundary condition
# Using discharge instead... distribute discharge over mouth of wetted inlet
# boundary, proportional to depth^2 (arbitrary, but probably somewhat accurate
# Note, this is different than a "boundary condition" as you'd normally specify.
# I've found tht specifying an inlet flow boundary condition tends to crash Delft3D
# When I start with an empty channel.  This seems to work just fine.

if (inlet=="south") {
inlet.x = (1:NX)[Wdepth[,2] > 0]
if (length(inlet.x) == 0) {inlet.y=(1:NX)[depth[2,]==max(depth[2,])]}
inlet.y = rep(2, length(inlet.x))}

if (inlet=="north") {
inlet.x = (1:NX)[Wdepth[,(NY-1)] > 0]
if (length(inlet.x) == 0) {inlet.y=(1:NX)[depth[(NY-1),]==max(depth[(NY-1),])]}
inlet.y = rep((NY-1), length(inlet.x))}

if (inlet=="west") {
inlet.y = (1:NY)[Wdepth[2,] > 0]
if (length(inlet.y) == 0) {inlet.y=(1:NY)[depth[2,]==max(depth[2,])]}
inlet.x = rep(2, length(inlet.y))}

if (inlet=="east") {
inlet.y = (1:NY)[Wdepth[(NX-1),] > 0]
# if all points dry, assign inlet as lowest point 
if (length(inlet.y) == 0) {inlet.y=(1:NY)[depth[(NX-1),]==max(depth[(NX-1),])]}

inlet.x = rep((NX-1), length(inlet.y))

}



inlet.y
inlet.x
inlet.depth = depth[cbind(inlet.x, inlet.y)]

###################################################
# Fix issue
# if inlet and outlet are on the same side, make sure discharge is zero
# at points closer to outlet
if (inlet == outlet) {
   d.to.outlet =  sqrt((GridX[cbind(inlet.x, inlet.y)] - outletX)^2 + (GridY[cbind(inlet.x,inlet.y)] - outletY)^2)
   d.to.inlet =  sqrt((GridX[cbind(inlet.x, inlet.y)] - inletX)^2 + (GridY[cbind(inlet.x,inlet.y)] - inletY)^2)
   inlet.depth[d.to.outlet < d.to.inlet] = 0
 }

inlet.y = inlet.y[inlet.depth > 0]
inlet.x = inlet.x[inlet.depth > 0]
inlet.depth = inlet.depth[inlet.depth > 0]

#######################################################


##Plot points on the latest plot, just to see where the inlet discharge ended up.
#points(
#GridX[cbind(inlet.x, inlet.y)]
#,GridY[cbind(inlet.x,inlet.y)], col="green", pch=19)


#make inlet discharge per cell proportional to depth^2 
# change to make dishcarge proportional to depth (not squared)
inlet.discharge = discharge* inlet.depth / sum(inlet.depth)

# Checks....
#sum(inlet.discharge)
#plot(GridX[cbind(inlet.x, inlet.y)], inlet.discharge)
#inlet.discharge


#jpeg(paste(QA.folder,"Boundary_Conditions_",site.list$SiteID[k],site.list$Year[k],".jpg",sep=""), 6,6, units='in', res=600)
jpeg(paste(QA.folder,"Boundary_Conditions.jpg",sep=""), 6,6, units='in', res=600)

################################################################
### Make plot for paper #####
################################################################
plot(data$X, data$Y, col="gray", pch=19, cex=.1,
 ylim=c(min(data$Y), min(data$Y)+1.2*(max(data$Y)-min(data$Y))),
 xlab="Easting (m)", ylab="Northing (m)",
 main=site.list$SiteID[k]
)

points(WSEDEM$ws.X, WSEDEM$ws.Y, col="cyan", pch=19, cex=.1)

lines(c(minX, maxX, maxX, minX, minX),
    c(minY, minY, maxY, maxY, minY), type="l", lt=3,lw=2, col="black")


legend("top", 
legend = c("Reach Extent", "Wetted Area", "CFD Grid Extent",
 "Thalweg", "Inlet Boundary", "Exit Boundary"), 
lt=c(1,1,3,1,1,1), lw=c(4,4,2,1,2,2), 
col=c("gray","cyan", "black", "black", "blue", "red"), bg="white",
bty="n",ncol=2, cex=.7)


lines(thalweg)

points(GridX[inlet.x, inlet.y], GridY[inlet.x, inlet.y],pch=19, cex=.7,
col="blue")


lines(GridX[outflow.x.idx, outflow.y.idx], GridY[outflow.x.idx, outflow.y.idx],
col="red", lw=4)
dev.off()
#####
# Done making plot for paper



######################## write input files

# test.src file
cat("", file="test.src")
for (i in 1:length(inlet.discharge)){
cat("discharge",i,"           Y      ", inlet.x[i], inlet.y[i], " 0    N","\n", file="test.src", append=T)
}


cat("", file="test.dis")
for (i in 1:length(inlet.discharge)){
# test.dis file
cat("table-name           'Discharge :",i,"'","\n", file="test.dis", append=T)
cat("contents             'regular   '", "\n", file="test.dis", append=T)
cat("location             'discharge",i,"    '","\n", file="test.dis", append=T)
cat("time-function        'non-equidistant'", "\n", file="test.dis", append=T)
cat("reference-time       20130101", "\n", file="test.dis", append=T)
cat("time-unit            'minutes'", "\n", file="test.dis", append=T)
cat("interpolation        'linear'", "\n", file="test.dis", append=T)
cat("parameter            'time                '                     unit '[min]'", "\n", file="test.dis", append=T)
cat("parameter            'flux/discharge rate '                     unit '[m3/s]'", "\n", file="test.dis", append=T)
cat("records-in-table     2", "\n", file="test.dis", append=T)
cat( 0.0000000e+000 , inlet.discharge[i],"\n", file="test.dis", append=T)
cat( simtime,   inlet.discharge[i], file="test.dis", "\n",append=T)

}
########################################################################
# Write Observations file (not sure why I need this.. not used, but
# needed to run model


cat("(5,5)                    5      5", "\n", file= "test.obs")
cat("(8,8)                    8      8", file= "test.obs", append=T)
##########################################################################





#############################################################
##### Write .grd file (grid file)
##### This just makes a rectangular grid in which the\
##### stream channel fits.  It does not fit a "best"
##### curvilinear orthogonal grid to the stream channel,
##### so the computation might end up being very inefficient
#############################################################


M=NX
N=NY


X = GridX
Y= GridY

# subtract minimum.  Need to pass smaller number (fewer digits) as
# grid input, to limit no. of significant digits.  Will output offset
# (and other info) to reference file, for use in post-processing

min(GridX)
min(GridY)

X = GridX-min(GridX)
Y= GridY-min(GridY)

# Write this file containing the offsets between true x and y locations
# and x and y locations used in the model.  Because the X and Y points
# have too many digits, I can't pass the entire X and Y locations
# to Delft 3D, or the last digits will be lost, and chaos ensues.  A
# pain in the butt, but that's how it needs to be for now....

write.csv(data.frame("X.offset"=min(GridX), "Y.offset"=min(GridY),
"inlet"=inlet, "outlet"=outlet),"offset.csv")

########################################################

inlet
write(
"*
* Test and Development Grid File
*", file = "test.grd")

write("Coordinate System = Cartesian", file = "test.grd",append=T)
write(paste(nrow(X), ncol(X)), file = "test.grd",append=T)
write("0 0 0",  file = "test.grd",append=T)


# Can't have list of points more than 5 per line....
# Do some data gymnastics here...
rows = trunc(nrow(Y)/5-.0000001)+1

r=1
# x's in columns, y's in rows...


for (j in 1:ncol(Y)) {
 cat("ETA=", j, " ", file=  "test.grd",append=T) 
for (r in 1:rows) {
  cat(round(X[((r-1)*5+1):min(nrow(Y),(r*5)),j],2),"\n", file=  "test.grd",append=T)
}} 



for (j in 1:ncol(Y)) {
 cat("ETA=", j," ", file=  "test.grd",append=T) 
for (r in 1:rows) {
  cat(round(Y[((r-1)*5+1):min(nrow(Y),(r*5)),j],2),"\n", file=  "test.grd",append=T)
}}



#################################3
# Bathymetry
# Define "B-Data".
# Note that Delft 3D requires one more dimension in X and Y than grid because geometry for X and
# Y is at the "centers" and Z is that the "nodes" of computationl cells.
# To be honest I ignore that and assign Z at the X and Y points.  Which means my results may actually be
# translated DX/2 in both the X and Y directions - but I don't really care if the entire set of results is
# shifted by, typically, 5 cmm or less in X and Y.  (Maybe I should, but for now, I'm ignoring that).

M= NX
N= NY

# Bathymetry data for input file has to be in an array of size M+1, N+1
Bdata = array(rep(-999.0, ((M+1)*(N+1))), c((M+1), (N+1)))
Bdata[1:M, 1:N] = depth
Bdata[M+1, 1:N] = Bdata[M, 1:N]
Bdata[,N+1] = Bdata[,N]

max(Bdata)
dim(Bdata)

# Can't have list of points more than 5 per line....
# Do some data gymnastics here...
Brows = trunc((M+1)/5-.0000001)+1
cat("", file = "test.dep")
for (j in 1:(N+1)){
for (r in 1:Brows) {
  cat(Bdata[((r-1)*5+1):min(nrow(Bdata),(r*5)),j],"\n", 
file=  "test.dep",append=T)
}
}

#############################################
# Define dry points (will be a lot for rectangular grids)
# Dry points are un-used computational points that will never have water in them.  
# Since I use a rectangular computational grid, but the actual sites are far from rectangular,
# There's usually a lot of "dry" points.  Defining them cuts down on computational time significantly.
# I don't have to define them tough....

cat("", file="test.dry")
for (m in 2:(M-1)){

# check to see if entire side is dry
if  ((sum(Bdata[m,2:(N-1)]>0) * 1) == 0){
 cat(m,"2", m," ", N-1,"\n", file="test.dry", append=T)
} else {

Bdata[3,]==0
if ((Bdata[m,2]==0) & (Bdata[m,1]==0)) {
 # bottom dry line
    cat(m,"2 ",m, " ", (match(F, Bdata[m,]==0)-1),"\n", file="test.dry", append=T)
  }

 # top dry line
if ((Bdata[m,(N-1)]==0) & (Bdata[m,N]==0))  {
 # bottom dry line
    cat(m," ",N-(match(F, Bdata[m,N:1]==0))+2," ",m," ", (N-1),"\n", file="test.dry", append=T)
  }
}
}


##########################################
# Define the computational grid enclosure
# Delft 3D requires an "enclosure" file.  Could get fancy, but for now, I just draw the
# rectangle of the computational grid.



cat(1, " ", 1, "  begin external polygon \n", file = "test.enc")
cat(M, " ", 1, "\n", file= "test.enc", append=T)
cat(M, " ", N, "\n", file= "test.enc", append=T)
cat(1, " ", N, "\n", file= "test.enc", append=T)
cat(1, " ", 1, "  end external polygon \n", file = "test.enc", append=T)



###################################################################################
############ Write the mdf file ####################
# The mdf file is the "main" input file for Delft 3D. It lists all the other input
# files, and also provides inputs for time step, DX and DY, HEV, etc; and tells
# Delft 3D which output files to write and when, etc.  To create this script, I just
# copied an mdf file created from the GUI, and re-wrote it here, changing the appropriate
# inputs as needed.  Note that all the input files are named "test.___" for all sites.  They
# are all written to the working directory, which is unique to each site - that's how different
# sites are differentiated from each other.

dt = .0025
if (simtime > 100) {dt = .005}
if (simtime > 200) {dt = .01}
if (simtime > 300) {dt = .025}
if (simtime > 500) {dt = .025}


# Write outputs for every 1 minute of simulated time, up to a maximum of 100 outputs.  
# It's sometimes useful to see solution progress, rather than just write last time step.
# Sometimes, if it's clear the solution has converged before the last time step, it's OK
# to stop the simulation early, and results will have been written.
# Note: There's a lot here that's un-needed, as far as what get's written to the output, tracking points,
# etc.  Some clean up to do.  I left it all in as it's probbly easier to cut stuff later than add stuff later.

output.interval = max(1, round(simtime / 5))
output.interval = max(1,round((simtime-5)/20))
#output.interval =.1

HEV = site.list$HEV[k]

# Write the mdf File!
cat("Ident  = #Delft3D-FLOW 3.43.05.22651#
Commnt =                  
Filcco = #test.grd#
Anglat =  0.0000000e+000
Grdang =  0.0000000e+000
Filgrd = #test.enc#
MNKmax = ",NX+1, NY+1, 1," 
Thick  =  1.0000000e+002
Commnt =                  
Fildep = #test.dep#
Commnt =                  
Commnt =                 no. dry points: 1437
Fildry = #test.dry#
Commnt =                 no. thin dams: 0
Commnt =                  
Itdate = #2013-01-01#
Tunit  = #M#
Tstart =  0.0000000e+000
Tstop  =  ", simtime,"
Dt     = ",dt,"
Tzone  = 0
Commnt =                  
Sub1   = #    #
Sub2   = #   #
Commnt =                  
Wnsvwp = #N#
Wndint = #Y#
Commnt =       
Zeta0  = ",outflow.ws.level,"
Commnt =                  
Commnt =                 no. open boundaries: 1
Filbnd = #test.bnd#
FilbcT = #test.bct#
Commnt =                  
Ag     =  9.8100000e+000
Rhow   =  1.0000000e+003
Tempw  =  1.5000000e+001
Salw   =  3.1000000e+001
Wstres =  6.3000000e-004  0.0000000e+000  7.2300000e-003  1.0000000e+002  7.2300000e-003  1.0000000e+002
Rhoa   =  1.0000000e+000
Betac  =  5.0000000e-001
Equili = #N#
Ktemp  = 0
Fclou  =  0.0000000e+000
Sarea  =  0.0000000e+000
Temint = #Y#
Commnt =                  
Roumet = #W#
Ccofu  =  ", Ccofu,"
Ccofv  =  ", Ccofv,"
Xlo    =  0.0000000e+000
Vicouv =  ", HEV, "
Dicouv =  1.0000000e+001
Htur2d = #N#
Irov   = 0
Commnt =                  
Iter   =      2
Dryflp = #YES#
Dpsopt = #MAX#
Dpuopt = #MEAN#
Dryflc =  1.0000000e-002
Dco    = -9.9900000e+002
Tlfsmo =  5.0000000e+000
ThetQH =  0.0000000e+000
Forfuv = #Y#
Forfww = #N#
Sigcor = #N#
Trasol = #Cyclic-method#
Momsol = #Cyclic#
Commnt =                  
Commnt =                 no. discharges:", length(inlet.discharge),"
Filsrc = #test.src#
Fildis = #test.dis#
Commnt =                 no. observation points: 2
Filsta = #test.obs#
Commnt =                 no. drogues: 0
Commnt =                  
Commnt =                  
Commnt =                 no. cross sections: 0
Commnt =                  
SMhydr = #YYYYY#     
SMderv = #YYYYYY#    
SMproc = #YYYYYYYYYY#
PMhydr = #YYYYYY#    
PMderv = #YYY#       
PMproc = #YYYYYYYYYY#
SHhydr = #YYYY#      
SHderv = #YYYYY#     
SHproc = #YYYYYYYYYY#
SHflux = #YYYY#      
PHhydr = #YYYYYY#    
PHderv = #YYY#       
PHproc = #YYYYYYYYYY#
PHflux = #YYYY#      
Online = #N#
Waqmod = #N#
Flmap  =  0.0000000e+000 ",output.interval," ",simtime,"
Flhis  =  0.0000000e+000 ",output.interval," ",simtime,"
Flpp   =  0.0000000e+000 ",output.interval," ",simtime,"
Flrst  = 1
Commnt =                  
Commnt =   ",
file= "test.mdf")

detach(thalweg)       
        



########################################################################################
# Write xml file, which is needed to run from command line prompt

cat("<?xml version=\"1.0\" encoding=\"iso-8859-1\"?>
<deltaresHydro xmlns=\"http://schemas.deltares.nl/deltaresHydro\" xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\" xsi:schemaLocation=\"http://schemas.deltares.nl/deltaresHydro http://content.oss.deltares.nl/schemas/d_hydro-1.00.xsd\">
    <control>
        <sequence>
            <start>myNameFlow</start>
        </sequence>
    </control>
    <flow2D3D name=\"myNameFlow\">
        <library>flow2d3d</library>
        <mdfFile>test.mdf</mdfFile>
    </flow2D3D>
    <delftOnline>
        <enabled>true</enabled>
        <urlFile>test.url</urlFile>
        <waitOnStart>false</waitOnStart>
        <clientControl>true</clientControl>    <!-- client allowed to start, step, stop, terminate -->
        <clientWrite>false</clientWrite>    <!-- client allowed to modify data -->
    </delftOnline>
</deltaresHydro>",
file= "test.xml")
#####################################################################


################
# Create MetaData file
# Surface Roughness

names(site.list)
meta.data = list(
"VisitID"=site.list$VisitID[k],
"D84"=site.list$D84[k],
"Roughness Input" = Ccofu,
"HEV" = HEV,
"TrimLength"=slop,
"Discharge" = discharge, #modeled discharge
"Measured_Discharge" = site.list$Measured.Discharge[k],
"Year"=site.list$Year[k],
"WatershedName"=site.list$WatershedName[k],
"Exit_BC"=outflow.ws.level,
"Left.Reference" =  ref.left,
"Top.Reference" = ref.top,
"Comp.Grid.Spacing" = DX,
"Build.Input.File.R.Version" = Build.Input.File.R.Version,
"Delft3D.Version" = Delft3D.Version,
"Operator" = Operator,
"Pre.Processing.Date.Time" = Sys.time(),
"rbt.version" = rbt.version
)

names(site.list)
write.csv(meta.data,paste(site.list$D3D.Input.Folder[k],"/Meta.Data.csv", sep=""))


WorkingDir
}

# Done!
# Files should be ready to run in Delft 3D.


setwd(savedwd)
