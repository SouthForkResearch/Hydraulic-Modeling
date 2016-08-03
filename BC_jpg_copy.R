
site.list$Results.Folder
BC_Dir=("C://Matt-SFR Files//Hydraulic Modeling//BC Copies")
dir.create(BC_Dir)

i=1
for (i in 1:length(site.list$Results.Folder)){
files=dir(site.list$Results.Folder[i])
if ("Boundary_Conditions.jpg" %in% files) {

copyfile=paste(site.list$Results.Folder[i], "//Boundary_Conditions.jpg", sep="")
pastefile = paste(BC_Dir,"//Boundary_Conditions_",i,".jpg", sep="")

file.copy(copyfile, pastefile, overwrite=T, copy.mode = TRUE, copy.date = FALSE)
}
}

