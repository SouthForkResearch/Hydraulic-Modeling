
site.list$Results.Folder
BC_Dir=("C://Matt-SFR Files//Hydraulic Modeling//BC Copies")
dir.create(BC_Dir)

names(site.list)
i=1
for (i in 1:length(site.list$Results.Folder)){
files=dir(site.list$Results.Folder[i])
if (site.list$Model[i] == "Yes") {
if ("Boundary_Conditions.jpg" %in% files) {

copyfile1=paste(site.list$Results.Folder[i], "//Boundary_Conditions.jpg", sep="")
copyfile2=paste(site.list$Results.Folder[i], "//Depth.Error.jpg", sep="")
copyfile3=paste(site.list$Results.Folder[i], "//Velocity.Magnitude.jpg", sep="")
copyfile4=paste(site.list$Results.Folder[i], "//Depth.jpg", sep="")


pastefile1 = paste(BC_Dir,"//Boundary_Conditions_",i,".jpg", sep="")
pastefile2 = paste(BC_Dir,"//Depth_Error_",i,".jpg", sep="")
pastefile3 =  paste(BC_Dir,"//Velocity_Magnitude_",i,".jpg", sep="")
pastefile4 =  paste(BC_Dir,"//Depth_",i,".jpg", sep="")


file.copy(copyfile1, pastefile1, overwrite=T, copy.mode = TRUE, copy.date = FALSE)
file.copy(copyfile2, pastefile2, overwrite=T, copy.mode = TRUE, copy.date = FALSE)
file.copy(copyfile3, pastefile3, overwrite=T, copy.mode = TRUE, copy.date = FALSE)
file.copy(copyfile4, pastefile4, overwrite=T, copy.mode = TRUE, copy.date = FALSE)

}
}
}

