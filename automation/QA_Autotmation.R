# First, download QA files as follows - or with some variant of....
# Note - make sure to change watershed in both sides of sync command (if including watershed)

# aws s3 sync s3://sfr-champdata/QA/2011/Methow "C:\Matt-SFR Files\Hydraulic Modeling\QA\champ data from bucket\2011\Methow" --exclude "*" --include "*\Depth.Error.jpg" --include "*\Boundary_Conditions.jpg" --exclude "*\Results\*" --exclude "*\QA_Reject\*" --exclude "*\ResultsZZ\*"
# aws s3 sync s3://sfr-champdata/QA/2011/UpperGrandeRonde "C:\Matt-SFR Files\Hydraulic Modeling\QA\champ data from bucket\2011\UpperGrandeRonde" --exclude "*" --include "*\Depth.Error.jpg" --include "*\Boundary_Conditions.jpg" --exclude "*\Results\*" --exclude "*\ResultsZZ\*" --exclude "*\QA_Reject\*"
# aws s3 sync s3://sfr-champdata/QA/2011/JohnDay "C:\Matt-SFR Files\Hydraulic Modeling\QA\champ data from bucket\2011\JohnDay" --exclude "*" --include "*\Depth.Error.jpg" --include "*\Boundary_Conditions.jpg" --exclude "*\Results\*" --exclude "*\ResultsZZ\*" --exclude "*\QA_Reject\*"
# aws s3 sync s3://sfr-champdata/QA/2011/Entiat "C:\Matt-SFR Files\Hydraulic Modeling\QA\champ data from bucket\2011\Entiat" --exclude "*" --include "*\Depth.Error.jpg" --include "*\Boundary_Conditions.jpg" --exclude "*\Results\*" --exclude "*\ResultsZZ\*" --exclude "*\QA_Reject\*"
# aws s3 sync s3://sfr-champdata/QA/2011/Wenatchee "C:\Matt-SFR Files\Hydraulic Modeling\QA\champ data from bucket\2011\Wenatchee" --exclude "*" --include "*\Depth.Error.jpg" --include "*\Boundary_Conditions.jpg" --exclude "*\Results\*" --exclude "*\ResultsZZ\*" --exclude "*\QA_Reject\*"
# aws s3 sync s3://sfr-champdata/QA/2011/Tucannon "C:\Matt-SFR Files\Hydraulic Modeling\QA\champ data from bucket\2011\Tucannon" --exclude "*" --include "*\Depth.Error.jpg" --include "*\Boundary_Conditions.jpg" --exclude "*\Results\*" --exclude "*\ResultsZZ\*" --exclude "*\QA_Reject\*"

# aws s3 sync s3://sfr-champdata/QA/2011/YankeeFork "C:\Matt-SFR Files\Hydraulic Modeling\QA\champ data from bucket\2011\YankeeFork" --exclude "*" --include "*\Depth.Error.jpg" --include "*\Boundary_Conditions.jpg" --exclude "*\Results\*" --exclude "*\ResultsZZ\*" --exclude "*\QA_Reject\*"
# aws s3 sync s3://sfr-champdata/QA/2011/Lemhi "C:\Matt-SFR Files\Hydraulic Modeling\QA\champ data from bucket\2011\Lemhi" --exclude "*" --include "*\Depth.Error.jpg" --include "*\Boundary_Conditions.jpg" --exclude "*\Results\*" --exclude "*\ResultsZZ\*" --exclude "*\QA_Reject\*"
# aws s3 sync s3://sfr-champdata/QA/2011/SouthForkSalmon "C:\Matt-SFR Files\Hydraulic Modeling\QA\champ data from bucket\2011\SouthForkSalmon" --exclude "*" --include "*\Depth.Error.jpg" --include "*\Boundary_Conditions.jpg" --exclude "*\Results\*" --exclude "*\ResultsZZ\*" --exclude "*\QA_Reject\*"

# aws s3 sync s3://sfr-champdata/QA/2011 "C:\Matt-SFR Files\Hydraulic Modeling\QA\champ data from bucket\2011" --exclude "*" --include "*\Depth.Error.jpg" --include "*\Boundary_Conditions.jpg" --exclude "*\Results\*"  --exclude "*\ResultsZZ\*"--exclude "*\QA_Reject\*"
# aws s3 sync s3://sfr-champdata/QA/2011 "C:\Matt-SFR Files\Hydraulic Modeling\QA\champ data from bucket\2011" --exclude "*" --include "*\Depth.Error.jpg" --include "*\Boundary_Conditions.jpg" --exclude "*\Results\*"  --exclude "*\ResultsZZ\*"--exclude "*\ResultsZZ\*" --exclude "*\QA_Reject\*"

# aws s3 sync s3://sfr-champdata/QA/2012 "C:\Matt-SFR Files\Hydraulic Modeling\QA\champ data from bucket\2012" --exclude "*" --include "*\Depth.Error.jpg" --include "*\Boundary_Conditions.jpg" --exclude "*\Results\*" --exclude "*\ResultsZZ\*" --exclude "*\QA_Reject\*"
# aws s3 sync s3://sfr-champdata/QA/2013 "C:\Matt-SFR Files\Hydraulic Modeling\QA\champ data from bucket\2013" --exclude "*" --include "*\Depth.Error.jpg" --include "*\Boundary_Conditions.jpg" --exclude "*\Results\*" --exclude "*\ResultsZZ\*" --exclude "*\QA_Reject\*"
# aws s3 sync s3://sfr-champdata/QA/2014 "C:\Matt-SFR Files\Hydraulic Modeling\QA\champ data from bucket\2014" --exclude "*" --include "*\Depth.Error.jpg" --include "*\Boundary_Conditions.jpg" --exclude "*\Results\*" --exclude "*\ResultsZZ\*" --exclude "*\QA_Reject\*"
# aws s3 sync s3://sfr-champdata/QA/2015 "C:\Matt-SFR Files\Hydraulic Modeling\QA\champ data from bucket\2015" --exclude "*" --include "*\Depth.Error.jpg" --include "*\Boundary_Conditions.jpg" --exclude "*\Results\*" --exclude "*\ResultsZZ\*" --exclude "*\QA_Reject\*"
# aws s3 sync s3://sfr-champdata/QA/2016 "C:\Matt-SFR Files\Hydraulic Modeling\QA\champ data from bucket\2016" --exclude "*" --include "*\Depth.Error.jpg" --include "*\Boundary_Conditions.jpg" --exclude "*\Results\*" --exclude "*\ResultsZZ\*" --exclude "*\QA_Reject\*"

# aws s3 sync s3://sfr-champdata/QA "C:\Matt-SFR Files\Hydraulic Modeling\QA\champ data from bucket" --exclude "*" --include "*\Depth.Error.jpg" --include "*\Boundary_Conditions.jpg" --exclude "*\Results\*" --exclude "*\ResultsZZ\*" --exclude "*\QA_Reject\*"


# Load this function, then manually enter "QA()" in the console.  For reasons I don't understand, the "readline()" function won't stop and
# wait for the user to enter data if it's not run either from a command line, or inside a function that's called from the command line.  


QA = function() {
library(jpeg)

root_dir = "c:/Matt-SFR Files/Hydraulic Modeling/R Code to Build Input Files/R-Code/automation"

dir()
savedwd =  "C:/Matt-SFR Files/Hydraulic Modeling/QA/champ data from bucket/"
setwd(savedwd)


dirs = list.files(recursive=T, pattern = "Boundary_Conditions.jpg")
dirs = sub("Boundary_Conditions.jpg", "", dirs)
dirs
print("length dirs =")
print(length(dirs))
dirs[1]
VisitIDs = rep("", length(dirs))
JobString = rep("", length(dirs))
RM_old_results = rep("", length(dirs))
RM_old_rejects = rep("", length(dirs))
Mv_results = rep("", length(dirs))
Mv_rejects = rep("", length(dirs))
Mv_command = rep("", length(dirs))

dirs

i=1
for (i in 1:length(dirs)){
A2=strsplit(dirs[i], "VISIT_")[[1]][2]
A1=strsplit(dirs[i], "VISIT_")[[1]][1]
A1
A2
VisitIDs[i] = (strsplit(A2,"/")[[1]][1])
VisitIDs[i]
VisitIDs

B=strsplit(dirs[i], "QA/")[[1]][2]
B=strsplit(B,"/")[[1]][1]
B
JobString[i] = B


RM_old_results[i] = paste("aws s3 rm --recursive s3://sfr-champdata/QA/",A1,"VISIT_",VisitIDs[i],"/Hydro/Results/", B, sep="")
RM_old_rejects[i] = paste("aws s3 rm --recursive s3://sfr-champdata/QA/",A1,"VISIT_",VisitIDs[i],"/Hydro/QA_Reject/", B, sep="")

Mv_results[i]=paste("aws s3 mv --recursive s3://sfr-champdata/QA/",A1,"VISIT_",VisitIDs[i],"/Hydro/QA/",B,
" s3://sfr-champdata/QA/",A1,"VISIT_",VisitIDs[i],"/Hydro/Results/",B, sep="")
Mv_results[i]

Mv_rejects[i]=paste("aws s3 mv --recursive s3://sfr-champdata/QA/",A1,"VISIT_",VisitIDs[i],"/Hydro/QA/",B,
" s3://sfr-champdata/QA/",A1,"VISIT_",VisitIDs[i],"/Hydro/QA_Reject/",B, sep="")

}
RM_old_results
RM_old_rejects
Mv_rejects
Mv_results

JobString
VisitIDs=as.numeric(VisitIDs)

QA_status = rep("fail", length(VisitIDs))


# Note: there are 844 visitIDs to check.  I may want to do these in
# smaller chunks

low = 1
high = length(dirs)
for (i in low:high){
#for (i in 1:length(VisitIDs)){
print(paste("i =",i,"VisitID =", VisitIDs[i]))
setwd(savedwd)
savedwd
dirs[i]
setwd(dirs[i])


jfile = readJPEG("Boundary_Conditions.jpg")
plot(1:2, type='n')
rasterImage(jfile, 1, 1, 2, 2)
print("QA OK?")
yQA_1= readline()

jfile = readJPEG("Depth.Error.jpg")
plot(1:2, type='n')
rasterImage(jfile, 1, 1, 2, 2)

print("QA OK?")
yQA_2= readline()


if ((yQA_1=="y") & (yQA_2=="y")) {QA_status[i] = "pass"}
if (QA_status[i] == "pass") {Mv_command[i] = Mv_results[i]} else {Mv_command[i] = Mv_rejects[i]}


setwd(savedwd)

}

## Now delete local folders so I never try to move these again!
do.call(file.remove, list(list.files(dirs, full.names = TRUE)))


VIDs = rep(VisitIDs[low:high], 3)
Jobs = rep(JobString[low:high], 3)
Bat = c(RM_old_results[low:high], RM_old_rejects[low:high], Mv_command[low:high])
results=data.frame(VIDs, Jobs, Bat)

#results=data.frame("VisitID"=VisitIDs, "Job"=JobString, QA_status, RM_old_results, RM_old_rejects, Mv_command)
#results=results[low:high,]


name=paste(root_dir, "/QA_results_",low,"_",high,".csv",sep="")
write.csv(results, name)
} # end of QA function

