# CHaMP Hydraulic Modeling of Porous Plates

CHaMP hydraulic modeling efforts have been expanded to exploit the porous plate feature of the Delft3D software in order to roughly approximate the effect of porous structures in flows such as large wood jams and beaver dams.  As of July 2017, this is a manual process, requiring a separate set of input files that describe the location and porosity of each wood structure, and a different pre-processing R script must be used.  Limited validation works done to date suggests that including porous structures in the hydro model significantly improves the agreement between modeled and measured results, as compared to simply ignoring the porous structures in the hydraulic models.

Because we're modeling in 2 dimensions, any porous structures modeled as assumed to be uniform in the x-y plane and, essentially, infinitely tall.  Water can flow through them and/or around them, but not over them.  Structures such as fallen logs, where water flows not only around but also below and/or over the top of, are generally not well suited for porous structure modeling.  



## Input files

Two extra input comma separated value (.csv) files are required to specify the location and porosity value of each jam modeled.  These must then be placed in the folder with the other, default input files (DEM.csv, WSEDEM.csv, and Thalweg.csv).  The first must be called "Jam_Locations.csv"  It must contain the x and y locations, on a 10 cm grid (the same grid as the DEM.csv file), for all points at each porous structure to be modeled.  The "Value" field is then used as an identifier to label the individual porous structure.  An example of a portion of an input file is shown below.  For more information on generating the Jam_Locations.csv input file from field surveys or for simulated restoration structures, contact Matt Nahorniak, Kelly Whitehead, or Andrew Hill.

### Jam_Locations.csv

| X         | Y          | Value |
| --------- | ---------- | ----- |
| 436033.95 | 5003011.55 | 4     |
| 436034.05 | 5003011.55 | 4     |
| 436034.15 | 5003011.55 | 4     |
| 436034.25 | 5003011.55 | 4     |
| 436034.35 | 5003011.55 | 4     |
| 436034.45 | 5003011.55 | 4     |
| 436034.55 | 5003011.55 | 4     |
| 436034.65 | 5003011.55 | 4     |
| 436034.75 | 5003011.55 | 4     |
| 436034.85 | 5003011.55 | 4     |
| 436034.95 | 5003011.55 | 4     |
| 436035.05 | 5003011.55 | 4     |
| 436031.75 | 5003011.45 | 4     |
| 436031.85 | 5003011.45 | 4     |



The 2nd input file must be called "Jam_Porosity.csv".  This file is simply a list of each jam to be modeled, and the approximate porosity of that jam.  porosity is defined, roughly, as the 1 minus the proportion of the volume of the jam that is filled by solid material.  Contact Matt Nahorniak, Roby Ventres Pake, or Andrew Hill for more information of field collection of porosity estimates and field protocols.



### Jam_Porosity.csv

| Jam  | Porosity_Pct |
| ---- | ------------ |
| 1    | 70           |
| 2    | 80           |
| 3    | 50           |
| 4    | 90           |
| 5    | 70           |
| 6    | 80           |
| 7    | 70           |
| 8    | 60           |
| 9    | 70           |



An alternate pre-processing script called "Build_Input_Files_Porous_Plates.R" is used instead of the default "Build_Input_Files.R" script.  This script will read the extra csv input files and generate additional Delft3D input files required to recognize and run the model with the specified porous structures.  This file will also generate an additional QA file called "porous_plates.jpg" that shows the locations of the user specified porous structures.  This file should always be checked to ensure the modeled porous plates are in the expected locations.

The remainder of the hydro modeling steps for porous structure modeling are exactly the same as default CHaMP hydraulic modeling, and the same set of results files are generated. 

## Comparisons with and without porous structures

Often users may wish to compare hydraulic model results with and without porous structures.  For example, users may wish to model the effects of simulated restoration and compare these to an un-restored state.  The same set of input files can be used to model both cases - simply run the model with porous structures as described above, and model the no-porous structure model using the default pre-processing file "Build_Input_Files.R", which will ignore all the porous structure inputs.  Upon completion of each model, users will want to move the results files to a separate location, so that each model run doesn't overwrite results from the previous run.

[HYDRAULIC MODEL HOME](README.md)
