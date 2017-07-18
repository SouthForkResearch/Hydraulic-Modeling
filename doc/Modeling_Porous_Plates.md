# CHaMP Hydraulic Modeling of Porous Plates

CHaMP hydraulic modeling efforts have been expanded to exploit the porous plate feature of the Delft3D software in order to roughly approximate the effect of porous structures in flows such as large wood jams and beaver dams.  As of July 2017, this is a manual process, requiring a separate set of input files that describe the location and porosity of each wood structure, and a different pre-processing R script must be used.  Limited validation works done to date suggests that including porous structures in the hydro model significantly improves the agreement between modeled and measured results, as compared to simply ignoring the porous structures in the hydraulic models.



describe assumptions, what sort of jams can be modeled



## Input files

Two extra input comma separated value (.csv) files are required to specify the location and porosity value of each jam modeled.  These must then be placed in the folder with the other, default input files (DEM.csv, WSEDEM.csv, and Thalweg.csv).  The first must be called "Jam_Locations.csv"  It must contain the x and y locations, on a 10 cm grid (the same grid as the DEM.csv file), for all points at each porous structure to be modeled.  The "Value" field is then used as an identifier to label the individual porous structure.  An example of a portion of an input file is shown below.  For more information on generating the Jam_Locations.csv input file from field surveys or for simulated restoration structures, contact Matt Nahorniak, Kelly Whitehead, or Andrew Hill.

Jam_Locations.csv

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



alternate pre-processing script

extra QA plot

results

