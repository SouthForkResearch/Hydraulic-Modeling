# CHaMP Hydraulic Modeling at Non-Default Flow Rates

Modeling non-default flow rates requires only one additional step beyond that used for default CHaMP hydraulic modeling: specifying the alternate flow rate in the file "CFD_Site_List.csv".  As part of the default modeling process, this file is generated, and may contain one row or multiple rows.  For each discharge to be modeled at each site/VisitID, specify the alternate discharge in the "Modeled.Discharge" column.  The units for this input are cubic meters per second (cms)

Running at an alternate discharge will necessarily create a boundary condition error at the downstream boundary of the modeled flow.  Generally we simply accept and ignore this, as these errors tend to be small for high gradient flows.  However, there are instances where the modeler may have information on the change in water depth at or near the exit boundary between the measured discharge and the modeled discharge.  If this information is known, it can be specified in the column "DeltaBC".  These units are in meters.

![cfd_sitelist](C:\Matt-SFR Files\Hydraulic Modeling\R Code to Build Input Files\R-Code\doc\cfd_sitelist.jpg)

Note that folders containing results from non-default flow rates with begin with an "M" for "modeled" instead of an "S" for surveyed.

Also note that the normal QA process for default hydraulic modeling includes examination of the "Depth Error" plot that compares surveyed depth and modeled depth.  In the case on modeling non-measured flow rates, we should not expect measured depth to be equal do modeled depth, and should therefore ignore the depth error QA plots.

