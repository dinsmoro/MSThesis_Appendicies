**This research was supported under NSF Grant AGS 12-41407 to The Pennsylvania State University.**


**MAIN C:** Code to Create and Analyze Synthetic GPS-TEC Data: **test_SyntheticNoiseANDWave.m**

This code creates synthetic GPS-TEC data for given latitude/longitude/time coordinates and then analyzes that synthetic data. It was developed on MATLAB R2016B.

This code relies on the code in **MSThesis_Appendicies/MSThesis_B_GPS-TEC_Dataset_Analysis/FUN_dataViewer.m** to be run for a latitude range of 35 to 50 arcdeg and a longitude range of -85 to -60 arcdeg. The resulting *pplat*, *pplong*, *time*, and *timeUnique* variables should then be saved into a combined file named **sTEC_supporting_values.mat**. 

This code has a setting for geomagnetic coordinates, but support for the geomagnetic coordinates is not fully implemented.
