**This research was supported under NSF Grant AGS 12-41407 to The Pennsylvania State University.**


**MAIN B**: Code to Read, Prepare, and Analyze GPS-TEC Data: **FUN_dataViewer.m**

This code imports GPS-TEC data and prepares it for analysis in a variety of ways. Flags at the beginning of the code control whether the radius averaging (a point-approximation) or the averaging slices algorithms are to be used. The GPS-TEC dataset is not directly downloadable. It was developed on MATLAB R2016b. For older MATLAB versions an alternate HDF5 file reading code is included but commented out that should work if the currently implemented HDF5 file reading code causes issues.

For the ISR dataset support, the code in **MSThesis_Appendicies/MSThesis_A_ISR_Dataset_Analysis/Ross_test.m** must be run and the variables *MISA_SNR*, *MISA_SNR_bp*, *MISA_height*, *MISA_time*, *MISA_vel*, *Zenith_SNR*, *Zenith_SNR_bp*, *Zenith_height*, *Zenith_time*, and *Zenith_vel* must be saved with the same names with the **.mat** file ending.

For the geomagnetic latitude and longitude support, this code relies on an external function set named *IGRF Magnetic Field* by Kip Knight on the MATLAB File Exchange. It can be found at: https://www.mathworks.com/matlabcentral/fileexchange/45606-igrf-magnetic-field?focused=3820319&tab=function. The two functions from that function set that are needed are IGRF.m and GEODIP.m. It is recommended to generate geomagnetic coordinates for the entire world (-90 to 90 arcdeg latitude and -180 to 180 arcdeg longitude) with the command to create a set of pplat and pplong (*[pplat,pplong] = sFUN_geoToGeomag(dateOfData,magF,pplat,pplong);*) and then save the geomagnetic pplat and pplong as **magcoord_pplat_126127128.mat** and **magcoord_pplong_126127128.mat**. Creating the geomagnetic coordinates takes much longer than loading them from a pre-prepared file.

For the AMPERE/JEDI Heating data, there is no support as the datasets are preliminary and not available.


*Supporting functions are also provided so the code can be run.*


**FUN B-1:** Function to Convert Geographic to Geomagnetic Coordinates: **sFUN_geoToGeomag.m**

This function converts geographic coordinates to geomagnetic coordinates using an external function set for the conversion process. It was developed on MATLAB R2016B.


**FUN B-2:** Function to get Kp Index Data: **sFUN_KpIndexGET.m**

This function gets Kp index data for the requested days. Requires internet access. It was developed on MATLAB R2016B.


**FUN B-3:** Function to get OMNI Data: **sFUN_OMNIGET.m**

This function gets OMNI data for the requested days. Requires internet access. It was developed on MATLAB R2016B.


**FUN B-4:** Function to Convert Date to Day Number: **sFUN_dateToDayNum.m**

This function converts a given date to day number (the cumulative number of days into a year the date is). It was developed on MATLAB R2016B.


**FUN B-5:** Function to Convert Date to Number of Days in the Month: **sFUN_dateToNumDaysInTheMonth.m**

This function converts a given date to the number of days in the month. It was developed on MATLAB R2016B.


**FUN B-6:** Function to Convert Day Number to its Date: **sFUN_dayNumber_to_Date.m**

This function converts day number (the cumulative number of days into a year the date is) to the current date. It was developed on MATLAB R2016B.


**FUN B-7:** Function to Convert Day Number to Dates for Multiple Day Numbers: **sFUN_dayNumber_to_Date_MULTIPLE.m**

This function converts a consecutive range of day numbers (the cumulative number of days into a year the date is) to their respective dates. It was developed on MATLAB R2016B.


**FUN B-8:** Function to Convert Month Number to Month Name: **sFUN_monthNum_to_word.m**

This function converts a month number to the month name, e.g. 2 becomes February, Feb, or Feb. depending on options. It was developed on MATLAB R2016B.
