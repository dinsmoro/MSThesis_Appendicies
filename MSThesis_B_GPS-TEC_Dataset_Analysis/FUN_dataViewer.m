clc
clear variables
close all

%% Parallel Computing
warning('off','parallel:convenience:RunningJobsExist'); %my coding takes no prisoners
CPU_Threads = 4; %Number of worker threads to use
myCluster = parcluster('local'); %work around to force the CPU-threads requested
myCluster.NumWorkers = CPU_Threads; %forcing it
saveProfile(myCluster);
poolobj = gcp('nocreate'); % If no pool, do not create new one.
if isempty(poolobj)
    parpool(CPU_Threads); %Start new matlab pool
else
    if(poolobj.NumWorkers ~= CPU_Threads) %Check to see if open w/ desired worker num
        delete(gcp('nocreate')); %closes it if one is open w/ undersired size
        parpool(CPU_Threads); %Start new matlab pool
    end
end
clear poolobj myCluster; %keep the excess variables down

%% Turn on or off various Data Crunching and Visualizers

%convert the lat/longitude measurements to magnetic lat/long measurements
FLG_geomagneticCoords = 0; %1 for on, 0 for off

%Replace TEC data with random data on a Gaussian curve that's around the real data
FLG_randomData = 0; %1 for on, 0 for off

%plot data number visualizations
FLG_dataNumber = 0; %1 for on, 0 for off

%plot solar/geomagnetic index stuff
FLG_solarGeoStuff = 0; %1 for on, 0 for off
%Load JEDI data and prep for plotting
FLG_JEDI_DATA = 0;%1 for on, 0 for off
%Load AMPERE data and prep for plotting (JEDI data is heavily computed AMPERE)
FLG_AMPERE_DATA = 1; %1 for on, 0 for off

%gather stec at requested points for point analysis
FLG_gatherDataAtPoint = 0; %1 for on, 0 for off
%filter sTEC gathered at requested points with a smoothing system
FLG_gatherDataAtPoint_Filter = 0; %1 for on, 0 for off
%plot filtered sTEC results for all requested points
FLG_gatherDataAtPoint_Filter_Plot = 0; %1 for on, 0 for off
%Lomb-Scargle (alt FFT) on filtered sTEC results for all requested points
FLG_gatherDataAtPoint_Filter_Scargle = 0; %1 for on, 0 for off
%band pass filtered sTEC gathered at requested points (does nothing because was revealed sTEC was already low passed apparently)
FLG_gatherDataAtPoint_Filter_BP = 0; %1 for on, 0 for off
%Lomb-Scargle (alt FFT) on filtered bandpassed sTEC results for all requested points
FLG_gatherDataAtPoint_Filter_BP_Scargle = 0; %1 for on, 0 for off


FLG_gatherDataAtPoint_Filter_TimeMatch_n_Smooth =0; %1 for on, 0 for off

FLG_gatherDataAtPoint_Filter_TimeMatch_n_Smooth_ISRCorr = 0; %1 for on, 0 for off

FLG_gatherDataAtPoint_Filter_TimeMatch_n_Smooth_BP = 0; %1 for on, 0 for off

FLG_gatherDataAtPoint_Filter_TimeMatch_n_Smooth_Scargle = 0; %1 for on, 0 for off

FLG_gatherDataAtPoint_Filter_TimeMatch_n_Smooth_BP_PlotsScargle = 0; %1 for on, 0 for off
%Walking Lomb-Scargle (alt FFT) covering X hour time chunks on filtered sTEC results with a truncated time span to match ISR time span
FLG_gatherDataAtPoint_Filter_ISRTimeMatch_Scargle_Walking = 0; %1 for on, 0 for off - MANY PLOTS

%Zenith load ISR data for plotting later and prep it up in various ways
FLG_ISRdata_Zenith = 1;%1 for on, 0 for off
%Zenith plot direct ISR to sTEC filtered gathered at same point and more
FLG_ISRdata_Zenith_Plot = 0;%1 for on, 0 for off
%Zenith plot direct ISR to sTEC filtered BP gathered at same point and more
FLG_ISRdata_Zenith_Plot_BP = 0; %1 for on, 0 for off
%Zenith correlation between ISR and sTEC filtered gathered at same point
FLG_ISRdata_Zenith_Corr = 0; %1 for on, 0 for off
%Zenith ISR data Lomb-Scargled (alt FFT)
FLG_ISRdata_Zenith_Scargle = 0;%1 for on, 0 for off

%MISA load ISR data for plotting later and prep it up in various ways
FLG_ISRdata_MISA = 1;%1 for on, 0 for off
%MISA plot direct ISR to sTEC gathered at same point and more
FLG_ISRdata_MISA_Plot = 0;%1 for on, 0 for off
%MISA plot direct ISR to sTEC filtered BP gathered at same point and more
FLG_ISRdata_MISA_Plot_BP = 0; %1 for on, 0 for off
%MISA correlation between ISR and sTEC filtered gathered at same point
FLG_ISRdata_MISA_Corr = 0; %1 for on, 0 for off
%MISA ISR data Lomb-Scargled (alt FFT)
FLG_ISRdata_MISA_Scargle = 0; %1 for on, 0 for off



%AVG along Longitude sTEC plotted vs Latitude and MISA ISR
FLG_AVG_LONG = 0; %1 for on, 0 for off
%AVG along Latitude sTEC plotted vs Longitude and MISA ISR
FLG_AVG_LAT = 0; %1 for on, 0 for off
%AVG along Longitude sTEC Lomb-Scargled (alt FFT) at Millstone latitude
FLG_AVG_LONG_Scargle = 0;% 1 for on, 0 for off
%AVG along Longitude sTEC ISR Time Matched and BandPassed then Scargled @ Millstone Latitude
FLG_AVG_LONG_TimeMatch_n_BP_Scargle = 0;% 1 for on, 0 for off

%ANY ANGLE INFO
FLG_AVG_ANYANGLE = 0; %1 for on, 0 for off
FLG_AVG_ANYANGLE_PLOT = 0; %1 for on, 0 for off
FLG_AVG_ANYANGLE_PLOT_TIMECUTOUT = 0; %1 for on, 0 for off. Uses the time cutout range defined below
FLG_AVG_ANYANGLE_ZENITHORMISA = 0; %0 for plot vs Zenith ISR, 1 for plot vs MISA ISR
FLG_AVG_ANYANGLE_TimeMatch_n_BP_Scargle = 0; %1 for on, 0 for off. Scargles the data
FLG_AVG_ANYANGLE_Scargle_FFT = 0; %0 for Scargle, 1 for FFT
FLG_AVG_ANYANGLE_FULLTIME = 1; %0 for match to Zenith/MISA (norm), 1 for full time (will break stuff later I bet)

%NOT WHAT IT DOES more for ISR comparisons ancillary
FLG_AVG_LONG_TimeMatch_n_BP_Velocity = 0; %1 for on, 0 for off

%enable movie creation
FLG_enable_movieCreation = 1; %1 for on, 0 for off 


%% Prep Stuff
fileName = ['tid_new_126_2013.h5';'tid_new_127_2013.h5';'tid_new_128_2013.h5']; %file names of data
fileName = ['tid_126_2013_w241_n01_e30.h5';'tid_127_2013_w241_n01_e30.h5';'tid_128_2013_w241_n01_e30.h5']; %file names of data
% fileName = 'tid_new_128_2013.h5';
% FORMAT OF DATA
%1 - time of day (days, UT I think??)
%2 - pierce point lat
%3 - pierce point long
%4 - TID estimate (dTEC est?)
%5 - vTEC from GPS site (calc'd)
%6 - site (guess site #?)
% h5disp('tid_127_2013.h5');
tecDataChunk = 10000000; %max chunk size (prevent memory bleedout by matlab)
tecDataLimPercent = 0.05; %0.05 = 5%, cut out times with very low data content (less than 5% of the mean data content #)

dateRange = [2013,126;2013,128]; %date range 
dateRange_zeroHr = [2013,127]; %the day for 0 hour to occur - counting in hours for ease, but need reference

plot_Freq_Lim = 120; %min, time to limit the freq analysis plots to show
Xaxisvar_SCARGLE = 0:10:plot_Freq_Lim; %create a periodogram x axis tick

plotLatRange = [35,50]; %latitude limit for plotting
%-90 to 90 is world, 35 to 50 is good for USA East Coast
plotLongRange = [-85,-60]; %longitude limit for plotting
%-180 to 180 is world, -85 to -60 is good for USA East Coast
pointRadius = 50; %km, radius around the point to search for nearby points
plot_Scatter_Point_Size = 325; %arb. scatter pt size to make it big and easy to see

pointAltitude = 300; %km, altitude where most e-'s are (F region max - assumed)

time_cutout_range = [18,35]; %hrs, cut-out this time period for comparison between the ISR & GPS data at the best possible time to compare (night)

% plotLatRange = [-90,90]; %latitude limit for plotting
%-90 to 90 is world, 35 to 50 is good for USA East Coast
% plotLongRange = [-180,180]; %longitude limit for plotting
%-180 to 180 is world, -85 to -60 is good for USA East Coast

% plotLatRange = [30,50]; %latitude limit for USA West Coast
% plotLongRange = [-130,-110]; %longitude limit for USA West Coast

% plotLatRange = [30,50]; %latitude limit for plotting % %-90 to 90 is world, 35 to 50 is good for USA East Coast
% plotLongRange = [-125,-60]; %longitude limit for plotting usa

% plotLatRange = [30,75]; %latitude limit for USA & ALASKA & CANADA & GREENLAND KINDA
% plotLongRange = [-165,-30]; %longitude limit for USA & ALASKA & CANADA & GREENLAND KINDA

 plotLatRange = [0,75]; %latitude limit for North America & Europe
 plotLongRange = [-180,45]; %longitude limit for North America & Europe

% plotLatRange = [30,75]; %latitude limit for Europe
% plotLongRange = [-15,40]; %longitude limit for Europe

% plotLatRange = [30,46]; %good for Japan close up
% plotLongRange = [130,146]; %good for Japan close up

% plotLatRange = [20,55]; %good for Japan a bit farther away
% plotLongRange = [115,160]; %good for Japan a bit farther away

Re = 6371.0; %km, Earth mean Radius
%from http://nssdc.gsfc.nasa.gov/planetary/factsheet/earthfact.html

%Location of Millstone Hill
latMillstone = 42.6233; %Deg North
longMillstone = -71.4882; %Deg East (71.4882 deg West)
% latMillstone = 42.6233 - (100/Re)*(180/pi); %Deg North
% plotLatRange = [latMillstone-2,latMillstone+2]; %new latitude limit for plotting
% plotLongRange = [longMillstone-2,longMillstone+2]; %new longitude limit for plotting
% longMillstone = -74.0; %deg east, philly
% longMillstone = 23+10/60; %deg west, Kalix Sweeden (over south africa)
% longMillstone = 140+52/60; %deg west, Sendai Japan (over Eastern Japan)
% longMillstone = 119+44/60+49/3600; %deg west, Marble Bar Australia (over East China/West Australia)
% latMillstone = mean(plotLatRange); %center
% longMillstone = mean(plotLongRange); %center
% latMillstone = 38; %put at 38 deg N
% longMillstone = 37.5; %put at 54 deg E

MillstoneMISA_azimuth = 168.5; %deg Azimuth from North towards East (assuming geo - none specified)
MillstoneMISA_elev = 66.26; %deg elevation
MillstoneZenith_elev = 88; %deg elevation

latMillstoneMISA = latMillstone + sind(MillstoneMISA_azimuth+90)*(((pointAltitude/tand(MillstoneMISA_elev))/Re)*180/pi); %deg North (see spreadsheet titled ISR Angle Calc)
longMillstoneMISA = longMillstone - cosd(MillstoneMISA_azimuth+90)*(((pointAltitude/tand(MillstoneMISA_elev))/Re)*180/pi); %deg East at MISA (Azimuth is clockwise - cosine gets a negative)


%FOR AVG'D ALONG LONG/LAT PLOTS
%avg'd long plot changers
pplatN = 200; %number of bits to split pplat into (plus one below, for counting purposes)
pplongAvgDeg = 2; %deg, degree to average up and down around millstone (2low, 8.5 is maxish before it gets way avg'd) (NOT TOTAL WIDTH)
%avg'd lat plot changers
pplongN = 200; %number of bits to split pplat into (plus one below, for counting purposes)
pplatAvgDeg = 2; %deg, degree to average up and down around millstone (NOT TOTAL WIDTH)

%any angle
avg_anyAngle = 90; %deg, user defined angle
avg_anyAngle_Width = 360; %arcdeg, total width - not an angle in this instance
avg_anyAngle_N = 200; %number of chunks to split the range into 
avg_anyAngle_45vsLatLong = 0; %0 for longitude on xaxis on a 45 degree angle (or multiple of it), 1 for latitude on xaxis
FLG_memSavr = 1; %1 on, 0 off, clears stuff will probably break next code
avg_anyAngle_Zoom = 5; %+/-# arcdeg around the Millstone Hill beam of choice, zoom for the time-cut plot

%JEDI Settings
JEDI_highpass_Freq = 1/2; %1/hr, high-pass filter cutoff freqeuncy
JEDI_nOrder = 42; %window order

%AMPERE Settings
fileName_AMPERE = ['JHeat_data_out_20130506';'JHeat_data_out_20130507';'JHeat_data_out_20130508']; %file names of data
gif_AMPERE_jouleHeating_cap = [1,5]; %erg/(cm^2*sec), low and high limits for plotting joule heating - taken from example video (seems Joule Heating low end limit is 1 - nothing less than 1)
AMPERE_jouleHeating_pos = 3; %set the position of the joule heating data
AMPERE_data_num = 3; %number of data entries per line
%1 Pedersen Conductance / 2 Hall Conductance / 3 Joule Heat (ergs/(cm^2*sec)) / 4 Electric Potential / 5 Field-Algined Current

latSpacing = 1; %deg, spacing of latitude points
% latPoints = [latMillstone+latSpacing,latMillstone,latMillstone-latSpacing]; %Deg North
% longPoints = [longMillstone,longMillstone,longMillstone]; %Deg East
latPoints = [latMillstone,latMillstoneMISA]; %Deg North
longPoints = [longMillstone,longMillstoneMISA]; %Deg East
% latMillstone = latMillstoneMISA; %force calls to go to MISA locale
% longMillstone = longMillstoneMISA; %force calls to go to MISA locale


dataReject = 2; %ratio that is multiplied by variance and used to determine how much data to eject
dataRejectMax = 4*dataReject; %maximum multiplier to reach
dataRejectLimit = 25; %percentage, limit on how much data can be jetisoned
dataRejectOrig = dataReject; %record original value
dataRejectLimitOrig = dataRejectLimit; %record original value

figNum = 1; %starting figure number
%starts at 2 b/c 1 reserved for movie (antiquated)

font_Weight = 'bold'; %sets font weight for all plots
font_Size = 18; %sets font size for all plots

gif_Type = 5; %see below
%0 = moving data points (Fastest, confusing to actually see what up)
%1 = stationary data points
%2 = stationary data points + Zenith ISR overlay
%3 = stationary data points + MISA ISR overlay
%4 = stationary data points + AMPERE data on same plot with time average to AMPERE data (every 10 min)
%5 = stationary data points + AMPERE data on same plot (no time average for TEC)
gif_TimeAvg = 6; %min, how many minutes to average into one picture (only used with gif_Type chosen appropriately)
gif_MP4 = 1; %1 for MP4, 0 for gif
gif_Grid_Lat_Spaces = 150; %how many spaces to break the latitude range up into (only for stationary data points)
gif_Grid_Long_Spaces = 300; %how many spaces to break the longitude range up into (only for stationary data points)
%below good for smaller zones (Japan)
% gif_Grid_Lat_Spaces = 40; %how many spaces to break the latitude range up into (only for stationary data points)
% gif_Grid_Long_Spaces = 80; %how many spaces to break the longitude range up into (only for stationary data points)
%keep to only square, interpolation seems to have an issue with non-square
%(to get around, set extrapolation to 'linear' or 'nearest' from 'none'
%-it's the last value in scatteredInterp or whateve rit's called)

gif_plotProjection = 0; %1 turns on mapping projections (uses X/Y coords and mapping toolbox), 0 directly plots things (plots as lat/long, only Mercator as that's just flat plotting)
%Note: X and Y labels do not appear with mapping projections
gif_plotProjection_Mercator = 0; %1 turns on the Mercator projection (looks same as no projection mode - but with mapping coords (e.g. they're X/Y coords not lat/long), 
    %0 is whatever projection MATLAB chooses as best
gif_ContinentFill = 1; %1 turns on coloring land/water, 0 only does continent outlines
gif_ContientColor = [0.5 0.7 0.5]; %sets the land color (if gif_ContinentFill = 1)
gif_ContinentWaterColor = [176/255,196/255,222/255]; %sets the water color (if gif_ContinentFill = 1)
gif_Millstone_Marker_Color = [255/255,36/255,0/255]; %sets the color of the Millstone Hill ISR Zenith marker on maps
gif_Millstone_Marker = '*'; %sets the marker used to denote the Millstone Hill ISR Zenith marker on maps

gif_DesiredMaxRunTime = 120; %sec, desired max run time (to keep video limited) - if calc'd over this it will bump from 30 FPS to 60 FPS
gif_FLGdisableFPSShift = 1; %use this to disable the jump to 60 FPS if it makes it go too fast (0 allow 60 FPS jump, 1 disabled)
gif_DesiredFrameTime = 0.35; %sec, desired time between frame for MP4 (how long a frame is on the screen)
gif_Delay = 0.05; %s, time between gif frames

gif_sTEC_cap = 0.5; %TECU, limit (+ and -) for plotting the gif and other plots too

gif_Save_Locale = 'C:\Users\Ode to Hagis\Folder of Folders\College\Penn State\Research\Semester 5\New New Big Data\Plots'; %location of gif
% gif_Save_Locale = 'C:\Users\Razzle Dazzle\Folder of Folders\College\Penn State\Research\Semester 5\New New Big Data\Plots';

gif_Name = 'sTEC_movie'; %name of movie to make

gif_Scatter_Point_Size = 20; %arb. size to make points
gif_Millstone_Marker_Size = 20; %arb. size to make a marker (bigger than the scatter pt size btw)
gif_Scatter_Point_Size_AMPERE = 50; %arb. size to make points

gif_Grid_Div = 4^1; %number of blocks to divide a 1 deg by 1 deg square into
%gotta be a power of 4 at this point don't remember why

dateOfData = decyear(2013,5,7); %get date needed (manually, approx)

% warning('off','dataViewer:FLG') %use to turn off flag warnings

%% Error Catch
if( sum(size(latPoints)) ~= sum(size(longPoints)) ) %Error if latitude and longitude vectors do not match size
    error(['PRERUN ERROR CHECK: Lat/long point sizes do not match. Lat size: ',num2str(size(latPoints)),' Long size: ',num2str(size(longPoints))]);
end
if( exist(gif_Save_Locale,'dir') == 0 ) %Errror if gif save locale folder does not exist
   error(['PRERUN ERROR CHECK: Folder does not exist: ',gif_Save_Locale]); 
end
%Check for dependancies that don't exist

%ERROR FOR DEPENDANT DISABLED
% if( FLG_gatherDataAtPoint == 0 && FLG_gatherDataAtPoint_Filter == 1 )
%     error('dataViewer:FLG','FLG_gatherDataAtPoint set off (0) while FLG_gatherDataAtPoint_Filter set on (1) but depends on turned off functionality.');
% end  
% %ERROR FOR DEPENDANT DISABLED
% if( FLG_gatherDataAtPoint_Filter == 0 && FLG_gatherDataAtPoint_Filter_Plot == 1 )
%     error('dataViewer:FLG','FLG_gatherDataAtPoint_Filter set off (0) while FLG_gatherDataAtPoint_Filter_Plot set on (1) but depends on turned off functionality.');
% end
% %ERROR FOR DEPENDANT DISABLED
% if( FLG_gatherDataAtPoint_Filter == 0 && FLG_gatherDataAtPoint_Filter_Scargle == 1 )
%     error('dataViewer:FLG','FLG_gatherDataAtPoint_Filter set off (0) while FLG_gatherDataAtPoint_Filter_Scargle set on (1) but depends on turned off functionality.');
% end
% %ERROR FOR DEPENDANT DISABLED
% if( FLG_gatherDataAtPoint_Filter == 0 && FLG_gatherDataAtPoint_Filter_BP == 1 )
%     error('dataViewer:FLG','FLG_gatherDataAtPoint_Filter set off (0) while FLG_gatherDataAtPoint_Filter_BP set on (1) but depends on turned off functionality.');
% end
% %ERROR FOR DEPENDANT DISABLED
% if( FLG_ISRdata_Zenith == 0 && FLG_ISRdata_Zenith_Plot == 1 || FLG_gatherDataAtPoint_Filter == 0 && FLG_ISRdata_Zenith_Plot == 1 )
%     error('dataViewer:FLG','FLG_ISRdata_Zenith set off (0) OR FLG_gatherDataAtPoint_Filter set off (0) while FLG_ISRdata_Zenith_Plot set on (1) but depends on turned off functionality.');
% end
% %ERROR FOR DEPENDANT DISABLED
% if( FLG_ISRdata_Zenith == 0 && FLG_ISRdata_Zenith_Plot_BP == 1 || FLG_gatherDataAtPoint_Filter_BP == 0 && FLG_ISRdata_Zenith_Plot_BP == 1 )
%     error('dataViewer:FLG','FLG_ISRdata_Zenith set off (0) OR FLG_gatherDataAtPoint_Filter_BP set off (0) while FLG_ISRdata_Zenith_Plot_BP set on (1) but depends on turned off functionality.');
% end
% %ERROR FOR DEPENDANT DISABLED
% if( FLG_ISRdata_MISA == 0 && FLG_ISRdata_MISA_Plot == 1 || FLG_gatherDataAtPoint_Filter == 0 && FLG_ISRdata_MISA_Plot == 1 )
%     error('dataViewer:FLG','FLG_ISRdata_MISA set off (0) OR FLG_gatherDataAtPoint_Filter set off (0) while FLG_ISRdata_MISA_Plot set on (1) but depends on turned off functionality.');
% end
% %ERROR FOR DEPENDANT DISABLED
% if( FLG_ISRdata_MISA == 0 && FLG_ISRdata_MISA_Plot_BP == 1 || FLG_gatherDataAtPoint_Filter_BP == 0 && FLG_ISRdata_MISA_Plot_BP == 1 )
%     error('dataViewer:FLG','FLG_ISRdata_MISA set off (0) OR FLG_gatherDataAtPoint_Filter_BP set off (0) while FLG_ISRdata_MISA_Plot_BP set on (1) but depends on turned off functionality.');
% end
% %ERROR FOR DEPENDANT DISABLED
% if( FLG_ISRdata_Zenith == 0 && FLG_ISRdata_Zenith_Corr == 1 || FLG_gatherDataAtPoint_Filter_BP == 0 && FLG_ISRdata_Zenith_Corr == 1 )
%     error('dataViewer:FLG','FLG_ISRdata_Zenith set off (0) OR FLG_gatherDataAtPoint_Filter_BP set off (0) while FLG_ISRdata_Zenith_Corr set on (1) but depends on turned off functionality.');
% end
% %ERROR FOR DEPENDANT DISABLED
% if( FLG_ISRdata_MISA == 0 && FLG_ISRdata_MISA_Corr == 1 || FLG_gatherDataAtPoint_Filter == 0 && FLG_ISRdata_MISA_Corr == 1 )
%     error('dataViewer:FLG','FLG_ISRdata_MISA set off (0) OR FLG_gatherDataAtPoint_Filter set off (0) while FLG_ISRdata_MISA_Corr set on (1) but depends on turned off functionality.');
% end
% %ERROR FOR DEPENDANT DISABLED
% if( FLG_ISRdata_MISA == 0 && FLG_AVG_LONG == 1 )
%     error('dataViewer:FLG','FLG_ISRdata_MISA set off (0) while FLG_AVG_LONG set on (1) but depends on turned off functionality.');
% end
% if( FLG_ISRdata_MISA_Corr == 0 && FLG_AVG_LONG == 1 )
%     error('dataViewer:FLG','FLG_ISRdata_MISA_Corr set off (0) while FLG_AVG_LONG set on (1) but depends on turned off functionality.');
% end
% %ERROR FOR DEPENDANT DISABLED
% if( FLG_ISRdata_MISA == 0 && FLG_AVG_LAT == 1 )
%     error('dataViewer:FLG','FLG_ISRdata_MISA set off (0) while FLG_AVG_LAT set on (1) but depends on turned off functionality.');
% end
% %ERROR FOR DEPENDANT DISABLED
% if( FLG_ISRdata_MISA == 0 && FLG_AVG_ANYANGLE == 1 )
%     error('dataViewer:FLG','FLG_ISRdata_MISA set off (0) while FLG_AVG_ANYANGL set on (1) but depends on turned off functionality.');
% end
% %ERROR FOR DEPENDANT DISABLED
% if( FLG_AVG_ANYANGLE == 0 && FLG_AVG_ANYANGLE_PLOT == 1 )
%     error('dataViewer:FLG','FLG_AVG_ANYANGLE set off (0) while FLG_AVG_ANYANGLE_PLOT set on (1) but depends on turned off functionality.');
% end
% %ERROR FOR DEPENDANT DISABLED
% if( FLG_ISRdata_MISA_Corr == 0 && FLG_AVG_LONG_Scargle == 1 )
%     error('dataViewer:FLG','FLG_ISRdata_MISA_Corr set off (0) while FLG_AVG_LONG_Scargle set on (1) but depends on turned off functionality.');
% end
% %ERROR FOR DEPENDANT DISABLED
% if( FLG_ISRdata_Zenith == 0 && FLG_ISRdata_Zenith_Scargle == 1 )
%     error('dataViewer:FLG','FLG_ISRdata_Zenith set off (0) while FLG_ISRdata_Zenith_Scargle set on (1) but depends on turned off functionality.');
% end
% %ERROR FOR DEPENDANT DISABLED
% if( FLG_ISRdata_MISA == 0 && FLG_ISRdata_MISA_Scargle == 1 )
%     error('dataViewer:FLG','FLG_ISRdata_MISA set off (0) while FLG_ISRdata_MISA_Scargle set on (1) but depends on turned off functionality.');
% end
% %ERROR FOR DEPENDANT DISABLED
% if( FLG_gatherDataAtPoint_Filter_BP == 0 && FLG_gatherDataAtPoint_Filter_BP_Scargle == 1 )
%     error('dataViewer:FLG','FLG_gatherDataAtPoint_Filter_BP set off (0) while FLG_gatherDataAtPoint_Filter_BP_Scargle set on (1) but depends on turned off functionality.');
% end
% %ERROR FOR DEPENDANT DISABLED
% if( FLG_gatherDataAtPoint_Filter == 0 && FLG_gatherDataAtPoint_Filter_TimeMatch_n_Smooth == 1 )
%     error('dataViewer:FLG','FLG_gatherDataAtPoint_Filter set off (0) while FLG_gatherDataAtPoint_Filter_TimeMatch_n_Smooth set on (1) but depends on turned off functionality.');
% end
% if( FLG_ISRdata_Zenith_Corr == 0 && FLG_gatherDataAtPoint_Filter_TimeMatch_n_Smooth == 1 )
%     error('dataViewer:FLG','FLG_ISRdata_Zenith_Corr set off (0) while FLG_gatherDataAtPoint_Filter_TimeMatch_n_Smooth set on (1) but depends on turned off functionality.');
% end 
% if( FLG_ISRdata_MISA_Corr == 0 && FLG_gatherDataAtPoint_Filter_TimeMatch_n_Smooth == 1 )
%     error('dataViewer:FLG','FLG_ISRdata_MISA_Corr set off (0) while FLG_gatherDataAtPoint_Filter_TimeMatch_n_Smooth set on (1) but depends on turned off functionality.');
% end 
% %ERROR FOR DEPENDANT DISABLED
% if( FLG_gatherDataAtPoint_Filter_TimeMatch_n_Smooth== 0 && FLG_gatherDataAtPoint_Filter_TimeMatch_n_Smooth_ISRCorr == 1 )
%     error('dataViewer:FLG','FLG_gatherDataAtPoint_Filter_TimeMatch_n_Smooth set off (0) while FLG_gatherDataAtPoint_Filter_TimeMatch_n_Smooth_ISRCorr set on (1) but depends on turned off functionality.');
% end
% %ERROR FOR DEPENDANT DISABLED
% if( FLG_AVG_LONG == 0 && FLG_AVG_LONG_TimeMatch_n_BP_Scargle == 1 )
%     error('dataViewer:FLG','FLG_AVG_LONG set off (0) while FLG_AVG_LONG_TimeMatch_n_BP_Scargle set on (1) but depends on turned off functionality.');
% end
% %ERROR FOR DEPENDANT DISABLED
% if( FLG_AVG_LONG == 0 && FLG_AVG_LONG_TimeMatch_n_BP_Scargle == 1 )
%     error('dataViewer:FLG','FLG_AVG_LONG set off (0) while FLG_AVG_LONG_TimeMatch_n_BP_Scargle set on (1) but depends on turned off functionality.');
% end
% %ERROR FOR DEPENDANT DISABLED
% if( FLG_gatherDataAtPoint_Filter_TimeMatch_n_Smooth == 0 && FLG_gatherDataAtPoint_Filter_TimeMatch_n_Smooth == 1 )
%     error('dataViewer:FLG','FLG_gatherDataAtPoint_Filter_TimeMatch_n_Smooth set off (0) while FLG_gatherDataAtPoint_Filter_TimeMatch_n_Smooth set on (1) but depends on turned off functionality.');
% end
% %ERROR FOR DEPENDANT DISABLED
% if( FLG_gatherDataAtPoint_Filter_TimeMatch_n_Smooth == 0 && FLG_gatherDataAtPoint_Filter_TimeMatch_n_Smooth_BP == 1 )
%     error('dataViewer:FLG','FLG_gatherDataAtPoint_Filter_TimeMatch_n_Smooth set off (0) while FLG_gatherDataAtPoint_Filter_TimeMatch_n_Smooth_BP set on (1) but depends on turned off functionality.');
% end
% %ERROR FOR DEPENDANT DISABLED
% if( FLG_gatherDataAtPoint_Filter_TimeMatch_n_Smooth_BP == 0 && FLG_gatherDataAtPoint_Filter_TimeMatch_n_Smooth_BP_PlotsScargle == 1 )
%     error('dataViewer:FLG','FLG_gatherDataAtPoint_Filter_TimeMatch_n_Smooth_BP set off (0) while FLG_gatherDataAtPoint_Filter_TimeMatch_n_Smooth_BP_PlotsScargle set on (1) but depends on turned off functionality.');
% end
% %ERROR FOR DEPENDANT DISABLED
% if( FLG_gatherDataAtPoint_Filter_TimeMatch_n_Smooth == 0 && FLG_gatherDataAtPoint_Filter_ISRTimeMatch_Scargle_Walking == 1 )
%     error('dataViewer:FLG','FLG_gatherDataAtPoint_Filter_TimeMatch_n_Smooth set off (0) while FLG_gatherDataAtPoint_Filter_ISRTimeMatch_Scargle_Walking set on (1) but depends on turned off functionality.');
% end

%forcing of super memory saving fix
if( ( FLG_memSavr == 1 ) && min(plotLatRange) == -90 && max(plotLatRange) == 90 && min(plotLongRange) == -180 && max(plotLongRange) == 180 )
    FLG_memSavr = 2; %force upgrade to 2 which doesn't do any data range culling since the whole world is chosen
end

%if mag conversion on, prep the system to call the IGRF conversion library
if( FLG_geomagneticCoords == 1 )
    magF = IGRF(); %prep using the external conversion library
    
end

%if gif_Type == 4 but FLG_AMPERE_DATA == 0 (gif_Type == 4 plots AMPERE data too)
if( gif_Type == 4 && FLG_AMPERE_DATA == 0 )
    warning(['PRERUN ERROR CHECK: gif_Type == 3 to plot TEC + AMPERE data, but FLG_AMPERE_DATA == 0. Forcing FLG_AMPERE_DATA == 1']); 
    FLG_AMPERE_DATA = 1; %force it on so we can do it
end

%plot help with autotick calculating
plotLongRange_autoTick = (max(plotLongRange) - min(plotLongRange))/25; %tries to split the longitude range into 25 parts (based off of 360/15+1)
if( plotLongRange_autoTick > 10 )
    plotLongRange_autoTick = 15; %sets the tick setting to 15 arcdegrees per tick
elseif( plotLongRange_autoTick > 5 )
    plotLongRange_autoTick = 10; %sets the tick setting to 10 arcdegrees per tick
elseif( plotLongRange_autoTick > 2 )
    plotLongRange_autoTick = 5; %sets the tick setting to 5 arcdegrees per tick
elseif( plotLongRange_autoTick > 1 )
    plotLongRange_autoTick = 2; %sets the tick setting to 5 arcdegrees per tick
elseif( plotLongRange_autoTick >= 0.6 ) %0.6 because 15/25 = 0.6, so there will be enough 1 arcdeg ticks
    plotLongRange_autoTick = 1; %sets the tick setting to 1 arcdegree per tick
else
    plotLongRange_autoTick = (max(plotLongRange) - min(plotLongRange))/15; %just goes for it if it's a super tiny range
end
plotLatRange_autoTick = (max(plotLatRange) - min(plotLatRange))/13; %tries to split the latitude range into 13 parts (based off of 180/15+1)
if( plotLatRange_autoTick > 10 )
    plotLatRange_autoTick = 15; %sets the tick setting to 15 arcdegrees per tick
elseif( plotLatRange_autoTick > 5 )
    plotLatRange_autoTick = 10; %sets the tick setting to 10 arcdegrees per tick
elseif( plotLatRange_autoTick > 2 )
    plotLatRange_autoTick = 5; %sets the tick setting to 5 arcdegrees per tick
elseif( plotLatRange_autoTick > 1 )
    plotLatRange_autoTick = 2; %sets the tick setting to 2 arcdegrees per tick
elseif( plotLatRange_autoTick > 0.75 ) %0.75 because 10/13 = 0.76something and it sounded good for enough 1 arcdeg ticks
    plotLatRange_autoTick = 1; %sets the tick setting to 1 arcdegree per tick
else
    plotLatRange_autoTick = (max(plotLatRange) - min(plotLatRange))/15; %just goes for it if it's a super tiny range
end
plotLatRange_autoTick_Crunched = (max(plotLatRange) - min(plotLatRange))/7; %tries to split the latitude range into 7 parts (based off of seems only 7 fits on a 2 subplot plot)
if( plotLatRange_autoTick_Crunched > 10 )
    plotLatRange_autoTick_Crunched = 15; %sets the tick setting to 15 arcdegrees per tick
elseif( plotLatRange_autoTick_Crunched > 5 )
    plotLatRange_autoTick_Crunched = 10; %sets the tick setting to 10 arcdegrees per tick
elseif( plotLatRange_autoTick_Crunched > 2 )
    plotLatRange_autoTick_Crunched = 5; %sets the tick setting to 5 arcdegrees per tick
elseif( plotLatRange_autoTick_Crunched > 1 )
    plotLatRange_autoTick_Crunched = 2; %sets the tick setting to 2 arcdegrees per tick
elseif( plotLatRange_autoTick_Crunched > 0.75 ) %0.75 because 10/13 = 0.76something and it sounded good for enough 1 arcdeg ticks
    plotLatRange_autoTick_Crunched = 1; %sets the tick setting to 1 arcdegree per tick
else
    plotLatRange_autoTick_Crunched = (max(plotLatRange) - min(plotLatRange))/15; %just goes for it if it's a super tiny range
end

%AVG any angle error fix - it can't do exactly 0 or exactly 90 degrees
if( avg_anyAngle == 0 || avg_anyAngle == 270 )
    avg_anyAngle = avg_anyAngle + 0.0001; %deg, adjust so not exactly 0 or 270 (to be symmetrically consistent)
elseif( avg_anyAngle == 90 || 180 )
    avg_anyAngle = avg_anyAngle - 0.0001; %deg, adjust so not exactly 90 or 180 (to be symmetrically consistent)
end
%AVG any angle also can't do exactly the full width (so a width that takes the whole plot area)
if( (round(avg_anyAngle) == 0 || round(avg_anyAngle) == 180) && (avg_anyAngle_Width >= abs(plotLatRange(1) - plotLatRange(2))) )
    if( avg_anyAngle_Width > abs(plotLatRange(1) - plotLatRange(2)) )
%         error('FUN_dataViewer:syntaxCheck',['Avg any angle width of ',num2str(avg_anyAngle_Width),' arcdeg is larger than the latitudinal plot area of ',num2str(abs(plotLatRange(1) - plotLatRange(2))),' arcdeg.']);
        warning('FUN_dataViewer:syntaxCheck',['Avg any angle width of ',num2str(avg_anyAngle_Width),' arcdeg is larger than the latitudinal plot area of ',num2str(abs(plotLatRange(1) - plotLatRange(2))),' arcdeg. Reducing to be size of latitudinal plot area.']);
        avg_anyAngle_Width = abs(plotLatRange(1) - plotLatRange(2)) - 0.001; %arcdeg, adjust so not exactly plot width, but basically since overstepped the range already
    else
        avg_anyAngle_Width = avg_anyAngle_Width - 0.001; %arcdeg, adjust so not exactly plot width
    end
elseif( (round(avg_anyAngle) == 90 || round(avg_anyAngle) == 270) && (avg_anyAngle_Width >= abs(plotLongRange(1) - plotLongRange(2))) )
    if( avg_anyAngle_Width > abs(plotLongRange(1) - plotLongRange(2)) )
%         error('FUN_dataViewer:syntaxCheck',['Avg any angle width of ',num2str(avg_anyAngle_Width),' arcdeg is larger than the longitudinal plot area of ',num2str(abs(plotLongRange(1) - plotLongRange(2))),' arcdeg.']);
        warning('FUN_dataViewer:syntaxCheck',['Avg any angle width of ',num2str(avg_anyAngle_Width),' arcdeg is larger than the longitudinal plot area of ',num2str(abs(plotLongRange(1) - plotLongRange(2))),' arcdeg. Reducing to be size of longitudinal plot area.']);
        avg_anyAngle_Width = abs(plotLongRange(1) - plotLongRange(2)) - 0.001; %arcdeg, adjust so not exactly plot width, but just under it since overstepped the range already
    else
        avg_anyAngle_Width = avg_anyAngle_Width - 0.001; %arcdeg, adjust so not exactly plot width
    end
end


%% Read the Data File
%OLDER VERSIONS OF MATLAB MUST SPLIT IT UP (SLOWER) - left for support
% for( fileCntr = 1:size(fileName,1) )
% 
%     fileInfo = h5info(fileName(fileCntr,:),'/tecs'); %get the info on the data
%     tecDataSize = fileInfo.Dataspace.MaxSize(2); %get how many data points there are
%     
%     if(fileCntr == 1) %this will init our variables
%         for(i = 1:ceil(tecDataSize/tecDataChunk) ) %loop to prevent RAM issues (read it in chunks, not all at once)
% 
%             if(i == ceil(tecDataSize/tecDataChunk))
%                 %means last i, which will not be a nice ending number
%                 tecDataTemp = double(h5read(fileName(fileCntr,:),'/tecs',[1,1+(i-1)*tecDataChunk],[6,tecDataSize-(1+(i-1)*tecDataChunk)+1]));
%                 %read an amount of data (to the end of file)
%             else
%                 tecDataTemp = double(h5read(fileName(fileCntr,:),'/tecs',[1,1+(i-1)*tecDataChunk],[6,tecDataChunk]));
%                 %read an amount of data
%             end
%             
%             if(FLG_geomagneticCoords == 0 ) %if mag coords, del later
%                 %remove data that is out of desired range (save that RAM)
%                 k = (tecDataTemp(3,:) > max(plotLongRange)) | (tecDataTemp(3,:) < min(plotLongRange)) ...
%                         | (tecDataTemp(2,:) > max(plotLatRange)) | (tecDataTemp(2,:) < min(plotLatRange));
%                 tecDataTemp(:,k) = []; %del it
%             end
%             
% %             if( FLG_geomagneticCoords == 1 ) %convert data to geomagnetic coords if desired
% %                 [tecDataTemp(2,:),tecDataTemp(3,:)] = sFUN_geoToGeomag(dateOfData,magF,tecDataTemp(2,:),tecDataTemp(3,:)); %convert from geographic to geomagnetic coordinates
% %             end %doing as we read data bit by bit because full math conversion will cause mem to run out with my 16 gigs omg
% 
%             if( i ~= 1 ) %tack on data (or create data variable)
%                 time = [time,tecDataTemp(1,:)]; %create desired data vectors
%                 pplat = [pplat,tecDataTemp(2,:)]; %create desired data vectors
%                 pplong = [pplong,tecDataTemp(3,:)]; %create desired data vectors
%                 sTEC = [sTEC,tecDataTemp(4,:)]; %create desired data vectors
%             else
%                 time = tecDataTemp(1,:); %create desired data vectors
%                 pplat = tecDataTemp(2,:); %create desired data vectors
%                 pplong = tecDataTemp(3,:); %create desired data vectors
%                 sTEC = tecDataTemp(4,:); %create desired data vectors
%             end
%         end
% %         clear tecDataTemp
% 
%     else %this will tack on the data
%         for(i = 1:ceil(tecDataSize/tecDataChunk) ) %loop to prevent RAM issues (read it in chunks, not all at once)
% 
%             if(i == ceil(tecDataSize/tecDataChunk))
%                 %means last i, which will not be a nice ending number
%                 tecDataTemp = double(h5read(fileName(fileCntr,:),'/tecs',[1,1+(i-1)*tecDataChunk],[6,tecDataSize-(1+(i-1)*tecDataChunk)+1]));
%                 %read an amount of data (to the end of file)
%             else
%                 tecDataTemp = double(h5read(fileName(fileCntr,:),'/tecs',[1,1+(i-1)*tecDataChunk],[6,tecDataChunk]));
%                 %read an amount of data
%             end
%             
%             if( FLG_geomagneticCoords == 0 ) %if mag coords, cut later
%                 %remove data that is out of desired range (save that RAM)
%                 k = (tecDataTemp(3,:) > max(plotLongRange)) | (tecDataTemp(3,:) < min(plotLongRange)) ...
%                         | (tecDataTemp(2,:) > max(plotLatRange)) | (tecDataTemp(2,:) < min(plotLatRange));
%                 tecDataTemp(:,k) = []; %del it
%             end
%             
% %             if( FLG_geomagneticCoords == 1 ) %convert data to geomagnetic coords if desired
% %                 tic
% %                 [tecDataTemp(2,:),tecDataTemp(3,:)] = sFUN_geoToGeomag(dateOfData,magF,tecDataTemp(2,:),tecDataTemp(3,:)); %convert from geographic to geomagnetic coordinates
% %                 toc
% %             end %doing as we read data bit by bit because full math conversion will cause mem to run out with my 16 gigs omg
% 
%             time = [time,tecDataTemp(1,:)]; %create desired data vectors
%             pplat = [pplat,tecDataTemp(2,:)]; %create desired data vectors
%             pplong = [pplong,tecDataTemp(3,:)]; %create desired data vectors
%             sTEC = [sTEC,tecDataTemp(4,:)]; %create desired data vectors
%         end
% %         clear tecDataTemp
% 
%     end
% 
%     fprintf(['File Name: %s\nOrig. Data Size: 6 x %d\nKept Data Size: 6 x %d\n',...
%         '%% Of Data Kept: %%%.4f\nLat @ Millstone: %.3f deg\tLong @ Millstone: %.3f deg\tRadius Around Pt: %.2f km\n\n'],...
%         fileName(fileCntr,:),tecDataSize,size(sTEC,2),size(sTEC,2)/tecDataSize*100,latMillstone,longMillstone,pointRadius);
% end
% %do some mem saving moves
% clear tecDataTemp  %clear unused var
% time = single(time); %convert time to a single floating point precision to save mem
% sTEC = single(sTEC); %math doesn't need precision, if it does on specific case can re-enable double for calc, but these numbers have few decimal places
% pplat = single(pplat);
% pplong = single(pplong);
% %done doing mem saving moves
%*************STILL FILE READING CODE**************
%NEWER (not sure when) VERSIONS OF MATLAB CAN DO WHOLE FILE:
timeHoldr = cell(size(fileName,1),1); %preallocate
pplatHoldr = cell(size(fileName,1),1); %preallocate
pplongHoldr = cell(size(fileName,1),1); %preallocate
sTECHoldr = cell(size(fileName,1),1); %preallocate

fprintf('Beginning Data Import\n\n');
tic
for( fileCntr = 1:size(fileName,1) )

    fileInfo = h5info(fileName(fileCntr,:),'/tecs'); %get the info on the data
    tecDataSize = fileInfo.Dataspace.MaxSize(2); %get how many data points there are
    
    tecDataTemp = single(h5read(fileName(fileCntr,:),'/tecs'));
    %read an amount of data (to the end of file)
    
    if(FLG_geomagneticCoords == 0 ) %if mag coords, del later
        %remove data that is out of desired range (save that RAM)
        k = (tecDataTemp(3,:) > max(plotLongRange)) | (tecDataTemp(3,:) < min(plotLongRange)) ...
                | (tecDataTemp(2,:) > max(plotLatRange)) | (tecDataTemp(2,:) < min(plotLatRange));
        tecDataTemp(:,k) = []; %del it
    end
    
    timeHoldr{fileCntr} = tecDataTemp(1,:); %pull out time stuff
    pplatHoldr{fileCntr} = tecDataTemp(2,:); %pull out time stuff
    pplongHoldr{fileCntr} = tecDataTemp(3,:); %pull out time stuff
    sTECHoldr{fileCntr} = tecDataTemp(4,:); %pull out time stuff
    
    fprintf(['File Name: %s\nOrig. Data Size: 6 x %d\nKept Data Size: 4 x %d\n',...
        '%% Of Data Kept: %%%.4f\n\n'],...
        fileName(fileCntr,:),tecDataSize,length(tecDataTemp(1,:)),length(tecDataTemp(1,:))/tecDataSize*100);
end
clear tecDataTemp k fileInfo tecDataSize  %clear unused var

%Unpack the holder variables
for(i = 1:size(fileName,1))
    if ( i == 1 )
        time = timeHoldr{i}; %save the data
        timeHoldr(i) = {double.empty(0)}; %clear it out to save memory, but don't delete the spot for counting purposes
        pplat = pplatHoldr{i}; 
        pplatHoldr(i) = {double.empty(0)}; %clear it out to save memory, but don't delete the spot for counting purposes
        pplong = pplongHoldr{i};
        pplongHoldr(i) = {double.empty(0)}; %clear it out to save memory, but don't delete the spot for counting purposes
        sTEC = sTECHoldr{i};
        sTECHoldr(i) = {double.empty(0)}; %clear it out to save memory, but don't delete the spot for counting purposes
    else
        time = [time,timeHoldr{i}]; %tack on
        timeHoldr(i) = {double.empty(0)}; %clear it out to save memory, but don't delete the spot for counting purposes
        pplat = [pplat,pplatHoldr{i}]; %tack on
        pplatHoldr(i) = {double.empty(0)}; %clear it out to save memory, but don't delete the spot for counting purposes
        pplong = [pplong,pplongHoldr{i}]; %tack on
        pplongHoldr(i) = {double.empty(0)}; %clear it out to save memory, but don't delete the spot for counting purposes
        sTEC = [sTEC,sTECHoldr{i}]; %tack on
        sTECHoldr(i) = {double.empty(0)}; %clear it out to save memory, but don't delete the spot for counting purposes
    end
end
tocTemp = toc; %save toc
fprintf('Data Import took %.2f sec / %.2f min\nData size: %d\n\n',tocTemp,tocTemp/60,length(time)); %report
clear timeHoldr pplatHoldr pplongHoldr sTECHoldr tocTemp %clear unused var

time = time + 1; %days, adjust because the data comes starting at 0 days not 1 days

timeUnique = unique(time); %days, gathers unique times

if( FLG_geomagneticCoords == 1 ) %convert data to geomagnetic coords if desired
%     [pplat,pplong] = sFUN_geoToGeomag(dateOfData,magF,pplat,pplong); %convert from geographic to geomagnetic coordinates
    load magcoord_pplat_126127128 %faster to just load these up that are pre calc'd
    load magcoord_pplong_126127128 %faster to just load these up that are pre calc'd
    k = (pplong > max(plotLongRange)) | (pplong < min(plotLongRange)) ...
                    | (pplat > max(plotLatRange)) | (pplat < min(plotLatRange));
    pplat(k) = []; %delete excess
    pplong(k) = []; %delete excess
    time(k) = []; %delete also to match
    sTEC(k) = []; %delete also to match
end %doing as we read data bit by bit because full math conversion will cause mem to run out with my 16 gigs omg

% find( min(abs(timeUnique - (126+1/24))) == abs(timeUnique - (126+1/24))) %hold over to investigate every hour major sharp data spike

%Cut off time stamps with very little data in the range selected
tecDataAvgNum = length(sTEC)/length(timeUnique); %average data per time
tecDataLim = round(tecDataLimPercent*tecDataAvgNum); %min number before eject time set
[k,~] = histcounts(time,timeUnique); %get the number of occurance for each set
k = k < tecDataLim; %get the time uniques to delete for lack of data
k = ismember(time,timeUnique(k)); %find where time matches times to delete
%now, it does NOT check the very last timeUnique slot, so manual is go
k2 = time == timeUnique(end);
if( sum(k) < tecDataLim ) %if number of TEC points in last slot is less than required, delete last slot too
    k = k | k2; %tack on k2 if the data there needs to be deleted
end
clear k2
time(k) = []; %delete the low data stuff
pplat(k) = []; %delete
pplong(k) = []; %delete
sTEC(k) = []; %delete
%now, it does NOT check the very last timeUnique slot, so manual is go
k = time == timeUnique(end);
if( sum(k) < tecDataLim ) %if number of TEC points in last slot is less than required, delete last slot too
    time(k) = []; %delete the low data stuff
    pplat(k) = []; %delete
    pplong(k) = []; %delete
    sTEC(k) = []; %delete
end

timeUnique = unique(time); %days, gathers unique times again after cuttin

%GOAL HERE IS TO INVESTIGATE DATA NUMBER SEEN PER PERIOD AND MAKE SURE IT DOES NOT AFFECT DATA
if( FLG_dataNumber == 1 )
    
    %THIS IS TO VIEW DATA AVG RANGE BEING TAKEN
    coastLines = load('coast');%loads coast line into structure, vars are lat and long - can interfere with previous code so hidden in structure
    coastLines_lat = coastLines.lat; %breaks this data out
    coastLines_long = coastLines.long;
    clear coastLines
    if( FLG_geomagneticCoords == 1 ) %convert data to geomagnetic coords if desired
        [coastLines_lat,coastLines_long] = sFUN_geoToGeomag(dateOfData,magF,coastLines_lat,coastLines_long); %convert from geographic to geomagnetic coordinates
    end %doing as we read data bit by bit because full math conversion will cause mem to run out with my 16 gigs omg
    jk = find(coastLines_lat < min(plotLatRange) ); %find less than min, remove
    coastLines_lat(jk) = []; %remove
    coastLines_long(jk) = []; %remove
    jk = find(coastLines_lat > max(plotLatRange) ); %find more than max, remove
    coastLines_lat(jk) = []; %remove
    coastLines_long(jk) = []; %remove
    jk = find(coastLines_long < min(plotLongRange) ); %find less than min, remove
    coastLines_lat(jk) = []; %remove
    coastLines_long(jk) = []; %remove
    jk = find(coastLines_long > max(plotLongRange) ); %find more than max, remove
    coastLines_lat(jk) = []; %remove
    coastLines_long(jk) = []; %remove
    coastLines_long( coastLines_long == 180 ) = 180-.01; %fix some weird thing where 180 flips to -180

    %Corral the data to the right place    
    k = find( time == timeUnique(1)); %gets during a time period
    %prep figure #
    % figure('units','normalized','outerposition',[0 0 1 1]); %maximizes figure window
    figure(figNum);
    figNum = figNum+1;
    scatter(pplong(k),pplat(k),plot_Scatter_Point_Size,sTEC(k),'.');
    caxis([-gif_sTEC_cap, gif_sTEC_cap]);
    %         colormap(flipud(colormap('hot'))) %heat map where 0 (low) is white
    colormap('jet') %heat map
    h = colorbar; %shows the color bar
    ylabel(h, 'delta-sTEC (TECU)','fontweight',font_Weight, 'FontSize',font_Size); %labels color bar
    hold on;
    %Now drawing line of interest
    plot(longMillstone,latMillstone,'*','Color',[255/255,36/255,0/255],'MarkerSize',25,'LineWidth',1.9); %plots a point with a red big *
    hold on;
    %Plot lines of contients
    geoshow(coastLines_lat,coastLines_long,'Color','k') %add continental outlines in black  
    xlabel('Longitude (arcdeg)'); %old: ,'fontweight','bold','FontSize',12 for all
    ylabel('Latitude (arcdeg)'); %sped up significantly when merged above
    string_Title = ['delta-sTEC at Day ',num2str(round(timeUnique(1),2)),'  - Lat [',num2str(min(plotLatRange)),'x',num2str(max(plotLatRange)),'] Long [',num2str(min(plotLongRange)),'x',num2str(max(plotLongRange)),']']; %create mecha title
    if( FLG_geomagneticCoords == 1)
        string_Title = [string_Title,' [Mag Coord Aligned]']; %tack on mag coordinate warning
    end
    title(string_Title);

    xticks( min(plotLongRange):plotLongRange_autoTick:max(plotLongRange) ); %creates x ticks automagically
    yticks( min(plotLatRange):plotLatRange_autoTick:max(plotLatRange) ); %creates y ticks automagically
    set(gca,'fontweight',font_Weight, 'FontSize',font_Size);
    
    [k,~] = histcounts(time,timeUnique); %get the number of occurance for each set

    %prep figure #
%     figure('units','normalized','outerposition',[0 0 1 1]); %maximizes figure window (caused graphical issues on Windows 10 + AMD GPU or Intel GPU)
    figure(figNum);
    figNum = figNum+1;
    subplot(2,1,1); %subplot time!
    plot(timeUnique(1:end-2),k(1:end-1));
    xlabel('Time (days) from Day 126 to 128, 2013 (May 6 to May 8)');
    ylabel('Number of Data Points');
    title(['Number of Data Points vs Time - Lat [',num2str(min(plotLatRange)),'x',num2str(max(plotLatRange)),'] Long [',num2str(min(plotLongRange)),'x',num2str(max(plotLongRange)),']']);
    set(gca,'fontweight',font_Weight, 'FontSize',font_Size);
%     curtick = get(gca, 'YTick');
%     set(gca, 'YTickLabel', cellstr(num2str(curtick(:))));

    %SCARLGE OPTION
    Ross_Lombscargle_optimized([(timeUnique(1:end-1)'-floor(mean(time))).*24 ,k']); 
    fid = fopen('power_freq.txt', 'r');
        f_nm = textscan(fid, '%f %f %f %f %f', 'headerlines', 15);
        P_Zenith_SNR=f_nm{:,2};
        F_Zenith_SNR=f_nm{:,1};
        T_Zenith_SNR=f_nm{:,5}*60; % in minutes
        gs_Zenith_SNR=f_nm{:,4};
        max_Zenith_SNR=T_Zenith_SNR(P_Zenith_SNR==max(P_Zenith_SNR));
        K=length(P_Zenith_SNR);
        gf_Zenith_SNR=f_nm{:,3};
    fclose('all');

    %prep figure #
    %     figure('units','normalized','outerposition',[0 0 1 1]); %maximizes figure window
%     figure(figNum);
%     figNum = figNum+1;
    subplot(2,1,2); %subplot time!
    plot(T_Zenith_SNR,P_Zenith_SNR,'b');
    xlim([0 plot_Freq_Lim]);
    hold on;
    % plot(T_SNR_08apr08_10_or,P_SNR_08apr08_10_or,'r');
    % hold on;
    plot(T_Zenith_SNR,gf_Zenith_SNR,'LineWidth',1.5,'color',[0.5 0.5 0.5]);
    xlim([0 plot_Freq_Lim]);  
    xlabel('Periods (min)');
    ylabel('Normalized Power');
    title(['Periodogram of Data Number Sighted from Sats on Day 126 to Day 128 of 2013']);
    set(gca,'fontweight',font_Weight, 'FontSize',font_Size);
    set(gca,'xtick',Xaxisvar_SCARGLE);

    %FFT OPTION
%     totalcount = size(k',1); % # of samples taken
%     timeDelta = mean(timeUnique(2:end)-timeUnique(1:end-1))*24; %hrs already
%     deltatime = timeDelta; %sample interval in time - units hours *)
%     deltafreq = 1/(totalcount*deltatime); %frequency resolution per above sampling interval in time
%     frequencysamples = ((0:(totalcount-1))*deltafreq)';
%     % deltafreq = 1/((2^nextpow2(totalcount)*FFT_coeff)*deltatime);
%     % frequencysamples = ((0:2^nextpow2(totalcount)*FFT_coeff-1)*deltafreq)';
%     %frequencysamples = linspace(0,ceil((totalcount-1)*deltafreq),totalcount)';
%     % spectrummax = frequencysamples(totalcount/2 - 1);
%     % timedomain = zeros(totalcount,1);
% 
%     frequencydomain = fft(k); %,2^nextpow2(totalcount)*FFT_coeff
%     powerspectrum = sqrt(abs(frequencydomain.*conj(frequencydomain)))./totalcount; %properly normalized now
% 
%     gf=(1-(.05/length(powerspectrum))^(1/(length(powerspectrum)-1))); %for gf and gs see pg 934-935 of the journal paper "SPECTRUM: SPECTRAL ANALYSIS OF UNEVENLY SPACED PALEOCLIMATIC TIME SERIES" by SCHULZ and STATTEGGER
% %         gs=0.4*gf; %Computers & Geosciences Vol. 23, No. 9, pp. 929-945, 1997 
% 
%     subplot(2,1,2);
%     plot(1./frequencysamples*60,powerspectrum)
%     hold on;
%     plot(1./frequencysamples*60,repmat(gf,size(powerspectrum)),'LineWidth',1.5,'color',[0.5 0.5 0.5]);
%     xlabel('Period (min)');
%     ylabel('Normalized Power');
%     xlim([0 plot_Freq_Lim]);  
%     set(gca,'xtick',Xaxisvar_SCARGLE);
%     title(['DFT of Data Number Sighted from Sats on Day 126 to Day 128 of 2013']);
%     set(gca,'fontweight',font_Weight, 'FontSize',font_Size);
    % axis 'auto y'
    
    
    %prep figure #
%     figure('units','normalized','outerposition',[0 0 1 1]); %maximizes figure window (caused graphical issues on Windows 10 + AMD GPU or Intel GPU)
    figure(figNum);
    figNum = figNum+1;
    subplot(1,2,1)
    hWorld = plot(Re*cos(0:pi/1000:2*pi), Re*sin(0:pi/1000:2*pi) - Re,'LineWidth',1.5,'Color',[35/255,142/255,35/255]);
    hold on
    hSky = plot((Re+pointAltitude)*cos(0:pi/1000:2*pi), (Re+pointAltitude)*sin(0:pi/1000:2*pi) - Re,'LineStyle','--','LineWidth',1.5,'Color',[159/255,182/255,205/255]);
    hold on
    hAngle = plot(35*cosd(0:0.5:MillstoneMISA_elev), 35*sind(0:0.5:MillstoneMISA_elev),'LineWidth',1.5,'Color','c');
    hold on
    text(35*cosd(MillstoneMISA_elev/2)+2,35*sind(MillstoneMISA_elev/2)+2,[num2str(MillstoneMISA_elev),'^{\circ} Elevation'],'Color','k','fontweight','bold','FontSize',18); %plots text that shows angle
    hZenith = plot( zeros(1,length(0:1:pointAltitude)) , (0:1:pointAltitude),'LineWidth',1.5,'Color','k');
    hold on
    hMISA = plot( (0:1:pointAltitude).*sind(90-MillstoneMISA_elev)./sind(MillstoneMISA_elev) , (0:1:pointAltitude),'LineWidth',1.5,'Color','b');
    hold on
    text(pointAltitude.*sind(90-MillstoneMISA_elev)./sind(MillstoneMISA_elev)/2-25,pointAltitude-10,[num2str(round(pointAltitude.*sind(90-MillstoneMISA_elev)./sind(MillstoneMISA_elev))),' km Separation'],...
        'Color','k','fontweight','bold','FontSize',18); %plots text that shows horizontal distance
    hMill = plot(2*cos(0:pi/50:2*pi), 2*sin(0:pi/50:2*pi),'LineWidth',10,'Color',[148/255,0/255,211/255]);
    hold off
    axis( [ -25, 200, -25, 350] ); %axis limits
    legend([hZenith, hMISA],'Zenith Beam','MISA Beam');
    title(['Millstone Hill ISR Beam Orientation at ',num2str(pointAltitude),' km']);
    xlabel('Distance from Millstone Hill ISR (km)');
    ylabel('Altitude (km)');
    set(gca,'fontweight',font_Weight, 'FontSize',font_Size);

%prep figure #
    % figure('units','normalized','outerposition',[0 0 1 1]); %maximizes figure window
%     figure(figNum);
%     figNum = figNum+1;
    subplot(1,2,2)
    %Now drawing line of interest
    hZenithGeo = plot(longMillstone,latMillstone,'*','Color','k','MarkerSize',17); %plots a point with a black *
    hold on;
    hMISAGeo = plot(longMillstoneMISA,latMillstoneMISA,'*','Color','b','MarkerSize',17); %plots a point with a blue *
    hold on;
    plot( linspace(longMillstone,longMillstoneMISA,avg_anyAngle_N) ,linspace(latMillstone,latMillstoneMISA,avg_anyAngle_N),'LineStyle','--','LineWidth',1.5,'Color','r'); %plots a point with a blue *
    hold on;
        text( (max([longMillstoneMISA,longMillstone])-min([longMillstoneMISA,longMillstone]))/2+min([longMillstoneMISA,longMillstone])+.07 , (max([latMillstoneMISA,latMillstone])-min([latMillstoneMISA,latMillstone]))/2+min([latMillstoneMISA,latMillstone])+.07 ,[num2str(round(pointAltitude.*sind(90-MillstoneMISA_elev)./sind(MillstoneMISA_elev))),' km Separation'],...
        'Color','k','fontweight','bold','FontSize',18); %plots text that shows horizontal distance
    hold on;
    %Plot lines of contients
    geoshow(coastLines_lat,coastLines_long,'Color','k') %add continental outlines in black  
    hold off
    xlabel('Longitude (arcdeg)'); %old: ,'fontweight','bold','FontSize',12 for all
    ylabel('Latitude (arcdeg)'); %sped up significantly when merged above
    string_Title = ['Geographic Beam Orientation with Beams at ',num2str(pointAltitude),' km']; %create mecha title
    if( FLG_geomagneticCoords == 1)
        string_Title = [string_Title,' [Mag Coord Aligned]']; %tack on mag coordinate warning
    end
    title(string_Title);
    legend([hZenithGeo, hMISAGeo],'Zenith Beam','MISA Beam');
    axis( [ min([longMillstoneMISA,longMillstone])-3, max([longMillstoneMISA,longMillstone])+3, min([latMillstoneMISA,latMillstone])-3, max([latMillstoneMISA,latMillstone])+3] ); %axis limits 
%     xticks( min(plotLongRange):plotLongRange_autoTick:max(plotLongRange) ); %creates x ticks automagically
%     yticks( min(plotLatRange):plotLatRange_autoTick:max(plotLatRange) ); %creates y ticks automagically
    set(gca,'fontweight',font_Weight, 'FontSize',font_Size);


%[148/255,0/255,211/255]    c
%text(2,8,'A Simple Plot','Color','red','FontSize',14)
end

%adjust geographic coords to geomagnetic coords for given constants
if( FLG_geomagneticCoords == 1 ) %convert data to geomagnetic coords if desired
    [latMillstone,longMillstone] = sFUN_geoToGeomag(dateOfData,magF,latMillstone,longMillstone); %convert from geographic to geomagnetic coordinates
    [latMillstoneMISA,longMillstoneMISA] = sFUN_geoToGeomag(dateOfData,magF,latMillstoneMISA,longMillstoneMISA); %convert from geographic to geomagnetic coordinates
    [latPoints,longPoints] = sFUN_geoToGeomag(dateOfData,magF,latPoints,longPoints); %convert from geographic to geomagnetic coordinates
%     [plotLatRange,plotLongRange] = sFUN_geoToGeomag(dateOfData,magF,plotLatRange,plotLongRange); %convert from geographic to geomagnetic coordinates
end

%random noise background parameters
if( FLG_randomData == 1 )
    noise_Background_Mean = 0; %delta_sTEC, avg of noise background
    noise_Background_STDEV = 1/3.85; %delta_sTEC, standard dev of noise background, GOAL: at more than 99% of sTEC is between +/-1 delta_STEC by #*STDEV
    noise_Background_STDEV = 0.3168; %delta_sTEC, standard dev of noise background taken from all sTEC data
    sTEC = normrnd(noise_Background_Mean,noise_Background_STDEV,size(sTEC)); %fill in random noise into the sTEC
end

%Remove asterisk/black line from plots if plots aren't in the zone of Millstone Hill (aka. Japan)
if( latMillstone < min(plotLatRange) || latMillstone > max(plotLatRange) || longMillstone < min(plotLongRange) || longMillstone > max(plotLongRange) )
    latMillstone = NaN;
    longMillstone = NaN;
    latMillstoneMISA = NaN;
    longMillstoneMISA = NaN; %NaN out so it won't plot anything
    
end


%% Create date range hours for full day time limit
dateRange_zeroHr_YrMonDay = sFUN_dayNumber_to_Date(dateRange_zeroHr); %get the Yr/Mon/Day of the zero hour day
dateRange_zeroHr_YrMonDay_monWord = sFUN_monthNum_to_word(dateRange_zeroHr_YrMonDay(2),1,1); %1's request abbreviation with no dots (Jan, Feb, Mar...)
dateRange_zeroHr_YrMonDay_string = [dateRange_zeroHr_YrMonDay_monWord,' ',num2str(dateRange_zeroHr_YrMonDay(3))]; %put in the month
if( dateRange_zeroHr_YrMonDay(3) == 1 || dateRange_zeroHr_YrMonDay(3) == 21 )
    dateRange_zeroHr_YrMonDay_string = [dateRange_zeroHr_YrMonDay_string,'st, ']; %appropriate abbrevs for beauty
elseif( dateRange_zeroHr_YrMonDay(3) == 2 || dateRange_zeroHr_YrMonDay(3) == 22 )
    dateRange_zeroHr_YrMonDay_string = [dateRange_zeroHr_YrMonDay_string,'nd, ']; %appropriate abbrevs for beauty
elseif( dateRange_zeroHr_YrMonDay(3) == 3 || dateRange_zeroHr_YrMonDay(3) == 33 )
    dateRange_zeroHr_YrMonDay_string = [dateRange_zeroHr_YrMonDay_string,'rd, ']; %appropriate abbrevs for beauty
else
    dateRange_zeroHr_YrMonDay_string = [dateRange_zeroHr_YrMonDay_string,'th, ']; %appropriate abbrevs for beauty
end
dateRange_zeroHr_MonDay_string = dateRange_zeroHr_YrMonDay_string(1:strfind(dateRange_zeroHr_YrMonDay_string,',')-1);
dateRange_zeroHr_YrMonDay_string = [dateRange_zeroHr_YrMonDay_string,num2str(dateRange_zeroHr_YrMonDay(1))]; %create a human-readable date string

dateRange_full = sFUN_dayNumber_to_Date_MULTIPLE(dateRange,1); %get the full date range
dateRange_full_zeroHrMatch = ( dateRange_full(:,1) == dateRange_zeroHr(:,1) ) & ( dateRange_full(:,2) == dateRange_zeroHr(:,2) ); %get when zero matches the full data
if( ( sum(dateRange_full_zeroHrMatch) == 0 ) || ( sum(dateRange_full_zeroHrMatch) > 1 ) )
    error(['Declared Zero Hour Date (using day number in yr) of ',num2str(dateRange_zeroHr(1)),'/',num2str(dateRange_zeroHr(2)),' is not within the date range ',num2str(dateRange(1,1)),'/',num2str(dateRange(1,2)),' to ',num2str(dateRange(2,1)),'/',num2str(dateRange(2,2)),'.']);
end
cntr = 0; %prep it up
FLG_break = 0; %prep the break flag
for( i = 1:size(dateRange_full,1) )
    if( dateRange_full_zeroHrMatch(i) == 1 ) %look for when it matches
        FLG_break = 1; %break flag on!
    end
    if( FLG_break == 0 )
        cntr = cntr + 1; %increment
    end
end
dateRange_numDays = size(dateRange_full,1); %days, record number of days dealing with
dateRange_zeroHr_hrOffset = cntr*24; %put in the hour offset that makes a 0 to 72 hour range into the right -24 to 48 hour range (ALWAYS SUBTRACT)
dateRange_zeroHr_hrBounds = [-cntr*24,dateRange_numDays*24 - dateRange_zeroHr_hrOffset]; %create hour date range and put zero hour at desired time
dateRange_zeroHr_hrs = (0:24:dateRange_numDays*24) - dateRange_zeroHr_hrOffset; %hr, the hours that split the day wrt the zero hour day

%% Compare with ISR Data - ZENITH - DATA IMPORT
if( FLG_ISRdata_Zenith == 1 )

load('Zenith_height.mat')
load('Zenith_SNR.mat')
load('Zenith_SNR_bp.mat')
load('Zenith_vel.mat')
load('Zenith_time.mat') %time was adjusted so 0 hour is on May 7th, we will adjust to be same as this date system
Zenith_timeOffset = dateRange_zeroHr(2); %days, offset applied to get the dates the same between these two May 7th
Zenith_time = Zenith_timeOffset + Zenith_time./24; %days, convert to same date system
Zenith_avgRange = 25; %km, average up and down that amount

Zenith_threeHun = pointAltitude; %km, center goal altitude (defined above)
Zenith_threeHunIndex = find(min(abs(Zenith_height - Zenith_threeHun)) == abs(Zenith_height - Zenith_threeHun)); %Zenith height index that corresponds to 300 km (47 used, 299 km, now dynamic)

Zenith_avgRangeIndex = [find( min(abs( (Zenith_height(Zenith_threeHunIndex)-Zenith_avgRange) - Zenith_height)) == abs(Zenith_height(Zenith_threeHunIndex)-Zenith_avgRange - Zenith_height) )...
    ; find( min(abs(Zenith_height(Zenith_threeHunIndex)+Zenith_avgRange - Zenith_height)) == abs(Zenith_height(Zenith_threeHunIndex)+Zenith_avgRange - Zenith_height) )]; 
%upper and lower indexes 

for(i = min(Zenith_avgRangeIndex):max(Zenith_avgRangeIndex) )
    Zenith_SNR_bp(:,i) = Zenith_SNR_bp(:,i) - mean(Zenith_SNR_bp(:,i)); %remove the mean
end
Zenith_SNR_threeHun_AVGD = mean(Zenith_SNR_bp(:,Zenith_avgRangeIndex(1):Zenith_avgRangeIndex(2)),2); %average to try to mimic ionosphere TEC readings

time_min = (min(timeUnique)); %days, UT I think
time_max = (max(timeUnique)); %days, UT I think

jm = find( min(abs(time_min - Zenith_time)) == abs(time_min - Zenith_time) ); %get starting Zenith time that corresponds to here
jn = find( min(abs(time_max - Zenith_time)) == abs(time_max - Zenith_time) ); %get ending Zenith time that corresponds to here

Zenith_time = Zenith_time(jm:jn); %snip out data that is not comparable
Zenith_SNR = Zenith_SNR(jm:jn,:); %snip out data that is not comparable
Zenith_SNR_bp = Zenith_SNR_bp(jm:jn,:); %snip out data that is not comparable
Zenith_SNR_threeHun_AVGD = Zenith_SNR_threeHun_AVGD(jm:jn,:); %snip out data that is not comparable

tj = find( min( abs(min(Zenith_time) - timeUnique) ) == abs(min(Zenith_time) - timeUnique) ); %find min timeUnique that fits ISR time
tk = find( min( abs(max(Zenith_time) - timeUnique) ) == abs(max(Zenith_time) - timeUnique) ); %find max timeUnique that fits ISR time
timeUnique_limISR_Zenith = timeUnique(tj:tk); %cut out time that matches ISR time

Zenith_SNR_threeHun_AVGD_expanded = interp1(Zenith_time,Zenith_SNR_threeHun_AVGD,timeUnique_limISR_Zenith,'spline'); %dB, expand Zenith_SNR_AVGD to be the same size as the sTEC_combined_limISR var

if( mod(round((min(Zenith_time)-dateRange_zeroHr(2))*24),2) ~= 0 )
    plotTime_ISR_autoTick = round((min(Zenith_time)-dateRange_zeroHr(2))*24)-1; %adjust to be even
else
    plotTime_ISR_autoTick = round((min(Zenith_time)-dateRange_zeroHr(2))*24); %leave as-is
end
plotTime_ISR_autoTick = (plotTime_ISR_autoTick:2:round((max(Zenith_time)-dateRange_zeroHr(2))*24))';

end %*********END OF FLG_ISRdata_Zenith*********


%% Compare with ISR Data - MISA - DATA IMPORT
if( FLG_ISRdata_MISA == 1 )

load('MISA_height.mat')
load('MISA_SNR.mat')
load('MISA_SNR_bp.mat')
load('MISA_vel.mat')
load('MISA_time.mat') %time was adjusted so 0 hour is on May 7th, we will adjust to be same as this date system

MISA_timeOffset = dateRange_zeroHr(2); %days, offset applied to get the dates the same between these two May 7th
MISA_time = MISA_timeOffset + MISA_time./24; %days, convert to same date system
MISA_avgRange = 25; %km, average up and down that amount
MISA_threeHun = pointAltitude; %km, center goal altitude (defined above)

MISA_threeHunIndex = find(min(abs(MISA_height - MISA_threeHun)) == abs(MISA_height - MISA_threeHun)); %MISA height index that corresponds to 300 km (47 used, 282 km, now dynamic)
MISA_avgRangeIndex = [find( min(abs( (MISA_height(MISA_threeHunIndex)-MISA_avgRange) - MISA_height)) == abs(MISA_height(MISA_threeHunIndex)-MISA_avgRange - MISA_height) )...
    ; find( min(abs(MISA_height(MISA_threeHunIndex)+MISA_avgRange - MISA_height)) == abs(MISA_height(MISA_threeHunIndex)+MISA_avgRange - MISA_height) )]; 
%upper and lower indexes 

MISA_SNR_threeHun_AVGD = mean(MISA_SNR_bp(:,MISA_avgRangeIndex(1):MISA_avgRangeIndex(2)),2); %average to try to mimic ionosphere TEC readings

time_min = (min(timeUnique)); %days, UT I think
time_max = (max(timeUnique)); %days, UT I think

jm = find( min(abs(time_min - MISA_time)) == abs(time_min - MISA_time) ); %get starting MISA time that corresponds to here
jn = find( min(abs(time_max - MISA_time)) == abs(time_max - MISA_time) ); %get ending MISA time that corresponds to here

MISA_time = MISA_time(jm:jn); %snip out data that is not comparable
MISA_SNR = MISA_SNR(jm:jn,:); %snip out data that is not comparable
MISA_SNR_bp = MISA_SNR_bp(jm:jn,:); %snip out data that is not comparable
MISA_SNR_threeHun_AVGD = MISA_SNR_threeHun_AVGD(jm:jn,:); %snip out data that is not comparable

tj = find( min( abs(min(MISA_time) - timeUnique) ) == abs(min(MISA_time) - timeUnique) ); %find min timeUnique that fits ISR time
tk = find( min( abs(max(MISA_time) - timeUnique) ) == abs(max(MISA_time) - timeUnique) ); %find max timeUnique that fits ISR time
timeUnique_limISR_MISA = timeUnique(tj:tk); %cut out time that matches ISR time

MISA_SNR_threeHun_AVGD_expanded = interp1(MISA_time,MISA_SNR_threeHun_AVGD,timeUnique_limISR_MISA,'spline'); %dB, expand MISA_SNR_AVGD to be the same size as the sTEC_combined_limISR var

if( mod(round((min(MISA_time)-dateRange_zeroHr(2))*24),2) ~= 0 )
    plotTime_ISR_autoTick = round((min(MISA_time)-dateRange_zeroHr(2))*24)-1; %adjust to be even
else
    plotTime_ISR_autoTick = round((min(MISA_time)-dateRange_zeroHr(2))*24); %leave as-is
end
plotTime_ISR_autoTick = (plotTime_ISR_autoTick:2:round((max(MISA_time)-dateRange_zeroHr(2))*24))';

end %*********END OF FLG_ISRdata_MISA*********


%% Import and Prep JEDI Data
if( FLG_JEDI_DATA == 1)
    
% Import Joule Heating Data
JEDI_data = importdata('20130505_to_20130509_data_v2'); %import data

JEDI_date = JEDI_data(:,1); %date in numerical format (20130505 = 2013, 05, 05)
JEDI_time = JEDI_data(:,2); %hr, time in UT
JEDI_energyFlux = JEDI_data(:,3); % integrated energy ?ux from precipitation
JEDI_jouleHeating = JEDI_data(:,4); % integrated Joule heat input 
clear JEDI_data %cleanup

% Adjust the time variable so it is continuous and matches current format
dateRangeZeroHr_YrMonDay = sFUN_dayNumber_to_Date(dateRange_zeroHr); %change zero hour day number to yr/mon/day version
temp_monHoldr = num2str(dateRangeZeroHr_YrMonDay(2)); %pull out the month
if( length(temp_monHoldr) == 1 )
    temp_monHoldr = ['0',temp_monHoldr]; %tack on a 0 so like 05 is the 5th month
end
temp_dayHoldr = num2str(dateRangeZeroHr_YrMonDay(3)); %pull out the day
if( length(temp_dayHoldr) == 1 )
    temp_dayHoldr = ['0',temp_dayHoldr]; %tack on a 0 so like 05 is the 5th day
end
JEDI_dateRangeZeroHr_YrMonDay = str2double([num2str(dateRangeZeroHr_YrMonDay(1)),temp_monHoldr,temp_dayHoldr]); %create a JEDI-friendly string then turn it into a number again bam bam bam
JEDI_time = JEDI_time + (JEDI_date-JEDI_dateRangeZeroHr_YrMonDay).*24; %hr, adjusts time so it does not repeat
clear temp_monHoldr temp_dayHoldr

% Clip out desired date range WRT full date range
tk = find( min(abs( min(dateRange_zeroHr_hrBounds) - JEDI_time)) == abs( min(dateRange_zeroHr_hrBounds) - JEDI_time) );
tl = find( min(abs( max(dateRange_zeroHr_hrBounds) - JEDI_time)) == abs( max(dateRange_zeroHr_hrBounds) - JEDI_time) );

JEDI_time_full = JEDI_time(tk:tl); %clip
JEDI_energyFlux_full = JEDI_energyFlux(tk:tl); %clip
JEDI_jouleHeating_full = JEDI_jouleHeating(tk:tl); %clip

% Clip out Desired Date Range WRT Zenith/MISA
tk = find( min(abs( (min(timeUnique_limISR_MISA)-floor(mean(timeUnique_limISR_MISA))).*24 - JEDI_time)) == abs( (min(timeUnique_limISR_MISA)-floor(mean(timeUnique_limISR_MISA))).*24 - JEDI_time) );
tl = find( min(abs( (max(timeUnique_limISR_MISA)-floor(mean(timeUnique_limISR_MISA))).*24 - JEDI_time)) == abs( (max(timeUnique_limISR_MISA)-floor(mean(timeUnique_limISR_MISA))).*24 - JEDI_time) );

JEDI_time = JEDI_time(tk:tl); %clip
JEDI_energyFlux = JEDI_energyFlux(tk:tl); %clip
JEDI_jouleHeating = JEDI_jouleHeating(tk:tl); %clip

% Filter the Data
JEDI_time_Delta = mean(JEDI_time(2:end) - JEDI_time(1:end-1)); %hr, calc time delta

JEDI_fp=JEDI_highpass_Freq; %1/hr, high-pass cutoff freq
JEDI_n = JEDI_nOrder; %set window order

JEDI_f= 1/JEDI_time_Delta; %1/hr, the sampling frequency, based off of the time delta calc'd

JEDI_wp=2*JEDI_fp/JEDI_f; % Normalizing the frequencies (Matlab is 0 to 1)

%Calculation of filter coefficients
% [b,a]=fir1(n,wp,'high'); %This takes the order and lower cut off
%Uses the default hamming window
% frequency ORIG
JEDI_W = hann(JEDI_n+1);
[JEDI_b,JEDI_a] = fir1(JEDI_n,JEDI_wp,'high',JEDI_W);
%Applys Hanning Window for shiggles and giggles

JEDI_jouleHeating_full_bp(:,:) = filtfilt(JEDI_b,JEDI_a,JEDI_jouleHeating_full(:,:)); %Applies the filter
JEDI_energyFlux_full_bp(:,:) = filtfilt(JEDI_b,JEDI_a,JEDI_energyFlux_full(:,:)); %Applies the filter
JEDI_jouleHeating_bp(:,:) = filtfilt(JEDI_b,JEDI_a,JEDI_jouleHeating(:,:)); %Appies the filter
JEDI_energyFlux_bp(:,:) = filtfilt(JEDI_b,JEDI_a,JEDI_energyFlux(:,:)); %Applies the filter
    
end


%% Load and prep AMPERE data
%JEDI data is heavily processed AMPERE btw

if( FLG_AMPERE_DATA == 1 )
    
    %info on how files were created:
%        for zt=0,239 do begin
%         printf,5,utout(zt)
%         for zlat=0,49 do begin
%           printf,5,clat(zlat)
%           for zlong=0,119 do begin
%             printf,5,qsigp(zt,zlong,zlat),qsigh(zt,zlong,zlat),qjh(zt,zlong,zlat),qpot(zt,zlong,zlat),qjp(zt,zlong,zlat)
%           endfor
%         endfor
%       endear
% In other words, each block of data starts with the UT, then a geographic latitude beginning with 40 degrees, and then a block of 120 lines, each containing six quantities:
% Pedersen Conductance, Hall Conductance, Joule Heat in ergs/cm2-sec, electric potential, field-aligned current
% The 120 lines correspond to geographic longitudes at 3-degree intervals beginning with 0 degrees.

%Code to read year from file name, decided unneeded
AMPERE_date = zeros(size(fileName_AMPERE,1),2,'int16'); %good till yr 32,767
for( i = 1:size(fileName_AMPERE,1) )
    jk = strfind(fileName_AMPERE(i,:),'_'); %get all the _'s
    tempYr = str2double(fileName_AMPERE(i, (jk(end)+1):((jk(end)+1)+3) )); %yr, get the year
    tempMon = str2double(fileName_AMPERE(i, (jk(end)+1+4):((jk(end)+1)+5) )); %mon, get the month
    tempDay = str2double(fileName_AMPERE(i, (jk(end)+1+6):end )); %day, get the day
    AMPERE_date(i,1) = tempYr; %yr, save year
    AMPERE_date(i,2) = sFUN_dateToDayNum([tempYr,tempMon,tempDay]); %dayNum, save day number (what we use for ez-er math)
end

AMPERE_timeStep = 0.1; %hr, set the time step to split hours into
AMPERE_time_local = repmat(single(0:AMPERE_timeStep:(24-AMPERE_timeStep))',size(fileName_AMPERE,1),1); %hr, AMPERE data time step that matches data file (for checking)
AMPERE_time_range = single(0:AMPERE_timeStep:(size(fileName_AMPERE,1)*24))'; %hr, AMPERE data is every 0.1 hr
AMPERE_time_range(end) = []; %hr, remove last entry b/c AMPERE counts 0 to 23.9 so the last day's entry of 24/48/72 shouldn't exist
AMPERE_time_range = AMPERE_time_range - (find(dateRange_full_zeroHrMatch == 1)-1)*24; %hr, AMPERE data 0 hr on the zero hour date given above (e.g. May 7th, 2013)
AMPERE_lat_range = single((40:1:89)'); %degc, lat, ampere's lat range
AMPERE_long_range = single(circshift((0:3:357)'-180,round(length(0:3:357)/2))); %degc, long, ampere's long range (starts at 0 long, increases, wraps around to -180, etc)
AMPERE_time_len = length(AMPERE_time_range); %length of time range (how many unique times there are)
AMPERE_lat_len = length(AMPERE_lat_range); %length of lat range (how many unique lats there are)
AMPERE_long_len = length(AMPERE_long_range); %length of long range (how many unique longs there are)
AMPERE_data_len = AMPERE_time_len*AMPERE_lat_len*AMPERE_long_len; %length of data vector (time/lat/long repeated in varying ways to match data)
AMPERE_time = zeros(AMPERE_data_len,1,'single'); %hr, time that corresponds to data points (50*120=6,000 hr -24 then 6,000 hr -23.9 etc.)
for( i = 1:length(AMPERE_time_range) ) %loop to make time data
    AMPERE_time( 1+(i-1)*AMPERE_lat_len*AMPERE_long_len:i*AMPERE_lat_len*AMPERE_long_len ) = repmat(AMPERE_time_range(i),AMPERE_lat_len*AMPERE_long_len,1); %hr, replicate each AMPERE_time_range entry 6,000 times (50lat*120long) to match data cadence
end
AMPERE_lat_local = repmat(AMPERE_lat_range,AMPERE_time_len,1); %degc, AMPERE's lat range that matches file readout for checking
AMPERE_lat = zeros(AMPERE_lat_len*AMPERE_long_len,1,'single'); %degc, lat that corresponds to data points prep (40 120 times then 41 120 times etc.)
for( i = 1:length(AMPERE_lat_range) ) %loop to make lat data
   	AMPERE_lat((i-1)*AMPERE_long_len+1:i*AMPERE_long_len) = repmat(AMPERE_lat_range(i),AMPERE_long_len,1); %degc, put in (40 120 times then 41 120 times etc.)
end
AMPERE_lat = repmat(AMPERE_lat,AMPERE_time_len,1); %degc, replicate to cover all of the time steps
AMPERE_long = repmat(AMPERE_long_range,AMPERE_time_len*AMPERE_lat_len,1); %degc, long that corresponds to data points

%key AMPERE data descrip here
AMPERE_data = zeros(AMPERE_data_len,AMPERE_data_num,'single'); %1 Pedersen Conductance / 2 Hall Conductance / 3 Joule Heat (ergs/(cm^2*sec)) / 4 Electric Potential / 5 Field-Algined Current

%load saved AMPERE to speed up
AMPERE_loaded = zeros(size(fileName_AMPERE,1),1,'uint8'); %prep up, 0 or 1
for( i = 1:size(fileName_AMPERE,1) ) %check for each file and load if needed
    fileMATName = ['AMPERE_',num2str(AMPERE_date(i,1)),'_',num2str(AMPERE_date(i,2)),'.mat']; %create file name it would be called
    k = exist(fileMATName,'file'); %see if it is there
    if( k == 2 )
        load(fileMATName); %load saved file
        AMPERE_data( 1+(i-1)*length(0:AMPERE_timeStep:(24-AMPERE_timeStep))*AMPERE_lat_len*AMPERE_long_len:i*length(0:AMPERE_timeStep:(24-AMPERE_timeStep))*AMPERE_lat_len*AMPERE_long_len , :) = AMPERE_data_temp; %save the data in one big file!
        AMPERE_loaded(i) = 1; %note it was loaded
    end
end
clear AMPERE_data_temp
% AMPERE_loaded = 0; %for testing, makes it read anew
%else read AMPERE files to get the data needed
if( sum(AMPERE_loaded) ~= size(fileName_AMPERE,1) ) %if there's days missing we gotta read it in, gonna take about 220 sec on my comp for up to 4 files at once
    fprintf('Some or None of "cached" pre-read AMPERE files found; beginning to read AMPERE data - est. 4 min on fast comp\n');
    tic
    AMPERE_data_par = cell(size(fileName_AMPERE,1),1);  %each file gets its own data, need cell for parfor thx matlab
    for( i = 1:size(fileName_AMPERE,1) )
        AMPERE_data_par{i} = zeros(length(0:AMPERE_timeStep:(24-AMPERE_timeStep))*AMPERE_lat_len*AMPERE_long_len,AMPERE_data_num,'single'); %prep up
    end
    
    parfor( i = 1:size(fileName_AMPERE,1) ) %open each file and roll through them
        %Import raw AMPERE data
        fid = fopen(fileName_AMPERE(i,:),'r'); %open file pointer
            AMPERE_raws = textscan(fid, '%s','delimiter','\n'); %multiple headers of varying sizes make reading this real annoying
        fclose(fid); %close up shop
        AMPERE_raws = AMPERE_raws{1,1}; %pull out the data to make it easier to work with

        cntr_time = 0; %prep counter for time (240*# of files)
        cntr_lat = AMPERE_lat_len+1; %prep counter for lats (50 per time) - starts off 50+1 for logic
        cntr_long = 0; %prep counter for longs (120 per lat)
        for(j = 1:length(AMPERE_raws) ) %run through each line
            k = strfind(AMPERE_raws{j,1},'.'); %get place of each .
            %Time & Lat header info have 1 . per line and time is before lat when it occurs
            %long inputs have 5 . per line
            if( length(k) == 1 &&  (cntr_lat >= AMPERE_lat_len) ) %if this is true it is a time entry
                cntr_time = cntr_time + 1; %increment
                cntr_lat = 0; %set to 0 to begin the fun
                % str2double(AMPERE_raws{j,1}) ~= AMPERE_time_local(cntr_time) had issues with differing levels of precision
    %             if(  ismembertol(str2double(AMPERE_raws{j,1}),AMPERE_time_local(cntr_time)) ~= 1 )
    %                 error(['Expected time of ',num2str(AMPERE_time_local(cntr_time)),' is not same as read line of ',AMPERE_raws{j,1},' on line #',num2str(j),' of file ',fileName_AMPERE(i,:)]); %Report an error in data mismatch expected vs read
    %             end
    %         DISABLED ALL CHECKS TO IMPROVE SPEED FARTHER

            elseif( length(k) == 1 &&  (cntr_lat < AMPERE_lat_len) ) %if this is true it is a lat entry
                cntr_lat = cntr_lat + 1; %increment lat cntr
                % str2double(AMPERE_raws{j,1}) ~= AMPERE_lat_local(cntr_lat) might have same issue as time had with differing levels of precision
    %             if( ismembertol(str2double(AMPERE_raws{j,1}),AMPERE_lat_local(cntr_lat)) ~= 1 )
    %                 error(['Expected lat of ',num2str(AMPERE_lat_local(cntr_time)),' is not same as read line of ',AMPERE_raws{j,1},' on line #',num2str(j),' of file ',fileName_AMPERE(i,:)]); %Report an error in data mismatch expected vs read
    %             end

            else %otherwise it is a long entry
                cntr_long = cntr_long + 1; %increment
    %             C = strsplit(AMPERE_raws{j,1},' '); %split up the 5 numbers
    %             if( length(C) ~= 5 )
    %                 error(['Expected number of data pts for a data line is ',length(C),' not 5. The read line is ',AMPERE_raws{j,1},' on line #',num2str(j),' of file ',fileName_AMPERE(i,:)]); %Report an error in data mismatch expected vs read
    %             end
                AMPERE_data_par{i}(cntr_long,:) = str2double(strsplit(AMPERE_raws{j,1},' ')); %place the 5 data points in            
            end
        end
    end
    for( i = 1:size(fileName_AMPERE,1) )
        AMPERE_data( 1+(i-1)*length(0:AMPERE_timeStep:(24-AMPERE_timeStep))*AMPERE_lat_len*AMPERE_long_len:i*length(0:AMPERE_timeStep:(24-AMPERE_timeStep))*AMPERE_lat_len*AMPERE_long_len , :) = AMPERE_data_par{i}(:,:); %read in the data
        %now save the data to use again
        AMPERE_data_temp = AMPERE_data_par{i}(:,:); %read in the data to temp var to save
        fileMATName = ['AMPERE_',num2str(AMPERE_date(i,1)),'_',num2str(AMPERE_date(i,2)),'.mat'];
        save(fileMATName,'AMPERE_data_temp','-v7.3'); %saves the mega variable
        %v7 is used because of reasonable file size. (7.3 used b/c just numbers yay!)
        %v7.3 (which is really HDF5) results in files going from 4 MB (v7) to 665 MB (v7.3)
        %Looks like it is due to my use of cell arrays here - will apply v7.3 as the final product of just numbers
    end
    
    clear AMPERE_data_par AMPERE_data_temp
    fprintf('Time to read AMPERE data:\tRuntime: %.3f sec\n',toc);
end

%Clear out AMPERE data less than the minimum
k = min(gif_AMPERE_jouleHeating_cap) > AMPERE_data(:,AMPERE_jouleHeating_pos); %find entries less than the min plotting number (clear it up)
AMPERE_data(k,:) = []; %delete entries that are less than the min plotting number (clear it up)
AMPERE_time(k) = []; %delete entries that are less than the min plotting number (clear it up)
AMPERE_lat(k) = []; %delete entries that are less than the min plotting number (clear it up)
AMPERE_long(k) = []; %delete entries that are less than the min plotting number (clear it up)

end %*********END OF FLG_AMPERE_DATA*********


%% Get Solar Activity/Geomagnetic Stuff
if(FLG_solarGeoStuff == 1 )
    if( FLG_ISRdata_Zenith == 1 )
        ISR_time_min_Zenith = (Zenith_time(1)-round(mean(Zenith_timeOffset)))*24; %UT hr, min time to compare to
        ISR_time_max_Zenith = (Zenith_time(end)-round(mean(Zenith_timeOffset)))*24; %UT hr, min time to compare to    
    end
    
    %===========================Kp Index===================================
    %Kp INDEX HERE
    [Kp_dateRange_Full_Month,Kp_output] = sFUN_KpIndexGET(dateRange,1); %get Kp index for the date range
    Kp_time = ((3:3:dateRange_numDays*24)-dateRange_zeroHr_hrOffset); %time for Kp measurements (0-3 is assumed to be relevant for hr 3, etc...)
    %Kp PLOT HERE
    %prep figure #
    % figure('units','normalized','outerposition',[0 0 1 1]); %maximizes figure window (caused graphical issues on Windows 10 + AMD GPU or Intel GPU)
    figure(figNum);
    figNum = figNum+1;
    Kp_output_Plot = repelem(reshape(Kp_output',[size(Kp_output,1)*size(Kp_output,2),1]),2); %replicate teh Kp for plotting
    Kp_time_Plot = repelem(Kp_time,2); %replicate the hours for plotting purposes
    Kp_time_Plot = [min(Kp_time_Plot)-3,Kp_time_Plot(1:end-1)]; %readjust so the hour range matches what is real (e.g. Kp lasts 3 hr, so 0-3 hr is same Kp value)
    plot( Kp_time_Plot,Kp_output_Plot ,'LineWidth',3 )
    if( FLG_ISRdata_Zenith == 1 )
        hold on;
        plot( repmat(ISR_time_min_Zenith,size(0:.5:(max(Kp_output_Plot)+.5))) , 0:.5:(max(Kp_output_Plot)+.5),'LineWidth',1.75,'Color','r'); %plot red lines showing ISR data time
        hold on;
        plot( repmat(ISR_time_max_Zenith,size(0:.5:(max(Kp_output_Plot)+.5))) , 0:.5:(max(Kp_output_Plot)+.5),'LineWidth',1.75,'Color','r'); %plot red lines showing ISR data time
        hold off;
    end
%     bar( Kp_time-1,reshape(Kp_output',[size(Kp_output,1)*size(Kp_output,2),1]) ) %bar graph attempt...
    xlabel('Time (hr)'); %old: ,'fontweight','bold','FontSize',12 for all
    ylabel('Kp Index, 3 Hr Resolution'); %sped up significantly when merged above
    title(['Kp Index for ',num2str(Kp_dateRange_Full_Month(1,2)),'/',num2str(Kp_dateRange_Full_Month(1,3)),...
        '/',num2str(Kp_dateRange_Full_Month(1,1)),' to ',num2str(Kp_dateRange_Full_Month(end,2)),...
        '/',num2str(Kp_dateRange_Full_Month(end,3)),'/',num2str(Kp_dateRange_Full_Month(end,1)),...
        ' (M/D/Y) with 0 Hr on ',dateRange_zeroHr_YrMonDay_string]);
    if( mod(min(Kp_time_Plot),2) == 0 )
        Kp_time_axis_min = min(Kp_time_Plot); %is even, good to go
    else
        Kp_time_axis_min = min(Kp_time_Plot)+1; %is odd, make even
    end
    
    if( mod(max(Kp_time_Plot),2) == 0 )
        Kp_time_axis_max = max(Kp_time_Plot); %is even, good to go
    else
        Kp_time_axis_max = max(Kp_time_Plot)-1; %is odd, make even
    end
    Xaxisvar = Kp_time_axis_min:4:Kp_time_axis_max;
    set(gca,'xtick',Xaxisvar);
    axis([min(Xaxisvar) max(Xaxisvar) 0 (max(Kp_output_Plot)+.5)])
    set(gca,'fontweight',font_Weight, 'FontSize',font_Size);
    
    %=========================OMNI=========================================
    [OMNI_dateRange_Full_Month,OMNI_output] = sFUN_OMNIGET(dateRange,1); %get OMNI index for the date range
    % 1 - yr
    % 2 - day#
    % 3 - hr
    % 4 - min
    % 5 - Bz GSE (nT)
    % 6 - Flow speed (km/s)
    % 7 - Flow pressure (nPa)
    % 8 - AE Index (nT)
    % 9 - SYM/H Index (nT)
    OMNI_time = OMNI_output(:,2) + OMNI_output(:,3)./24 + OMNI_output(:,4)./(24*60); %get the time in day#, ignore year atm
    OMNI_time_hr = (OMNI_time - dateRange_zeroHr(2)).*24; %0 hr on May 7th, 2013 (Day 127) custom hour numbering
    if( mod(round(min(OMNI_time_hr)),2) == 0 )
        OMNI_time_hr_axis_min = round(min(OMNI_time_hr)); %is even, good to go
    else
        OMNI_time_hr_axis_min = round(min(OMNI_time_hr))+1; %is odd, make even
    end
    if( mod(round(max(OMNI_time_hr)),2) == 0 )
        OMNI_time_hr_axis_max = round(max(OMNI_time_hr)); %is even, good to go
    else
        OMNI_time_hr_axis_max = round(max(OMNI_time_hr))-1; %is odd, make even
    end
    %OMNI PLOT HERE
    %prep figure #set(gca,'xticklabel',[])
    % figure('units','normalized','outerposition',[0 0 1 1]); %maximizes figure window (caused graphical issues on Windows 10 + AMD GPU or Intel GPU)
    figure(figNum);
    figNum = figNum+1;
    Xaxisvar = OMNI_time_hr_axis_min:4:OMNI_time_hr_axis_max;
    ax1 = subplot(5,1,1);
    plot( OMNI_time_hr,OMNI_output(:,5) );
    if( FLG_ISRdata_Zenith == 1 )
        hold on;
        yVarrr = linspace(min(OMNI_output(:,5)),max(OMNI_output(:,5)),10);
        plot( repmat(ISR_time_min_Zenith,size(yVarrr)) , yVarrr,'LineWidth',1.75,'Color','r'); %plot red lines showing ISR data time
        hold on;
        plot( repmat(ISR_time_max_Zenith,size(yVarrr)) , yVarrr,'LineWidth',1.75,'Color','r'); %plot red lines showing ISR data time
        hold off;
    end
    ylabel('Bz GSE (nT)'); %sped up significantly when merged above
    title(['OMNI-Sourced Data (1 min resolution) for ',num2str(OMNI_dateRange_Full_Month(1,2)),'/',num2str(OMNI_dateRange_Full_Month(1,3)),...
        '/',num2str(OMNI_dateRange_Full_Month(1,1)),' to ',num2str(OMNI_dateRange_Full_Month(end,2)),...
        '/',num2str(OMNI_dateRange_Full_Month(end,3)),'/',num2str(OMNI_dateRange_Full_Month(end,1)),...
        ' (M/D/Y) with 0 Hr on ',dateRange_zeroHr_YrMonDay_string]);
    set(ax1,'xtick',Xaxisvar);
    set(ax1,'xticklabels',[])
    axis([min(Xaxisvar) max(Xaxisvar) min(OMNI_output(:,5)) max(OMNI_output(:,5))])
    set(gca,'fontweight',font_Weight, 'FontSize',font_Size);
    grid on
    ax2 = subplot(5,1,2);
    plot( OMNI_time_hr,OMNI_output(:,6) );
    if( FLG_ISRdata_Zenith == 1 )
        hold on;
        yVarrr = linspace(min(OMNI_output(:,6)),max(OMNI_output(:,6)),10);
        plot( repmat(ISR_time_min_Zenith,size(yVarrr)) , yVarrr,'LineWidth',1.75,'Color','r'); %plot red lines showing ISR data time
        hold on;
        plot( repmat(ISR_time_max_Zenith,size(yVarrr)) , yVarrr,'LineWidth',1.75,'Color','r'); %plot red lines showing ISR data time
        hold off;
    end
    ylabel('Vsw (km/s)'); %sped up significantly when merged above
    set(ax2,'xtick',Xaxisvar);
    set(ax2,'xticklabels',[])
    axis([min(Xaxisvar) max(Xaxisvar) min(OMNI_output(:,6)) max(OMNI_output(:,6))])
    set(gca,'fontweight',font_Weight, 'FontSize',font_Size);
    grid on
    ax3 = subplot(5,1,3);
    plot( OMNI_time_hr,OMNI_output(:,7) );
    if( FLG_ISRdata_Zenith == 1 )
        hold on;
        yVarrr = linspace(min(OMNI_output(:,7)),max(OMNI_output(:,7)),10);
        plot( repmat(ISR_time_min_Zenith,size(yVarrr)) , yVarrr,'LineWidth',1.75,'Color','r'); %plot red lines showing ISR data time
        hold on;
        plot( repmat(ISR_time_max_Zenith,size(yVarrr)) , yVarrr,'LineWidth',1.75,'Color','r'); %plot red lines showing ISR data time
        hold off;
    end
    ylabel('Psw (nPa)'); %sped up significantly when merged above
    set(ax3,'xtick',Xaxisvar);
    set(ax3,'xticklabels',[])
    axis([min(Xaxisvar) max(Xaxisvar) min(OMNI_output(:,7)) max(OMNI_output(:,7))])
    set(gca,'fontweight',font_Weight, 'FontSize',font_Size);
    grid on
    ax4 = subplot(5,1,4);
    plot( OMNI_time_hr,OMNI_output(:,8) );
    if( FLG_ISRdata_Zenith == 1 )
        hold on;
        yVarrr = linspace(min(OMNI_output(:,8)),max(OMNI_output(:,8)),10);
        plot( repmat(ISR_time_min_Zenith,size(yVarrr)) , yVarrr,'LineWidth',1.75,'Color','r'); %plot red lines showing ISR data time
        hold on;
        plot( repmat(ISR_time_max_Zenith,size(yVarrr)) , yVarrr,'LineWidth',1.75,'Color','r'); %plot red lines showing ISR data time
        hold off;
    end
    ylabel('AE (nT)'); %sped up significantly when merged above
    set(ax4,'xtick',Xaxisvar);
    set(ax4,'xticklabels',[])
    axis([min(Xaxisvar) max(Xaxisvar) min(OMNI_output(:,8)) max(OMNI_output(:,8))])
    set(gca,'fontweight',font_Weight, 'FontSize',font_Size);
    grid on
    ax5 = subplot(5,1,5);
    plot( OMNI_time_hr,OMNI_output(:,9) );
    if( FLG_ISRdata_Zenith == 1 )
        hold on;
        yVarrr = linspace(min(OMNI_output(:,9)),max(OMNI_output(:,9)),10);
        plot( repmat(ISR_time_min_Zenith,size(yVarrr)) , yVarrr,'LineWidth',1.75,'Color','r'); %plot red lines showing ISR data time
        hold on;
        plot( repmat(ISR_time_max_Zenith,size(yVarrr)) , yVarrr,'LineWidth',1.75,'Color','r'); %plot red lines showing ISR data time
        hold off;
    end
    ylabel('SYM/H (nT)'); %sped up significantly when merged above
    xlabel('Time (hr)'); %old: ,'fontweight','bold','FontSize',12 for all
    set(ax5,'xtick',Xaxisvar);
    axis([min(Xaxisvar) max(Xaxisvar) min(OMNI_output(:,9)) max(OMNI_output(:,9))])
    set(gca,'fontweight',font_Weight, 'FontSize',font_Size);
    grid on
    
    
    
    bs=(1/2);   % highpass cuttoff frequency (not sure what to make of it)
    %I think it is related to lf above (which is 1/2 )

    n=42; % order of the Hamming window used in the custom function
    % c = 3.32*pi; %Hamming constant
    % M = n/2;
    % bandwidth = c/M;
    %The above was included to investigate the bandwidth - 1/2. Related to
    %bs/lf?

    fp=bs; % stop band stoppin freq (or pass band passing freq, depending on how you look at it)
    timeDelta = mean(OMNI_time_hr(2:end)-OMNI_time_hr(1:end-1)); %hrs already
    f= 1/(timeDelta); % the sampling frequency, based off of the time delta calc'd

    wp=2*fp/f; % Normalizing the frequencies (Matlab is 0 to 1)

    %Calculation of filter coefficients
    % [b,a]=fir1(n,wp,'high'); %This takes the order and lower cut off
    %Uses the default hamming window
    % frequency ORIG
    W = hann(n+1);
    [b,a] = fir1(n,wp,'high',W);
    %Applys Hanning Window for shiggles and giggles


    OMNI_output_Bz_BP = filtfilt(b,a,OMNI_output(:,5)); %Appies the filter
    OMNI_output_AE_BP = filtfilt(b,a,OMNI_output(:,8)); %Appies the filter

    Ross_Lombscargle_optimized([ OMNI_time_hr , OMNI_output_Bz_BP ]); 
    fid = fopen('power_freq.txt', 'r');
        f_nm = textscan(fid, '%f %f %f %f %f', 'headerlines', 15);
        P_Zenith_SNR=f_nm{:,2};
        F_Zenith_SNR=f_nm{:,1};
        T_Zenith_SNR=f_nm{:,5}*60; % in minutes
        gs_Zenith_SNR=f_nm{:,4};
        max_Zenith_SNR=T_Zenith_SNR(P_Zenith_SNR==max(P_Zenith_SNR));
        K=length(P_Zenith_SNR);
        gf_Zenith_SNR=f_nm{:,3};
    fclose('all');

    %prep figure #
    % figure('units','normalized','outerposition',[0 0 1 1]); %maximizes figure window
    figure(figNum);
    figNum = figNum+1;

    plot(T_Zenith_SNR,P_Zenith_SNR,'b');
    xlim([0 plot_Freq_Lim]);
    hold on;
    % plot(T_SNR_08apr08_10_or,P_SNR_08apr08_10_or,'r');
    % hold on;
    plot(T_Zenith_SNR,gf_Zenith_SNR,'LineWidth',1.5,'color',[0.5 0.5 0.5]);
    xlim([0 plot_Freq_Lim]);  
    xlabel('Periods (min)');
    ylabel('Normalized Power');
    string_Title = ['Bz GSE High-pass Cutoff Freq of ',num2str(round(1/bs)),' Lomb-Scargle Periodogram -- 0 Hr on ',dateRange_zeroHr_YrMonDay_string];
    if( FLG_geomagneticCoords == 1)
        string_Title = [string_Title,' [Mag Coord Aligned]']; %tack on mag coordinate warning
    end
    title(string_Title);
    set(gca,'fontweight',font_Weight, 'FontSize',font_Size);
    set(gca,'xtick',Xaxisvar_SCARGLE);
    
    
    %===============JEDI data in-addition to OMNI plot=====================
    if( FLG_JEDI_DATA == 1 )
        %OMNI & JEDI PLOT HERE
        %prep figure #set(gca,'xticklabel',[])
        % figure('units','normalized','outerposition',[0 0 1 1]); %maximizes figure window (caused graphical issues on Windows 10 + AMD GPU or Intel GPU)
        figure(figNum);
        figNum = figNum+1;
        Xaxisvar = OMNI_time_hr_axis_min:4:OMNI_time_hr_axis_max;
        ax1 = subplot(6,1,1);
        plot( OMNI_time_hr,OMNI_output(:,5) );
        if( FLG_ISRdata_Zenith == 1 )
            hold on;
            yVarrr = linspace(min(OMNI_output(:,5)),max(OMNI_output(:,5)),10);
            plot( repmat(ISR_time_min_Zenith,size(yVarrr)) , yVarrr,'LineWidth',1.75,'Color','r'); %plot red lines showing ISR data time
            hold on;
            plot( repmat(ISR_time_max_Zenith,size(yVarrr)) , yVarrr,'LineWidth',1.75,'Color','r'); %plot red lines showing ISR data time
            hold off;
        end
        ylabel('Bz GSE (nT)'); %sped up significantly when merged above
        title(['OMNI (1 min resolution) & JEDI (',num2str(round(JEDI_time_Delta*60)),' min res.) Data for ',num2str(OMNI_dateRange_Full_Month(1,2)),'/',num2str(OMNI_dateRange_Full_Month(1,3)),...
            '/',num2str(OMNI_dateRange_Full_Month(1,1)),' to ',num2str(OMNI_dateRange_Full_Month(end,2)),...
            '/',num2str(OMNI_dateRange_Full_Month(end,3)),'/',num2str(OMNI_dateRange_Full_Month(end,1)),...
            ' (M/D/Y) with 0 Hr on ',dateRange_zeroHr_YrMonDay_string]);
        set(ax1,'xtick',Xaxisvar);
        set(ax1,'xticklabels',[])
        axis([min(Xaxisvar) max(Xaxisvar) min(OMNI_output(:,5)) max(OMNI_output(:,5))])
        set(gca,'fontweight',font_Weight, 'FontSize',font_Size);
        grid on
        ax2 = subplot(6,1,2);
        plot( OMNI_time_hr,OMNI_output(:,6) );
        if( FLG_ISRdata_Zenith == 1 )
            hold on;
            yVarrr = linspace(min(OMNI_output(:,6)),max(OMNI_output(:,6)),10);
            plot( repmat(ISR_time_min_Zenith,size(yVarrr)) , yVarrr,'LineWidth',1.75,'Color','r'); %plot red lines showing ISR data time
            hold on;
            plot( repmat(ISR_time_max_Zenith,size(yVarrr)) , yVarrr,'LineWidth',1.75,'Color','r'); %plot red lines showing ISR data time
            hold off;
        end
        ylabel('Vsw (km/s)'); %sped up significantly when merged above
        set(ax2,'xtick',Xaxisvar);
        set(ax2,'xticklabels',[])
        axis([min(Xaxisvar) max(Xaxisvar) min(OMNI_output(:,6)) max(OMNI_output(:,6))])
        set(gca,'fontweight',font_Weight, 'FontSize',font_Size);
        grid on
        ax3 = subplot(6,1,3);
        plot( OMNI_time_hr,OMNI_output(:,7) );
        if( FLG_ISRdata_Zenith == 1 )
            hold on;
            yVarrr = linspace(min(OMNI_output(:,7)),max(OMNI_output(:,7)),10);
            plot( repmat(ISR_time_min_Zenith,size(yVarrr)) , yVarrr,'LineWidth',1.75,'Color','r'); %plot red lines showing ISR data time
            hold on;
            plot( repmat(ISR_time_max_Zenith,size(yVarrr)) , yVarrr,'LineWidth',1.75,'Color','r'); %plot red lines showing ISR data time
            hold off;
        end
        ylabel('Psw (nPa)'); %sped up significantly when merged above
        set(ax3,'xtick',Xaxisvar);
        set(ax3,'xticklabels',[])
        axis([min(Xaxisvar) max(Xaxisvar) min(OMNI_output(:,7)) max(OMNI_output(:,7))])
        set(gca,'fontweight',font_Weight, 'FontSize',font_Size);
        grid on
        ax4 = subplot(6,1,4);
        plot( OMNI_time_hr,OMNI_output(:,8) );
        if( FLG_ISRdata_Zenith == 1 )
            hold on;
            yVarrr = linspace(min(OMNI_output(:,8)),max(OMNI_output(:,8)),10);
            plot( repmat(ISR_time_min_Zenith,size(yVarrr)) , yVarrr,'LineWidth',1.75,'Color','r'); %plot red lines showing ISR data time
            hold on;
            plot( repmat(ISR_time_max_Zenith,size(yVarrr)) , yVarrr,'LineWidth',1.75,'Color','r'); %plot red lines showing ISR data time
            hold off;
        end
        ylabel('AE (nT)'); %sped up significantly when merged above
        set(ax4,'xtick',Xaxisvar);
        set(ax4,'xticklabels',[])
        axis([min(Xaxisvar) max(Xaxisvar) min(OMNI_output(:,8)) max(OMNI_output(:,8))])
        set(gca,'fontweight',font_Weight, 'FontSize',font_Size);
        grid on
        ax5 = subplot(6,1,5);
        plot( OMNI_time_hr,OMNI_output(:,9) );
        if( FLG_ISRdata_Zenith == 1 )
            hold on;
            yVarrr = linspace(min(OMNI_output(:,9)),max(OMNI_output(:,9)),10);
            plot( repmat(ISR_time_min_Zenith,size(yVarrr)) , yVarrr,'LineWidth',1.75,'Color','r'); %plot red lines showing ISR data time
            hold on;
            plot( repmat(ISR_time_max_Zenith,size(yVarrr)) , yVarrr,'LineWidth',1.75,'Color','r'); %plot red lines showing ISR data time
            hold off;
        end
        ylabel('SYM/H (nT)'); %sped up significantly when merged above
        set(ax5,'xtick',Xaxisvar);
        set(ax5,'xticklabels',[])
        axis([min(Xaxisvar) max(Xaxisvar) min(OMNI_output(:,9)) max(OMNI_output(:,9))])
        set(gca,'fontweight',font_Weight, 'FontSize',font_Size);
        grid on
        ax6 = subplot(6,1,6); %JEDI BIT
        plot(JEDI_time_full,JEDI_jouleHeating_full);
        if( FLG_ISRdata_Zenith == 1 )
            hold on;
            yVarrr = linspace(min(JEDI_jouleHeating_full),max(JEDI_jouleHeating_full),10);
            plot( repmat(ISR_time_min_Zenith,size(yVarrr)) , yVarrr,'LineWidth',1.75,'Color','r'); %plot red lines showing ISR data time
            hold on;
            plot( repmat(ISR_time_max_Zenith,size(yVarrr)) , yVarrr,'LineWidth',1.75,'Color','r'); %plot red lines showing ISR data time
            hold off;
        end
        ylabel('JouleHeat'); %sped up significantly when merged above
        xlabel('Time (hr)'); %old: ,'fontweight','bold','FontSize',12 for all
        set(ax6,'xtick',Xaxisvar);
        axis([min(Xaxisvar) max(Xaxisvar) min(JEDI_jouleHeating_full) max(JEDI_jouleHeating_full)])
        set(gca,'fontweight',font_Weight, 'FontSize',font_Size8);
        grid on
        
        %OMNI & JEDI PLOT just AE and Joule Heating HERE
        %prep figure #set(gca,'xticklabel',[])
        % figure('units','normalized','outerposition',[0 0 1 1]); %maximizes figure window (caused graphical issues on Windows 10 + AMD GPU or Intel GPU)
        figure(figNum);
        figNum = figNum+1;
        Xaxisvar = OMNI_time_hr_axis_min:4:OMNI_time_hr_axis_max;
        ax1 = subplot(3,1,1);
        plot( OMNI_time_hr,OMNI_output(:,5) );
        if( FLG_ISRdata_Zenith == 1 )
            hold on;
            yVarrr = linspace(min(OMNI_output(:,5)),max(OMNI_output(:,5)),10);
            plot( repmat(ISR_time_min_Zenith,size(yVarrr)) , yVarrr,'LineWidth',1.75,'Color','r'); %plot red lines showing ISR data time
            hold on;
            plot( repmat(ISR_time_max_Zenith,size(yVarrr)) , yVarrr,'LineWidth',1.75,'Color','r'); %plot red lines showing ISR data time
            hold off;
        end
        ylabel('Bz GSE (nT)'); %sped up significantly when merged above
        title(['OMNI (1 min resolution) & JEDI (',num2str(round(JEDI_time_Delta*60)),' min res.) Data for ',num2str(OMNI_dateRange_Full_Month(1,2)),'/',num2str(OMNI_dateRange_Full_Month(1,3)),...
            '/',num2str(OMNI_dateRange_Full_Month(1,1)),' to ',num2str(OMNI_dateRange_Full_Month(end,2)),...
            '/',num2str(OMNI_dateRange_Full_Month(end,3)),'/',num2str(OMNI_dateRange_Full_Month(end,1)),...
            ' (M/D/Y) with 0 Hr on ',dateRange_zeroHr_YrMonDay_string]);
        set(ax1,'xtick',Xaxisvar);
        set(ax1,'xticklabels',[]);
        axis([min(Xaxisvar) max(Xaxisvar) min(OMNI_output(:,5)) max(OMNI_output(:,5))])
        set(gca,'fontweight',font_Weight, 'FontSize',font_Size);
        grid on
        ax2 = subplot(3,1,2);
        plot( OMNI_time_hr,OMNI_output(:,8) );
        if( FLG_ISRdata_Zenith == 1 )
            hold on;
            yVarrr = linspace(min(OMNI_output(:,8)),max(OMNI_output(:,8)),10);
            plot( repmat(ISR_time_min_Zenith,size(yVarrr)) , yVarrr,'LineWidth',1.75,'Color','r'); %plot red lines showing ISR data time
            hold on;
            plot( repmat(ISR_time_max_Zenith,size(yVarrr)) , yVarrr,'LineWidth',1.75,'Color','r'); %plot red lines showing ISR data time
            hold off;
        end
        ylabel('AE (nT)'); %sped up significantly when merged above
        set(ax2,'xtick',Xaxisvar);
        set(ax2,'xticklabels',[]);
        axis([min(Xaxisvar) max(Xaxisvar) min(OMNI_output(:,8)) max(OMNI_output(:,8))])
        set(gca,'fontweight',font_Weight, 'FontSize',font_Size);
        grid on
        ax3 = subplot(3,1,3); %JEDI BIT
        plot(JEDI_time_full,JEDI_jouleHeating_full);
        if( FLG_ISRdata_Zenith == 1 )
            hold on;
            yVarrr = linspace(min(JEDI_jouleHeating_full),max(JEDI_jouleHeating_full),10);
            plot( repmat(ISR_time_min_Zenith,size(yVarrr)) , yVarrr,'LineWidth',1.75,'Color','r'); %plot red lines showing ISR data time
            hold on;
            plot( repmat(ISR_time_max_Zenith,size(yVarrr)) , yVarrr,'LineWidth',1.75,'Color','r'); %plot red lines showing ISR data time
            hold off;
        end
        ylabel('Joule Heating (?)'); %sped up significantly when merged above
        xlabel('Time (hr)'); %old: ,'fontweight','bold','FontSize',12 for all
        set(ax3,'xtick',Xaxisvar);
        axis([min(Xaxisvar) max(Xaxisvar) min(JEDI_jouleHeating_full) max(JEDI_jouleHeating_full)])
        set(gca,'fontweight',font_Weight, 'FontSize',font_Size);
        grid on
        
        
        %OMNI & JEDI PLOT just AE and Joule Heating HERE THAT ARE HIGH-PASSED
        %prep figure #set(gca,'xticklabel',[])
        % figure('units','normalized','outerposition',[0 0 1 1]); %maximizes figure window (caused graphical issues on Windows 10 + AMD GPU or Intel GPU)
        figure(figNum);
        figNum = figNum+1;
        Xaxisvar = OMNI_time_hr_axis_min:4:OMNI_time_hr_axis_max;
        ax1 = subplot(3,1,1);
        plot( OMNI_time_hr,OMNI_output_Bz_BP );
        if( FLG_ISRdata_Zenith == 1 )
            hold on;
            yVarrr = linspace(min(OMNI_output_Bz_BP),max(OMNI_output_Bz_BP),10);
            plot( repmat(ISR_time_min_Zenith,size(yVarrr)) , yVarrr,'LineWidth',1.75,'Color','r'); %plot red lines showing ISR data time
            hold on;
            plot( repmat(ISR_time_max_Zenith,size(yVarrr)) , yVarrr,'LineWidth',1.75,'Color','r'); %plot red lines showing ISR data time
            hold off;
        end
        ylabel('Bz GSE HP (nT)'); %sped up significantly when merged above
        title(['OMNI (1 min resolution) & JEDI (',num2str(round(JEDI_time_Delta*60)),' min res.) ',num2str(round(1/JEDI_highpass_Freq,2)),' Hr High-Passed Data for ',num2str(OMNI_dateRange_Full_Month(1,2)),'/',num2str(OMNI_dateRange_Full_Month(1,3)),...
            '/',num2str(OMNI_dateRange_Full_Month(1,1)),' to ',num2str(OMNI_dateRange_Full_Month(end,2)),...
            '/',num2str(OMNI_dateRange_Full_Month(end,3)),'/',num2str(OMNI_dateRange_Full_Month(end,1)),...
            ' (M/D/Y) with 0 Hr on ',dateRange_zeroHr_YrMonDay_string]);
        set(ax1,'xtick',Xaxisvar);
        set(ax1,'xticklabels',[]);
        axis([min(Xaxisvar) max(Xaxisvar) min(OMNI_output_Bz_BP) max(OMNI_output_Bz_BP)])
        set(gca,'fontweight',font_Weight, 'FontSize',font_Size);
        grid on
        ax2 = subplot(3,1,2);
        plot( OMNI_time_hr,OMNI_output_AE_BP );
        if( FLG_ISRdata_Zenith == 1 )
            hold on;
            yVarrr = linspace(min(OMNI_output_AE_BP),max(OMNI_output_AE_BP),10);
            plot( repmat(ISR_time_min_Zenith,size(yVarrr)) , yVarrr,'LineWidth',1.75,'Color','r'); %plot red lines showing ISR data time
            hold on;
            plot( repmat(ISR_time_max_Zenith,size(yVarrr)) , yVarrr,'LineWidth',1.75,'Color','r'); %plot red lines showing ISR data time
            hold off;
        end
        ylabel('AE HP (nT)'); %sped up significantly when merged above
        set(ax2,'xtick',Xaxisvar);
        set(ax2,'xticklabels',[]);
        axis([min(Xaxisvar) max(Xaxisvar) min(OMNI_output_AE_BP) max(OMNI_output_AE_BP)])
        set(gca,'fontweight',font_Weight, 'FontSize',font_Size);
        grid on
        ax3 = subplot(3,1,3); %JEDI BIT
        plot(JEDI_time_full,JEDI_jouleHeating_full_bp);
        if( FLG_ISRdata_Zenith == 1 )
            hold on;
            yVarrr = linspace(min(JEDI_jouleHeating_full_bp),max(JEDI_jouleHeating_full_bp),10);
            plot( repmat(ISR_time_min_Zenith,size(yVarrr)) , yVarrr,'LineWidth',1.75,'Color','r'); %plot red lines showing ISR data time
            hold on;
            plot( repmat(ISR_time_max_Zenith,size(yVarrr)) , yVarrr,'LineWidth',1.75,'Color','r'); %plot red lines showing ISR data time
            hold off;
        end
        ylabel('Joule Heating HP (?)'); %sped up significantly when merged above
        xlabel('Time (hr)'); %old: ,'fontweight','bold','FontSize',12 for all
        set(ax3,'xtick',Xaxisvar);
        axis([min(Xaxisvar) max(Xaxisvar) min(JEDI_jouleHeating_full_bp) max(JEDI_jouleHeating_full_bp)])
        set(gca,'fontweight',font_Weight, 'FontSize',font_Size);
        grid on

        
%         Xaxisvar = [ (round(min(JEDI_time))-mod(round((min(JEDI_time))),2)):2:...
%             (round(max(JEDI_time))-mod(round((max(JEDI_time))),2)) ];
%         set(gca,'xtick',Xaxisvar);
%         xlim([min(JEDI_time) max(JEDI_time)]);
    
        %****************FREQ OF OMNI & JEDI*********************
        Ross_Lombscargle_optimized([ OMNI_time_hr , OMNI_output_Bz_BP ]); 
        fid = fopen('power_freq.txt', 'r');
            f_nm = textscan(fid, '%f %f %f %f %f', 'headerlines', 15);
            P_Zenith_SNR=f_nm{:,2};
            F_Zenith_SNR=f_nm{:,1};
            T_Zenith_SNR=f_nm{:,5}*60; % in minutes
            gs_Zenith_SNR=f_nm{:,4};
            max_Zenith_SNR=T_Zenith_SNR(P_Zenith_SNR==max(P_Zenith_SNR));
            K=length(P_Zenith_SNR);
            gf_Zenith_SNR=f_nm{:,3};
        fclose('all');

        %prep figure #
        % figure('units','normalized','outerposition',[0 0 1 1]); %maximizes figure window
        figure(figNum);
        figNum = figNum+1;
        subplot(3,1,1);
        plot(T_Zenith_SNR,P_Zenith_SNR,'b');
        xlim([0 plot_Freq_Lim]);
        hold on;
        % plot(T_SNR_08apr08_10_or,P_SNR_08apr08_10_or,'r');
        % hold on;
        plot(T_Zenith_SNR,gf_Zenith_SNR,'LineWidth',1.5,'color',[0.5 0.5 0.5]);
        xlim([0 plot_Freq_Lim]);  
        xlabel('Periods (min)');
        ylabel('Normalized Power');
        string_Title = ['OMNI Bz GSE High-pass Cutoff Period of ',num2str(round(1/bs)),' Hr Lomb-Scargle Periodogram'];
        if( FLG_geomagneticCoords == 1)
            string_Title = [string_Title,' [Mag Coord Aligned]']; %tack on mag coordinate warning
        end
        title(string_Title);
        set(gca,'fontweight',font_Weight, 'FontSize',font_Size);
        set(gca,'xtick',Xaxisvar_SCARGLE);
        hold off;
        
        Ross_Lombscargle_optimized([ OMNI_time_hr , OMNI_output_AE_BP ]); 
        fid = fopen('power_freq.txt', 'r');
            f_nm = textscan(fid, '%f %f %f %f %f', 'headerlines', 15);
            P_Zenith_SNR=f_nm{:,2};
            F_Zenith_SNR=f_nm{:,1};
            T_Zenith_SNR=f_nm{:,5}*60; % in minutes
            gs_Zenith_SNR=f_nm{:,4};
            max_Zenith_SNR=T_Zenith_SNR(P_Zenith_SNR==max(P_Zenith_SNR));
            K=length(P_Zenith_SNR);
            gf_Zenith_SNR=f_nm{:,3};
        fclose('all');
        
        subplot(3,1,2); %put down the AE now
        plot(T_Zenith_SNR,P_Zenith_SNR,'b');
        xlim([0 plot_Freq_Lim]);
        hold on;
        % plot(T_SNR_08apr08_10_or,P_SNR_08apr08_10_or,'r');
        % hold on;
        plot(T_Zenith_SNR,gf_Zenith_SNR,'LineWidth',1.5,'color',[0.5 0.5 0.5]);
        xlim([0 plot_Freq_Lim]);  
        xlabel('Periods (min)');
        ylabel('Normalized Power');
        string_Title = ['OMNI AE Index High-pass Cutoff Period of ',num2str(round(1/bs)),' Hr Lomb-Scargle Periodogram'];
        if( FLG_geomagneticCoords == 1)
            string_Title = [string_Title,' [Mag Coord Aligned]']; %tack on mag coordinate warning
        end
        title(string_Title);
        set(gca,'fontweight',font_Weight, 'FontSize',font_Size);
        set(gca,'xtick',Xaxisvar_SCARGLE); 
        hold off; 
        
        Ross_Lombscargle_optimized([ JEDI_time_full , JEDI_jouleHeating_full_bp ]); 
        fid = fopen('power_freq.txt', 'r');
            f_nm = textscan(fid, '%f %f %f %f %f', 'headerlines', 15);
            P_Zenith_SNR=f_nm{:,2};
            F_Zenith_SNR=f_nm{:,1};
            T_Zenith_SNR=f_nm{:,5}*60; % in minutes
            gs_Zenith_SNR=f_nm{:,4};
            max_Zenith_SNR=T_Zenith_SNR(P_Zenith_SNR==max(P_Zenith_SNR));
            K=length(P_Zenith_SNR);
            gf_Zenith_SNR=f_nm{:,3};
        fclose('all');

        subplot(3,1,3);
        plot(T_Zenith_SNR,P_Zenith_SNR,'b');
        xlim([0 plot_Freq_Lim]);
        hold on;
        % plot(T_SNR_08apr08_10_or,P_SNR_08apr08_10_or,'r');
        % hold on;
        plot(T_Zenith_SNR,gf_Zenith_SNR,'LineWidth',1.5,'color',[0.5 0.5 0.5]);
        xlim([0 plot_Freq_Lim]);  
        xlabel('Periods (min)');
        ylabel('Normalized Power');
        string_Title = ['JEDI Joule Heating High-pass Cutoff Period of ',num2str(round(1/JEDI_highpass_Freq)),' Hr Lomb-Scargle Periodogram'];
        if( FLG_geomagneticCoords == 1)
            string_Title = [string_Title,' [Mag Coord Aligned]']; %tack on mag coordinate warning
        end
        title(string_Title);
        set(gca,'fontweight',font_Weight, 'FontSize',font_Size);
        set(gca,'xtick',Xaxisvar_SCARGLE); 
        hold off;        
    end
    
    
    %===============Basic Day/Night Time Plot==============================
    if( FLG_ISRdata_Zenith == 1 )
        Xaxisvar_min = (Zenith_time(1)-round(mean(Zenith_timeOffset)))*24; %UT hr, min time to compare to
        Xaxisvar_max = (Zenith_time(end)-round(mean(Zenith_timeOffset)))*24; %UT hr, max time to compare to

        %FIRST BATTLE: REGION WHERE LAT/LONG IS   
        url = ['http://api.geonames.org/timezone?lat=',num2str(latMillstone),'&lng=',num2str(longMillstone),'&username=razzluhdzuul']; %create link for lat/long
        siteData = xmlread(url); %DL the returned data
        dstOffset = siteData.getElementsByTagName('dstOffset'); %get the daylight savings time offset
        utOffset = siteData.getElementsByTagName('gmtOffset'); %get the UT offset
        dayNite_DSToffset = char(dstOffset.item(0).getTextContent());
        dayNite_UToffset = char(utOffset.item(0).getTextContent());
        dayNite_timeZoneID = siteData.getElementsByTagName('timezoneId'); %get the time zone ID
        dayNite_timeZoneID = char(dayNite_timeZoneID.item(0).getTextContent()); %get a name that is usable

        dateRange_dates = sFUN_dayNumber_to_Date_MULTIPLE(dateRange); %convert day# range to yr/mon/day range with all inbetween filled in
        dateRange_extended = sFUN_dayNumber_to_Date_MULTIPLE(dateRange,1); %convert day# range to yr/day# flushed out
        dayNite_DSTactive = isdst(datetime(dateRange_dates(1,1),dateRange_dates(1,2),dateRange_dates(:,3),'TimeZone',dayNite_timeZoneID)); %get DST active or not for each day
        
        %SECOND STEP: CALC LOCAL SUNRISE/SUNSET TIMES
        %will adjust for DST later - calcs done in UT/GMT
        % based on calc steps in https://www.mathworks.com/examples/matlab/community/21093-estimating-sunrise-and-sunset
        long_corrected = 4*(longMillstone - 15*str2double(dayNite_UToffset)); %calc corrected longitude, for sunrise/sunset time

        dayNite_B = 360*(dateRange_extended(:,2) - 81)/365; %some sort of angle based on days and stuff
        dayNite_EoT_corrected = 9.87*sind(2*dayNite_B) - 7.53*cosd(dayNite_B) - 1.5*sind(dayNite_B); %eq for Time Correction
        dayNite_solar_corrected = long_corrected + dayNite_EoT_corrected; %min, solar time correction - for noon

        dayNite_solar_declination = asind(sind(23.45)*sind(360*(dateRange_extended(:,2) - 81)/365)); %deg, solar declination

        dayNite_sunrise = 12 - acosd(-tand(latMillstone)*tand(dayNite_solar_declination))/15 - dayNite_solar_corrected/60; %hr, sunrise time
        dayNite_sunset = 12 + acosd(-tand(latMillstone)*tand(dayNite_solar_declination))/15 - dayNite_solar_corrected/60; %hr, sunrise time

        %Daylight savings local time fix
        for( i = 1:size(dateRange_dates,1) )
            %Run through the days
            dayNite_DSToffset_delta = str2double(dayNite_DSToffset) - str2double(dayNite_UToffset); %hr, calc the delta

            if( dayNite_DSTactive(i) == 1 ) %daylight savings time is active
                dayNite_sunrise(i) = dayNite_sunrise(i) + dayNite_DSToffset_delta; %apply the offset so it's local time
                dayNite_sunset(i) = dayNite_sunset(i) + dayNite_DSToffset_delta; %apply the offset so it's local time
            else

            end
        end
        
        %THIRD STEP: PREP FOR PLOTTING BY ALIGNING TIMES, MAKING PLOT VARIABLES
        if( dayNite_DSTactive(1) == 1 )
            Xaxisvar_min = Xaxisvar_min + str2double(dayNite_DSToffset); %hr local, UT time of -12 conv. to local
        else
            Xaxisvar_min = Xaxisvar_min + str2double(dayNite_UToffset); %hr local, UT time of -12 conv. to local
        end

        if( dayNite_DSTactive(end) == 1 )
            Xaxisvar_max = Xaxisvar_max + str2double(dayNite_DSToffset); %hr local, UT time of -12 conv. to local
        else
            Xaxisvar_max = Xaxisvar_max + str2double(dayNite_UToffset); %hr local, UT time of -12 conv. to local
        end

        for( i = 1:size(dateRange_dates,1) )
            dayNite_sunrise(i) = dayNite_sunrise(i) + (dateRange_extended(i,2) - round(mean(dateRange(:,2))))*24; %hr, adjust so 0 time is middle of date range
            dayNite_sunset(i) = dayNite_sunset(i) + (dateRange_extended(i,2) - round(mean(dateRange(:,2))))*24; %hr, adjust so 0 time is middle of date range
        end


        Xaxisvar = floor(Xaxisvar_min):2:ceil(Xaxisvar_max); %hrs, axis vars to use

        xTime = sort([dayNite_sunrise;dayNite_sunset])'; %hr, xvar to plot against
        yDayNite = repmat([1,0],[1,length(xTime)/2]); %day == 1, night == 0
        yDayNite( xTime < Xaxisvar_min) = []; %delete less than the min time req in the day/nite vector
        xTime( xTime < Xaxisvar_min) = []; %delete less than the min time req
        yDayNite( xTime > Xaxisvar_max) = []; %delete more than the max time req in the day/nite vector
        xTime( xTime > Xaxisvar_max) = []; %delete more than the max time req
        xTime = [Xaxisvar_min,xTime,Xaxisvar_max]; %tack on the start and ending times
        yDayNite = [1,yDayNite,0]; %tack on a day/night to correspond to the start and ending times which are probably day at the start and night at the end though idk
        %I think it works logically thru magic and the fact days have sunrise and sunset you know shh it cool no it doesn't work perfectly but it'll do here
        %THIS IS DEF NOT 100% SOLID ALL THE TIME but works for right now watch this
        xTime = repelem(xTime,2); %for plotting purposes (so it looks like a flat line)
        yDayNite = repelem(yDayNite,2); %for plotting purposes (so it looks like a flat line)
        yDayNite = circshift(yDayNite,1); %for plotting purposes, need to make the 1st thing be the same as the 2nd thing
        yDayNite(1) = yDayNite(2); %fix this up
        
        %FOURTH STEP: ACTUALLY PLOTTING
        dayNite_timeZoneID = strrep(dayNite_timeZoneID,'_',' ');
        %prep figure #
        % figure('units','normalized','outerposition',[0 0 1 1]); %maximizes figure window
        figure(figNum);
        figNum = figNum+1;
        plot(xTime,yDayNite,'k','LineWidth',4);
        hold on;
    %     title('Local Day/Night at Millstone Hill');
        set(gca,'yticklabel',[])
        if( dayNite_DSTactive(1) == 1 )
            strang = [dayNite_DSToffset,' (Daylight Savings) in the ',dayNite_timeZoneID,' Time Zone'];
        else
             strang = [dayNite_UToffset,' in the ',dayNite_timeZoneID,' Time Zone'];
        end
        if( yDayNite(end) == yDayNite(end-1) )
            %hacky offset so day or night isn't printed off the plot (or skipped)
            offset = 3;
        else
            offset = 1;
        end
        for( i = 2:2:length(yDayNite)-offset )
            if( yDayNite(i) == 1 )
                text( (xTime(i+1)-xTime(i))/2+xTime(i)-1.3,1-.15,'Day',...
                    'Color','k','fontweight','bold','FontSize',40); %plots text that shows day or night
                hold on;
            else
                text( (xTime(i+1)-xTime(i))/2+xTime(i)-1.9,.25,'Night',...
                    'Color','k','fontweight','bold','FontSize',40); %plots text that shows day or night
                hold on;
            end
        end

        xlabel(['Local Time, ',strang,' (hr)']);
        xlim([Xaxisvar_min , Xaxisvar_max]); %limit to time desired
        set(gca,'xtick',Xaxisvar);
        set(gca, 'Ticklength', [0 0])
        set(gca,'fontweight',font_Weight, 'FontSize',font_Size);
    end %end of Day/Night IF Zenith data is enabled

end


%% Gather Data at a Point
if(FLG_gatherDataAtPoint == 1)

%Select the point that is the main one to compare to
gg_ZenithSite = 1:length(latPoints);
if( length(unique(latPoints)) == length(latPoints) )
    gg_ZenithSite = gg_ZenithSite(min(abs(latPoints - latMillstone)) == abs(latPoints - latMillstone)); %get which lat is closer to Millstone
else
    gg_ZenithSite = gg_ZenithSite(min(abs(longPoints - longMillstone)) == abs(longPoints - longMillstone)); %get which long is closer to Millstone
end

%This is the MISA site in event it is needed
gg_MISASite = 1:length(latPoints);
if( length(unique(latPoints)) == length(latPoints) )
    gg_MISASite = gg_MISASite(min(abs(latPoints - latMillstone)) == abs(latPoints - latMillstone)); %get which lat is closer to Millstone
else
    gg_MISASite = gg_MISASite(min(abs(longPoints - longMillstone)) == abs(longPoints - longMillstone)); %get which long is closer to Millstone
end

tic
pointRadiusAngular = (pointRadius/Re)*180/pi; %deg, (based on s=R*theta) angular radius around the point allowed

k = zeros(size(pplat)); %preallocate this
tempRadius = zeros(size(latPoints,2),length(pplat)); %record the radius calc'd
for(gg = 1:size(latPoints,2) ) %gather data for multiple points
    tempRadius(gg,:) = sqrt((pplat - latPoints(gg)).^2 + (pplong - longPoints(gg)).^2); %record the radius calc'd b/c why not save calc's
    k = k | (pointRadiusAngular > tempRadius(gg,:));
end
pplat_relev = pplat(k); %use results of points within the needed areas to get relevant data
pplong_relev = pplong(k);
sTEC_relev = sTEC(k);
time_relev = time(k);
tempRadius = tempRadius(:,k); %this is still useful promise :) despite covering multiple points

sTEC_holder = cell(2,1); %preallocate this 
parfor(gg = 1:size(latPoints,2) ) %gather data for multiple points
    
    
    j = time_relev == timeUnique(1); %find the time that connect to the time we are at

    k = pointRadiusAngular > tempRadius(gg,j); %get data within the radius
%     k = pointRadiusAngular > sqrt((pplat_relev(j) - latPoints(gg)).^2 + (pplong_relev(j) - longPoints(gg)).^2); %get data within the radius
    %WAS: pointRadiusAngular > sqrt((pplat(j) - latPoints(gg)).^2 + (pplong(j) - longPoints(gg)).^2)
    %to check for latitude points within the radius - calc'd it all at once before

    tempsTEC = sTEC_relev(j);
    tempLat = pplat_relev(j);
    tempLong = pplong_relev(j);
    sTEC_holder{gg,:,:} = [tempsTEC(k);tempLat(k);tempLong(k);repmat(timeUnique(1),1,size(tempLong(k),2))]';
    %sTEC is in 1
    %lat is in 2
    %long is in 3
    %timeUnique index is in 4


    for( t = 2:size(timeUnique,2) ) %loops through time

        j = time_relev == timeUnique(t); %find the time that connect to the time we are at

        k = pointRadiusAngular > tempRadius(gg,j); %get data within the radius
%         k = pointRadiusAngular > sqrt((pplat_relev(j) - latPoints(gg)).^2 + (pplong_relev(j) - longPoints(gg)).^2); %get data within the radius
        %WAS: pointRadiusAngular > sqrt((pplat(j) - latPoints(gg)).^2 + (pplong(j) - longPoints(gg)).^2)
        %to check for latitude points within the radius - calc'd it all at once before

        tempsTEC = sTEC_relev(j); %preps my temporary variables
        tempLat = pplat_relev(j);
        tempLong = pplong_relev(j);
        
        sTEC_holder{gg,:,:} = [sTEC_holder{gg,:,:};[tempsTEC(k);tempLat(k);tempLong(k);repmat(timeUnique(t),1,size(tempLong(k),2))]' ]; %tacs on new data

    end
    
end

clear pplat_relev pplong_relev sTEC_relev time_relev tempRadius

%OLD LONG WAY for reference, top way confirmed to match perfectly
% finAnnounce_percentToUpdateAt = 10; % every % to update info at
% finAnnounce_div = round(100/finAnnounce_percentToUpdateAt); %calc divisor to use
% for(gg = 1:size(latPoints,2) ) %gather data for multiple points
% 
%     j = time == timeUnique(1); %find the time that connect to the time we are at
% 
%     k = pointRadiusAngular > sqrt((pplat(j) - latPoints(gg)).^2 + (pplong(j) - longPoints(gg)).^2);
% 
%     tempsTEC = sTEC(j);
%     sTEC_holderTEMP(:,1) = tempsTEC(k); %sTEC is in 1
% 
%     tempLat = pplat(j);
%     sTEC_holderTEMP(:,2) = tempLat(k); %lat is in 2
% 
%     tempLong = pplong(j);
%     sTEC_holderTEMP(:,3) = tempLong(k); %long is in 3
% 
%     sTEC_holderTEMP(:,4) = repmat(timeUnique(1),1,size(tempLong(k),2)); %timeUnique index is in 4
% 
% 
%     parfor( t = 2:size(timeUnique,2) ) %loops through time
% 
%         j = time == timeUnique(t); %find the time that connect to the time we are at
% 
%         k = pointRadiusAngular > sqrt((pplat(j) - latPoints(gg)).^2 + (pplong(j) - longPoints(gg)).^2);
% 
%         tempsTEC = sTEC(j); %preps my temporary variables
%         tempLat = pplat(j);
%         tempLong = pplong(j);
% 
%         sTEC_holderTEMP = [sTEC_holderTEMP;[tempsTEC(k);tempLat(k);tempLong(k);repmat(timeUnique(t),1,size(tempLong(k),2))]' ]; %tacs on new data
%         
%     end
%     
%     if( mod(gg,round(length(latPoints)/finAnnounce_div)) == 0 )
%         fprintf('PROGRESS: Averaging TEC point zones %d/%d\t%%Complete: %.1f%%\tRuntime: %.3f sec\n',gg,length(latPoints),100*gg/length(latPoints),toc); %announce job, % complete, and time so far
%     end
%     
%     sTEC_holder{gg,:,:} = sTEC_holderTEMP;
%     clear sTEC_holderTEMP
%     
%     %Generic progress reporter
%     
% end
tocSave = toc;
fprintf('GATHER DATA AT A POINT: Time to complete Averaging %d TEC Point Radius (%.2f km) Zones:\t%.2f sec / %.2f min\n\n',length(latPoints),pointRadius,tocSave,tocSave/60); %report


%THIS IS TO VIEW DATA AVG RANGE BEING TAKEN
coastLines = load('coast');%loads coast line into structure, vars are lat and long - can interfere with previous code so hidden in structure
coastLines_lat = coastLines.lat; %breaks this data out
coastLines_long = coastLines.long;
clear coastLines
if( FLG_geomagneticCoords == 1 ) %convert data to geomagnetic coords if desired
    [coastLines_lat,coastLines_long] = sFUN_geoToGeomag(dateOfData,magF,coastLines_lat,coastLines_long); %convert from geographic to geomagnetic coordinates
end %doing as we read data bit by bit because full math conversion will cause mem to run out with my 16 gigs omg
jk = find(coastLines_lat < min(plotLatRange) ); %find less than min, remove
coastLines_lat(jk) = []; %remove
coastLines_long(jk) = []; %remove
jk = find(coastLines_lat > max(plotLatRange) ); %find more than max, remove
coastLines_lat(jk) = []; %remove
coastLines_long(jk) = []; %remove
jk = find(coastLines_long < min(plotLongRange) ); %find less than min, remove
coastLines_lat(jk) = []; %remove
coastLines_long(jk) = []; %remove
jk = find(coastLines_long > max(plotLongRange) ); %find more than max, remove
coastLines_lat(jk) = []; %remove
coastLines_long(jk) = []; %remove
coastLines_long( coastLines_long == 180 ) = 180-.01; %fix some weird thing where 180 flips to -180

%Corral the data to the right place    
k = find( time == timeUnique(1)); %gets during a time period
%prep figure #
% figure('units','normalized','outerposition',[0 0 1 1]); %maximizes figure window
figure(figNum);
figNum = figNum+1;
scatter(pplong(k),pplat(k),plot_Scatter_Point_Size,sTEC(k),'.');
caxis([-gif_sTEC_cap, gif_sTEC_cap]);
%         colormap(flipud(colormap('hot'))) %heat map where 0 (low) is white
colormap('jet') %heat map
h = colorbar; %shows the color bar
ylabel(h, 'delta-sTEC (TECU)','fontweight',font_Weight, 'FontSize',font_Size); %labels color bar
hold on;
%Now drawing line of interest
plot(longMillstone,latMillstone,'*','Color',[255/255,36/255,0/255],'MarkerSize',25,'LineWidth',1.9); %plots a point with a red big *
hold on;
%Plot lines of contients
geoshow(coastLines_lat,coastLines_long,'Color','k') %add continental outlines in black  
xlabel('Longitude (arcdeg)'); %old: ,'fontweight','bold','FontSize',12 for all
ylabel('Latitude (arcdeg)'); %sped up significantly when merged above
string_Title = ['Radius of ',num2str(pointRadius),' km around Millstone Hill ISR']; %create mecha title
if( FLG_geomagneticCoords == 1)
    string_Title = [string_Title,' [Mag Coord Aligned]']; %tack on mag coordinate warning
end
title(string_Title);
axis([(longMillstone-(pointRadius/Re)*180/pi*17), (longMillstone+pointRadius/Re*180/pi*17), (latMillstone-pointRadius/Re*180/pi*10), (latMillstone+pointRadius/Re*180/pi*10)])
% xticks( min(plotLongRange):plotLongRange_autoTick:max(plotLongRange) ); %creates x ticks automagically
% yticks( min(plotLatRange):plotLatRange_autoTick:max(plotLatRange) ); %creates y ticks automagically
hold on
h = plot((pointRadius/Re)*180/pi*cos(0:pi/50:2*pi) + longMillstone, (pointRadius/Re)*180/pi*sin(0:pi/50:2*pi) + latMillstone,'LineWidth',4,'Color',[148/255,0/255,211/255]);
hold off
set(gca,'fontweight',font_Weight, 'FontSize',font_Size);

end %*********END OF FLG_gatherDataAtPoint*********


%% Gather Data at a Point - Do a little Median Filter & Combine Data
if( FLG_gatherDataAtPoint_Filter == 1 )  
    
sTEC_combined = zeros(size(latPoints,2),size(timeUnique,2)); %preallocate (this has missing time gaps - most likely)
sTEC_combinedInterped = zeros(size(latPoints,2),size(timeUnique,2)); %preallocate (this has no missing time gaps through interpolation)
for(gg = 1:size(latPoints,2))
    
    for( i = 1:length(timeUnique) ) %combine data at same time into one value at one time

        jp = find( timeUnique(i) == sTEC_holder{gg,1}(:,4) ); %find the data points at the same time
        
        temp = sTEC_holder{gg,1}(jp,1); %save data into a temp variable
        
        tempMean = mean( sTEC_holder{gg,1}(jp,1) ); %get the mean
        
        tempVar = var( sTEC_holder{gg,1}(jp,1) ); %get variance
        
        %loop to prevent too much data being nixed
        while( (length(temp) - sum(((tempMean+dataReject*tempVar) > temp) & ((tempMean-dataReject*tempVar) < temp) ) )*100/length(temp) > dataRejectLimit )
            dataReject = dataReject*1.05; %increase the data rejection ratio to include more data
            if( dataReject > dataRejectMax )
                dataRejectLimit = 100; %if this goes on for a bit then it is halted with this manuever
            end
%             (length(temp) - sum(((tempMean+dataReject*tempVar) > temp) & ((tempMean-dataReject*tempVar) < temp) ) )*100/length(temp)
%             dataReject
%             dataRejectLimit
        end
        
        if( dataReject > dataRejectMax ) %if this occured leave the data, it is too sparse or too varied to deal with effectively
            
        else
            temp = temp( ((tempMean+dataReject*tempVar) > temp) & ((tempMean-dataReject*tempVar) < temp) ); %delete some extraneous data
            %this is normal operation
        end
        
        dataReject = dataRejectOrig; %reset dataReject
        dataRejectLimit = dataRejectLimitOrig; %reset dataRejectLimit

        sTEC_combined(gg,i) = mean( temp ); %TECU, average data at same time into one value
    end
    
    jk = find(isnan(sTEC_combined(gg,:)) == 1); %find NaNs existing
    jl = find(isnan(sTEC_combined(gg,:)) ~= 1); %find NaNs not existing
    sTEC_combinedInterped(gg,:) = sTEC_combined(gg,:); %save the O.G. data
    sTEC_combinedInterped(gg,jk) = interp1(jl,sTEC_combined(gg,jl),jk,'linear','extrap'); %fill in the gaps with the spline fill in
    fprintf('At lat %.3f deg/long %.3f deg %%%.3f of data was NaN and interpolated over.\n',latPoints(gg),longPoints(gg),length(jk)/length(sTEC_combined(gg,:))*100); %report % data filled in
    
end
clear temp tempMean tempVar
fprintf('\n'); %add a space for good measure

timeUnique_combined_noNaN = timeUnique( ~isnan(sTEC_combined(gg_ZenithSite,:)) ); %get some time that's lost the NaN times
sTEC_combined( : , isnan(sTEC_combined(gg_ZenithSite,:)) ) = []; %delete the NaNs

end %*********END OF FLG_gatherDataAtPoint_Filter*********


%% Gather Data at a Point - Plot the Data Up
if( FLG_gatherDataAtPoint_Filter_Plot == 1 )

%prep figure #
% figure('units','normalized','outerposition',[0 0 1 1]); %maximizes figure window (caused graphical issues on Windows 10 + AMD GPU or Intel GPU)
figure(figNum);
figNum = figNum+1;

for(gg = 1:size(latPoints,2) ) %gather data for multiple points
    
    subplot(size(latPoints,2),1,gg)
    h = plot( (timeUnique_combined_noNaN-floor(mean(time))).*24 ,sTEC_combined(gg,:) );
    if(gg == 1)
        set(h(1),'Color','g')
    elseif(gg == 2)
        set(h(1),'Color','b')
    elseif(gg == 3)
        set(h(1),'Color','r')
    elseif(gg == 4)
        set(h(1),'Color','c')
    elseif(gg == 5)
        set(h(1),'Color','k')
    end
    xlabel('Time (hr)'); %old: ,'fontweight','bold','FontSize',12 for all
    ylabel('Delta Slant TEC (TECU)'); %sped up significantly when merged above
    title(['Delta sTEC at Lat ',num2str(latPoints(:,gg)),' deg/ Long ',num2str(longPoints(:,gg)),'deg on Day ',num2str(floor(mean(time))),', 2013']);
    axis([-inf,inf,-1.5,1.5])
end


% for(gg = 1:size(latPoints,2) ) %gather data for multiple points
% 
% 
%     plot( (sTEC_holder{gg,1}(:,4)-floor(mean(time))).*24 ,sTEC_holder{gg,1}(:,1),'.');
%     xlabel('Time (hr)'); %old: ,'fontweight','bold','FontSize',12 for all
%     ylabel('Delta Slant TEC (TECU)'); %sped up significantly when merged above
%     title(['Delta sTEC at ',num2str(latPoints(gg)),' deg/',num2str(longPoints(gg)),'deg on Day ',num2str(floor(mean(time))),', 2013']);
%     axis([0,24,min(sTEC_holder{gg,1}(:,1))-3,max(sTEC_holder{gg,1}(:,1))+3])
% end

%prep figure #
% figure('units','normalized','outerposition',[0 0 1 1]); %maximizes figure window
figure(figNum);
figNum = figNum+1;
for(gg = 1:size(latPoints,2) ) %gather data for multiple points
    
    h = plot( (timeUnique_combined_noNaN-floor(mean(time))).*24 ,sTEC_combined(gg,:) ); %plot all at once on this
    hold on;
    if(gg == 1)
        set(h(1),'Color','g')
    elseif(gg == 2)
        set(h(1),'Color','b')
    elseif(gg == 3)
        set(h(1),'Color','r')
    elseif(gg == 4)
        set(h(1),'Color','c')
    elseif(gg == 5)
        set(h(1),'Color','k')
    end
    
    legendString{gg,1} = ['Lat ',num2str(latPoints(:,gg)),' / Long ',num2str(longPoints(:,gg))]; %create legend string

end
xlabel('Time (hr)'); %old: ,'fontweight','bold','FontSize',12 for all
ylabel('Delta Slant TEC (TECU)'); %sped up significantly when merged above
title(['Delta sTEC at +/-',num2str(latSpacing),' deg from ',num2str(latMillstone),' deg/',num2str(longMillstone),'deg on Day ',num2str(floor(mean(time))),', 2013']);
axis([-inf,inf,-1.5,1.5])
legend(legendString);

end %*********END OF FLG_gatherDataAtPoint_Filter_Plot*********


%% Gather Data at a Point - IT'S SCARGLIN TIM
if( FLG_gatherDataAtPoint_Filter_Scargle == 1 )

tic
for(gg = 1:size(latPoints,2) ) %gather data for multiple points
%     Ross_Lombscargle([(timeUnique'-floor(mean(time))).*24 ,sTEC_combined(gg,:)']); 
    Ross_Lombscargle_optimized([(timeUnique_combined_noNaN'-floor(mean(time))).*24 ,sTEC_combined(gg,:)']); 
    fid = fopen('power_freq.txt', 'r');
        f_nm = textscan(fid, '%f %f %f %f %f', 'headerlines', 15);
        P_Zenith_SNR=f_nm{:,2};
        F_Zenith_SNR=f_nm{:,1};
        T_Zenith_SNR=f_nm{:,5}*60; % in minutes
        gs_Zenith_SNR=f_nm{:,4};
        max_Zenith_SNR=T_Zenith_SNR(P_Zenith_SNR==max(P_Zenith_SNR));
        K=length(P_Zenith_SNR);
        gf_Zenith_SNR=f_nm{:,3};
    fclose('all');

    %prep figure #
%     figure('units','normalized','outerposition',[0 0 1 1]); %maximizes figure window
    figure(figNum);
    figNum = figNum+1;
    
    plot(T_Zenith_SNR,P_Zenith_SNR,'b');
    xlim([0 plot_Freq_Lim]);
    hold on;
    % plot(T_SNR_08apr08_10_or,P_SNR_08apr08_10_or,'r');
    % hold on;
    plot(T_Zenith_SNR,gf_Zenith_SNR,'LineWidth',1.5,'color',[0.5 0.5 0.5]);
    xlim([0 plot_Freq_Lim]);  
    xlabel('Periods (min)');
    ylabel('Normalized Power');
    title(['Delta sTEC at deg from ',num2str(latPoints(gg)),' deg lat/',num2str(longPoints(gg)),'deg long on Day ',num2str(floor(mean(time))),', 2013 - SCARGLED']);
    set(gca,'fontweight',font_Weight, 'FontSize',font_Size);
    set(gca,'xtick',Xaxisvar_SCARGLE);
end
toc

end %*********END OF FLG_gatherDataAtPoint_Filter_Scargle*********


%% Gather Data At Point - BAND PASS FILTER STEC DATA
if( FLG_gatherDataAtPoint_Filter_BP == 1 )

% sTEC_time_unique = unique(sTEC_holder{jk,1}(:,4)); %get unique times for sTEC
% sTEC_combined = zeros(size(sTEC_time_unique,1),1); %preallocate
% for( i = 1:size(sTEC_time_unique,1) ) %combine data at same time into one value at one time
%     
%     jp = find( sTEC_time_unique(i) == sTEC_holder{jk,1}(:,4) ); %find the data points at the same time
% 
%     sTEC_combined(i) = mean(  sTEC_holder{jk,1}(jp,1) ); %TECU, average data at same time into one value
% end

% Zenith_SNR_bp_expanded = interp1(Zenith_time,Zenith_SNR_bp(:,Zenith_threeHunIndex),sTEC_combined,'spline'); %dB, expand Zenith_SNR_bp to be the same size as the sTEC_combined var
%this doesn't work btw look in next section for proper setup

%These are unused, but accentuate the periods desired
%Original file used highpass_fir.m, it has been merged
% lp=1; % Lower period (in hrs)
% hp=2; % Higher period (in hrs)
% lf=(1/hp);  % Lowpass frequency corner (1/hr)
% hf=(1/lp);  % Highpass frequency corner (1/hr)

bs=(1/2);   % highpass cuttoff frequency (not sure what to make of it)
%I think it is related to lf above (which is 1/2 )

n=42; % order of the Hamming window used in the custom function
% c = 3.32*pi; %Hamming constant
% M = n/2;
% bandwidth = c/M;
%The above was included to investigate the bandwidth - 1/2. Related to
%bs/lf?

fp=bs; % stop band stoppin freq (or pass band passing freq, depending on how you look at it)
timeDelta = timeUnique(50)-timeUnique(49); %hrs already
f= 1/(timeDelta); % the sampling frequency, based off of the time delta calc'd

wp=2*fp/f; % Normalizing the frequencies (Matlab is 0 to 1)

%Calculation of filter coefficients
% [b,a]=fir1(n,wp,'high'); %This takes the order and lower cut off
%Uses the default hamming window
% frequency ORIG
W = hann(n+1);
[b,a] = fir1(n,wp,'high',W);
%Applys Hanning Window for shiggles and giggles

sTEC_combined_bp = zeros(size(sTEC_combined)); %preallocate
for(jk = 1:size(latPoints,2) ) %gather data for multiple points
    sTEC_combined_bp(jk,:) = filtfilt(b,a,sTEC_combined(jk,:)); %Appies the filter
end

end %*********END OF FLG_gatherDataAtPoint_Filter_BP*********


%% Plot Zenith ISR and then sTEC filtered to Compare
if( FLG_ISRdata_Zenith_Plot == 1)

% figure('units','normalized','outerposition',[0 0 1 1]); %maximizes figure window
figure(figNum);
figNum = figNum+1;
h = pcolor((Zenith_time-floor(mean(Zenith_time))).*24,Zenith_height,Zenith_SNR_bp');
set(h, 'EdgeColor', 'none'); %remove grid lines (which make it look hilarious)
shading interp;
lighting phong;
colorbar;
caxis([-0.1 0.1]);
colormap('gray');
ylim([90 700]);
% ylim([90 500]);
% xlim([0 70]);
xlabel('Time in UT - 0 Hr at 0 UT May 7th');
ylabel('Height (km)');
title('May 6-7-8 2013, Zenith SNR 2 Hr Highpass Filtered');
set(gca,'fontweight',font_Weight, 'FontSize',font_Size);
% Xaxisvar = [-12, -8, -4, 0, 4, 8, 12, 16, 20, 24 28, 32, 36, 40, 44, 48]; %shows off important magnitude points
Xaxisvar = [-12,-8,-4,0, 4, 8, 12, 16, 20, 24 28, 32, 36, 40, 44, 48, 52, 56, 60, 64, 68, 72]; %shows off important magnitude points
set(gca,'xtick',Xaxisvar);

% figure('units','normalized','outerposition',[0 0 1 1]); %maximizes figure window
figure(figNum);
figNum = figNum+1;
h = plot( (timeUnique_combined_noNaN-floor(mean(time))).*24 ,sTEC_combined(gg_ZenithSite,:) ); %plot all at once on this
hold on;
set(h(1),'Color','g')
h = plot( (Zenith_time-floor(mean(Zenith_time))).*24 ,Zenith_SNR(:,Zenith_threeHunIndex) ); %plot all at once on this
set(h(1),'Color','b')
xlabel('Time (hr)'); %old: ,'fontweight','bold','FontSize',12 for all
ylabel('Delta Slant TEC (TECU) or Zenith SNR bp (dB)'); %sped up significantly when merged above
title(['Delta sTEC AND Zenith SNR at ',num2str(latMillstone),' deg lat/',num2str(longMillstone),'deg long on Day ',num2str(floor(mean(time))),', 2013']);
axis([-inf,inf,-inf,inf])
legend('sTEC high density data','Zenith SNR');

% figure('units','normalized','outerposition',[0 0 1 1]); %maximizes figure window
figure(figNum);
figNum = figNum+1;
h = plot( (timeUnique_combined_noNaN-floor(mean(time))).*24 ,sTEC_combined(gg_ZenithSite,:) ); %plot all at once on this
hold on;
set(h(1),'Color','g')
h = plot( (Zenith_time-floor(mean(Zenith_time))).*24 ,Zenith_SNR_bp(:,Zenith_threeHunIndex) ); %plot all at once on this
set(h(1),'Color','b')
xlabel('Time (hr)'); %old: ,'fontweight','bold','FontSize',12 for all
ylabel('Delta Slant TEC (TECU) or Zenith SNR bp (dB)'); %sped up significantly when merged above
title(['Delta sTEC AND Zenith SNR BP at ',num2str(latMillstone),' deg lat/',num2str(longMillstone),'deg long on Day ',num2str(floor(mean(time))),', 2013']);
axis([-inf,inf,-inf,inf])
legend('sTEC high density data','Zenith SNR band-passed');

end %*********END OF FLG_ISRdata_Zenith_Plot*********


%% Plot MISA ISR and then sTEC filtered to Compare
if( FLG_ISRdata_MISA_Plot == 1 )

% figure('units','normalized','outerposition',[0 0 1 1]); %maximizes figure window
figure(figNum);
figNum = figNum+1;
h = pcolor((MISA_time-floor(mean(MISA_time))).*24,MISA_height,MISA_SNR_bp');
set(h, 'EdgeColor', 'none'); %remove grid lines (which make it look hilarious)
shading interp;
lighting phong;
colorbar;
caxis([-0.1 0.1]);
colormap('gray');
ylim([90 700]);
% ylim([90 500]);
% xlim([0 70]);
xlabel('Time in UT - 0 Hr at 0 UT May 7th');
ylabel('Height (km)');
title('May 6-7-8 2013, MISA SNR 2 Hr Highpass Filtered');
set(gca,'fontweight',font_Weight, 'FontSize',font_Size);
% Xaxisvar = [-12, -8, -4, 0, 4, 8, 12, 16, 20, 24 28, 32, 36, 40, 44, 48]; %shows off important magnitude points
Xaxisvar = [-12,-8,-4,0, 4, 8, 12, 16, 20, 24 28, 32, 36, 40, 44, 48, 52, 56, 60, 64, 68, 72]; %shows off important magnitude points
set(gca,'xtick',Xaxisvar);

% figure('units','normalized','outerposition',[0 0 1 1]); %maximizes figure window
figure(figNum);
figNum = figNum+1;
h = plot( (timeUnique-floor(mean(time))).*24 ,sTEC_combined(gg_ZenithSite,:) ); %plot all at once on this
hold on;
set(h(1),'Color','g')
h = plot( (MISA_time-floor(mean(MISA_time))).*24 ,MISA_SNR(:,MISA_threeHunIndex) ); %plot all at once on this
set(h(1),'Color','b')
xlabel('Time (hr)'); %old: ,'fontweight','bold','FontSize',12 for all
ylabel('Delta Slant TEC (TECU) or MISA SNR bp (dB)'); %sped up significantly when merged above
title(['Delta sTEC AND MISA SNR at ',num2str(latMillstone),' deg lat/',num2str(longMillstone),'deg long on Day ',num2str(floor(mean(time))),', 2013']);
axis([-inf,inf,-inf,inf])
legend('sTEC high density data','MISA SNR');

% figure('units','normalized','outerposition',[0 0 1 1]); %maximizes figure window
figure(figNum);
figNum = figNum+1;
h = plot( (timeUnique-floor(mean(time))).*24 ,sTEC_combined(gg_ZenithSite,:) ); %plot all at once on this
hold on;
set(h(1),'Color','g')
h = plot( (MISA_time-floor(mean(MISA_time))).*24 ,MISA_SNR_bp(:,MISA_threeHunIndex) ); %plot all at once on this
set(h(1),'Color','b')
xlabel('Time (hr)'); %old: ,'fontweight','bold','FontSize',12 for all
ylabel('Delta Slant TEC (TECU) or MISA SNR bp (dB)'); %sped up significantly when merged above
title(['Delta sTEC AND MISA SNR BP at ',num2str(latMillstone),' deg lat/',num2str(longMillstone),'deg long on Day ',num2str(floor(mean(time))),', 2013']);
axis([-inf,inf,-inf,inf])
legend('sTEC high density data','MISA SNR band-passed');

end %*********END OF FLG_ISRdata_MISA_Plot*********


%% Plot Zenith ISR and then sTEC filtered BP to Compare
if( FLG_ISRdata_Zenith_Plot_BP == 1)
    
% figure('units','normalized','outerposition',[0 0 1 1]); %maximizes figure window
figure(figNum);
figNum = figNum+1;
h = plot( (timeUnique_combined_noNaN-dateRange_zeroHr(2)).*24 ,sTEC_combined_bp(gg_ZenithSite,:) ); %plot all at once on this
hold on;
set(h(1),'Color','g')
h = plot( (Zenith_time - dateRange_zeroHr(2)).*24 ,Zenith_SNR_bp(:,Zenith_threeHunIndex) ); %plot all at once on this
set(h(1),'Color','b')
xlabel('Time (hr)'); %old: ,'fontweight','bold','FontSize',12 for all
ylabel('Delta Slant TEC bp (TECU) or Zenith SNR bp (dB)'); %sped up significantly when merged above
title(['Delta sTEC BP AND Zenith SNR BP at ',num2str(latMillstone),' deg lat/',num2str(longMillstone),'deg long on Day ',num2str(dateRange_zeroHr(2)),', ',num2str(dateRange_zeroHr(1))]);
axis([-inf,inf,-inf,inf])
legend('sTEC high density data band-passed','Zenith SNR band-passed');

% figure('units','normalized','outerposition',[0 0 1 1]); %maximizes figure window
figure(figNum);
figNum = figNum+1;
h = plot( (timeUnique_combined_noNaN - dateRange_zeroHr(2)).*24 ,sTEC_combined_bp(gg_ZenithSite,:) ); %plot all at once on this
hold on;
set(h(1),'Color','g')
h = plot( (Zenith_time - dateRange_zeroHr(2)).*24 ,Zenith_SNR_threeHun_AVGD ); %plot all at once on this
set(h(1),'Color','b')
xlabel('Time (hr)'); %old: ,'fontweight','bold','FontSize',12 for all
ylabel('Delta Slant TEC bp (TECU) or Zenith SNR bp (dB)'); %sped up significantly when merged above
title(['Delta sTEC BP AND Zenith SNR BP AVG''d +/-',num2str(Zenith_avgRange),' km around ',num2str(Zenith_threeHun),' km at ',num2str(latMillstone),' deg lat/',num2str(longMillstone),'deg long on Day ',num2str(dateRange_zeroHr(2)),', ',num2str(dateRange_zeroHr(1))]);
axis([-inf,inf,-inf,inf])
legend('sTEC high density data band-passed','Zenith SNR band-passed AVG''d');

end %*********END OF FLG_ISRdata_Zenith_Plot_BP*********


%% Plot MISA ISR and then sTEC filtered BP to Compare
if( FLG_ISRdata_MISA_Plot_BP == 1)
    
% figure('units','normalized','outerposition',[0 0 1 1]); %maximizes figure window
figure(figNum);
figNum = figNum+1;
h = plot( (timeUnique-floor(mean(time))).*24 ,sTEC_combined_bp(gg_ZenithSite,:) ); %plot all at once on this
hold on;
set(h(1),'Color','g')
h = plot( (MISA_time-floor(mean(MISA_time))).*24 ,MISA_SNR_bp(:,MISA_threeHunIndex) ); %plot all at once on this
set(h(1),'Color','b')
xlabel('Time (hr)'); %old: ,'fontweight','bold','FontSize',12 for all
ylabel('Delta Slant TEC bp (TECU) or MISA SNR bp (dB)'); %sped up significantly when merged above
title(['Delta sTEC BP AND MISA SNR BP at ',num2str(latMillstone),' deg lat/',num2str(longMillstone),'deg long on Day ',num2str(floor(mean(time))),', 2013']);
axis([-inf,inf,-inf,inf])
legend('sTEC high density data band-passed','MISA SNR band-passed');

% figure('units','normalized','outerposition',[0 0 1 1]); %maximizes figure window
figure(figNum);
figNum = figNum+1;
h = plot( (timeUnique-floor(mean(time))).*24 ,sTEC_combined_bp(gg_ZenithSite,:) ); %plot all at once on this
hold on;
set(h(1),'Color','g')
h = plot( (MISA_time-floor(mean(MISA_time))).*24 ,MISA_SNR_threeHun_AVGD ); %plot all at once on this
set(h(1),'Color','b')
xlabel('Time (hr)'); %old: ,'fontweight','bold','FontSize',12 for all
ylabel('Delta Slant TEC bp (TECU) or MISA SNR bp (dB)'); %sped up significantly when merged above
title(['Delta sTEC BP AND MISA SNR BP AVG''d +/-',num2str(MISA_avgRange),' km around ',num2str(MISA_threeHun),' km at ',num2str(latMillstone),' deg lat/',num2str(longMillstone),'deg long on Day ',num2str(floor(mean(time))),', 2013']);
axis([-inf,inf,-inf,inf])
legend('sTEC high density data band-passed','MISA SNR band-passed AVG''d');

end %*********END OF FLG_ISRdata_MISA_Plot_BP*********


%% Do some Correlation Comparison Dealio - ZENITH
if( FLG_ISRdata_Zenith_Corr == 1)

tj = find( min( abs(min(Zenith_time) - timeUnique) ) == abs(min(Zenith_time) - timeUnique) ); %find min timeUnique that fits ISR time
tk = find( min( abs(max(Zenith_time) - timeUnique) ) == abs(max(Zenith_time) - timeUnique) ); %find max timeUnique that fits ISR time
sTEC_combined_limISR_Zenith = sTEC_combinedInterped(gg_ZenithSite,tj:tk); %cut out time that matches ISR time
    
R_Zenith = corrcoef(sTEC_combined_limISR_Zenith,Zenith_SNR_threeHun_AVGD_expanded);
fprintf('Correlation coeff between Zenith ISR and sTEC over Millstonefull time: %f\n',R_Zenith(2));

% figure('units','normalized','outerposition',[0 0 1 1]); %maximizes figure window
figure(figNum);
figNum = figNum+1;
h = plot( (timeUnique_limISR_Zenith-floor(mean(timeUnique_limISR_Zenith))).*24 ,sTEC_combined_limISR_Zenith ); %plot all at once on this
hold on;
set(h(1),'Color','g')
h = plot( (timeUnique_limISR_Zenith-floor(mean(timeUnique_limISR_Zenith))).*24 ,Zenith_SNR_threeHun_AVGD_expanded ); %plot all at once on this
set(h(1),'Color','b')
xlabel('Time (hr)'); %old: ,'fontweight','bold','FontSize',12 for all
ylabel('Delta Slant TEC bp (TECU) or Zenith SNR bp (dB)'); %sped up significantly when merged above
title(['Delta sTEC BP AND Zenith SNR BP AVG''d w/ Interp. +/-',num2str(Zenith_avgRange),' km around ',num2str(Zenith_threeHun),' km at ',num2str(latMillstone),' deg lat/',num2str(longMillstone),'deg long on Day ',num2str(floor(mean(time))),', 2013']);
axis([-inf,inf,-inf,inf])
legend('sTEC high density data band-passed','Zenith SNR band-passed AVG''d');

end %*********END OF FLG_ISRdata_Zenith_Corr*********


%% Do some Correlation Comparison Dealio - MISA
if( FLG_ISRdata_MISA_Corr == 1)

tj = find( min( abs(min(MISA_time) - timeUnique) ) == abs(min(MISA_time) - timeUnique) ); %find min timeUnique that fits ISR time
tk = find( min( abs(max(MISA_time) - timeUnique) ) == abs(max(MISA_time) - timeUnique) ); %find max timeUnique that fits ISR time
sTEC_combined_limISR_MISA = sTEC_combinedInterped(gg_MISASite,tj:tk); %cut out time that matches ISR time

R_MISA = corrcoef(sTEC_combined_limISR_MISA,MISA_SNR_threeHun_AVGD_expanded);
fprintf('Correlation coeff between MISA ISR and sTEC over Millstone: %f\n\n',R_MISA(2));

% figure('units','normalized','outerposition',[0 0 1 1]); %maximizes figure window
figure(figNum);
figNum = figNum+1;
h = plot( (timeUnique_limISR_MISA-floor(mean(timeUnique_limISR_MISA))).*24 ,sTEC_combined_limISR_MISA ); %plot all at once on this
hold on;
set(h(1),'Color','g')
h = plot( (timeUnique_limISR_MISA-floor(mean(timeUnique_limISR_MISA))).*24 ,MISA_SNR_threeHun_AVGD_expanded ); %plot all at once on this
set(h(1),'Color','b')
xlabel('Time (hr)'); %old: ,'fontweight','bold','FontSize',12 for all
ylabel('Delta Slant TEC bp (TECU) or MISA SNR bp (dB)'); %sped up significantly when merged above
title(['Delta sTEC BP AND MISA SNR BP AVG''d w/ Interp. +/-',num2str(MISA_avgRange),' km around ',num2str(MISA_threeHun),' km at ',num2str(latMillstone),' deg lat/',num2str(longMillstone),'deg long on Day ',num2str(floor(mean(time))),', 2013']);
axis([-inf,inf,-inf,inf])
legend('sTEC high density data band-passed','MISA SNR band-passed');

end %*********END OF FLG_ISRdata_MISA_Corr*********


%% Average Longitude, Plot Time vs Latitude
if( FLG_AVG_LONG == 1 )

pplatChunks = linspace(min(pplat),max(pplat),pplatN+1); %bits of latitude to average within
pplatChunksPlot = linspace( (pplatChunks(2) + pplatChunks(1))/2 , (pplatChunks(length(pplatChunks)) + pplatChunks(length(pplatChunks)-1))/2 , pplatN); %for plotting, the mean between each point

pplongRange = [longMillstone-pplongAvgDeg,longMillstone+pplongAvgDeg]; %deg, range to average between

% pplongChunked = zeros(pplatN,length(timeUnique)); %preallocate
sTECChunked_longavg = zeros(length(timeUnique_limISR_MISA),pplatN); %preallocate
parfor( i = 1:length(timeUnique_limISR_MISA) )
    %Corral the data to the right place    
    k = time == timeUnique_limISR_MISA(i); %gets during a time period
    
    temp_sTEC = sTEC(k); %record the sTEC for the time
    temp_pplat = pplat(k); %record the pierce-point latitude
    temp_pplong = pplong(k); %record the pierce-point longitude
    kk = temp_pplong >= min(pplongRange) & temp_pplong <= max(pplongRange); %get the range of long to keep
    temp_sTEC = temp_sTEC(kk); %these get cut too to match the counting
    temp_pplat = temp_pplat(kk); %these get cut too to match the counting
%     temp_pplong = temp_pplong(kk); %cut out only the long wanted

    for(j = 1: pplatN )
        %average sTEC for a latitude range chunk within the given longitude
        %range
        kl = (pplatChunks(j) <= temp_pplat) & (pplatChunks(j+1) >= temp_pplat) ; %get pplats in range
        
        sTECChunked_longavg(i,j) = mean(temp_sTEC(kl)); %average the sTEC in this lat band with the given long range
    end        
end
clear temp_sTEC temp_pplat temp_pplong

%prep figure #
% figure('units','normalized','outerposition',[0 0 1 1]); %maximizes figure window
figure(figNum);
figNum = figNum+1;
%Prep plot space
ax1 = subplot(2,1,1); %upper bit
h1 = pcolor( (MISA_time-floor(mean(MISA_time))).*24 ,MISA_height,MISA_SNR_bp');
set(h1, 'EdgeColor', 'none'); %remove grid lines (which make it look hilarious)
shading interp;
lighting phong;
colorbar;
caxis([-0.1 0.1]);
colormap(ax1,'gray');
ylim([90 700]);
% ylim([90 500]);
% xlim([0 70]);
xlabel('Time in UT - 0 Hr at 0 UT May 7th');
ylabel('Height (km)');
title('May 6-7-8 2013, MISA SNR 2 Hr Highpass Filtered');
set(gca,'fontweight',font_Weight, 'FontSize',font_Size);
% Xaxisvar = [-12, -8, -4, 0, 4, 8, 12, 16, 20, 24 28, 32, 36, 40, 44, 48]; %shows off important magnitude points
%         Xaxisvar = [-12,-8,-4,0, 4, 8, 12, 16, 20, 24 28, 32, 36, 40, 44, 48, 52, 56, 60, 64, 68, 72]; %shows off important magnitude points
Xaxisvar = round((min(MISA_time)-floor(mean(MISA_time)))*24):1:round((max(MISA_time)--floor(mean(MISA_time)))*24);
set(gca,'xtick',Xaxisvar);
% 
ax2 = subplot(2,1,2); %lower bit
h = pcolor( (timeUnique_limISR_MISA -floor(mean(timeUnique))).*24 , pplatChunksPlot , sTECChunked_longavg');% pseudocolor plot "stretched" to the grid
set(h, 'EdgeColor', 'none'); %remove grid lines (which make it look hilarious)
%             shading interp; %stuff for interpolation
%             lighting phong; %stuff for interpolation
caxis([-gif_sTEC_cap, gif_sTEC_cap]);
%         colormap(flipud(colormap('hot'))) %heat map where 0 (low) is white
colormap(ax2,'jet') %heat map
% colormap('jet');
colorbar; %shows the color bar

hold on;

%Now drawing line of interest
line([ min((timeUnique_limISR_MISA -floor(mean(timeUnique))).*24) , max((timeUnique_limISR_MISA -floor(mean(timeUnique))).*24) ],[latMillstone,latMillstone],'Color','k','LineWidth',.25);  %plots a point with a blue *

xlabel('Time, 0 Hr on Day 127 2016 (hr)'); %old: ,'fontweight','bold','FontSize',12 for all
ylabel('Latitude (deg)'); %sped up significantly when merged above
title(['TEC Averaged from Longitudes ',num2str(max(pplongRange)),' deg to ',num2str(min(pplongRange)),' deg -- 0 Hr on May 7th, 2013']);
hold off;

end %*********END OF FLG_AVG_LONG*********


%% Average Latitude, Plot Time vs Longitude
if( FLG_AVG_LAT == 1 )

pplongChunks = linspace(min(pplong),max(pplong),pplongN+1); %bits of latitude to average within
pplongChunksPlot = linspace( (pplongChunks(2) + pplongChunks(1))/2 , (pplongChunks(length(pplongChunks)) + pplongChunks(length(pplongChunks)-1))/2 , pplongN); %for plotting, the mean between each point

pplatRange = [latMillstone-pplatAvgDeg,latMillstone+pplatAvgDeg]; %deg, range to average between

% sTEC_RANDTEST = 2.*(rand(size(sTEC))-.5); %1 to -1 random noise

% pplongChunked = zeros(pplatN,length(timeUnique)); %preallocate
sTECChunked_latavg = zeros(length(timeUnique_limISR_MISA),pplongN); %preallocate
parfor( i = 1:length(timeUnique_limISR_MISA) )
    %Corral the data to the right place    
    k = time == timeUnique_limISR_MISA(i); %gets during a time period
    
    temp_sTEC = sTEC(k); %record the sTEC for the time
    temp_pplat = pplat(k); %record the pierce-point latitude
    temp_pplong = pplong(k); %record the pierce-point longitude
    kk = temp_pplat >= min(pplatRange) & temp_pplat <= max(pplatRange); %get the range of long to keep
    temp_sTEC = temp_sTEC(kk); %these get cut too to match the counting
    temp_pplong = temp_pplong(kk); %these get cut too to match the counting
%     temp_pplat = temp_pplat(kk); %cut out only the lat wanted

    for(j = 1: pplongN )
        %average sTEC for a latitude range chunk within the given longitude
        %range
        kl = (pplongChunks(j) <= temp_pplong) & (pplongChunks(j+1) >= temp_pplong) ; %get pplats in range
        
        sTECChunked_latavg(i,j) = mean(temp_sTEC(kl)); %average the sTEC in this lat band with the given long range
    end        
end
clear temp_sTEC temp_pplat temp_pplong

%prep figure #
% figure('units','normalized','outerposition',[0 0 1 1]); %maximizes figure window
figure(figNum);
figNum = figNum+1;
%Prep plot space
ax1 = subplot(2,1,1); %upper bit
h1 = pcolor( (MISA_time-floor(mean(MISA_time))).*24 ,MISA_height,MISA_SNR_bp');
set(h1, 'EdgeColor', 'none'); %remove grid lines (which make it look hilarious)
shading interp;
lighting phong;
colorbar;
caxis([-0.1 0.1]);
colormap(ax1,'gray');
ylim([90 700]);
% ylim([90 500]);
% xlim([0 70]);
xlabel('Time in UT - 0 Hr at 0 UT May 7th');
ylabel('Height (km)');
title('May 6-7-8 2013, MISA SNR 2 Hr Highpass Filtered');
set(gca,'fontweight',font_Weight, 'FontSize',font_Size);
% Xaxisvar = [-12, -8, -4, 0, 4, 8, 12, 16, 20, 24 28, 32, 36, 40, 44, 48]; %shows off important magnitude points
%         Xaxisvar = [-12,-8,-4,0, 4, 8, 12, 16, 20, 24 28, 32, 36, 40, 44, 48, 52, 56, 60, 64, 68, 72]; %shows off important magnitude points
Xaxisvar = round((min(MISA_time)-floor(mean(MISA_time)))*24):1:round((max(MISA_time)--floor(mean(MISA_time)))*24);
set(gca,'xtick',Xaxisvar);
% 
% %add on a vertical line that shows where we are at in time
% 
ax2 = subplot(2,1,2); %lower bit
h = pcolor( (timeUnique_limISR_MISA -floor(mean(timeUnique))).*24 , pplongChunksPlot , sTECChunked_latavg');% pseudocolor plot "stretched" to the grid
set(h, 'EdgeColor', 'none'); %remove grid lines (which make it look hilarious)
%             shading interp; %stuff for interpolation
%             lighting phong; %stuff for interpolation
caxis([-gif_sTEC_cap, gif_sTEC_cap]);
%         colormap(flipud(colormap('hot'))) %heat map where 0 (low) is white
colormap(ax2,'jet');
colorbar; %shows the color bar

hold on;
%Now drawing line of interest
line([ min((timeUnique_limISR_MISA -floor(mean(timeUnique))).*24) , max((timeUnique_limISR_MISA -floor(mean(timeUnique))).*24) ],[longMillstone,longMillstone],'Color','k','LineWidth',.25);  %plots a point with a blue *

xlabel('Time, 0 Hr on Day 127 2016 (hr)'); %old: ,'fontweight','bold','FontSize',12 for all
ylabel('Longitude (deg)'); %sped up significantly when merged above
title(['TEC Averaged from Latitudes ',num2str(max(pplatRange)),' deg to ',num2str(min(pplatRange)),' deg -- 0 Hr on May 7th, 2013']);
hold off;


%SECOND IS TIME ON Y AXIS
%prep figure #
% figure('units','normalized','outerposition',[0 0 1 1]); %maximizes figure window
figure(figNum);
figNum = figNum+1;
%Prep plot space
ax1 = subplot(1,2,1); %upper bit
h1 = pcolor( MISA_height, (MISA_time-floor(mean(MISA_time))).*24 , MISA_SNR_bp);
set(h1, 'EdgeColor', 'none'); %remove grid lines (which make it look hilarious)
shading interp;
lighting phong;
colorbar;
caxis([-0.1 0.1]);
colormap(ax1,'gray');
xlim([90 700]);
% ylim([90 500]);
% xlim([0 70]);
ylabel('Time in UT - 0 Hr at 0 UT May 7th');
xlabel('Height (km)');
title('May 6-7-8 2013, MISA SNR 2 Hr Highpass Filtered');
set(gca,'fontweight',font_Weight, 'FontSize',font_Size);
% Xaxisvar = [-12, -8, -4, 0, 4, 8, 12, 16, 20, 24 28, 32, 36, 40, 44, 48]; %shows off important magnitude points
%         Xaxisvar = [-12,-8,-4,0, 4, 8, 12, 16, 20, 24 28, 32, 36, 40, 44, 48, 52, 56, 60, 64, 68, 72]; %shows off important magnitude points
Yaxisvar = round((min(MISA_time)-floor(mean(MISA_time)))*24):1:round((max(MISA_time)--floor(mean(MISA_time)))*24);
set(gca,'ytick',Yaxisvar);
% 
ax2 = subplot(1,2,2); %lower bit
h = pcolor( pplongChunksPlot ,(timeUnique_limISR_MISA -floor(mean(timeUnique))).*24 , sTECChunked_latavg );% pseudocolor plot "stretched" to the grid
set(h, 'EdgeColor', 'none'); %remove grid lines (which make it look hilarious)
%             shading interp; %stuff for interpolation
%             lighting phong; %stuff for interpolation
caxis([-gif_sTEC_cap, gif_sTEC_cap]);
%         colormap(flipud(colormap('hot'))) %heat map where 0 (low) is white
colormap(ax2,'jet') %heat map
colorbar; %shows the color bar

hold on;

%Now drawing line of interest
% line([ min((timeUnique -floor(mean(timeUnique))).*24) , max((timeUnique -floor(mean(timeUnique))).*24) ],[longMillstone,longMillstone],'Color','k','LineWidth',.25);  %plots a point with a blue *
line([longMillstone,longMillstone],[ min((timeUnique_limISR_MISA -floor(mean(timeUnique))).*24) , max((timeUnique_limISR_MISA -floor(mean(timeUnique))).*24) ],'Color','k','LineWidth',.25);  %plots a point with a blue *

ylabel('Time, 0 Hr on Day 127 2016 (hr)'); %old: ,'fontweight','bold','FontSize',12 for all
xlabel('Longitude (deg)'); %sped up significantly when merged above
title(['TEC Averaged from Latitudes ',num2str(max(pplatRange)),' deg to ',num2str(min(pplatRange)),' deg -- 0 Hr on May 7th, 2013']);
hold off;

end %*********END OF FLG_AVG_LAT*********


%% Average Any Angle, user defined width
if( FLG_AVG_ANYANGLE == 1 )
    
disp('TIME TO RUN ANY ANGLE ALG'); 
tic
avg_anyAngle = avg_anyAngle*pi/180; %rad, convert to radians

%
avg_anyAngle_slope = tan(avg_anyAngle); %get slope of line required
%this conversion is for y=LATITUDE x=LONGITUDE line action

%idea here is if 5 arcdeg width, half it and go 2.5 arcdeg(x-y at this point - c denotes) 90 deg (real angle) to the req
%angle - y=Latitude x=Longitude solving for intercept with 2 points and
%slope
avg_anyAngle_upLine_int = avg_anyAngle_Width/2*sin(avg_anyAngle + pi/2) + mean(plotLatRange)  ...
    - avg_anyAngle_slope*(avg_anyAngle_Width/2*cos(avg_anyAngle + pi/2) + mean(plotLongRange)); %get intercept of upper line
%upper and lower lines are parallel
%idea here is if 5 arcdeg width, half it and go 2.5 arcdeg(x-y at this point - c denotes) -90 deg (real angle) to the req
%angle - y=Latitude x=Longitude solving for intercept with 2 points and
%slope
avg_anyAngle_loLine_int = avg_anyAngle_Width/2*sin(avg_anyAngle - pi/2) + mean(plotLatRange)  ...
    - avg_anyAngle_slope*(avg_anyAngle_Width/2*cos(avg_anyAngle - pi/2) + mean(plotLongRange)); %get intercept of lower line

% avg_anyAngle_LimMaxLong = max(plotLongRange); %arcdeg, get max longitude
% avg_anyAngle_LimMinLong = min(plotLongRange); %arcdeg, get min longitude

avg_anyAngle_LatLim_upLine = avg_anyAngle_slope*sort(plotLongRange) + avg_anyAngle_upLine_int; %arcdeg,
%latitude range from largest possible longitude range
avg_anyAngle_LatLim_upLine( avg_anyAngle_LatLim_upLine > max(plotLatRange)) = max(plotLatRange); %arcdeg, limit lat to current range
avg_anyAngle_LatLim_upLine( avg_anyAngle_LatLim_upLine < min(plotLatRange)) = min(plotLatRange); %arcdeg, limit lat to current range

% avg_anyAngle_LongLim_upLine( avg_anyAngle_LongLim_upLine > max(plotLongRange)) = max(plotLongRange); %arcdeg, limit long to current range
% avg_anyAngle_LongLim_upLine( avg_anyAngle_LongLim_upLine < min(plotLongRange)) = min(plotLongRange); %arcdeg, limit long to current range
%theo not needed

% avg_anyAngle_LatLim_loLine = avg_anyAngle_slope*plotLongRange + avg_anyAngle_loLine_int; %arcdeg,
% %latitude range from largest possible longitude range
% avg_anyAngle_LatLim_loLine( avg_anyAngle_LatLim_loLine > max(plotLatRange)) = max(plotLatRange); %arcdeg, limit lat to current range
% avg_anyAngle_LatLim_loLine( avg_anyAngle_LatLim_loLine < min(plotLatRange)) = min(plotLatRange); %arcdeg, limit lat to current range
% avg_anyAngle_LongLim_loLine = (avg_anyAngle_LatLim_loLine - avg_anyAngle_loLine_int)/avg_anyAngle_slope; %arcdeg, get longitudes that match

%Project upper line pts to lower pts
avg_anyAngle_LatLim_loLine = [min(avg_anyAngle_LatLim_upLine) + avg_anyAngle_Width*sin(avg_anyAngle - pi/2) , max(avg_anyAngle_LatLim_upLine) + avg_anyAngle_Width*sin(avg_anyAngle - pi/2)]; %arcdeg, calc lower pts based on upper
avg_anyAngle_LatLim_loLine( avg_anyAngle_LatLim_loLine > max(plotLatRange)) = max(plotLatRange); %arcdeg, limit lat to current range
avg_anyAngle_LatLim_loLine( avg_anyAngle_LatLim_loLine < min(plotLatRange)) = min(plotLatRange); %arcdeg, limit lat to current range
%redefine upper pts based on possibly adjusted lower pts
avg_anyAngle_LatLim_upLine = [min(avg_anyAngle_LatLim_loLine) + avg_anyAngle_Width*sin(avg_anyAngle + pi/2) , max(avg_anyAngle_LatLim_loLine) + avg_anyAngle_Width*sin(avg_anyAngle + pi/2)]; %arcdeg, calc upper pts based on lower
avg_anyAngle_LongLim_upLine = (avg_anyAngle_LatLim_upLine - avg_anyAngle_upLine_int)/avg_anyAngle_slope; %arcdeg, get longitudes that match
avg_anyAngle_LongLim_loLine = (avg_anyAngle_LatLim_loLine - avg_anyAngle_loLine_int)/avg_anyAngle_slope; %arcdeg, get longitudes that match

avg_anyAngle_Range = [ [avg_anyAngle_LatLim_upLine' ; avg_anyAngle_LatLim_loLine'] , [avg_anyAngle_LongLim_upLine' ; avg_anyAngle_LongLim_loLine'] ]; %arcdeg, record pts for use
% 
% avg_anyAngle_Range(3,:) = [avg_anyAngle_Range(1,1) + avg_anyAngle_Width*sin(avg_anyAngle - pi/2), avg_anyAngle_Range(1,2) + avg_anyAngle_Width*cos(avg_anyAngle - pi/2)]; %arcdeg, record max up pts
% avg_anyAngle_Range(4,:) = [avg_anyAngle_Range(2,1) + avg_anyAngle_Width*sin(avg_anyAngle - pi/2), avg_anyAngle_Range(2,2) + avg_anyAngle_Width*cos(avg_anyAngle - pi/2)]; %arcdeg, record max up pts
% 
% %Chec
% if( (max(avg_anyAngle_LatLim_upLine) + avg_anyAngle_Width*sin(avg_anyAngle - pi/2)) < max(plotLatRange) && (max(avg_anyAngle_LongLim_upLine) + avg_anyAngle_Width*cos(avg_anyAngle - pi/2)) < max(plotLongRange) )
% 
% elseif( (max(avg_anyAngle_LatLim_loLine) + avg_anyAngle_Width*sin(avg_anyAngle + pi/2)) < max(plotLatRange) && (max(avg_anyAngle_LongLim_loLine) + avg_anyAngle_Width*cos(avg_anyAngle + pi/2)) < max(plotLongRange) )
% %check lo now
%     avg_anyAngle_LineLim = 2; %1 for up, 2 for lo
% end
% if( avg_anyAngle_LineLim == 0 )
%     error('Any Angle Avg did not find a suitable limit.')
% end
% %both prob won't happen so I won't go there
% 
% avg_anyAngle_Range = zeros(4,2); %preallocate, [LAT, LONG] setup OR [Y, X]
% if( avg_anyAngle_LineLim == 1 )
%     %up was limit, get pts from it
%     avg_anyAngle_Range(1,:) = [min(avg_anyAngle_LatLim_upLine), min(avg_anyAngle_LongLim_upLine)]; %arcdeg, record min up pts
%     avg_anyAngle_Range(2,:) = [max(avg_anyAngle_LatLim_upLine), max(avg_anyAngle_LongLim_upLine)]; %arcdeg, record max up pts
%     avg_anyAngle_Range(3,:) = [avg_anyAngle_Range(1,1) + avg_anyAngle_Width*sin(avg_anyAngle - pi/2), avg_anyAngle_Range(1,2) + avg_anyAngle_Width*cos(avg_anyAngle - pi/2)]; %arcdeg, record max up pts
%     avg_anyAngle_Range(4,:) = [avg_anyAngle_Range(2,1) + avg_anyAngle_Width*sin(avg_anyAngle - pi/2), avg_anyAngle_Range(2,2) + avg_anyAngle_Width*cos(avg_anyAngle - pi/2)]; %arcdeg, record max up pts
% else
%     %lo was limit, get pts from it 
%     avg_anyAngle_Range = avg_anyAngle_LatLim_loLine; %arcdeg, get box pts to get data from
%     avg_anyAngle_LongRange = avg_anyAngle_LongLim_loLine; %arcdeg, get box pts to get data from
% end


% avg_anyAngle_Range_Length = sqrt( (avg_anyAngle_Range(1,2) - avg_anyAngle_Range(2,2))^2 + (avg_anyAngle_Range(1,1) - avg_anyAngle_Range(2,1))^2 ); %arcdeg, length of line
%other one should be same lenght bc science 

avg_anyAngle_Range_Chunks_Long_up = linspace( min(avg_anyAngle_Range(1:2,2)) , max(avg_anyAngle_Range(1:2,2)) , avg_anyAngle_N + 1)'; %arcdeg, chunks in longitude to go between
%chose longitude because 0 deg will be stable - 90 deg would never be with
%my hella math *wasn't stable at 0 anyway lol*
avg_anyAngle_Range_Chunks_Long_lo = linspace( min(avg_anyAngle_Range(3:4,2)) , max(avg_anyAngle_Range(3:4,2)) , avg_anyAngle_N + 1)'; %arcdeg, chunks in longitude to go between
% avg_anyAngle_Range_Chunks_Lat = linspace( min(avg_anyAngle_Range(1:2,1)) , max(avg_anyAngle_Range(1:2,1)) , avg_anyAngle_N + 1)'; %arcdeg, chunks in latitude to go between

avg_anyAngle_Range_Chunks_Long_Plot = linspace( (avg_anyAngle_Range_Chunks_Long_up(1) + avg_anyAngle_Range_Chunks_Long_lo(1))/2 ,...
    (avg_anyAngle_Range_Chunks_Long_up(end) + avg_anyAngle_Range_Chunks_Long_lo(end))/2 , avg_anyAngle_N)'; %for plotting, the mean between each point
%for plotting using up

if( avg_anyAngle_45vsLatLong == 1 && (avg_anyAngle == 45*pi/180 || avg_anyAngle == 135*pi/180) ) %override mechanism to allow LATITUDE when normally longitude is defaulted to
    
    avg_anyAngle_Range_Chunks_Long_Plot_Name = 'Latitude'; %name of longitude
    
    avg_anyAngle_Range_Chunks_Long_Plot_int = mean(plotLatRange)  ...
        - avg_anyAngle_slope*mean(plotLongRange); %get intercept of upper line
    
    avg_anyAngle_Range_Chunks_Long_Plot = avg_anyAngle_slope.*avg_anyAngle_Range_Chunks_Long_Plot + avg_anyAngle_Range_Chunks_Long_Plot_int;
    %for plotting rocking it up
    
elseif( avg_anyAngle <= 45*pi/180 || avg_anyAngle >= 135*pi/180 ) %actually LONGITUDE on the axis

    avg_anyAngle_Range_Chunks_Long_Plot_Name = 'Longitude'; %name of longitude

else %otherwise LATITUDE on the axis
    avg_anyAngle_Range_Chunks_Long_Plot_Name = 'Latitude'; %name of longitude
    
    avg_anyAngle_Range_Chunks_Long_Plot_int = mean(plotLatRange)  ...
        - avg_anyAngle_slope*mean(plotLongRange); %get intercept of upper line
    
    avg_anyAngle_Range_Chunks_Long_Plot = avg_anyAngle_slope.*avg_anyAngle_Range_Chunks_Long_Plot + avg_anyAngle_Range_Chunks_Long_Plot_int;
    %for plotting rocking it up
end

% pplatRange = [latMillstone-pplatAvgDeg,latMillstone+pplatAvgDeg]; %deg, range to average between
% pplongChunked = zeros(pplatN,length(timeUnique)); %preallocate
temp_Longs_up = zeros(avg_anyAngle_N,2); %preallocate
temp_Lats_up = zeros(avg_anyAngle_N,2); %preallocate
temp_Longs_lo = zeros(avg_anyAngle_N,2); %preallocate
temp_Lats_lo = zeros(avg_anyAngle_N,2); %preallocate
for(j = 1:avg_anyAngle_N ) %preallocate and fill
    temp_Longs_up(j,:) = [avg_anyAngle_Range_Chunks_Long_up(j) , avg_anyAngle_Range_Chunks_Long_up(j+1)]; %arcdeg, get longitudes needed upper line
    temp_Longs_lo(j,:) = fliplr([avg_anyAngle_Range_Chunks_Long_lo(j) , avg_anyAngle_Range_Chunks_Long_lo(j+1)]); %arcdeg, get longitudes needed low
    temp_Lats_up(j,:) = avg_anyAngle_slope*temp_Longs_up(j,:) + avg_anyAngle_upLine_int; %arcdeg, get latitudes needed up
    temp_Lats_lo(j,:) = avg_anyAngle_slope*temp_Longs_lo(j,:) + avg_anyAngle_loLine_int; %arcdeg, get latitudes needed lower line
end

temp_Long_List = zeros(avg_anyAngle_N,5); %preallocate, 5 pts to define shape (start and end are same - really 4 for square thing)
temp_Lat_List = zeros(avg_anyAngle_N,5); %preallocate, 5 pts to define shape (start and end are same - really 4 for square thing)
for(j = 1:avg_anyAngle_N)
    temp_Long_List(j,:) = [temp_Longs_up(j,:),temp_Longs_lo(j,:),temp_Longs_up(j,1)]; %so this isn't done dynamtically
    temp_Lat_List(j,:) = [temp_Lats_up(j,:),temp_Lats_lo(j,:),temp_Lats_up(j,1)];
end


%THIS IS TO VIEW DATA AVG RANGE BEING TAKEN
coastLines = load('coast');%loads coast line into structure, vars are lat and long - can interfere with previous code so hidden in structure
coastLines_lat = coastLines.lat; %breaks this data out
coastLines_long = coastLines.long;
clear coastLines
if( FLG_geomagneticCoords == 1 ) %convert data to geomagnetic coords if desired
    [coastLines_lat,coastLines_long] = sFUN_geoToGeomag(dateOfData,magF,coastLines_lat,coastLines_long); %convert from geographic to geomagnetic coordinates
end %doing as we read data bit by bit because full math conversion will cause mem to run out with my 16 gigs omg
jk = find(coastLines_lat < min(plotLatRange) ); %find less than min, remove
coastLines_lat(jk) = []; %remove
coastLines_long(jk) = []; %remove
jk = find(coastLines_lat > max(plotLatRange) ); %find more than max, remove
coastLines_lat(jk) = []; %remove
coastLines_long(jk) = []; %remove
jk = find(coastLines_long < min(plotLongRange) ); %find less than min, remove
coastLines_lat(jk) = []; %remove
coastLines_long(jk) = []; %remove
jk = find(coastLines_long > max(plotLongRange) ); %find more than max, remove
coastLines_lat(jk) = []; %remove
coastLines_long(jk) = []; %remove
coastLines_long( coastLines_long == 180 ) = 180-.01; %fix some weird thing where 180 flips to -180

%Corral the data to the right place    
k = find( time == timeUnique_limISR_MISA(1)); %gets during a time period
%prep figure #
% figure('units','normalized','outerposition',[0 0 1 1]); %maximizes figure window
figure(figNum);
figNum = figNum+1;
scatter(pplong(k),pplat(k),plot_Scatter_Point_Size,sTEC(k),'.');
caxis([-gif_sTEC_cap, gif_sTEC_cap]);
%         colormap(flipud(colormap('hot'))) %heat map where 0 (low) is white
colormap('jet') %heat map
h = colorbar; %shows the color bar
ylabel(h, 'delta-sTEC (TECU)','fontweight',font_Weight, 'FontSize',font_Size); %labels color bar
hold on;
%Now drawing line of interest
plot(longMillstone,latMillstone,'*','Color',[255/255,36/255,0/255],'MarkerSize',25,'LineWidth',1.9); %plots a point with a red big *
hold on;
%Plot lines of contients
geoshow(coastLines_lat,coastLines_long,'Color','k') %add continental outlines in black  
xlabel('Longitude (arcdeg)'); %old: ,'fontweight','bold','FontSize',12 for all
ylabel('Latitude (arcdeg)'); %sped up significantly when merged above
string_Title = ['delta-sTEC Averaging with Angle = ',num2str(round(avg_anyAngle*180/pi,2)),' deg, Avg Total Width = ',num2str(round(avg_anyAngle_Width,2)),...
    ' arcdeg, Avg Step # = ',num2str(avg_anyAngle_N),...
    ' Avg per Step Width = ',num2str(round(sqrt( (temp_Longs_up(1,1) - temp_Longs_up(1,2))^2 + (temp_Lats_up(1,1) - temp_Lats_up(1,2))^2 ),2)),' arcdeg']; %create mecha title
if( FLG_geomagneticCoords == 1)
    string_Title = [string_Title,' [Mag Coord Aligned]']; %tack on mag coordinate warning
end
title(string_Title);

xticks( round(floor(min(plotLongRange)):plotLongRange_autoTick:ceil(max(plotLongRange)),2) ); %creates x ticks automagically
yticks( round(floor(min(plotLatRange)):plotLatRange_autoTick:ceil(max(plotLatRange)),2) ); %creates y ticks automagically
set(gca,'fontweight',font_Weight, 'FontSize',font_Size);
hold on;
%draw a box where the data will be gathered
line( [avg_anyAngle_Range(1:2,2) ; flipud(avg_anyAngle_Range(3:4,2)) ; avg_anyAngle_Range(1,2)],... %X longitude arcdeg
    [avg_anyAngle_Range(1:2,1) ; flipud(avg_anyAngle_Range(3:4,1)) ; avg_anyAngle_Range(1,1)],... %Y latitude arcdeg
    'Color','m','LineWidth',4)
for( i = 10:10:avg_anyAngle_N )
    line( [ temp_Long_List(i,2:3) ],... %X longitude arcdeg
    [ temp_Lat_List(i,2:3) ],... %Y latitude arcdeg
    'Color','m','LineWidth',1)
end
hold off;


% time_par = parallel.pool.Constant(time( keepr ));
% sTEC_par = parallel.pool.Constant(sTEC( keepr ));
% pplat_par = parallel.pool.Constant(pplat( keepr ));
% pplong_par = parallel.pool.Constant(pplong( keepr ));
% temp_Long_List_par = parallel.pool.Constant(temp_Long_List);
% temp_Lat_List_par = parallel.pool.Constant(temp_Lat_List);
% time_par =  time(keepr);
% time_par = time;
% clear time
% sTEC_par =  sTEC(keepr);
% sTEC_par  = sTEC;
% clear sTEC
% pplat_par = pplat(keepr);
% pplat_par = pplat;
% clear pplat
% pplong_par = pplong(keepr);
% pplong_par = pplong;
% clear pplong keepr
if( FLG_memSavr == 1 )
    % keepr = pplong < max(max(temp_Long_List)) &  pplong > min(min(temp_Long_List)) & pplat < max(max(temp_Lat_List)) & pplat > min(min(temp_Lat_List)); %limit memory sent to parallel workers
    keepr = inpolygon(pplong,pplat, [avg_anyAngle_Range(1:2,2) ; flipud(avg_anyAngle_Range(3:4,2)) ; avg_anyAngle_Range(1,2)] , [avg_anyAngle_Range(1:2,1) ; flipud(avg_anyAngle_Range(3:4,1)) ; avg_anyAngle_Range(1,1)]);
    %advanced limiter to greatly reduce memory sent to parallel workers based on if data is in polygon of avging 
    time_par =  time(keepr);
    clear time
    sTEC_par =  sTEC(keepr);
    clear sTEC
    pplat_par = pplat(keepr);
    clear pplat
    pplong_par = pplong(keepr);
    clear pplong keepr
elseif( FLG_memSavr == 2) %special mode where keepr is skipped because it is the whole world - must choose specifically
    time_par = time;
    clear time
    sTEC_par = sTEC;
    clear sTEC
    pplat_par = pplat;
    clear pplat
    pplong_par = pplong;
    clear pplong
else
    % keepr = pplong < max(max(temp_Long_List)) &  pplong > min(min(temp_Long_List)) & pplat < max(max(temp_Lat_List)) & pplat > min(min(temp_Lat_List)); %limit memory sent to parallel workers
    keepr = inpolygon(pplong,pplat, [avg_anyAngle_Range(1:2,2) ; flipud(avg_anyAngle_Range(3:4,2)) ; avg_anyAngle_Range(1,2)] , [avg_anyAngle_Range(1:2,1) ; flipud(avg_anyAngle_Range(3:4,1)) ; avg_anyAngle_Range(1,1)]);
    %advanced limiter to greatly reduce memory sent to parallel workers based on if data is in polygon of avging 
    time_par =  time(keepr);
    sTEC_par =  sTEC(keepr);
    pplat_par = pplat(keepr);
    pplong_par = pplong(keepr);
end

%Used to evaluate if there's space for parallel stuff or not
s_time = whos('time_par'); %get the sizes of the variables
s_sTEC = whos('sTEC_par');
s_pplat = whos('pplat_par');
s_pplong = whos('pplong_par');
memSizePredicted = (s_time.bytes/(1024^3) + s_sTEC.bytes/(1024^3) + s_pplat.bytes/(1024^3) + s_pplong.bytes/(1024^3) )*(CPU_Threads)*1.4;
%CPU threads and *1.4 for safety because does it go cray
%     memSizeCurrent = (s_time.bytes/(1024^3) + s_sTEC.bytes/(1024^3) + s_pplat.bytes/(1024^3) + s_pplong.bytes/(1024^3) );

[~,systemview] = memory; %checks the system memory for sizes and stuff
%     systemview.PhysicalMemory.Available/(1024^3)
%     systemview.PhysicalMemory.Total/(1024^3)
memSizeAvail = systemview.PhysicalMemory.Available/(1024^3);
clear s_time s_sTEC s_pplat s_pplong systemview %keep excess vars down

%did not work, didn't feel like figuring out why - "using exist to check for a variable" error yet not using exist
if( memSizeAvail > memSizePredicted )
    %parallel can work and will do it - this keeps copies from moving around too much (can't avoid copy I HAVE TRIED MANY TIMES)
    %SPMD may avoid copy but I ain't
    time_par  = parallel.pool.Constant(time_par); %copies to a paralell constant
    sTEC_par  = parallel.pool.Constant(sTEC_par); %copies to a paralell constant
    pplat_par  = parallel.pool.Constant(pplat_par); %copies to a paralell constant
    pplong_par = parallel.pool.Constant(pplong_par); %copies to a paralell constant
end
% %no turning back, will break things and can't get data outa this

if( FLG_AVG_ANYANGLE_FULLTIME == 1 ) %RUNS FOR THE FULL TIME BREAKS ASSUMPTIONS LATERRRRRRRRR
    %NO CHECK FOR Zenith OR MISA BECAUSE FULL TIME DOESN'T CARE
    sTECChunked_anyAngleAvg = zeros(length(timeUnique),avg_anyAngle_N); %preallocate
    % finAnnounce_cntr = 1; %prep for progress reporter
    finAnnounce_percentToUpdateAt = 10; % every % to update info at
    finAnnounce_div = round(100/finAnnounce_percentToUpdateAt); %calc divisor to use

    if( memSizeAvail < memSizePredicted )
        %if memory available isn't big enough for the dataset, don't parallel
        for( i = 1:length(timeUnique) )
            %Corral the data to the right place  
            k = time_par == timeUnique(i); %gets during a time period

            temp_sTEC = sTEC_par(k); %record the sTEC for the time
            temp_pplat = pplat_par(k); %record the pierce-point latitude
            temp_pplong = pplong_par(k); %record the pierce-point longitude

        %     kk = temp_pplat >= min(pplatRange) & temp_pplat <= max(pplatRange); %get the range of long to keep
        %     temp_sTEC = temp_sTEC(kk); %these get cut too to match the counting
        %     temp_pplong = temp_pplong(kk); %these get cut too to match the counting
        %     temp_pplat = temp_pplat(kk); %cut out only the lat wanted
            for(j = 1: avg_anyAngle_N )
                %average sTEC for a range chunk on an angle
        %         tLong = temp_Long_List_par(j,:);
        %         tLat = temp_Lat_List_par(j,:);
                kl = inpolygon(temp_pplong,temp_pplat,temp_Long_List(j,:),temp_Lat_List(j,:)); %get pts inside the area defined
                %Gets pts inside the range

        %         figure %Plots bits in and bits out
        %         plot([temp_Longs_up(j,:),temp_Longs_lo(j,:),temp_Longs_up(j,1)],[temp_Lats_up(j,:),temp_Lats_lo(j,:),temp_Lats_up(j,1)]) % polygon
        %         axis equal
        %         hold on
        %         plot(temp_pplong(kl),temp_pplat(kl),'r+') % points inside
        %         plot(temp_pplong(~kl),temp_pplat(~kl),'bo') % points outside
        %         hold off

                sTECChunked_anyAngleAvg(i,j) = mean(temp_sTEC(kl)); %average the sTEC in this lat band with the given long range
            end      

            %Generic progress reporter
            if( mod(i,round(length(timeUnique)/finAnnounce_div)) == 0 )
                fprintf('PROGRESS: Averaging TEC\t%%Complete: %.1f%%\tRuntime: %.3f sec\n',finAnnounce_percentToUpdateAt*i/660,toc); %announce job, % complete, and time so far
            end

        end %end of parallel loop
    else
        %PARALLEL TIME CAN FIT IT
        parfor( i = 1:length(timeUnique) )
            %Corral the data to the right place  
            k = time_par.Value == timeUnique(i); %gets during a time period

            temp_sTEC = sTEC_par.Value(k); %record the sTEC for the time
            temp_pplat = pplat_par.Value(k); %record the pierce-point latitude
            temp_pplong = pplong_par.Value(k); %record the pierce-point longitude

        %     kk = temp_pplat >= min(pplatRange) & temp_pplat <= max(pplatRange); %get the range of long to keep
        %     temp_sTEC = temp_sTEC(kk); %these get cut too to match the counting
        %     temp_pplong = temp_pplong(kk); %these get cut too to match the counting
        %     temp_pplat = temp_pplat(kk); %cut out only the lat wanted
            for(j = 1: avg_anyAngle_N )
                %average sTEC for a range chunk on an angle
        %         tLong = temp_Long_List_par(j,:);
        %         tLat = temp_Lat_List_par(j,:);
                kl = inpolygon(temp_pplong,temp_pplat,temp_Long_List(j,:),temp_Lat_List(j,:)); %get pts inside the area defined
                %Gets pts inside the range

        %         figure %Plots bits in and bits out
        %         plot([temp_Longs_up(j,:),temp_Longs_lo(j,:),temp_Longs_up(j,1)],[temp_Lats_up(j,:),temp_Lats_lo(j,:),temp_Lats_up(j,1)]) % polygon
        %         axis equal
        %         hold on
        %         plot(temp_pplong(kl),temp_pplat(kl),'r+') % points inside
        %         plot(temp_pplong(~kl),temp_pplat(~kl),'bo') % points outside
        %         hold off

                sTECChunked_anyAngleAvg(i,j) = mean(temp_sTEC(kl)); %average the sTEC in this lat band with the given long range
            end      

        end %end of parallel loop
    end %end of parallel or no parallel
    
else % ELSE MATCH TIME TO Zenith OR MISA
    
    if( FLG_AVG_ANYANGLE_ZENITHORMISA == 0 ) %if 0, Zenith ISR is plotted
        sTECChunked_anyAngleAvg = zeros(length(timeUnique_limISR_Zenith),avg_anyAngle_N); %preallocate
        % finAnnounce_cntr = 1; %prep for progress reporter
        finAnnounce_percentToUpdateAt = 10; % every % to update info at
        finAnnounce_div = round(100/finAnnounce_percentToUpdateAt); %calc divisor to use

        if( memSizeAvail < memSizePredicted )
            %if memory available isn't big enough for the dataset, don't parallel
            for( i = 1:length(timeUnique_limISR_Zenith) )
                %Corral the data to the right place  
                k = time_par == timeUnique_limISR_Zenith(i); %gets during a time period

                temp_sTEC = sTEC_par(k); %record the sTEC for the time
                temp_pplat = pplat_par(k); %record the pierce-point latitude
                temp_pplong = pplong_par(k); %record the pierce-point longitude

            %     kk = temp_pplat >= min(pplatRange) & temp_pplat <= max(pplatRange); %get the range of long to keep
            %     temp_sTEC = temp_sTEC(kk); %these get cut too to match the counting
            %     temp_pplong = temp_pplong(kk); %these get cut too to match the counting
            %     temp_pplat = temp_pplat(kk); %cut out only the lat wanted
                for(j = 1: avg_anyAngle_N )
                    %average sTEC for a range chunk on an angle
            %         tLong = temp_Long_List_par(j,:);
            %         tLat = temp_Lat_List_par(j,:);
                    kl = inpolygon(temp_pplong,temp_pplat,temp_Long_List(j,:),temp_Lat_List(j,:)); %get pts inside the area defined
                    %Gets pts inside the range

            %         figure %Plots bits in and bits out
            %         plot([temp_Longs_up(j,:),temp_Longs_lo(j,:),temp_Longs_up(j,1)],[temp_Lats_up(j,:),temp_Lats_lo(j,:),temp_Lats_up(j,1)]) % polygon
            %         axis equal
            %         hold on
            %         plot(temp_pplong(kl),temp_pplat(kl),'r+') % points inside
            %         plot(temp_pplong(~kl),temp_pplat(~kl),'bo') % points outside
            %         hold off

                    sTECChunked_anyAngleAvg(i,j) = mean(temp_sTEC(kl)); %average the sTEC in this lat band with the given long range
                end      

                %Generic progress reporter
                if( mod(i,round(length(timeUnique_limISR_Zenith)/finAnnounce_div)) == 0 )
                    fprintf('PROGRESS: Averaging TEC\t%%Complete: %.1f%%\tRuntime: %.3f sec\n',finAnnounce_percentToUpdateAt*i/660,toc); %announce job, % complete, and time so far
                end

            end %end of parallel loop
        else
            %PARALLEL TIME CAN FIT IT
            parfor( i = 1:length(timeUnique_limISR_Zenith) )
                %Corral the data to the right place  
                k = time_par.Value == timeUnique_limISR_Zenith(i); %gets during a time period

                temp_sTEC = sTEC_par.Value(k); %record the sTEC for the time
                temp_pplat = pplat_par.Value(k); %record the pierce-point latitude
                temp_pplong = pplong_par.Value(k); %record the pierce-point longitude

            %     kk = temp_pplat >= min(pplatRange) & temp_pplat <= max(pplatRange); %get the range of long to keep
            %     temp_sTEC = temp_sTEC(kk); %these get cut too to match the counting
            %     temp_pplong = temp_pplong(kk); %these get cut too to match the counting
            %     temp_pplat = temp_pplat(kk); %cut out only the lat wanted
                for(j = 1: avg_anyAngle_N )
                    %average sTEC for a range chunk on an angle
            %         tLong = temp_Long_List_par(j,:);
            %         tLat = temp_Lat_List_par(j,:);
                    kl = inpolygon(temp_pplong,temp_pplat,temp_Long_List(j,:),temp_Lat_List(j,:)); %get pts inside the area defined
                    %Gets pts inside the range

            %         figure %Plots bits in and bits out
            %         plot([temp_Longs_up(j,:),temp_Longs_lo(j,:),temp_Longs_up(j,1)],[temp_Lats_up(j,:),temp_Lats_lo(j,:),temp_Lats_up(j,1)]) % polygon
            %         axis equal
            %         hold on
            %         plot(temp_pplong(kl),temp_pplat(kl),'r+') % points inside
            %         plot(temp_pplong(~kl),temp_pplat(~kl),'bo') % points outside
            %         hold off

                    sTECChunked_anyAngleAvg(i,j) = mean(temp_sTEC(kl)); %average the sTEC in this lat band with the given long range
                end      

            end %end of parallel loop
        end %end of parallel or no parallel
    else %if 1, MISA ISR is plotted
        sTECChunked_anyAngleAvg = zeros(length(timeUnique_limISR_MISA),avg_anyAngle_N); %preallocate
        % finAnnounce_cntr = 1; %prep for progress reporter
        finAnnounce_percentToUpdateAt = 10; % every % to update info at
        finAnnounce_div = round(100/finAnnounce_percentToUpdateAt); %calc divisor to use

        %CAN BE PAR - DISABLE FOR MORE MEM RANGE
        if( memSizeAvail < memSizePredicted )
            %if memory available isn't big enough for the dataset, don't parallel
            for( i = 1:length(timeUnique_limISR_MISA) )
                %Corral the data to the right place  
                k = time_par == timeUnique_limISR_MISA(i); %gets during a time period

                temp_sTEC = sTEC_par(k); %record the sTEC for the time
                temp_pplat = pplat_par(k); %record the pierce-point latitude
                temp_pplong = pplong_par(k); %record the pierce-point longitude

            %     kk = temp_pplat >= min(pplatRange) & temp_pplat <= max(pplatRange); %get the range of long to keep
            %     temp_sTEC = temp_sTEC(kk); %these get cut too to match the counting
            %     temp_pplong = temp_pplong(kk); %these get cut too to match the counting
            %     temp_pplat = temp_pplat(kk); %cut out only the lat wanted
                for(j = 1: avg_anyAngle_N )
                    %average sTEC for a range chunk on an angle
            %         tLong = temp_Long_List_par(j,:);
            %         tLat = temp_Lat_List_par(j,:);
                    kl = inpolygon(temp_pplong,temp_pplat,temp_Long_List(j,:),temp_Lat_List(j,:)); %get pts inside the area defined
                    %Gets pts inside the range

            %         figure %Plots bits in and bits out
            %         plot([temp_Longs_up(j,:),temp_Longs_lo(j,:),temp_Longs_up(j,1)],[temp_Lats_up(j,:),temp_Lats_lo(j,:),temp_Lats_up(j,1)]) % polygon
            %         axis equal
            %         hold on
            %         plot(temp_pplong(kl),temp_pplat(kl),'r+') % points inside
            %         plot(temp_pplong(~kl),temp_pplat(~kl),'bo') % points outside
            %         hold off

                    sTECChunked_anyAngleAvg(i,j) = mean(temp_sTEC(kl)); %average the sTEC in this lat band with the given long range
                end      

                %Generic progress reporter
                if( mod(i,round(length(timeUnique_limISR_MISA)/finAnnounce_div)) == 0 )
                    fprintf('PROGRESS: Averaging TEC\t%%Complete: %.1f%%\tRuntime: %.3f sec\n',finAnnounce_percentToUpdateAt*i/(length(timeUnique_limISR_MISA)/10),toc); %announce job, % complete, and time so far
                end

            end
        else %CAN PARALLEL LETS DO THIS
            parfor( i = 1:length(timeUnique_limISR_MISA) )
                %Corral the data to the right place  
                k = time_par.Value == timeUnique_limISR_MISA(i); %gets during a time period

                temp_sTEC = sTEC_par.Value(k); %record the sTEC for the time
                temp_pplat = pplat_par.Value(k); %record the pierce-point latitude
                temp_pplong = pplong_par.Value(k); %record the pierce-point longitude

            %     kk = temp_pplat >= min(pplatRange) & temp_pplat <= max(pplatRange); %get the range of long to keep
            %     temp_sTEC = temp_sTEC(kk); %these get cut too to match the counting
            %     temp_pplong = temp_pplong(kk); %these get cut too to match the counting
            %     temp_pplat = temp_pplat(kk); %cut out only the lat wanted
                for(j = 1: avg_anyAngle_N )
                    %average sTEC for a range chunk on an angle
            %         tLong = temp_Long_List_par(j,:);
            %         tLat = temp_Lat_List_par(j,:);
                    kl = inpolygon(temp_pplong,temp_pplat,temp_Long_List(j,:),temp_Lat_List(j,:)); %get pts inside the area defined
                    %Gets pts inside the range

            %         figure %Plots bits in and bits out
            %         plot([temp_Longs_up(j,:),temp_Longs_lo(j,:),temp_Longs_up(j,1)],[temp_Lats_up(j,:),temp_Lats_lo(j,:),temp_Lats_up(j,1)]) % polygon
            %         axis equal
            %         hold on
            %         plot(temp_pplong(kl),temp_pplat(kl),'r+') % points inside
            %         plot(temp_pplong(~kl),temp_pplat(~kl),'bo') % points outside
            %         hold off

                    sTECChunked_anyAngleAvg(i,j) = mean(temp_sTEC(kl)); %average the sTEC in this lat band with the given long range
                end      

            end %end of parallel loop
        end %end of parallel or no parallel check
    end

end %END OF FULL TIME OR MATCH TO Zenith OR MISA

if( memSizeAvail > memSizePredicted )
    fprintf('PROGRESS for PARALLEL: Averaging TEC\t%%Complete: %%100\tRuntime: %.3f sec\n',toc);
end

clear temp_sTEC temp_pplat temp_pplong temp_Long_List temp_Lat_List keepr temp_Longs_up temp_Longs_lo temp_Lats_up temp_Lats_lo
toc

end %*********END OF FLG_AVG_ANYANGL*********


%% Average Any Angle, PLOT IT UP
if( FLG_AVG_ANYANGLE_PLOT == 1 )
    
    if( FLG_AVG_ANYANGLE_FULLTIME == 1 )
        
        %Plot just the TEC
        %prep figure #
        % figure('units','normalized','outerposition',[0 0 1 1]); %maximizes figure window
        figure(figNum);
        figNum = figNum+1;
        h = pcolor( (timeUnique - dateRange_zeroHr(2)).*24 , avg_anyAngle_Range_Chunks_Long_Plot , sTECChunked_anyAngleAvg');% pseudocolor plot "stretched" to the grid
        set(h, 'EdgeColor', 'none'); %remove grid lines (which make it look hilarious)
        %             shading interp; %stuff for interpolation
        %             lighting phong; %stuff for interpolation
        caxis([-gif_sTEC_cap, gif_sTEC_cap]);
        %         colormap(flipud(colormap('hot'))) %heat map where 0 (low) is white
        colormap('jet');
        h2 = colorbar; %shows the color bar
        ylabel(h2, 'delta-sTEC (TECU)','fontweight',font_Weight, 'FontSize',font_Size); %labels color bar
        Xaxisvar = [ (round((min(timeUnique)-dateRange_zeroHr(2))*24) - mod(round((min(timeUnique)-dateRange_zeroHr(2))*24),2)) :2:...
                (round((max(timeUnique)-dateRange_zeroHr(2))*24) - mod(round((max(timeUnique)-dateRange_zeroHr(2))*24),2)) ];
        set(gca,'TickDir','out'); %put tick marks on outside of grid to help connect to time elements
        set(gca,'box','off'); %keeps tick marks from being on top and right side too
        set(gca,'xtick',Xaxisvar);

        avg_anyAngle_Range_Chunks_Long_Plot_autoTick = (ceil(max(avg_anyAngle_Range_Chunks_Long_Plot)) - floor(min(avg_anyAngle_Range_Chunks_Long_Plot)))/13; %tries to split the latitude range into 13 parts (based off of 180/15+1)
        if( avg_anyAngle_Range_Chunks_Long_Plot_autoTick > 25 )
            avg_anyAngle_Range_Chunks_Long_Plot_autoTick = 30; %sets the tick setting to 15 arcdegrees per tick
        elseif( avg_anyAngle_Range_Chunks_Long_Plot_autoTick > 10 )
            avg_anyAngle_Range_Chunks_Long_Plot_autoTick = 15; %sets the tick setting to 15 arcdegrees per tick
        elseif( avg_anyAngle_Range_Chunks_Long_Plot_autoTick > 5 )
            avg_anyAngle_Range_Chunks_Long_Plot_autoTick = 10; %sets the tick setting to 10 arcdegrees per tick
        elseif( avg_anyAngle_Range_Chunks_Long_Plot_autoTick > 2 )
            avg_anyAngle_Range_Chunks_Long_Plot_autoTick = 5; %sets the tick setting to 5 arcdegrees per tick
        elseif( avg_anyAngle_Range_Chunks_Long_Plot_autoTick > 1 )
            avg_anyAngle_Range_Chunks_Long_Plot_autoTick = 2; %sets the tick setting to 5 arcdegrees per tick
        elseif( avg_anyAngle_Range_Chunks_Long_Plot_autoTick >= 0.6 ) %0.6 because 15/25 = 0.6, so there will be enough 1 arcdeg ticks
            avg_anyAngle_Range_Chunks_Long_Plot_autoTick = 1; %sets the tick setting to 1 arcdegree per tick
        else
            avg_anyAngle_Range_Chunks_Long_Plot_autoTick = (max(plotLongRange) - min(plotLongRange))/13; %just goes for it if it's a super tiny range
        end
        Yaxisvar = round(floor(min(avg_anyAngle_Range_Chunks_Long_Plot)):avg_anyAngle_Range_Chunks_Long_Plot_autoTick:ceil(max(avg_anyAngle_Range_Chunks_Long_Plot)),2); %creates y ticks automagically
        set(gca,'ytick',Yaxisvar);

        hold on;
        %Now drawing line of interest
        if( strcmp(avg_anyAngle_Range_Chunks_Long_Plot_Name,'Longitude') == 1 ) %if true, longitude
            line([ min((timeUnique_limISR_Zenith -floor(mean(timeUnique))).*24) , max((timeUnique_limISR_Zenith -floor(mean(timeUnique))).*24) ],[longMillstone,longMillstone],'Color','k','LineWidth',.25);  %plots a point with a black line
        else %else latitude
            line([ min((timeUnique_limISR_Zenith -floor(mean(timeUnique))).*24) , max((timeUnique_limISR_Zenith -floor(mean(timeUnique))).*24) ],[latMillstone,latMillstone],'Color','k','LineWidth',.25);  %plots a point with a black line
        end

        xlabel(['Time in UT - 0 Hr on Day ',num2str(dateRange_zeroHr(2)),', ',num2str(dateRange_zeroHr(1)),' (hr)']); %old: ,'fontweight','bold','FontSize',12 for all
        ylabel([avg_anyAngle_Range_Chunks_Long_Plot_Name,' (arcdeg)']); %sped up significantly when merged above
        string_Title = ['TEC Averaged on Angle of ',num2str(round(avg_anyAngle*180/pi,2)),' deg and Width of ',num2str(round(avg_anyAngle_Width,2)),' arcdeg, Line Shows ',avg_anyAngle_Range_Chunks_Long_Plot_Name,' of Millstone Hill Zenith Beam']; %create title string
        if( FLG_geomagneticCoords == 1)
            string_Title = [string_Title,' [Mag Coord Aligned]']; %tack on mag coordinate warning
        end
        title(string_Title);
        set(gca,'fontweight',font_Weight, 'FontSize',font_Size);
        hold off;
        
    else %OTHERWISE PLOT VS Zenith/MISA

        if( FLG_AVG_ANYANGLE_ZENITHORMISA == 0 ) %Zenith option
            %Plot just the TEC
            %prep figure #
            % figure('units','normalized','outerposition',[0 0 1 1]); %maximizes figure window
            figure(figNum);
            figNum = figNum+1;
            h = pcolor( (timeUnique_limISR_Zenith - dateRange_zeroHr(2)).*24 , avg_anyAngle_Range_Chunks_Long_Plot , sTECChunked_anyAngleAvg');% pseudocolor plot "stretched" to the grid
            set(h, 'EdgeColor', 'none'); %remove grid lines (which make it look hilarious)
            %             shading interp; %stuff for interpolation
            %             lighting phong; %stuff for interpolation
            caxis([-gif_sTEC_cap, gif_sTEC_cap]);
            %         colormap(flipud(colormap('hot'))) %heat map where 0 (low) is white
            colormap('jet');
            h2 = colorbar; %shows the color bar
            ylabel(h2, 'delta-sTEC (TECU)','fontweight',font_Weight, 'FontSize',font_Size); %labels color bar
            Xaxisvar = [ (round((min(Zenith_time)-floor(mean(Zenith_time)))*24)-mod(round((min(Zenith_time)-floor(mean(Zenith_time)))*24),2)):2:...
                (round((max(Zenith_time)-floor(mean(Zenith_time)))*24)-mod(round((max(Zenith_time)-floor(mean(Zenith_time)))*24),2)) ];
            set(gca,'TickDir','out'); %put tick marks on outside of grid to help connect to time elements
            set(gca,'box','off'); %keeps tick marks from being on top and right side too
            set(gca,'xtick',Xaxisvar);

            avg_anyAngle_Range_Chunks_Long_Plot_autoTick = (ceil(max(avg_anyAngle_Range_Chunks_Long_Plot)) - floor(min(avg_anyAngle_Range_Chunks_Long_Plot)))/13; %tries to split the latitude range into 13 parts (based off of 180/15+1)
            if( avg_anyAngle_Range_Chunks_Long_Plot_autoTick > 25 )
                avg_anyAngle_Range_Chunks_Long_Plot_autoTick = 30; %sets the tick setting to 15 arcdegrees per tick
            elseif( avg_anyAngle_Range_Chunks_Long_Plot_autoTick > 10 )
                avg_anyAngle_Range_Chunks_Long_Plot_autoTick = 15; %sets the tick setting to 15 arcdegrees per tick
            elseif( avg_anyAngle_Range_Chunks_Long_Plot_autoTick > 5 )
                avg_anyAngle_Range_Chunks_Long_Plot_autoTick = 10; %sets the tick setting to 10 arcdegrees per tick
            elseif( avg_anyAngle_Range_Chunks_Long_Plot_autoTick > 2 )
                avg_anyAngle_Range_Chunks_Long_Plot_autoTick = 5; %sets the tick setting to 5 arcdegrees per tick
            elseif( avg_anyAngle_Range_Chunks_Long_Plot_autoTick > 1 )
                avg_anyAngle_Range_Chunks_Long_Plot_autoTick = 2; %sets the tick setting to 5 arcdegrees per tick
            elseif( avg_anyAngle_Range_Chunks_Long_Plot_autoTick >= 0.6 ) %0.6 because 15/25 = 0.6, so there will be enough 1 arcdeg ticks
                avg_anyAngle_Range_Chunks_Long_Plot_autoTick = 1; %sets the tick setting to 1 arcdegree per tick
            else
                avg_anyAngle_Range_Chunks_Long_Plot_autoTick = (max(plotLongRange) - min(plotLongRange))/13; %just goes for it if it's a super tiny range
            end
            Yaxisvar = round(floor(min(avg_anyAngle_Range_Chunks_Long_Plot)):avg_anyAngle_Range_Chunks_Long_Plot_autoTick:ceil(max(avg_anyAngle_Range_Chunks_Long_Plot)),2); %creates y ticks automagically
            set(gca,'ytick',Yaxisvar);

            hold on;
            %Now drawing line of interest
            if( strcmp(avg_anyAngle_Range_Chunks_Long_Plot_Name,'Longitude') == 1 ) %if true, longitude
                line([ min((timeUnique_limISR_Zenith -floor(mean(timeUnique))).*24) , max((timeUnique_limISR_Zenith -floor(mean(timeUnique))).*24) ],[longMillstone,longMillstone],'Color','k','LineWidth',.25);  %plots a point with a black line
            else %else latitude
                line([ min((timeUnique_limISR_Zenith -floor(mean(timeUnique))).*24) , max((timeUnique_limISR_Zenith -floor(mean(timeUnique))).*24) ],[latMillstone,latMillstone],'Color','k','LineWidth',.25);  %plots a point with a black line
            end

            xlabel(['Time in UT - 0 Hr on Day ',num2str(dateRange_zeroHr(2)),', ',num2str(dateRange_zeroHr(1)),' (hr)']); %old: ,'fontweight','bold','FontSize',12 for all
            ylabel([avg_anyAngle_Range_Chunks_Long_Plot_Name,' (arcdeg)']); %sped up significantly when merged above
            string_Title = ['TEC Averaged on Angle of ',num2str(round(avg_anyAngle*180/pi,2)),' deg and Width of ',num2str(round(avg_anyAngle_Width,2)),' arcdeg, Line Shows ',avg_anyAngle_Range_Chunks_Long_Plot_Name,' of Millstone Hill Zenith Beam']; %create title string
            if( FLG_geomagneticCoords == 1)
                string_Title = [string_Title,' [Mag Coord Aligned]']; %tack on mag coordinate warning
            end
            title(string_Title);
            set(gca,'fontweight',font_Weight, 'FontSize',font_Size);
            hold off;

            %PLOT WITH Zenith
            %prep figure #
            % figure('units','normalized','outerposition',[0 0 1 1]); %maximizes figure window
            figure(figNum);
            figNum = figNum+1;
            %Prep plot space
            ax1 = subplot(2,1,1); %upper bit
            h1 = pcolor( (Zenith_time - dateRange_zeroHr(2)).*24 ,Zenith_height,Zenith_SNR_bp');
            set(h1, 'EdgeColor', 'none'); %remove grid lines (which make it look hilarious)
            shading interp;
            lighting phong;
            colorbar;
            caxis([-0.1 0.1]);
            colormap(ax1,'gray');
            ylim([90 700]);
            % ylim([90 500]);
            % xlim([0 70]);
    %         xlabel('Time in UT - 0 Hr at 0 UT May 7th');
            ylabel('Height (km)');
            title(['May 6-7-8 2013, Zenith SNR 2 Hr Highpass Filtered, Latitude = ',num2str(round(latMillstone,2)),' arcdeg and Longitude = ',num2str(round(longMillstone,2)),' arcdeg']);
            % set(gca, 'fontweight','bold', 'FontSize',12);
            % Xaxisvar = [-12, -8, -4, 0, 4, 8, 12, 16, 20, 24 28, 32, 36, 40, 44, 48]; %shows off important magnitude points
            %         Xaxisvar = [-12,-8,-4,0, 4, 8, 12, 16, 20, 24 28, 32, 36, 40, 44, 48, 52, 56, 60, 64, 68, 72]; %shows off important magnitude points
            Xaxisvar = [ (round((min(Zenith_time)-floor(mean(Zenith_time)))*24)-mod(round((min(Zenith_time)-floor(mean(Zenith_time)))*24),2)):2:...
                (round((max(Zenith_time)-floor(mean(Zenith_time)))*24)-mod(round((max(Zenith_time)-floor(mean(Zenith_time)))*24),2)) ];
            set(gca,'TickDir','out'); %put tick marks on outside of grid to help connect to time elements
            set(gca,'box','off'); %keeps tick marks from being on top and right side too
            set(gca,'xtick',Xaxisvar);
            set(gca,'fontweight',font_Weight, 'FontSize',font_Size);
            % 
            % %add on a vertical line that shows where we are at in time
            % 
            ax2 = subplot(2,1,2); %lower bit
            h = pcolor( (timeUnique_limISR_Zenith - dateRange_zeroHr(2)).*24 , avg_anyAngle_Range_Chunks_Long_Plot , sTECChunked_anyAngleAvg');% pseudocolor plot "stretched" to the grid
            set(h, 'EdgeColor', 'none'); %remove grid lines (which make it look hilarious)
            %             shading interp; %stuff for interpolation
            %             lighting phong; %stuff for interpolation
            caxis([-gif_sTEC_cap, gif_sTEC_cap]);
            %         colormap(flipud(colormap('hot'))) %heat map where 0 (low) is white
            colormap(ax2,'jet');
            h2 = colorbar; %shows the color bar
            ylabel(h2, 'delta-sTEC (TECU)','fontweight',font_Weight, 'FontSize',font_Size); %labels color bar
            set(gca,'TickDir','out'); %put tick marks on outside of grid to help connect to time elements
            set(gca,'box','off'); %keeps tick marks from being on top and right side too
            set(gca,'xtick',Xaxisvar);

            avg_anyAngle_Range_Chunks_Long_Plot_autoTick = (ceil(max(avg_anyAngle_Range_Chunks_Long_Plot)) - floor(min(avg_anyAngle_Range_Chunks_Long_Plot)))/13; %tries to split the latitude range into 13 parts (based off of 180/15+1)
            if( avg_anyAngle_Range_Chunks_Long_Plot_autoTick > 25 )
                avg_anyAngle_Range_Chunks_Long_Plot_autoTick = 30; %sets the tick setting to 15 arcdegrees per tick
            elseif( avg_anyAngle_Range_Chunks_Long_Plot_autoTick > 10 )
                avg_anyAngle_Range_Chunks_Long_Plot_autoTick = 15; %sets the tick setting to 15 arcdegrees per tick
            elseif( avg_anyAngle_Range_Chunks_Long_Plot_autoTick > 5 )
                avg_anyAngle_Range_Chunks_Long_Plot_autoTick = 10; %sets the tick setting to 10 arcdegrees per tick
            elseif( avg_anyAngle_Range_Chunks_Long_Plot_autoTick > 2 )
                avg_anyAngle_Range_Chunks_Long_Plot_autoTick = 5; %sets the tick setting to 5 arcdegrees per tick
            elseif( avg_anyAngle_Range_Chunks_Long_Plot_autoTick > 1 )
                avg_anyAngle_Range_Chunks_Long_Plot_autoTick = 2; %sets the tick setting to 5 arcdegrees per tick
            elseif( avg_anyAngle_Range_Chunks_Long_Plot_autoTick >= 0.6 ) %0.6 because 15/25 = 0.6, so there will be enough 1 arcdeg ticks
                avg_anyAngle_Range_Chunks_Long_Plot_autoTick = 1; %sets the tick setting to 1 arcdegree per tick
            else
                avg_anyAngle_Range_Chunks_Long_Plot_autoTick = (max(plotLongRange) - min(plotLongRange))/13; %just goes for it if it's a super tiny range
            end
            Yaxisvar = round(floor(min(avg_anyAngle_Range_Chunks_Long_Plot)):avg_anyAngle_Range_Chunks_Long_Plot_autoTick:ceil(max(avg_anyAngle_Range_Chunks_Long_Plot)),2); %creates y ticks automagically
            set(gca,'ytick',Yaxisvar);

            hold on;
            %Now drawing line of interest
            if( strcmp(avg_anyAngle_Range_Chunks_Long_Plot_Name,'Longitude') == 1 ) %if true, longitude
                line([ min((timeUnique_limISR_Zenith -floor(mean(timeUnique))).*24) , max((timeUnique_limISR_Zenith -floor(mean(timeUnique))).*24) ],[longMillstone,longMillstone],'Color','k','LineWidth',.25);  %plots a point with a black line
            else %else latitude
                line([ min((timeUnique_limISR_Zenith -floor(mean(timeUnique))).*24) , max((timeUnique_limISR_Zenith -floor(mean(timeUnique))).*24) ],[latMillstone,latMillstone],'Color','k','LineWidth',.25);  %plots a point with a black line
            end

            xlabel(['Time in UT - 0 Hr on Day ',num2str(dateRange_zeroHr(2)),', ',num2str(dateRange_zeroHr(1)),' (hr)']); %old: ,'fontweight','bold','FontSize',12 for all
            ylabel([avg_anyAngle_Range_Chunks_Long_Plot_Name,' (arcdeg)']); %sped up significantly when merged above
            string_Title = ['TEC Averaged on Angle of ',num2str(round(avg_anyAngle*180/pi,2)),' deg and Width of ',num2str(round(avg_anyAngle_Width,2)),' arcdeg, Line Shows ',avg_anyAngle_Range_Chunks_Long_Plot_Name,' of Millstone Hill Zenith Beam']; %create title string
            if( FLG_geomagneticCoords == 1)
                string_Title = [string_Title,' [Mag Coord Aligned]']; %tack on mag coordinate warning
            end
            title(string_Title);
            set(gca,'fontweight',font_Weight, 'FontSize',font_Size);
            hold off;


            %prep figure #
            % figure('units','normalized','outerposition',[0 0 1 1]); %maximizes figure window
            figure(figNum);
            figNum = figNum+1;
            %Prep plot space
            ax1 = subplot(2,1,1); %upper bit
            h = pcolor( (timeUnique_limISR_Zenith - dateRange_zeroHr(2)).*24 , avg_anyAngle_Range_Chunks_Long_Plot , sTECChunked_anyAngleAvg');% pseudocolor plot "stretched" to the grid
            set(h, 'EdgeColor', 'none'); %remove grid lines (which make it look hilarious)
            %             shading interp; %stuff for interpolation
            %             lighting phong; %stuff for interpolation
            caxis([-gif_sTEC_cap, gif_sTEC_cap]);
            %         colormap(flipud(colormap('hot'))) %heat map where 0 (low) is white
            colormap(ax1,'jet');
            h2 = colorbar; %shows the color bar
            ylabel(h2, 'delta-sTEC (TECU)','fontweight',font_Weight, 'FontSize',font_Size); %labels color bar
            set(gca,'TickDir','out'); %put tick marks on outside of grid to help connect to time elements
            set(gca,'box','off'); %keeps tick marks from being on top and right side too
            set(gca,'xtick',Xaxisvar);

            avg_anyAngle_Range_Chunks_Long_Plot_autoTick = (ceil(max(avg_anyAngle_Range_Chunks_Long_Plot)) - floor(min(avg_anyAngle_Range_Chunks_Long_Plot)))/13; %tries to split the latitude range into 13 parts (based off of 180/15+1)
            if( avg_anyAngle_Range_Chunks_Long_Plot_autoTick > 25 )
                avg_anyAngle_Range_Chunks_Long_Plot_autoTick = 30; %sets the tick setting to 15 arcdegrees per tick
            elseif( avg_anyAngle_Range_Chunks_Long_Plot_autoTick > 10 )
                avg_anyAngle_Range_Chunks_Long_Plot_autoTick = 15; %sets the tick setting to 15 arcdegrees per tick
            elseif( avg_anyAngle_Range_Chunks_Long_Plot_autoTick > 5 )
                avg_anyAngle_Range_Chunks_Long_Plot_autoTick = 10; %sets the tick setting to 10 arcdegrees per tick
            elseif( avg_anyAngle_Range_Chunks_Long_Plot_autoTick > 2 )
                avg_anyAngle_Range_Chunks_Long_Plot_autoTick = 5; %sets the tick setting to 5 arcdegrees per tick
            elseif( avg_anyAngle_Range_Chunks_Long_Plot_autoTick > 1 )
                avg_anyAngle_Range_Chunks_Long_Plot_autoTick = 2; %sets the tick setting to 5 arcdegrees per tick
            elseif( avg_anyAngle_Range_Chunks_Long_Plot_autoTick >= 0.6 ) %0.6 because 15/25 = 0.6, so there will be enough 1 arcdeg ticks
                avg_anyAngle_Range_Chunks_Long_Plot_autoTick = 1; %sets the tick setting to 1 arcdegree per tick
            else
                avg_anyAngle_Range_Chunks_Long_Plot_autoTick = (max(plotLongRange) - min(plotLongRange))/13; %just goes for it if it's a super tiny range
            end
            Yaxisvar = round(floor(min(avg_anyAngle_Range_Chunks_Long_Plot)):avg_anyAngle_Range_Chunks_Long_Plot_autoTick:ceil(max(avg_anyAngle_Range_Chunks_Long_Plot)),2); %creates y ticks automagically
            set(gca,'ytick',Yaxisvar);

            hold on;
            %Now drawing line of interest
            if( strcmp(avg_anyAngle_Range_Chunks_Long_Plot_Name,'Longitude') == 1 ) %if true, longitude
                line([ min((timeUnique_limISR_Zenith -floor(mean(timeUnique))).*24) , max((timeUnique_limISR_Zenith -floor(mean(timeUnique))).*24) ],[longMillstone,longMillstone],'Color','k','LineWidth',.25);  %plots a point with a black line
            else %else latitude
                line([ min((timeUnique_limISR_Zenith -floor(mean(timeUnique))).*24) , max((timeUnique_limISR_Zenith -floor(mean(timeUnique))).*24) ],[latMillstone,latMillstone],'Color','k','LineWidth',.25);  %plots a point with a black line
            end

            ylabel([avg_anyAngle_Range_Chunks_Long_Plot_Name,' (arcdeg)']); %sped up significantly when merged above
            string_Title = ['TEC Averaged on Angle of ',num2str(round(avg_anyAngle*180/pi,2)),' deg and Width of ',num2str(round(avg_anyAngle_Width,2)),' arcdeg, Line Shows ',avg_anyAngle_Range_Chunks_Long_Plot_Name,' of Millstone Hill Zenith Beam']; %create title string
            if( FLG_geomagneticCoords == 1)
                string_Title = [string_Title,' [Mag Coord Aligned]']; %tack on mag coordinate warning
            end
            title(string_Title);
            set(gca,'fontweight',font_Weight, 'FontSize',font_Size);
            hold off;
            % 
            ax2 = subplot(2,1,2); %lower bit     
            h1 = pcolor( (Zenith_time - dateRange_zeroHr(2)).*24 ,Zenith_height,Zenith_SNR_bp');
            set(h1, 'EdgeColor', 'none'); %remove grid lines (which make it look hilarious)
            shading interp;
            lighting phong;
            colorbar;
            caxis([-0.1 0.1]);
            colormap(ax2,'gray');
            ylim([90 700]);
            % ylim([90 500]);
            % xlim([0 70]);
    %         xlabel('Time in UT - 0 Hr at 0 UT May 7th');
            xlabel(['Time in UT - 0 Hr on Day ',num2str(dateRange_zeroHr(2)),', ',num2str(dateRange_zeroHr(1)),' (hr)']); %old: ,'fontweight','bold','FontSize',12 for all
            ylabel('Height (km)');
            title(['May 6-7-8 2013, Zenith SNR 2 Hr Highpass Filtered, Latitude = ',num2str(round(latMillstone,2)),' arcdeg and Longitude = ',num2str(round(longMillstone,2)),' arcdeg']);
            % set(gca, 'fontweight','bold', 'FontSize',12);
            % Xaxisvar = [-12, -8, -4, 0, 4, 8, 12, 16, 20, 24 28, 32, 36, 40, 44, 48]; %shows off important magnitude points
            %         Xaxisvar = [-12,-8,-4,0, 4, 8, 12, 16, 20, 24 28, 32, 36, 40, 44, 48, 52, 56, 60, 64, 68, 72]; %shows off important magnitude points
            Xaxisvar = [ (round((min(Zenith_time)-floor(mean(Zenith_time)))*24)-mod(round((min(Zenith_time)-floor(mean(Zenith_time)))*24),2)):2:...
                (round((max(Zenith_time)-floor(mean(Zenith_time)))*24)-mod(round((max(Zenith_time)-floor(mean(Zenith_time)))*24),2)) ];
            set(gca,'TickDir','out'); %put tick marks on outside of grid to help connect to time elements
            set(gca,'box','off'); %keeps tick marks from being on top and right side too
            set(gca,'xtick',Xaxisvar);
            set(gca,'fontweight',font_Weight, 'FontSize',font_Size);

        else % if 1, MISA ISR is plotted

            %Plot just the TEC
            %prep figure #
            % figure('units','normalized','outerposition',[0 0 1 1]); %maximizes figure window
            figure(figNum);
            figNum = figNum+1;
            h = pcolor( (timeUnique_limISR_MISA - dateRange_zeroHr(2)).*24 , avg_anyAngle_Range_Chunks_Long_Plot , sTECChunked_anyAngleAvg');% pseudocolor plot "stretched" to the grid
            set(h, 'EdgeColor', 'none'); %remove grid lines (which make it look hilarious)
            %             shading interp; %stuff for interpolation
            %             lighting phong; %stuff for interpolation
            caxis([-gif_sTEC_cap, gif_sTEC_cap]);
            %         colormap(flipud(colormap('hot'))) %heat map where 0 (low) is white
            colormap('jet');
            h2 = colorbar; %shows the color bar
            ylabel(h2, 'delta-sTEC (TECU)','fontweight',font_Weight, 'FontSize',font_Size); %labels color bar
            Xaxisvar = [ (round((min(Zenith_time)-floor(mean(Zenith_time)))*24)-mod(round((min(Zenith_time)-floor(mean(Zenith_time)))*24),2)):2:...
                (round((max(Zenith_time)-floor(mean(Zenith_time)))*24)-mod(round((max(Zenith_time)-floor(mean(Zenith_time)))*24),2)) ];
            set(gca,'TickDir','out'); %put tick marks on outside of grid to help connect to time elements
            set(gca,'box','off'); %keeps tick marks from being on top and right side too
            set(gca,'xtick',Xaxisvar);

            avg_anyAngle_Range_Chunks_Long_Plot_autoTick = (ceil(max(avg_anyAngle_Range_Chunks_Long_Plot)) - floor(min(avg_anyAngle_Range_Chunks_Long_Plot)))/13; %tries to split the latitude range into 13 parts (based off of 180/15+1)
            if( avg_anyAngle_Range_Chunks_Long_Plot_autoTick > 10 )
                avg_anyAngle_Range_Chunks_Long_Plot_autoTick = 15; %sets the tick setting to 15 arcdegrees per tick
            elseif( avg_anyAngle_Range_Chunks_Long_Plot_autoTick > 5 )
                avg_anyAngle_Range_Chunks_Long_Plot_autoTick = 10; %sets the tick setting to 10 arcdegrees per tick
            elseif( avg_anyAngle_Range_Chunks_Long_Plot_autoTick > 2 )
                avg_anyAngle_Range_Chunks_Long_Plot_autoTick = 5; %sets the tick setting to 5 arcdegrees per tick
            elseif( avg_anyAngle_Range_Chunks_Long_Plot_autoTick > 1 )
                avg_anyAngle_Range_Chunks_Long_Plot_autoTick = 2; %sets the tick setting to 5 arcdegrees per tick
            elseif( avg_anyAngle_Range_Chunks_Long_Plot_autoTick >= 0.6 ) %0.6 because 15/25 = 0.6, so there will be enough 1 arcdeg ticks
                avg_anyAngle_Range_Chunks_Long_Plot_autoTick = 1; %sets the tick setting to 1 arcdegree per tick
            else
                avg_anyAngle_Range_Chunks_Long_Plot_autoTick = (max(plotLongRange) - min(plotLongRange))/15; %just goes for it if it's a super tiny range
            end
            Yaxisvar = round(floor(min(avg_anyAngle_Range_Chunks_Long_Plot)):avg_anyAngle_Range_Chunks_Long_Plot_autoTick:ceil(max(avg_anyAngle_Range_Chunks_Long_Plot)),2); %creates y ticks automagically
            set(gca,'ytick',Yaxisvar);

            hold on;
            %Now drawing line of interest
            if( strcmp(avg_anyAngle_Range_Chunks_Long_Plot_Name,'Longitude') == 1 ) %if true, longitude
                line([ min((timeUnique_limISR_MISA -floor(mean(timeUnique))).*24) , max((timeUnique_limISR_MISA -floor(mean(timeUnique))).*24) ],[longMillstoneMISA,longMillstoneMISA],'Color','k','LineWidth',.25);  %plots a point with a black line
            else %else latitude
                line([ min((timeUnique_limISR_MISA -floor(mean(timeUnique))).*24) , max((timeUnique_limISR_MISA -floor(mean(timeUnique))).*24) ],[latMillstoneMISA,latMillstoneMISA],'Color','k','LineWidth',.25);  %plots a point with a black line
            end

            xlabel(['Time in UT - 0 Hr on Day ',num2str(dateRange_zeroHr(2)),', ',num2str(dateRange_zeroHr(1)),' (hr)']); %old: ,'fontweight','bold','FontSize',12 for all
            ylabel([avg_anyAngle_Range_Chunks_Long_Plot_Name,' (arcdeg)']); %sped up significantly when merged above
            string_Title = ['TEC Averaged on Angle of ',num2str(round(avg_anyAngle*180/pi,2)),' deg and Width of ',num2str(round(avg_anyAngle_Width,2)),' arcdeg, Line Shows ',avg_anyAngle_Range_Chunks_Long_Plot_Name,' of Millstone Hill MISA @ ',num2str(pointAltitude),' km']; %create title string
            if( FLG_geomagneticCoords == 1)
                string_Title = [string_Title,' [Mag Coord Aligned]']; %tack on mag coordinate warning
            end
            title(string_Title);
            set(gca,'fontweight',font_Weight, 'FontSize',font_Size);
            hold off;

            %PLOT WITH MISA
            %prep figure #
            % figure('units','normalized','outerposition',[0 0 1 1]); %maximizes figure window
            figure(figNum);
            figNum = figNum+1;
            %Prep plot space
            ax1 = subplot(2,1,1); %upper bit
            h1 = pcolor( (MISA_time - dateRange_zeroHr(2)).*24 ,MISA_height,MISA_SNR_bp');
            set(h1, 'EdgeColor', 'none'); %remove grid lines (which make it look hilarious)
            shading interp;
            lighting phong;
            colorbar;
            caxis([-0.1 0.1]);
            colormap(ax1,'gray');
            ylim([90 700]);
            % ylim([90 500]);
            % xlim([0 70]);
    %         xlabel('Time in UT - 0 Hr at 0 UT May 7th');
            ylabel('Height (km)');
            title(['May 6-7-8 2013, MISA SNR 2 Hr Highpass Filtered, Latitude = ',num2str(round(latMillstoneMISA,2)),' arcdeg and Longitude = ',num2str(round(longMillstoneMISA,2)),' arcdeg']);
            % set(gca, 'fontweight','bold', 'FontSize',12);
            % Xaxisvar = [-12, -8, -4, 0, 4, 8, 12, 16, 20, 24 28, 32, 36, 40, 44, 48]; %shows off important magnitude points
            %         Xaxisvar = [-12,-8,-4,0, 4, 8, 12, 16, 20, 24 28, 32, 36, 40, 44, 48, 52, 56, 60, 64, 68, 72]; %shows off important magnitude points
            Xaxisvar = [ (round((min(MISA_time)-floor(mean(MISA_time)))*24)-mod(round((min(MISA_time)-floor(mean(MISA_time)))*24),2)):2:...
                (round((max(MISA_time)-floor(mean(MISA_time)))*24)-mod(round((max(MISA_time)-floor(mean(MISA_time)))*24),2)) ];
            set(gca,'TickDir','out'); %put tick marks on outside of grid to help connect to time elements
            set(gca,'box','off'); %keeps tick marks from being on top and right side too
            set(gca,'xtick',Xaxisvar);
            set(gca,'fontweight',font_Weight, 'FontSize',font_Size);
            % 
            % %add on a vertical line that shows where we are at in time
            % 
            ax2 = subplot(2,1,2); %lower bit
            h = pcolor( (timeUnique_limISR_MISA - dateRange_zeroHr(2)).*24 , avg_anyAngle_Range_Chunks_Long_Plot , sTECChunked_anyAngleAvg');% pseudocolor plot "stretched" to the grid
            set(h, 'EdgeColor', 'none'); %remove grid lines (which make it look hilarious)
            %             shading interp; %stuff for interpolation
            %             lighting phong; %stuff for interpolation
            caxis([-gif_sTEC_cap, gif_sTEC_cap]);
            %         colormap(flipud(colormap('hot'))) %heat map where 0 (low) is white
            colormap(ax2,'jet');
            h2 = colorbar; %shows the color bar
            ylabel(h2, 'delta-sTEC (TECU)','fontweight',font_Weight, 'FontSize',font_Size); %labels color bar
            set(gca,'TickDir','out'); %put tick marks on outside of grid to help connect to time elements
            set(gca,'box','off'); %keeps tick marks from being on top and right side too
            set(gca,'xtick',Xaxisvar);

            avg_anyAngle_Range_Chunks_Long_Plot_autoTick = (ceil(max(avg_anyAngle_Range_Chunks_Long_Plot)) - floor(min(avg_anyAngle_Range_Chunks_Long_Plot)))/13; %tries to split the latitude range into 13 parts (based off of 180/15+1)
            if( avg_anyAngle_Range_Chunks_Long_Plot_autoTick > 10 )
                avg_anyAngle_Range_Chunks_Long_Plot_autoTick = 15; %sets the tick setting to 15 arcdegrees per tick
            elseif( avg_anyAngle_Range_Chunks_Long_Plot_autoTick > 5 )
                avg_anyAngle_Range_Chunks_Long_Plot_autoTick = 10; %sets the tick setting to 10 arcdegrees per tick
            elseif( avg_anyAngle_Range_Chunks_Long_Plot_autoTick > 2 )
                avg_anyAngle_Range_Chunks_Long_Plot_autoTick = 5; %sets the tick setting to 5 arcdegrees per tick
            elseif( avg_anyAngle_Range_Chunks_Long_Plot_autoTick > 1 )
                avg_anyAngle_Range_Chunks_Long_Plot_autoTick = 2; %sets the tick setting to 5 arcdegrees per tick
            elseif( avg_anyAngle_Range_Chunks_Long_Plot_autoTick >= 0.6 ) %0.6 because 15/25 = 0.6, so there will be enough 1 arcdeg ticks
                avg_anyAngle_Range_Chunks_Long_Plot_autoTick = 1; %sets the tick setting to 1 arcdegree per tick
            else
                avg_anyAngle_Range_Chunks_Long_Plot_autoTick = (max(plotLongRange) - min(plotLongRange))/15; %just goes for it if it's a super tiny range
            end
            Yaxisvar = round(floor(min(avg_anyAngle_Range_Chunks_Long_Plot)):avg_anyAngle_Range_Chunks_Long_Plot_autoTick:ceil(max(avg_anyAngle_Range_Chunks_Long_Plot)),2); %creates y ticks automagically
            set(gca,'ytick',Yaxisvar);

            hold on;
            %Now drawing line of interest
            if( strcmp(avg_anyAngle_Range_Chunks_Long_Plot_Name,'Longitude') == 1 ) %if true, longitude
                line([ min((timeUnique_limISR_MISA -floor(mean(timeUnique))).*24) , max((timeUnique_limISR_MISA -floor(mean(timeUnique))).*24) ],[longMillstoneMISA,longMillstoneMISA],'Color','k','LineWidth',.25);  %plots a point with a black line
            else %else latitude
                line([ min((timeUnique_limISR_MISA -floor(mean(timeUnique))).*24) , max((timeUnique_limISR_MISA -floor(mean(timeUnique))).*24) ],[latMillstoneMISA,latMillstoneMISA],'Color','k','LineWidth',.25);  %plots a point with a black line
            end

            xlabel(['Time in UT - 0 Hr on Day ',num2str(dateRange_zeroHr(2)),', ',num2str(dateRange_zeroHr(1)),' (hr)']); %old: ,'fontweight','bold','FontSize',12 for all
            ylabel([avg_anyAngle_Range_Chunks_Long_Plot_Name,' (arcdeg)']); %sped up significantly when merged above
            string_Title = ['TEC Averaged on Angle of ',num2str(round(avg_anyAngle*180/pi,2)),' deg and Width of ',num2str(round(avg_anyAngle_Width,2)),' arcdeg, Line Shows ',avg_anyAngle_Range_Chunks_Long_Plot_Name,' of Millstone Hill MISA @ ',num2str(pointAltitude),' km']; %create title string
            if( FLG_geomagneticCoords == 1)
                string_Title = [string_Title,' [Mag Coord Aligned]']; %tack on mag coordinate warning
            end
            title(string_Title);
            set(gca,'fontweight',font_Weight, 'FontSize',font_Size);
            hold off;

        end
    end %END OF FULL TIME OR Zenith/MISA PLOTTING

    %JEDI PLOT
    if( FLG_JEDI_DATA == 1)
        if( FLG_AVG_ANYANGLE_ZENITHORMISA == 0 ) %Zenith option
            %prep figure #
            % figure('units','normalized','outerposition',[0 0 1 1]); %maximizes figure window
            figure(figNum);
            figNum = figNum+1;
            %Prep plot space
            ax1 = subplot(2,1,1); %upper bit
            plot(JEDI_time,JEDI_jouleHeating_bp);
            % ylim([90 500]);
            % xlim([0 70]);
            xlabel('Time in UT - 0 Hr at 0 UT May 7th');
            ylabel('Joule Heating (?))');
            title(['May 6-7-8 2013, JEDI Joule Heating Integrated Auroral Zone, High-Passed at ',num2str(1/JEDI_highpass_Freq),' Hr Cutoff Period']);
            % set(gca, 'fontweight','bold', 'FontSize',12);
            % Xaxisvar = [-12, -8, -4, 0, 4, 8, 12, 16, 20, 24 28, 32, 36, 40, 44, 48]; %shows off important magnitude points
            %         Xaxisvar = [-12,-8,-4,0, 4, 8, 12, 16, 20, 24 28, 32, 36, 40, 44, 48, 52, 56, 60, 64, 68, 72]; %shows off important magnitude points
            Xaxisvar = [ (round(min(JEDI_time))-mod(round((min(JEDI_time))),2)):2:...
                (round(max(JEDI_time))-mod(round((max(JEDI_time))),2)) ];
            set(gca,'xtick',Xaxisvar);
            xlim([min(JEDI_time) max(JEDI_time)]);
            % 
            % %add on a vertical line that shows where we are at in time
            % 
            ax2 = subplot(2,1,2); %lower bit
            h = pcolor( (timeUnique_limISR_Zenith - dateRange_zeroHr(2)).*24 , avg_anyAngle_Range_Chunks_Long_Plot , sTECChunked_anyAngleAvg');% pseudocolor plot "stretched" to the grid
            set(h, 'EdgeColor', 'none'); %remove grid lines (which make it look hilarious)
            %             shading interp; %stuff for interpolation
            %             lighting phong; %stuff for interpolation
            caxis([-gif_sTEC_cap, gif_sTEC_cap]);
            %         colormap(flipud(colormap('hot'))) %heat map where 0 (low) is white
            colormap(ax2,'jet');
            h2 = colorbar; %shows the color bar
            ylabel(h2, 'delta-sTEC (TECU)','fontweight',font_Weight, 'FontSize',font_Size); %labels color bar
            set(gca,'TickDir','out'); %put tick marks on outside of grid to help connect to time elements
            set(gca,'box','off'); %keeps tick marks from being on top and right side too
            set(gca,'xtick',Xaxisvar);

            avg_anyAngle_Range_Chunks_Long_Plot_autoTick = (ceil(max(avg_anyAngle_Range_Chunks_Long_Plot)) - floor(min(avg_anyAngle_Range_Chunks_Long_Plot)))/13; %tries to split the latitude range into 13 parts (based off of 180/15+1)
            if( avg_anyAngle_Range_Chunks_Long_Plot_autoTick > 10 )
                avg_anyAngle_Range_Chunks_Long_Plot_autoTick = 15; %sets the tick setting to 15 arcdegrees per tick
            elseif( avg_anyAngle_Range_Chunks_Long_Plot_autoTick > 5 )
                avg_anyAngle_Range_Chunks_Long_Plot_autoTick = 10; %sets the tick setting to 10 arcdegrees per tick
            elseif( avg_anyAngle_Range_Chunks_Long_Plot_autoTick > 2 )
                avg_anyAngle_Range_Chunks_Long_Plot_autoTick = 5; %sets the tick setting to 5 arcdegrees per tick
            elseif( avg_anyAngle_Range_Chunks_Long_Plot_autoTick > 1 )
                avg_anyAngle_Range_Chunks_Long_Plot_autoTick = 2; %sets the tick setting to 5 arcdegrees per tick
            elseif( avg_anyAngle_Range_Chunks_Long_Plot_autoTick >= 0.6 ) %0.6 because 15/25 = 0.6, so there will be enough 1 arcdeg ticks
                avg_anyAngle_Range_Chunks_Long_Plot_autoTick = 1; %sets the tick setting to 1 arcdegree per tick
            else
                avg_anyAngle_Range_Chunks_Long_Plot_autoTick = (max(plotLongRange) - min(plotLongRange))/15; %just goes for it if it's a super tiny range
            end
            Yaxisvar = round(floor(min(avg_anyAngle_Range_Chunks_Long_Plot)):avg_anyAngle_Range_Chunks_Long_Plot_autoTick:ceil(max(avg_anyAngle_Range_Chunks_Long_Plot)),2); %creates y ticks automagically
            set(gca,'ytick',Yaxisvar);

            hold on;
            %Now drawing line of interest
            if( strcmp(avg_anyAngle_Range_Chunks_Long_Plot_Name,'Longitude') == 1 ) %if true, longitude
                line([ min((timeUnique_limISR_MISA -floor(mean(timeUnique))).*24) , max((timeUnique_limISR_MISA -floor(mean(timeUnique))).*24) ],[longMillstoneMISA,longMillstoneMISA],'Color','k','LineWidth',.25);  %plots a point with a black line
            else %else latitude
                line([ min((timeUnique_limISR_MISA -floor(mean(timeUnique))).*24) , max((timeUnique_limISR_MISA -floor(mean(timeUnique))).*24) ],[latMillstoneMISA,latMillstoneMISA],'Color','k','LineWidth',.25);  %plots a point with a black line
            end

            xlabel(['Time in UT - 0 Hr on Day ',num2str(dateRange_zeroHr(2)),', ',num2str(dateRange_zeroHr(1)),' (hr)']); %old: ,'fontweight','bold','FontSize',12 for all
            ylabel([avg_anyAngle_Range_Chunks_Long_Plot_Name,' (deg)']); %sped up significantly when merged above
            string_Title = ['TEC Averaged on Angle of ',num2str(round(avg_anyAngle*180/pi,2)),' deg and Width of ',num2str(round(avg_anyAngle_Width,2)),' arcdeg, Line Shows ',avg_anyAngle_Range_Chunks_Long_Plot_Name,' of Millstone Hill MISA @ ',num2str(pointAltitude),' km']; %create title string
            if( FLG_geomagneticCoords == 1)
                string_Title = [string_Title,' [Mag Coord Aligned]']; %tack on mag coordinate warning
            end
            title(string_Title);
            set(gca,'fontweight',font_Weight, 'FontSize',font_Size);
            hold off;
            
        else %1 for MISA option
            
            %prep figure #
            % figure('units','normalized','outerposition',[0 0 1 1]); %maximizes figure window
            figure(figNum);
            figNum = figNum+1;
            %Prep plot space
            ax1 = subplot(2,1,1); %upper bit
            plot(JEDI_time,JEDI_jouleHeating_bp);
            % ylim([90 500]);
            % xlim([0 70]);
            xlabel('Time in UT - 0 Hr at 0 UT May 7th');
            ylabel('Joule Heating (?))');
            title(['May 6-7-8 2013, JEDI Joule Heating Integrated Auroral Zone, High-Passed at ',num2str(1/JEDI_highpass_Freq),' Hr Cutoff Period']);
            % set(gca, 'fontweight','bold', 'FontSize',12);
            % Xaxisvar = [-12, -8, -4, 0, 4, 8, 12, 16, 20, 24 28, 32, 36, 40, 44, 48]; %shows off important magnitude points
            %         Xaxisvar = [-12,-8,-4,0, 4, 8, 12, 16, 20, 24 28, 32, 36, 40, 44, 48, 52, 56, 60, 64, 68, 72]; %shows off important magnitude points
            Xaxisvar = [ (round(min(JEDI_time))-mod(round((min(JEDI_time))),2)):2:...
                (round(max(JEDI_time))-mod(round((max(JEDI_time))),2)) ];
            set(gca,'xtick',Xaxisvar);
            xlim([min(JEDI_time) max(JEDI_time)]);
            % 
            % %add on a vertical line that shows where we are at in time
            % 
            ax2 = subplot(2,1,2); %lower bit
            h = pcolor( (timeUnique_limISR_MISA - dateRange_zeroHr(2)).*24 , avg_anyAngle_Range_Chunks_Long_Plot , sTECChunked_anyAngleAvg');% pseudocolor plot "stretched" to the grid
            set(h, 'EdgeColor', 'none'); %remove grid lines (which make it look hilarious)
            %             shading interp; %stuff for interpolation
            %             lighting phong; %stuff for interpolation
            caxis([-gif_sTEC_cap, gif_sTEC_cap]);
            %         colormap(flipud(colormap('hot'))) %heat map where 0 (low) is white
            colormap(ax2,'jet');
            h2 = colorbar; %shows the color bar
            ylabel(h2, 'delta-sTEC (TECU)','fontweight',font_Weight, 'FontSize',font_Size); %labels color bar
            set(gca,'TickDir','out'); %put tick marks on outside of grid to help connect to time elements
            set(gca,'box','off'); %keeps tick marks from being on top and right side too
            set(gca,'xtick',Xaxisvar);

            avg_anyAngle_Range_Chunks_Long_Plot_autoTick = (ceil(max(avg_anyAngle_Range_Chunks_Long_Plot)) - floor(min(avg_anyAngle_Range_Chunks_Long_Plot)))/13; %tries to split the latitude range into 13 parts (based off of 180/15+1)
            if( avg_anyAngle_Range_Chunks_Long_Plot_autoTick > 10 )
                avg_anyAngle_Range_Chunks_Long_Plot_autoTick = 15; %sets the tick setting to 15 arcdegrees per tick
            elseif( avg_anyAngle_Range_Chunks_Long_Plot_autoTick > 5 )
                avg_anyAngle_Range_Chunks_Long_Plot_autoTick = 10; %sets the tick setting to 10 arcdegrees per tick
            elseif( avg_anyAngle_Range_Chunks_Long_Plot_autoTick > 2 )
                avg_anyAngle_Range_Chunks_Long_Plot_autoTick = 5; %sets the tick setting to 5 arcdegrees per tick
            elseif( avg_anyAngle_Range_Chunks_Long_Plot_autoTick > 1 )
                avg_anyAngle_Range_Chunks_Long_Plot_autoTick = 2; %sets the tick setting to 5 arcdegrees per tick
            elseif( avg_anyAngle_Range_Chunks_Long_Plot_autoTick >= 0.6 ) %0.6 because 15/25 = 0.6, so there will be enough 1 arcdeg ticks
                avg_anyAngle_Range_Chunks_Long_Plot_autoTick = 1; %sets the tick setting to 1 arcdegree per tick
            else
                avg_anyAngle_Range_Chunks_Long_Plot_autoTick = (max(plotLongRange) - min(plotLongRange))/15; %just goes for it if it's a super tiny range
            end
            Yaxisvar = round(floor(min(avg_anyAngle_Range_Chunks_Long_Plot)):avg_anyAngle_Range_Chunks_Long_Plot_autoTick:ceil(max(avg_anyAngle_Range_Chunks_Long_Plot)),2); %creates y ticks automagically
            set(gca,'ytick',Yaxisvar);

            hold on;
            %Now drawing line of interest
            if( strcmp(avg_anyAngle_Range_Chunks_Long_Plot_Name,'Longitude') == 1 ) %if true, longitude
                line([ min((timeUnique_limISR_MISA -floor(mean(timeUnique))).*24) , max((timeUnique_limISR_MISA -floor(mean(timeUnique))).*24) ],[longMillstoneMISA,longMillstoneMISA],'Color','k','LineWidth',.25);  %plots a point with a black line
            else %else latitude
                line([ min((timeUnique_limISR_MISA -floor(mean(timeUnique))).*24) , max((timeUnique_limISR_MISA -floor(mean(timeUnique))).*24) ],[latMillstoneMISA,latMillstoneMISA],'Color','k','LineWidth',.25);  %plots a point with a black line
            end

            xlabel(['Time in UT - 0 Hr on Day ',num2str(dateRange_zeroHr(2)),', ',num2str(dateRange_zeroHr(1)),' (hr)']); %old: ,'fontweight','bold','FontSize',12 for all
            ylabel([avg_anyAngle_Range_Chunks_Long_Plot_Name,' (deg)']); %sped up significantly when merged above
            string_Title = ['TEC Averaged on Angle of ',num2str(round(avg_anyAngle*180/pi,2)),' deg and Width of ',num2str(round(avg_anyAngle_Width,2)),' arcdeg, Line Shows ',avg_anyAngle_Range_Chunks_Long_Plot_Name,' of Millstone Hill MISA @ ',num2str(pointAltitude),' km']; %create title string
            if( FLG_geomagneticCoords == 1)
                string_Title = [string_Title,' [Mag Coord Aligned]']; %tack on mag coordinate warning
            end
            title(string_Title);
            set(gca,'fontweight',font_Weight, 'FontSize',font_Size);
            hold off;
        end
    end
    
    %Kp Index Plot & OMNI Bz GSE
    if( FLG_solarGeoStuff == 1 )
        %Kp Time Limited PLOT HERE
        Kp_time_limit_lo = find( Kp_time > min((timeUnique_limISR_Zenith -floor(mean(timeUnique))).*24) ,1); %get the nearest min time to the ISR stuff that is within the range
%         Kp_time_limit_hi = find(min(abs(Kp_time - max((timeUnique_limISR_Zenith -floor(mean(timeUnique))).*24))) == abs(Kp_time - max((timeUnique_limISR_Zenith -floor(mean(timeUnique))).*24))); %get the nearest max time to the ISR stuff
        Kp_time_limit_hi = find( Kp_time < max((timeUnique_limISR_Zenith -floor(mean(timeUnique))).*24) ,1,'last'); %get the nearest max time to the ISR stuff that is within the range
        %prep figure #
        % figure('units','normalized','outerposition',[0 0 1 1]); %maximizes figure window (caused graphical issues on Windows 10 + AMD GPU or Intel GPU)
        Kp_time_limISR_Zenith = Kp_time(Kp_time_limit_lo:Kp_time_limit_hi); %cut out the limited time
        Kp_output_limISR_Zenith = reshape(Kp_output',[size(Kp_output,1)*size(Kp_output,2),1]);
        Kp_output_limISR_Zenith = Kp_output_limISR_Zenith(Kp_time_limit_lo:Kp_time_limit_hi); %cut out the data that matches the time limit
        figure(figNum);
        figNum = figNum+1;
        ax1 = subplot(2,1,1); %upper bit
        plot( Kp_time_limISR_Zenith , Kp_output_limISR_Zenith )
        xlabel('Time (hr)'); %old: ,'fontweight','bold','FontSize',12 for all
        ylabel('Kp Index, 3 Hr Resolution'); %sped up significantly when merged above
        title(['Kp Index for ',num2str(Kp_dateRange_Full_Month(1,2)),'/',num2str(Kp_dateRange_Full_Month(1,3)),...
            '/',num2str(Kp_dateRange_Full_Month(1,1)),' to ',num2str(Kp_dateRange_Full_Month(end,2)),...
            '/',num2str(Kp_dateRange_Full_Month(end,3)),'/',num2str(Kp_dateRange_Full_Month(end,1)),...
            ' (M/D/Y) with 0 Hr on May 7th, 2013']);
        Xaxisvar = [ (round((min(Zenith_time)-floor(mean(Zenith_time)))*24)-mod(round((min(Zenith_time)-floor(mean(Zenith_time)))*24),2)):2:...
            (round((max(Zenith_time)-floor(mean(Zenith_time)))*24)-mod(round((max(Zenith_time)-floor(mean(Zenith_time)))*24),2)) ];
        set(gca,'xtick',Xaxisvar);
        set(gca,'fontweight',font_Weight, 'FontSize',font_Size);
        % 
        % %add on a vertical line that shows where we are at in time
        % 
        ax2 = subplot(2,1,2); %lower bit
        h = pcolor( (timeUnique_limISR_Zenith - dateRange_zeroHr(2)).*24 , avg_anyAngle_Range_Chunks_Long_Plot , sTECChunked_anyAngleAvg');% pseudocolor plot "stretched" to the grid
        set(h, 'EdgeColor', 'none'); %remove grid lines (which make it look hilarious)
        %             shading interp; %stuff for interpolation
        %             lighting phong; %stuff for interpolation
        caxis([-gif_sTEC_cap, gif_sTEC_cap]);
        %         colormap(flipud(colormap('hot'))) %heat map where 0 (low) is white
        colormap(ax2,'jet');
        h2 = colorbar; %shows the color bar
        ylabel(h2, 'delta-sTEC (TECU)','fontweight',font_Weight, 'FontSize',font_Size); %labels color bar
        set(gca,'TickDir','out'); %put tick marks on outside of grid to help connect to time elements
        set(gca,'box','off'); %keeps tick marks from being on top and right side too
        set(gca,'xtick',Xaxisvar);

        avg_anyAngle_Range_Chunks_Long_Plot_autoTick = (ceil(max(avg_anyAngle_Range_Chunks_Long_Plot)) - floor(min(avg_anyAngle_Range_Chunks_Long_Plot)))/13; %tries to split the latitude range into 13 parts (based off of 180/15+1)
        if( avg_anyAngle_Range_Chunks_Long_Plot_autoTick > 25 )
            avg_anyAngle_Range_Chunks_Long_Plot_autoTick = 30; %sets the tick setting to 15 arcdegrees per tick
        elseif( avg_anyAngle_Range_Chunks_Long_Plot_autoTick > 10 )
            avg_anyAngle_Range_Chunks_Long_Plot_autoTick = 15; %sets the tick setting to 15 arcdegrees per tick
        elseif( avg_anyAngle_Range_Chunks_Long_Plot_autoTick > 5 )
            avg_anyAngle_Range_Chunks_Long_Plot_autoTick = 10; %sets the tick setting to 10 arcdegrees per tick
        elseif( avg_anyAngle_Range_Chunks_Long_Plot_autoTick > 2 )
            avg_anyAngle_Range_Chunks_Long_Plot_autoTick = 5; %sets the tick setting to 5 arcdegrees per tick
        elseif( avg_anyAngle_Range_Chunks_Long_Plot_autoTick > 1 )
            avg_anyAngle_Range_Chunks_Long_Plot_autoTick = 2; %sets the tick setting to 5 arcdegrees per tick
        elseif( avg_anyAngle_Range_Chunks_Long_Plot_autoTick >= 0.6 ) %0.6 because 15/25 = 0.6, so there will be enough 1 arcdeg ticks
            avg_anyAngle_Range_Chunks_Long_Plot_autoTick = 1; %sets the tick setting to 1 arcdegree per tick
        else
            avg_anyAngle_Range_Chunks_Long_Plot_autoTick = (max(plotLongRange) - min(plotLongRange))/13; %just goes for it if it's a super tiny range
        end
        Yaxisvar = round(floor(min(avg_anyAngle_Range_Chunks_Long_Plot)):avg_anyAngle_Range_Chunks_Long_Plot_autoTick:ceil(max(avg_anyAngle_Range_Chunks_Long_Plot)),2); %creates y ticks automagically
        set(gca,'ytick',Yaxisvar);

        hold on;
        %Now drawing line of interest
        if( strcmp(avg_anyAngle_Range_Chunks_Long_Plot_Name,'Longitude') == 1 ) %if true, longitude
            line([ min((timeUnique_limISR_Zenith -floor(mean(timeUnique))).*24) , max((timeUnique_limISR_Zenith -floor(mean(timeUnique))).*24) ],[longMillstone,longMillstone],'Color','k','LineWidth',.25);  %plots a point with a black line
        else %else latitude
            line([ min((timeUnique_limISR_Zenith -floor(mean(timeUnique))).*24) , max((timeUnique_limISR_Zenith -floor(mean(timeUnique))).*24) ],[latMillstone,latMillstone],'Color','k','LineWidth',.25);  %plots a point with a black line
        end

        xlabel(['Time in UT - 0 Hr on Day ',num2str(dateRange_zeroHr(2)),', ',num2str(dateRange_zeroHr(1)),' (hr)']); %old: ,'fontweight','bold','FontSize',12 for all
        ylabel([avg_anyAngle_Range_Chunks_Long_Plot_Name,' (arcdeg)']); %sped up significantly when merged above
        string_Title = ['TEC Averaged on Angle of ',num2str(round(avg_anyAngle*180/pi,2)),' deg and Width of ',num2str(round(avg_anyAngle_Width,2)),' arcdeg, Line Shows ',avg_anyAngle_Range_Chunks_Long_Plot_Name,' of Millstone Hill Zenith Beam']; %create title string
        if( FLG_geomagneticCoords == 1)
            string_Title = [string_Title,' [Mag Coord Aligned]']; %tack on mag coordinate warning
        end
        title(string_Title);
        set(gca,'fontweight',font_Weight, 'FontSize',font_Size);
        hold off;

        
        
        
        
        %OMNI Time Limited PLOT HERE
        OMNI_time_hr_limit_lo = find( OMNI_time_hr > min((timeUnique_limISR_Zenith -floor(mean(timeUnique))).*24) ,1); %get the nearest min time to the ISR stuff that is within the range
%         Kp_time_limit_hi = find(min(abs(Kp_time - max((timeUnique_limISR_Zenith -floor(mean(timeUnique))).*24))) == abs(Kp_time - max((timeUnique_limISR_Zenith -floor(mean(timeUnique))).*24))); %get the nearest max time to the ISR stuff
        OMNI_time_hr_limit_hi = find( OMNI_time_hr < max((timeUnique_limISR_Zenith -floor(mean(timeUnique))).*24) ,1,'last'); %get the nearest max time to the ISR stuff that is within the range
        %prep figure #
        % figure('units','normalized','outerposition',[0 0 1 1]); %maximizes figure window (caused graphical issues on Windows 10 + AMD GPU or Intel GPU)
        OMNI_time_hr_limISR_Zenith = OMNI_time_hr(OMNI_time_hr_limit_lo:OMNI_time_hr_limit_hi); %cut out the limited time
        OMNI_output_Bz_BP_limISR_Zenith = OMNI_output_Bz_BP(OMNI_time_hr_limit_lo:OMNI_time_hr_limit_hi); %cut out the data that matches the time limit
        figure(figNum);
        figNum = figNum+1;
        ax1 = subplot(2,1,1); %upper bit
        plot( OMNI_time_hr_limISR_Zenith , OMNI_output_Bz_BP_limISR_Zenith )
        xlabel('Time (hr)'); %old: ,'fontweight','bold','FontSize',12 for all
        ylabel('OMNI Bz GSE HP (nT)'); %sped up significantly when merged above
        title(['OMNI Bz GSE High-passed w/ Cutoff Freq ',num2str(round(1/bs)),' hr for ',num2str(OMNI_dateRange_Full_Month(1,2)),'/',num2str(OMNI_dateRange_Full_Month(1,3)),...
            '/',num2str(OMNI_dateRange_Full_Month(1,1)),' to ',num2str(OMNI_dateRange_Full_Month(end,2)),...
            '/',num2str(OMNI_dateRange_Full_Month(end,3)),'/',num2str(OMNI_dateRange_Full_Month(end,1)),...
            ' (M/D/Y) with 0 Hr on May 7th, 2013']);
        Xaxisvar = [ (round((min(Zenith_time)-floor(mean(Zenith_time)))*24)-mod(round((min(Zenith_time)-floor(mean(Zenith_time)))*24),2)):2:...
            (round((max(Zenith_time)-floor(mean(Zenith_time)))*24)-mod(round((max(Zenith_time)-floor(mean(Zenith_time)))*24),2)) ];
        set(gca,'xtick',Xaxisvar);
        set(gca,'fontweight',font_Weight, 'FontSize',font_Size);
        % 
        % %add on a vertical line that shows where we are at in time
        % 
        ax2 = subplot(2,1,2); %lower bit
        h = pcolor( (timeUnique_limISR_Zenith - dateRange_zeroHr(2)).*24 , avg_anyAngle_Range_Chunks_Long_Plot , sTECChunked_anyAngleAvg');% pseudocolor plot "stretched" to the grid
        set(h, 'EdgeColor', 'none'); %remove grid lines (which make it look hilarious)
        %             shading interp; %stuff for interpolation
        %             lighting phong; %stuff for interpolation
        caxis([-gif_sTEC_cap, gif_sTEC_cap]);
        %         colormap(flipud(colormap('hot'))) %heat map where 0 (low) is white
        colormap(ax2,'jet');
        h2 = colorbar; %shows the color bar
        ylabel(h2, 'delta-sTEC (TECU)','fontweight',font_Weight, 'FontSize',font_Size); %labels color bar
        set(gca,'TickDir','out'); %put tick marks on outside of grid to help connect to time elements
        set(gca,'box','off'); %keeps tick marks from being on top and right side too
        set(gca,'xtick',Xaxisvar);

        avg_anyAngle_Range_Chunks_Long_Plot_autoTick = (ceil(max(avg_anyAngle_Range_Chunks_Long_Plot)) - floor(min(avg_anyAngle_Range_Chunks_Long_Plot)))/13; %tries to split the latitude range into 13 parts (based off of 180/15+1)
        if( avg_anyAngle_Range_Chunks_Long_Plot_autoTick > 25 )
            avg_anyAngle_Range_Chunks_Long_Plot_autoTick = 30; %sets the tick setting to 15 arcdegrees per tick
        elseif( avg_anyAngle_Range_Chunks_Long_Plot_autoTick > 10 )
            avg_anyAngle_Range_Chunks_Long_Plot_autoTick = 15; %sets the tick setting to 15 arcdegrees per tick
        elseif( avg_anyAngle_Range_Chunks_Long_Plot_autoTick > 5 )
            avg_anyAngle_Range_Chunks_Long_Plot_autoTick = 10; %sets the tick setting to 10 arcdegrees per tick
        elseif( avg_anyAngle_Range_Chunks_Long_Plot_autoTick > 2 )
            avg_anyAngle_Range_Chunks_Long_Plot_autoTick = 5; %sets the tick setting to 5 arcdegrees per tick
        elseif( avg_anyAngle_Range_Chunks_Long_Plot_autoTick > 1 )
            avg_anyAngle_Range_Chunks_Long_Plot_autoTick = 2; %sets the tick setting to 5 arcdegrees per tick
        elseif( avg_anyAngle_Range_Chunks_Long_Plot_autoTick >= 0.6 ) %0.6 because 15/25 = 0.6, so there will be enough 1 arcdeg ticks
            avg_anyAngle_Range_Chunks_Long_Plot_autoTick = 1; %sets the tick setting to 1 arcdegree per tick
        else
            avg_anyAngle_Range_Chunks_Long_Plot_autoTick = (max(plotLongRange) - min(plotLongRange))/13; %just goes for it if it's a super tiny range
        end
        Yaxisvar = round(floor(min(avg_anyAngle_Range_Chunks_Long_Plot)):avg_anyAngle_Range_Chunks_Long_Plot_autoTick:ceil(max(avg_anyAngle_Range_Chunks_Long_Plot)),2); %creates y ticks automagically
        set(gca,'ytick',Yaxisvar);

        hold on;
        %Now drawing line of interest
        if( strcmp(avg_anyAngle_Range_Chunks_Long_Plot_Name,'Longitude') == 1 ) %if true, longitude
            line([ min((timeUnique_limISR_Zenith -floor(mean(timeUnique))).*24) , max((timeUnique_limISR_Zenith -floor(mean(timeUnique))).*24) ],[longMillstone,longMillstone],'Color','k','LineWidth',.25);  %plots a point with a black line
        else %else latitude
            line([ min((timeUnique_limISR_Zenith -floor(mean(timeUnique))).*24) , max((timeUnique_limISR_Zenith -floor(mean(timeUnique))).*24) ],[latMillstone,latMillstone],'Color','k','LineWidth',.25);  %plots a point with a black line
        end

        xlabel(['Time in UT - 0 Hr on Day ',num2str(dateRange_zeroHr(2)),', ',num2str(dateRange_zeroHr(1)),' (hr)']); %old: ,'fontweight','bold','FontSize',12 for all
        ylabel([avg_anyAngle_Range_Chunks_Long_Plot_Name,' (arcdeg)']); %sped up significantly when merged above
        string_Title = ['TEC Averaged on Angle of ',num2str(round(avg_anyAngle*180/pi,2)),' deg and Width of ',num2str(round(avg_anyAngle_Width,2)),' arcdeg, Line Shows ',avg_anyAngle_Range_Chunks_Long_Plot_Name,' of Millstone Hill Zenith Beam']; %create title string
        if( FLG_geomagneticCoords == 1)
            string_Title = [string_Title,' [Mag Coord Aligned]']; %tack on mag coordinate warning
        end
        title(string_Title);
        set(gca,'fontweight',font_Weight, 'FontSize',font_Size);
        hold off;

        
        % OMNI Bz GSE Scarglin'
        Ross_Lombscargle_optimized([ OMNI_time_hr_limISR_Zenith , OMNI_output_Bz_BP_limISR_Zenith ]); 
        fid = fopen('power_freq.txt', 'r');
            f_nm = textscan(fid, '%f %f %f %f %f', 'headerlines', 15);
            P_Zenith_SNR=f_nm{:,2};
            F_Zenith_SNR=f_nm{:,1};
            T_Zenith_SNR=f_nm{:,5}*60; % in minutes
            gs_Zenith_SNR=f_nm{:,4};
            max_Zenith_SNR=T_Zenith_SNR(P_Zenith_SNR==max(P_Zenith_SNR));
            K=length(P_Zenith_SNR);
            gf_Zenith_SNR=f_nm{:,3};
        fclose('all');

        %prep figure #
        % figure('units','normalized','outerposition',[0 0 1 1]); %maximizes figure window
        figure(figNum);
        figNum = figNum+1;

        plot(T_Zenith_SNR,P_Zenith_SNR,'b');
        xlim([0 plot_Freq_Lim]);
        hold on;
        % plot(T_SNR_08apr08_10_or,P_SNR_08apr08_10_or,'r');
        % hold on;
        plot(T_Zenith_SNR,gf_Zenith_SNR,'LineWidth',1.5,'color',[0.5 0.5 0.5]);
        xlim([0 plot_Freq_Lim]);  
        xlabel('Periods (min)');
        ylabel('Normalized Power');
        string_Title = ['Bz GSE High-pass Cutoff Freq of ',num2str(round(1/bs)),' Lomb-Scargle Periodogram -- 0 Hr on May 7th, 2013'];
        if( FLG_geomagneticCoords == 1)
            string_Title = [string_Title,' [Mag Coord Aligned]']; %tack on mag coordinate warning
        end
        title(string_Title);
        set(gca,'fontweight',font_Weight, 'FontSize',font_Size);
        set(gca,'xtick',Xaxisvar_SCARGLE);
    end
    

end %*********END OF FLG_AVG_ANYANGLE_PLOT*********


%% Average Any Angle, PLOT SELECT TIME UP
if( FLG_AVG_ANYANGLE_PLOT_TIMECUTOUT == 1 )
    
    %THIS BIT IS FOR FINAL PLOT OF STEC vs EXPANDED ZENITH vs EXPANDED MISA
    tmin = find( min( abs(min(time_cutout_range) - (timeUnique_limISR_MISA - floor(mean(timeUnique))).*24 ) ) == abs(min(time_cutout_range) - (timeUnique_limISR_MISA- floor(mean(timeUnique))).*24 ) ); %find min timeUnique that fits time limit
    tmax = find( min( abs(max(time_cutout_range) - (timeUnique_limISR_MISA -  floor(mean(timeUnique))).*24 ) ) == abs(max(time_cutout_range) - (timeUnique_limISR_MISA- floor(mean(timeUnique))).*24 ) ); %find max timeUnique that fits time limit
    timeUnique_limTIME_anyAngleAvg = timeUnique_limISR_MISA(tmin:tmax); %cut out time that matches time limits (Zenith or MISA have been expanded to this time cadence)
    sTECChunked_anyAngleAvg_limTIME = sTECChunked_anyAngleAvg(tmin:tmax,:); %cut out the time that matches the TEC time
    %THIS BIT IS FOR FINAL PLOT OF STEC vs EXPANDED ZENITH vs EXPANDED Zenith
    tmin = find( min( abs(min(time_cutout_range) - (timeUnique_limISR_Zenith - floor(mean(timeUnique))).*24 ) ) == abs(min(time_cutout_range) - (timeUnique_limISR_Zenith- floor(mean(timeUnique))).*24 ) ); %find min timeUnique that fits time limit
    tmax = find( min( abs(max(time_cutout_range) - (timeUnique_limISR_Zenith -  floor(mean(timeUnique))).*24 ) ) == abs(max(time_cutout_range) - (timeUnique_limISR_Zenith- floor(mean(timeUnique))).*24 ) ); %find max timeUnique that fits time limit
    timeUnique_limTIME_Zenith = timeUnique_limISR_Zenith(tmin:tmax); %cut out time that matches time limits (Zenith or MISA have been expanded to this time cadence)
    %THIS BIT IS FOR FINAL PLOT OF STEC vs EXPANDED ZENITH vs EXPANDED MISA
    tmin = find( min( abs(min(time_cutout_range) - (timeUnique_limISR_MISA - floor(mean(timeUnique))).*24 ) ) == abs(min(time_cutout_range) - (timeUnique_limISR_MISA- floor(mean(timeUnique))).*24 ) ); %find min timeUnique that fits time limit
    tmax = find( min( abs(max(time_cutout_range) - (timeUnique_limISR_MISA -  floor(mean(timeUnique))).*24 ) ) == abs(max(time_cutout_range) - (timeUnique_limISR_MISA- floor(mean(timeUnique))).*24 ) ); %find max timeUnique that fits time limit
    timeUnique_limTIME_MISA = timeUnique_limISR_MISA(tmin:tmax); %cut out time that matches time limits (Zenith or MISA have been expanded to this time cadence)
    
    %THIS BIT IS FOR ZENITH SNR PLOTS (w/ height)
    tmin = find( min( abs(min(time_cutout_range) - (Zenith_time - floor(mean(Zenith_time))).*24 ) ) == abs(min(time_cutout_range) - (Zenith_time- floor(mean(Zenith_time))).*24 ) ); %find min timeUnique that fits time limit
    tmax = find( min( abs(max(time_cutout_range) - (Zenith_time -  floor(mean(Zenith_time))).*24 ) ) == abs(max(time_cutout_range) - (Zenith_time- floor(mean(Zenith_time))).*24 ) ); %find max timeUnique that fits time limit
    Zenith_time_limTIME = Zenith_time(tmin:tmax); %cut out time that matches time limits (Zenith or MISA have been expanded to this time cadence)
    Zenith_SNR_bp_limTIME = Zenith_SNR_bp(tmin:tmax,:); %cut out time that matches ISR time

    %THIS BIT IS FOR MISA SNR PLOTS (w/ height)
    tmin = find( min( abs(min(time_cutout_range) - (MISA_time - floor(mean(MISA_time))).*24 ) ) == abs(min(time_cutout_range) - (MISA_time- floor(mean(MISA_time))).*24 ) ); %find min timeUnique that fits time limit
    tmax = find( min( abs(max(time_cutout_range) - (MISA_time -  floor(mean(MISA_time))).*24 ) ) == abs(max(time_cutout_range) - (MISA_time- floor(mean(MISA_time))).*24 ) ); %find max timeUnique that fits time limit
    MISA_time_limTIME = MISA_time(tmin:tmax); %cut out time that matches time limits (Zenith or MISA have been expanded to this time cadence)
    MISA_SNR_bp_limTIME = MISA_SNR_bp(tmin:tmax,:); %cut out time that matches ISR time
    
    if( FLG_AVG_ANYANGLE_ZENITHORMISA == 0 ) %Zenith option
        %PLOT W/ ZOOM
        if( isnan(latMillstone) ~= 1) %only run if latMillstone is relevant
            %prep figure #
            % figure('units','normalized','outerposition',[0 0 1 1]); %maximizes figure window
            figure(figNum);
            figNum = figNum+1;
            %Prep plot space
            ax1 = subplot(2,1,1); %upper bit
            h1 = pcolor( (Zenith_time_limTIME-floor(mean(Zenith_time))).*24 ,Zenith_height,Zenith_SNR_bp_limTIME');
            set(h1, 'EdgeColor', 'none'); %remove grid lines (which make it look hilarious)
            shading interp;
            lighting phong;
            colorbar;
            caxis([-0.1 0.1]);
            colormap(ax1,'gray');
            ylim([90 700]);
            % ylim([90 500]);
            % xlim([0 70]);
    %         xlabel('Time in UT - 0 Hr at 0 UT May 7th');
            ylabel('Height (km)');
            title(['Time Cutout Range:',num2str(min(time_cutout_range)),':',num2str(max(time_cutout_range)),', Zenith SNR 2 Hr Highpass Filtered, Latitude = ',num2str(round(latMillstone,2)),' arcdeg and Longitude = ',num2str(round(longMillstone,2)),' arcdeg']);
            % set(gca, 'fontweight','bold', 'FontSize',12);
            % Xaxisvar = [-12, -8, -4, 0, 4, 8, 12, 16, 20, 24 28, 32, 36, 40, 44, 48]; %shows off important magnitude points
            %         Xaxisvar = [-12,-8,-4,0, 4, 8, 12, 16, 20, 24 28, 32, 36, 40, 44, 48, 52, 56, 60, 64, 68, 72]; %shows off important magnitude points
            Xaxisvar = [ (round((min(Zenith_time_limTIME)-floor(mean(Zenith_time)))*24)-mod(round((min(Zenith_time_limTIME)-floor(mean(Zenith_time)))*24),2)):2:...
                (round((max(Zenith_time_limTIME)-floor(mean(Zenith_time)))*24)-mod(round((max(Zenith_time_limTIME)-floor(mean(Zenith_time)))*24),2)) ];
            set(gca,'TickDir','out'); %put tick marks on outside of grid to help connect to time elements
            set(gca,'box','off'); %keeps tick marks from being on top and right side too
            set(gca,'xtick',Xaxisvar);
            set(gca,'fontweight',font_Weight, 'FontSize',font_Size);
            % 
            % %add on a vertical line that shows where we are at in time
            % 
            ax2 = subplot(2,1,2); %lower bit
            if( strcmp(avg_anyAngle_Range_Chunks_Long_Plot_Name,'Longitude') == 1 ) %if Longitude, use Zenith or MISA Longitude value
                avg_anyAngle_Zoom_Range = (longMillstone + avg_anyAngle_Zoom) >= avg_anyAngle_Range_Chunks_Long_Plot & (longMillstone - avg_anyAngle_Zoom) <= avg_anyAngle_Range_Chunks_Long_Plot; %gets the zoom range 
            else
                 avg_anyAngle_Zoom_Range = (latMillstone + avg_anyAngle_Zoom) >= avg_anyAngle_Range_Chunks_Long_Plot & (latMillstone - avg_anyAngle_Zoom) <= avg_anyAngle_Range_Chunks_Long_Plot; %gets the zoom range 
            end
            h = pcolor( (timeUnique_limTIME_anyAngleAvg -floor(mean(timeUnique))).*24 , avg_anyAngle_Range_Chunks_Long_Plot( avg_anyAngle_Zoom_Range ) , sTECChunked_anyAngleAvg_limTIME(:,avg_anyAngle_Zoom_Range)');% pseudocolor plot "stretched" to the grid
            set(h, 'EdgeColor', 'none'); %remove grid lines (which make it look hilarious)
            %             shading interp; %stuff for interpolation
            %             lighting phong; %stuff for interpolation
            caxis([-gif_sTEC_cap, gif_sTEC_cap]);
            %         colormap(flipud(colormap('hot'))) %heat map where 0 (low) is white
            colormap(ax2,'jet');
            h2 = colorbar; %shows the color bar
            ylabel(h2, 'delta-sTEC (TECU)','fontweight',font_Weight, 'FontSize',font_Size); %labels color bar
            set(gca,'xtick',Xaxisvar);

            avg_anyAngle_Range_Chunks_Long_Plot_autoTick = (ceil(max(avg_anyAngle_Range_Chunks_Long_Plot(avg_anyAngle_Zoom_Range))) - floor(min(avg_anyAngle_Range_Chunks_Long_Plot(avg_anyAngle_Zoom_Range))))/13; %tries to split the latitude range into 13 parts (based off of 180/15+1)
            if( avg_anyAngle_Range_Chunks_Long_Plot_autoTick > 10 )
                avg_anyAngle_Range_Chunks_Long_Plot_autoTick = 15; %sets the tick setting to 15 arcdegrees per tick
            elseif( avg_anyAngle_Range_Chunks_Long_Plot_autoTick > 5 )
                avg_anyAngle_Range_Chunks_Long_Plot_autoTick = 10; %sets the tick setting to 10 arcdegrees per tick
            elseif( avg_anyAngle_Range_Chunks_Long_Plot_autoTick > 2 )
                avg_anyAngle_Range_Chunks_Long_Plot_autoTick = 5; %sets the tick setting to 5 arcdegrees per tick
            elseif( avg_anyAngle_Range_Chunks_Long_Plot_autoTick > 1 )
                avg_anyAngle_Range_Chunks_Long_Plot_autoTick = 2; %sets the tick setting to 5 arcdegrees per tick
            elseif( avg_anyAngle_Range_Chunks_Long_Plot_autoTick >= 0.6 ) %0.6 because 15/25 = 0.6, so there will be enough 1 arcdeg ticks
                avg_anyAngle_Range_Chunks_Long_Plot_autoTick = 1; %sets the tick setting to 1 arcdegree per tick
            else
                avg_anyAngle_Range_Chunks_Long_Plot_autoTick = (max(avg_anyAngle_Range_Chunks_Long_Plot(avg_anyAngle_Zoom_Range)) - min(avg_anyAngle_Range_Chunks_Long_Plot(avg_anyAngle_Zoom_Range)))/13; %just goes for it if it's a super tiny range
            end
            Yaxisvar = round(floor(min(avg_anyAngle_Range_Chunks_Long_Plot(avg_anyAngle_Zoom_Range))):avg_anyAngle_Range_Chunks_Long_Plot_autoTick:ceil(max(avg_anyAngle_Range_Chunks_Long_Plot(avg_anyAngle_Zoom_Range))),2); %creates y ticks automagically
            set(gca,'ytick',Yaxisvar);

            hold on;
            %Now drawing line of interest
            if( strcmp(avg_anyAngle_Range_Chunks_Long_Plot_Name,'Longitude') == 1 ) %if true, longitude
                line([ min((timeUnique_limISR_Zenith -floor(mean(timeUnique))).*24) , max((timeUnique_limISR_Zenith -floor(mean(timeUnique))).*24) ],[longMillstone,longMillstone],'Color','k','LineWidth',.25);  %plots a point with a black line
            else %else latitude
                line([ min((timeUnique_limISR_Zenith -floor(mean(timeUnique))).*24) , max((timeUnique_limISR_Zenith -floor(mean(timeUnique))).*24) ],[latMillstone,latMillstone],'Color','k','LineWidth',.25);  %plots a point with a black line
            end

            xlabel(['Time in UT - 0 Hr on Day ',num2str(dateRange_zeroHr(2)),', ',num2str(dateRange_zeroHr(1)),' (hr)']); %old: ,'fontweight','bold','FontSize',12 for all
            ylabel([avg_anyAngle_Range_Chunks_Long_Plot_Name,' (arcdeg)']); %sped up significantly when merged above
            string_Title = ['TEC Averaged on Angle of ',num2str(round(avg_anyAngle*180/pi,2)),' deg and Width of ',num2str(round(avg_anyAngle_Width,2)),' arcdeg, Line Shows ',avg_anyAngle_Range_Chunks_Long_Plot_Name,' of Millstone Hill Zenith Beam with Zoom of +/- ',num2str(avg_anyAngle_Zoom),'deg around']; %create title string
            if( FLG_geomagneticCoords == 1)
                string_Title = [string_Title,' [Mag Coord Aligned]']; %tack on mag coordinate warning
            end
            title(string_Title);
            set(gca,'fontweight',font_Weight, 'FontSize',font_Size);
            hold off;
        end
        
        %PLOT W/O ZOOM
        %prep figure #
        % figure('units','normalized','outerposition',[0 0 1 1]); %maximizes figure window
        figure(figNum);
        figNum = figNum+1;
        %Prep plot space
        ax1 = subplot(2,1,1); %upper bit
        h1 = pcolor( (Zenith_time_limTIME-floor(mean(Zenith_time))).*24 ,Zenith_height,Zenith_SNR_bp_limTIME');
        set(h1, 'EdgeColor', 'none'); %remove grid lines (which make it look hilarious)
        shading interp;
        lighting phong;
        colorbar;
        caxis([-0.1 0.1]);
        colormap(ax1,'gray');
        ylim([90 700]);
        % ylim([90 500]);
        % xlim([0 70]);
%         xlabel('Time in UT - 0 Hr at 0 UT May 7th');
        ylabel('Height (km)');
        title(['Time Cutout Range:',num2str(min(time_cutout_range)),':',num2str(max(time_cutout_range)),', Zenith SNR 2 Hr Highpass Filtered, Latitude = ',num2str(round(latMillstone,2)),' arcdeg and Longitude = ',num2str(round(longMillstone,2)),' arcdeg']);
        % set(gca, 'fontweight','bold', 'FontSize',12);
        % Xaxisvar = [-12, -8, -4, 0, 4, 8, 12, 16, 20, 24 28, 32, 36, 40, 44, 48]; %shows off important magnitude points
        %         Xaxisvar = [-12,-8,-4,0, 4, 8, 12, 16, 20, 24 28, 32, 36, 40, 44, 48, 52, 56, 60, 64, 68, 72]; %shows off important magnitude points
        Xaxisvar = [ (round((min(Zenith_time_limTIME)-floor(mean(Zenith_time)))*24)-mod(round((min(Zenith_time_limTIME)-floor(mean(Zenith_time)))*24),2)):2:...
            (round((max(Zenith_time_limTIME)-floor(mean(Zenith_time)))*24)-mod(round((max(Zenith_time_limTIME)-floor(mean(Zenith_time)))*24),2)) ];
        set(gca,'TickDir','out'); %put tick marks on outside of grid to help connect to time elements
        set(gca,'box','off'); %keeps tick marks from being on top and right side too
        set(gca,'xtick',Xaxisvar);
        set(gca,'fontweight',font_Weight, 'FontSize',font_Size);
        % 
        % %add on a vertical line that shows where we are at in time
        % 
        ax2 = subplot(2,1,2); %lower bit
        h = pcolor( (timeUnique_limTIME_anyAngleAvg -floor(mean(timeUnique))).*24 , avg_anyAngle_Range_Chunks_Long_Plot , sTECChunked_anyAngleAvg_limTIME(:,:)');% pseudocolor plot "stretched" to the grid
        set(h, 'EdgeColor', 'none'); %remove grid lines (which make it look hilarious)
        %             shading interp; %stuff for interpolation
        %             lighting phong; %stuff for interpolation
        caxis([-gif_sTEC_cap, gif_sTEC_cap]);
        %         colormap(flipud(colormap('hot'))) %heat map where 0 (low) is white
        colormap(ax2,'jet');
        h2 = colorbar; %shows the color bar
        ylabel(h2, 'delta-sTEC (TECU)','fontweight',font_Weight, 'FontSize',font_Size); %labels color bar
        set(gca,'xtick',Xaxisvar);

        avg_anyAngle_Range_Chunks_Long_Plot_autoTick = (ceil(max(avg_anyAngle_Range_Chunks_Long_Plot)) - floor(min(avg_anyAngle_Range_Chunks_Long_Plot)))/13; %tries to split the latitude range into 13 parts (based off of 180/15+1)
        if( avg_anyAngle_Range_Chunks_Long_Plot_autoTick > 25 )
            avg_anyAngle_Range_Chunks_Long_Plot_autoTick = 30; %sets the tick setting to 15 arcdegrees per tick
        elseif( avg_anyAngle_Range_Chunks_Long_Plot_autoTick > 10 )
            avg_anyAngle_Range_Chunks_Long_Plot_autoTick = 15; %sets the tick setting to 15 arcdegrees per tick
        elseif( avg_anyAngle_Range_Chunks_Long_Plot_autoTick > 5 )
            avg_anyAngle_Range_Chunks_Long_Plot_autoTick = 10; %sets the tick setting to 10 arcdegrees per tick
        elseif( avg_anyAngle_Range_Chunks_Long_Plot_autoTick > 2 )
            avg_anyAngle_Range_Chunks_Long_Plot_autoTick = 5; %sets the tick setting to 5 arcdegrees per tick
        elseif( avg_anyAngle_Range_Chunks_Long_Plot_autoTick > 1 )
            avg_anyAngle_Range_Chunks_Long_Plot_autoTick = 2; %sets the tick setting to 5 arcdegrees per tick
        elseif( avg_anyAngle_Range_Chunks_Long_Plot_autoTick >= 0.6 ) %0.6 because 15/25 = 0.6, so there will be enough 1 arcdeg ticks
            avg_anyAngle_Range_Chunks_Long_Plot_autoTick = 1; %sets the tick setting to 1 arcdegree per tick
        else
            avg_anyAngle_Range_Chunks_Long_Plot_autoTick = (max(plotLongRange) - min(plotLongRange))/13; %just goes for it if it's a super tiny range
        end
        Yaxisvar = round(floor(min(avg_anyAngle_Range_Chunks_Long_Plot)):avg_anyAngle_Range_Chunks_Long_Plot_autoTick:ceil(max(avg_anyAngle_Range_Chunks_Long_Plot)),2); %creates y ticks automagically
        set(gca,'ytick',Yaxisvar);

        hold on;
        %Now drawing line of interest
        if( strcmp(avg_anyAngle_Range_Chunks_Long_Plot_Name,'Longitude') == 1 ) %if true, longitude
            line([ min((timeUnique_limISR_Zenith -floor(mean(timeUnique))).*24) , max((timeUnique_limISR_Zenith -floor(mean(timeUnique))).*24) ],[longMillstone,longMillstone],'Color','k','LineWidth',.25);  %plots a point with a black line
        else %else latitude
            line([ min((timeUnique_limISR_Zenith -floor(mean(timeUnique))).*24) , max((timeUnique_limISR_Zenith -floor(mean(timeUnique))).*24) ],[latMillstone,latMillstone],'Color','k','LineWidth',.25);  %plots a point with a black line
        end

        xlabel(['Time in UT - 0 Hr on Day ',num2str(dateRange_zeroHr(2)),', ',num2str(dateRange_zeroHr(1)),' (hr)']); %old: ,'fontweight','bold','FontSize',12 for all
        ylabel([avg_anyAngle_Range_Chunks_Long_Plot_Name,' (arcdeg)']); %sped up significantly when merged above
        string_Title = ['TEC Averaged on Angle of ',num2str(round(avg_anyAngle*180/pi,2)),' deg and Width of ',num2str(round(avg_anyAngle_Width,2)),' arcdeg, Line Shows ',avg_anyAngle_Range_Chunks_Long_Plot_Name,' of Millstone Hill Zenith Beam']; %create title string
        if( FLG_geomagneticCoords == 1)
            string_Title = [string_Title,' [Mag Coord Aligned]']; %tack on mag coordinate warning
        end
        title(string_Title);
        set(gca,'fontweight',font_Weight, 'FontSize',font_Size);
        hold off;

    else % if 1, MISA ISR is plotted
        %PLOT W/ ZOOM
        if( isnan(latMillstone) ~= 1) %only run if latMillstone is relevant
            %prep figure #
            % figure('units','normalized','outerposition',[0 0 1 1]); %maximizes figure window
            figure(figNum);
            figNum = figNum+1;
            %Prep plot space
            ax1 = subplot(2,1,1); %upper bit
            h1 = pcolor( (MISA_time_limTIME-floor(mean(MISA_time))).*24 ,MISA_height,MISA_SNR_bp_limTIME');
            set(h1, 'EdgeColor', 'none'); %remove grid lines (which make it look hilarious)
            shading interp;
            lighting phong;
            colorbar;
            caxis([-0.1 0.1]);
            colormap(ax1,'gray');
            ylim([90 700]);
            % ylim([90 500]);
            % xlim([0 70]);
    %         xlabel('Time in UT - 0 Hr at 0 UT May 7th');
            ylabel('Height (km)');
            title(['Time Cutout Range:',num2str(min(time_cutout_range)),':',num2str(max(time_cutout_range)),', MISA SNR 2 Hr Highpass Filtered, Latitude = ',num2str(round(latMillstoneMISA,2)),' arcdeg and Longitude = ',num2str(round(longMillstoneMISA,2)),' arcdeg']);
            % set(gca, 'fontweight','bold', 'FontSize',12);
            % Xaxisvar = [-12, -8, -4, 0, 4, 8, 12, 16, 20, 24 28, 32, 36, 40, 44, 48]; %shows off important magnitude points
            %         Xaxisvar = [-12,-8,-4,0, 4, 8, 12, 16, 20, 24 28, 32, 36, 40, 44, 48, 52, 56, 60, 64, 68, 72]; %shows off important magnitude points
            Xaxisvar = [ (round((min(MISA_time_limTIME)-floor(mean(MISA_time)))*24)-mod(round((min(MISA_time_limTIME)-floor(mean(MISA_time)))*24),2)):2:...
                (round((max(MISA_time_limTIME)-floor(mean(MISA_time)))*24)-mod(round((max(MISA_time_limTIME)-floor(mean(MISA_time)))*24),2)) ];
            set(gca,'TickDir','out'); %put tick marks on outside of grid to help connect to time elements
            set(gca,'box','off'); %keeps tick marks from being on top and right side too
            set(gca,'xtick',Xaxisvar);
            set(gca,'fontweight',font_Weight, 'FontSize',font_Size);
            % 
            % %add on a vertical line that shows where we are at in time
            % 
            ax2 = subplot(2,1,2); %lower bit
            if( strcmp(avg_anyAngle_Range_Chunks_Long_Plot_Name,'Longitude') == 1 ) %if Longitude, use Zenith or MISA Longitude value
                avg_anyAngle_Zoom_Range = (longMillstoneMISA + avg_anyAngle_Zoom) >= avg_anyAngle_Range_Chunks_Long_Plot & (longMillstoneMISA - avg_anyAngle_Zoom) <= avg_anyAngle_Range_Chunks_Long_Plot; %gets the zoom range 
            else
                 avg_anyAngle_Zoom_Range = (latMillstoneMISA + avg_anyAngle_Zoom) >= avg_anyAngle_Range_Chunks_Long_Plot & (latMillstoneMISA - avg_anyAngle_Zoom) <= avg_anyAngle_Range_Chunks_Long_Plot; %gets the zoom range 
            end
            h = pcolor( (timeUnique_limTIME_anyAngleAvg -floor(mean(timeUnique))).*24 , avg_anyAngle_Range_Chunks_Long_Plot(avg_anyAngle_Zoom) , sTECChunked_anyAngleAvg_limTIME(:,avg_anyAngle_Zoom)');% pseudocolor plot "stretched" to the grid
            set(h, 'EdgeColor', 'none'); %remove grid lines (which make it look hilarious)
            %             shading interp; %stuff for interpolation
            %             lighting phong; %stuff for interpolation
            caxis([-gif_sTEC_cap, gif_sTEC_cap]);
            %         colormap(flipud(colormap('hot'))) %heat map where 0 (low) is white
            colormap(ax2,'jet');
            h2 = colorbar; %shows the color bar
            ylabel(h2, 'delta-sTEC (TECU)','fontweight',font_Weight, 'FontSize',font_Size); %labels color bar
            set(gca,'xtick',Xaxisvar);

            avg_anyAngle_Range_Chunks_Long_Plot_autoTick = (ceil(max(avg_anyAngle_Range_Chunks_Long_Plot(avg_anyAngle_Zoom))) - floor(min(avg_anyAngle_Range_Chunks_Long_Plot(avg_anyAngle_Zoom))))/13; %tries to split the latitude range into 13 parts (based off of 180/15+1)
            if( avg_anyAngle_Range_Chunks_Long_Plot_autoTick > 10 )
                avg_anyAngle_Range_Chunks_Long_Plot_autoTick = 15; %sets the tick setting to 15 arcdegrees per tick
            elseif( avg_anyAngle_Range_Chunks_Long_Plot_autoTick > 5 )
                avg_anyAngle_Range_Chunks_Long_Plot_autoTick = 10; %sets the tick setting to 10 arcdegrees per tick
            elseif( avg_anyAngle_Range_Chunks_Long_Plot_autoTick > 2 )
                avg_anyAngle_Range_Chunks_Long_Plot_autoTick = 5; %sets the tick setting to 5 arcdegrees per tick
            elseif( avg_anyAngle_Range_Chunks_Long_Plot_autoTick > 1 )
                avg_anyAngle_Range_Chunks_Long_Plot_autoTick = 2; %sets the tick setting to 5 arcdegrees per tick
            elseif( avg_anyAngle_Range_Chunks_Long_Plot_autoTick >= 0.6 ) %0.6 because 15/25 = 0.6, so there will be enough 1 arcdeg ticks
                avg_anyAngle_Range_Chunks_Long_Plot_autoTick = 1; %sets the tick setting to 1 arcdegree per tick
            else
                avg_anyAngle_Range_Chunks_Long_Plot_autoTick = (max(avg_anyAngle_Range_Chunks_Long_Plot(avg_anyAngle_Zoom_Range)) - min(avg_anyAngle_Range_Chunks_Long_Plot(avg_anyAngle_Zoom_Range)))/13; %just goes for it if it's a super tiny range
            end
            Yaxisvar = round(floor(min(avg_anyAngle_Range_Chunks_Long_Plot(avg_anyAngle_Zoom))):avg_anyAngle_Range_Chunks_Long_Plot_autoTick:ceil(max(avg_anyAngle_Range_Chunks_Long_Plot(avg_anyAngle_Zoom))),2); %creates y ticks automagically
            set(gca,'ytick',Yaxisvar);

            hold on;
            %Now drawing line of interest
            if( strcmp(avg_anyAngle_Range_Chunks_Long_Plot_Name,'Longitude') == 1 ) %if true, longitude
                line([ min((timeUnique_limISR_MISA -floor(mean(timeUnique))).*24) , max((timeUnique_limISR_MISA -floor(mean(timeUnique))).*24) ],[longMillstoneMISA,longMillstoneMISA],'Color','k','LineWidth',.25);  %plots a point with a black line
            else %else latitude
                line([ min((timeUnique_limISR_MISA -floor(mean(timeUnique))).*24) , max((timeUnique_limISR_MISA -floor(mean(timeUnique))).*24) ],[latMillstoneMISA,latMillstoneMISA],'Color','k','LineWidth',.25);  %plots a point with a black line
            end

            xlabel(['Time in UT - 0 Hr on Day ',num2str(dateRange_zeroHr(2)),', ',num2str(dateRange_zeroHr(1)),' (hr)']); %old: ,'fontweight','bold','FontSize',12 for all
            ylabel([avg_anyAngle_Range_Chunks_Long_Plot_Name,' (arcdeg)']); %sped up significantly when merged above
            string_Title = ['TEC Averaged on Angle of ',num2str(round(avg_anyAngle*180/pi,2)),' deg and Width of ',num2str(round(avg_anyAngle_Width,2)),' arcdeg, Line Shows ',avg_anyAngle_Range_Chunks_Long_Plot_Name,' of Millstone Hill MISA @ ',num2str(pointAltitude),' km with Zoom of +/- ',num2str(avg_anyAngle_Zoom),'deg around']; %create title string
            if( FLG_geomagneticCoords == 1)
                string_Title = [string_Title,' [Mag Coord Aligned]']; %tack on mag coordinate warning
            end
            title(string_Title);
            set(gca,'fontweight',font_Weight, 'FontSize',font_Size);
            hold off;
        end
        
        %PLOT W/O ZOOM
        %prep figure #
        % figure('units','normalized','outerposition',[0 0 1 1]); %maximizes figure window
        figure(figNum);
        figNum = figNum+1;
        %Prep plot space
        ax1 = subplot(2,1,1); %upper bit
        h1 = pcolor( (MISA_time_limTIME-floor(mean(MISA_time))).*24 ,MISA_height,MISA_SNR_bp_limTIME');
        set(h1, 'EdgeColor', 'none'); %remove grid lines (which make it look hilarious)
        shading interp;
        lighting phong;
        colorbar;
        caxis([-0.1 0.1]);
        colormap(ax1,'gray');
        ylim([90 700]);
        % ylim([90 500]);
        % xlim([0 70]);
%         xlabel('Time in UT - 0 Hr at 0 UT May 7th');
        ylabel('Height (km)');
        title(['Time Cutout Range:',num2str(min(time_cutout_range)),':',num2str(max(time_cutout_range)),', MISA SNR 2 Hr Highpass Filtered, Latitude = ',num2str(round(latMillstoneMISA,2)),' arcdeg and Longitude = ',num2str(round(longMillstoneMISA,2)),' arcdeg']);
        % set(gca, 'fontweight','bold', 'FontSize',12);
        % Xaxisvar = [-12, -8, -4, 0, 4, 8, 12, 16, 20, 24 28, 32, 36, 40, 44, 48]; %shows off important magnitude points
        %         Xaxisvar = [-12,-8,-4,0, 4, 8, 12, 16, 20, 24 28, 32, 36, 40, 44, 48, 52, 56, 60, 64, 68, 72]; %shows off important magnitude points
        Xaxisvar = [ (round((min(MISA_time_limTIME)-floor(mean(MISA_time)))*24)-mod(round((min(MISA_time_limTIME)-floor(mean(MISA_time)))*24),2)):2:...
            (round((max(MISA_time_limTIME)-floor(mean(MISA_time)))*24)-mod(round((max(MISA_time_limTIME)-floor(mean(MISA_time)))*24),2)) ];
        set(gca,'TickDir','out'); %put tick marks on outside of grid to help connect to time elements
        set(gca,'box','off'); %keeps tick marks from being on top and right side too
        set(gca,'xtick',Xaxisvar);
        set(gca,'fontweight',font_Weight, 'FontSize',font_Size);
        % 
        % %add on a vertical line that shows where we are at in time
        % 
        ax2 = subplot(2,1,2); %lower bit
        if( strcmp(avg_anyAngle_Range_Chunks_Long_Plot_Name,'Longitude') == 1 ) %if Longitude, use Zenith or MISA Longitude value
            avg_anyAngle_Zoom_Range = (longMillstoneMISA + avg_anyAngle_Zoom) >= avg_anyAngle_Range_Chunks_Long_Plot & (longMillstoneMISA - avg_anyAngle_Zoom) <= avg_anyAngle_Range_Chunks_Long_Plot; %gets the zoom range 
        else
             avg_anyAngle_Zoom_Range = (latMillstoneMISA + avg_anyAngle_Zoom) >= avg_anyAngle_Range_Chunks_Long_Plot & (latMillstoneMISA - avg_anyAngle_Zoom) <= avg_anyAngle_Range_Chunks_Long_Plot; %gets the zoom range 
        end
        h = pcolor( (timeUnique_limTIME_anyAngleAvg -floor(mean(timeUnique))).*24 , avg_anyAngle_Range_Chunks_Long_Plot , sTECChunked_anyAngleAvg_limTIME(:,:)');% pseudocolor plot "stretched" to the grid
        set(h, 'EdgeColor', 'none'); %remove grid lines (which make it look hilarious)
        %             shading interp; %stuff for interpolation
        %             lighting phong; %stuff for interpolation
        caxis([-gif_sTEC_cap, gif_sTEC_cap]);
        %         colormap(flipud(colormap('hot'))) %heat map where 0 (low) is white
        colormap(ax2,'jet');
        h2 = colorbar; %shows the color bar
        ylabel(h2, 'delta-sTEC (TECU)','fontweight',font_Weight, 'FontSize',font_Size); %labels color bar
        set(gca,'xtick',Xaxisvar);

        avg_anyAngle_Range_Chunks_Long_Plot_autoTick = (ceil(max(avg_anyAngle_Range_Chunks_Long_Plot)) - floor(min(avg_anyAngle_Range_Chunks_Long_Plot)))/13; %tries to split the latitude range into 13 parts (based off of 180/15+1)
        if( avg_anyAngle_Range_Chunks_Long_Plot_autoTick > 25 )
            avg_anyAngle_Range_Chunks_Long_Plot_autoTick = 30; %sets the tick setting to 15 arcdegrees per tick
        elseif( avg_anyAngle_Range_Chunks_Long_Plot_autoTick > 10 )
            avg_anyAngle_Range_Chunks_Long_Plot_autoTick = 15; %sets the tick setting to 15 arcdegrees per tick
        elseif( avg_anyAngle_Range_Chunks_Long_Plot_autoTick > 5 )
            avg_anyAngle_Range_Chunks_Long_Plot_autoTick = 10; %sets the tick setting to 10 arcdegrees per tick
        elseif( avg_anyAngle_Range_Chunks_Long_Plot_autoTick > 2 )
            avg_anyAngle_Range_Chunks_Long_Plot_autoTick = 5; %sets the tick setting to 5 arcdegrees per tick
        elseif( avg_anyAngle_Range_Chunks_Long_Plot_autoTick > 1 )
            avg_anyAngle_Range_Chunks_Long_Plot_autoTick = 2; %sets the tick setting to 5 arcdegrees per tick
        elseif( avg_anyAngle_Range_Chunks_Long_Plot_autoTick >= 0.6 ) %0.6 because 15/25 = 0.6, so there will be enough 1 arcdeg ticks
            avg_anyAngle_Range_Chunks_Long_Plot_autoTick = 1; %sets the tick setting to 1 arcdegree per tick
        else
            avg_anyAngle_Range_Chunks_Long_Plot_autoTick = (max(plotLongRange) - min(plotLongRange))/13; %just goes for it if it's a super tiny range
        end
        Yaxisvar = round(floor(min(avg_anyAngle_Range_Chunks_Long_Plot)):avg_anyAngle_Range_Chunks_Long_Plot_autoTick:ceil(max(avg_anyAngle_Range_Chunks_Long_Plot)),2); %creates y ticks automagically
        set(gca,'ytick',Yaxisvar);

        hold on;
        %Now drawing line of interest
        if( strcmp(avg_anyAngle_Range_Chunks_Long_Plot_Name,'Longitude') == 1 ) %if true, longitude
            line([ min((timeUnique_limISR_MISA -floor(mean(timeUnique))).*24) , max((timeUnique_limISR_MISA -floor(mean(timeUnique))).*24) ],[longMillstoneMISA,longMillstoneMISA],'Color','k','LineWidth',.25);  %plots a point with a black line
        else %else latitude
            line([ min((timeUnique_limISR_MISA -floor(mean(timeUnique))).*24) , max((timeUnique_limISR_MISA -floor(mean(timeUnique))).*24) ],[latMillstoneMISA,latMillstoneMISA],'Color','k','LineWidth',.25);  %plots a point with a black line
        end

        xlabel(['Time in UT - 0 Hr on Day ',num2str(dateRange_zeroHr(2)),', ',num2str(dateRange_zeroHr(1)),' (hr)']); %old: ,'fontweight','bold','FontSize',12 for all
        ylabel([avg_anyAngle_Range_Chunks_Long_Plot_Name,' (arcdeg)']); %sped up significantly when merged above
        string_Title = ['TEC Averaged on Angle of ',num2str(round(avg_anyAngle*180/pi,2)),' deg and Width of ',num2str(round(avg_anyAngle_Width,2)),' arcdeg, Line Shows ',avg_anyAngle_Range_Chunks_Long_Plot_Name,' of Millstone Hill MISA @ ',num2str(pointAltitude),' km']; %create title string
        if( FLG_geomagneticCoords == 1)
            string_Title = [string_Title,' [Mag Coord Aligned]']; %tack on mag coordinate warning
        end
        title(string_Title);
        set(gca,'fontweight',font_Weight, 'FontSize',font_Size);
        hold off;
        
    end
    
    
end %*********END OF FLG_AVG_ANYANGLE_PLOT_TIMECUTOUT*********


%% Average Any Angle, IT'S SCARGLIN TIME TIME LIMITED AND BAND PASSED
if( FLG_AVG_ANYANGLE_TimeMatch_n_BP_Scargle == 1 )
    
    if(FLG_AVG_ANYANGLE_Scargle_FFT == 0 ) %scargle option
    
        if( FLG_AVG_ANYANGLE_ZENITHORMISA == 0 ) %Zenith option
            if( strcmp(avg_anyAngle_Range_Chunks_Long_Plot_Name,'Longitude') == 1 ) %if true, longitude
                jk = find( min(abs( avg_anyAngle_Range_Chunks_Long_Plot - longMillstone )) == abs( avg_anyAngle_Range_Chunks_Long_Plot - longMillstone ) ); %get the location closest
            else %else vs latitude
                jk = find( min(abs( avg_anyAngle_Range_Chunks_Long_Plot - latMillstone )) == abs( avg_anyAngle_Range_Chunks_Long_Plot - latMillstone ) ); %get the location closest
            end
            jm = isnan(sTECChunked_anyAngleAvg(:,jk)); %get the NaN
            dataTime = (timeUnique_limISR_Zenith(~jm)' -floor(mean(timeUnique))).*24; %time data time w/o NaNs
%             sTECChunked_anyAngleAvg_atMillstone = interp1( dataTime , sTECChunked_anyAngleAvg(~jm,jk) , (timeUnique_limISR_Zenith' -floor(mean(timeUnique))).*24 ,'linear','extrap'); %interpolate extra data in
            
            sTECChunked_anyAngleAvg_atMillstone = sTECChunked_anyAngleAvg(~jm,jk); %ignore that interp, just ditch the NaNs
            
%             tmin = find( min( abs(min(time_cutout_range) - (timeUnique_limISR_Zenith - floor(mean(timeUnique))).*24 ) ) == abs(min(time_cutout_range) - (timeUnique_limISR_Zenith- floor(mean(timeUnique))).*24 ) ); %find min timeUnique that fits time limit
%             tmax = find( min( abs(max(time_cutout_range) - (timeUnique_limISR_Zenith -  floor(mean(timeUnique))).*24 ) ) == abs(max(time_cutout_range) - (timeUnique_limISR_Zenith- floor(mean(timeUnique))).*24 ) ); %find max timeUnique that fits time limit
%             timeUnique_limTIME_Zenith = timeUnique_limISR_Zenith(tmin:tmax); %cut out time that matches time limits (Zenith or MISA have been expanded to this time cadence)

            tmin = find( min( abs(min(time_cutout_range) - dataTime) ) == abs(min(time_cutout_range) - dataTime ) ); %find min timeUnique that fits time limit
            tmax = find( min( abs(max(time_cutout_range) - dataTime) ) == abs(max(time_cutout_range) - dataTime ) ); %find max timeUnique that fits time limit
            dataTime_limTIME = dataTime(tmin:tmax); %cut out time that matches time limits (Zenith or MISA have been expanded to this time cadence)
            sTECChunked_anyAngleAvg_atMillstone_limTIME = sTECChunked_anyAngleAvg_atMillstone(tmin:tmax); %time limited

            bs=(1/2);   % highpass cuttoff frequency (not sure what to make of it)
            %I think it is related to lf above (which is 1/2 )

            n=42; % order of the Hamming window used in the custom function
            % c = 3.32*pi; %Hamming constant
            % M = n/2;
            % bandwidth = c/M;
            %The above was included to investigate the bandwidth - 1/2. Related to
            %bs/lf?

            fp=bs; % stop band stoppin freq (or pass band passing freq, depending on how you look at it)
            timeDelta = mean(timeUnique_limTIME_Zenith(2:end)-timeUnique_limTIME_Zenith(1:end-1))*24; %hrs already
            f= 1/(timeDelta); % the sampling frequency, based off of the time delta calc'd

            wp=2*fp/f; % Normalizing the frequencies (Matlab is 0 to 1)

            %Calculation of filter coefficients
            % [b,a]=fir1(n,wp,'high'); %This takes the order and lower cut off
            %Uses the default hamming window
            % frequency ORIG
            W = hann(n+1);
            [b,a] = fir1(n,wp,'high',W);
            %Applys Hanning Window for shiggles and giggles


            sTECChunked_anyAngleAvg_atMillstone_bp = filtfilt(b,a,sTECChunked_anyAngleAvg_atMillstone); %Appies the filter

%             Ross_Lombscargle_optimized([(timeUnique_limISR_Zenith' -floor(mean(timeUnique))).*24 , sTECChunked_anyAngleAvg_atMillstone_bp ]); 
            Ross_Lombscargle_optimized([dataTime , sTECChunked_anyAngleAvg_atMillstone_bp ]); 
            fid = fopen('power_freq.txt', 'r');
                f_nm = textscan(fid, '%f %f %f %f %f', 'headerlines', 15);
                P_Zenith_SNR=f_nm{:,2};
                F_Zenith_SNR=f_nm{:,1};
                T_Zenith_SNR=f_nm{:,5}*60; % in minutes
                gs_Zenith_SNR=f_nm{:,4};
                max_Zenith_SNR=T_Zenith_SNR(P_Zenith_SNR==max(P_Zenith_SNR));
                K=length(P_Zenith_SNR);
                gf_Zenith_SNR=f_nm{:,3};
            fclose('all');

            %prep figure #
            % figure('units','normalized','outerposition',[0 0 1 1]); %maximizes figure window
            figure(figNum);
            figNum = figNum+1;

            plot(T_Zenith_SNR,P_Zenith_SNR,'b');
            xlim([0 plot_Freq_Lim]);
            hold on;
            % plot(T_SNR_08apr08_10_or,P_SNR_08apr08_10_or,'r');
            % hold on;
            plot(T_Zenith_SNR,gf_Zenith_SNR,'LineWidth',1.5,'color',[0.5 0.5 0.5]);
            xlim([0 plot_Freq_Lim]);  
            xlabel('Periods (min)');
            ylabel('Normalized Power');
            string_Title = ['Periodogram of TEC Avg''d on Angle of ',num2str(round(avg_anyAngle*180/pi,2)),' deg and Width of ',num2str(round(avg_anyAngle_Width,2)),...
                ' arcdeg, High-Passed @ ',num2str(round(1/bs,2)),' hr at Zenith ',avg_anyAngle_Range_Chunks_Long_Plot_Name];
            if( FLG_geomagneticCoords == 1)
                string_Title = [string_Title,' [Mag Coord Aligned]']; %tack on mag coordinate warning
            end
            title(string_Title);
            set(gca,'fontweight',font_Weight, 'FontSize',font_Size);
            set(gca,'xtick',Xaxisvar_SCARGLE);

            sTECChunked_anyAngleAvg_atMillstone_limTIME_bp = filtfilt(b,a,sTECChunked_anyAngleAvg_atMillstone_limTIME); %Appies the filter

%             Ross_Lombscargle_optimized([(timeUnique_limTIME_Zenith' -floor(mean(timeUnique))).*24 , sTECChunked_anyAngleAvg_atMillstone_limTIME_bp ]); 
            Ross_Lombscargle_optimized([dataTime_limTIME , sTECChunked_anyAngleAvg_atMillstone_limTIME_bp ]); 
            fid = fopen('power_freq.txt', 'r');
                f_nm = textscan(fid, '%f %f %f %f %f', 'headerlines', 15);
                P_Zenith_SNR=f_nm{:,2};
                F_Zenith_SNR=f_nm{:,1};
                T_Zenith_SNR=f_nm{:,5}*60; % in minutes
                gs_Zenith_SNR=f_nm{:,4};
                max_Zenith_SNR=T_Zenith_SNR(P_Zenith_SNR==max(P_Zenith_SNR));
                K=length(P_Zenith_SNR);
                gf_Zenith_SNR=f_nm{:,3};
            fclose('all');

            %prep figure #
            % figure('units','normalized','outerposition',[0 0 1 1]); %maximizes figure window
            figure(figNum);
            figNum = figNum+1;

            plot(T_Zenith_SNR,P_Zenith_SNR,'b');
            xlim([0 plot_Freq_Lim]);
            hold on;
            % plot(T_SNR_08apr08_10_or,P_SNR_08apr08_10_or,'r');
            % hold on;
            plot(T_Zenith_SNR,gf_Zenith_SNR,'LineWidth',1.5,'color',[0.5 0.5 0.5]);
            xlim([0 plot_Freq_Lim]);  
            xlabel('Periods (min)');
            ylabel('Normalized Power');
            string_Title = ['Periodogram of TEC Avg''d on Angle of ',num2str(round(avg_anyAngle*180/pi,2)),' deg and Width of ',num2str(round(avg_anyAngle_Width,2)),...
                ' arcdeg, High-Passed @ ',num2str(round(1/bs,2)),' hr & Time Lim ',num2str(min(time_cutout_range)),'-',num2str(max(time_cutout_range)),...
                ' hrs at Zenith ',avg_anyAngle_Range_Chunks_Long_Plot_Name];
            if( FLG_geomagneticCoords == 1)
                string_Title = [string_Title,' [Mag Coord Aligned]']; %tack on mag coordinate warning
            end
            title(string_Title);
            set(gca,'fontweight',font_Weight, 'FontSize',font_Size);
            set(gca,'xtick',Xaxisvar_SCARGLE);
            
            
            %=========COMPARE TO ZENITH AND MISA ISR===========
            % SCARGLE IT
            Ross_Lombscargle_optimized([dataTime_limTIME , sTECChunked_anyAngleAvg_atMillstone_limTIME_bp  ]); 
            fid = fopen('power_freq.txt', 'r');
                f_nm = textscan(fid, '%f %f %f %f %f', 'headerlines', 15);
                P_Zenith_SNR=f_nm{:,2};
                F_Zenith_SNR=f_nm{:,1};
                T_Zenith_SNR=f_nm{:,5}*60; % in minutes
                gs_Zenith_SNR=f_nm{:,4};
                max_Zenith_SNR=T_Zenith_SNR(P_Zenith_SNR==max(P_Zenith_SNR));
                K=length(P_Zenith_SNR);
                gf_Zenith_SNR=f_nm{:,3};
            fclose('all');

            %prep figure #
            % figure('units','normalized','outerposition',[0 0 1 1]); %maximizes figure window
            figure(figNum);
            figNum = figNum+1;
            subplot(3,1,1);
            plot(T_Zenith_SNR,P_Zenith_SNR,'b');
            xlim([0 plot_Freq_Lim]);
            hold on;
            % plot(T_SNR_08apr08_10_or,P_SNR_08apr08_10_or,'r');
            % hold on;
            plot(T_Zenith_SNR,gf_Zenith_SNR,'LineWidth',1.5,'color',[0.5 0.5 0.5]);
            xlim([0 plot_Freq_Lim]);  
            xlabel('Periods (min)');
            ylabel('Normalized Power');
            string_Title = ['Periodogram of TEC Avg''d on Angle of ',num2str(round(avg_anyAngle*180/pi,2)),' deg and Width of ',num2str(round(avg_anyAngle_Width,2)),...
                ' arcdeg, High-Passed @ ',num2str(round(1/bs,2)),' hr at Zenith ',avg_anyAngle_Range_Chunks_Long_Plot_Name, ' and Cut to Hr ',num2str(min(time_cutout_range)),' to ',num2str(max(time_cutout_range))];
            if( FLG_geomagneticCoords == 1)
                string_Title = [string_Title,' [Mag Coord Aligned]']; %tack on mag coordinate warning
            end
            title(string_Title);
            set(gca,'fontweight',font_Weight, 'FontSize',font_Size);
            set(gca,'xtick',Xaxisvar_SCARGLE);

            tmin_Zenith = find(Zenith_time_limTIME(1) == Zenith_time); %get time indicies to limit the Zenith or MISA ISR data
            tmax_Zenith = find(Zenith_time_limTIME(end) == Zenith_time);
            tmin_MISA = find(MISA_time_limTIME(1) == MISA_time);
            tmax_MISA = find(MISA_time_limTIME(end) == MISA_time);
            % SCARGLE IT (ZENITH!)
            Ross_Lombscargle_optimized([(Zenith_time_limTIME-floor(mean(Zenith_time))).*24 ,Zenith_SNR_threeHun_AVGD(tmin_Zenith:tmax_Zenith)  ]); 
            fid = fopen('power_freq.txt', 'r');
                f_nm = textscan(fid, '%f %f %f %f %f', 'headerlines', 15);
                P_Zenith_SNR=f_nm{:,2};
                F_Zenith_SNR=f_nm{:,1};
                T_Zenith_SNR=f_nm{:,5}*60; % in minutes
                gs_Zenith_SNR=f_nm{:,4};
                max_Zenith_SNR=T_Zenith_SNR(P_Zenith_SNR==max(P_Zenith_SNR));
                K=length(P_Zenith_SNR);
                gf_Zenith_SNR=f_nm{:,3};
            fclose('all');

            %prep 
            subplot(3,1,2);
            plot(T_Zenith_SNR,P_Zenith_SNR,'b');
            xlim([0 plot_Freq_Lim]);
            hold on;
            % plot(T_SNR_08apr08_10_or,P_SNR_08apr08_10_or,'r');
            % hold on;
            plot(T_Zenith_SNR,gf_Zenith_SNR,'LineWidth',1.5,'color',[0.5 0.5 0.5]);
            xlim([0 plot_Freq_Lim]);  
            xlabel('Periods (min)');
            ylabel('Normalized Power');
            title(['Zenith ISR avg''d +/-',num2str(Zenith_avgRange),' km around ',num2str(round(Zenith_height(Zenith_threeHunIndex))),' km - ',num2str(latMillstone),' deg lat/',num2str(longMillstone),'deg long - Hours ',num2str(min(time_cutout_range)),' to ',num2str(max(time_cutout_range)),' on Day ',num2str(floor(mean(Zenith_time))),', 2013 - Lomb-Scargle Freq']);
            set(gca,'xtick',Xaxisvar_SCARGLE);
            set(gca,'fontweight',font_Weight, 'FontSize',font_Size);

            % SCARGLE IT (MISA!)
            Ross_Lombscargle_optimized([(MISA_time_limTIME -floor(mean(MISA_time))).*24 ,MISA_SNR_threeHun_AVGD(tmin_MISA:tmax_MISA)  ]); 
            fid = fopen('power_freq.txt', 'r');
                f_nm = textscan(fid, '%f %f %f %f %f', 'headerlines', 15);
                P_Zenith_SNR=f_nm{:,2};
                F_Zenith_SNR=f_nm{:,1};
                T_Zenith_SNR=f_nm{:,5}*60; % in minutes
                gs_Zenith_SNR=f_nm{:,4};
                max_Zenith_SNR=T_Zenith_SNR(P_Zenith_SNR==max(P_Zenith_SNR));
                K=length(P_Zenith_SNR);
                gf_Zenith_SNR=f_nm{:,3};
            fclose('all');

            %prep
            subplot(3,1,3);
            plot(T_Zenith_SNR,P_Zenith_SNR,'b');
            xlim([0 plot_Freq_Lim]);
            hold on;
            % plot(T_SNR_08apr08_10_or,P_SNR_08apr08_10_or,'r');
            % hold on;
            plot(T_Zenith_SNR,gf_Zenith_SNR,'LineWidth',1.5,'color',[0.5 0.5 0.5]);
            xlim([0 plot_Freq_Lim]);  
            xlabel('Periods (min)');
            ylabel('Normalized Power');
            title(['MISA ISR avg''d +/-',num2str(MISA_avgRange),' km around ',num2str(round(MISA_height(MISA_threeHunIndex))),' km - ',num2str(latMillstoneMISA),' deg lat/',num2str(longMillstoneMISA),'deg long - Hours ',num2str(min(time_cutout_range)),' to ',num2str(max(time_cutout_range)),' on Day ',num2str(floor(mean(MISA_time))),', 2013 - Lomb-Scargle Freq']);
            set(gca,'xtick',Xaxisvar_SCARGLE);
            set(gca,'fontweight',font_Weight, 'FontSize',font_Size);

            
            %======PUT the AVG_anyangle on the same time cadence as ISR Zenith & MISA=======
            sTECChunked_anyAngleAvg_atMillstone_5min_Zenith = zeros(length(Zenith_time),1); %percolate
            Zenith_time_delta = mean(Zenith_time(2:end)-Zenith_time(1:end-1)); %days, delta of time between readings
            dataTime_dayConv = dataTime/24+floor(mean(Zenith_time)); %days, convert from hr to day
            for(i = 1:length(Zenith_time) )

                if(i == 1)
                    jold = find( min(abs((Zenith_time(i)-Zenith_time_delta) - dataTime_dayConv)) == abs((Zenith_time(i)-Zenith_time_delta) - dataTime_dayConv) ); %get the matching time 5 min prev. (since it's an integration up to the time stamp given)
                    jcurrent = find( min(abs((Zenith_time(i)) - dataTime_dayConv)) == abs((Zenith_time(i)) - dataTime_dayConv) ); %get the matching time

                    sTECChunked_anyAngleAvg_atMillstone_5min_Zenith(i) = mean(sTECChunked_anyAngleAvg_atMillstone(jold:jcurrent)); % average over the ISR time period
                else
                    jold = jcurrent; %set the old time
                    jcurrent = find( min(abs((Zenith_time(i)) - dataTime_dayConv)) == abs((Zenith_time(i)) - dataTime_dayConv) ); %get the matching time

                    sTECChunked_anyAngleAvg_atMillstone_5min_Zenith(i) = mean(sTECChunked_anyAngleAvg_atMillstone(jold+1:jcurrent)); % average over the ISR time period (increment jold by 1 so not using same data)
                end

            end
            
            jm = isnan(sTECChunked_anyAngleAvg_atMillstone_5min_Zenith); %find the NaNs
            if( sum(isnan(sTECChunked_anyAngleAvg_atMillstone_5min_Zenith)) > 0 )
                 sTECChunked_anyAngleAvg_atMillstone_5min_Zenith = interp1( dataTime(~jm) , sTECChunked_anyAngleAvg_atMillstone_5min_Zenith(~jm) , dataTime ,'linear','extrap'); %interpolate extra data in
            end

            sTECChunked_anyAngleAvg_atMillstone_5min_MISA = zeros(length(MISA_time),1);
            MISA_time_delta = mean(MISA_time(2:end)-MISA_time(1:end-1)); %days, delta of time between readings

            for(i = 1:length(MISA_time) )

                if(i == 1) %special case to start - prev 5 min
                    jold = find( min(abs((MISA_time(i)-MISA_time_delta) - dataTime_dayConv)) == abs((MISA_time(i)-MISA_time_delta) - dataTime_dayConv) ); %get the matching time 5 min prev.
                    jcurrent = find( min(abs((MISA_time(i)) - dataTime_dayConv)) == abs((MISA_time(i)) - dataTime_dayConv) ); %get the matching time

                    sTECChunked_anyAngleAvg_atMillstone_5min_MISA(i) = mean(sTECChunked_anyAngleAvg_atMillstone(jold:jcurrent)); % average over the ISR time period
                else
                    jold = jcurrent; %set the old time
                    jcurrent = find( min(abs((MISA_time(i)) - dataTime_dayConv)) == abs((MISA_time(i)) - dataTime_dayConv) ); %get the matching time

                    sTECChunked_anyAngleAvg_atMillstone_5min_MISA(i) = mean(sTECChunked_anyAngleAvg_atMillstone(jold+1:jcurrent)); % average over the ISR time period (increment jold by 1 so not using same data)
                end

            end
            
            jm = isnan(sTECChunked_anyAngleAvg_atMillstone_5min_MISA); %find the NaNs
            if( sum(isnan(sTECChunked_anyAngleAvg_atMillstone_5min_MISA)) > 0 )
                 sTECChunked_anyAngleAvg_atMillstone_5min_MISA = interp1( dataTime(~jm) , sTECChunked_anyAngleAvg_atMillstone_5min_MISA(~jm) , dataTime ,'linear','extrap'); %interpolate extra data in
            end
            
            
            sTECChunked_anyAngleAvg_atMillstone_5min_Zenith_bp = filtfilt(b,a,sTECChunked_anyAngleAvg_atMillstone_5min_Zenith); %Appies the filter
            sTECChunked_anyAngleAvg_atMillstone_5min_MISA_bp = filtfilt(b,a,sTECChunked_anyAngleAvg_atMillstone_5min_MISA); %Appies the filter

            R_Zenith = corrcoef(sTECChunked_anyAngleAvg_atMillstone_5min_Zenith_bp(tmin_Zenith:tmax_Zenith),Zenith_SNR_threeHun_AVGD(tmin_Zenith:tmax_Zenith));
            fprintf('Correlation coeff between Zenith ISR and anyAngleAVG-sTEC 5 Min Averaged over Millstone w/ High-pass at time cutout range %d to %d: %f\n',min(time_cutout_range),max(time_cutout_range),R_Zenith(2));

            R_MISA = corrcoef(sTECChunked_anyAngleAvg_atMillstone_5min_MISA_bp(tmin_MISA:tmax_MISA),MISA_SNR_threeHun_AVGD(tmin_MISA:tmax_MISA));
            fprintf('Correlation coeff between MISA ISR and anyAngleAVG-sTEC 5 Min Averaged over Millstone w/ High-pass at time cutout range %d to %d: %f\n\n',min(time_cutout_range),max(time_cutout_range),R_MISA(2));
            
            %COMPARE TO ZENITH AND MISA ISR
            % SCARGLE IT
            Ross_Lombscargle_optimized([(Zenith_time_limTIME-floor(mean(Zenith_time))).*24 , sTECChunked_anyAngleAvg_atMillstone_5min_Zenith_bp(tmin_Zenith:tmax_Zenith)  ]); 
            fid = fopen('power_freq.txt', 'r');
                f_nm = textscan(fid, '%f %f %f %f %f', 'headerlines', 15);
                P_Zenith_SNR=f_nm{:,2};
                F_Zenith_SNR=f_nm{:,1};
                T_Zenith_SNR=f_nm{:,5}*60; % in minutes
                gs_Zenith_SNR=f_nm{:,4};
                max_Zenith_SNR=T_Zenith_SNR(P_Zenith_SNR==max(P_Zenith_SNR));
                K=length(P_Zenith_SNR);
                gf_Zenith_SNR=f_nm{:,3};
            fclose('all');

            %prep figure #
            % figure('units','normalized','outerposition',[0 0 1 1]); %maximizes figure window
            figure(figNum);
            figNum = figNum+1;
            subplot(3,1,1);
            plot(T_Zenith_SNR,P_Zenith_SNR,'b');
            xlim([0 plot_Freq_Lim]);
            hold on;
            % plot(T_SNR_08apr08_10_or,P_SNR_08apr08_10_or,'r');
            % hold on;
            plot(T_Zenith_SNR,gf_Zenith_SNR,'LineWidth',1.5,'color',[0.5 0.5 0.5]);
            xlim([0 plot_Freq_Lim]);  
            xlabel('Periods (min)');
            ylabel('Normalized Power');
            string_Title = ['Periodogram of TEC Avg''d on Angle of ',num2str(round(avg_anyAngle*180/pi,2)),' deg and Width of ',num2str(round(avg_anyAngle_Width,2)),...
                ' arcdeg, High-Passed @ ',num2str(round(1/bs,2)),' hr at Zenith ',avg_anyAngle_Range_Chunks_Long_Plot_Name,' and on Zenith ISR Cadence Cut to Hr ',num2str(min(time_cutout_range)),' to ',num2str(max(time_cutout_range))];
            if( FLG_geomagneticCoords == 1)
                string_Title = [string_Title,' [Mag Coord Aligned]']; %tack on mag coordinate warning
            end
            title(string_Title);
            set(gca,'fontweight',font_Weight, 'FontSize',font_Size);
            set(gca,'xtick',Xaxisvar_SCARGLE);

            % SCARGLE IT (ZENITH!)
            Ross_Lombscargle_optimized([(Zenith_time_limTIME-floor(mean(Zenith_time))).*24 ,Zenith_SNR_threeHun_AVGD(tmin_Zenith:tmax_Zenith)  ]); 
            fid = fopen('power_freq.txt', 'r');
                f_nm = textscan(fid, '%f %f %f %f %f', 'headerlines', 15);
                P_Zenith_SNR=f_nm{:,2};
                F_Zenith_SNR=f_nm{:,1};
                T_Zenith_SNR=f_nm{:,5}*60; % in minutes
                gs_Zenith_SNR=f_nm{:,4};
                max_Zenith_SNR=T_Zenith_SNR(P_Zenith_SNR==max(P_Zenith_SNR));
                K=length(P_Zenith_SNR);
                gf_Zenith_SNR=f_nm{:,3};
            fclose('all');

            %prep 
            subplot(3,1,2);
            plot(T_Zenith_SNR,P_Zenith_SNR,'b');
            xlim([0 plot_Freq_Lim]);
            hold on;
            % plot(T_SNR_08apr08_10_or,P_SNR_08apr08_10_or,'r');
            % hold on;
            plot(T_Zenith_SNR,gf_Zenith_SNR,'LineWidth',1.5,'color',[0.5 0.5 0.5]);
            xlim([0 plot_Freq_Lim]);  
            xlabel('Periods (min)');
            ylabel('Normalized Power');
            title(['Zenith ISR avg''d +/-',num2str(Zenith_avgRange),' km around ',num2str(round(Zenith_height(Zenith_threeHunIndex))),' km - ',num2str(latMillstone),' deg lat/',num2str(longMillstone),'deg long - Hours ',num2str(min(time_cutout_range)),' to ',num2str(max(time_cutout_range)),' on Day ',num2str(floor(mean(Zenith_time))),', 2013 - Lomb-Scargle Freq']);
            set(gca,'xtick',Xaxisvar_SCARGLE);
            set(gca,'fontweight',font_Weight, 'FontSize',font_Size);

            % SCARGLE IT (MISA!)
            Ross_Lombscargle_optimized([(MISA_time_limTIME -floor(mean(MISA_time))).*24 ,MISA_SNR_threeHun_AVGD(tmin_MISA:tmax_MISA)  ]); 
            fid = fopen('power_freq.txt', 'r');
                f_nm = textscan(fid, '%f %f %f %f %f', 'headerlines', 15);
                P_Zenith_SNR=f_nm{:,2};
                F_Zenith_SNR=f_nm{:,1};
                T_Zenith_SNR=f_nm{:,5}*60; % in minutes
                gs_Zenith_SNR=f_nm{:,4};
                max_Zenith_SNR=T_Zenith_SNR(P_Zenith_SNR==max(P_Zenith_SNR));
                K=length(P_Zenith_SNR);
                gf_Zenith_SNR=f_nm{:,3};
            fclose('all');

            %prep
            subplot(3,1,3);
            plot(T_Zenith_SNR,P_Zenith_SNR,'b');
            xlim([0 plot_Freq_Lim]);
            hold on;
            % plot(T_SNR_08apr08_10_or,P_SNR_08apr08_10_or,'r');
            % hold on;
            plot(T_Zenith_SNR,gf_Zenith_SNR,'LineWidth',1.5,'color',[0.5 0.5 0.5]);
            xlim([0 plot_Freq_Lim]);  
            xlabel('Periods (min)');
            ylabel('Normalized Power');
            title(['MISA ISR avg''d +/-',num2str(MISA_avgRange),' km around ',num2str(round(MISA_height(MISA_threeHunIndex))),' km - ',num2str(latMillstoneMISA),' deg lat/',num2str(longMillstoneMISA),'deg long - Hours ',num2str(min(time_cutout_range)),' to ',num2str(max(time_cutout_range)),' on Day ',num2str(floor(mean(MISA_time))),', 2013 - Lomb-Scargle Freq']);
            set(gca,'xtick',Xaxisvar_SCARGLE);
            set(gca,'fontweight',font_Weight, 'FontSize',font_Size);

            %Testing phase difference of Zenith at 300 km using Cross Power 
            % Spectral Density=============================================
    % 
    %         hr_timeUnique_limISR_Zenith = (timeUnique_limISR_Zenith -floor(mean(timeUnique))).*24; %get the hour version
    %         Fs = 1/(hr_timeUnique_limISR_Zenith(178)-hr_timeUnique_limISR_Zenith(177)); %get the delta time (hr)
    %         nfft = 512; %even, size of vectors will be nfft/2 + 1;
    %         window = hamming(110);
    %         Cxy_Zenith = zeros(nfft/2+1,length(Zenith_height));
    %         F_Zenith = zeros(nfft/2+1,length(Zenith_height));
    %         Axy_Zenith = zeros(nfft/2+1,length(Zenith_height));
    %         Pxy_Zenith = zeros(nfft/2+1,length(Zenith_height));
    %         for i=1:size(sTECChunked_anyAngleAvg,2)
    %             [Cxy_Zenith(:,i),F_Zenith(:,i)] = cpsd(filtfilt(b,a,sTECChunked_anyAngleAvg(:,i)),filtfilt(b,a,sTECChunked_anyAngleAvg_atMillstone),window,100,nfft,Fs);
    %             %I just CAN'T get cpsd to work with a vector. DESPITE THE DOCUMENTATION
    %             %SAYING IT IS GOOD TO GO. It errors on things that are correct.......
    %             Axy_Zenith(:,i)=angle(Cxy_Zenith(:,i))*180/pi;
    %             Pxy_Zenith(:,i)=abs(Cxy_Zenith(:,i));
    %         end
    %         % 
    %         % gg = repmat(Zenith_SNR_47,1,length(Zenith_height)); %DOES NOT WORK BUT
    %         % % SHOULD
    %         % [Cxy_Zenith2,F_Zenith2]  = cpsd(abs(Zenith_SNR(1:end,1:end)),abs(gg(1:end,1:end)),window,100,nfft,Fs);
    % 
    %         F_Zenith = F_Zenith(:,end); %Puts F_Zenith at the end of the array used.
    % 
    % %         Pxy_log_Zenith(:,:) = log10(Pxy_Zenith(:,:)); %Takes the log
    % 
    %         F_Zenith(:) = 60./F_Zenith(:); %Converts some units... somewhere (sec -> min I think)
    %         
    %         
    %         figure(figNum);
    %         figNum = figNum+1;
    %         subplot(2,1,1);
    %         pcolor(F_Zenith,avg_anyAngle_Range_Chunks_Long_Plot,Axy_Zenith');
    %         shading interp;
    %         lighting phong;
    %         colorbar;
    %         xlim([0 plot_Freq_Lim]);
    % %         ylim([200 500]);
    %         ylabel([avg_anyAngle_Range_Chunks_Long_Plot_Name,' (arcdeg)']);
    %         xlabel('Period (min)');
    %         title(['Cross spectral phase differences between delta-sTEC AVG''d at each ',avg_anyAngle_Range_Chunks_Long_Plot_Name,' band and delta-sTEC AVG''d at Haystack Observatory']);
    %         set(gca,'fontweight',font_Weight, 'FontSize',font_Size);
    %         
    %         subplot(2,1,2);
    %         pcolor(F_Zenith,avg_anyAngle_Range_Chunks_Long_Plot,Pxy_Zenith');
    %         shading interp;
    %         lighting phong;
    %         colorbar;
    %         xlim([0 plot_Freq_Lim]);
    % %         ylim([200 500]);
    % %         caxis([0 0.1]);
    %         ylabel([avg_anyAngle_Range_Chunks_Long_Plot_Name,' (arcdeg)']);
    %         xlabel('Period (min)');
    %         title(['Cross spectral amplitude of delta-sTEC AVG''d at each ',avg_anyAngle_Range_Chunks_Long_Plot_Name,' band and delta-sTEC AVG''d at Haystack Observatory']);
    %         set(gca,'fontweight',font_Weight, 'FontSize',font_Size);

        else %MISA option!

            if( strcmp(avg_anyAngle_Range_Chunks_Long_Plot_Name,'Longitude') == 1 ) %if true, longitude
                jk = find( min(abs( avg_anyAngle_Range_Chunks_Long_Plot - longMillstoneMISA )) == abs( avg_anyAngle_Range_Chunks_Long_Plot - longMillstoneMISA ) ); %get the location closest
            else %else vs latitude
                jk = find( min(abs( avg_anyAngle_Range_Chunks_Long_Plot - latMillstoneMISA )) == abs( avg_anyAngle_Range_Chunks_Long_Plot - latMillstoneMISA ) ); %get the location closest
            end
            jm = isnan(sTECChunked_anyAngleAvg(:,jk)); %get the NaN
            dataTime = (timeUnique_limISR_MISA(~jm)' -floor(mean(timeUnique))).*24; %time of data that is not NaN
%             sTECChunked_anyAngleAvg_atMillstone = interp1( dataTime , sTECChunked_anyAngleAvg(~jm,jk) , (timeUnique_limISR_Zenith' -floor(mean(timeUnique))).*24 ,'linear','extrap'); %interpolate extra data in
            
            sTECChunked_anyAngleAvg_atMillstone = sTECChunked_anyAngleAvg(~jm,jk); %ignore that interp, just ditch the NaNs
            
%             tmin = find( min( abs(min(time_cutout_range) - (timeUnique_limISR_Zenith - floor(mean(timeUnique))).*24 ) ) == abs(min(time_cutout_range) - (timeUnique_limISR_Zenith- floor(mean(timeUnique))).*24 ) ); %find min timeUnique that fits time limit
%             tmax = find( min( abs(max(time_cutout_range) - (timeUnique_limISR_Zenith -  floor(mean(timeUnique))).*24 ) ) == abs(max(time_cutout_range) - (timeUnique_limISR_Zenith- floor(mean(timeUnique))).*24 ) ); %find max timeUnique that fits time limit
%             timeUnique_limTIME_Zenith = timeUnique_limISR_Zenith(tmin:tmax); %cut out time that matches time limits (Zenith or MISA have been expanded to this time cadence)

            tmin = find( min( abs(min(time_cutout_range) - dataTime) ) == abs(min(time_cutout_range) - dataTime ) ); %find min timeUnique that fits time limit
            tmax = find( min( abs(max(time_cutout_range) - dataTime) ) == abs(max(time_cutout_range) - dataTime ) ); %find max timeUnique that fits time limit
            dataTime_limTIME = dataTime(tmin:tmax); %cut out time that matches time limits (Zenith or MISA have been expanded to this time cadence)
            sTECChunked_anyAngleAvg_atMillstone_limTIME = sTECChunked_anyAngleAvg_atMillstone(tmin:tmax); %time limited

            bs=(1/2);   % highpass cuttoff frequency (not sure what to make of it)
            %I think it is related to lf above (which is 1/2 )

            n=42; % order of the Hamming window used in the custom function
            % c = 3.32*pi; %Hamming constant
            % M = n/2;
            % bandwidth = c/M;
            %The above was included to investigate the bandwidth - 1/2. Related to
            %bs/lf?

            fp=bs; % stop band stoppin freq (or pass band passing freq, depending on how you look at it)
            timeDelta = mean(timeUnique_limTIME_MISA(2:end)-timeUnique_limTIME_MISA(1:end-1))*24; %hrs already
            f= 1/(timeDelta); % the sampling frequency, based off of the time delta calc'd

            wp=2*fp/f; % Normalizing the frequencies (Matlab is 0 to 1)

            %Calculation of filter coefficients
            % [b,a]=fir1(n,wp,'high'); %This takes the order and lower cut off
            %Uses the default hamming window
            % frequency ORIG
            W = hann(n+1);
            [b,a] = fir1(n,wp,'high',W);
            %Applys Hanning Window for shiggles and giggles

            sTECChunked_anyAngleAvg_atMillstone_bp = filtfilt(b,a,sTECChunked_anyAngleAvg_atMillstone); %Appies the filter


%             Ross_Lombscargle_optimized([(timeUnique_limISR_MISA' -floor(mean(timeUnique))).*24 , sTECChunked_anyAngleAvg_atMillstone_bp ]); 
            Ross_Lombscargle_optimized([dataTime , sTECChunked_anyAngleAvg_atMillstone_bp ]); 
            fid = fopen('power_freq.txt', 'r');
                f_nm = textscan(fid, '%f %f %f %f %f', 'headerlines', 15);
                P_Zenith_SNR=f_nm{:,2};
                F_Zenith_SNR=f_nm{:,1};
                T_Zenith_SNR=f_nm{:,5}*60; % in minutes
                gs_Zenith_SNR=f_nm{:,4};
                max_Zenith_SNR=T_Zenith_SNR(P_Zenith_SNR==max(P_Zenith_SNR));
                K=length(P_Zenith_SNR);
                gf_Zenith_SNR=f_nm{:,3};
            fclose('all');

            %prep figure #
            % figure('units','normalized','outerposition',[0 0 1 1]); %maximizes figure window
            figure(figNum);
            figNum = figNum+1;

            plot(T_Zenith_SNR,P_Zenith_SNR,'b');
            xlim([0 plot_Freq_Lim]);
            hold on;
            % plot(T_SNR_08apr08_10_or,P_SNR_08apr08_10_or,'r');
            % hold on;
            plot(T_Zenith_SNR,gf_Zenith_SNR,'LineWidth',1.5,'color',[0.5 0.5 0.5]);
            xlim([0 plot_Freq_Lim]);  
            xlabel('Periods (min)');
            ylabel('Normalized Power');
            string_Title = ['Periodogram of TEC Avg''d on Angle of ',num2str(round(avg_anyAngle*180/pi,2)),' deg and Width of ',num2str(round(avg_anyAngle_Width,2)),...
                ' arcdeg, High-Passed @ ',num2str(round(1/bs,2)),' hr',...
                ' at MISA ',avg_anyAngle_Range_Chunks_Long_Plot_Name];
            if( FLG_geomagneticCoords == 1)
                string_Title = [string_Title,' [Mag Coord Aligned]']; %tack on mag coordinate warning
            end
            title(string_Title);
            set(gca,'fontweight',font_Weight, 'FontSize',font_Size);
            set(gca,'xtick',Xaxisvar_SCARGLE);

            sTECChunked_anyAngleAvg_atMillstone_limTIME_bp = filtfilt(b,a,sTECChunked_anyAngleAvg_atMillstone_limTIME); %Appies the filter

%             Ross_Lombscargle_optimized([(timeUnique_limTIME_MISA' -floor(mean(timeUnique))).*24 , sTECChunked_anyAngleAvg_atMillstone_limTIME_bp ]); 
            Ross_Lombscargle_optimized([dataTime_limTIME , sTECChunked_anyAngleAvg_atMillstone_limTIME_bp ]); 
            fid = fopen('power_freq.txt', 'r');
                f_nm = textscan(fid, '%f %f %f %f %f', 'headerlines', 15);
                P_Zenith_SNR=f_nm{:,2};
                F_Zenith_SNR=f_nm{:,1};
                T_Zenith_SNR=f_nm{:,5}*60; % in minutes
                gs_Zenith_SNR=f_nm{:,4};
                max_Zenith_SNR=T_Zenith_SNR(P_Zenith_SNR==max(P_Zenith_SNR));
                K=length(P_Zenith_SNR);
                gf_Zenith_SNR=f_nm{:,3};
            fclose('all');

            %prep figure #
            % figure('units','normalized','outerposition',[0 0 1 1]); %maximizes figure window
            figure(figNum);
            figNum = figNum+1;

            plot(T_Zenith_SNR,P_Zenith_SNR,'b');
            xlim([0 plot_Freq_Lim]);
            hold on;
            % plot(T_SNR_08apr08_10_or,P_SNR_08apr08_10_or,'r');
            % hold on;
            plot(T_Zenith_SNR,gf_Zenith_SNR,'LineWidth',1.5,'color',[0.5 0.5 0.5]);
            xlim([0 plot_Freq_Lim]);  
            xlabel('Periods (min)');
            ylabel('Normalized Power');
            string_Title = ['Periodogram of TEC Avg''d on Angle of ',num2str(round(avg_anyAngle*180/pi,2)),' deg and Width of ',num2str(round(avg_anyAngle_Width,2)),...
                ' arcdeg, High-Passed @ ',num2str(round(1/bs,2)),' hr & Time Lim ',num2str(min(time_cutout_range)),'-',num2str(max(time_cutout_range)),...
                ' hrs at MISA ',avg_anyAngle_Range_Chunks_Long_Plot_Name];
            if( FLG_geomagneticCoords == 1)
                string_Title = [string_Title,' [Mag Coord Aligned]']; %tack on mag coordinate warning
            end
            title(string_Title);
            set(gca,'fontweight',font_Weight, 'FontSize',font_Size);
            set(gca,'xtick',Xaxisvar_SCARGLE);
        end
    
    else % FFT Option
        
        if( FLG_AVG_ANYANGLE_ZENITHORMISA == 0 ) %Zenith option
            if( strcmp(avg_anyAngle_Range_Chunks_Long_Plot_Name,'Longitude') == 1 ) %if true, longitude
                jk = find( min(abs( avg_anyAngle_Range_Chunks_Long_Plot - longMillstone )) == abs( avg_anyAngle_Range_Chunks_Long_Plot - longMillstone ) ); %get the location closest
            else %else vs latitude
                jk = find( min(abs( avg_anyAngle_Range_Chunks_Long_Plot - latMillstone )) == abs( avg_anyAngle_Range_Chunks_Long_Plot - latMillstone ) ); %get the location closest
            end
            jm = isnan(sTECChunked_anyAngleAvg(:,jk)); %get the NaN
            dataTime = (timeUnique_limISR_Zenith(~jm)' -floor(mean(timeUnique))).*24; %time data time
            sTECChunked_anyAngleAvg_atMillstone = interp1( dataTime , sTECChunked_anyAngleAvg(~jm,jk) , (timeUnique_limISR_Zenith' -floor(mean(timeUnique))).*24 ,'linear','extrap'); %interpolate extra data in

            tmin = find( min( abs(min(time_cutout_range) - (timeUnique_limISR_Zenith - floor(mean(timeUnique))).*24 ) ) == abs(min(time_cutout_range) - (timeUnique_limISR_Zenith- floor(mean(timeUnique))).*24 ) ); %find min timeUnique that fits time limit
            tmax = find( min( abs(max(time_cutout_range) - (timeUnique_limISR_Zenith -  floor(mean(timeUnique))).*24 ) ) == abs(max(time_cutout_range) - (timeUnique_limISR_Zenith- floor(mean(timeUnique))).*24 ) ); %find max timeUnique that fits time limit
            timeUnique_limTIME_Zenith = timeUnique_limISR_Zenith(tmin:tmax); %cut out time that matches time limits (Zenith or MISA have been expanded to this time cadence)
            sTECChunked_anyAngleAvg_atMillstone_limTIME = sTECChunked_anyAngleAvg_atMillstone(tmin:tmax); %time limited

            bs=(1/2);   % highpass cuttoff frequency (not sure what to make of it)
            %I think it is related to lf above (which is 1/2 )

            n=42; % order of the Hamming window used in the custom function
            % c = 3.32*pi; %Hamming constant
            % M = n/2;
            % bandwidth = c/M;
            %The above was included to investigate the bandwidth - 1/2. Related to
            %bs/lf?

            fp=bs; % stop band stoppin freq (or pass band passing freq, depending on how you look at it)
            timeDelta = mean(timeUnique_limTIME_Zenith(2:end)-timeUnique_limTIME_Zenith(1:end-1))*24; %hrs already
            f= 1/(timeDelta); % the sampling frequency, based off of the time delta calc'd

            wp=2*fp/f; % Normalizing the frequencies (Matlab is 0 to 1)

            %Calculation of filter coefficients
            % [b,a]=fir1(n,wp,'high'); %This takes the order and lower cut off
            %Uses the default hamming window
            % frequency ORIG
            W = hann(n+1);
            [b,a] = fir1(n,wp,'high',W);
            %Applys Hanning Window for shiggles and giggles


            sTECChunked_anyAngleAvg_atMillstone_bp = filtfilt(b,a,sTECChunked_anyAngleAvg_atMillstone); %Appies the filter

            %prep figure #
            % figure('units','normalized','outerposition',[0 0 1 1]); %maximizes figure window
            totalcount = size(sTECChunked_anyAngleAvg_atMillstone_bp,1); % # of samples taken
            deltatime = timeDelta; %sample interval in time - units hours *)
            deltafreq = 1/(totalcount*deltatime); %frequency resolution per above sampling interval in time
            frequencysamples = ((0:(totalcount-1))*deltafreq)';
            % deltafreq = 1/((2^nextpow2(totalcount)*FFT_coeff)*deltatime);
            % frequencysamples = ((0:2^nextpow2(totalcount)*FFT_coeff-1)*deltafreq)';
            %frequencysamples = linspace(0,ceil((totalcount-1)*deltafreq),totalcount)';
            % spectrummax = frequencysamples(totalcount/2 - 1);
            % timedomain = zeros(totalcount,1);

            frequencydomain = fft(sTECChunked_anyAngleAvg_atMillstone_bp); %,2^nextpow2(totalcount)*FFT_coeff
            powerspectrum = sqrt(abs(frequencydomain.*conj(frequencydomain)))./totalcount; %properly normalized now
            
            gf=(1-(.05/length(powerspectrum))^(1/(length(powerspectrum)-1))); %for gf and gs see pg 934-935 of the journal paper "SPECTRUM: SPECTRAL ANALYSIS OF UNEVENLY SPACED PALEOCLIMATIC TIME SERIES" by SCHULZ and STATTEGGER
%         gs=0.4*gf; %Computers & Geosciences Vol. 23, No. 9, pp. 929-945, 1997 

            figure(figNum);
            figNum = figNum+1;
            plot(1./frequencysamples*60,powerspectrum)
            hold on;
            plot(1./frequencysamples*60,repmat(gf,size(powerspectrum)),'LineWidth',1.5,'color',[0.5 0.5 0.5]);
            string_Title = ['DFT of TEC Avg''d on Angle of ',num2str(round(avg_anyAngle*180/pi,2)),' deg and Width of ',num2str(round(avg_anyAngle_Width,2)),...
                ' arcdeg, High-passed @ ',num2str(round(1/bs,2)),' hr at Zenith ',avg_anyAngle_Range_Chunks_Long_Plot_Name];
            if( FLG_geomagneticCoords == 1)
                string_Title = [string_Title,' [Mag Coord Aligned]']; %tack on mag coordinate warning
            end
            title(string_Title);
            xlabel('Period (min)');
            ylabel('Normalized Power');
            axis([0, plot_Freq_Lim, -inf, inf])
            set(gca,'fontweight',font_Weight, 'FontSize',font_Size);
            % axis 'auto y'

            sTECChunked_anyAngleAvg_atMillstone_limTIME_bp = filtfilt(b,a,sTECChunked_anyAngleAvg_atMillstone_limTIME); %Appies the filter

            totalcount = size(sTECChunked_anyAngleAvg_atMillstone_limTIME_bp,1); % # of samples taken
            deltatime = timeDelta; %sample interval in time - units hours *)
            deltafreq = 1/(totalcount*deltatime); %frequency resolution per above sampling interval in time
            frequencysamples = ((0:(totalcount-1))*deltafreq)';
            % deltafreq = 1/((2^nextpow2(totalcount)*FFT_coeff)*deltatime);
            % frequencysamples = ((0:2^nextpow2(totalcount)*FFT_coeff-1)*deltafreq)';
            %frequencysamples = linspace(0,ceil((totalcount-1)*deltafreq),totalcount)';
            % spectrummax = frequencysamples(totalcount/2 - 1);
            % timedomain = zeros(totalcount,1);

            frequencydomain = fft(sTECChunked_anyAngleAvg_atMillstone_limTIME_bp); %,2^nextpow2(totalcount)*FFT_coeff
            powerspectrum = sqrt(abs(frequencydomain.*conj(frequencydomain)))./totalcount; %properly normalized now
            
            gf=(1-(.05/length(powerspectrum))^(1/(length(powerspectrum)-1))); %for gf and gs see pg 934-935 of the journal paper "SPECTRUM: SPECTRAL ANALYSIS OF UNEVENLY SPACED PALEOCLIMATIC TIME SERIES" by SCHULZ and STATTEGGER
%         gs=0.4*gf; %Computers & Geosciences Vol. 23, No. 9, pp. 929-945, 1997 

            figure(figNum);
            figNum = figNum+1;
            plot(1./frequencysamples*60,powerspectrum) 
            hold on;
            plot(1./frequencysamples*60,repmat(gf,size(powerspectrum)),'LineWidth',1.5,'color',[0.5 0.5 0.5]);
            string_Title = ['DFT of TEC Avg''d on Angle of ',num2str(round(avg_anyAngle*180/pi,2)),' deg and Width of ',num2str(round(avg_anyAngle_Width,2)),...
                ' arcdeg, High-passed @ ',num2str(round(1/bs,2)),' hr & Time Lim ',num2str(min(time_cutout_range)),'-',num2str(max(time_cutout_range)),...
                ' hrs at Zenith ',avg_anyAngle_Range_Chunks_Long_Plot_Name];
            if( FLG_geomagneticCoords == 1)
                string_Title = [string_Title,' [Mag Coord Aligned]']; %tack on mag coordinate warning
            end
            title(string_Title);
            xlabel('Period (min)');
            ylabel('Normalized Power');
            axis([0, plot_Freq_Lim, -inf, inf])
            set(gca,'fontweight',font_Weight, 'FontSize',font_Size);
            % axis 'auto y'



            %Testing phase difference of Zenith at 300 km using Cross Power 
            % Spectral Density=============================================
    % 
    %         hr_timeUnique_limISR_Zenith = (timeUnique_limISR_Zenith -floor(mean(timeUnique))).*24; %get the hour version
    %         Fs = 1/(hr_timeUnique_limISR_Zenith(178)-hr_timeUnique_limISR_Zenith(177)); %get the delta time (hr)
    %         nfft = 512; %even, size of vectors will be nfft/2 + 1;
    %         window = hamming(110);
    %         Cxy_Zenith = zeros(nfft/2+1,length(Zenith_height));
    %         F_Zenith = zeros(nfft/2+1,length(Zenith_height));
    %         Axy_Zenith = zeros(nfft/2+1,length(Zenith_height));
    %         Pxy_Zenith = zeros(nfft/2+1,length(Zenith_height));
    %         for i=1:size(sTECChunked_anyAngleAvg,2)
    %             [Cxy_Zenith(:,i),F_Zenith(:,i)] = cpsd(filtfilt(b,a,sTECChunked_anyAngleAvg(:,i)),filtfilt(b,a,sTECChunked_anyAngleAvg_atMillstone),window,100,nfft,Fs);
    %             %I just CAN'T get cpsd to work with a vector. DESPITE THE DOCUMENTATION
    %             %SAYING IT IS GOOD TO GO. It errors on things that are correct.......
    %             Axy_Zenith(:,i)=angle(Cxy_Zenith(:,i))*180/pi;
    %             Pxy_Zenith(:,i)=abs(Cxy_Zenith(:,i));
    %         end
    %         % 
    %         % gg = repmat(Zenith_SNR_47,1,length(Zenith_height)); %DOES NOT WORK BUT
    %         % % SHOULD
    %         % [Cxy_Zenith2,F_Zenith2]  = cpsd(abs(Zenith_SNR(1:end,1:end)),abs(gg(1:end,1:end)),window,100,nfft,Fs);
    % 
    %         F_Zenith = F_Zenith(:,end); %Puts F_Zenith at the end of the array used.
    % 
    % %         Pxy_log_Zenith(:,:) = log10(Pxy_Zenith(:,:)); %Takes the log
    % 
    %         F_Zenith(:) = 60./F_Zenith(:); %Converts some units... somewhere (sec -> min I think)
    %         
    %         
    %         figure(figNum);
    %         figNum = figNum+1;
    %         subplot(2,1,1);
    %         pcolor(F_Zenith,avg_anyAngle_Range_Chunks_Long_Plot,Axy_Zenith');
    %         shading interp;
    %         lighting phong;
    %         colorbar;
    %         xlim([0 plot_Freq_Lim]);
    % %         ylim([200 500]);
    %         ylabel([avg_anyAngle_Range_Chunks_Long_Plot_Name,' (arcdeg)']);
    %         xlabel('Period (min)');
    %         title(['Cross spectral phase differences between delta-sTEC AVG''d at each ',avg_anyAngle_Range_Chunks_Long_Plot_Name,' band and delta-sTEC AVG''d at Haystack Observatory']);
    %         set(gca,'fontweight',font_Weight, 'FontSize',font_Size);
    %         
    %         subplot(2,1,2);
    %         pcolor(F_Zenith,avg_anyAngle_Range_Chunks_Long_Plot,Pxy_Zenith');
    %         shading interp;
    %         lighting phong;
    %         colorbar;
    %         xlim([0 plot_Freq_Lim]);
    % %         ylim([200 500]);
    % %         caxis([0 0.1]);
    %         ylabel([avg_anyAngle_Range_Chunks_Long_Plot_Name,' (arcdeg)']);
    %         xlabel('Period (min)');
    %         title(['Cross spectral amplitude of delta-sTEC AVG''d at each ',avg_anyAngle_Range_Chunks_Long_Plot_Name,' band and delta-sTEC AVG''d at Haystack Observatory']);
    %         set(gca,'fontweight',font_Weight, 'FontSize',font_Size);

        else %MISA option!

            if( strcmp(avg_anyAngle_Range_Chunks_Long_Plot_Name,'Longitude') == 1 ) %if true, longitude
                jk = find( min(abs( avg_anyAngle_Range_Chunks_Long_Plot - longMillstoneMISA )) == abs( avg_anyAngle_Range_Chunks_Long_Plot - longMillstoneMISA ) ); %get the location closest
            else %else vs latitude
                jk = find( min(abs( avg_anyAngle_Range_Chunks_Long_Plot - latMillstoneMISA )) == abs( avg_anyAngle_Range_Chunks_Long_Plot - latMillstoneMISA ) ); %get the location closest
            end
            jm = isnan(sTECChunked_anyAngleAvg(:,jk)); %get the NaN
            dataTime = (timeUnique_limISR_MISA(~jm)' -floor(mean(timeUnique))).*24; %time of data that is not NaN
            sTECChunked_anyAngleAvg_atMillstone = interp1( dataTime , sTECChunked_anyAngleAvg(~jm,jk) , (timeUnique_limISR_MISA' -floor(mean(timeUnique))).*24 ,'linear','extrap'); %interpolate extra data in
            
            tmin = find( min( abs(min(time_cutout_range) - (timeUnique_limISR_MISA - floor(mean(timeUnique))).*24 ) ) == abs(min(time_cutout_range) - (timeUnique_limISR_MISA- floor(mean(timeUnique))).*24 ) ); %find min timeUnique that fits time limit
            tmax = find( min( abs(max(time_cutout_range) - (timeUnique_limISR_MISA -  floor(mean(timeUnique))).*24 ) ) == abs(max(time_cutout_range) - (timeUnique_limISR_MISA- floor(mean(timeUnique))).*24 ) ); %find max timeUnique that fits time limit
            timeUnique_limTIME_MISA = timeUnique_limISR_MISA(tmin:tmax); %cut out time that matches time limits (Zenith or MISA have been expanded to this time cadence)
            sTECChunked_anyAngleAvg_atMillstone_limTIME = sTECChunked_anyAngleAvg_atMillstone(tmin:tmax); %time limited

            bs=(1/2);   % highpass cuttoff frequency (not sure what to make of it)
            %I think it is related to lf above (which is 1/2 )

            n=42; % order of the Hamming window used in the custom function
            % c = 3.32*pi; %Hamming constant
            % M = n/2;
            % bandwidth = c/M;
            %The above was included to investigate the bandwidth - 1/2. Related to
            %bs/lf?

            fp=bs; % stop band stoppin freq (or pass band passing freq, depending on how you look at it)
            timeDelta = mean(timeUnique_limTIME_MISA(2:end)-timeUnique_limTIME_MISA(1:end-1))*24; %hrs already
            f= 1/(timeDelta); % the sampling frequency, based off of the time delta calc'd

            wp=2*fp/f; % Normalizing the frequencies (Matlab is 0 to 1)

            %Calculation of filter coefficients
            % [b,a]=fir1(n,wp,'high'); %This takes the order and lower cut off
            %Uses the default hamming window
            % frequency ORIG
            W = hann(n+1);
            [b,a] = fir1(n,wp,'high',W);
            %Applys Hanning Window for shiggles and giggles

            sTECChunked_anyAngleAvg_atMillstone_bp = filtfilt(b,a,sTECChunked_anyAngleAvg_atMillstone); %Appies the filter

            %prep figure #
            % figure('units','normalized','outerposition',[0 0 1 1]); %maximizes figure window
            totalcount = size(sTECChunked_anyAngleAvg_atMillstone_bp,1); % # of samples taken
            deltatime = timeDelta; %sample interval in time - units hours *)
            deltafreq = 1/(totalcount*deltatime); %frequency resolution per above sampling interval in time
            frequencysamples = ((0:(totalcount-1))*deltafreq)';
            % deltafreq = 1/((2^nextpow2(totalcount)*FFT_coeff)*deltatime);
            % frequencysamples = ((0:2^nextpow2(totalcount)*FFT_coeff-1)*deltafreq)';
            %frequencysamples = linspace(0,ceil((totalcount-1)*deltafreq),totalcount)';
            % spectrummax = frequencysamples(totalcount/2 - 1);
            % timedomain = zeros(totalcount,1);

            frequencydomain = fft(sTECChunked_anyAngleAvg_atMillstone_bp); %,2^nextpow2(totalcount)*FFT_coeff
            powerspectrum = sqrt(abs(frequencydomain.*conj(frequencydomain)))./totalcount; %properly normalized now
            
            gf=(1-(.05/length(powerspectrum))^(1/(length(powerspectrum)-1))); %for gf and gs see pg 934-935 of the journal paper "SPECTRUM: SPECTRAL ANALYSIS OF UNEVENLY SPACED PALEOCLIMATIC TIME SERIES" by SCHULZ and STATTEGGER
%         gs=0.4*gf; %Computers & Geosciences Vol. 23, No. 9, pp. 929-945, 1997 

            figure(figNum);
            figNum = figNum+1;
            plot(1./frequencysamples*60,powerspectrum)
            hold on;
            plot(1./frequencysamples*60,repmat(gf,size(powerspectrum)),'LineWidth',1.5,'color',[0.5 0.5 0.5]);
            string_Title = ['DFT of TEC Avg''d on Angle of ',num2str(round(avg_anyAngle*180/pi,2)),' deg and Width of ',num2str(round(avg_anyAngle_Width,2)),...
                ' arcdeg, High-passed @ ',num2str(round(1/bs,2)),' hr at MISA ',avg_anyAngle_Range_Chunks_Long_Plot_Name];
            if( FLG_geomagneticCoords == 1)
                string_Title = [string_Title,' [Mag Coord Aligned]']; %tack on mag coordinate warning
            end
            title(string_Title);
            xlabel('Period (min)');
            ylabel('Normalized Power');
            axis([0, plot_Freq_Lim, -inf, inf])
            set(gca,'fontweight',font_Weight, 'FontSize',font_Size);
            % axis 'auto y'

            sTECChunked_anyAngleAvg_atMillstone_limTIME_bp = filtfilt(b,a,sTECChunked_anyAngleAvg_atMillstone_limTIME); %Appies the filter

            totalcount = size(sTECChunked_anyAngleAvg_atMillstone_limTIME_bp,1); % # of samples taken
            deltatime = timeDelta; %sample interval in time - units hours *)
            deltafreq = 1/(totalcount*deltatime); %frequency resolution per above sampling interval in time
            frequencysamples = ((0:(totalcount-1))*deltafreq)';
            % deltafreq = 1/((2^nextpow2(totalcount)*FFT_coeff)*deltatime);
            % frequencysamples = ((0:2^nextpow2(totalcount)*FFT_coeff-1)*deltafreq)';
            %frequencysamples = linspace(0,ceil((totalcount-1)*deltafreq),totalcount)';
            % spectrummax = frequencysamples(totalcount/2 - 1);
            % timedomain = zeros(totalcount,1);

            frequencydomain = fft(sTECChunked_anyAngleAvg_atMillstone_limTIME_bp); %,2^nextpow2(totalcount)*FFT_coeff
            powerspectrum = sqrt(abs(frequencydomain.*conj(frequencydomain)))./totalcount; %properly normalized now
            
            gf=(1-(.05/length(powerspectrum))^(1/(length(powerspectrum)-1))); %for gf and gs see pg 934-935 of the journal paper "SPECTRUM: SPECTRAL ANALYSIS OF UNEVENLY SPACED PALEOCLIMATIC TIME SERIES" by SCHULZ and STATTEGGER
%         gs=0.4*gf; %Computers & Geosciences Vol. 23, No. 9, pp. 929-945, 1997 

            figure(figNum);
            figNum = figNum+1;
            plot(1./frequencysamples*60,powerspectrum)
            hold on;
            plot(1./frequencysamples*60,repmat(gf,size(powerspectrum)),'LineWidth',1.5,'color',[0.5 0.5 0.5]);
            string_Title = ['DFT of TEC Avg''d on Angle of ',num2str(round(avg_anyAngle*180/pi,2)),' deg and Width of ',num2str(round(avg_anyAngle_Width,2)),...
                ' arcdeg, High-passed @ ',num2str(round(1/bs,2)),' hr & Time Lim ',num2str(min(time_cutout_range)),'-',num2str(max(time_cutout_range)),...
                ' hrs at MISA ',avg_anyAngle_Range_Chunks_Long_Plot_Name];
            if( FLG_geomagneticCoords == 1)
                string_Title = [string_Title,' [Mag Coord Aligned]']; %tack on mag coordinate warning
            end
            title(string_Title);
            xlabel('Period (min)');
            ylabel('Normalized Power');
            axis([0, plot_Freq_Lim, -inf, inf])
            set(gca,'fontweight',font_Weight, 'FontSize',font_Size);
            % axis 'auto y'
        end
        
    end

end %*********END OF FLG_AVG_LONG_TimeMatch_n_BP_Scargle*********


%% IT'S SCARGLIN TIME LONG AVG'D (so time vs latitude is plot)
if( FLG_AVG_LONG_Scargle == 1 )

jk = find( min(abs( pplatChunksPlot - latMillstone )) == abs( pplatChunksPlot - latMillstone ) ); %get the location closest
jm = isnan(sTECChunked_longavg(:,jk)); %get the NaN
dataTime = (timeUnique_limISR_MISA(~jm)' -floor(mean(timeUnique))).*24; %time data time
sTECChunked_longavg_atMillstone = interp1( dataTime , sTECChunked_longavg(~jm,jk) , (timeUnique_limISR_MISA' -floor(mean(timeUnique))).*24 ,'linear','extrap'); %interpolate extra data in

Ross_Lombscargle_optimized([(timeUnique_limISR_MISA' -floor(mean(timeUnique))).*24 , sTECChunked_longavg_atMillstone ]); 
fid = fopen('power_freq.txt', 'r');
    f_nm = textscan(fid, '%f %f %f %f %f', 'headerlines', 15);
    P_Zenith_SNR=f_nm{:,2};
    F_Zenith_SNR=f_nm{:,1};
    T_Zenith_SNR=f_nm{:,5}*60; % in minutes
    gs_Zenith_SNR=f_nm{:,4};
    max_Zenith_SNR=T_Zenith_SNR(P_Zenith_SNR==max(P_Zenith_SNR));
    K=length(P_Zenith_SNR);
    gf_Zenith_SNR=f_nm{:,3};
fclose('all');

%prep figure #
% figure('units','normalized','outerposition',[0 0 1 1]); %maximizes figure window
figure(figNum);
figNum = figNum+1;

plot(T_Zenith_SNR,P_Zenith_SNR,'b');
xlim([0 plot_Freq_Lim]);
hold on;
% plot(T_SNR_08apr08_10_or,P_SNR_08apr08_10_or,'r');
% hold on;
plot(T_Zenith_SNR,gf_Zenith_SNR,'LineWidth',1.5,'color',[0.5 0.5 0.5]);
xlim([0 plot_Freq_Lim]);  
xlabel('Periods (min)');
ylabel('Normalized Power');
title(['Longitude Avg''d at Millstone Latitude of ',num2str(latMillstone),' deg lat with longitude AVD''d +/-',num2str(pplongAvgDeg),'deg on Day ',num2str(floor(mean(time))),', 2013 - SCARGLED']);
set(gca,'fontweight',font_Weight, 'FontSize',font_Size);
set(gca,'xtick',Xaxisvar_SCARGLE);

end %*********END OF FLG_AVG_LONG_Scargle*********


%% IT'S SCARGLIN TIME (ZENITH AVG'd) - ZENITH
if( FLG_ISRdata_Zenith_Scargle == 1 )

Ross_Lombscargle_optimized([(Zenith_time-floor(mean(Zenith_time))).*24, Zenith_SNR_threeHun_AVGD ]); 
fid = fopen('power_freq.txt', 'r');
    f_nm = textscan(fid, '%f %f %f %f %f', 'headerlines', 15);
    P_Zenith_SNR=f_nm{:,2};
    F_Zenith_SNR=f_nm{:,1};
    T_Zenith_SNR=f_nm{:,5}*60; % in minutes
    gs_Zenith_SNR=f_nm{:,4};
    max_Zenith_SNR=T_Zenith_SNR(P_Zenith_SNR==max(P_Zenith_SNR));
    K=length(P_Zenith_SNR);
    gf_Zenith_SNR=f_nm{:,3};
fclose('all');

%prep figure #
% figure('units','normalized','outerposition',[0 0 1 1]); %maximizes figure window
figure(figNum);
figNum = figNum+1;

plot(T_Zenith_SNR,P_Zenith_SNR,'b');
xlim([0 plot_Freq_Lim]);
hold on;
% plot(T_SNR_08apr08_10_or,P_SNR_08apr08_10_or,'r');
% hold on;
plot(T_Zenith_SNR,gf_Zenith_SNR,'LineWidth',1.5,'color',[0.5 0.5 0.5]);
xlim([0 plot_Freq_Lim]);  
xlabel('Periods (min)');
ylabel('Normalized Power');
title(['Zenith SNR BP at ',num2str(latMillstone),' deg lat/',num2str(longMillstone),'deg long on Day ',num2str(floor(mean(time))),', 2013 - SCARGLED']);
set(gca,'fontweight',font_Weight, 'FontSize',font_Size);
set(gca,'xtick',Xaxisvar_SCARGLE);

end %*********END OF FLG_ISRdata_Zenith_Scargle*********


%% IT'S SCARGLIN TIME (MISA AVG'd) - MISA
if( FLG_ISRdata_MISA_Scargle == 1 )

Ross_Lombscargle_optimized([(MISA_time-floor(mean(MISA_time))).*24, MISA_SNR_threeHun_AVGD ]); 
fid = fopen('power_freq.txt', 'r');
    f_nm = textscan(fid, '%f %f %f %f %f', 'headerlines', 15);
    P_Zenith_SNR=f_nm{:,2};
    F_Zenith_SNR=f_nm{:,1};
    T_Zenith_SNR=f_nm{:,5}*60; % in minutes
    gs_Zenith_SNR=f_nm{:,4};
    max_Zenith_SNR=T_Zenith_SNR(P_Zenith_SNR==max(P_Zenith_SNR));
    K=length(P_Zenith_SNR);
    gf_Zenith_SNR=f_nm{:,3};
fclose('all');

%prep figure #
% figure('units','normalized','outerposition',[0 0 1 1]); %maximizes figure window
figure(figNum);
figNum = figNum+1;

plot(T_Zenith_SNR,P_Zenith_SNR,'b');
xlim([0 plot_Freq_Lim]);
hold on;
% plot(T_SNR_08apr08_10_or,P_SNR_08apr08_10_or,'r');
% hold on;
plot(T_Zenith_SNR,gf_Zenith_SNR,'LineWidth',1.5,'color',[0.5 0.5 0.5]);
xlim([0 plot_Freq_Lim]);  
xlabel('Periods (min)');
ylabel('Normalized Power');
title(['MISA SNR BP at ',num2str(latMillstone),' deg lat/',num2str(longMillstone),'deg long on Day ',num2str(floor(mean(time))),', 2013 - SCARGLED']);
set(gca,'fontweight',font_Weight, 'FontSize',font_Size);
set(gca,'xtick',Xaxisvar_SCARGLE);

end %*********END OF FLG_ISRdata_MISA_Scargle*********


%% IT'S SCARGLIN TIME (sTEC filtered BP)
if( FLG_gatherDataAtPoint_Filter_BP_Scargle == 1 )

tic
for(gg = 1:size(latPoints,2) ) %gather data for multiple points
    Ross_Lombscargle_optimized([(timeUnique_combined_noNaN'-floor(mean(time))).*24 ,sTEC_combined_bp(gg,:)']); 
    fid = fopen('power_freq.txt', 'r');
        f_nm = textscan(fid, '%f %f %f %f %f', 'headerlines', 15);
        P_Zenith_SNR=f_nm{:,2};
        F_Zenith_SNR=f_nm{:,1};
        T_Zenith_SNR=f_nm{:,5}*60; % in minutes
        gs_Zenith_SNR=f_nm{:,4};
        max_Zenith_SNR=T_Zenith_SNR(P_Zenith_SNR==max(P_Zenith_SNR));
        K=length(P_Zenith_SNR);
        gf_Zenith_SNR=f_nm{:,3};
    fclose('all');

    %prep figure #
%     figure('units','normalized','outerposition',[0 0 1 1]); %maximizes figure window
    figure(figNum);
    figNum = figNum+1;
    
    plot(T_Zenith_SNR,P_Zenith_SNR,'b');
    xlim([0 plot_Freq_Lim]);
    hold on;
    % plot(T_SNR_08apr08_10_or,P_SNR_08apr08_10_or,'r');
    % hold on;
    plot(T_Zenith_SNR,gf_Zenith_SNR,'LineWidth',1.5,'color',[0.5 0.5 0.5]);
    xlim([0 plot_Freq_Lim]);  
    xlabel('Periods (min)');
    ylabel('Normalized Power');
    title(['Delta sTEC Band-Passed at deg from ',num2str(latPoints(gg)),' deg lat/',num2str(longPoints(gg)),'deg long on Day ',num2str(floor(mean(time))),', 2013 - SCARGLED']);
    set(gca,'fontweight',font_Weight, 'FontSize',font_Size);
    set(gca,'xtick',Xaxisvar_SCARGLE);
end
toc

end %*********END OF FLG_gatherDataAtPoint_Filter_BP_Scargle*********


%% ==============================Bonus FFT attempt=========================

% for(gg = 1:size(latPoints,2) ) %gather data for multiple points
%     totalcount = size(sTEC_combined_bp(gg,:),2); % # of samples taken
%     deltatime = timeDelta*24; %sample interval in time - units hours *)
%     deltafreq = 1/(totalcount*deltatime); %frequency resolution per above sampling interval in time
%     frequencysamples = ((0:(totalcount-1))*deltafreq)';
%     % frequencysamples = linspace(0,ceil((totalcount-1)*deltafreq),totalcount)';
%     % spectrummax = frequencysamples(totalcount/2 - 1);
%     timedomain = zeros(totalcount,1);
% 
%     frequencydomain = fft( sTEC_combined_bp(gg,:) );
%     powerspectrum = sqrt(abs(frequencydomain.*conj(frequencydomain)))/totalcount; %properly normalized now
% 
% 
%     %prep figure #
% %     figure('units','normalized','outerposition',[0 0 1 1]); %maximizes figure window
%     figure(figNum);
%     figNum = figNum+1;
%     plot(1./frequencysamples*60,powerspectrum)
%     title(['FFT of Delta sTEC Band-Passed at deg from ',num2str(latPoints(gg)),' deg lat/',num2str(longPoints(gg)),'deg long on Day ',num2str(floor(mean(time))),', 2013'],'FontSize',14);
%     xlabel('Period (min)','FontSize',14);
%     ylabel('Power','FontSize',14);
%     axis([0, plot_Freq_Lim, min(powerspectrum), max(powerspectrum)])
%     axis 'auto y'
% end


%% Compare Strong MSTID times (chosen time spanS), Smooth sTEC TOO!
if( FLG_gatherDataAtPoint_Filter_TimeMatch_n_Smooth == 1 )

gg = 1:length(latPoints);
gg = gg(min(abs(latPoints - latMillstone)) == abs(latPoints - latMillstone)); %get which lat is closer to Millstone

%THIS BIT IS FOR FINAL PLOT OF STEC vs EXPANDED ZENITH vs EXPANDED MISA
tmin = find( min( abs(min(time_cutout_range) - (timeUnique_combined_noNaN - dateRange_zeroHr(2)).*24 ) ) == abs(min(time_cutout_range) - (timeUnique_combined_noNaN- dateRange_zeroHr(2)).*24 ) ); %find min timeUnique that fits time limit
tmax = find( min( abs(max(time_cutout_range) - (timeUnique_combined_noNaN -  dateRange_zeroHr(2)).*24 ) ) == abs(max(time_cutout_range) - (timeUnique_combined_noNaN- dateRange_zeroHr(2)).*24 ) ); %find max timeUnique that fits time limit
timeUnique_limTIME = timeUnique_combined_noNaN(tmin:tmax); %cut out time that matches time limits (Zenith or MISA have been expanded to this time cadence)
sTEC_combined_limTIME = sTEC_combined(gg_ZenithSite,tmin:tmax); %cut out time that matches TEC time
sTEC_combined_limTIME_smoothed = smooth(sTEC_combined_limTIME)'; %smooth it out (hopefully match up better with below!)

tmin = find( min( abs(min(time_cutout_range) - (timeUnique - dateRange_zeroHr(2)).*24 ) ) == abs(min(time_cutout_range) - (timeUnique - dateRange_zeroHr(2)).*24 ) ); %find min timeUnique that fits time limit
tmax = find( min( abs(max(time_cutout_range) - (timeUnique -  dateRange_zeroHr(2)).*24 ) ) == abs(max(time_cutout_range) - (timeUnique- dateRange_zeroHr(2)).*24 ) ); %find max timeUnique that fits time limit
timeUnique_Interped_limTIME = timeUnique(tmin:tmax); %cut out time that matches time limits (Zenith or MISA have been expanded to this time cadence)
sTEC_combinedInterped_limTIME = sTEC_combined(gg_ZenithSite,tmin:tmax); %cut out time that matches TEC time
sTEC_combinedInterped_limTIME_smoothed = smooth(sTEC_combinedInterped_limTIME)'; %smooth it out (hopefully match up better with below!)

%THIS BIT IS FOR FINAL PLOT OF STEC vs EXPANDED ZENITH vs EXPANDED MISA
tmin = find( min( abs(min(time_cutout_range) - (timeUnique_limISR_Zenith - dateRange_zeroHr(2)).*24 ) ) == abs(min(time_cutout_range) - (timeUnique_limISR_Zenith- dateRange_zeroHr(2)).*24 ) ); %find min timeUnique that fits time limit
tmax = find( min( abs(max(time_cutout_range) - (timeUnique_limISR_Zenith -  dateRange_zeroHr(2)).*24 ) ) == abs(max(time_cutout_range) - (timeUnique_limISR_Zenith- dateRange_zeroHr(2)).*24 ) ); %find max timeUnique that fits time limit
timeUnique_limTIME_Zenith = timeUnique_limISR_Zenith(tmin:tmax); %cut out time that matches time limits (Zenith or MISA have been expanded to this time cadence)
Zenith_SNR_threeHun_AVGD_expanded_limTIME = Zenith_SNR_threeHun_AVGD_expanded(tmin:tmax); %cut out time that matches ISR time
%THIS BIT IS FOR FINAL PLOT OF STEC vs EXPANDED ZENITH vs EXPANDED MISA
tmin = find( min( abs(min(time_cutout_range) - (timeUnique_limISR_MISA - dateRange_zeroHr(2)).*24 ) ) == abs(min(time_cutout_range) - (timeUnique_limISR_MISA- dateRange_zeroHr(2)).*24 ) ); %find min timeUnique that fits time limit
tmax = find( min( abs(max(time_cutout_range) - (timeUnique_limISR_MISA -  dateRange_zeroHr(2)).*24 ) ) == abs(max(time_cutout_range) - (timeUnique_limISR_MISA- dateRange_zeroHr(2)).*24 ) ); %find max timeUnique that fits time limit
timeUnique_limTIME_MISA = timeUnique_limISR_MISA(tmin:tmax); %cut out time that matches time limits (Zenith or MISA have been expanded to this time cadence)
MISA_SNR_threeHun_AVGD_expanded_limTIME = MISA_SNR_threeHun_AVGD_expanded(tmin:tmax); %cut out time that matches ISR time

%THIS BIT IS FOR ZENITH SNR PLOTS (w/ height)
tmin = find( min( abs(min(time_cutout_range) - (Zenith_time - dateRange_zeroHr(2)).*24 ) ) == abs(min(time_cutout_range) - (Zenith_time- dateRange_zeroHr(2)).*24 ) ); %find min timeUnique that fits time limit
tmax = find( min( abs(max(time_cutout_range) - (Zenith_time -  dateRange_zeroHr(2)).*24 ) ) == abs(max(time_cutout_range) - (Zenith_time- dateRange_zeroHr(2)).*24 ) ); %find max timeUnique that fits time limit
Zenith_time_limTIME = Zenith_time(tmin:tmax); %cut out time that matches time limits (Zenith or MISA have been expanded to this time cadence)
Zenith_SNR_bp_limTIME = Zenith_SNR_bp(tmin:tmax,:); %cut out time that matches ISR time
Zenith_SNR_threeHun_AVGD_limTIME = Zenith_SNR_threeHun_AVGD(tmin:tmax,:); %cut out time that matches ISR time
Zenith_vel_limTIME = Zenith_vel(tmin:tmax,:); %cut out time that matches ISR time


%THIS BIT IS FOR MISA SNR PLOTS (w/ height)
tmin = find( min( abs(min(time_cutout_range) - (MISA_time - dateRange_zeroHr(2)).*24 ) ) == abs(min(time_cutout_range) - (MISA_time- dateRange_zeroHr(2)).*24 ) ); %find min timeUnique that fits time limit
tmax = find( min( abs(max(time_cutout_range) - (MISA_time -  dateRange_zeroHr(2)).*24 ) ) == abs(max(time_cutout_range) - (MISA_time- dateRange_zeroHr(2)).*24 ) ); %find max timeUnique that fits time limit
MISA_time_limTIME = MISA_time(tmin:tmax); %cut out time that matches time limits (Zenith or MISA have been expanded to this time cadence)
MISA_SNR_bp_limTIME = MISA_SNR_bp(tmin:tmax,:); %cut out time that matches ISR time
MISA_SNR_threeHun_AVGD_limTIME = MISA_SNR_threeHun_AVGD(tmin:tmax,:); %cut out time that matches ISR time
MISA_vel_limTIME = MISA_vel(tmin:tmax,:); %cut out time that matches ISR time

%THIS BIT IS FOR ZENITH SNR PLOTS (w/ height)

%ZENITH SNR PLOT
% figure('units','normalized','outerposition',[0 0 1 1]); %maximizes figure window
figure(figNum);
figNum = figNum+1;
subplot(3,1,1)
h = pcolor((Zenith_time_limTIME-floor(mean(Zenith_time))).*24,Zenith_height,Zenith_SNR_bp_limTIME');
set(h, 'EdgeColor', 'none'); %remove grid lines (which make it look hilarious)
shading interp;
lighting phong;
colorbar;
caxis([-0.1 0.1]);
colormap('gray');
ylim([90 700]);
% ylim([90 500]);
% xlim([0 70]);
% xlabel('Time in UT - 0 Hr at 0 UT May 7th','fontweight','bold','FontSize',12);
ylabel('Height (km)');
title(['Hours ',num2str(min(time_cutout_range)),' to ',num2str(max(time_cutout_range)),' (0 hr May 7) ',...
    'That Accentuate MSTIDs - Zenith SNR HP']);
set(gca,'fontweight',font_Weight, 'FontSize',font_Size);
Xaxisvar = floor(min((Zenith_time_limTIME-floor(mean(Zenith_time))).*24)):1:ceil(max((Zenith_time_limTIME-floor(mean(Zenith_time))).*24)); %1 step
% Xaxisvar = [-12, -8, -4, 0, 4, 8, 12, 16, 20, 24 28, 32, 36, 40, 44, 48]; %shows off important magnitude points
% Xaxisvar = [-12,-8,-4,0, 4, 8, 12, 16, 20, 24 28, 32, 36, 40, 44, 48, 52, 56, 60, 64, 68, 72]; %shows off important magnitude points
set(gca,'xtick',Xaxisvar);

%MISA SNR PLOT
subplot(3,1,2)
h = pcolor((MISA_time_limTIME-floor(mean(MISA_time))).*24,MISA_height,MISA_SNR_bp_limTIME');
set(h, 'EdgeColor', 'none'); %remove grid lines (which make it look hilarious)
shading interp;
lighting phong;
colorbar;
caxis([-0.1 0.1]);
colormap('gray');
ylim([90 700]);
% ylim([90 500]);
% xlim([0 70]);
xlabel('Time in UT - 0 Hr at 0 UT May 7th');
ylabel('Height (km)');
title('MISA SNR HP');
set(gca,'fontweight',font_Weight, 'FontSize',font_Size);
% Xaxisvar = [-12, -8, -4, 0, 4, 8, 12, 16, 20, 24 28, 32, 36, 40, 44, 48]; %shows off important magnitude points
% Xaxisvar = [-12,-8,-4,0, 4, 8, 12, 16, 20, 24 28, 32, 36, 40, 44, 48, 52, 56, 60, 64, 68, 72]; %shows off important magnitude points
set(gca,'xtick',Xaxisvar);

%STEC/ZENITH AVG'D/MISA AVG'D
subplot(3,1,3)
h = plot( (timeUnique_limTIME-floor(mean(timeUnique))).*24 ,sTEC_combined_limTIME_smoothed ); %plot all at once on this
hold on;
set(h(1),'Color','g')
h = plot( (timeUnique_limTIME_Zenith-floor(mean(timeUnique))).*24 ,Zenith_SNR_threeHun_AVGD_expanded_limTIME ); %plot all at once on this
set(h(1),'Color','b','LineStyle','--')
h = plot( (timeUnique_limTIME_MISA-floor(mean(timeUnique))).*24 ,MISA_SNR_threeHun_AVGD_expanded_limTIME ); %plot all at once on this
set(h(1),'Color','r','Marker','.')
xlabel('Time (hr)'); %old: ,'fontweight','bold','FontSize',12 for all
ylabel('TECU or SNR'); %sped up significantly when merged above
title(['Delta sTEC BP AND Zenith/MISA SNR BP AVG''d w/ Interp. +/-',num2str(Zenith_avgRange),' km around ',num2str(Zenith_threeHun),' km at ',num2str(latMillstone),' deg lat/',num2str(longMillstone),'deg long']);
axis([-inf,inf,-inf,inf])
legend('sTEC band-passed','Zenith SNR band-passed AVG''d','MISA SNR band-passed AVG''d');
set(gca,'fontweight',font_Weight, 'FontSize',font_Size);

R_Zenith = corrcoef(sTEC_combinedInterped_limTIME_smoothed,Zenith_SNR_threeHun_AVGD_expanded_limTIME);
fprintf('Correlation coeff between Zenith ISR and sTEC Smoothed over Millstone Cutout %d to %d hr: %f\n',min(time_cutout_range),max(time_cutout_range),R_Zenith(2));

R_MISA = corrcoef(sTEC_combinedInterped_limTIME_smoothed,MISA_SNR_threeHun_AVGD_expanded_limTIME);
fprintf('Correlation coeff between MISA ISR and sTEC Smoothed over Millstone Cutout %d to %d hr: %f\n\n',min(time_cutout_range),max(time_cutout_range),R_MISA(2));


% SCARGLE IT
Ross_Lombscargle_optimized([(timeUnique_limTIME'-floor(mean(timeUnique))).*24 ,sTEC_combined_limTIME_smoothed'  ]); 
fid = fopen('power_freq.txt', 'r');
    f_nm = textscan(fid, '%f %f %f %f %f', 'headerlines', 15);
    P_Zenith_SNR=f_nm{:,2};
    F_Zenith_SNR=f_nm{:,1};
    T_Zenith_SNR=f_nm{:,5}*60; % in minutes
    gs_Zenith_SNR=f_nm{:,4};
    max_Zenith_SNR=T_Zenith_SNR(P_Zenith_SNR==max(P_Zenith_SNR));
    K=length(P_Zenith_SNR);
    gf_Zenith_SNR=f_nm{:,3};
fclose('all');

%prep figure #
% figure('units','normalized','outerposition',[0 0 1 1]); %maximizes figure window
figure(figNum);
figNum = figNum+1;

plot(T_Zenith_SNR,P_Zenith_SNR,'b');
xlim([0 plot_Freq_Lim]);
hold on;
% plot(T_SNR_08apr08_10_or,P_SNR_08apr08_10_or,'r');
% hold on;
plot(T_Zenith_SNR,gf_Zenith_SNR,'LineWidth',1.5,'color',[0.5 0.5 0.5]);
xlim([0 plot_Freq_Lim]);  
xlabel('Periods (min)');
ylabel('Normalized Power');
title(['Spectrum - sTEC Time Limited Band-Passed Smoothed at ',num2str(latMillstone),' deg lat/',num2str(longMillstone),'deg long on Day ',num2str(floor(mean(time))),', 2013']);
set(gca,'fontweight',font_Weight, 'FontSize',font_Size);
set(gca,'xtick',Xaxisvar_SCARGLE);

end %*********END OF FLG_gatherDataAtPoint_Filter_TimeMatch_n_Smooth*********


%% Compare Strong MSTID times with Interp'd to 5 min on STEC to Match ISR
if( FLG_gatherDataAtPoint_Filter_TimeMatch_n_Smooth_ISRCorr == 1 )

sTEC_5min_Zenith = zeros(length(Zenith_time),1);
sTEC_5minInterped_Zenith  = zeros(length(Zenith_time),1);
Zenith_time_delta = Zenith_time(162)-Zenith_time(161); %days, delta of time between readings

for(i = 1:length(Zenith_time) ) %no nan time table is different
    
    if(i == 1)
        jold = find( min(abs((Zenith_time(i)-Zenith_time_delta) - timeUnique_combined_noNaN)) == abs((Zenith_time(i)-Zenith_time_delta) - timeUnique_combined_noNaN) ); %get the matching time 5 min prev. (since it's an integration up to the time stamp given)
        jcurrent = find( min(abs((Zenith_time(i)) - timeUnique_combined_noNaN)) == abs((Zenith_time(i)) - timeUnique_combined_noNaN) ); %get the matching time
        
        sTEC_5min_Zenith(i) = mean(sTEC_combined(gg_ZenithSite,jold:jcurrent)); % average over the ISR time period
    else
        jold = jcurrent; %set the old time
        jcurrent = find( min(abs((Zenith_time(i)) - timeUnique_combined_noNaN)) == abs((Zenith_time(i)) - timeUnique_combined_noNaN) ); %get the matching time
        
        sTEC_5min_Zenith(i) = mean(sTEC_combined(gg_ZenithSite,jold+1:jcurrent)); % average over the ISR time period (increment jold by 1 so not using same data)
    end
    
end

for(i = 1:length(Zenith_time) )
    
    if(i == 1)
        jold = find( min(abs((Zenith_time(i)-Zenith_time_delta) - timeUnique)) == abs((Zenith_time(i)-Zenith_time_delta) - timeUnique) ); %get the matching time 5 min prev. (since it's an integration up to the time stamp given)
        jcurrent = find( min(abs((Zenith_time(i)) - timeUnique)) == abs((Zenith_time(i)) - timeUnique) ); %get the matching time
        
        sTEC_5minInterped_Zenith(i) = mean(sTEC_combinedInterped(gg_ZenithSite,jold:jcurrent)); % average over the ISR time period
    else
        jold = jcurrent; %set the old time
        jcurrent = find( min(abs((Zenith_time(i)) - timeUnique)) == abs((Zenith_time(i)) - timeUnique) ); %get the matching time
        
        sTEC_5minInterped_Zenith(i) = mean(sTEC_combinedInterped(gg_ZenithSite,jold:jcurrent)); % average over the ISR time period
    end
    
end
Zenith_time_NoNaN = Zenith_time; %copy zenith time
Zenith_time_NoNaN( isnan(sTEC_5min_Zenith) ) = []; %delete the NaNs
sTEC_5min_Zenith( isnan(sTEC_5min_Zenith) ) = []; %delete the NaNs

sTEC_5min_MISA = zeros(length(MISA_time),1);
sTEC_5minInterped_MISA  = zeros(length(Zenith_time),1);
MISA_time_delta = MISA_time(162)-MISA_time(161); %days, delta of time between readings

for(i = 1:length(MISA_time) ) %no nan time table is different (MISA DOES NOT WORK FOR THIS -included in case I make it work)
    
    if(i == 1) %special case to start - prev 5 min
        jold = find( min(abs((MISA_time(i)-MISA_time_delta) - timeUnique_combined_noNaN)) == abs((MISA_time(i)-MISA_time_delta) - timeUnique_combined_noNaN) ); %get the matching time 5 min prev.
        jcurrent = find( min(abs((MISA_time(i)) - timeUnique_combined_noNaN)) == abs((MISA_time(i)) - timeUnique_combined_noNaN) ); %get the matching time
        
        sTEC_5min_MISA(i) = mean(sTEC_combined(gg_MISASite,jold:jcurrent)); % average over the ISR time period
    else
        jold = jcurrent; %set the old time
        jcurrent = find( min(abs((MISA_time(i)) - timeUnique_combined_noNaN)) == abs((MISA_time(i)) - timeUnique_combined_noNaN) ); %get the matching time
        
        sTEC_5min_MISA(i) = mean(sTEC_combined(gg_MISASite,jold+1:jcurrent)); % average over the ISR time period (increment jold by 1 so not using same data)
    end
    
end

for(i = 1:length(MISA_time) )
    
    if(i == 1) %special case to start - prev 5 min
        jold = find( min(abs((MISA_time(i)-MISA_time_delta) - timeUnique)) == abs((MISA_time(i)-MISA_time_delta) - timeUnique) ); %get the matching time 5 min prev.
        jcurrent = find( min(abs((MISA_time(i)) - timeUnique)) == abs((MISA_time(i)) - timeUnique) ); %get the matching time
        
        sTEC_5minInterped_MISA(i) = mean(sTEC_combinedInterped(gg_MISASite,jold:jcurrent)); % average over the ISR time period
    else
        jold = jcurrent; %set the old time
        jcurrent = find( min(abs((MISA_time(i)) - timeUnique)) == abs((MISA_time(i)) - timeUnique) ); %get the matching time
        
        sTEC_5minInterped_MISA(i) = mean(sTEC_combinedInterped(gg_MISASite,jold:jcurrent)); % average over the ISR time period
    end
    
end
MISA_time_NoNaN = MISA_time; %copy MISA time
MISA_time_NoNaN( isnan(sTEC_5min_MISA) ) = []; %delete the NaNs
sTEC_5min_MISA( isnan(sTEC_5min_MISA) ) = []; %delete the NaNs

tmin_Zenith = find(Zenith_time_limTIME(1) == Zenith_time); %get time indicies to limit the Zenith or MISA ISR data
tmax_Zenith = find(Zenith_time_limTIME(end) == Zenith_time);
tmin_MISA = find(MISA_time_limTIME(1) == MISA_time);
tmax_MISA = find(MISA_time_limTIME(end) == MISA_time);

R_Zenith = corrcoef(sTEC_5minInterped_Zenith(tmin_Zenith:tmax_Zenith),Zenith_SNR_threeHun_AVGD(tmin_Zenith:tmax_Zenith));
fprintf('Correlation coeff between Zenith ISR and sTEC 5 Min Averaged over Millstone Cutout %d to %d hr: %f\n',min(time_cutout_range),max(time_cutout_range),R_Zenith(2));

R_MISA = corrcoef(sTEC_5minInterped_MISA(tmin_MISA:tmax_MISA),MISA_SNR_threeHun_AVGD(tmin_MISA:tmax_MISA));
fprintf('Correlation coeff between MISA ISR and sTEC 5 Min Averaged over Millstone Cutout %d to %d hr: %f\n\n',min(time_cutout_range),max(time_cutout_range),R_MISA(2));

end %*********END OF FLG_gatherDataAtPoint_Filter_TimeMatch_n_Smooth_ISRCorr*********


%% IT'S SCARGLIN TIME LONG AVG'D (so time vs latitude is plot) TIME LIMITED AND BAND PASSED
if( FLG_AVG_LONG_TimeMatch_n_BP_Scargle == 1 )

jk = find( min(abs( pplatChunksPlot - latMillstone )) == abs( pplatChunksPlot - latMillstone ) ); %get the location closest
jm = isnan(sTECChunked_longavg(:,jk)); %get the NaN
dataTime = (timeUnique_limISR_MISA(~jm)' -floor(mean(timeUnique))).*24; %time data time
sTECChunked_longavg_atMillstone = interp1( dataTime , sTECChunked_longavg(~jm,jk) , (timeUnique_limISR_MISA' -floor(mean(timeUnique))).*24 ,'linear','extrap'); %interpolate extra data in

tmin = find( min( abs(min(time_cutout_range) - (timeUnique_limISR_MISA - floor(mean(timeUnique))).*24 ) ) == abs(min(time_cutout_range) - (timeUnique_limISR_MISA- floor(mean(timeUnique))).*24 ) ); %find min timeUnique that fits time limit
tmax = find( min( abs(max(time_cutout_range) - (timeUnique_limISR_MISA -  floor(mean(timeUnique))).*24 ) ) == abs(max(time_cutout_range) - (timeUnique_limISR_MISA- floor(mean(timeUnique))).*24 ) ); %find max timeUnique that fits time limit
timeUnique_limTIME_MISA = timeUnique_limISR_MISA(tmin:tmax); %cut out time that matches time limits (Zenith or MISA have been expanded to this time cadence)
sTECChunked_longavg_atMillstone_limTIME = sTECChunked_longavg_atMillstone(tmin:tmax); %time limited

bs=(1/2);   % highpass cuttoff frequency (not sure what to make of it)
%I think it is related to lf above (which is 1/2 )

n=42; % order of the Hamming window used in the custom function
% c = 3.32*pi; %Hamming constant
% M = n/2;
% bandwidth = c/M;
%The above was included to investigate the bandwidth - 1/2. Related to
%bs/lf?

fp=bs; % stop band stoppin freq (or pass band passing freq, depending on how you look at it)
timeDelta = timeUnique_limTIME_MISA(50)-timeUnique_limTIME_MISA(49); %hrs already
f= 1/(timeDelta); % the sampling frequency, based off of the time delta calc'd

wp=2*fp/f; % Normalizing the frequencies (Matlab is 0 to 1)

%Calculation of filter coefficients
% [b,a]=fir1(n,wp,'high'); %This takes the order and lower cut off
%Uses the default hamming window
% frequency ORIG
W = hann(n+1);
[b,a] = fir1(n,wp,'high',W);
%Applys Hanning Window for shiggles and giggles

sTECChunked_longavg_atMillstone_limTIME_bp = filtfilt(b,a,sTECChunked_longavg_atMillstone_limTIME); %Appies the filter


Ross_Lombscargle_optimized([(timeUnique_limTIME_MISA' -floor(mean(timeUnique))).*24 , sTECChunked_longavg_atMillstone_limTIME_bp ]); 
fid = fopen('power_freq.txt', 'r');
    f_nm = textscan(fid, '%f %f %f %f %f', 'headerlines', 15);
    P_Zenith_SNR=f_nm{:,2};
    F_Zenith_SNR=f_nm{:,1};
    T_Zenith_SNR=f_nm{:,5}*60; % in minutes
    gs_Zenith_SNR=f_nm{:,4};
    max_Zenith_SNR=T_Zenith_SNR(P_Zenith_SNR==max(P_Zenith_SNR));
    K=length(P_Zenith_SNR);
    gf_Zenith_SNR=f_nm{:,3};
fclose('all');

%prep figure #
% figure('units','normalized','outerposition',[0 0 1 1]); %maximizes figure window
figure(figNum);
figNum = figNum+1;

plot(T_Zenith_SNR,P_Zenith_SNR,'b');
xlim([0 plot_Freq_Lim]);
hold on;
% plot(T_SNR_08apr08_10_or,P_SNR_08apr08_10_or,'r');
% hold on;
plot(T_Zenith_SNR,gf_Zenith_SNR,'LineWidth',1.5,'color',[0.5 0.5 0.5]);
xlim([0 plot_Freq_Lim]);  
xlabel('Periods (min)');
ylabel('Normalized Power');
title(['Longitude Avg''d Band-Passed Time Limited to ',num2str(min(time_cutout_range)),' to ',num2str(max(time_cutout_range)),' hrs at Millstone Latitude of ',num2str(latMillstone),' deg lat with longitude AVD''d +/-',num2str(pplongAvgDeg),'deg on Day ',num2str(floor(mean(time))),', 2013 - SCARGLED']);
set(gca,'fontweight',font_Weight, 'FontSize',font_Size);
set(gca,'xtick',Xaxisvar_SCARGLE);

end %*********END OF FLG_AVG_LONG_TimeMatch_n_BP_Scargle*********


%% LONG AVG'D AND ISR VELOCITY AND PLOTTED ONTOP OF EACH OTHER
if( FLG_AVG_LONG_TimeMatch_n_BP_Velocity == 1 )
        
% MISA at 300km VS LONG AVG at same Lat plot
% figure('units','normalized','outerposition',[0 0 1 1]); %maximizes figure window
figure(figNum);
figNum = figNum+1;
h = plot( (timeUnique_limTIME_MISA' -floor(mean(timeUnique))).*24 ,sTECChunked_longavg_atMillstone_limTIME_bp ); %plot all at once on this
hold on;
set(h(1),'Color','g')
h = plot( (MISA_time_limTIME-floor(mean(MISA_time))).*24 ,MISA_SNR_threeHun_AVGD_limTIME ); %plot all at once on this
hold on;
set(h(1),'Color','b','LineStyle','--')
h = plot( (timeUnique_limTIME_MISA' -floor(mean(timeUnique))).*24 ,sTEC_combined_limTIME ); %plot all at once on this
hold on;
set(h(1),'Color','r','LineStyle','-.')
xlabel('Time (hr)'); %old: ,'fontweight','bold','FontSize',12 for all
ylabel('sTEC AVG''D in Long (TECU) or MISA SNR bp AVG''d (dB)'); %sped up significantly when merged above
title(['sTEC AVG''D ',num2str(pplongAvgDeg),' deg in Long and at ',num2str(latMillstone),' deg Latitude AND MISA SNR BP AVG''d w/ Interp. +/-',num2str(Zenith_avgRange),' km around ',num2str(Zenith_threeHun),' km at ',num2str(latMillstone),' deg lat/',num2str(longMillstone),'deg long']);
axis([-inf,inf,-inf,inf])
legend('sTEC AVG''d in Long','MISA SNR band-passed AVG''d','sTEC at MISA');
set(gca,'fontweight',font_Weight, 'FontSize',font_Size);

R_AVGLong_sTEC = corrcoef(sTEC_combined_limTIME,sTECChunked_longavg_atMillstone_limTIME_bp);
fprintf('Correlation coeff between Zenith ISR and sTEC 5 Min Averaged over Millstone: %f\n',R_AVGLong_sTEC(2));


% MISA SCARGLE PLOT FOR TIME LIMITED TO THE HOURS DESCRIBED
Ross_Lombscargle_optimized([(MISA_time_limTIME-floor(mean(MISA_time))).*24 ,MISA_SNR_threeHun_AVGD_limTIME ]); 
fid = fopen('power_freq.txt', 'r');
    f_nm = textscan(fid, '%f %f %f %f %f', 'headerlines', 15);
    P_Zenith_SNR=f_nm{:,2};
    F_Zenith_SNR=f_nm{:,1};
    T_Zenith_SNR=f_nm{:,5}*60; % in minutes
    gs_Zenith_SNR=f_nm{:,4};
    max_Zenith_SNR=T_Zenith_SNR(P_Zenith_SNR==max(P_Zenith_SNR));
    K=length(P_Zenith_SNR);
    gf_Zenith_SNR=f_nm{:,3};
fclose('all');

%prep figure #
% figure('units','normalized','outerposition',[0 0 1 1]); %maximizes figure window
figure(figNum);
figNum = figNum+1;

plot(T_Zenith_SNR,P_Zenith_SNR,'b');
xlim([0 plot_Freq_Lim]);
hold on;
% plot(T_SNR_08apr08_10_or,P_SNR_08apr08_10_or,'r');
% hold on;
plot(T_Zenith_SNR,gf_Zenith_SNR,'LineWidth',1.5,'color',[0.5 0.5 0.5]);
xlim([0 plot_Freq_Lim]);  
xlabel('Periods (min)');
ylabel('Normalized Power');
title(['ISR MISA Time Limited to ',num2str(min(time_cutout_range)),' to ',num2str(max(time_cutout_range)),' hrs on Day ',num2str(floor(mean(MISA_time_limTIME))),', 2013 - SCARGLED']);
set(gca,'fontweight',font_Weight, 'FontSize',font_Size);
set(gca,'xtick',Xaxisvar_SCARGLE);


% sTEC AT MISA POINT SCARGLE PLOT FOR TIME LIMITED TO THE HOURS DESCRIBED
Ross_Lombscargle_optimized([(timeUnique_limTIME_MISA' -floor(mean(timeUnique))).*24 ,sTEC_combined_limTIME' ]); 
fid = fopen('power_freq.txt', 'r');
    f_nm = textscan(fid, '%f %f %f %f %f', 'headerlines', 15);
    P_Zenith_SNR=f_nm{:,2};
    F_Zenith_SNR=f_nm{:,1};
    T_Zenith_SNR=f_nm{:,5}*60; % in minutes
    gs_Zenith_SNR=f_nm{:,4};
    max_Zenith_SNR=T_Zenith_SNR(P_Zenith_SNR==max(P_Zenith_SNR));
    K=length(P_Zenith_SNR);
    gf_Zenith_SNR=f_nm{:,3};
fclose('all');

%prep figure #
% figure('units','normalized','outerposition',[0 0 1 1]); %maximizes figure window
figure(figNum);
figNum = figNum+1;

plot(T_Zenith_SNR,P_Zenith_SNR,'b');
xlim([0 plot_Freq_Lim]);
hold on;
% plot(T_SNR_08apr08_10_or,P_SNR_08apr08_10_or,'r');
% hold on;
plot(T_Zenith_SNR,gf_Zenith_SNR,'LineWidth',1.5,'color',[0.5 0.5 0.5]);
xlim([0 plot_Freq_Lim]);  
xlabel('Periods (min)');
ylabel('Normalized Power');
title(['Spectrum - sTEC at MISA Time Limited to ',num2str(min(time_cutout_range)),' to ',num2str(max(time_cutout_range)),' hrs on Day ',num2str(floor(mean(MISA_time_limTIME))),', 2013']);
set(gca,'fontweight',font_Weight, 'FontSize',font_Size);
set(gca,'xtick',Xaxisvar_SCARGLE);

    
end %*********END OF FLG_AVG_LONG_TimeMatch_n_BP_Velocity*********


%% Strong MSTID times w/ sTEC on ISR Time Scale, Plots n Such
if( FLG_gatherDataAtPoint_Filter_TimeMatch_n_Smooth_Scargle == 1 )

%ZENITH SNR PLOT
% figure('units','normalized','outerposition',[0 0 1 1]); %maximizes figure window
figure(figNum);
figNum = figNum+1;
subplot(3,1,1)
h = pcolor((Zenith_time_limTIME-floor(mean(Zenith_time))).*24,Zenith_height,Zenith_SNR_bp_limTIME');
set(h, 'EdgeColor', 'none'); %remove grid lines (which make it look hilarious)
shading interp;
lighting phong;
colorbar;
caxis([-0.1 0.1]);
colormap('gray');
ylim([90 700]);
% ylim([90 500]);
% xlim([0 70]);
xlabel('Time in UT - 0 Hr at 0 UT May 7th');
ylabel('Height (km) FOR ZENITH');
title(['Hours ',num2str(min(time_cutout_range)),' to ',num2str(max(time_cutout_range)),' (0 hr May 7) ',...
    'That Accentuate MSTIDs - Zenith SNR BP']);
set(gca,'fontweight',font_Weight, 'FontSize',font_Size);
% Xaxisvar = [-12, -8, -4, 0, 4, 8, 12, 16, 20, 24 28, 32, 36, 40, 44, 48]; %shows off important magnitude points
% Xaxisvar = [-12,-8,-4,0, 4, 8, 12, 16, 20, 24 28, 32, 36, 40, 44, 48, 52, 56, 60, 64, 68, 72]; %shows off important magnitude points
set(gca,'xtick',Xaxisvar);

%MISA SNR PLOT
subplot(3,1,2)
h = pcolor((MISA_time_limTIME-floor(mean(MISA_time))).*24,MISA_height,MISA_SNR_bp_limTIME');
set(h, 'EdgeColor', 'none'); %remove grid lines (which make it look hilarious)
shading interp;
lighting phong;
colorbar;
caxis([-0.1 0.1]);
colormap('gray');
ylim([90 700]);
% ylim([90 500]);
% xlim([0 70]);
xlabel('Time in UT - 0 Hr at 0 UT May 7th');
ylabel('Height (km) FOR MISA');
title('MISA SNR BP');
set(gca,'fontweight',font_Weight, 'FontSize',font_Size);
% Xaxisvar = [-12, -8, -4, 0, 4, 8, 12, 16, 20, 24 28, 32, 36, 40, 44, 48]; %shows off important magnitude points
% Xaxisvar = [-12,-8,-4,0, 4, 8, 12, 16, 20, 24 28, 32, 36, 40, 44, 48, 52, 56, 60, 64, 68, 72]; %shows off important magnitude points
set(gca,'xtick',Xaxisvar);

tmin_Zenith_NoNaN = find(Zenith_time_limTIME(1) == Zenith_time_NoNaN); %get time indicies to limit the Zenith or MISA ISR data
tmax_Zenith_NoNaN = find(Zenith_time_limTIME(end) == Zenith_time_NoNaN);
%STEC/ZENITH AVG'D/MISA AVG'D
subplot(3,1,3)
h = plot( (Zenith_time_NoNaN(tmin_Zenith_NoNaN:tmax_Zenith_NoNaN)-floor(mean(Zenith_time))).*24 ,sTEC_5min_Zenith(tmin_Zenith_NoNaN:tmax_Zenith_NoNaN) ); %plot all at once on this
hold on;
set(h(1),'Color','g')
h = plot( (Zenith_time_limTIME-floor(mean(Zenith_time))).*24 ,Zenith_SNR_threeHun_AVGD(tmin_Zenith:tmax_Zenith) ); %plot all at once on this
hold on;
set(h(1),'Color','b','LineStyle','--')
h = plot( (MISA_time_limTIME -floor(mean(MISA_time))).*24 ,MISA_SNR_threeHun_AVGD(tmin_MISA:tmax_MISA) ); %plot all at once on this
set(h(1),'Color','r','Marker','.')
xlabel('Time (hr)'); %old: ,'fontweight','bold','FontSize',12 for all
ylabel('Delta sTEC bp (TECU) or Zenith/MISA SNR bp AVG''d (dB)'); %sped up significantly when merged above
title(['Delta sTEC BP ISR match AVG''d AND Zenith/MISA SNR BP AVG''d w/ Interp. +/-',num2str(Zenith_avgRange),' km around ',num2str(Zenith_threeHun),' km at ',num2str(latMillstone),' deg lat/',num2str(longMillstone),'deg long']);
axis([-inf,inf,-inf,inf])
legend('sTEC band-passed ISR matched AVG''d','Zenith SNR band-passed AVG''d','MISA SNR band-passed AVG''d');
set(gca,'fontweight',font_Weight, 'FontSize',font_Size);

% SCARGLE IT
Ross_Lombscargle_optimized([(Zenith_time_limTIME-floor(mean(Zenith_time))).*24 ,sTEC_5min_Zenith(tmin_Zenith:tmax_Zenith)  ]); 
fid = fopen('power_freq.txt', 'r');
    f_nm = textscan(fid, '%f %f %f %f %f', 'headerlines', 15);
    P_Zenith_SNR=f_nm{:,2};
    F_Zenith_SNR=f_nm{:,1};
    T_Zenith_SNR=f_nm{:,5}*60; % in minutes
    gs_Zenith_SNR=f_nm{:,4};
    max_Zenith_SNR=T_Zenith_SNR(P_Zenith_SNR==max(P_Zenith_SNR));
    K=length(P_Zenith_SNR);
    gf_Zenith_SNR=f_nm{:,3};
fclose('all');

%prep figure #
% figure('units','normalized','outerposition',[0 0 1 1]); %maximizes figure window
figure(figNum);
figNum = figNum+1;

plot(T_Zenith_SNR,P_Zenith_SNR,'b');
xlim([0 plot_Freq_Lim]);
hold on;
% plot(T_SNR_08apr08_10_or,P_SNR_08apr08_10_or,'r');
% hold on;
plot(T_Zenith_SNR,gf_Zenith_SNR,'LineWidth',1.5,'color',[0.5 0.5 0.5]);
xlim([0 plot_Freq_Lim]);  
xlabel('Periods (min)');
ylabel('Normalized Power');
title(['Spectrum - sTEC on ISR Interval at ',num2str(latMillstone),' deg lat/',num2str(longMillstone),'deg long on Day ',num2str(floor(mean(time))),', 2013']);
set(gca,'fontweight',font_Weight, 'FontSize',font_Size);
set(gca,'xtick',Xaxisvar_SCARGLE);

end %*********END OF FLG_gatherDataAtPoint_Filter_TimeMatch_n_Smooth_Scargle*********


%% FILTER STEC 5 MIN MATCHED TO ISR DATA
if( FLG_gatherDataAtPoint_Filter_TimeMatch_n_Smooth_BP == 1 )

%These are unused, but accentuate the periods desired
%Original file used highpass_fir.m, it has been merged
% lp=1; % Lower period (in hrs)
% hp=2; % Higher period (in hrs)
% lf=(1/hp);  % Lowpass frequency corner (1/hr)
% hf=(1/lp);  % Highpass frequency corner (1/hr)

bs=(1/2);   % highpass cuttoff frequency (not sure what to make of it)
%I think it is related to lf above (which is 1/2 )

n=42; % order of the Hamming window used in the custom function
% c = 3.32*pi; %Hamming constant
% M = n/2;
% bandwidth = c/M;
%The above was included to investigate the bandwidth - 1/2. Related to
%bs/lf?

fp=bs; % stop band stoppin freq (or pass band passing freq, depending on how you look at it)
timeDelta = (Zenith_time_limTIME(50)-floor(mean(Zenith_time))).*24-(Zenith_time_limTIME(49)-floor(mean(Zenith_time))).*24; %hrs forced
f= 1/(timeDelta); % the sampling frequency, based off of the time delta calc'd

wp=2*fp/f; % Normalizing the frequencies (Matlab is 0 to 1)

%Calculation of filter coefficients
% [b,a]=fir1(n,wp,'high'); %This takes the order and lower cut off
%Uses the default hamming window
% frequency ORIG
W = hann(n+1);
[b,a] = fir1(n,wp,'high',W);
%Applys Hanning Window for shiggles and giggles

sTEC_5min_Zenith_bp = filtfilt(b,a,sTEC_5min_Zenith); %Appies the filter
sTEC_5minInterped_Zenith_bp = filtfilt(b,a,sTEC_5minInterped_Zenith); %Appies the filter

sTEC_5min_MISA_bp = filtfilt(b,a,sTEC_5min_MISA); %Appies the filter
sTEC_5minInterped_MISA_bp = filtfilt(b,a,sTEC_5minInterped_MISA); %Appies the filter

end %*********END OF FLG_gatherDataAtPoint_Filter_TimeMatch_n_Smooth_BP*********


%% Strong MSTID times w/ sTEC on ISR Time Scale, Plots n Such, BANDPASSED STEC
if( FLG_gatherDataAtPoint_Filter_TimeMatch_n_Smooth_BP_PlotsScargle == 1 )

time_cutout_range_Override = time_cutout_range; %save
time_cutout_range = [22,33]; %stuff for great comparison
tmin_Zenith_limTIME = find(min(abs(min(time_cutout_range) - (Zenith_time-dateRange_zeroHr(2))*24)) == abs(min(time_cutout_range) - (Zenith_time-dateRange_zeroHr(2))*24) ); %get time indicies to limit the Zenith or MISA ISR data
tmax_Zenith_limTIME = find(min(abs(max(time_cutout_range) - (Zenith_time-dateRange_zeroHr(2))*24)) == abs(max(time_cutout_range) - (Zenith_time-dateRange_zeroHr(2))*24) ); %get time indicies to limit the Zenith or MISA ISR data
tmin_MISA_limTIME = find(min(abs(min(time_cutout_range) - (MISA_time-dateRange_zeroHr(2))*24)) == abs(min(time_cutout_range) - (MISA_time-dateRange_zeroHr(2))*24) ); %get time indicies to limit the Zenith or MISA ISR data
tmax_MISA_limTIME = find(min(abs(max(time_cutout_range) - (MISA_time-dateRange_zeroHr(2))*24)) == abs(max(time_cutout_range) - (MISA_time-dateRange_zeroHr(2))*24) ); %get time indicies to limit the Zenith or MISA ISR data

R_Zenith = corrcoef(sTEC_5minInterped_Zenith_bp(tmin_Zenith_limTIME:tmax_Zenith_limTIME),Zenith_SNR_threeHun_AVGD(tmin_Zenith_limTIME:tmax_Zenith_limTIME));
fprintf('Correlation coeff between Zenith ISR and sTEC 5 Min Averaged Band-Passed over Millstone Cutout %d to %d hr (hard coded time): %f\n',min(time_cutout_range),max(time_cutout_range),R_Zenith(2));

R_MISA = corrcoef(sTEC_5minInterped_MISA_bp(tmin_MISA_limTIME:tmax_MISA_limTIME),MISA_SNR_threeHun_AVGD(tmin_MISA_limTIME:tmax_MISA_limTIME));
fprintf('Correlation coeff between MISA ISR and sTEC 5 Min Averaged Band-Passed over Millstone Cutout %d to %d hr (hard coded time): %f\n\n',min(time_cutout_range),max(time_cutout_range),R_MISA(2));
time_cutout_range = time_cutout_range_Override; %revert

R_Zenith = corrcoef(sTEC_5minInterped_Zenith_bp(tmin_Zenith:tmax_Zenith),Zenith_SNR_threeHun_AVGD(tmin_Zenith:tmax_Zenith));
fprintf('Correlation coeff between Zenith ISR and sTEC 5 Min Averaged Band-Passed over Millstone Cutout %d to %d hr: %f\n',min(time_cutout_range),max(time_cutout_range),R_Zenith(2));

R_MISA = corrcoef(sTEC_5minInterped_MISA_bp(tmin_MISA:tmax_MISA),MISA_SNR_threeHun_AVGD(tmin_MISA:tmax_MISA));
fprintf('Correlation coeff between MISA ISR and sTEC 5 Min Averaged Band-Passed over Millstone Cutout %d to %d hr: %f\n\n',min(time_cutout_range),max(time_cutout_range),R_MISA(2));

%FFT COMPARISON
% totalcount = length(sTEC_5minInterped_Zenith_bp(tmin_Zenith:tmax_Zenith)); % # of samples taken
% deltatime = Zenith_time_delta*24; %sample interval in time - units hours *)
% deltafreq = 1/(totalcount*deltatime); %frequency resolution per above sampling interval in time
% frequencysamples = ((0:(totalcount-1))*deltafreq)';
% % frequencysamples = linspace(0,ceil((totalcount-1)*deltafreq),totalcount)';
% % spectrummax = frequencysamples(totalcount/2 - 1);
% timedomain = zeros(totalcount,1);
% 
% frequencydomain = fft( sTEC_5minInterped_Zenith_bp(tmin_Zenith:tmax_Zenith) );
% powerspectrum = sqrt(abs(frequencydomain.*conj(frequencydomain)))/totalcount; %properly normalized now
% 
% %prep figure #
% %     figure('units','normalized','outerposition',[0 0 1 1]); %maximizes figure window
% figure(figNum);
% figNum = figNum+1;
% subplot(2,1,1);
% plot(1./frequencysamples*60,powerspectrum)
% title(['FFT of sTEC - Hours ',num2str(min(time_cutout_range)),' to ',num2str(max(time_cutout_range)),' ',...
%     'That Accentuate MSTIDs']);
% axis([0, plot_Freq_Lim, min(powerspectrum), max(powerspectrum)])
% % axis 'auto y'
% set(gca,'fontweight',font_Weight, 'FontSize',font_Size);
% 
% totalcount = length(Zenith_SNR_threeHun_AVGD(tmin_Zenith:tmax_Zenith)); % # of samples taken
% deltatime = Zenith_time_delta*24; %sample interval in time - units hours *)
% deltafreq = 1/(totalcount*deltatime); %frequency resolution per above sampling interval in time
% frequencysamples = ((0:(totalcount-1))*deltafreq)';
% % frequencysamples = linspace(0,ceil((totalcount-1)*deltafreq),totalcount)';
% % spectrummax = frequencysamples(totalcount/2 - 1);
% timedomain = zeros(totalcount,1);
% 
% frequencydomain = fft( Zenith_SNR_threeHun_AVGD(tmin_Zenith:tmax_Zenith) );
% powerspectrum = sqrt(abs(frequencydomain.*conj(frequencydomain)))/totalcount; %properly normalized now
% 
% subplot(2,1,2);
% plot(1./frequencysamples*60,powerspectrum)
% title(['FFT of Zenith HP SNR - Hours ',num2str(min(time_cutout_range)),' to ',num2str(max(time_cutout_range)),' ',...
%     'That Accentuate MSTIDs']);
% axis([0, plot_Freq_Lim, min(powerspectrum), max(powerspectrum)])
% % axis 'auto y'
% set(gca,'fontweight',font_Weight, 'FontSize',font_Size);

R_Zenith = corrcoef(sTEC_5minInterped_Zenith_bp,Zenith_SNR_threeHun_AVGD);
fprintf('Correlation coeff between Zenith ISR and sTEC 5 Min Averaged Band-Passed over Millstone full time: %f\n',R_Zenith(2));

R_MISA = corrcoef(sTEC_5minInterped_MISA_bp,MISA_SNR_threeHun_AVGD);
fprintf('Correlation coeff between MISA ISR and sTEC 5 Min Averaged Band-Passed over Millstone Cutout full time: %f\n\n',R_MISA(2));

%ZENITH SNR PLOT 
% figure('units','normalized','outerposition',[0 0 1 1]); %maximizes figure window
figure(figNum);
figNum = figNum+1;
%STEC/ZENITH AVG'D/MISA AVG'D
subplot(3,1,1)
h = plot( (Zenith_time_limTIME-floor(mean(Zenith_time))).*24 ,Zenith_SNR_threeHun_AVGD(tmin_Zenith:tmax_Zenith) ); %plot all at once on this
hold on;
set(h(1),'Color','b','LineStyle','--','LineWidth',1.25)
h = plot( (MISA_time_limTIME -floor(mean(MISA_time))).*24 ,MISA_SNR_threeHun_AVGD(tmin_MISA:tmax_MISA) ); %plot all at once on this
set(h(1),'Color','g','LineWidth',1.25)
hold on;
h = plot( (Zenith_time_NoNaN(tmin_Zenith_NoNaN:tmax_Zenith_NoNaN)-floor(mean(Zenith_time))).*24 ,sTEC_5min_Zenith_bp(tmin_Zenith_NoNaN:tmax_Zenith_NoNaN) ); %plot all at once on this
set(h(1),'Color','r','Marker','.','LineWidth',1.25,'MarkerSize',12)
ylabel('Ampltiude (TECU/Unitless)'); %sped up significantly when merged above
title(['delta-sTEC ISR Interval Matched AVG''d High-passed AND Zenith/MISA SNR HP AVG''d +/-',num2str(Zenith_avgRange),' km around ',num2str(Zenith_threeHun),' km at ',num2str(latMillstone),' deg lat/',num2str(longMillstone),'deg long']);
axis([-inf,inf,-inf,inf])
legend({'Zenith SNR High-passed AVG''d','MISA SNR High-passed AVG''d','delta-sTEC ISR Interval Matched AVG''d High-passed'});
set(gca,'fontweight',font_Weight, 'FontSize',font_Size);
set(gca,'xtick',Xaxisvar);

subplot(3,1,2)
h = pcolor((Zenith_time_limTIME-floor(mean(Zenith_time))).*24,Zenith_height,Zenith_SNR_bp_limTIME');
set(h, 'EdgeColor', 'none'); %remove grid lines (which make it look hilarious)
shading interp;
lighting phong;
% colorbar;
caxis([-0.1 0.1]);
colormap('gray');
ylim([90 700]);
% ylim([90 500]);
% xlim([0 70]);
ylabel('Height (km)');
title(['Zenith SNR High-pass - Hours ',num2str(min(time_cutout_range)),' to ',num2str(max(time_cutout_range))]);
set(gca,'fontweight',font_Weight, 'FontSize',font_Size);
Xaxisvar = min(time_cutout_range):1:max(time_cutout_range); %shows off important magnitude points
% Xaxisvar = [-12, -8, -4, 0, 4, 8, 12, 16, 20, 24 28, 32, 36, 40, 44, 48]; %shows off important magnitude points
% Xaxisvar = [-12,-8,-4,0, 4, 8, 12, 16, 20, 24 28, 32, 36, 40, 44, 48, 52, 56, 60, 64, 68, 72]; %shows off important magnitude points
set(gca,'xtick',Xaxisvar);

%MISA SNR PLOT
subplot(3,1,3)
h = pcolor((MISA_time_limTIME-floor(mean(MISA_time))).*24,MISA_height,MISA_SNR_bp_limTIME');
set(h, 'EdgeColor', 'none'); %remove grid lines (which make it look hilarious)
shading interp;
lighting phong;
% colorbar;
caxis([-0.1 0.1]);
colormap('gray');
ylim([90 700]);
% ylim([90 500]);
% xlim([0 70]);
xlabel(['Time in UT - 0 Hr at 0 UT ',dateRange_zeroHr_YrMonDay_string]); %old: ,'fontweight','bold','FontSize',12 for all
ylabel('Height (km)');
title('MISA SNR High-pass');
set(gca,'fontweight',font_Weight, 'FontSize',font_Size);
% Xaxisvar = [-12, -8, -4, 0, 4, 8, 12, 16, 20, 24 28, 32, 36, 40, 44, 48]; %shows off important magnitude points
% Xaxisvar = [-12,-8,-4,0, 4, 8, 12, 16, 20, 24 28, 32, 36, 40, 44, 48, 52, 56, 60, 64, 68, 72]; %shows off important magnitude points
set(gca,'xtick',Xaxisvar);


%-------PLOT USED IN PAPERS, SAME AS ABOVE BUT NO MISA AND STEC ON TOP-----
%ZENITH SNR PLOT
% figure('units','normalized','outerposition',[0 0 1 1]); %maximizes figure window
figure(figNum);
figNum = figNum+1;
%STEC/ZENITH AVG'D/MISA AVG'D
subplot(2,1,1)
h = plot( (Zenith_time_limTIME-floor(mean(Zenith_time))).*24 ,Zenith_SNR_threeHun_AVGD(tmin_Zenith:tmax_Zenith) ); %plot all at once on this
hold on;
set(h(1),'Color','b','LineStyle','--','LineWidth',1.25)
h = plot( (MISA_time_limTIME -floor(mean(MISA_time))).*24 ,MISA_SNR_threeHun_AVGD(tmin_MISA:tmax_MISA) ); %plot all at once on this
hold on;
set(h(1),'Color','g','LineWidth',1.25)
h = plot( (Zenith_time_NoNaN(tmin_Zenith_NoNaN:tmax_Zenith_NoNaN)-floor(mean(Zenith_time))).*24 ,sTEC_5min_Zenith_bp(tmin_Zenith_NoNaN:tmax_Zenith_NoNaN) ); %plot all at once on this
set(h(1),'Color','r','Marker','.','LineWidth',1.25,'MarkerSize',12)
ylabel('Amplitude (TECU/Unitless)'); %sped up significantly when merged above
title(['delta-sTEC ISR Interval Matched AVG''d High-passed AND Zenith/MISA SNR HP AVG''d +/-',num2str(Zenith_avgRange),' km around ',num2str(Zenith_threeHun),' km at ',num2str(latMillstone),' deg lat/',num2str(longMillstone),'deg long']);
axis([-inf,inf,-inf,inf])
legend({'Zenith SNR High-passed AVG''d','MISA SNR High-passed AVG''d','delta-sTEC ISR Interval Matched AVG''d High-passed'},'Location','southeast');
set(gca,'fontweight',font_Weight, 'FontSize',font_Size);
set(gca,'xtick',Xaxisvar);

subplot(2,1,2)
h = pcolor((Zenith_time_limTIME-dateRange_zeroHr(2)).*24,Zenith_height,Zenith_SNR_bp_limTIME');
set(h, 'EdgeColor', 'none'); %remove grid lines (which make it look hilarious)
shading interp;
lighting phong;
% colorbar;
caxis([-0.1 0.1]);
colormap('gray');
ylim([90 700]);
% ylim([90 500]);
% xlim([0 70]);
xlabel(['Time in UT - 0 Hr at 0 UT ',dateRange_zeroHr_MonDay_string]); %old: ,'fontweight','bold','FontSize',12 for all
ylabel('Height (km)');
title(['Zenith SNR High-pass - Hours ',num2str(min(time_cutout_range)),' to ',num2str(max(time_cutout_range)),' ',...
    'That Accentuate MSTIDs']);
set(gca,'fontweight',font_Weight, 'FontSize',font_Size);
Xaxisvar = min(time_cutout_range):1:max(time_cutout_range); %shows off important magnitude points
% Xaxisvar = [-12, -8, -4, 0, 4, 8, 12, 16, 20, 24 28, 32, 36, 40, 44, 48]; %shows off important magnitude points
% Xaxisvar = [-12,-8,-4,0, 4, 8, 12, 16, 20, 24 28, 32, 36, 40, 44, 48, 52, 56, 60, 64, 68, 72]; %shows off important magnitude points
set(gca,'xtick',Xaxisvar);


%ION VELOCITY
% figure('units','normalized','outerposition',[0 0 1 1]); %maximizes figure window
figure(figNum);
figNum = figNum+1;
subplot(2,1,1);
h = pcolor((Zenith_time_limTIME-dateRange_zeroHr(2)).*24,Zenith_height,Zenith_vel_limTIME');
set(h, 'EdgeColor', 'none'); %remove grid lines (which make it look hilarious)
shading interp;
lighting phong;
% colorbar;
caxis([-20 20]);
colormap('gray');
ylim([90 700]);
% xlabel('Time in EDT (UT-4 hrs)');
ylabel('Height (km)');
% tstr1 = strcat(sprintf('%s', date_str),' (Zenith Velocity)');
title(['Zenith LOS Ion Velocity - Hours ',num2str(min(time_cutout_range)),' to ',num2str(max(time_cutout_range))]);
set(gca,'xtick',Xaxisvar);
set(gca, 'fontweight','bold', 'FontSize',18);

subplot(2,1,2);
h = pcolor((MISA_time_limTIME-dateRange_zeroHr(2)).*24,MISA_height,MISA_vel_limTIME');
set(h, 'EdgeColor', 'none'); %remove grid lines (which make it look hilarious)
shading interp;
lighting phong;
% colorbar;
colormap('gray');
caxis([-20 20]);
ylim([90 700]);
xlabel(['Time in UT - 0 Hr at 0 UT ',dateRange_zeroHr_MonDay_string]); %old: ,'fontweight','bold','FontSize',12 for all
ylabel('Height (km)');
% tstr1 = strcat(sprintf('%s', date_str),' (MISA Velocity)');
title('MISA LOS Ion Velocity');
set(gca,'xtick',Xaxisvar);
set(gca, 'fontweight','bold', 'FontSize',18);


%TEC AND ISR but no RTI
% figure('units','normalized','outerposition',[0 0 1 1]); %maximizes figure window
figure(figNum);
figNum = figNum+1;
%STEC/ZENITH AVG'D/MISA AVG'D
h = plot( (Zenith_time_limTIME-floor(mean(Zenith_time))).*24 ,Zenith_SNR_threeHun_AVGD(tmin_Zenith:tmax_Zenith) ); %plot all at once on this
hold on;
set(h(1),'Color','b','LineStyle','--','LineWidth',1.25)
h = plot( (MISA_time_limTIME -floor(mean(MISA_time))).*24 ,MISA_SNR_threeHun_AVGD(tmin_MISA:tmax_MISA) ); %plot all at once on this
hold on;
set(h(1),'Color','g','LineWidth',1.25)
h = plot( (Zenith_time_NoNaN(tmin_Zenith_NoNaN:tmax_Zenith_NoNaN)-floor(mean(Zenith_time))).*24 ,sTEC_5min_Zenith_bp(tmin_Zenith_NoNaN:tmax_Zenith_NoNaN) ); %plot all at once on this
set(h(1),'Color','r','Marker','.','LineWidth',1.25,'MarkerSize',12)
xlabel(['Time in UT - 0 Hr at 0 UT ',dateRange_zeroHr_MonDay_string]); %old: ,'fontweight','bold','FontSize',12 for all
ylabel('Amplitude (TECU/Unitless)'); %sped up significantly when merged above
title(['delta-sTEC ISR Interval Matched AVG''d High-passed AND Zenith/MISA SNR HP AVG''d +/-',num2str(Zenith_avgRange),' km around ',num2str(Zenith_threeHun),' km at ',num2str(latMillstone),' deg lat/',num2str(longMillstone),'deg long']);
axis([-inf,inf,-inf,inf])
legend({'Zenith SNR High-passed AVG''d','MISA SNR High-passed AVG''d','delta-sTEC ISR Interval Matched AVG''d High-passed'},'Location','southeast');
set(gca,'fontweight',font_Weight, 'FontSize',font_Size);
set(gca,'xtick',Xaxisvar);

%NO TEC JUST ISR
figure(figNum);
figNum = figNum+1;
h = plot( (Zenith_time_limTIME-floor(mean(Zenith_time))).*24, Zenith_SNR_bp(tmin_Zenith:tmax_Zenith,Zenith_threeHunIndex) ); %plot all at once on this
hold on;
set(h(1),'Color','b','LineStyle','--','LineWidth',1.25)
h = plot( (MISA_time_limTIME -floor(mean(MISA_time))).*24 ,MISA_SNR_bp(tmin_MISA:tmax_MISA,MISA_threeHunIndex) ); %plot all at once on this
set(h(1),'Color','g','LineWidth',1.25)
xlabel(['Time in UT - 0 Hr at 0 UT ',dateRange_zeroHr_MonDay_string]); %old: ,'fontweight','bold','FontSize',12 for all
ylabel('SNR'); %sped up significantly when merged above
title(['Zenith/MISA SNR HP at ',num2str(pointAltitude),' km Altitude']);
axis([-inf,inf,-inf,inf])
legend('Zenith SNR High-passed','MISA SNR High-passed');
set(gca,'xtick',Xaxisvar);
set(gca,'fontweight',font_Weight, 'FontSize',font_Size);


%------------------SCARGLE POINT TEC/ZENITH/MISA---------------------------
Ross_Lombscargle_optimized([(Zenith_time_NoNaN(tmin_Zenith_NoNaN:tmax_Zenith_NoNaN)-floor(mean(Zenith_time))).*24 ,sTEC_5min_Zenith_bp(tmin_Zenith_NoNaN:tmax_Zenith_NoNaN)  ]); 
fid = fopen('power_freq.txt', 'r');
    f_nm = textscan(fid, '%f %f %f %f %f', 'headerlines', 15);
    P_Zenith_SNR=f_nm{:,2};
    F_Zenith_SNR=f_nm{:,1};
    T_Zenith_SNR=f_nm{:,5}*60; % in minutes
    gs_Zenith_SNR=f_nm{:,4};
    max_Zenith_SNR=T_Zenith_SNR(P_Zenith_SNR==max(P_Zenith_SNR));
    K=length(P_Zenith_SNR);
    gf_Zenith_SNR=f_nm{:,3};
fclose('all');

%prep figure #
% figure('units','normalized','outerposition',[0 0 1 1]); %maximizes figure window
figure(figNum);
figNum = figNum+1;
subplot(3,1,1);
plot(T_Zenith_SNR,P_Zenith_SNR,'b');
xlim([0 plot_Freq_Lim]);
hold on;
% plot(T_SNR_08apr08_10_or,P_SNR_08apr08_10_or,'r');
% hold on;
plot(T_Zenith_SNR,gf_Zenith_SNR,'LineWidth',1.5,'color',[0.5 0.5 0.5]);
xlim([0 plot_Freq_Lim]);  
xlabel('Periods (min)');
ylabel('Normalized Power');
title(['delta-sTEC avg''d to ISR Interval (',num2str(round(Zenith_time_delta*24*60,1)),' min) & Cut to ISR Time & High-passed - ',num2str(latMillstone),' deg lat/',num2str(longMillstone),'deg long - sTEC point radius of ',num2str(pointRadius),' km - Hours ',num2str(min(time_cutout_range)),' to ',num2str(max(time_cutout_range)),' on Day ',num2str(floor(mean(time))),', 2013 - Lomb-Scargle Freq']);
set(gca,'xtick',Xaxisvar_SCARGLE);
set(gca,'fontweight',font_Weight, 'FontSize',font_Size);


% SCARGLE IT (ZENITH!)
Ross_Lombscargle_optimized([(Zenith_time_limTIME-floor(mean(Zenith_time))).*24 ,Zenith_SNR_threeHun_AVGD(tmin_Zenith:tmax_Zenith)  ]); 
fid = fopen('power_freq.txt', 'r');
    f_nm = textscan(fid, '%f %f %f %f %f', 'headerlines', 15);
    P_Zenith_SNR=f_nm{:,2};
    F_Zenith_SNR=f_nm{:,1};
    T_Zenith_SNR=f_nm{:,5}*60; % in minutes
    gs_Zenith_SNR=f_nm{:,4};
    max_Zenith_SNR=T_Zenith_SNR(P_Zenith_SNR==max(P_Zenith_SNR));
    K=length(P_Zenith_SNR);
    gf_Zenith_SNR=f_nm{:,3};
fclose('all');

%prep 
subplot(3,1,2);
plot(T_Zenith_SNR,P_Zenith_SNR,'b');
xlim([0 plot_Freq_Lim]);
hold on;
% plot(T_SNR_08apr08_10_or,P_SNR_08apr08_10_or,'r');
% hold on;
plot(T_Zenith_SNR,gf_Zenith_SNR,'LineWidth',1.5,'color',[0.5 0.5 0.5]);
xlim([0 plot_Freq_Lim]);  
xlabel('Periods (min)');
ylabel('Normalized Power');
title(['Zenith ISR avg''d +/-',num2str(Zenith_avgRange),' km around ',num2str(round(Zenith_height(Zenith_threeHunIndex))),' km - ',num2str(latMillstone),' deg lat/',num2str(longMillstone),'deg long - Hours ',num2str(min(time_cutout_range)),' to ',num2str(max(time_cutout_range)),' on Day ',num2str(floor(mean(time))),', 2013 - Lomb-Scargle Freq']);
set(gca,'xtick',Xaxisvar_SCARGLE);
set(gca,'fontweight',font_Weight, 'FontSize',font_Size);

% SCARGLE IT (MISA!)
Ross_Lombscargle_optimized([(MISA_time_limTIME -floor(mean(MISA_time))).*24 ,MISA_SNR_threeHun_AVGD(tmin_MISA:tmax_MISA)  ]); 
fid = fopen('power_freq.txt', 'r');
    f_nm = textscan(fid, '%f %f %f %f %f', 'headerlines', 15);
    P_Zenith_SNR=f_nm{:,2};
    F_Zenith_SNR=f_nm{:,1};
    T_Zenith_SNR=f_nm{:,5}*60; % in minutes
    gs_Zenith_SNR=f_nm{:,4};
    max_Zenith_SNR=T_Zenith_SNR(P_Zenith_SNR==max(P_Zenith_SNR));
    K=length(P_Zenith_SNR);
    gf_Zenith_SNR=f_nm{:,3};
fclose('all');

%prep
subplot(3,1,3);
plot(T_Zenith_SNR,P_Zenith_SNR,'b');
xlim([0 plot_Freq_Lim]);
hold on;
% plot(T_SNR_08apr08_10_or,P_SNR_08apr08_10_or,'r');
% hold on;
plot(T_Zenith_SNR,gf_Zenith_SNR,'LineWidth',1.5,'color',[0.5 0.5 0.5]);
xlim([0 plot_Freq_Lim]);  
xlabel('Periods (min)');
ylabel('Normalized Power');
title(['MISA ISR avg''d +/-',num2str(MISA_avgRange),' km around ',num2str(round(MISA_height(MISA_threeHunIndex))),' km - ',num2str(latMillstoneMISA),' deg lat/',num2str(longMillstoneMISA),'deg long - Hours ',num2str(min(time_cutout_range)),' to ',num2str(max(time_cutout_range)),' on Day ',num2str(floor(mean(time))),', 2013 - Lomb-Scargle Freq']);
set(gca,'xtick',Xaxisvar_SCARGLE);
set(gca,'fontweight',font_Weight, 'FontSize',font_Size);

end %*********END OF FLG_gatherDataAtPoint_Filter_TimeMatch_n_Smooth_BP*********


%% Walking sTEC Scargle (4 hr) based on sTEC filtered and cut time span to match ISR time span
if( FLG_gatherDataAtPoint_Filter_ISRTimeMatch_Scargle_Walking == 1)

sTEC_timeDelta = (timeUnique_limTIME(162)-timeUnique_limTIME(161)).*24; %hr
sTEC_hrChunkTime = 4; %hr, time to split the chunks

sTEC_countTo4hr = floor(sTEC_hrChunkTime/sTEC_timeDelta); %4 is hour, this is a count

sTEC_time_size = length(timeUnique_limTIME); %since I use it a lot, var it
sTEC_4hrChunks = floor(sTEC_time_size/sTEC_countTo4hr); %# 4 hr chunks in dataset


for(i = 1:sTEC_4hrChunks)

    if(i ~= sTEC_4hrChunks)
        Ross_Lombscargle_optimized( [ (timeUnique_limTIME(((i-1)*sTEC_countTo4hr+1):(i*sTEC_countTo4hr))-floor(mean(timeUnique))).*24 ; sTEC_combined_limTIME(((i-1)*sTEC_countTo4hr+1):(i*sTEC_countTo4hr)) ]');
    else
        Ross_Lombscargle_optimized([ (timeUnique_limTIME(((i-1)*sTEC_countTo4hr+1):(sTEC_time_size))-floor(mean(timeUnique))).*24 ; sTEC_combined_limTIME(((i-1)*sTEC_countTo4hr+1):(sTEC_time_size)) ]' );
    end
    fid = fopen('power_freq.txt', 'r');
        f_nm = textscan(fid, '%f %f %f %f %f', 'headerlines', 15);
        P_Zenith_SNR=f_nm{:,2};
        F_Zenith_SNR=f_nm{:,1};
        T_Zenith_SNR=f_nm{:,5}*60; % in minutes
        gs_Zenith_SNR=f_nm{:,4};
        max_Zenith_SNR=T_Zenith_SNR(P_Zenith_SNR==max(P_Zenith_SNR));
        K=length(P_Zenith_SNR);
        gf_Zenith_SNR=f_nm{:,3};
    fclose('all');
    
%     P_Zenith_SNR_chunk(i,:) = P_Zenith_SNR;
%     T_Zenith_SNR_chunk(i,:) = T_Zenith_SNR;
%     gf_Zenith_SNR_chunk(i,:) = gf_Zenith_SNR;
    
%     figure('units','normalized','outerposition',[0 0 1 1]); %maximizes figure window
    figure(figNum);
    figNum = figNum+1;
    plot(T_Zenith_SNR,P_Zenith_SNR,'b');
    xlim([0 plot_Freq_Lim]);
    hold on;
    % plot(T_SNR_08apr08_10_or,P_SNR_08apr08_10_or,'r');
    % hold on;
    plot(T_Zenith_SNR,gf_Zenith_SNR,'LineWidth',1.5,'color',[0.5 0.5 0.5]);
    xlim([0 plot_Freq_Lim]);  
    xlabel('Periods (min)');
    ylabel('Normalized Power');
    if(i ~= sTEC_4hrChunks)
        tstr1 = strcat(sprintf('sTEC - Time Limited, %.2f', (i-1)*sTEC_hrChunkTime+(timeUnique_limTIME(1)-floor(mean(timeUnique))).*24 ),' to');
        %Zenith_time(1) adds the initial offset since this doesn't start at
        %0
        tstr1 = strcat(tstr1,sprintf(' %.2f', (i)*sTEC_hrChunkTime+(timeUnique_limTIME(1)-floor(mean(timeUnique))).*24  ));
        title(tstr1);
    else
        tstr1 = strcat(sprintf('sTEC - Time Limited, %.2f', (i-1)*sTEC_hrChunkTime+(timeUnique_limTIME(1)-floor(mean(timeUnique))).*24  ),' to');
        tstr1 = strcat(tstr1,sprintf(' %.2f', (timeUnique_limTIME(end)-floor(mean(timeUnique))).*24 ));
        title(tstr1);
    end
    
    set(gca,'fontweight',font_Weight, 'FontSize',font_Size);
    set(gca,'xtick',Xaxisvar_SCARGLE);
end

end %*********END OF FLG_gatherDataAtPoint_Filter_ISRTimeMatch_Scargle_Walking*********

