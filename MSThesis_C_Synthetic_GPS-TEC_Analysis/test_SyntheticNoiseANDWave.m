clc
clear variables
close all

%This research was supported under NSF Grant AGS 12-41407 to The Pennsylvania State University.

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
clear poolobj CPU_Threads myCluster; %keep the excess variables down


%% CONSTANTS
Re = 6371.0; %km, Earth mean Radius
%from http://nssdc.gsfc.nasa.gov/planetary/factsheet/earthfact.html

worldLatRange = [-90,90]; %arcdeg, latitude range maximum
worldLongRange = [-180,180]; %arcdeg, longitude range maximum

data = load('sTEC_supporting_values.mat');
pplat = data.pplat; %copy out
pplong = data.pplong; %copy out
time = data.time; %copy out
timeUnique = data.timeUnique; %copy out
clear data %save mems

FLG_AVG_ANYANGLE_TimeMatch_n_BP_Scargle = 1; %1 for on, 0 for off

FLG_AVG_ANYANGLE_Scargle_FFT = 0; %0 for scargle, 1 for FFT

FLG_AVG_ANYANGLE_ZENITHORMISA = 1; %0 for Zenith, 1 for MISA (only MISA read in now, so go MISA bro)

%% USER DEFINED
N = 300; %splits lat/long ranges into chunks

plotLatRange = [35,50]; %latitude limit for plotting
%-90 to 90 is world, 35 to 50 is good for USA East Coast
plotLongRange = [-85,-60]; %longitude limit for plotting
%-180 to 180 is world, -85 to -60 is good for USA East Coast
pointAltitude = 300; %km, altitude where most e-'s are (F region max - assumed)

%Location of Millstone Hill
latMillstone = 42.6233; %Deg North
longMillstone = -71.4882; %Deg East (71.4882 deg West)
MillstoneMISA_azimuth = 168.5; %deg Azimuth from North (assuming geo - none specified)
MillstoneMISA_elev = 66.26; %deg elevation

latMillstoneMISA = latMillstone + sind(MillstoneMISA_azimuth+90)*(((pointAltitude/tand(MillstoneMISA_elev))/Re)*180/pi); %deg North (see spreadsheet titled ISR Angle Calc)
longMillstoneMISA = longMillstone - cosd(MillstoneMISA_azimuth+90)*(((pointAltitude/tand(MillstoneMISA_elev))/Re)*180/pi); %deg East at MISA (Azimuth is clockwise - cosine gets a negative)

time_cutout_range = [18,35]; %hrs, cut-out this time period for comparison between the ISR & GPS data at the best possible time to compare (night)
plot_Freq_Lim = 120; %min, for scarglin and FFT
Xaxisvar_SCARGLE = 0:10:plot_Freq_Lim; %create a periodogram x axis tick

%any angle avging
avg_anyAngle = 89.999; %deg, user defined angle
avg_anyAngle = 0.0001; %deg, user defined angle
avg_anyAngle = 135; %deg, user defined angle
avg_anyAngle_Width = 4; %arcdeg, width - not an angle in this instance
avg_anyAngle_N = 200; %number of chunks to split the range into 
avg_anyAngle_45vsLatLong = 0; %0 for longitude on xaxis on a 45 degree angle (or multiple of it), 1 for latitude on xaxis
FLG_geomagneticCoords = 0; %set geo coords off
FLG_memSavr = 1; %1 on, 0 off, clears stuff will probably break next code

%synthetic wave parameters
wave_Period = 1; %hr, period of wave
wave_WaveLength = 200; %km, wavelength of wave given in paper EPS_2008_Seker
% wave_WaveLength = 700; %km, wavelength of wave
% wave_WaveLength = 1.21*2800; %km, wavelength of wave
wave_Angle = 135; %deg, angle of wave direction
wave_Phase = 0; %deg, 0 to 360, phase of wave
wave_Amp = 0.2; %delta_sTEC, amplitude of wave
%2ND WAVE
wave_Period(2) = 1; %hr, period of wave
wave_WaveLength(2) = 200; %km, wavelength of wave given in paper EPS_2008_Seker
wave_WaveLength(2) = 700; %km, wavelength of wave
% wave_WaveLength(2) = 1.21*2800; %km, wavelength of wave
wave_Angle(2) = 225; %deg, angle of wave direction
wave_Phase(2) = 0; %deg, 0 to 360, phase of wave
wave_Amp(2) = 0.0; %delta_sTEC, amplitude of wave


%random noise background parameters
noise_Background_Mean = 0; %delta_sTEC, avg of noise background
noise_Background_STDEV = 1/3.85; %delta_sTEC, standard dev of noise background, GOAL: at more than 99% of sTEC is between +/-1 delta_STEC by #*STDEV
noise_Background_STDEV = 0.3168; %delta_sTEC, standard dev of noise background taken from all sTEC data
gif_sTEC_cap = 0.5; %TECU, limit (+ and -) for plotting the gif and other plots too
gif_Scatter_Point_Size = 325; %arb. size to make points


%% Create autoticks for plotting & other fixes

%forcing of super memory saving fix
if( ( FLG_memSavr == 1 ) && min(plotLatRange) == -90 && max(plotLatRange) == 90 && min(plotLongRange) == -180 && max(plotLongRange) == 180 )
    FLG_memSavr = 2; %force upgrade to 2 which doesn't do any data range culling since the whole world is chosen
end

%if mag conversion on, prep the system to call the IGRF conversion library
if( FLG_geomagneticCoords == 1 )
    magF = IGRF(); %prep using the external conversion library
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


%% CALC N CONVERT AS NEEDED
worldLatRange_lin = linspace(min(worldLatRange),max(worldLatRange),N); %arcdeg, get a vector of pts along this
worldLongRange_lin = linspace(min(worldLongRange),max(worldLongRange),N); %arcdeg, get a vector of pts along this

[worldLongRange_mesh, worldLatRange_mesh] = meshgrid(worldLongRange_lin,worldLatRange_lin); %arcdeg, meshgrid of X (long) and Y (lat) for plotting and stuff

wave_Angle = wave_Angle.*pi/180; %rad, convert
wave_Phase = wave_Phase.*pi/180; %rad, convert

wave_WaveLength_km = wave_WaveLength; %km, record
wave_WaveLength = wave_WaveLength./Re*(180/pi); %arcdeg, convert to arcdeg convention from km arc
wave_Freq = 1./wave_Period; %1/hr, freq of wave
wave_Speed = wave_WaveLength./wave_Period; %arcdeg/hr
wave_WaveNumber = 2.*pi./wave_WaveLength; %rad/arcdeg, wave number
wave_FreqAngular = 2.*pi.*wave_Freq; %rad/hr, angular freq


%% TEST WAVEFORM MADE

grating = flipud( wave_Amp(1)*sin( wave_WaveNumber(1)*cos(wave_Angle(1))*worldLongRange_mesh + wave_WaveNumber(1)*sin(wave_Angle(1))*worldLatRange_mesh + wave_FreqAngular(1)*0 + wave_Phase(1)) ...
+ wave_Amp(2)*sin( wave_WaveNumber(2)*cos(wave_Angle(2)).*worldLongRange_mesh + wave_WaveNumber(2)*sin(wave_Angle(2)).*worldLatRange_mesh + wave_FreqAngular(2).*0 + wave_Phase(2)) ...    
);

figure;
% imagesc( grating, [-1 1] );
h = pcolor( worldLongRange_lin , worldLatRange_lin , fliplr(grating));% pseudocolor plot "stretched" to the grid
set(h, 'EdgeColor', 'none'); %remove grid lines (which make it look hilarious)
%             shading interp; %stuff for interpolation
%             lighting phong; %stuff for interpolation
caxis([-gif_sTEC_cap, gif_sTEC_cap]);
colormap gray(256);
string_Title = 'Synthetic Wave - ';
if(wave_Amp(2) == 0)
    string_Title = [string_Title,' WaveParams:',num2str(wave_Angle(1)*180/pi),'deg,',num2str(wave_WaveLength_km(1)),'km,',num2str(wave_Period(1)),'hr'];
else
    string_Title = [string_Title,' WaveParams: W1:',num2str(wave_Angle(1)*180/pi),'deg,',num2str(wave_WaveLength_km(1)),'km,',num2str(wave_Period(1)),'hr & W2:',num2str(wave_Angle(2)*180/pi),'deg,',num2str(wave_WaveLength_km(2)),'km,',num2str(wave_Period(2)),'hr'];
end
title(string_Title);
xlabel('Longitude (arcdeg)');
ylabel('Latitude (arcdeg)');
h = colorbar; %shows the color bar
ylabel(h, 'delta-sTEC (TECU)','fontweight','bold','FontSize',18); %labels color bar
xticks( -180:15:180 ); %creates x ticks -180 to 180
yticks( -90:15:90 ); %creates y ticks -90 to 90
set(gca, 'fontweight','bold', 'FontSize',18);

figure;
% imagesc( grating + normrnd(noise_Background_Mean,noise_Background_STDEV,size(grating)) , [-1 1] );
h = pcolor( worldLongRange_lin , worldLatRange_lin , fliplr(grating + normrnd(noise_Background_Mean,noise_Background_STDEV,size(grating))));% pseudocolor plot "stretched" to the grid
set(h, 'EdgeColor', 'none'); %remove grid lines (which make it look hilarious)
%             shading interp; %stuff for interpolation
%             lighting phong; %stuff for interpolation
caxis([-gif_sTEC_cap, gif_sTEC_cap]);
colormap gray(256);
string_Title = 'Synthetic Wave with Noise - ';
if(wave_Amp(2) == 0)
    string_Title = [string_Title,' WaveParams:',num2str(wave_Angle(1)*180/pi),'deg,',num2str(wave_WaveLength_km(1)),'km,',num2str(wave_Period(1)),'hr'];
else
    string_Title = [string_Title,' WaveParams: W1:',num2str(wave_Angle(1)*180/pi),'deg,',num2str(wave_WaveLength_km(1)),'km,',num2str(wave_Period(1)),'hr & W2:',num2str(wave_Angle(2)*180/pi),'deg,',num2str(wave_WaveLength_km(2)),'km,',num2str(wave_Period(2)),'hr'];
end
string_Title = [string_Title, ' - NoiseParams: Normal,Mean:',num2str(round(noise_Background_Mean,2)),',Stdev:',num2str(round(noise_Background_STDEV,2)),''];
title(string_Title);
xlabel('Longitude (arcdeg)');
ylabel('Latitude (arcdeg)');
h = colorbar; %shows the color bar
ylabel(h, 'delta-sTEC (TECU)','fontweight','bold','FontSize',18); %labels color bar
xticks( -180:15:180 ); %creates x ticks -180 to 180
yticks( -90:15:90 ); %creates y ticks -90 to 90
set(gca, 'fontweight','bold', 'FontSize',18);


%% CREATE REAL-WORLD-LIKE DATA

sTEC = wave_Amp(1)*sin( wave_WaveNumber(1)*cos(wave_Angle(1)).*pplong + wave_WaveNumber(1)*sin(wave_Angle(1)).*pplat + wave_FreqAngular(1).*(time - floor(mean(time))).*24 + wave_Phase(1)) ...
+ wave_Amp(2)*sin( wave_WaveNumber(2)*cos(wave_Angle(2)).*pplong + wave_WaveNumber(2)*sin(wave_Angle(2)).*pplat + wave_FreqAngular(2).*(time - floor(mean(time))).*24 + wave_Phase(2)) ...    
+ normrnd(noise_Background_Mean,noise_Background_STDEV,size(pplat)); %delta_sTEC, calculate synthetic sTEC


%% IMPORT MISA NEEDED INFO

load('MISA_height.mat')
load('MISA_SNR.mat')
load('MISA_SNR_bp.mat')
load('MISA_time.mat') %time was adjusted so 0 hour is on May 7th, we will adjust to be same as this date system
MISA_threeHunIndex = 47; %MISA height index that corresponds to 300 km
MISA_timeOffset = 127; %days, offset applied to get the dates the same between these two May 7th
MISA_time = MISA_timeOffset + MISA_time./24; %days, convert to same date system

time_min = (min(timeUnique)); %days, UTC I think
time_max = (max(timeUnique)); %days, UTC I think

jm = find( min(abs(time_min - MISA_time)) == abs(time_min - MISA_time) ); %get starting MISA time that corresponds to here
jn = find( min(abs(time_max - MISA_time)) == abs(time_max - MISA_time) ); %get ending MISA time that corresponds to here

MISA_time = MISA_time(jm:jn); %snip out data that is not comparable
MISA_SNR = MISA_SNR(jm:jn,:); %snip out data that is not comparable
MISA_SNR_bp = MISA_SNR_bp(jm:jn,:); %snip out data that is not comparable

tj = find( min( abs(min(MISA_time) - timeUnique) ) == abs(min(MISA_time) - timeUnique) ); %find min timeUnique that fits ISR time
tk = find( min( abs(max(MISA_time) - timeUnique) ) == abs(max(MISA_time) - timeUnique) ); %find max timeUnique that fits ISR time
timeUnique_limISR_MISA = timeUnique(tj:tk); %cut out time that matches ISR time


%% CALCULATE ANY ANGLE

disp('TIME TO RUN ANY ANGLE ALG'); 
tic
avg_anyAngle = avg_anyAngle*pi/180; %rad, convert to radians

%
avg_anyAngle_slope = tan(avg_anyAngle); %get slope of line required
%this conversion is for y=LATITUDE x=LONGITUDE line action

%idea here is if 5 degc width, half it and go 2.5 degc(x-y at this point - c denotes) 90 deg (real angle) to the req
%angle - y=Latitude x=Longitude solving for intercept with 2 points and
%slope
avg_anyAngle_upLine_int = avg_anyAngle_Width/2*sin(avg_anyAngle + pi/2) + mean(plotLatRange)  ...
    - avg_anyAngle_slope*(avg_anyAngle_Width/2*cos(avg_anyAngle + pi/2) + mean(plotLongRange)); %get intercept of upper line
%upper and lower lines are parallel
%idea here is if 5 degc width, half it and go 2.5 degc(x-y at this point - c denotes) -90 deg (real angle) to the req
%angle - y=Latitude x=Longitude solving for intercept with 2 points and
%slope
avg_anyAngle_loLine_int = avg_anyAngle_Width/2*sin(avg_anyAngle - pi/2) + mean(plotLatRange)  ...
    - avg_anyAngle_slope*(avg_anyAngle_Width/2*cos(avg_anyAngle - pi/2) + mean(plotLongRange)); %get intercept of lower line

% avg_anyAngle_LimMaxLong = max(plotLongRange); %degc, get max longitude
% avg_anyAngle_LimMinLong = min(plotLongRange); %degc, get min longitude

avg_anyAngle_LatLim_upLine = avg_anyAngle_slope*sort(plotLongRange) + avg_anyAngle_upLine_int; %degc,
%latitude range from largest possible longitude range
avg_anyAngle_LatLim_upLine( avg_anyAngle_LatLim_upLine > max(plotLatRange)) = max(plotLatRange); %degc, limit lat to current range
avg_anyAngle_LatLim_upLine( avg_anyAngle_LatLim_upLine < min(plotLatRange)) = min(plotLatRange); %degc, limit lat to current range

% avg_anyAngle_LongLim_upLine( avg_anyAngle_LongLim_upLine > max(plotLongRange)) = max(plotLongRange); %degc, limit long to current range
% avg_anyAngle_LongLim_upLine( avg_anyAngle_LongLim_upLine < min(plotLongRange)) = min(plotLongRange); %degc, limit long to current range
%theo not needed

% avg_anyAngle_LatLim_loLine = avg_anyAngle_slope*plotLongRange + avg_anyAngle_loLine_int; %degc,
% %latitude range from largest possible longitude range
% avg_anyAngle_LatLim_loLine( avg_anyAngle_LatLim_loLine > max(plotLatRange)) = max(plotLatRange); %degc, limit lat to current range
% avg_anyAngle_LatLim_loLine( avg_anyAngle_LatLim_loLine < min(plotLatRange)) = min(plotLatRange); %degc, limit lat to current range
% avg_anyAngle_LongLim_loLine = (avg_anyAngle_LatLim_loLine - avg_anyAngle_loLine_int)/avg_anyAngle_slope; %degc, get longitudes that match

%Project upper line pts to lower pts
avg_anyAngle_LatLim_loLine = [min(avg_anyAngle_LatLim_upLine) + avg_anyAngle_Width*sin(avg_anyAngle - pi/2) , max(avg_anyAngle_LatLim_upLine) + avg_anyAngle_Width*sin(avg_anyAngle - pi/2)]; %degc, calc lower pts based on upper
avg_anyAngle_LatLim_loLine( avg_anyAngle_LatLim_loLine > max(plotLatRange)) = max(plotLatRange); %degc, limit lat to current range
avg_anyAngle_LatLim_loLine( avg_anyAngle_LatLim_loLine < min(plotLatRange)) = min(plotLatRange); %degc, limit lat to current range
%redefine upper pts based on possibly adjusted lower pts
avg_anyAngle_LatLim_upLine = [min(avg_anyAngle_LatLim_loLine) + avg_anyAngle_Width*sin(avg_anyAngle + pi/2) , max(avg_anyAngle_LatLim_loLine) + avg_anyAngle_Width*sin(avg_anyAngle + pi/2)]; %degc, calc upper pts based on lower
avg_anyAngle_LongLim_upLine = (avg_anyAngle_LatLim_upLine - avg_anyAngle_upLine_int)/avg_anyAngle_slope; %degc, get longitudes that match
avg_anyAngle_LongLim_loLine = (avg_anyAngle_LatLim_loLine - avg_anyAngle_loLine_int)/avg_anyAngle_slope; %degc, get longitudes that match

avg_anyAngle_Range = [ [avg_anyAngle_LatLim_upLine' ; avg_anyAngle_LatLim_loLine'] , [avg_anyAngle_LongLim_upLine' ; avg_anyAngle_LongLim_loLine'] ]; %degc, record pts for use
% 
% avg_anyAngle_Range(3,:) = [avg_anyAngle_Range(1,1) + avg_anyAngle_Width*sin(avg_anyAngle - pi/2), avg_anyAngle_Range(1,2) + avg_anyAngle_Width*cos(avg_anyAngle - pi/2)]; %degc, record max up pts
% avg_anyAngle_Range(4,:) = [avg_anyAngle_Range(2,1) + avg_anyAngle_Width*sin(avg_anyAngle - pi/2), avg_anyAngle_Range(2,2) + avg_anyAngle_Width*cos(avg_anyAngle - pi/2)]; %degc, record max up pts
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
%     avg_anyAngle_Range(1,:) = [min(avg_anyAngle_LatLim_upLine), min(avg_anyAngle_LongLim_upLine)]; %degc, record min up pts
%     avg_anyAngle_Range(2,:) = [max(avg_anyAngle_LatLim_upLine), max(avg_anyAngle_LongLim_upLine)]; %degc, record max up pts
%     avg_anyAngle_Range(3,:) = [avg_anyAngle_Range(1,1) + avg_anyAngle_Width*sin(avg_anyAngle - pi/2), avg_anyAngle_Range(1,2) + avg_anyAngle_Width*cos(avg_anyAngle - pi/2)]; %degc, record max up pts
%     avg_anyAngle_Range(4,:) = [avg_anyAngle_Range(2,1) + avg_anyAngle_Width*sin(avg_anyAngle - pi/2), avg_anyAngle_Range(2,2) + avg_anyAngle_Width*cos(avg_anyAngle - pi/2)]; %degc, record max up pts
% else
%     %lo was limit, get pts from it 
%     avg_anyAngle_Range = avg_anyAngle_LatLim_loLine; %degc, get box pts to get data from
%     avg_anyAngle_LongRange = avg_anyAngle_LongLim_loLine; %degc, get box pts to get data from
% end


% avg_anyAngle_Range_Length = sqrt( (avg_anyAngle_Range(1,2) - avg_anyAngle_Range(2,2))^2 + (avg_anyAngle_Range(1,1) - avg_anyAngle_Range(2,1))^2 ); %degc, length of line
%other one should be same lenght bc science 

avg_anyAngle_Range_Chunks_Long_up = linspace( min(avg_anyAngle_Range(1:2,2)) , max(avg_anyAngle_Range(1:2,2)) , avg_anyAngle_N + 1)'; %degc, chunks in longitude to go between
%chose longitude because 0 deg will be stable - 90 deg would never be with
%my hella math *wasn't stable at 0 anyway lol*
avg_anyAngle_Range_Chunks_Long_lo = linspace( min(avg_anyAngle_Range(3:4,2)) , max(avg_anyAngle_Range(3:4,2)) , avg_anyAngle_N + 1)'; %degc, chunks in longitude to go between
% avg_anyAngle_Range_Chunks_Lat = linspace( min(avg_anyAngle_Range(1:2,1)) , max(avg_anyAngle_Range(1:2,1)) , avg_anyAngle_N + 1)'; %degc, chunks in latitude to go between

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
    temp_Longs_up(j,:) = [avg_anyAngle_Range_Chunks_Long_up(j) , avg_anyAngle_Range_Chunks_Long_up(j+1)]; %degc, get longitudes needed upper line
    temp_Longs_lo(j,:) = fliplr([avg_anyAngle_Range_Chunks_Long_lo(j) , avg_anyAngle_Range_Chunks_Long_lo(j+1)]); %degc, get longitudes needed low
    temp_Lats_up(j,:) = avg_anyAngle_slope*temp_Longs_up(j,:) + avg_anyAngle_upLine_int; %degc, get latitudes needed up
    temp_Lats_lo(j,:) = avg_anyAngle_slope*temp_Longs_lo(j,:) + avg_anyAngle_loLine_int; %degc, get latitudes needed lower line
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

%Corral the data to the right place    
k = find( time == timeUnique_limISR_MISA(1)); %gets during a time period
%prep figure #
% figure('units','normalized','outerposition',[0 0 1 1]); %maximizes figure window
figure;
scatter(pplong(k),pplat(k),gif_Scatter_Point_Size,sTEC(k),'.');
caxis([-gif_sTEC_cap, gif_sTEC_cap]);
%         colormap(flipud(colormap('hot'))) %heat map where 0 (low) is white
colormap('jet') %heat map
h = colorbar; %shows the color bar
ylabel(h, 'delta-sTEC (TECU)','fontweight','bold','FontSize',18); %labels color bar
hold on;
%Now drawing line of interest
plot(longMillstone,latMillstone,'*','Color',[255/255,36/255,0/255],'MarkerSize',25,'LineWidth',1.9); %plots a point with a red big *
hold on;
%Plot lines of contients
geoshow(coastLines_lat,coastLines_long,'Color','k') %add continental outlines in black  
xlabel('Longitude (arcdeg)'); %old: ,'fontweight','bold','FontSize',12 for all
ylabel('Latitude (arcdeg)'); %sped up significantly when merged above
string_Title = ['Synth delta-sTEC Averaging with Angle = ',num2str(round(avg_anyAngle*180/pi,2)),' deg, Avg Total Width = ',num2str(round(avg_anyAngle_Width,2)),...
    ' arcdeg, Avg Step # = ',num2str(avg_anyAngle_N),...
    ' Avg per Step Width = ',num2str(round(sqrt( (temp_Longs_up(1,1) - temp_Longs_up(1,2))^2 + (temp_Lats_up(1,1) - temp_Lats_up(1,2))^2 ),2)),' arcdeg']; %create mecha title
if( FLG_geomagneticCoords == 1)
    string_Title = [string_Title,' [Mag Coord Aligned]']; %tack on mag coordinate warning
end
title(string_Title);

xticks( min(plotLongRange):plotLongRange_autoTick:max(plotLongRange) ); %creates x ticks automagically
yticks( min(plotLatRange):plotLatRange_autoTick:max(plotLatRange) ); %creates y ticks automagically
set(gca, 'fontweight','bold', 'FontSize',18);
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


sTECChunked_anyAngleAvg = zeros(length(timeUnique_limISR_MISA),avg_anyAngle_N); %preallocate


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

% finAnnounce_cntr = 1; %prep for progress reporter
finAnnounce_percentToUpdateAt = 10; % every % to update info at
finAnnounce_div = round(100/finAnnounce_percentToUpdateAt); %calc divisor to use
%CAN BE PAR - DISABLE FOR MORE MEM RANGE
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
        fprintf('PROGRESS: Averaging TEC\t%%Complete: %.1f%%\tRuntime: %.3f sec\n',finAnnounce_percentToUpdateAt*i/660,toc); %announce job, % complete, and time so far
    end
    
end
clear temp_sTEC temp_pplat temp_pplong temp_Long_List temp_Lat_List keepr temp_Longs_up temp_Longs_lo temp_Lats_up temp_Lats_lo
toc


%% PLOT THE SYNTHETIC DATA LIKE REGULAR DATA

%prep figure #
% figure('units','normalized','outerposition',[0 0 1 1]); %maximizes figure window
figure;
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
xlabel('Time in UT (hr) - 0 Hr at 0 UT May 7th, 2013');
ylabel('Height (km)');
title(['May 6-7-8 2013, MISA SNR 2 Hr High-pass Filtered, Latitude = ',num2str(round(latMillstoneMISA,2)),' arcdeg and Longitude = ',num2str(round(longMillstoneMISA,2)),' arcdeg'],'fontweight','bold','FontSize',12);
% set(gca, 'fontweight','bold', 'FontSize',12);
% Xaxisvar = [-12, -8, -4, 0, 4, 8, 12, 16, 20, 24 28, 32, 36, 40, 44, 48]; %shows off important magnitude points
%         Xaxisvar = [-12,-8,-4,0, 4, 8, 12, 16, 20, 24 28, 32, 36, 40, 44, 48, 52, 56, 60, 64, 68, 72]; %shows off important magnitude points
Xaxisvar = [ (round((min(MISA_time)-floor(mean(MISA_time)))*24)-mod(round((min(MISA_time)-floor(mean(MISA_time)))*24),2)):2:...
    (round((max(MISA_time)-floor(mean(MISA_time)))*24)-mod(round((max(MISA_time)-floor(mean(MISA_time)))*24),2)) ];
set(gca,'xtick',Xaxisvar);
set(gca, 'fontweight','bold', 'FontSize',18);
% 
% %add on a vertical line that shows where we are at in time
% 
ax2 = subplot(2,1,2); %lower bit
h = pcolor( (timeUnique_limISR_MISA -floor(mean(timeUnique))).*24 , avg_anyAngle_Range_Chunks_Long_Plot , sTECChunked_anyAngleAvg');% pseudocolor plot "stretched" to the grid
set(h, 'EdgeColor', 'none'); %remove grid lines (which make it look hilarious)
%             shading interp; %stuff for interpolation
%             lighting phong; %stuff for interpolation
caxis([-gif_sTEC_cap, gif_sTEC_cap]);
%         colormap(flipud(colormap('hot'))) %heat map where 0 (low) is white
colormap(ax2,'jet');
h2 = colorbar; %shows the color bar
ylabel(h2, 'delta-sTEC (TECU)','fontweight','bold','FontSize',18); %labels color bar
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

Yaxisvar = floor(min(avg_anyAngle_Range_Chunks_Long_Plot)):avg_anyAngle_Range_Chunks_Long_Plot_autoTick:ceil(max(avg_anyAngle_Range_Chunks_Long_Plot)); %creates y ticks automagically
set(gca,'ytick',Yaxisvar);

hold on;
%Now drawing line of interest
if( strcmp(avg_anyAngle_Range_Chunks_Long_Plot_Name,'Longitude') == 1 ) %if true, longitude
    line([ min((timeUnique_limISR_MISA -floor(mean(timeUnique))).*24) , max((timeUnique_limISR_MISA -floor(mean(timeUnique))).*24) ],[longMillstoneMISA,longMillstoneMISA],'Color','k','LineWidth',.25);  %plots a point with a black line
else %else latitude
    line([ min((timeUnique_limISR_MISA -floor(mean(timeUnique))).*24) , max((timeUnique_limISR_MISA -floor(mean(timeUnique))).*24) ],[latMillstoneMISA,latMillstoneMISA],'Color','k','LineWidth',.25);  %plots a point with a black line
end

xlabel('Time in UT (hr) - 0 Hr at 0 UT May 7th, 2013'); %old: ,'fontweight','bold','FontSize',12 for all
ylabel([avg_anyAngle_Range_Chunks_Long_Plot_Name,' (arcdeg)']); %sped up significantly when merged above
string_Title = ['Synth delta-sTEC Averaged on Angle of ',num2str(round(avg_anyAngle*180/pi,2)),' deg and Width of ',num2str(round(avg_anyAngle_Width,2)),' arcdeg']; %create title string
if(wave_Amp(2) == 0)
    string_Title = [string_Title,', WaveParams:',num2str(wave_Angle(1)*180/pi),'deg,',num2str(wave_WaveLength_km(1)),'km,',num2str(wave_Period(1)),'hr'];
else
    string_Title = [string_Title,', WaveParams: W1:',num2str(wave_Angle(1)*180/pi),'deg,',num2str(wave_WaveLength_km(1)),'km,',num2str(wave_Period(1)),'hr & W2:',num2str(wave_Angle(2)*180/pi),'deg,',num2str(wave_WaveLength_km(2)),'km,',num2str(wave_Period(2)),'hr'];
end
if( FLG_geomagneticCoords == 1)
    string_Title = [string_Title,' [Mag Coord Aligned]']; %tack on mag coordinate warning
end
title(string_Title);
set(gca, 'fontweight','bold', 'FontSize',18);
hold off;


%% Scargle or FFT codes!

if( FLG_AVG_ANYANGLE_TimeMatch_n_BP_Scargle == 1 )
    
    tj = find( min( abs(min(MISA_time) - timeUnique) ) == abs(min(MISA_time) - timeUnique) ); %find min timeUnique that fits ISR time
    tk = find( min( abs(max(MISA_time) - timeUnique) ) == abs(max(MISA_time) - timeUnique) ); %find max timeUnique that fits ISR time
    timeUnique_limISR_MISA = timeUnique(tj:tk); %cut out time that matches ISR time
    
    tmin = find( min( abs(min(time_cutout_range) - (timeUnique_limISR_MISA - floor(mean(timeUnique))).*24 ) ) == abs(min(time_cutout_range) - (timeUnique_limISR_MISA- floor(mean(timeUnique))).*24 ) ); %find min timeUnique that fits time limit
    tmax = find( min( abs(max(time_cutout_range) - (timeUnique_limISR_MISA -  floor(mean(timeUnique))).*24 ) ) == abs(max(time_cutout_range) - (timeUnique_limISR_MISA- floor(mean(timeUnique))).*24 ) ); %find max timeUnique that fits time limit
    timeUnique_limTIME_MISA = timeUnique_limISR_MISA(tmin:tmax); %cut out time that matches time limits (Zenith or MISA have been expanded to this time cadence)
    
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
            figure;
            plot(T_Zenith_SNR,P_Zenith_SNR,'b');
            xlim([0 plot_Freq_Lim]);
            hold on;
            % plot(T_SNR_08apr08_10_or,P_SNR_08apr08_10_or,'r');
            % hold on;
            plot(T_Zenith_SNR,gf_Zenith_SNR,'LineWidth',1.5,'color',[0.5 0.5 0.5]);
            xlim([0 plot_Freq_Lim]);  
            xlabel('Periods (min)','fontweight','bold','FontSize',18);
            ylabel('Normalized Power','fontweight','bold','FontSize',18);
            string_Title = ['Periodogram of TEC Avg''d on Angle of ',num2str(round(avg_anyAngle*180/pi,2)),' deg and Width of ',num2str(round(avg_anyAngle_Width,2)),...
                ' arcdeg, High-Passed @ ',num2str(round(1/bs,2)),' hr at Zenith ',avg_anyAngle_Range_Chunks_Long_Plot_Name];
            if( FLG_geomagneticCoords == 1)
                string_Title = [string_Title,' [Mag Coord Aligned]']; %tack on mag coordinate warning
            end
            title(string_Title);
            set(gca, 'fontweight','bold', 'FontSize',18);
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
            figure;
            plot(T_Zenith_SNR,P_Zenith_SNR,'b');
            xlim([0 plot_Freq_Lim]);
            hold on;
            % plot(T_SNR_08apr08_10_or,P_SNR_08apr08_10_or,'r');
            % hold on;
            plot(T_Zenith_SNR,gf_Zenith_SNR,'LineWidth',1.5,'color',[0.5 0.5 0.5]);
            xlim([0 plot_Freq_Lim]);  
            xlabel('Periods (min)','fontweight','bold','FontSize',18);
            ylabel('Normalized Power','fontweight','bold','FontSize',18);
            string_Title = ['Periodogram of TEC Avg''d on Angle of ',num2str(round(avg_anyAngle*180/pi,2)),' deg and Width of ',num2str(round(avg_anyAngle_Width,2)),...
                ' arcdeg, High-Passed @ ',num2str(round(1/bs,2)),' hr & Time Lim ',num2str(min(time_cutout_range)),'-',num2str(max(time_cutout_range)),...
                ' hrs at Zenith ',avg_anyAngle_Range_Chunks_Long_Plot_Name];
            if( FLG_geomagneticCoords == 1)
                string_Title = [string_Title,' [Mag Coord Aligned]']; %tack on mag coordinate warning
            end
            title(string_Title);
            set(gca, 'fontweight','bold', 'FontSize',18);
            set(gca,'xtick',Xaxisvar_SCARGLE);



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
    %         set(gca, 'fontweight','bold', 'FontSize',18);
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
    %         set(gca, 'fontweight','bold', 'FontSize',18);

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
            figure;
            plot(T_Zenith_SNR,P_Zenith_SNR,'b');
            xlim([0 plot_Freq_Lim]);
            hold on;
            % plot(T_SNR_08apr08_10_or,P_SNR_08apr08_10_or,'r');
            % hold on;
            plot(T_Zenith_SNR,gf_Zenith_SNR,'LineWidth',1.5,'color',[0.5 0.5 0.5]);
            xlim([0 plot_Freq_Lim]);  
            xlabel('Periods (min)','fontweight','bold','FontSize',18);
            ylabel('Normalized Power','fontweight','bold','FontSize',18);
            string_Title = ['Periodogram of TEC Avg''d on Angle of ',num2str(round(avg_anyAngle*180/pi,2)),' deg and Width of ',num2str(round(avg_anyAngle_Width,2)),...
                ' arcdeg, High-Passed @ ',num2str(round(1/bs,2)),' hr',...
                ' at MISA ',avg_anyAngle_Range_Chunks_Long_Plot_Name];
            if( FLG_geomagneticCoords == 1)
                string_Title = [string_Title,' [Mag Coord Aligned]']; %tack on mag coordinate warning
            end
            title(string_Title);
            set(gca, 'fontweight','bold', 'FontSize',18);
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
            figure;
            plot(T_Zenith_SNR,P_Zenith_SNR,'b');
            xlim([0 plot_Freq_Lim]);
            hold on;
            % plot(T_SNR_08apr08_10_or,P_SNR_08apr08_10_or,'r');
            % hold on;
            plot(T_Zenith_SNR,gf_Zenith_SNR,'LineWidth',1.5,'color',[0.5 0.5 0.5]);
            xlim([0 plot_Freq_Lim]);  
            xlabel('Periods (min)','fontweight','bold','FontSize',18);
            ylabel('Normalized Power','fontweight','bold','FontSize',18);
            string_Title = ['Periodogram of TEC Avg''d on Angle of ',num2str(round(avg_anyAngle*180/pi,2)),' deg and Width of ',num2str(round(avg_anyAngle_Width,2)),...
                ' arcdeg, High-Passed @ ',num2str(round(1/bs,2)),' hr & Time Lim ',num2str(min(time_cutout_range)),'-',num2str(max(time_cutout_range)),...
                ' hrs at MISA ',avg_anyAngle_Range_Chunks_Long_Plot_Name];
            if( FLG_geomagneticCoords == 1)
                string_Title = [string_Title,' [Mag Coord Aligned]']; %tack on mag coordinate warning
            end
            title(string_Title);
            set(gca, 'fontweight','bold', 'FontSize',18);
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

            figure;
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
            set(gca, 'fontweight','bold', 'FontSize',18);
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

            figure;
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
            set(gca, 'fontweight','bold', 'FontSize',18);
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
    %         set(gca, 'fontweight','bold', 'FontSize',18);
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
    %         set(gca, 'fontweight','bold', 'FontSize',18);

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

            figure;
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
            set(gca, 'fontweight','bold', 'FontSize',18);
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

            figure;
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
            set(gca, 'fontweight','bold', 'FontSize',18);
            % axis 'auto y'
        end
        
    end

end %*********END OF FLG_AVG_LONG_TimeMatch_n_BP_Scargle*********