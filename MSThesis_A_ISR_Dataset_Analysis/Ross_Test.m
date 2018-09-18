clear variables;
close all;
clc; 

%This research was supported under NSF Grant AGS 12-41407 to The Pennsylvania State University.

%% =====================PARALLEL COMPUTING=================================
% warning('off','parallel:convenience:RunningJobsExist'); %my coding takes no prisoners
% CPUs = 4; %Number of worker threads to use (4 max cause wut)
% if matlabpool('size') ~= 0 % checking to see if my pool is already open
%     if(matlabpool('size') ~= CPUs) %Check to see if open w/ desired worker num
%         matlabpool close; %closes it if one is open w/ undersired size
%         matlabpool(CPUs); %Start new matlab pool
%     end
% else
%     matlabpool(CPUs); %Start new matlab pool
% end

%% ======================SETTINGS & PREP===================================

%SETTINGS
plot_Freq_Lim = 120; %min, limits the scargle'd freq plots
bs = 2; %hr, high-pass cutoff period (any higher period is removed - 2 is usual)
walkingScargleHour = 4; %hr, time to split the span into
plot_SNR_Lim = 0.1; %sets SNR limit for SNR plots

Zenith_FFT_Lim = 0.025; %limits the amplitude for FFT plotting in 3D
MISA_FFT_Lim = 0.035; %limits the amplitude for FFT plotting in 3D
FFT_coeff = 1; %mult FFT # by for more FFT action *unused atm*
Zenith_Scargle_Lim = 0.007; %limits the amplitude for scargle plotting in 3D
MISA_Scargle_Lim = 0.012;  %limits the amplitude for scargle plotting in 3D

%This should set the 0 at May 7th 2013, since time is in UT theoretically
%Data begins at 11.5 hr something May 6th UT
timeOffset = 24; %hrs, offsets the times by 24 hours in certain areas for usefulness

%% ====================DISENGAGE PLOTS HERE================================

FLG_Scargle_OR_FFT = 0; %0 for Scargle, 1 for FFT
FLG_Walking_Zenith = 0; %1 for on, 0 for off
FLG_Walking_MISA = 0; %1 for on, 0 for off
FLG_Walking_CherryPick_Zenith = 0; %1 for on, 0 for off
FLG_Walking_CherryPick_MISA = 0; %doesn't exist yet
FLG_Scargle_AllAlt = 0; %1 for on, 0 for off - FFT version on at all times b/c it fast


%% =====================READING THE DATA===================================
% filename ='mlh130506a.000.hdf5';
% filename ='mlh130506b.000.hdf5';
% filename ='mlh130506c.000.hdf5';
filename ='mlh130506g.001.hdf5'; %this one is read in prev. versions of files so I went with it
% filename ='mlh130506g.002.hdf5'; %got this from madrigal site - claims to be final version

fileinfo  = h5info(filename);

fns=textscan(filename, '%3s %2d %2d %2d %10s'); %rips apart the filename to get the date (for mlh130506g.001.hdf5 it is mlh 13 05 06 g.001.hdf5)
year=double(2000+fns{1,2}); %Year from 2nd piece of file name (with same example, 13 or 2013)
month=double(fns{1,3}); %Month from 3rd piece of file name (with same example, 05)
day=double(fns{1,4}); %Day from 4th piece of file name (with same example, 06)
date=[year, month, day, 00, 00, 00]; %makes the date the year, month, day from the file name and says time is 0 hr 0 min 0 sec

range = double(h5read(filename,'/Data/Array Layout/2D Parameters/range')); %Reads into the h5 file.
time = double(h5read(filename,'/Data/Array Layout/timestamps')); %The file has the structure given.
SNR = double(h5read(filename,'/Data/Array Layout/2D Parameters/snp3')); %The double enforces that each value is a double snp3
Ti = double(h5read(filename,'/Data/Array Layout/2D Parameters/ti')); %Some values will be imported as integers else
el1 = double(h5read(filename,'/Data/Array Layout/1D Parameters/el1'));
az1 = double(h5read(filename,'/Data/Array Layout/1D Parameters/az1'));
el2 = double(h5read(filename,'/Data/Array Layout/1D Parameters/el2'));
vel = double(h5read(filename,'/Data/Array Layout/2D Parameters/vo'));
dopplar = double(h5read(filename,'/Data/Array Layout/2D Parameters/vdopp'));
% mdtyp = double(h5read(filename,'/Data/Array Layout/1D Parameters/mdtyp')); %this tells you what radar pulses were used - 115 for single pulse, 97 for alternating code

% There is a offset in time. The starting time is the total number of
% seconds counted from 01 Jan 1970.
% This experiment started on day month year. Therefore the difference between
% day month year at 00:00:00 and 01 Jan 1970 at 00:00:00 will give the offset.
% Julian day is calculated for both dates. The difference between Julian days
% gives the total number days (taking care of all leap years) and hence the
% offset time in seconds...
format long g
offset=(date2jd(year,month,day,00,00,00)-date2jd(1970,01,01,00,00,00))*24*3600; %Uses a fn someone wrote to convert date to Julian date, then finds the time difference in sec
%Julian time is used incidentally, as this function just gives the seconds
%in the date - so that they can be subtracted directly.
time=(time-offset)/3600; %removes offset, converts to hours

%Excess hours (more than 24 hours) need to be converted to days
day1 = day + floor(double(time(end))/24);
date1=[year, month, day1, 00, 00, 00];

% time_sm=smooth(time,2); %This is commented out, as unused in later code

%Clearing out excess variables
clear day1
clear offset
clear day
clear month
clear year
clear fns
clear fileinfo
clear filename

%% ========================Replacing NAN Values with 0======================
SNR(isnan(SNR)) = 0; %Replaces NAN with 0
SNR = SNR'; %Reorientates

range(isnan(range)) = 0; %Replaces NAN with 0
range = range'; %Reorientates

vel(isnan(vel)) = 0; %Replaces NAN with 0
vel = vel'; %Reorientates

dopplar(isnan(dopplar)) = 0; %Replaces NAN with 0
dopplar = dopplar'; %Reorientates

%% =================SPLITTING DATA INTO ZENITH OR MISA======================
% SNR, Range and time have both Zenith and MISA information. Based on the
% elevation angle (el1), Zenith or MISA SNR, Range and time are decoupled.
%==========================================================================

k = find( el1(:)<85 ); %Zenith is greater than 85 for eli, so it finds less than and removes it.
Zenith_SNR = SNR; %Zenith SNR copied
Zenith_SNR(k,:) = 0; %0's non-Zenith values
Zenith_range = range; %Zenith range
Zenith_range(k,:) = 0;
Zenith_time = time; %Zenith time
Zenith_time(k,:) = 0;
Zenith_el = el1; %Zenith elevation angle
Zenith_el(k,:) = 0;
Zenith_vel = vel; %Zenith velocity
Zenith_vel(k,:) = 0;
Zenith_dopplar = dopplar; %Zenith dopplar
Zenith_dopplar(k,:) = 0;

k = find( el1(:)>85 ); %MISA is less than 85 for eli, so it finds greater than and removes it.
MISA_SNR = SNR; %MISA SNR copied
MISA_SNR(k,:) = 0; %Removes non-MISA values
MISA_range = range; %MISA range
MISA_range(k,:) = 0;
MISA_time = time; %MISA time
MISA_time(k,:) = 0;
MISA_el = el1; %MISA elevation angle
MISA_el(k,:) = 0;
MISA_vel = vel; %MISA velocity
MISA_vel(k,:) = 0;
MISA_dopplar = dopplar; %MISA dopplar
MISA_dopplar(k,:) = 0;

clear SNR; %Memory clear, time is used later on so left
clear range; %This is to help compress memory in the event of large data sets
clear Ti;
clear el1;
clear el2;
clear vel;
clear dopplar;


%% ===============Clear out when Zenith SNR is 0===========================
k = find( sum(Zenith_SNR(:,:),2)== 0 ); %Finds when sum of Zenith_SNR row is 0
Zenith_SNR(k,:) = []; %Clearing 0's out!
Zenith_range(k,:) = [];
Zenith_vel(k,:) = [];
Zenith_dopplar(k,:) = [];
Zenith_time(k,:) = [];
Zenith_el(k,:) = []; %Clearing 0's

%% ==============Clear out when Zenith range is 0==========================
[row, ~]=size(Zenith_range); %Gets size of range, as everything is linearized during removal
Zenith_SNR(Zenith_range == 0)=[]; %Removes
Zenith_SNR=reshape(Zenith_SNR,row,[]); %Reconstitutes


Zenith_vel(Zenith_range == 0)=[];
Zenith_vel=reshape(Zenith_vel,row,[]);

Zenith_dopplar(Zenith_range == 0)=[];
Zenith_dopplar=reshape(Zenith_dopplar,row,[]);

Zenith_range(Zenith_range == 0)=[]; %This is last because the rest depend on it.
Zenith_range=reshape(Zenith_range,row,[]);

%% ===============Clear out when MISA SNR is 0=============================
k = find( sum(MISA_SNR(:,:),2)== 0 ); %Finds when sum of MISA_SNR column is 0
MISA_SNR(k,:) = []; %Clearing 0's out!
MISA_range(k,:) = [];
MISA_vel(k,:) = [];
MISA_dopplar(k,:) = [];
MISA_time(k,:) = [];
MISA_el(k,:) = []; %Clearing 0's

%% ==============Clear out when MISA range is 0============================
[row, ~]=size(MISA_range); %Gets size of range, as everything is linearized during removal
MISA_SNR(MISA_range == 0)=[]; %Removes
MISA_SNR=reshape(MISA_SNR,row,[]); %Reconstitutes

MISA_vel(MISA_range == 0)=[];
MISA_vel=reshape(MISA_vel,row,[]);

MISA_dopplar(MISA_range == 0)=[];
MISA_dopplar=reshape(MISA_dopplar,row,[]);

MISA_range(MISA_range == 0)=[]; %This is last because the rest depend on it.
MISA_range=reshape(MISA_range,row,[]);

%% ============Convert Time to EDT (UT-4 hrs)==============================
% Zenith_time=Zenith_time-4;
% MISA_time=MISA_time-4;
% Disabled currently

%% ===============Zenith and MISA log of SNR Calcs=========================
Zenith_SNR(Zenith_SNR <= 0) = 0.1; %Removes negative numbers/0's and replaces with 0.1
Zenith_SNR_log = log10(Zenith_SNR);

MISA_SNR(MISA_SNR <= 0) = 0.1; %Removes negative numbers/0's and replaces with 0.1
MISA_SNR_log = log10(MISA_SNR);

%% ==============Zenith and MISA Height Calcs (from range + angle)=========
Zenith_height(:,1) = (Zenith_range(1,:))'./((1.+((tand(90.-Zenith_el(1:length(Zenith_range(1,:)),1))).^2)).^0.5);

MISA_height(:,1) =(MISA_range(1,:))'./((1.+((tand(90.-MISA_el(1:length(MISA_range(1,:)),1))).^2)).^0.5);


%% ===================Median Filter (kinda) to Remove Noise================

%user input here yo
% cutOut_time_points = [-12,-4.445,0.7761,10.24,16.11,20.35,24.38,34.06,40.8,43.73]; %hr, approximate time to cut out for, corresponds to two heights
% cutOut_height_points_lower = [640,608,698,501,555,577,703.3,501.1,550,600]; %km, approximate LOWER height to cut out for, corresponds to each time
% cutOut_height_points_upper = [1000,1000,1000,1000,1000,1000,1000,1000,1000,1000]; %km, approximate UPPER height to cut out for, corresponds to each time
cutOut_time_points = [-1.726,1.32,2.733,4.474,5.996,10.67,12.09]; %hr, approximate time to cut out for, corresponds to two heights
cutOut_height_points_lower = [1,1,1,1,1,1,1]; %km, approximate LOWER height to cut out for, corresponds to each time
cutOut_height_points_upper = [200,253.9,325.8,240.4,420.2,411.2,177.5]; %km, approximate UPPER height to cut out for, corresponds to each time
%error detection here yo
if( isequal(length(cutOut_time_points),length(cutOut_height_points_lower)) ~= isequal(length(cutOut_height_points_lower),length(cutOut_height_points_upper)) )
    rounded = round(mean([length(cutOut_time_points),length(cutOut_height_points_lower),length(cutOut_height_points_upper)]));
    if(length(cutOut_time_points) ~= rounded)
        error(['cutOut_time_points size is wrong: ',num2str(length(cutOut_time_points)),' while average size is ',num2str(rounded),'.']);
    elseif(length(cutOut_height_points_lower) ~= rounded)
        error(['cutOut_height_points_lower size is wrong: ',num2str(length(cutOut_height_points_lower)),' while average size is ',num2str(rounded),'.']);
    else
        error(['cutOut_height_points_upper size is wrong: ',num2str(length(cutOut_height_points_upper)),' while average size is ',num2str(rounded),'.']);
    end
end

for(i = 1:length(cutOut_time_points))
    cutOut_time_points(i) = Zenith_time( min(abs(cutOut_time_points(i)+timeOffset - Zenith_time)) == abs(cutOut_time_points(i)+timeOffset - Zenith_time) ) - timeOffset;
    %adjust times to be closest to the real data time points
    
    cutOut_height_points_lower(i) = Zenith_height( min(abs(cutOut_height_points_lower(i) - Zenith_height)) == abs(cutOut_height_points_lower(i) - Zenith_height) );
    %adjust heights to be closest to the real data time points

    cutOut_height_points_upper(i) = Zenith_height( min(abs(cutOut_height_points_upper(i) - Zenith_height)) == abs(cutOut_height_points_upper(i) - Zenith_height) );
    %adjust heights to be closest to the real data time points
end

%plot a shot before adjustment
%Zenith plot
Zenith_SNR_mean = mean(Zenith_SNR,1); %get the mean of each range (for the entire time)
figure;
subplot(2,2,1)
plot(Zenith_height,Zenith_SNR_mean);
xlabel('Height (km)','fontweight','bold','FontSize',12);
ylabel('Mean SNR (dB?)','fontweight','bold','FontSize',12);
set(gca, 'fontweight','bold', 'FontSize',12);
title('Mean SNR for Zenith at all Ranges','fontweight','bold','FontSize',12);
Yaxisvar = [0,2,4,6];
set(gca,'ytick',Yaxisvar);

subplot(2,2,3)
pcolor(Zenith_time-timeOffset,Zenith_height,Zenith_SNR');
shading interp;
lighting phong;
colorbar;
% caxis([-plot_SNR_Lim plot_SNR_Lim]);
colormap('gray');
ylabel('Height (km)','fontweight','bold','FontSize',12);
xlabel('Time in UT - 0 Hr at 0 UT May 7th','fontweight','bold','FontSize',12);
title('Zenith SNR pre-noise-removal','fontweight','bold','FontSize',12);
set(gca, 'fontweight','bold', 'FontSize',12);
Xaxisvar = [-12,-8,-4,0, 4, 8, 12, 16, 20, 24 28, 32, 36, 40, 44, 48, 52, 56, 60, 64, 68, 72]; %shows off important magnitude points
set(gca,'xtick',Xaxisvar);
hold on;
plot(cutOut_time_points,cutOut_height_points_lower,'r');
hold on;
plot(cutOut_time_points,cutOut_height_points_upper,'r');
hold on;
plot([cutOut_time_points(1),cutOut_time_points(1)],[cutOut_height_points_lower(1),cutOut_height_points_upper(1)],'r');
hold on;
plot([cutOut_time_points(end),cutOut_time_points(end)],[cutOut_height_points_lower(end),cutOut_height_points_upper(1)],'r');

%MISA plot
MISA_SNR_mean = mean(MISA_SNR,1); %get the mean of each range (for the entire time)
subplot(2,2,2)
plot(MISA_height,MISA_SNR_mean);
xlabel('Height (km)','fontweight','bold','FontSize',12);
ylabel('Mean SNR (dB?)','fontweight','bold','FontSize',12);
set(gca, 'fontweight','bold', 'FontSize',12);
title('Mean SNR for MISA at all Ranges','fontweight','bold','FontSize',12);
Yaxisvar = [0,2,4,6];
set(gca,'ytick',Yaxisvar);

subplot(2,2,4)
pcolor(MISA_time-timeOffset,MISA_height,MISA_SNR');
shading interp;
lighting phong;
colorbar;
% caxis([-plot_SNR_Lim plot_SNR_Lim]);
colormap('gray');
ylabel('Height (km)','fontweight','bold','FontSize',12);
xlabel('Time in UT - 0 Hr at 0 UT May 7th','fontweight','bold','FontSize',12);
title('MISA SNR pre-noise-removal','fontweight','bold','FontSize',12);
set(gca, 'fontweight','bold', 'FontSize',12);
Xaxisvar = [-12,-8,-4,0, 4, 8, 12, 16, 20, 24 28, 32, 36, 40, 44, 48, 52, 56, 60, 64, 68, 72]; %shows off important magnitude points
set(gca,'xtick',Xaxisvar);
hold on;
plot(cutOut_time_points,cutOut_height_points_lower,'r');
hold on;
plot(cutOut_time_points,cutOut_height_points_upper,'r');
hold on;
plot([cutOut_time_points(1),cutOut_time_points(1)],[cutOut_height_points_lower(1),cutOut_height_points_upper(1)],'r');
hold on;
plot([cutOut_time_points(end),cutOut_time_points(end)],[cutOut_height_points_lower(end),cutOut_height_points_upper(1)],'r');


%Build the cut out lines
cutOut_height_points_lengthOfRuns = (length(cutOut_height_points_lower)); %recorded because it wanted it so
for( i = 2:cutOut_height_points_lengthOfRuns )
    cutOut_line_lower_slope = (cutOut_height_points_lower(i) - cutOut_height_points_lower(i-1))./(cutOut_time_points(i) - cutOut_time_points(i-1)); %slope between two points
    cutOut_line_lower_intercept = cutOut_height_points_lower(i) - cutOut_line_lower_slope*cutOut_time_points(i); %intercept to make the line equation go
    
    cutOut_time_start = find( min(abs(cutOut_time_points(i-1)+timeOffset - Zenith_time)) == abs(cutOut_time_points(i-1)+timeOffset - Zenith_time) );
    %search to find closest Zenith_time index to the desired cut out time start
    cutOut_time_end = find( min(abs(cutOut_time_points(i)+timeOffset - Zenith_time)) == abs(cutOut_time_points(i)+timeOffset - Zenith_time) );
    %search to find closest Zenith_time index to the desired cut out time end
    if(i == 2)
        cutOut_line_lower = cutOut_line_lower_slope.*(Zenith_time(cutOut_time_start:cutOut_time_end-1)-timeOffset) + cutOut_line_lower_intercept; %km, heights for the cut out line at each time period
    elseif( i == cutOut_height_points_lengthOfRuns )
        cutOut_line_lower = [cutOut_line_lower;cutOut_line_lower_slope.*(Zenith_time(cutOut_time_start:cutOut_time_end)-timeOffset) + cutOut_line_lower_intercept]; %km, heights for the cut out line at each time period
    else
        cutOut_line_lower = [cutOut_line_lower;cutOut_line_lower_slope.*(Zenith_time(cutOut_time_start:cutOut_time_end-1)-timeOffset) + cutOut_line_lower_intercept]; %km, heights for the cut out line at each time period
    end
    
    cutOut_line_upper_slope = (cutOut_height_points_upper(i) - cutOut_height_points_upper(i-1))./(cutOut_time_points(i) - cutOut_time_points(i-1)); %slope between two points
    cutOut_line_upper_intercept = cutOut_height_points_upper(i) - cutOut_line_upper_slope*cutOut_time_points(i); %intercept to make the line equation go
    
    cutOut_time_start = find( min(abs(cutOut_time_points(i-1)+timeOffset - Zenith_time)) == abs(cutOut_time_points(i-1)+timeOffset - Zenith_time) );
    %search to find closest Zenith_time index to the desired cut out time start
    cutOut_time_end = find( min(abs(cutOut_time_points(i)+timeOffset - Zenith_time)) == abs(cutOut_time_points(i)+timeOffset - Zenith_time) );
    %search to find closest Zenith_time index to the desired cut out time end
    if(i == 2)
        cutOut_line_upper = cutOut_line_upper_slope.*(Zenith_time(cutOut_time_start:cutOut_time_end-1)-timeOffset) + cutOut_line_upper_intercept; %km, heights for the cut out line at each time period
    elseif( i == cutOut_height_points_lengthOfRuns )
        cutOut_line_upper = [cutOut_line_upper;cutOut_line_upper_slope.*(Zenith_time(cutOut_time_start:cutOut_time_end)-timeOffset) + cutOut_line_upper_intercept]; %km, heights for the cut out line at each time period
    else
        cutOut_line_upper = [cutOut_line_upper;cutOut_line_upper_slope.*(Zenith_time(cutOut_time_start:cutOut_time_end-1)-timeOffset) + cutOut_line_upper_intercept]; %km, heights for the cut out line at each time period
    end
end

%goal: remove mean SNR from ~800 km as it is transient noise
for(i = 1:size(Zenith_time,1))
    
     Zenith_SNR(i,:) = Zenith_SNR(i,:)-mean( Zenith_SNR(i,150:180) ); %dB, remove noise amplitude at ~800 km (160) from all data
     
     MISA_SNR(i,:) = MISA_SNR(i,:)-mean( MISA_SNR(i,150:180) ); %dB, remove noise amplitude at ~800 km (160) from all data
     %note: MISA has less height ranges and its heights are slightly
     %different   
    
end

%goal: set cutout SNR to NaN for ez detection
Zenith_SNR_cutOut = Zenith_SNR; %prep it up
MISA_SNR_cutOut = MISA_SNR; %prep it up
for(i = 1:size(cutOut_line_lower,1) ) %go for how long the cut out period demands
    
    cutOut_minHeight = find( min(abs(cutOut_line_lower(i)-Zenith_height)) == abs(cutOut_line_lower(i)-Zenith_height)); % find min height involved
    cutOut_maxHeight = find( min(abs(cutOut_line_upper(i)-Zenith_height)) == abs(cutOut_line_upper(i)-Zenith_height)); % find max height involved
    
    Zenith_SNR_cutOut(i,cutOut_minHeight:cutOut_maxHeight) = NaN; %set the stuff desired to be adjusted to NaN for detection
   
%     for(j = cutOut_minHeight:cutOut_maxHeight)
%     Zenith_SNR(i,j) = mean(Zenith_SNR(:,round(rand(1)*5)+85)).*(rand(1)*2); %put some random data into the area to be cut out
%     %data is seeded at Zenith_height 85 which is 469 km, approx 475 km
%     end

    MISA_SNR_cutOut(i,cutOut_minHeight:cutOut_maxHeight) = NaN; %set the stuff desired to be adjusted to NaN for detection

end

% cutOut_minHeight = find( min(abs(min(cutOut_line_lower)-Zenith_height)) == abs(min(cutOut_line_lower)-Zenith_height)); % find min height involved
% cutOut_minHeight = cutOut_minHeight(1);
% cutOut_maxHeight = find( min(abs(max(cutOut_line_upper)-Zenith_height)) == abs(max(cutOut_line_upper)-Zenith_height)); % find max height involved
% cutOut_maxHeight = cutOut_maxHeight(1);
% % error('n');
% %goal: go height by height and interpolate stuff
% cntr = 10; %prep
% if( cutOut_minHeight < 11 )
%     cntr = cutOut_minHeight - 1;
%     if(cntr < 1)
%        cntr = 0; 
%     end
% end
% for(i = cutOut_minHeight:cutOut_maxHeight )
%     
%     jk = find( isnan(Zenith_SNR_cutOut(:,i)) == 1 ); %find the NaN's I just put in
%     jm = find( isnan(Zenith_SNR_cutOut(:,i)) ~= 1 ); %find the not NaN's I just didn't put in
%     
%     cutOut_dataBloque = Zenith_SNR(:,i-cntr:i); %cut out the data to duck with
%     
% %     if( abs(length(jk) - length(Zenith_time)) > 5 )
% %         Zenith_SNR(jk,i) = interp1(Zenith_time,Zenith_SNR(:,i),Zenith_time(jk),'spline');
% %     end
% 
% %     Zenith_SNR(jk,i) = interp2(Zenith_time,Zenith_height(i-cntr:i),cutOut_dataBloque',Zenith_time(jk),repmat(Zenith_height(cutOut_minHeight),length(jk),1),'spline');
%     %2D interpolation for missing data that utilizes
%     
% 
%     
%     if( abs(length(jk) - length(Zenith_time)) > 10 )
%         tmp = interp1(Zenith_time(jm),Zenith_SNR(jm,i),Zenith_time(jk),'linear');
%         if( sum(isnan(tmp)) < 0 )
%             Zenith_SNR(jk,i) = tmp;
%         else
% %             Zenith_SNR(jk,i) = interp2(Zenith_time,Zenith_height(i-cntr:i),cutOut_dataBloque',Zenith_time(jk),repmat(Zenith_height(cutOut_minHeight),length(jk),1),'spline');
%         end
%     end
% 
%     if( i - cntr+1 > 1 )
%         cntr = cntr + 1; %custom counter increment
%     end
%     
% end
% % cutOut_dataBloque = Zenith_SNR(:,(cutOut_minHeight-1):cutOut_maxHeight); %cut out the data to duck with
% tmp = interp2(Zenith_time,Zenith_height,Zenith_SNR,Zenith_time,Zenith_height(((cutOut_minHeight-1):cutOut_maxHeight)),'spline');
% Zenith_SNR(:,((cutOut_minHeight-1):cutOut_maxHeight)) 

%plot a shot after adjustment
Zenith_SNR_mean = mean(Zenith_SNR,1); %get the mean of each range (for the entire time)
figure;
subplot(2,2,1)
plot(Zenith_height,Zenith_SNR_mean);
xlabel('Height (km)','fontweight','bold','FontSize',12);
ylabel('Mean SNR (dB?)','fontweight','bold','FontSize',12);
set(gca, 'fontweight','bold', 'FontSize',12);
title('Mean SNR for Zenith at all Ranges, post noise adjustment','fontweight','bold','FontSize',12);
Yaxisvar = [0,2,4,6];
set(gca,'ytick',Yaxisvar);

subplot(2,2,3)
pcolor(Zenith_time-timeOffset,Zenith_height,Zenith_SNR');
shading interp;
lighting phong;
colorbar;
% caxis([-plot_SNR_Lim plot_SNR_Lim]);
colormap('gray');
ylabel('Height (km)','fontweight','bold','FontSize',12);
xlabel('Time in UT - 0 Hr at 0 UT May 7th','fontweight','bold','FontSize',12);
title('Zenith SNR post-noise-removal','fontweight','bold','FontSize',12);
set(gca, 'fontweight','bold', 'FontSize',12);
Xaxisvar = [-12,-8,-4,0, 4, 8, 12, 16, 20, 24 28, 32, 36, 40, 44, 48, 52, 56, 60, 64, 68, 72]; %shows off important magnitude points
set(gca,'xtick',Xaxisvar);
hold on;
plot(cutOut_time_points,cutOut_height_points_lower,'r');
hold on;
plot(cutOut_time_points,cutOut_height_points_upper,'r');
hold on;
plot([cutOut_time_points(1),cutOut_time_points(1)],[cutOut_height_points_lower(1),cutOut_height_points_upper(1)],'r');
hold on;
plot([cutOut_time_points(end),cutOut_time_points(end)],[cutOut_height_points_lower(end),cutOut_height_points_upper(1)],'r');

%MISA plot
MISA_SNR_mean = mean(MISA_SNR,1); %get the mean of each range (for the entire time)
subplot(2,2,2)
plot(MISA_height,MISA_SNR_mean);
xlabel('Height (km)','fontweight','bold','FontSize',12);
ylabel('Mean SNR (dB?)','fontweight','bold','FontSize',12);
set(gca, 'fontweight','bold', 'FontSize',12);
title('Mean SNR for MISA at all Ranges, post noise adjustment','fontweight','bold','FontSize',12);
Yaxisvar = [0,2,4,6];
set(gca,'ytick',Yaxisvar);

subplot(2,2,4)
pcolor(MISA_time-timeOffset,MISA_height,MISA_SNR');
shading interp;
lighting phong;
colorbar;
% caxis([-plot_SNR_Lim plot_SNR_Lim]);
colormap('gray');
ylabel('Height (km)','fontweight','bold','FontSize',12);
xlabel('Time in UT - 0 Hr at 0 UT May 7th','fontweight','bold','FontSize',12);
title('MISA SNR post-noise-removal','fontweight','bold','FontSize',12);
set(gca, 'fontweight','bold', 'FontSize',12);
Xaxisvar = [-12,-8,-4,0, 4, 8, 12, 16, 20, 24 28, 32, 36, 40, 44, 48, 52, 56, 60, 64, 68, 72]; %shows off important magnitude points
set(gca,'xtick',Xaxisvar);
hold on;
plot(cutOut_time_points,cutOut_height_points_lower,'r');
hold on;
plot(cutOut_time_points,cutOut_height_points_upper,'r');
hold on;
plot([cutOut_time_points(1),cutOut_time_points(1)],[cutOut_height_points_lower(1),cutOut_height_points_upper(1)],'r');
hold on;
plot([cutOut_time_points(end),cutOut_time_points(end)],[cutOut_height_points_lower(end),cutOut_height_points_upper(1)],'r');


%plot a shot after adjustment (only SNR bit)
figure;
subplot(2,1,1)
h = pcolor(Zenith_time-timeOffset,Zenith_height,Zenith_SNR');
% set(h, 'EdgeColor', 'none');
shading interp;
lighting phong;
colorbar;
caxis([0 6]);
colormap('gray');
ylabel('Height (km)');
xlabel('Time in UT - 0 Hr at 0 UT May 7th');
title('Zenith SNR post-noise-removal');
set(gca, 'fontweight','bold', 'FontSize',12);
Xaxisvar = -12:2:44; %shows off important magnitude points
set(gca,'xtick',Xaxisvar);
hold on;
plot(cutOut_time_points,cutOut_height_points_lower,'r');
hold on;
plot(cutOut_time_points,cutOut_height_points_upper,'r');
hold on;
plot([cutOut_time_points(1),cutOut_time_points(1)],[cutOut_height_points_lower(1),cutOut_height_points_upper(1)],'r');
hold on;
plot([cutOut_time_points(end),cutOut_time_points(end)],[cutOut_height_points_lower(end),cutOut_height_points_upper(1)],'r');
set(gca, 'fontweight','bold', 'FontSize',18);

%MISA plot

subplot(2,1,2)
h = pcolor(MISA_time-timeOffset,MISA_height,MISA_SNR');
% set(h, 'EdgeColor', 'none');
shading interp;
lighting phong;
colorbar;
caxis([0 6]);
colormap('gray');
ylabel('Height (km)');
xlabel('Time in UT - 0 Hr at 0 UT May 7th');
title('MISA SNR post-noise-removal');
set(gca, 'fontweight','bold', 'FontSize',12);
Xaxisvar = -12:2:44; %shows off important magnitude points
set(gca,'xtick',Xaxisvar);
hold on;
plot(cutOut_time_points,cutOut_height_points_lower,'r');
hold on;
plot(cutOut_time_points,cutOut_height_points_upper,'r');
hold on;
plot([cutOut_time_points(1),cutOut_time_points(1)],[cutOut_height_points_lower(1),cutOut_height_points_upper(1)],'r');
hold on;
plot([cutOut_time_points(end),cutOut_time_points(end)],[cutOut_height_points_lower(end),cutOut_height_points_upper(1)],'r');
set(gca, 'fontweight','bold', 'FontSize',18);


%% ===============High Pass Filtering to Excentuate Hourly Periods=========
%These are unused, but accentuate the periods desired
%Original file used highpass_fir.m, it has been merged
% lp=1; % Lower period (in hrs)
% hp=2; % Higher period (in hrs)
% lf=(1/hp);  % Lowpass frequency corner (1/hr)
% hf=(1/lp);  % Highpass frequency corner (1/hr)


bs=(1/bs);   % highpass cuttoff frequency (not sure what to make of it)
%I think it is related to lf above (which is 1/2 )
ls=(1/(1/2)); %low-pass cutoff frequency

Zenith_delt=Zenith_time(200)-Zenith_time(199); %Calculates a time delta

%% ===============Highpass filtering on original Zenith SNR================

n=42; % order of the Hamming window used in the custom function
% c = 3.32*pi; %Hamming constant
% M = n/2;
% bandwidth = c/M;
%The above was included to investigate the bandwidth - 1/2. Related to
%bs/lf?

fp=bs; % stop band stoppin freq (or pass band passing freq, depending on how you look at it)

f= 1/Zenith_delt; % the sampling frequency, based off of the time delta calc'd

wp=2*fp/f; % Normalizing the frequencies (Matlab is 0 to 1)
wl = 2*ls/f; %norm this
%Calculation of filter coefficients
% [b,a]=fir1(n,wp,'high'); %This takes the order and lower cut off
%Uses the default hamming window
% frequency ORIG
W = hann(n+1);
% [b,a] = fir1(n,[wp,wl],W); %option for band-pass (30 min to 2 hr)
[b,a]=fir1(n,wp,'high',W); %just high-pass (2 hr and lower OK)
%Applys Hanning Window for shiggles and giggles

Zenith_SNR_bp(:,:) = filtfilt(b,a,Zenith_SNR(:,:)); %Appies the filter
% Zenith_SNR_bp = Zenith_SNR; %No filter option to see what it really looks
% like


MISA_delt=MISA_time(200)-MISA_time(199);

%% ===============Highpass filtering on original MISA SNR==================
n=42; % order of the Hamming window;

% if(MISA_delt ~= Zenith_delt) %If they're not the same, it re-calcs the window
    fp=bs; % stop band stoppin freq (or pass band passing freq, depending on how you look at it)

    f= 1/MISA_delt; % the sampling frequency, based off of the time delta calc'd
    %Should be same as Zenith

    wp=2*fp./f; % Normalizing the frequencies

    %Calculation of filter coefficients
    [b,a]=fir1(n,wp,'high'); %This takes the order and lower cut off frequency
    %Uses the default hamming window
    %Should be the same as Zenith
% end

MISA_SNR_bp(:,:) = filtfilt(b,a,MISA_SNR(:,:)); %Appies the filter

%For testing with Mathematica - Mathematica visually same but averages on
%bp diff. Orig Zenith_SNR are identical.
% mean(MISA_SNR_bp(200,:))
% 
% mean(MISA_SNR_bp(320,:))
% 
% mean(MISA_SNR(200,:))
% 
% mean(MISA_SNR(320,:))
% 
% break;

clear a %Memory cleanup
clear b
% clear bs
clear wp
clear f
clear fp
clear n


%% ====================Optional Time Cutout Settings================
% timeCutout = [18, 35]+timeOffset; %choose UT time cutout range
% %REALLY BREAKS STUFF SORRY
% 
% Zenith_timeCutout_minMax = [ find(min(abs(Zenith_time - min(timeCutout))) == abs(Zenith_time - min(timeCutout))) ...
%      find(min(abs(Zenith_time - max(timeCutout))) == abs(Zenith_time - max(timeCutout))) ]; %get the time points for Zenith that match
% MISA_timeCutout_minMax = [ find(min(abs(MISA_time - min(timeCutout))) == abs(MISA_time - min(timeCutout))) ...
%      find(min(abs(MISA_time - max(timeCutout))) == abs(MISA_time - max(timeCutout))) ]; %get the time points for Zenith that match
% 
% Zenith_SNR = Zenith_SNR(min(Zenith_timeCutout_minMax):max(Zenith_timeCutout_minMax),:); %tonedown Zenith_SNR
% Zenith_SNR_bp = Zenith_SNR_bp(min(Zenith_timeCutout_minMax):max(Zenith_timeCutout_minMax),:); %tonedown Zenith_SNR_bp
% Zenith_time = Zenith_time(min(Zenith_timeCutout_minMax):max(Zenith_timeCutout_minMax),:); %tonedown Zenith_time
% 
% MISA_SNR = MISA_SNR(min(MISA_timeCutout_minMax):max(MISA_timeCutout_minMax),:); %tonedown MISA_SNR
% MISA_SNR_bp = MISA_SNR_bp(min(MISA_timeCutout_minMax):max(MISA_timeCutout_minMax),:); %tonedown MISA_SNR_bp
% MISA_time = MISA_time(min(MISA_timeCutout_minMax):max(MISA_timeCutout_minMax),:); %tonedown MISA_time


%% =============== FFT Anlalysis and calculating phase angle ==============
nfft = 1024; % Length of FFT
Zenith_f = (1/Zenith_delt)*(0:nfft/2-1)/nfft; %Zenith Frequency, pre-calculated
MISA_f = (1/MISA_delt)*(0:nfft/2-1)/nfft; %MISA Frequency, pre-calculated

% Zenith_f_max = zeros(1,length(Zenith_height));
% Zenith_T_max = zeros(length(Zenith_height),1);
% Zenith_angle_max = zeros(length(Zenith_height),1);
% for i=1:length(Zenith_height)
%     Zenith_X = fft(Zenith_SNR_bp(:,i),nfft);
%     Zenith_X = Zenith_X(1:nfft/2); % FFT is symmetric, throw away second half
%     Zenith_P = abs(Zenith_X);
%     Zenith_f_max(i) = Zenith_f(Zenith_P==max(Zenith_P));
%     Zenith_T_max(i,:)=60/Zenith_f_max(i); % Period corresponding to Max power;
%     Zenith_angle=angle(Zenith_X);
%     Zenith_angle_max(i,:)=(Zenith_angle(Zenith_P==max(Zenith_P)))*180/pi; % Phase angle in degrees;
% end
%This may currently be more efficient... WHOOPS!

Zenith_X = fft(Zenith_SNR_bp(:,:),nfft);
Zenith_X = Zenith_X(1:nfft/2,:); % FFT is symmetric, throw away second half
Zenith_P = abs(Zenith_X); %Power!
Zenith_f_max = (Zenith_f*ismember(Zenith_P,max(Zenith_P,[],1)) )'; %Finds the maximum frequency for each instance
Zenith_T_max = 60./Zenith_f_max; % Period corresponding to maximum power;
Zenith_angle = angle(Zenith_X); %Angle!
Zenith_angle_max = (Zenith_angle(ismember(Zenith_P,max(Zenith_P,[],1))) ).*180/pi; % Phase angle in degrees;

clear Zenith_f; %Memory control
clear Zenith_X;
clear Zenith_P;
clear Zenith_angle;
clear Zenith_delt %Deleted here, not above

MISA_X = fft(MISA_SNR_bp(:,:),nfft);
MISA_X = MISA_X(1:nfft/2,:); % FFT is symmetric, throw away second half
MISA_P = abs(MISA_X); %Power!
MISA_f_max = (MISA_f*ismember(MISA_P,max(MISA_P,[],1)) )'; %Finds the maximum frequency for each instance
MISA_T_max = 60./MISA_f_max; % Period corresponding to maximum power;
MISA_angle = angle(MISA_X); %Angle!
MISA_angle_max = (MISA_angle(ismember(MISA_P,max(MISA_P,[],1))) ).*180/pi; % Phase angle in degrees;

clear nfft;
clear MISA_f; %Memory control
clear MISA_X;
clear MISA_P;
clear MISA_angle;
clear MISA_delt %Deleted here, not above

%% ===============Commented out Time Constraints===========================
% Constraining time between 25 to 52 hrs
% Zenith_time=Zenith_time(151:399);
% Zenith_SNR=Zenith_SNR(151:399,:);
% Zenith_SNR_bp=Zenith_SNR_bp(151:399,:);
% Zenith_SNR_log=Zenith_SNR_log(151:399,:);
% Zenith_vel=Zenith_vel(151:399,:);
% Zenith_dopplar=Zenith_dopplar(151:399,:);

% MISA_time=MISA_time(252:473);
% MISA_SNR=MISA_SNR(252:473,:);
% MISA_SNR_bp=MISA_SNR_bp(252:473,:);
% MISA_SNR_log=MISA_SNR_log(252:473,:);
% MISA_vel=MISA_vel(252:473,:);
% MISA_dopplar=MISA_dopplar(252:473,:);

%% =============== Spectral Estimation of Filtered Zenith SNR at 300 km ===
% 47 corresponds to 300 km
% 25 corresponds to 200 km
% 3 corresponds to 100 km
% 70 corresponds to 400 km
% 92 corresponds to 500 km (all approximate)
filtHeightZ = find( min(abs(Zenith_height - 300)) == abs(Zenith_height - 300) ); %See 'table' above (300 km was original)

% Removing the mean from Zenith SNR time series at 300 km....
Zenith_SNR_bp(:,filtHeightZ)=Zenith_SNR_bp(:,filtHeightZ)-mean(Zenith_SNR_bp(:,filtHeightZ));

% Spectral estimation of filtered Zenith SNR at 300 km....
Zenith_data=[Zenith_time, Zenith_SNR_bp(:,filtHeightZ)];

%Uses Lombscargle instead of Fourier as it is not necessarially periodic
Ross_Lombscargle_optimized(Zenith_data); 
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

%% =============== Spectral Estimation of Filtered MISA SNR at 330 km =====
%330 km MISA approximates to 300 km Zenith. Looking at height, very rough
%guesstimation
%51 corresponds to 330 km

filtHeightM = find( min(abs(MISA_height - 330)) == abs(MISA_height - 330) ); %range at which filtering for(see 'table' above-330 orig)

% Removing the mean from MISA SNR time series at 330 km....
MISA_SNR_bp(:,filtHeightM)=MISA_SNR_bp(:,filtHeightM)-mean(MISA_SNR_bp(:,filtHeightM));

% Spectral estimation of filtered MISA SNR at 330 km (corresponds to 300 km
% for Zenith beam)....
MISA_data=[MISA_time, MISA_SNR_bp(:,filtHeightM)];

Ross_Lombscargle_optimized(MISA_data);
fid = fopen('power_freq.txt', 'r');
    f_nm = textscan(fid, '%f %f %f %f %f', 'headerlines', 15);
    P_MISA_SNR=f_nm{:,2};
    F_MISA_SNR=f_nm{:,1};
    T_MISA_SNR=f_nm{:,5}*60; % in minutes
    gs_MISA_SNR=f_nm{:,4};
    max_MISA_SNR=T_MISA_SNR(P_MISA_SNR==max(P_MISA_SNR));
    K=length(P_MISA_SNR);
    gf_MISA_SNR=f_nm{:,3};
fclose('all');


%% =============== Fixing Dates and Time it Seems =========================
% Creating date string...
if time(end)<24
   date_str=datestr(date,'dd mmm yyyy');
elseif time(end)>24 
   date_str=strcat(sprintf('%s',datestr(date,'dd- ')), datestr(date1,'dd mmm yyyy'));
end

% Observing the discontinuity in time....
a = zeros(1,length(time)-1);

for i=1:length(time)-1
    a(i)=time(i+1)-time(i);  
end
figure
plot(time(1:end-1),a);
title(['This seems to be showing discontinuities in time?'],'FontSize',14);


%% =========Testing phase difference of Zenith at 300 km using Cross Power 
% Spectral Density=========================================================

% t=Zenith_time;
Fs = 1/0.1086; %hard-coded but it's the Zenith delta time in hr

Zenith_SNR_47 = Zenith_SNR(:,filtHeightZ);
nfft = 512; %even, size of vectors will be nfft/2 + 1;
window = hamming(110);
Cxy_Zenith = zeros(nfft/2+1,length(Zenith_height));
F_Zenith = zeros(nfft/2+1,length(Zenith_height));
Axy_Zenith = zeros(nfft/2+1,length(Zenith_height));
Pxy_Zenith = zeros(nfft/2+1,length(Zenith_height));
for i=1:length(Zenith_height)
    [Cxy_Zenith(:,i),F_Zenith(:,i)] = cpsd(Zenith_SNR(:,i),Zenith_SNR_47,window,100,nfft,Fs);
    %I just CAN'T get cpsd to work with a vector. DESPITE THE DOCUMENTATION
    %SAYING IT IS GOOD TO GO. It errors on things that are correct.......
    Axy_Zenith(:,i)=angle(Cxy_Zenith(:,i))*180/pi;
    Pxy_Zenith(:,i)=abs(Cxy_Zenith(:,i));
end
% 
% gg = repmat(Zenith_SNR_47,1,length(Zenith_height)); %DOES NOT WORK BUT
% % SHOULD
% [Cxy_Zenith2,F_Zenith2]  = cpsd(abs(Zenith_SNR(1:end,1:end)),abs(gg(1:end,1:end)),window,100,nfft,Fs);

F_Zenith = F_Zenith(:,end); %Puts F_Zenith at the end of the array used.

Pxy_log_Zenith(:,:) = log10(Pxy_Zenith(:,:)); %Takes the log

F_Zenith(:) = 60./F_Zenith(:); %Converts some units... somewhere (sec -> min I think)


%% =========Testing phase difference of MISA at 300 km using Cross Power 
% Spectral Density====================================================

Fs = 1/0.1086;
MISA_SNR_51 = MISA_SNR(:,filtHeightM); %Makes a constant vector for parallel for loops
nfft = 512; %even, size of vectors will be nfft/2 + 1;
Cxy_MISA = zeros(nfft/2+1,length(MISA_height));
F_MISA = zeros(nfft/2+1,length(MISA_height));
Axy_MISA = zeros(nfft/2+1,length(MISA_height));
Pxy_MISA = zeros(nfft/2+1,length(MISA_height));
for i=1:length(MISA_height)
    [Cxy_MISA(:,i),F_MISA(:,i)] = cpsd(MISA_SNR(:,i),MISA_SNR_51,hamming(110),100,nfft,Fs);
    Axy_MISA(:,i)=angle(Cxy_MISA(:,i))*180/pi;
    Pxy_MISA(:,i)=abs(Cxy_MISA(:,i));
end
F_MISA = F_MISA(:,end);

Pxy_log_MISA(:,:) = log10(Pxy_MISA(:,:)); %accomplishes same thing

F_MISA(:) = 60./F_MISA(:);  %Accomplishes same thing


%% ==================2D Cross Correlation Between MISA and Zenith==========

xcorr2D = xcorr2(Zenith_SNR_bp,MISA_SNR_bp);
figure;
image(xcorr2D);
title('2D Correlation between Zenith SNR bp and MISA SNR bp');

figure;
plot(diag(xcorr2D));
title('Diagonal from 2D Correlation between Zenith SNR bp and MISA SNR bp');

xcorr1D = xcorr(Zenith_SNR_47,MISA_SNR_51);
figure;
plot(xcorr1D);
title('Correlation between Zenith SNR at 300 km and MISA SNR at 300 km');

corrcoef1D = corrcoef(interp1(Zenith_time,Zenith_SNR_47,MISA_time,'linear','extrap'),MISA_SNR_51);
corrcoef1D


%% =================PLOT PREP CONSTANTS COMMONLY USED======================

Zenith_time=Zenith_time-timeOffset; %hr, apply offset to force Zenith_time to be 0 hr at May 7th UT
MISA_time=MISA_time-timeOffset; %hr, apply offset to force MISA_time to be 0 hr at May 7th UT

Zenith_timeDelta = Zenith_time(148)-Zenith_time(147); %hr
MISA_timeDelta = MISA_time(148)-MISA_time(147); %hr


%% ==================Possible Plots (there are MANY!)======================

figure;
subplot(2,1,1);
pcolor(F_Zenith,Zenith_height,Axy_Zenith');
shading interp;
lighting phong;
colorbar;
xlim([0 plot_Freq_Lim]);
ylim([200 500]);
ylabel('Height (km)','fontweight','bold','FontSize',12);
xlabel('Period (min)','fontweight','bold','FontSize',12);
set(gca, 'fontweight','bold', 'FontSize',12);
title('Cross spectral phase differences between original Zenith SNR time series at various heights and Zenith SNR at 300 km','fontweight','bold','FontSize',12);


subplot(2,1,2);
pcolor(F_Zenith,Zenith_height,Pxy_Zenith');
shading interp;
lighting phong;
colorbar;
xlim([0 plot_Freq_Lim]);
ylim([200 500]);
caxis([0 0.1]);
ylabel('Height (km)','fontweight','bold','FontSize',12);
xlabel('Period (min)','fontweight','bold','FontSize',12);
set(gca, 'fontweight','bold', 'FontSize',12);
title('Cross spectral amplitude of original Zenith SNR time series at various heights and Zenith SNR at 300 km','fontweight','bold','FontSize',12);



figure;
subplot(2,1,1);
pcolor(F_MISA,MISA_height,Axy_MISA');
shading interp;
lighting phong;
colorbar;
xlim([0 plot_Freq_Lim]);
ylim([200 500]);
ylabel('Height (km)','fontweight','bold','FontSize',12);
xlabel('Period (min)','fontweight','bold','FontSize',12);
set(gca, 'fontweight','bold', 'FontSize',12);
title('Cross spectral phase differences between original MISA SNR time series at various heights and MISA SNR at 300 km','fontweight','bold','FontSize',12);


subplot(2,1,2);
pcolor(F_MISA,MISA_height,Pxy_MISA');
shading interp;
lighting phong;
colorbar;
xlim([0 plot_Freq_Lim]);
ylim([200 500]);
caxis([0 0.1]);
ylabel('Height (km)','fontweight','bold','FontSize',12);
xlabel('Period (min)','fontweight','bold','FontSize',12);
set(gca, 'fontweight','bold', 'FontSize',12);
title('Cross spectral amplitude of original MISA SNR time series at various heights and MISA SNR at 300 km','fontweight','bold','FontSize',12);



% x=Zenith_SNR_bp(1:505,47);
% y=MISA_SNR_bp(1:505,51);
% t=(Zenith_time(1:505)+MISA_time(1:505))/2;
% Fs = 1/0.1086;
% a=30;
% b=60;
% c=90;
% d=24*60;
% 
% figure;
% plot(t,x);
% hold on;
% plot(t,y,'k');
% xlabel('Time (hrs)'); 


% [Pxy,F] = mscohere(x,y,hamming(100),80,512,Fs);
% figure;
% for i=1:length(F)
%     F(i)=60/F(i);
% end
% plot(F,Pxy,'linewidth',2); 
% xlim([0 plot_Freq_Lim]);
% set(gca,'xtick',[a b c]);
% title('Magnitude-squared Coherence');
% xlabel('Period (min)'); 
% grid on;
% 
% [Cxy,F] = cpsd(x,y,hamming(100),80,512,Fs);
% for i=1:length(F)
%     F(i)=60/F(i);
% end
% figure;
% plot(F,(angle(Cxy))*180/pi,'linewidth',2); 
% xlim([0 plot_Freq_Lim]);
% title('Cross Spectrum Phase');
% xlabel('Period (min)');
% set(gca,'xtick',[a b c]);
% set(gca,'ytick',[-180 -90 -45 0 45 90 180]);
% grid on;
% 
% figure;
% plot(F,abs(Cxy),'linewidth',2); 
% xlim([0 plot_Freq_Lim]);
% set(gca,'xtick',[a b c]);
% title('Cross Spectrum Amplitude');
% xlabel('Period (min)'); grid on;
% 
% figure0=figure;
% pcolor(Zenith_time,Zenith_height,Zenith_SNR');
% shading interp;
% lighting phong;
% colorbar;
% caxis([0 10]);
% colormap('jet');
% ylim([90 700]);
% % ylim([90 500]);
% % xlim([0 70]);

%==========================================================================

% 
% 
% 
% 
% 
% 
% 
% 
% figure0=figure;
% subplot(2,1,1);
% pcolor(Zenith_time,Zenith_height,Zenith_SNR');
% shading interp;
% lighting phong;
% colorbar;
% caxis([0 10]);
% colormap('jet');
% ylim([90 700]);
% % ylim([90 500]);
% % xlim([0 70]);
% 
% ylabel('Height (km)','fontweight','bold','FontSize',12);
% tstr1 = strcat(sprintf('%s', date_str),' (Zenith SNR)');
% title(tstr1,'fontweight','bold','FontSize',12);
% set(gca, 'fontweight','bold', 'FontSize',12);
% 
% subplot(2,1,2);
% pcolor(MISA_time,MISA_height,MISA_SNR');
% shading interp;
% lighting phong;
% colorbar;
% caxis([0 10]);
% ylim([90 700]);
% xlabel('Time in EDT (UT-4 hrs)','fontweight','bold','FontSize',12);
% ylabel('Height (km)','fontweight','bold','FontSize',12);
% tstr1 = strcat(sprintf('%s', date_str),' (MISA SNR)');
% title(tstr1,'fontweight','bold','FontSize',12);
% set(gca, 'fontweight','bold', 'FontSize',12);
% 
% 
% 
% 
% figure1=figure;
% subplot(2,1,1);
% pcolor(Zenith_time,Zenith_height,Zenith_SNR_log');
% shading interp;
% lighting phong;
% colorbar;
% caxis([-1 1.5]);
% colormap('jet');
% ylim([90 700]);
% % ylim([90 500]);
% % xlim([0 70]);
% 
% ylabel('Height (km)','fontweight','bold','FontSize',12);
% tstr1 = strcat(sprintf('%s', date_str),' (Zenith log SNR)');
% title(tstr1,'fontweight','bold','FontSize',12);
% set(gca, 'fontweight','bold', 'FontSize',12);
% 
% subplot(2,1,2);
% pcolor(MISA_time,MISA_height,MISA_SNR_log');
% shading interp;
% lighting phong;
% colorbar;
% caxis([-1 1.5]);
% ylim([90 700]);
% xlabel('Time in EDT (UT-4 hrs)','fontweight','bold','FontSize',12);
% ylabel('Height (km)','fontweight','bold','FontSize',12);
% tstr1 = strcat(sprintf('%s', date_str),' (MISA log SNR)');
% title(tstr1,'fontweight','bold','FontSize',12);
% set(gca, 'fontweight','bold', 'FontSize',12);
% 
% 
% 
% 
% 

figure;

subplot(2,1,1);
pcolor(Zenith_time,Zenith_height,Zenith_SNR');
shading interp;
lighting phong;
colorbar;
caxis([-plot_SNR_Lim plot_SNR_Lim]);
colormap('gray');
ylim([90 700]);
% ylim([90 500]);
% xlim([0 70]);

ylabel('Height (km)','fontweight','bold','FontSize',12);
tstr1 = strcat(sprintf('%s', date_str),' (Zenith SNR No Filters)');
title(tstr1,'fontweight','bold','FontSize',12);
set(gca, 'fontweight','bold', 'FontSize',12);
% Xaxisvar = [-12, -8, -4, 0, 4, 8, 12, 16, 20, 24 28, 32, 36, 40, 44, 48]; %shows off important magnitude points
Xaxisvar = [-12,-8,-4,0, 4, 8, 12, 16, 20, 24 28, 32, 36, 40, 44, 48, 52, 56, 60, 64, 68, 72]; %shows off important magnitude points
set(gca,'xtick',Xaxisvar);

subplot(2,1,2);
pcolor(MISA_time,MISA_height,MISA_SNR');
shading interp;
lighting phong;
colorbar;
caxis([-plot_SNR_Lim plot_SNR_Lim]);
ylim([90 700]);
% xlabel('Time in EDT (UT-4 hrs)','fontweight','bold','FontSize',12);
xlabel('Time in UT - 0 Hr at 0 UT May 7th','fontweight','bold','FontSize',12);
ylabel('Height (km)','fontweight','bold','FontSize',12);
tstr1 = strcat(sprintf('%s', date_str),' (MISA SNR No Filters)');
title(tstr1,'fontweight','bold','FontSize',12);
set(gca, 'fontweight','bold', 'FontSize',12);
% Xaxisvar = [-12, -8, -4, 0, 4, 8, 12, 16, 20, 24 28, 32, 36, 40, 44, 48]; %shows off important magnitude points
Xaxisvar = [-12,-8,-4,0, 4, 8, 12, 16, 20, 24 28, 32, 36, 40, 44, 48, 52, 56, 60, 64, 68, 72]; %shows off important magnitude points
set(gca,'xtick',Xaxisvar);

figure2=figure;
subplot(2,1,1);
pcolor(Zenith_time,Zenith_height,Zenith_SNR_bp');
shading interp;
lighting phong;
colorbar;
caxis([-plot_SNR_Lim plot_SNR_Lim]);
colormap('gray');
ylim([90 700]);
% ylim([90 500]);
% xlim([0 70]);
xlabel(date_str)
ylabel('Height (km)');
title(['Zenith SNR with High-pass Filter ',num2str(round(1/bs)),' hr Cutoff']);
set(gca, 'fontweight','bold', 'FontSize',18);
% Xaxisvar = [-12, -8, -4, 0, 4, 8, 12, 16, 20, 24 28, 32, 36, 40, 44, 48]; %shows off important magnitude points
% Xaxisvar = [-12,-8,-4,0, 4, 8, 12, 16, 20, 24 28, 32, 36, 40, 44, 48, 52, 56, 60, 64, 68, 72]; %shows off important magnitude points
Xaxisvar = -12:2:44; %important times 
set(gca,'xtick',Xaxisvar);
set(gca, 'fontweight','bold', 'FontSize',18);

subplot(2,1,2);
pcolor(MISA_time,MISA_height,MISA_SNR_bp');
shading interp;
lighting phong;
colorbar;
caxis([-plot_SNR_Lim plot_SNR_Lim]);
ylim([90 700]);
% xlabel('Time in EDT (UT-4 hrs)','fontweight','bold','FontSize',12);
xlabel('Time in UT - 0 Hr at 0 UT May 7th');
ylabel('Height (km)');
title(['MISA SNR with High-pass Filter ',num2str(round(1/bs)),' hr Cutoff']);
% Xaxisvar = [-12, -8, -4, 0, 4, 8, 12, 16, 20, 24 28, 32, 36, 40, 44, 48]; %shows off important magnitude points
% Xaxisvar = [-12,-8,-4,0, 4, 8, 12, 16, 20, 24 28, 32, 36, 40, 44, 48, 52, 56, 60, 64, 68, 72]; %shows off important magnitude points
Xaxisvar = -12:2:44; %important times 
set(gca,'xtick',Xaxisvar);
set(gca, 'fontweight','bold', 'FontSize',18);

% 
% ZENITH 300 KM PLOT
figure3=figure;
subplot(3,2,1);
plot(Zenith_time,Zenith_SNR(:,filtHeightZ));
% hold on;
% plot(time_08apr08_10,SNR_08apr08_10_sg(:,50),'r');
% xlabel('Time in EDT (UT-4 hrs)','fontweight','bold','FontSize',11);
% xlabel('Time in UT - 0 Hr at 0 UT May 7th');
ylabel('SNR at 300 km');
title('Zenith Beam');
axis([ -inf inf -inf inf]);
set(gca,'xtick',Xaxisvar);
set(gca, 'fontweight','bold', 'FontSize',18);

subplot(3,2,3);
plot(Zenith_time,Zenith_SNR_bp(:,filtHeightZ));
% xlabel('Time in EDT (UT-4 hrs)','fontweight','bold','FontSize',12);
xlabel('Time in UT - 0 Hr at 0 UT May 7th');
ylabel('High-pass SNR at 300 km');
axis([ -inf inf -inf inf]);
set(gca,'xtick',Xaxisvar);
set(gca, 'fontweight','bold', 'FontSize',18);

subplot(3,2,5);
plot(T_Zenith_SNR,P_Zenith_SNR,'b');
xlim([0 plot_Freq_Lim]);
hold on;
% plot(T_SNR_08apr08_10_or,P_SNR_08apr08_10_or,'r');
% hold on;
plot(T_Zenith_SNR,gf_Zenith_SNR,'LineWidth',1.5,'color',[0.5 0.5 0.5]);
xlim([0 plot_Freq_Lim]);  
xlabel('Periods (min)');
ylabel('Normalized Power');
set(gca, 'fontweight','bold', 'FontSize',18);

%MISA 300 KM PLOT
% figure4=figure;
subplot(3,2,2);
plot(MISA_time,MISA_SNR(:,filtHeightM));
% hold on;
% plot(time_08apr08_10,SNR_08apr08_10_sg(:,50),'r');
% xlabel('Time in EDT (UT-4 hrs)','fontweight','bold','FontSize',11);
xlabel('Time in UT - 0 Hr at 0 UT May 7th');
ylabel('SNR at 300 km');
title('MISA Beam');
axis([ -inf inf -inf inf]);
set(gca,'xtick',Xaxisvar);
set(gca, 'fontweight','bold', 'FontSize',18);

subplot(3,2,4);
plot(MISA_time,MISA_SNR_bp(:,filtHeightM));
% xlabel('Time in EDT (UT-4 hrs)','fontweight','bold','FontSize',11);
xlabel('Time in UT - 0 Hr at 0 UT May 7th');
ylabel('High-pass SNR at 300 km');
set(gca,'xtick',Xaxisvar);
axis([ -inf inf -inf inf]);
set(gca, 'fontweight','bold', 'FontSize',18);

subplot(3,2,6);
plot(T_MISA_SNR,P_MISA_SNR,'b');
xlim([0 plot_Freq_Lim]);
hold on;
% plot(T_SNR_08apr08_10_or,P_SNR_08apr08_10_or,'r');
% hold on;
plot(T_MISA_SNR,gf_MISA_SNR,'LineWidth',1.5,'color',[0.5 0.5 0.5]);
xlim([0 plot_Freq_Lim]);  
xlabel('Periods (min)');
ylabel('Normalized Power');
set(gca, 'fontweight','bold', 'FontSize',18);

%% Bonus Walking Lomb-Scargle - ZENITH
if( FLG_Walking_Zenith == 1 )
    
Zenith_hrChunkTime = walkingScargleHour ; %hr, time to split the chunks

Zenith_countTo4hr = floor(Zenith_hrChunkTime/Zenith_timeDelta); %4 is hour, this is a count

Zenith_time_size = size(Zenith_time,1); %since I use it a lot, var it
Zenith_4hrChunks = floor(Zenith_time_size/Zenith_countTo4hr); %# 4 hr chunks in dataset


for(i = 1:Zenith_4hrChunks)
    
    if( FLG_Scargle_OR_FFT == 0 )
        
        if(i ~= Zenith_4hrChunks)
            Ross_Lombscargle_optimized(Zenith_data(((i-1)*Zenith_countTo4hr+1):(i*Zenith_countTo4hr),:));
        else
            Ross_Lombscargle_optimized(Zenith_data(((i-1)*Zenith_countTo4hr+1):(Zenith_time_size),:));
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

        figure;
        plot(T_Zenith_SNR,P_Zenith_SNR,'b');
        xlim([0 plot_Freq_Lim]);
        hold on;
        % plot(T_SNR_08apr08_10_or,P_SNR_08apr08_10_or,'r');
        % hold on;
        plot(T_Zenith_SNR,gf_Zenith_SNR,'LineWidth',1.5,'color',[0.5 0.5 0.5]);
        xlim([0 plot_Freq_Lim]);  
        xlabel('Periods (min)');
        ylabel('Normalized Power');
        if(i ~= Zenith_4hrChunks)
            tstr1 = strcat(sprintf('Zenith %.2f', (i-1)*Zenith_hrChunkTime+Zenith_time(1) ),' to');
            %Zenith_time(1) adds the initial offset since this doesn't start at
            %0
            tstr1 = strcat(tstr1,sprintf(' %.2f', (i)*Zenith_hrChunkTime+Zenith_time(1)));
            title(tstr1);
        else
            tstr1 = strcat(sprintf('Zenith %.2f', (i-1)*Zenith_hrChunkTime+Zenith_time(1)),' to');
            tstr1 = strcat(tstr1,sprintf(' %.2f', Zenith_time(end)));
            title(tstr1);
        end
        set(gca, 'fontweight','bold', 'FontSize',18);
    
    elseif( FLG_Scargle_OR_FFT == 1 )
   
        totalcount = size(Zenith_data(((i-1)*Zenith_countTo4hr+1):(i*Zenith_countTo4hr),:),1); % # of samples taken
        deltatime = Zenith_timeDelta; %sample interval in time - units hours *)
        deltafreq = 1/(totalcount*deltatime); %frequency resolution per above sampling interval in time
        frequencysamples = ((0:(totalcount-1))*deltafreq)';
        % frequencysamples = linspace(0,ceil((totalcount-1)*deltafreq),totalcount)';
        % spectrummax = frequencysamples(totalcount/2 - 1);
%         timedomain = zeros(totalcount,1);

        frequencydomain = fft(Zenith_data(((i-1)*Zenith_countTo4hr+1):(i*Zenith_countTo4hr),2));
        powerspectrum = sqrt(abs(frequencydomain.*conj(frequencydomain)))/totalcount; %properly normalized now
        
        gf=(1-(.05/length(powerspectrum))^(1/(length(powerspectrum)-1))); %for gf and gs see pg 934-935 of the journal paper "SPECTRUM: SPECTRAL ANALYSIS OF UNEVENLY SPACED PALEOCLIMATIC TIME SERIES" by SCHULZ and STATTEGGER
%         gs=0.4*gf; %Computers & Geosciences Vol. 23, No. 9, pp. 929-945, 1997 
        
        figure;
        plot(1./frequencysamples*60,powerspectrum)
        hold on;
        plot(1./frequencysamples*60,repmat(gf,size(powerspectrum)),'LineWidth',1.5,'color',[0.5 0.5 0.5]);
        if(i ~= Zenith_4hrChunks)
            tstr1 = strcat(sprintf('Zenith %.2f', (i-1)*Zenith_hrChunkTime+Zenith_time(1) ),' to');
            %Zenith_time(1) adds the initial offset since this doesn't start at
            %0
            tstr1 = strcat(tstr1,sprintf(' %.2f FFT', (i)*Zenith_hrChunkTime+Zenith_time(1)));
            title(tstr1);
        else
            tstr1 = strcat(sprintf('Zenith %.2f', (i-1)*Zenith_hrChunkTime+Zenith_time(1)),' to');
            tstr1 = strcat(tstr1,sprintf(' %.2f FFT', Zenith_time(end)));
            title(tstr1);
        end
        xlabel('Period (min)');
        ylabel('Power');
        axis([0, plot_Freq_Lim, min(powerspectrum), max(powerspectrum)])
        axis 'auto y'
        set(gca, 'fontweight','bold', 'FontSize',18);
    
    end
end


end


%% =====================Bonus FFT attempt - ZENITH=========================

totalcount = size(Zenith_data,1); % # of samples taken
deltatime = Zenith_timeDelta; %sample interval in time - units hours *)
deltafreq = 1/(totalcount*deltatime); %frequency resolution per above sampling interval in time
frequencysamples = ((0:(totalcount-1))*deltafreq)';
% frequencysamples = linspace(0,ceil((totalcount-1)*deltafreq),totalcount)';
% spectrummax = frequencysamples(totalcount/2 - 1);
timedomain = zeros(totalcount,1);

frequencydomain = fft(Zenith_data(:,2));
powerspectrum = sqrt(abs(frequencydomain.*conj(frequencydomain)))/totalcount; %properly normalized now

gf=(1-(.05/length(powerspectrum))^(1/(length(powerspectrum)-1))); %for gf and gs see pg 934-935 of the journal paper "SPECTRUM: SPECTRAL ANALYSIS OF UNEVENLY SPACED PALEOCLIMATIC TIME SERIES" by SCHULZ and STATTEGGER
%         gs=0.4*gf; %Computers & Geosciences Vol. 23, No. 9, pp. 929-945, 1997 

figure;
plot(1./frequencysamples*60,powerspectrum)
hold on;
plot(1./frequencysamples*60,repmat(gf,size(powerspectrum)),'LineWidth',1.5,'color',[0.5 0.5 0.5]);
title('FFT of Zenith SNR at 300 km','FontSize',14);
xlabel('Period (min)','FontSize',14);
ylabel('Power','FontSize',14);
axis([0, plot_Freq_Lim, min(powerspectrum), max(powerspectrum)])
axis 'auto y'


%% Bonus Walking Lomb-Scargle - MISA
if( FLG_Walking_MISA == 1 )

MISA_hrChunkTime = walkingScargleHour ; %hr, time to split the chunks

MISA_countTo4hr = floor(MISA_hrChunkTime/MISA_timeDelta); %4 is hour, this is a count

MISA_time_size = size(MISA_time,1); %since I use it a lot, var it
MISA_4hrChunks = floor(MISA_time_size/MISA_countTo4hr); %# 4 hr chunks in dataset


for(i = 1:MISA_4hrChunks)
    
    if( FLG_Scargle_OR_FFT == 0 )
        
    if(i ~= MISA_4hrChunks)
        Ross_Lombscargle_optimized(MISA_data(((i-1)*MISA_countTo4hr+1):(i*MISA_countTo4hr),:));
    else
        Ross_Lombscargle_optimized(MISA_data(((i-1)*MISA_countTo4hr+1):(MISA_time_size),:));
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
    
    figure;
    plot(T_Zenith_SNR,P_Zenith_SNR,'b');
    xlim([0 plot_Freq_Lim]);
    hold on;
    % plot(T_SNR_08apr08_10_or,P_SNR_08apr08_10_or,'r');
    % hold on;
    plot(T_Zenith_SNR,gf_Zenith_SNR,'LineWidth',1.5,'color',[0.5 0.5 0.5]);
    xlim([0 plot_Freq_Lim]);  
    xlabel('Periods (min)','fontweight','bold','FontSize',11);
    ylabel('Normalized Power','fontweight','bold','FontSize',11);
    if(i ~= MISA_4hrChunks)
        tstr1 = strcat(sprintf('MISA %.2f', (i-1)*MISA_hrChunkTime+MISA_time(1) ),' to');
        %Zenith_time(1) adds the initial offset since this doesn't start at
        %0
        tstr1 = strcat(tstr1,sprintf(' %.2f', (i)*MISA_hrChunkTime+MISA_time(1)));
        title(tstr1,'fontweight','bold','FontSize',12);
    else
        tstr1 = strcat(sprintf('MISA %.2f', (i-1)*MISA_hrChunkTime+MISA_time(1)),' to');
        tstr1 = strcat(tstr1,sprintf(' %.2f', MISA_time(end)));
        title(tstr1,'fontweight','bold','FontSize',12);
    end
    
    set(gca, 'fontweight','bold', 'FontSize',11);
    
    elseif( FLG_Scargle_OR_FFT == 1 )
        
        totalcount = size(MISA_data(((i-1)*MISA_countTo4hr+1):(i*MISA_countTo4hr),:),1); % # of samples taken
        deltatime = MISA_timeDelta; %sample interval in time - units hours *)
        deltafreq = 1/(totalcount*deltatime); %frequency resolution per above sampling interval in time
        frequencysamples = ((0:(totalcount-1))*deltafreq)';
        % frequencysamples = linspace(0,ceil((totalcount-1)*deltafreq),totalcount)';
        % spectrummax = frequencysamples(totalcount/2 - 1);
%         timedomain = zeros(totalcount,1);

        frequencydomain = fft(MISA_data(((i-1)*MISA_countTo4hr+1):(i*MISA_countTo4hr),2));
        powerspectrum = sqrt(abs(frequencydomain.*conj(frequencydomain)))/totalcount; %properly normalized now
        gf=(1-(.05/length(powerspectrum))^(1/(length(powerspectrum)-1))); %for gf and gs see pg 934-935 of the journal paper "SPECTRUM: SPECTRAL ANALYSIS OF UNEVENLY SPACED PALEOCLIMATIC TIME SERIES" by SCHULZ and STATTEGGER
        %         gs=0.4*gf; %Computers & Geosciences Vol. 23, No. 9, pp. 929-945, 1997 
        

        figure;
        plot(1./frequencysamples*60,powerspectrum)
        hold on;
        plot(1./frequencysamples*60,repmat(gf,size(powerspectrum)),'LineWidth',1.5,'color',[0.5 0.5 0.5]);
        if(i ~= MISA_4hrChunks)
            tstr1 = strcat(sprintf('MISA %.2f', (i-1)*MISA_hrChunkTime+MISA_time(1) ),' to');
            %Zenith_time(1) adds the initial offset since this doesn't start at
            %0
            tstr1 = strcat(tstr1,sprintf(' %.2f FFT', (i)*MISA_hrChunkTime+MISA_time(1)));
            title(tstr1);
        else
            tstr1 = strcat(sprintf('MISA %.2f', (i-1)*MISA_hrChunkTime+MISA_time(1)),' to');
            tstr1 = strcat(tstr1,sprintf(' %.2f FFT', MISA_time(end)));
            title(tstr1);
        end
        xlabel('Period (min)');
        ylabel('Power');
        axis([0, plot_Freq_Lim, min(powerspectrum), max(powerspectrum)])
        axis 'auto y'
        set(gca, 'fontweight','bold', 'FontSize',18);        
    end
end

end


%% =====================Bonus FFT attempt - MISA===========================

totalcount = size(MISA_data,1); % # of samples taken
deltatime = MISA_timeDelta; %sample interval in time - units hours *)
deltafreq = 1/(totalcount*deltatime); %frequency resolution per above sampling interval in time
frequencysamples = ((0:(totalcount-1))*deltafreq)';
% frequencysamples = linspace(0,ceil((totalcount-1)*deltafreq),totalcount)';
% spectrummax = frequencysamples(totalcount/2 - 1);
% timedomain = zeros(totalcount,1);

frequencydomain = fft(MISA_data(:,2));
powerspectrum = sqrt(abs(frequencydomain.*conj(frequencydomain)))/totalcount; %properly normalized now

gf=(1-(.05/length(powerspectrum))^(1/(length(powerspectrum)-1))); %for gf and gs see pg 934-935 of the journal paper "SPECTRUM: SPECTRAL ANALYSIS OF UNEVENLY SPACED PALEOCLIMATIC TIME SERIES" by SCHULZ and STATTEGGER
%         gs=0.4*gf; %Computers & Geosciences Vol. 23, No. 9, pp. 929-945, 1997 

figure;
plot(1./frequencysamples*60,powerspectrum)
hold on;
plot(1./frequencysamples*60,repmat(gf,size(powerspectrum)),'LineWidth',1.5,'color',[0.5 0.5 0.5]);
title('FFT of MISA SNR at 300 km','FontSize',14);
xlabel('Period (min)','FontSize',14);
ylabel('Power','FontSize',14);
axis([0, plot_Freq_Lim, min(powerspectrum), max(powerspectrum)])
axis 'auto y'


%% Unused stuff I guess

% 
% 
% 
% figure4=figure;
% subplot(3,1,1);
% plot(MISA_time,MISA_SNR(:,51));
% % hold on;
% % plot(time_08apr08_10,SNR_08apr08_10_sg(:,75),'r');
% set(gca, 'fontweight','bold', 'FontSize',11);
% xlabel('Time in EDT (UT-4 hrs)','fontweight','bold','FontSize',11);
% ylabel('MISA SNR at 330 km','fontweight','bold','FontSize',11);
% 
% subplot(3,1,2);
% plot(MISA_time,MISA_SNR_bp(:,51));
% xlabel('Time in EDT (UT-4 hrs)','fontweight','bold','FontSize',11);
% set(gca, 'fontweight','bold', 'FontSize',11);
% ylabel('MISA Filtered SNR at 330 km','fontweight','bold','FontSize',11);
% 
% subplot(3,1,3);
% plot(T_MISA_SNR,P_MISA_SNR,'b');
% xlim([0 plot_Freq_Lim]);
% hold on;
% % plot(T_SNR_08apr08_10_or,P_SNR_08apr08_10_or,'r');
% % hold on;
% plot(T_MISA_SNR,gf_MISA_SNR,'LineWidth',1.5,'color',[0.5 0.5 0.5]);
% xlim([0 plot_Freq_Lim]);  
% xlabel('Periods (min)','fontweight','bold','FontSize',11);
% ylabel('Normalized Power','fontweight','bold','FontSize',11);
% set(gca, 'fontweight','bold', 'FontSize',11);
% 
% 
% 
% figure5=figure;
% subplot(2,1,1);
% pcolor(Zenith_time,Zenith_height,Zenith_vel');
% shading interp;
% lighting phong;
% colorbar;
% caxis([-100 100]);
% colormap('jet');
% ylim([100 600]);
% % ylim([90 500]);
% % xlim([0 70]);
% 
% ylabel('Range (km)','fontweight','bold','FontSize',12);
% tstr1 = strcat(sprintf('%s', date_str),' (Zenith LOS Velocity)');
% title(tstr1,'fontweight','bold','FontSize',12);
% set(gca, 'fontweight','bold', 'FontSize',12);
% 
% subplot(2,1,2);
% pcolor(MISA_time,MISA_height,MISA_vel');
% shading interp;
% lighting phong;
% colorbar;
% caxis([-100 100]);
% ylim([100 600]);
% xlabel('Time in EDT (UT-4 hrs)','fontweight','bold','FontSize',12);
% ylabel('Range (km)','fontweight','bold','FontSize',12);
% tstr1 = strcat(sprintf('%s', date_str),' (MISA LOS Velocity)');
% title(tstr1,'fontweight','bold','FontSize',12);
% set(gca, 'fontweight','bold', 'FontSize',12);
% 
% 
% 
% figure6=figure;
% subplot(2,1,1);
% pcolor(Zenith_time,Zenith_height,Zenith_dopplar');
% shading interp;
% lighting phong;
% colorbar;
% caxis([-100 100]);
% colormap('jet');
% ylim([100 600]);
% % ylim([90 500]);
% % xlim([0 70]);
% 
% ylabel('Range (km)','fontweight','bold','FontSize',12);
% tstr1 = strcat(sprintf('%s', date_str),' (Zenith LOS Dopplar)');
% title(tstr1,'fontweight','bold','FontSize',12);
% set(gca, 'fontweight','bold', 'FontSize',12);
% 
% subplot(2,1,2);
% pcolor(MISA_time,MISA_height,MISA_dopplar');
% shading interp;
% lighting phong;
% colorbar;
% caxis([-100 100]);
% ylim([100 600]);
% xlabel('Time in EDT (UT-4 hrs)','fontweight','bold','FontSize',12);
% ylabel('Range (km)','fontweight','bold','FontSize',12);
% tstr1 = strcat(sprintf('%s', date_str),' (MISA LOS Dopplar)');
% title(tstr1,'fontweight','bold','FontSize',12);
% set(gca, 'fontweight','bold', 'FontSize',12);
% 
% 
% figure;
% plot(Zenith_T_max,Zenith_height,'*');
% ylim([90 700]);xlim([0 plot_Freq_Lim]);
% xlabel('Dominant Period (min)');
% ylabel('Height (km)');
% 
% 
% figure;
% plot(Zenith_angle_max,Zenith_height,'*');
% ylim([90 700]);xlim([0 plot_Freq_Lim]);
% xlabel('Dominant Period Phase (degree)');
% ylabel('Height (km)');
% 
% % figure;
% % subplot(3,1,1);
% % plot(time,SNR_ring(47,:));
% % % hold on;
% % % plot(time,SNR_sg(47,:),'r');
% % set(gca, 'fontweight','bold', 'FontSize',11);
% % xlabel('Time in EDT (UT-4 hrs)','fontweight','bold','FontSize',11);
% % ylabel('Modelled Ionosphere at 300 km','fontweight','bold','FontSize',11);
% % 
% % subplot(3,1,2);
% % plot(time,SNR_ring_bp(47,:));
% % xlabel('Time in EDT (UT-4 hrs)','fontweight','bold','FontSize',11);
% % set(gca, 'fontweight','bold', 'FontSize',11);
% % ylabel('Residuals','fontweight','bold','FontSize',11);
% % 
% % subplot(3,1,3);
% % plot(T_SNR_ring,P_SNR_ring,'b');
% % xlim([0 plot_Freq_Lim]);
% % hold on;
% % plot(T_SNR_ring,gf_ring,'LineWidth',1.5,'color',[0.5 0.5 0.5]);
% % xlim([0 plot_Freq_Lim]);  
% % xlabel('Periods (min)','fontweight','bold','FontSize',11);
% % ylabel('Normalized Power','fontweight','bold','FontSize',11);
% % set(gca, 'fontweight','bold', 'FontSize',11);


%% "Cherry Picked" ISR Data since MSTIDs seem to vary in height for Zenith

% cherryPick_time_points = [-12,-2.5,-0.5,8.8,11.43,22.9,23.8,33,36,44]; %hr, approximate time to chery pick for, corresponds to each height
% cherryPick_height_points = [240,240,300,300,240,240,300,300,240,240]; %km, approximate height to cherry pick for, corresponds to each time
cherryPick_time_points_Orig = [-12,-2,-2,9.8,9.8,23.07,23.07,35.25,35.25,44]; %hr, approximate time to chery pick for, corresponds to each height
cherryPick_height_points_Orig = [240,240,300,300,240,240,300,300,240,240]; %km, approximate height to cherry pick for, corresponds to each time
%error detection here yo
if( isequal(length(cherryPick_time_points_Orig),length(cherryPick_height_points_Orig)) ~= 1 )
    rounded = round(mean([length(cherryPick_time_points_Orig),length(cherryPick_height_points_Orig)]));
    if(length(cherryPick_time_points_Orig) ~= rounded)
        error(['cherryPick_time_points size is wrong: ',num2str(length(cherryPick_time_points_Orig)),' while average size is ',num2str(rounded),'.']);
    else
        error(['cherryPick_height_points size is wrong: ',num2str(length(cherryPick_height_points_Orig)),' while average size is ',num2str(rounded),'.']);
    end
end


cherryPick_time_points = cherryPick_time_points_Orig; %Zenith prep adjusting
cherryPick_height_points = cherryPick_height_points_Orig; %Zenith prep adjusting

for(i = 1:length(cherryPick_time_points))
    cherryPick_time_points(i) = Zenith_time( min(abs(cherryPick_time_points(i) - Zenith_time)) == abs(cherryPick_time_points(i) - Zenith_time) );
    %adjust times to be closest to the real data time points
    
    cherryPick_height_points(i) = Zenith_height( min(abs(cherryPick_height_points(i) - Zenith_height)) == abs(cherryPick_height_points(i) - Zenith_height) );
    %adjust heights to be closest to the real data time points
end

figure;
subplot(2,1,1);
pcolor(Zenith_time,Zenith_height,Zenith_SNR');
shading interp;
lighting phong;
colorbar;
% caxis([-plot_SNR_Lim plot_SNR_Lim]);
colormap('gray');
ylabel('Height (km)');
xlabel('Time in UT - 0 Hr at 0 UT May 7th');
title('Zenith SNR pre-noise-removal');
Xaxisvar = [-12,-8,-4,0, 4, 8, 12, 16, 20, 24 28, 32, 36, 40, 44, 48, 52, 56, 60, 64, 68, 72]; %shows off important magnitude points
set(gca,'xtick',Xaxisvar);
hold on;
plot(cherryPick_time_points,cherryPick_height_points,'r');
set(gca, 'fontweight','bold', 'FontSize',18);

for(i = 1:size(cherryPick_time_points,2) )
    if(i ~= size(cherryPick_time_points,2) && (cherryPick_time_points(i) == cherryPick_time_points(i+1)) )
        cherryPick_time_points(i) =  Zenith_time( find(min(abs(cherryPick_time_points(i) - Zenith_time)) == abs(cherryPick_time_points(i) - Zenith_time)) - 1 );
        %fix two heights at the same time - causes issues for interpolating
    end
end
cherryPick_height_expanded = interp1(cherryPick_time_points,cherryPick_height_points,Zenith_time); %interpolate the height between the requested values to match size of Zenith_time

cherryPick_Zenith_SNR_bp = zeros(size(cherryPick_height_expanded,1),1); %preallocate
cherryPick_Zenith_SNR = zeros(size(cherryPick_height_expanded,1),1); %preallocate
for(i = 1:size(cherryPick_height_expanded,1) )
    
    cherryPick_Zenith_SNR(i) = Zenith_SNR(i, min(abs(cherryPick_height_expanded(i) - Zenith_height)) == abs(cherryPick_height_expanded(i) - Zenith_height) );
    %pick out Zenith_SNR values that correspond to the desired height
    
    cherryPick_Zenith_SNR_bp(i) = Zenith_SNR_bp(i, min(abs(cherryPick_height_expanded(i) - Zenith_height)) == abs(cherryPick_height_expanded(i) - Zenith_height) );
    %pick out Zenith_SNR_bp values that correspond to the desired height
end


%% "Cherry Picked" ISR Data since MSTIDs seem to vary in height for Zenith

cherryPick_time_points = cherryPick_time_points_Orig; %restore for MISA adjusting
cherryPick_height_points = cherryPick_height_points_Orig; %restore for MISA adjusting

for(i = 1:length(cherryPick_time_points))
    cherryPick_time_points(i) = MISA_time( min(abs(cherryPick_time_points(i) - MISA_time)) == abs(cherryPick_time_points(i) - MISA_time) );
    %adjust times to be closest to the real data time points
    
    cherryPick_height_points(i) = MISA_height( min(abs(cherryPick_height_points(i) - MISA_height)) == abs(cherryPick_height_points(i) - MISA_height) );
    %adjust heights to be closest to the real data time points
end

subplot(2,1,2);
pcolor(MISA_time,MISA_height,MISA_SNR');
shading interp;
lighting phong;
colorbar;
% caxis([-plot_SNR_Lim plot_SNR_Lim]);
colormap('gray');
ylabel('Height (km)');
xlabel('Time in UT - 0 Hr at 0 UT May 7th');
title('MISA SNR pre-noise-removal');
Xaxisvar = [-12,-8,-4,0, 4, 8, 12, 16, 20, 24 28, 32, 36, 40, 44, 48, 52, 56, 60, 64, 68, 72]; %shows off important magnitude points
set(gca,'xtick',Xaxisvar);
hold on;
plot(cherryPick_time_points,cherryPick_height_points,'r');
set(gca, 'fontweight','bold', 'FontSize',18);

for(i = 1:size(cherryPick_time_points,2) )
    if(i ~= size(cherryPick_time_points,2) && (cherryPick_time_points(i) == cherryPick_time_points(i+1)) )
        cherryPick_time_points(i) =  MISA_time( find(min(abs(cherryPick_time_points(i) - MISA_time)) == abs(cherryPick_time_points(i) - MISA_time)) - 1 );
        %fix two heights at the same time - causes issues for interpolating
    end
end
cherryPick_height_expanded = interp1(cherryPick_time_points,cherryPick_height_points,MISA_time); %interpolate the height between the requested values to match size of MISA_time

cherryPick_MISA_SNR_bp = zeros(size(cherryPick_height_expanded,1),1); %preallocate
cherryPick_MISA_SNR = zeros(size(cherryPick_height_expanded,1),1); %preallocate
for(i = 1:size(cherryPick_height_expanded,1) )
    
    cherryPick_MISA_SNR(i) = MISA_SNR(i, min(abs(cherryPick_height_expanded(i) - MISA_height)) == abs(cherryPick_height_expanded(i) - MISA_height) );
    %pick out Zenith_SNR values that correspond to the desired height
    
    cherryPick_MISA_SNR_bp(i) = MISA_SNR_bp(i, min(abs(cherryPick_height_expanded(i) - MISA_height)) == abs(cherryPick_height_expanded(i) - MISA_height) );
    %pick out Zenith_SNR_bp values that correspond to the desired height
end


figure;
subplot(2,1,1);
pcolor(Zenith_time,Zenith_height,Zenith_SNR_bp');
shading interp;
lighting phong;
colorbar;
caxis([-plot_SNR_Lim plot_SNR_Lim]);
colormap('gray');
ylabel('Height (km)');
xlabel('Time in UT - 0 Hr at 0 UT May 7th');
title(['Zenith SNR High-passed - ',num2str(round(1/bs)),' Hr Cutoff Period']);
Xaxisvar = [-12,-8,-4,0, 4, 8, 12, 16, 20, 24 28, 32, 36, 40, 44, 48, 52, 56, 60, 64, 68, 72]; %shows off important magnitude points
set(gca,'xtick',Xaxisvar);
hold on;
plot(cherryPick_time_points,cherryPick_height_points,'r');
set(gca, 'fontweight','bold', 'FontSize',18);

subplot(2,1,2);
pcolor(MISA_time,MISA_height,MISA_SNR_bp');
shading interp;
lighting phong;
colorbar;
caxis([-plot_SNR_Lim plot_SNR_Lim]);
colormap('gray');
ylabel('Height (km)');
xlabel('Time in UT - 0 Hr at 0 UT May 7th');
title(['MISA SNR High-passed - ',num2str(round(1/bs)),' Hr Cutoff Period']);
Xaxisvar = [-12,-8,-4,0, 4, 8, 12, 16, 20, 24 28, 32, 36, 40, 44, 48, 52, 56, 60, 64, 68, 72]; %shows off important magnitude points
set(gca,'xtick',Xaxisvar);
hold on;
plot(cherryPick_time_points,cherryPick_height_points,'r');
set(gca, 'fontweight','bold', 'FontSize',18);


%% ==== Spectral Estimation of Filtered Zenith SNR at VARIABLE altitude ===
% 47 corresponds to 300 km
% 25 corresponds to 200 km
% 3 corresponds to 100 km
% 70 corresponds to 400 km
% 92 corresponds to 500 km (all approximate)

% Spectral estimation of filtered Zenith SNR at variable km....
Zenith_data=[Zenith_time, cherryPick_Zenith_SNR_bp];

%Uses Lombscargle instead of Fourier as it is not necessarially periodic
Ross_Lombscargle_optimized(Zenith_data); 
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


% ZENITH 300 KM PLOT
figure;
subplot(3,2,1);
plot(Zenith_time,cherryPick_Zenith_SNR);
axis([ -inf inf -inf inf]);
% hold on;
% plot(time_08apr08_10,SNR_08apr08_10_sg(:,50),'r');
xlabel('Time in EDT (UT-4 hrs)');
ylabel('Zenith SNR at variable km');
title('Zenith Cherry Picked');
set(gca, 'fontweight','bold', 'FontSize',18);


subplot(3,2,3);
plot(Zenith_time,cherryPick_Zenith_SNR_bp);
axis([ -inf inf -inf inf]);
xlabel('Time in EDT (UT-4 hrs)');
ylabel('Zenith Filtered SNR at variable km');
set(gca, 'fontweight','bold', 'FontSize',18);

subplot(3,2,5);
plot(T_Zenith_SNR,P_Zenith_SNR,'b');
xlim([0 plot_Freq_Lim]);
hold on;
% plot(T_SNR_08apr08_10_or,P_SNR_08apr08_10_or,'r');
% hold on;
plot(T_Zenith_SNR,gf_Zenith_SNR,'LineWidth',1.5,'color',[0.5 0.5 0.5]);
xlim([0 plot_Freq_Lim]);  
xlabel('Periods (min)');
ylabel('Normalized Power');
set(gca, 'fontweight','bold', 'FontSize',18);


%% ==== Spectral Estimation of Filtered MISA SNR at VARIABLE altitude ===
% 47 corresponds to 300 km
% 25 corresponds to 200 km
% 3 corresponds to 100 km
% 70 corresponds to 400 km
% 92 corresponds to 500 km (all approximate)

% Spectral estimation of filtered Zenith SNR at variable km....
MISA_data=[MISA_time, cherryPick_MISA_SNR_bp];

%Uses Lombscargle instead of Fourier as it is not necessarially periodic
Ross_Lombscargle_optimized(MISA_data); 
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


% MISA 300 KM PLOT
subplot(3,2,2);
plot(MISA_time,cherryPick_MISA_SNR);
axis([ -inf inf -inf inf]);
% hold on;
% plot(time_08apr08_10,SNR_08apr08_10_sg(:,50),'r');
xlabel('Time in EDT (UT-4 hrs)');
ylabel('MISA SNR at variable km');
title('MISA Cherry Picked');
set(gca, 'fontweight','bold', 'FontSize',18);

subplot(3,2,4);
plot(MISA_time,cherryPick_MISA_SNR_bp);
axis([ -inf inf -inf inf]);
xlabel('Time in EDT (UT-4 hrs)');
ylabel('MISA Filtered SNR at variable km');
set(gca, 'fontweight','bold', 'FontSize',18);

subplot(3,2,6);
plot(T_Zenith_SNR,P_Zenith_SNR,'b');
xlim([0 plot_Freq_Lim]);
hold on;
% plot(T_SNR_08apr08_10_or,P_SNR_08apr08_10_or,'r');
% hold on;
plot(T_Zenith_SNR,gf_Zenith_SNR,'LineWidth',1.5,'color',[0.5 0.5 0.5]);
xlim([0 plot_Freq_Lim]);  
xlabel('Periods (min)');
ylabel('Normalized Power');
set(gca, 'fontweight','bold', 'FontSize',18);


%% =======Bonus Walking Lomb-Scargle for Cherry Picked Altitudes===========
if(FLG_Walking_CherryPick_Zenith == 1) %turn on walking scargle
    
Zenith_timeDelta = Zenith_time(162)-Zenith_time(161); %hr
Zenith_hrChunkTime = walkingScargleHour; %hr, time to split the chunks

Zenith_countTo4hr = floor(Zenith_hrChunkTime/Zenith_timeDelta); %4 is hour, this is a count

Zenith_time_size = size(Zenith_time,1); %since I use it a lot, var it
Zenith_4hrChunks = floor(Zenith_time_size/Zenith_countTo4hr); %# 4 hr chunks in dataset


for(i = 1:Zenith_4hrChunks)

    if(i ~= Zenith_4hrChunks)
        Ross_Lombscargle_optimized(Zenith_data(((i-1)*Zenith_countTo4hr+1):(i*Zenith_countTo4hr),:));
    else
        Ross_Lombscargle_optimized(Zenith_data(((i-1)*Zenith_countTo4hr+1):(Zenith_time_size),:));
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
    
    figure;
    plot(T_Zenith_SNR,P_Zenith_SNR,'b');
    xlim([0 plot_Freq_Lim]);
    hold on;
    % plot(T_SNR_08apr08_10_or,P_SNR_08apr08_10_or,'r');
    % hold on;
    plot(T_Zenith_SNR,gf_Zenith_SNR,'LineWidth',1.5,'color',[0.5 0.5 0.5]);
    xlim([0 plot_Freq_Lim]);  
    xlabel('Periods (min)','fontweight','bold','FontSize',11);
    ylabel('Normalized Power','fontweight','bold','FontSize',11);
    if(i ~= Zenith_4hrChunks)
        tstr1 = strcat(sprintf('%.2f', (i-1)*Zenith_hrChunkTime+Zenith_time(1) ),' to');
        %Zenith_time(1) adds the initial offset since this doesn't start at
        %0
        tstr1 = strcat(tstr1,sprintf(' %.2f - Cherry Picked Lomb-Scargle', (i)*Zenith_hrChunkTime+Zenith_time(1)));
        title(tstr1,'fontweight','bold','FontSize',12);
    else
        tstr1 = strcat(sprintf('%.2f', (i-1)*Zenith_hrChunkTime+Zenith_time(1)),' to');
        tstr1 = strcat(tstr1,sprintf(' %.2f - Cherry Picked Lomb-Scargle', Zenith_time(end)));
        title(tstr1,'fontweight','bold','FontSize',12);
    end
    
    set(gca, 'fontweight','bold', 'FontSize',11);
end

end

%% Zenith and MISA Velocity

% Filter, optional - does nothing...
% n=42; % order of the Hamming window used in the custom function
% % c = 3.32*pi; %Hamming constant
% % M = n/2;
% % bandwidth = c/M;
% %The above was included to investigate the bandwidth - 1/2. Related to
% %bs/lf?
% fp=bs; % stop band stoppin freq (or pass band passing freq, depending on how you look at it)
% f= 1/Zenith_timeDelta; % the sampling frequency, based off of the time delta calc'd
% wp=2*fp/f; % Normalizing the frequencies (Matlab is 0 to 1)
% %Calculation of filter coefficients
% % [b,a]=fir1(n,wp,'high'); %This takes the order and lower cut off
% %Uses the default hamming window
% % frequency ORIG
% W = hann(n+1);
% % [b,a] = fir1(n,[wp,wl],W); %option for band-pass (30 min to 2 hr)
% [b,a]=fir1(n,wp,'high',W); %just high-pass (2 hr and lower OK)
% Zenith_vel(:,:) = filtfilt(b,a,Zenith_vel(:,:)); %Appies the filter


figure;
subplot(2,1,1);
pcolor(Zenith_time,Zenith_height,Zenith_vel');
shading interp;
lighting phong;
colorbar;
caxis([-100 100]);
colormap('gray');
ylim([90 700]);
% xlabel('Time in EDT (UT-4 hrs)');
xlabel('Time in UT - 0 Hr at 0 UT May 7th');
ylabel('Height (km)');
% tstr1 = strcat(sprintf('%s', date_str),' (Zenith Velocity)');
title('Zenith LOS Ion Velocity');
set(gca, 'fontweight','bold', 'FontSize',18);
% Xaxisvar = [-12, -8, -4, 0, 4, 8, 12, 16, 20, 24 28, 32, 36, 40, 44, 48]; %shows off important magnitude points
% Xaxisvar = [-12,-8,-4,0, 4, 8, 12, 16, 20, 24 28, 32, 36, 40, 44, 48, 52, 56, 60, 64, 68, 72]; %shows off important magnitude points
Xaxisvar = [-12:2:44]; %important times
set(gca,'xtick',Xaxisvar);
set(gca, 'fontweight','bold', 'FontSize',18);



subplot(2,1,2);
pcolor(MISA_time,MISA_height,MISA_vel');
shading interp;
lighting phong;
colorbar;
colormap('gray');
caxis([-100 100]);
ylim([90 700]);
% xlabel('Time in EDT (UT-4 hrs)','fontweight','bold','FontSize',12);
xlabel('Time in UT - 0 Hr at 0 UT May 7th');
ylabel('Height (km)');
% tstr1 = strcat(sprintf('%s', date_str),' (MISA Velocity)');
title('MISA LOS Ion Velocity');
set(gca, 'fontweight','bold', 'FontSize',18);
% Xaxisvar = [-12, -8, -4, 0, 4, 8, 12, 16, 20, 24 28, 32, 36, 40, 44, 48]; %shows off important magnitude points
% Xaxisvar = [-12,-8,-4,0, 4, 8, 12, 16, 20, 24 28, 32, 36, 40, 44, 48, 52, 56, 60, 64, 68, 72]; %shows off important magnitude points
Xaxisvar = [-12:2:44]; %important times
set(gca,'xtick',Xaxisvar);
set(gca, 'fontweight','bold', 'FontSize',18);


figure;
subplot(2,1,1);
plot(Zenith_time,Zenith_vel(:,filtHeightZ));
xlabel('Time in UT - 0 Hr at 0 UT May 7th');
ylabel('Zenith LOS Ion Velocity at 300 km');
title('LOS Ion Velocity');
set(gca, 'fontweight','bold', 'FontSize',18);
set(gca,'xtick',Xaxisvar);
xlim([-inf inf]);


subplot(2,1,2);
plot(MISA_time,MISA_vel(:,filtHeightM));
xlabel('Time in UT - 0 Hr at 0 UT May 7th');
ylabel('MISA LOS Ion Velocity at 300 km');
title('LOS Ion Velocity');
set(gca, 'fontweight','bold', 'FontSize',18);
set(gca,'xtick',Xaxisvar);
xlim([-inf inf]);


vel_Lim = 20;
%======== Zenith Vel no filt=======
totalcount = size(Zenith_vel,1); % # of samples taken
deltatime = Zenith_timeDelta; %sample interval in time - units hours *)
deltafreq = 1/(totalcount*deltatime); %frequency resolution per above sampling interval in time
frequencysamples = ((0:(totalcount-1))*deltafreq)';
% deltafreq = 1/((2^nextpow2(totalcount)*FFT_coeff)*deltatime);
% frequencysamples = ((0:2^nextpow2(totalcount)*FFT_coeff-1)*deltafreq)';
% frequencysamples = linspace(0,ceil((totalcount-1)*deltafreq),totalcount)';
% spectrummax = frequencysamples(totalcount/2 - 1);
% timedomain = zeros(totalcount,1);

frequencydomain = fft(Zenith_vel); %,2^nextpow2(totalcount)*FFT_coeff
powerspectrum = sqrt(abs(frequencydomain.*conj(frequencydomain)))./totalcount; %properly normalized now

figure;
subplot(2,1,1);
% plot(1./frequencysamples*60,powerspectrum)
pcolor(1./frequencysamples*60,Zenith_height,powerspectrum')
shading interp;
lighting phong;
colorbar;
caxis([0 vel_Lim]);
title('DFT of Zenith Vel');
xlabel('Period (min)');
ylabel('Height (km)');
% axis([0, plot_Freq_Lim, -inf, inf])
ylim([90 700]);
xlim([0, plot_Freq_Lim]);
set(gca, 'fontweight','bold', 'FontSize',18);
% axis 'auto y'

%======== MISA Vel no filt=======
totalcount = size(MISA_vel,1); % # of samples taken
deltatime = MISA_timeDelta; %sample interval in time - units hours *)
deltafreq = 1/(totalcount*deltatime); %frequency resolution per above sampling interval in time
frequencysamples = ((0:(totalcount-1))*deltafreq)';
% deltafreq = 1/((2^nextpow2(totalcount)*FFT_coeff)*deltatime);
% frequencysamples = ((0:2^nextpow2(totalcount)*FFT_coeff-1)*deltafreq)';
% frequencysamples = linspace(0,ceil((totalcount-1)*deltafreq),totalcount)';
% spectrummax = frequencysamples(totalcount/2 - 1);
% timedomain = zeros(totalcount,1);

frequencydomain = fft(MISA_vel); %,2^nextpow2(totalcount)*FFT_coeff
powerspectrum = sqrt(abs(frequencydomain.*conj(frequencydomain)))./totalcount; %properly normalized now

subplot(2,1,2);
% plot(1./frequencysamples*60,powerspectrum)
pcolor(1./frequencysamples*60,MISA_height,powerspectrum')
shading interp;
lighting phong;
colorbar;
caxis([0 vel_Lim]);
title('DFT of MISA Vel');
xlabel('Period (min)');
ylabel('Height (km)');
axis([0, plot_Freq_Lim, 90, 700])
set(gca, 'fontweight','bold', 'FontSize',18);
% axis 'auto y'


%% ====================Optional Time Cutout Settings================
% timeCutout = [18, 35]; %choose UT time cutout range
% %REALLY BREAKS STUFF SORRY
% 
% Zenith_timeCutout_minMax = [ find(min(abs(Zenith_time - min(timeCutout))) == abs(Zenith_time - min(timeCutout))) ...
%      find(min(abs(Zenith_time - max(timeCutout))) == abs(Zenith_time - max(timeCutout))) ]; %get the time points for Zenith that match
% MISA_timeCutout_minMax = [ find(min(abs(MISA_time - min(timeCutout))) == abs(MISA_time - min(timeCutout))) ...
%      find(min(abs(MISA_time - max(timeCutout))) == abs(MISA_time - max(timeCutout))) ]; %get the time points for Zenith that match
% 
% Zenith_SNR = Zenith_SNR(min(Zenith_timeCutout_minMax):max(Zenith_timeCutout_minMax),:); %tonedown Zenith_SNR
% Zenith_SNR_bp = Zenith_SNR_bp(min(Zenith_timeCutout_minMax):max(Zenith_timeCutout_minMax),:); %tonedown Zenith_SNR_bp
% Zenith_time = Zenith_time(min(Zenith_timeCutout_minMax):max(Zenith_timeCutout_minMax),:); %tonedown Zenith_time
% 
% MISA_SNR = MISA_SNR(min(MISA_timeCutout_minMax):max(MISA_timeCutout_minMax),:); %tonedown MISA_SNR
% MISA_SNR_bp = MISA_SNR_bp(min(MISA_timeCutout_minMax):max(MISA_timeCutout_minMax),:); %tonedown MISA_SNR_bp
% MISA_time = MISA_time(min(MISA_timeCutout_minMax):max(MISA_timeCutout_minMax),:); %tonedown MISA_time


%% Zenith and MISA Scargle ALL ALTITUDES

if( FLG_Scargle_AllAlt == 1 )
%=========== Zenith SNR no filt ================

for(i = 1:length(Zenith_height) )
    Ross_Lombscargle_optimized([ Zenith_time , Zenith_SNR(:,i) ]);
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
    if( i == 1 )
        powerspectrum = zeros(length(P_Zenith_SNR),length(Zenith_height)); %preallocate
    end
    
    powerspectrum(:,i) = P_Zenith_SNR; %record the power spectrum
end
    
figure;
subplot(2,2,1);
pcolor(T_Zenith_SNR,Zenith_height,powerspectrum')
shading interp;
lighting phong;
colorbar;
caxis([0 Zenith_Scargle_Lim]);
title('Periodogram of Zenith SNR');
xlabel('Period (min)');
ylabel('Height (km)');
axis([0, plot_Freq_Lim, 90, 700])
set(gca, 'fontweight','bold', 'FontSize',18);
% axis 'auto y'

%=========== MISA SNR no filt ================

for(i = 1:length(MISA_height) )
    Ross_Lombscargle_optimized([ MISA_time , MISA_SNR(:,i) ]);
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
    if( i == 1 )
        powerspectrum = zeros(length(P_Zenith_SNR),length(MISA_height)); %preallocate
    end
    
    powerspectrum(:,i) = P_Zenith_SNR; %record the power spectrum
end
    
subplot(2,2,3);
pcolor(T_Zenith_SNR,MISA_height,powerspectrum')
shading interp;
lighting phong;
colorbar;
caxis([0 MISA_Scargle_Lim]);
title('Periodogram of MISA SNR');
xlabel('Period (min)');
ylabel('Height (km)');
axis([0, plot_Freq_Lim, 90, 700])
set(gca, 'fontweight','bold', 'FontSize',18);
% axis 'auto y'

%=========== Zenith SNR no filt ================

for(i = 1:length(Zenith_height) )
    Ross_Lombscargle_optimized([ Zenith_time , Zenith_SNR_bp(:,i) ]);
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
    if( i == 1 )
        powerspectrum = zeros(length(P_Zenith_SNR),length(Zenith_height)); %preallocate
    end
    
    powerspectrum(:,i) = P_Zenith_SNR; %record the power spectrum
end
    
subplot(2,2,2);
pcolor(T_Zenith_SNR,Zenith_height,powerspectrum')
shading interp;
lighting phong;
colorbar;
caxis([0 Zenith_Scargle_Lim]);
title('Periodogram of Zenith SNR High-passed');
xlabel('Period (min)');
ylabel('Height (km)');
axis([0, plot_Freq_Lim, 90, 700])
set(gca, 'fontweight','bold', 'FontSize',18);
% axis 'auto y'

%=========== MISA SNR BP ================

for(i = 1:length(MISA_height) )
    Ross_Lombscargle_optimized([ MISA_time , MISA_SNR_bp(:,i) ]);
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
    if( i == 1 )
        powerspectrum = zeros(length(P_Zenith_SNR),length(MISA_height)); %preallocate
    end
    
    powerspectrum(:,i) = P_Zenith_SNR; %record the power spectrum
end
    
subplot(2,2,4);
pcolor(T_Zenith_SNR,MISA_height,powerspectrum')
shading interp;
lighting phong;
colorbar;
caxis([0 MISA_Scargle_Lim]);
title('Periodogram of MISA SNR High-passed');
xlabel('Period (min)');
ylabel('Height (km)');
axis([0, plot_Freq_Lim, 90, 700])
set(gca, 'fontweight','bold', 'FontSize',18);
% axis 'auto y'

end


%% Zenith and MISA FFT ALL ALTITUDES
% cutOut_time = [11,23]; %hr, time range to cut out
% cutOut_time_min_Zenith = find( min(abs(Zenith_time - min(cutOut_time))) == abs(Zenith_time - min(cutOut_time))); %get the min/max index closest to the time req
% cutOut_time_max_Zenith = find( min(abs(Zenith_time - max(cutOut_time))) == abs(Zenith_time - max(cutOut_time)));
% cutOut_time_min_MISA = find( min(abs(MISA_time - min(cutOut_time))) == abs(MISA_time - min(cutOut_time)));
% cutOut_time_max_MISA = find( min(abs(MISA_time - max(cutOut_time))) == abs(MISA_time - max(cutOut_time)));
% 
% Zenith_SNR = Zenith_SNR(cutOut_time_min_Zenith:cutOut_time_max_Zenith,:); %cut out directly yolo
% Zenith_SNR_bp = Zenith_SNR_bp(cutOut_time_min_Zenith:cutOut_time_max_Zenith,:); %cut out directly yolo
% MISA_SNR = MISA_SNR(cutOut_time_min_MISA:cutOut_time_max_MISA,:); %cut out directly yolo
% MISA_SNR_bp = MISA_SNR_bp(cutOut_time_min_MISA:cutOut_time_max_MISA,:); %cut out directly yolo

%======== Zenith SNR no filt=======
totalcount = size(Zenith_SNR,1); % # of samples taken
deltatime = Zenith_timeDelta; %sample interval in time - units hours *)
deltafreq = 1/(totalcount*deltatime); %frequency resolution per above sampling interval in time
frequencysamples = ((0:(totalcount-1))*deltafreq)';
% deltafreq = 1/((2^nextpow2(totalcount)*FFT_coeff)*deltatime);
% frequencysamples = ((0:2^nextpow2(totalcount)*FFT_coeff-1)*deltafreq)';
% frequencysamples = linspace(0,ceil((totalcount-1)*deltafreq),totalcount)';
% spectrummax = frequencysamples(totalcount/2 - 1);
% timedomain = zeros(totalcount,1);

frequencydomain = fft(Zenith_SNR); %,2^nextpow2(totalcount)*FFT_coeff
powerspectrum = sqrt(abs(frequencydomain.*conj(frequencydomain)))./totalcount; %properly normalized now

figure;
subplot(2,2,1);
% plot(1./frequencysamples*60,powerspectrum)
pcolor(1./frequencysamples*60,Zenith_height,powerspectrum')
shading interp;
lighting phong;
colorbar;
caxis([0 Zenith_FFT_Lim]);
title('DFT of Zenith SNR');
xlabel('Period (min)');
ylabel('Height (km)');
% axis([0, plot_Freq_Lim, -inf, inf])
ylim([90 700]);
xlim([0, plot_Freq_Lim]);
set(gca, 'fontweight','bold', 'FontSize',18);
% axis 'auto y'

%======== MISA SNR no filt=======
totalcount = size(MISA_SNR,1); % # of samples taken
deltatime = MISA_timeDelta; %sample interval in time - units hours *)
deltafreq = 1/(totalcount*deltatime); %frequency resolution per above sampling interval in time
frequencysamples = ((0:(totalcount-1))*deltafreq)';
% deltafreq = 1/((2^nextpow2(totalcount)*FFT_coeff)*deltatime);
% frequencysamples = ((0:2^nextpow2(totalcount)*FFT_coeff-1)*deltafreq)';
% frequencysamples = linspace(0,ceil((totalcount-1)*deltafreq),totalcount)';
% spectrummax = frequencysamples(totalcount/2 - 1);
% timedomain = zeros(totalcount,1);

frequencydomain = fft(MISA_SNR); %,2^nextpow2(totalcount)*FFT_coeff
powerspectrum = sqrt(abs(frequencydomain.*conj(frequencydomain)))./totalcount; %properly normalized now

subplot(2,2,3);
% plot(1./frequencysamples*60,powerspectrum)
pcolor(1./frequencysamples*60,MISA_height,powerspectrum')
shading interp;
lighting phong;
colorbar;
caxis([0 MISA_FFT_Lim]);
title('DFT of MISA SNR');
xlabel('Period (min)');
ylabel('Height (km)');
axis([0, plot_Freq_Lim, 90, 700])
set(gca, 'fontweight','bold', 'FontSize',18);
% axis 'auto y'

%========ZENITH BP=======
totalcount = size(Zenith_SNR_bp,1); % # of samples taken
deltatime = Zenith_timeDelta; %sample interval in time - units hours *)
deltafreq = 1/(totalcount*deltatime); %frequency resolution per above sampling interval in time
frequencysamples = ((0:(totalcount-1))*deltafreq)';
% deltafreq = 1/((2^nextpow2(totalcount)*FFT_coeff)*deltatime);
% frequencysamples = ((0:2^nextpow2(totalcount)*FFT_coeff-1)*deltafreq)';
% frequencysamples = linspace(0,ceil((totalcount-1)*deltafreq),totalcount)';
% spectrummax = frequencysamples(totalcount/2 - 1);
% timedomain = zeros(totalcount,1);

frequencydomain = fft(Zenith_SNR_bp); %,2^nextpow2(totalcount)*FFT_coeff
powerspectrum = sqrt(abs(frequencydomain.*conj(frequencydomain)))./totalcount; %properly normalized now

subplot(2,2,2);
% plot(1./frequencysamples*60,powerspectrum)
pcolor(1./frequencysamples*60,Zenith_height,powerspectrum')
shading interp;
lighting phong;
colorbar;
caxis([0 Zenith_FFT_Lim]);
title('DFT of Zenith SNR High-passed');
xlabel('Period (min)');
ylabel('Height (km)');
axis([0, plot_Freq_Lim, 90, 700])
set(gca, 'fontweight','bold', 'FontSize',18);
% axis 'auto y'

%======= MISA BP =============
totalcount = size(MISA_SNR_bp,1); % # of samples taken
deltatime = MISA_timeDelta; %sample interval in time - units hours *)
deltafreq = 1/(totalcount*deltatime); %frequency resolution per above sampling interval in time
frequencysamples = ((0:(totalcount-1))*deltafreq)';
% deltafreq = 1/((2^nextpow2(totalcount)*FFT_coeff)*deltatime);
% frequencysamples = ((0:2^nextpow2(totalcount)*FFT_coeff-1)*deltafreq)';
%frequencysamples = linspace(0,ceil((totalcount-1)*deltafreq),totalcount)';
% spectrummax = frequencysamples(totalcount/2 - 1);
% timedomain = zeros(totalcount,1);

frequencydomain = fft(MISA_SNR_bp); %,2^nextpow2(totalcount)*FFT_coeff
powerspectrum = sqrt(abs(frequencydomain.*conj(frequencydomain)))./totalcount; %properly normalized now

subplot(2,2,4);
pcolor(1./frequencysamples*60,MISA_height,powerspectrum')
shading interp;
lighting phong;
colorbar;
caxis([0 MISA_FFT_Lim]);
ylim([215 700]);
title('DFT of MISA SNR High-passed');
xlabel('Period (min)');
ylabel('Height (km)');
axis([0, plot_Freq_Lim, 90, 700])
set(gca, 'fontweight','bold', 'FontSize',18);
% axis 'auto y'