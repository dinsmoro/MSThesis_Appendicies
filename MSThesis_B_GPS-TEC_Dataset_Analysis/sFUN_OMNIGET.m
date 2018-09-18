function [dateRange_Full,OMNI_output] = sFUN_OMNIGET(dateRange,monthOutputOption)
%GOAL: Get OMNI data
%expecting: [year,day#min;year,day#max] numerical format. So, say [2013,126;2013,128] for 2013 Day 126 to Day 128. No chars pls
%expecting: monthOutputOption is either 0 or 1. 1 for month form, 0 for day form. 0 day form is default.
%NOTE: OMNI is every 1 minute.
%NOTE: Reads OMNI data, which is 525600x46 (minutes in a year x 46 parameters reported, including 4 time parameters)
%NOTE: ftp://spdf.gsfc.nasa.gov/pub/data/omni/high_res_omni/hroformat.txt describes OMNI data setup

% FOR TESTING
% clear variables
% dateRange = [2013,126;2013,128]; %date range 
% FOR TESTING
ftp_client = ftp('spdf.gsfc.nasa.gov'); %ftp site (totes derpy, urlread can't deal with a big file)
cd(ftp_client, '/pub/data/omni/high_res_omni/'); %folder of ftp site desired
OMNI_columns_wanted = [1,2,3,4,17,22,28,38,42]; %columns wanted
% 1 - yr
% 2 - day#
% 3 - hr
% 4 - min
% 17 - Bz GSE (nT)
% 22 - Flow speed (km/s)
% 28 - Flow pressure (nPa)
% 38 - AE Index (nT)
% 42 - SYM/H Index (nT)



dateRange_Full = dateRange(2,1) - dateRange(1,1); %yr difference
if( dateRange_Full == 0 )
    dateRange_Full = dateRange(2,2) - dateRange(1,2); %day difference
    dateRange_Full = [repmat(dateRange(1,1),dateRange_Full+1,1),(dateRange(1,2):1:dateRange(2,2))']; %create full date range from min/max days
    dateRange_yearRange = dateRange(1,1); %set the year range as one year
else
%    error('So didn''t put in multi-year support yet sorry my b reign it in lad&/orlass'); %I tried
   dateRange_yearRange = (dateRange(1,1):1:dateRange(2,1))'; %get the full year range from min to max
   for(i = 1:length(dateRange_yearRange))
       %Leap Year Detection
        if( mod(dateRange_yearRange(i,1),4) == 0) %leap year
            %LEAP YEAR! Possibly.
            if((mod(dateRange_yearRange(i,1),100) == 0) && (mod(dateRange_yearRange(i,1),400) ~= 0))
                %NO LEAP YEAR
                %Leap Year Skipped Detected - next will be 2100
                if( length(dateRange_Full) < 2 ) %see if date range is being used a temp var or not
                    dateRange_Full = [repmat(dateRange_yearRange(i,1),365,1),(1:1:365)']; %create the date range if not it yet
                else
                    dateRange_Full = [dateRange_Full;[repmat(dateRange_yearRange(i,1),365,1),(1:1:365)']]; %if exists, tack on
                end
            else
                %Leap Year Confirmed (2000,2004,2008,2012,2016,2020...)
                if( length(dateRange_Full) < 2 ) %see if date range is being used a temp var or not
                    dateRange_Full = [repmat(dateRange_yearRange(i,1),366,1),(1:1:366)']; %create the date range if not it yet
                else
                    dateRange_Full = [dateRange_Full;[repmat(dateRange_yearRange(i,1),366,1),(1:1:366)']]; %if exists, tack on
                end
            end
        else %no leap year if this
            %no leap year
            if( length(dateRange_Full) < 2 ) %see if date range is being used a temp var or not
                dateRange_Full = [repmat(dateRange_yearRange(i,1),365,1),(1:1:365)']; %create the date range if not it yet
            else
                dateRange_Full = [dateRange_Full;[repmat(dateRange_yearRange(i,1),365,1),(1:1:365)']]; %if exists, tack on
            end
        end 
   end %end for loop for per year
   
   dateRange_Full_min = find(dateRange_Full(:,1) == dateRange(1,1) & dateRange_Full(:,2) == dateRange(1,2)); %get the min day to start on
   dateRange_Full_max = find(dateRange_Full(:,1) == dateRange(2,1) & dateRange_Full(:,2) == dateRange(2,2)); %get the max day to start on
   dateRange_Full = dateRange_Full(dateRange_Full_min:dateRange_Full_max,:); %cut out the extra
end

dateRange_Full_Month = zeros(size(dateRange_Full,1),3); %preallocate
for( i = 1:size(dateRange_Full,1)) %convert from year/day# format to year/month/day format
    dateRange_Full_Month(i,:) = sFUN_dayNumber_to_Date(dateRange_Full(i,:)); %get out the year/month/day from year/day# format
end
% for( yarr = 1:length(dateRange_yearRange))
%     temp = unique( dateRange_Full_Month( (dateRange_yearRange(yarr,1) == dateRange_Full_Month(:,1)), 2 ) ); %get number of unique months
%     if( exist('dateRange_yearRange_numberOfMonths','var') == 0 )
%         dateRange_yearRange_numberOfMonths = [repmat(dateRange_yearRange(yarr,1),length(temp),1),temp]; %create the vector
%     else
%         dateRange_yearRange_numberOfMonths = [dateRange_yearRange_numberOfMonths;[repmat(dateRange_yearRange(yarr,1),length(temp),1),temp]]; %create the vector
%     end
% end


%OMNI data cadence is 1 per minute
%download all that OMNI data
OMNI_output_raw = zeros(size(dateRange_Full,1)*1440,46); %preallocate, 46 total entries for the OMNI data
cntr = 1; %preallocate
for( i = 1:size(dateRange_yearRange,1) ) %steps through each year
    siteURL = ['omni_min',num2str(dateRange_yearRange(i,1)),'.asc']; %put in the 2 month digits (padded w/ 0 if needed)
    
    mget(ftp_client,siteURL); %download the right file
    OMNI_data = dlmread(siteURL); %read the file in
    
    dateRange_Full_CurrentYear = dateRange_Full(:,1) == dateRange_yearRange(i,1) ; %get the bits that are for the current year
    
    dateRange_Full_AllDays = ismember(OMNI_data(:,2),dateRange_Full(dateRange_Full_CurrentYear,2)); %get the days of the current year needed
    
    OMNI_output_raw(cntr:cntr+(sum(dateRange_Full_CurrentYear)*1440)-1,:) = OMNI_data(dateRange_Full_AllDays,:); %copy in all relevant data
    
    cntr = cntr + sum(dateRange_Full_CurrentYear)*1440; %increment the counter
    delete(siteURL); %delete the file when done
    
end

%Weed out bad numbers and such (denoted by varying amounts of 9's)
%Columns 7, 8, and 9 seem to be solid markers for data being bad - not sure for sure but doin it
OMNI_output_badData = OMNI_output_raw(:,7) == 999 | OMNI_output_raw(:,8) == 999 | OMNI_output_raw(:,9) == 999; %find the bad data, it's denoted by 999
OMNI_output_raw = OMNI_output_raw(~OMNI_output_badData,:); %remove the bad data
OMNI_output = OMNI_output_raw(:,OMNI_columns_wanted); %get the data desired only, call it a day

close(ftp_client); %close out the ftp when done

%final output format
if( exist('monthOutputOption','var') == 0 )
    %if doesn't exist, don't do anything
else
    %if does exist, check what it is
    if( monthOutputOption == 0 )
        %if 0, don't do anything
    else
        %if 1, output in month format
        dateRange_Full = dateRange_Full_Month; %just overwrite the output with the month form
    end
end


end
