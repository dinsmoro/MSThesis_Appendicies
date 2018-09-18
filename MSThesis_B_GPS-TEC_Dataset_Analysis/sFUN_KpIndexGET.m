function [dateRange_Full,Kp_output] = sFUN_KpIndexGET(dateRange,monthOutputOption)
%GOAL: Get 3 hour Kp indexes and stuff
%expecting: [year,day#min;year,day#max] numerical format. So, say [2013,126;2013,128] for 2013 Day 126 to Day 128. No chars pls
%expecting: monthOutputOption is either 0 or 1. 1 for month form, 0 for day form. 0 day form is default.
%NOTE: Kp is every 3 hours, 0-3, 3-6, 6-9, 9-12, 12-15, 15-18, 18-21, 21-24 UT


% FOR TESTING
% dateRange = [2013,126;2013,128]; %date range 
% FOR TESTING
siteURL_base = 'ftp://ftp.gfz-potsdam.de/pub/home/obs/kp-ap/tab/kp'; %base for Kp index, will tack on dates needed

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
for( yarr = 1:length(dateRange_yearRange))
    temp = unique( dateRange_Full_Month( (dateRange_yearRange(yarr,1) == dateRange_Full_Month(:,1)), 2 ) ); %get number of unique months
    if( exist('dateRange_yearRange_numberOfMonths','var') == 0 )
        dateRange_yearRange_numberOfMonths = [repmat(dateRange_yearRange(yarr,1),length(temp),1),temp]; %create the vector
    else
        dateRange_yearRange_numberOfMonths = [dateRange_yearRange_numberOfMonths;[repmat(dateRange_yearRange(yarr,1),length(temp),1),temp]]; %create the vector
    end
end



Kp_output = zeros(size(dateRange_Full_Month,1),8); %preallocate (8 kp's per day)
cntr = 1; %preallocate
for( i = 1:size(dateRange_yearRange_numberOfMonths,1) )
    siteURL_prep = num2str(dateRange_yearRange_numberOfMonths(i,1)); %get the year
    siteURL = [siteURL_base,siteURL_prep(3:end)]; %put in the last 2 year digits
    temp_YearMonthBase = siteURL_prep(3:end); %create a temp var that holds the year and month info
    siteURL_prep = num2str(dateRange_yearRange_numberOfMonths(i,2)); %get the month
    if( length(siteURL_prep) == 1 ) %if 1, tack on a 0
        siteURL_prep = ['0',siteURL_prep]; %tack it on
    elseif( length(siteURL_prep) > 2 ) %shouldn't happen
        error(['Month length greater than 2 chars, reported month is ',siteURL_prep]);
    end
    siteURL = [siteURL,siteURL_prep,'.tab']; %put in the 2 month digits (padded w/ 0 if needed)
    temp_YearMonthBase = [temp_YearMonthBase,siteURL_prep];
    Kp_data = urlread(siteURL); %download the data needed
    
    numDaysInMonth = dateRange_Full_Month( (dateRange_yearRange_numberOfMonths(i,1) == dateRange_Full_Month(:,1)) & (dateRange_yearRange_numberOfMonths(i,2) == dateRange_Full_Month(:,2)), 3 ); %long but figures out what days are in the month and year desired
    
    for(j = 1:length(numDaysInMonth)) %run through each day
        temp = num2str(numDaysInMonth(j)); %get the day
        if( length(temp) == 1 ) %if 1, tack on a 0
            temp = ['0',temp]; %tack it on
        elseif( length(temp) > 2 ) %shouldn't happen
            error(['Day length greater than 2 chars, reported day is ',temp]);
        end
        temp_YearMonthDay = [temp_YearMonthBase,temp]; %create the Year/Month/Day string
        locale = strfind(Kp_data,temp_YearMonthDay); %get the locale data starts
        Kp_data_temp = strsplit(Kp_data(1,locale+8:locale+31)); %get the Kp data for the day in question
        for( k = 1:8)
            temp = Kp_data_temp{1,k}; %get the string data
            if( temp(2) == 'o' )
                Kp_output(cntr,k) = str2double(temp(1)); %write in the value read
            elseif( temp(2) == '-' )
                Kp_output(cntr,k) = str2double(temp(1)) - 1/3; %write in the value read
            elseif( temp(2) == '+' )
                Kp_output(cntr,k) = str2double(temp(1)) + 1/3; %write in the value read
            end
        end
        cntr = cntr + 1; %increment counter to next day line
        
    end

end

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
