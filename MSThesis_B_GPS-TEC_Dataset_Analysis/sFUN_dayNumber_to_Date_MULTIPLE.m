function [dateRange_Full_Month] = sFUN_dayNumber_to_Date_MULTIPLE(dateRange,dayNumRequest)
%GOAL: Spit out Yr/Month/Day for a given Year/Day# to Year/Day# range
%expecting: [year,day#min;year,day#max] numerical format. So, say [2013,126;2013,128] for 2013 Day 126 to Day 128. No chars pls
%expecting: dayRequest exists (pref 1 I guess whatever), then yr/day# output will occur. Otherwise normal yr/mon/day output YOU CAN CALL THIS W/O dayRequest existing yolo coding
%output: [year,month,day_min;...;year,month,day_max], should do multi-year


dateRange_Full = dateRange(2,1) - dateRange(1,1); %yr difference
if( dateRange_Full == 0 )
    dateRange_Full = dateRange(2,2) - dateRange(1,2); %day difference
    dateRange_Full = [repmat(dateRange(1,1),dateRange_Full+1,1),(dateRange(1,2):1:dateRange(2,2))']; %create full date range from min/max days
%     dateRange_yearRange = dateRange(1,1); %set the year range as one year
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

if( exist('dayNumRequest','var') == 0 ) %if dayNumRequest doesn't exist, turn this into yr/mon/day, but if it does then keep yr/day# format (bonus option!)
    dateRange_Full_Month = zeros(size(dateRange_Full,1),3); %preallocate
    for( i = 1:size(dateRange_Full,1)) %convert from year/day# format to year/month/day format
        dateRange_Full_Month(i,:) = sFUN_dayNumber_to_Date(dateRange_Full(i,:)); %get out the year/month/day from year/day# format
    end
else
    dateRange_Full_Month = dateRange_Full; %copy it over for output
end

end