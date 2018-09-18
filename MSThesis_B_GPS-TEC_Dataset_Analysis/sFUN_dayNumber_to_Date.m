function date = sFUN_dayNumber_to_Date(dateGiven) 
%Order to input is 1st is year, 2nd is month, 3rd is day. Example below:
% dateGiven{1,1} = '2013'; %year in UTC (UTC really not a req btw)
% dateGiven{1,2} = '127'; %day number in UTC

if( iscellstr(dateGiven) == 1 ) %checks for cells
    dateGiven = str2double(dateGiven); %converts from string to double (1 is year, 2 day #)
    
    if( (size(dateGiven,1) == 1) && (size(dateGiven,2) == 2) ) %checks for correct orientation
        
    elseif( (size(dateGiven,1) == 2) && (size(dateGiven,2) == 1) ) %checks for incorrect orientation, corrects
        dateGiven = dateGiven';
        
    else %errors out if size is not correct and uncorrectable
        error(['IN: ',mfilename,' - Numbers provided but not correct format of 2 values with 1 holding year, 2 holding date #.']);
    end
    
elseif( isnumeric(dateGiven) == 1) %checks for integers (really numbers I think but whatev)
    if( (size(dateGiven,1) == 1) && (size(dateGiven,2) == 2) ) %checks for correct orientation
        
    elseif( (size(dateGiven,1) == 2) && (size(dateGiven,2) == 1) ) %checks for incorrect orientation, corrects
        dateGiven = dateGiven';
        
    else %errors out if size is not correct and uncorrectable
        error(['IN: ',mfilename,' - Numbers provided but not correct format of 2 values with 1 holding year, 2 holding date #.']);
    end
    
else %errors out because not stuff I can deal with
    dateGivenType = whos('dateGiven'); %gets that data class
    error(['IN: ',mfilename,' - Unsupported data class provided. Provided: ',dateGivenType.class,'. Must be a cell containing integers or integers, both with 2 values.']);
end


%% Prep

MonthDays = [31; 28; 31; 30; 31; 30; 31; 31; 30; 31; 30; 31]; %preps number of days in a month
MonthDays_Leap = [31; 29; 31; 30; 31; 30; 31; 31; 30; 31; 30; 31]; %preps number of days in a month



%% Leap Year Detection
if( mod(dateGiven(1,1),4) == 0) %leap year
    %% Leap Year Skipped Detected - next will be 2100
    if((mod(dateGiven(1,1),100) == 0) && (mod(dateGiven(1,1),400) ~= 0))
        %Same alg as no leap year used here
        
        tempDayCntr = dateGiven(2); %get the current day number
        tempMonCntr = 0; %counts the months
        while( tempDayCntr > 0 && tempMonCntr <= 12 ) %makes sure we don't go to 13 months - quits when days left are 0 or negative which is what happens with correct inputs
            tempMonCntr = tempMonCntr + 1; %increment month used
            tempDayCntr = tempDayCntr - MonthDays(tempMonCntr); %days, remove days in month from year
        end
        if( tempDayCntr > 0 ) %check for errors
            error(['Day count greater than 0 and days left in year: ',num2str(tempDayCntr),' days left. Please check that input of ',num2str(dateGiven(1)),' year and ',num2str(dateGiven(2)),' day number.']);
        end
        tempDayCntr = tempDayCntr + MonthDays(tempMonCntr); %days, calculate the day in the month we found
        
    else
    %% Leap Year Confirmed (2000,2004,2008,2012,2016,2020...)
        
        tempDayCntr = dateGiven(2); %get the current day number
        tempMonCntr = 0; %counts the months
        while( tempDayCntr > 0 && tempMonCntr <= 12 ) %makes sure we don't go to 13 months - quits when days left are 0 or negative which is what happens with correct inputs
            tempMonCntr = tempMonCntr + 1; %increment month used
            tempDayCntr = tempDayCntr - MonthDays_Leap(tempMonCntr); %days, remove days in month from year
        end
        if( tempDayCntr > 0 ) %check for errors
            error(['Day count greater than 0 and days left in year: ',num2str(tempDayCntr),' days left. Please check that input of ',num2str(dateGiven(1)),' year and ',num2str(dateGiven(2)),' day number.']);
        end
        tempDayCntr = tempDayCntr + MonthDays_Leap(tempMonCntr); %days, calculate the day in the month we found
    
    end
%% No leap year detected
else %no leap year if this
    
    tempDayCntr = dateGiven(2); %get the current day number
    tempMonCntr = 0; %counts the months
    while( tempDayCntr > 0 && tempMonCntr <= 12 ) %makes sure we don't go to 13 months - quits when days left are 0 or negative which is what happens with correct inputs
        tempMonCntr = tempMonCntr + 1; %increment month used
        tempDayCntr = tempDayCntr - MonthDays(tempMonCntr); %days, remove days in month from year
    end
    if( tempDayCntr > 0 ) %check for errors
        error(['Day count greater than 0 and days left in year: ',num2str(tempDayCntr),' days left. Please check that input of ',num2str(dateGiven(1)),' year and ',num2str(dateGiven(2)),' day number.']);
    end
    tempDayCntr = tempDayCntr + MonthDays(tempMonCntr); %days, calculate the day in the month we found
    
end

date = [dateGiven(1),tempMonCntr,tempDayCntr]; %pack up the date and ship out

end