function dayNum = rinex_sFUN_dateToDayNum_Razzle(dateGiven) 
%Order to input is 1st is year, 2nd is month, 3rd is day. Example below:
% dateGiven{1,1} = '2013'; %start year in UTC (UTC really not a req btw)
% dateGiven{1,2} = '05'; %start month in UTC
% dateGiven{1,3} = '06'; %start day in UTC

if( iscellstr(dateGiven) == 1 ) %checks for cells
    dateGiven = str2double(dateGiven); %converts from string to double (1 is year, 2 is month, 3 is day)
    
    if( (size(dateGiven,1) == 1) && (size(dateGiven,2) == 3) ) %checks for correct orientation
        
    elseif( (size(dateGiven,1) == 3) && (size(dateGiven,2) == 1) ) %checks for incorrect orientation, corrects
        dateGiven = dateGiven';
        
    else %errors out if size is not correct and uncorrectable
        error(['IN: ',mfilename,' - Numbers provided but not correct format of 3 values with 1 holding year, 2 holding month, 3 holding date.']);
    end
    
elseif( isnumeric(dateGiven) == 1) %checks for integers
    if( (size(dateGiven,1) == 1) && (size(dateGiven,2) == 3) ) %checks for correct orientation
        
    elseif( (size(dateGiven,1) == 3) && (size(dateGiven,2) == 1) ) %checks for incorrect orientation, corrects
        dateGiven = dateGiven';
        
    else %errors out if size is not correct and uncorrectable
        error(['IN: ',mfilename,' - Numbers provided but not correct format of 3 values with 1 holding year, 2 holding month, 3 holding date.']);
    end
    
else %errors out because not stuff I can deal with
    dateGivenType = whos('dateGiven'); %gets that data class
    error(['IN: ',mfilename,' - Unsupported data class provided. Provided: ',dateGivenType.class,'. Must be a cell containing integers or integers, both with 3 values.']);
end



%% Leap Year Detection
if( mod(dateGiven(1,1),4) == 0) %leap year
    %% Leap Year Skipped Detected - next will be 2100
    if((mod(dateGiven(1,1),100) == 0) && (mod(dateGiven(1,1),400) ~= 0))
        %Same alg as no leap year used here
        if( dateGiven(1,2) == 1) %Jan 31
        %Nothing happens, days count up here
            tempMonth = 0; %no contribution
        elseif( dateGiven(1,2) == 2) %Feb 28
            tempMonth = 31; %jan's contribution
        elseif( dateGiven(1,2) == 3) %March 31
            tempMonth = (31+28); %prev month's contribution
        elseif( dateGiven(1,2) == 4) %April 30
            tempMonth = (31+28+31); %prev month's contribution
        elseif( dateGiven(1,2) == 5) %May 31
            tempMonth = (31+28+31+30); %prev month's contribution
        elseif( dateGiven(1,2) == 6) %June 30
            tempMonth = (31+28+31+30+31); %prev month's contribution
        elseif( dateGiven(1,2) == 7) %July 31
            tempMonth = (31+28+31+30+31+30); %prev month's contribution
        elseif( dateGiven(1,2) == 8) %August 31
            tempMonth = (31+28+31+30+31+30+31); %prev month's contribution
        elseif( dateGiven(1,2) == 9) %September 30
            tempMonth = (31+28+31+30+31+30+31+31); %prev month's contribution
        elseif( dateGiven(1,2) == 10) %Oct 31
            tempMonth = (31+28+31+30+31+30+31+31+30); %prev month's contribution
        elseif( dateGiven(1,2) == 11) %Nov 30
            tempMonth = (31+28+31+30+31+30+31+31+30+31); %prev month's contribution
        elseif( dateGiven(1,2) == 12) %Dec 31
            tempMonth = (31+28+31+30+31+30+31+31+30+31+30); %prev month's contribution
        end
    else
    %% Leap Year Confirmed (2000,2004,2008,2012,2016,2020...)
        if(dateGiven(1,2) == 1) %Jan 31
            %Nothing happens, days count up here
            tempMonth = 0; %no contribution
        elseif( dateGiven(1,2) == 2) %Feb 29 - leap year
            tempMonth = 31; %jan's contribution
        elseif( dateGiven(1,2) == 3) %March 31
            tempMonth = (31+29); %prev month's contribution
        elseif( dateGiven(1,2) == 4) %April 30
            tempMonth = (31+29+31); %prev month's contribution
        elseif( dateGiven(1,2) == 5) %May 31
            tempMonth = (31+29+31+30); %prev month's contribution
        elseif( dateGiven(1,2) == 6) %June 30
            tempMonth = (31+29+31+30+31); %prev month's contribution
        elseif( dateGiven(1,2) == 7) %July 31
            tempMonth = (31+29+31+30+31+30); %prev month's contribution
        elseif( dateGiven(1,2) == 8) %August 31
            tempMonth = (31+29+31+30+31+30+31); %prev month's contribution
        elseif( dateGiven(1,2) == 9) %September 30
            tempMonth = (31+29+31+30+31+30+31+31); %prev month's contribution
        elseif( dateGiven(1,2) == 10) %Oct 31
            tempMonth = (31+29+31+30+31+30+31+31+30); %prev month's contribution
        elseif( dateGiven(1,2) == 11) %Nov 30
            tempMonth = (31+29+31+30+31+30+31+31+30+31); %prev month's contribution
        elseif( dateGiven(1,2) == 12) %Dec 31
            tempMonth = (31+29+31+30+31+30+31+31+30+31+30); %prev month's contribution
        end
    end
%% No leap year detected
else %no leap year if this
    if( dateGiven(1,2) == 1) %Jan 31
        %Nothing happens, days count up here
        tempMonth = 0; %no contribution
    elseif( dateGiven(1,2) == 2) %Feb 28
        tempMonth = 31; %jan's contribution
    elseif( dateGiven(1,2) == 3) %March 31
        tempMonth = (31+28); %prev month's contribution
    elseif( dateGiven(1,2) == 4) %April 30
        tempMonth = (31+28+31); %prev month's contribution
    elseif( dateGiven(1,2) == 5) %May 31
        tempMonth = (31+28+31+30); %prev month's contribution
    elseif( dateGiven(1,2) == 6) %June 30
        tempMonth = (31+28+31+30+31); %prev month's contribution
    elseif( dateGiven(1,2) == 7) %July 31
        tempMonth = (31+28+31+30+31+30); %prev month's contribution
    elseif( dateGiven(1,2) == 8) %August 31
        tempMonth = (31+28+31+30+31+30+31); %prev month's contribution
    elseif( dateGiven(1,2) == 9) %September 30
        tempMonth = (31+28+31+30+31+30+31+31); %prev month's contribution
    elseif( dateGiven(1,2) == 10) %Oct 31
        tempMonth = (31+28+31+30+31+30+31+31+30); %prev month's contribution
    elseif( dateGiven(1,2) == 11) %Nov 30
        tempMonth = (31+28+31+30+31+30+31+31+30+31); %prev month's contribution
    elseif( dateGiven(1,2) == 12) %Dec 31
        tempMonth = (31+28+31+30+31+30+31+31+30+31+30); %prev month's contribution
    end
end

dayNum = tempMonth + dateGiven(1,3); %Adds the days up to the month plus the days in the month


end