function monthWord = sFUN_monthNum_to_word(monthNumber,abbrev,dotOrNot) 
%Input: a number between 1 and 12
%Output: month word
%abbrev input (can be anything) requests abbreviated month word when applicable (e.g. May don't get abbrev but January -> Jan.
%dotOrNot input (can be anything) requests no dot if it's there (e.g. Jan. -> Jan) does nothing if abbrev not activated already
%threeMax input (can be anything) requests three letters max - overrides everything else before it

if( ndims(monthNumber) ~= 2 || size(monthNumber,1) ~= 1 || size(monthNumber,2) ~= 1 )
    %make some adjustments, make it (1,1)
    monthNumber = monthNumber(1,1);
end

if( isnumeric(monthNumber) == 1 )
    %good to go
elseif( ischar(monthNumber) == 1 )
    monthNumber = str2double(monthNumber); %convert to number
elseif( iscell(monthNumber) == 1 )
    if( isnumeric(monthNumber{1,1}) == 1 )
        monthNumber = monthNumber{1,1}; %pull out the number
    elseif( ischar(monthNumber{1,1}) == 1 )
        monthNumber = str2double(monthNumber{1,1}); %pull out the character and convert to number
    else
        monthNumber
        error('ERROR in sFUN_monthNum_to_word: input cell does not hold a number or string so failure, input is printed above.');
    end
else
    monthNumber
    error('ERROR in sFUN_monthNum_to_word: input is not a number/character/cell holding num/char so failure, input is printed above.');
end

if( exist('abbrev','var') == 0 ) %if abbrev doesn't exist, output whole name
    if( monthNumber == 1 )
        monthWord = 'January';
    elseif( monthNumber == 2 )
        monthWord = 'February';
    elseif( monthNumber == 3 )
        monthWord = 'March';
    elseif( monthNumber == 4 )
        monthWord = 'April';
    elseif( monthNumber == 5 )
        monthWord = 'May';
    elseif( monthNumber == 6 )
        monthWord = 'June';
    elseif( monthNumber == 7 )
        monthWord = 'July';
    elseif( monthNumber == 8 )
        monthWord = 'August';
    elseif( monthNumber == 9 )
        monthWord = 'September';
    elseif( monthNumber == 10 )
        monthWord = 'October';
    elseif( monthNumber == 11 )
        monthWord = 'November';
    elseif( monthNumber == 12 )
        monthWord = 'December';
    else
        error(['ERROR in sFUN_monthNum_to_word: ',num2str(monthNum),' IS NOT A MONTH NUMBER WOW']);
    end
else %if abbrev did exist, output abbreviation
    if( exist('dotOrNot','var') == 0 ) %if dotOrNot doesn't exist, put on a dot!
        if( monthNumber == 1 )
            monthWord = 'Jan.';
        elseif( monthNumber == 2 )
            monthWord = 'Feb.';
        elseif( monthNumber == 3 )
            monthWord = 'Mar.';
        elseif( monthNumber == 4 )
            monthWord = 'Apr.';
        elseif( monthNumber == 5 )
            monthWord = 'May';
        elseif( monthNumber == 6 )
            monthWord = 'June';
        elseif( monthNumber == 7 )
            monthWord = 'July';
        elseif( monthNumber == 8 )
            monthWord = 'Aug.';
        elseif( monthNumber == 9 )
            monthWord = 'Sept.';
        elseif( monthNumber == 10 )
            monthWord = 'Oct.';
        elseif( monthNumber == 11 )
            monthWord = 'Nov.';
        elseif( monthNumber == 12 )
            monthWord = 'Dec.';
        else
            error(['ERROR in sFUN_monthNum_to_word: ',num2str(monthNum),' IS NOT A MONTH NUMBER WOW']);
        end
    else %if dotOrNot did exist, no dots
        if( monthNumber == 1 )
            monthWord = 'Jan';
        elseif( monthNumber == 2 )
            monthWord = 'Feb';
        elseif( monthNumber == 3 )
            monthWord = 'Mar';
        elseif( monthNumber == 4 )
            monthWord = 'Apr';
        elseif( monthNumber == 5 )
            monthWord = 'May';
        elseif( monthNumber == 6 )
            monthWord = 'June';
        elseif( monthNumber == 7 )
            monthWord = 'July';
        elseif( monthNumber == 8 )
            monthWord = 'Aug';
        elseif( monthNumber == 9 )
            monthWord = 'Sept';
        elseif( monthNumber == 10 )
            monthWord = 'Oct';
        elseif( monthNumber == 11 )
            monthWord = 'Nov';
        elseif( monthNumber == 12 )
            monthWord = 'Dec';
        else
            error(['ERROR in sFUN_monthNum_to_word: ',num2str(monthNum),' IS NOT A MONTH NUMBER WOW']);
        end
    end
end


end