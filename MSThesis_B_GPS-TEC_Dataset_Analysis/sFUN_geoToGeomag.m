function [mag_lat, mag_long] = sFUN_geoToGeomag(dateNeeded,magF,lat,long)
%Converts Geographic Coordinates to Geomagnetic Coordinates
%uses external stuff https://www.mathworks.com/matlabcentral/fileexchange/45606-igrf-magnetic-field?focused=3820319&tab=function

mag_lat = zeros(size(lat)); %preallocate
mag_long = zeros(size(long)); %preallocate
parfor( i = 1:length(lat) ) %convert all lat/long sets given
    
    if( isnan(lat(i)) == 0 || isnan(long(i)) == 0 ) %catch for NaN lat/long hitting this (from coastline conversion mostly)
        [ mag_lat(i), mag_long(i), ~ ] = magF.GEODIP( dateNeeded , lat(i) , long(i) , 1 ); %1 at end is for saying "geomagnetic coordinates"
    else
        mag_lat(i) = lat(i); %work around for NaNs that break the drawing up (coastlines)
        mag_long(i) = long(i); %work around for NaNs that break the drawing up (coastlines)
    end

end

mag_long = mag_long - 180; %deg, convert from 360 to 0 to 180 to -180 convention

end %END OF FUNCTION