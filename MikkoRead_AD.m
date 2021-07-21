function [ Lightcurve, Npoints ] = MikkoRead_v2(filename,NUflag,magflag)
% Read in lightcurve file called filename
%
% - NUflag=1 -- Assumes 8-column format with uncertainties set to 0.01 by default
% - NUflag=0 -- Assumes 9-column format with uncertainties in column 9
%
% - magflag=0 -- Assumes intensities in Column 2 
% - magflag=1 -- Assumes magnitudes in Column 2
%
% Default (omitting the flags) assumes 9-column format with intensities in
% column 2.

% Slight change to original: also outputs number of points in each
% lightcurve
% For use in AD scripts reading in many, single point LCs


fid = fopen(filename,'r');
line = fgets(fid);
NoOfCurves = sscanf(line,'%d');
i=1;
index = 1;
Lightcurve = [];
while i<=NoOfCurves
    
    line=fgets(fid);
    NoOfPoints = sscanf(line,'%d %d');
    Npoints(i) = NoOfPoints(1,1);
    j=1;
    
    while j<= NoOfPoints(1)
        
        line = fgets(fid);
        
        if exist('NUflag','var') && (NUflag==1)
            temp(1:8) = sscanf(line,'%f %f %f %f %f %f %f %f');
            LCP(j,11) = 0.01;
        else 
            temp(1:9) = sscanf(line,'%f %f %f %f %f %f %f %f %f');
            LCP(j,11) = temp(9);
        end    

        LCP(j,1:8) = temp(1:8);
        LCP(j,9) = i;
        
        if exist('magflag','var') && (magflag==1)
            LCP(j,10) = LCP(j,2);
        else
            LCP(j,10) = -2.5*log10(LCP(j,2))+5;
        end
        
        j=j+1;
        index = index+1;
    end
    i=i+1;
    % Center magnitudes on 0
    LCP(:,10) = LCP(:,10) - mean(LCP(:,10));
    Lightcurve = [Lightcurve;LCP];
    clear LCP;

end

Npoints = Npoints'

end



