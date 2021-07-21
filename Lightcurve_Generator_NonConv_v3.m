function [results] = Lightcurve_Generator_NonConv_v3(Fn,FNA,shadows,LCSfileVariable,T0,lambda,beta,P,omega,B0,gF,hwidth,rough,x,yorp,dataset,lightcurveNo,phi,bestscalingFixed,ylims,dir)
% Lightcurve generator for producing publication quality lightcurves 
% for non-convex models

% Line 824- edited to stop fix plot
%
% Input patameters:
% - Fn -- matrix of facet normals
% - FNA -- matrix of facet areas
% - shadows -- container structure for the pre-computed shadowing geomertry
% information;
% - LCSfileVariable - the 10-column matrix holding lightcurves loaded to
% memory
% - T0 -- model epoch -- time at wchich the X-axis crosses the
% plane-of-sky, given in [days]
% - lambda, beta -- ecliptic coordinates of the rotation pole, given in 
% [degrees];
% - P -- siderial rotation period at T0, given in [h];
% - omega --  Single-particle scattering albedo for Hapke scattering model,
% it is used as c, or the mixing parameter for the combination of
% Lommel-Seelinger and Lambertian scattering;
% - B0 -- Opposition surge amplitude (used for Hapke only);
% - gF -- Asymmetry parameter of the single-particle phase function (used 
% for Hapke only);
% - hwidth -- Opposition surge width (used for Hapke only);
% - rough -- Macroscopic roughness(used for Hapke only);
% - x -- flag used to indicate which scattering model to use: 
%   - 3 for Lambertian;
%   - 4 for Lommel-Seelinger (this is the most used option);
%   - 5 for Hapke;
%   - 8 for a combination of Lommel-Seelinger + omega * Lambertian; 
% - yorp -- the YORP spin-up factor in [rad/d^2];
% - dataset -- the fist part of output file name;
% - lightcurveNo -- ID of the lightcurve to be plotted, use 0 to plot all
% the lightcurves
% - phi -- rotational phase offset to be added to lightcurves
% - bestscalingFixed -- a vertical offset to be added to the lightcurve
% points
% - ylims -- plotting limits for the "_fix" plot
%
% Output:
% * results container sctructure with fields:
% - phasedCurve: the matrix has size 
% (total number of lightcurve points + 102 * number of lightcurves) * 8 --
% the artificial lightcurve matrix with brightness 
% calculated for rotational phases corresponding to all the points in the
% observed lightcure + 102 evenly spaced points across a full rotation, the
% rows for each lightcurve are sorted by rotational phase. Each row has 8
% columns: 
% 1. t-T0 - time since T0 (in days);
% 2. plotphase - rotation phase as a fraction of full rotation (accounting 
% for YORP) ;
% 3. SL_illumination - magnitude as calculated using Lambertian scattering, centered around 0 
% 4. SLS_illumination - magnitude as calculated using Lommel-Seelinger
% scattering, centered around 0
% 5. scatteringmodel - magnitude as calculated using Hapke scattering,
% centered around 0
% (nonsense if x!=8 as Hapke scattering parameters might be simply set to 0)
% 6. phase - rotation phase as a fraction of full rotation, but with yorp=0
% and using starting period
% 7. lightcurveNo 
% 8. SLLS_illumination - magnitude as calulated using combination of L and
% LS, centered around 0
% 
% - angle2sun: matrix size (number of facets) x (total number od lightcurve points)
% the cosine of the angle between each facet normal and the direction to the Sun 
% for each lightcurve point 
% - angle2earth: matrix size (number of facets) x (total number od lightcurve points)
% the cosine of the angle between each facet normal and the direction to
% the observer for each lightcurve point
% - Lambert: matrix size (number of facets) x (total number od lightcurve points)
% contribution to the lightcurve from each facet, calulated assuming Lambertian light scattering 
% - LommelSeelinger: matrix size (number of facets) x (total number od lightcurve points)
% contribution to the lightcurve from each facet, calulated assuming Lommel-Seelinger light scattering
% - Hapke: matrix size (number of facets) x (total number od lightcurve points)
% contribution to the lightcurve from each facet, calulated assuming Hapke light scattering 
% - observer: matrix size (total number of lightcurve points) x  3 
% the asteroid-Earth vector in asteroid-centric coordinates (viewing
% vector)
% - sun: matrix size (total number of lightcurve points) x  3
% the asteroid-Sun vector in asteroid-centric coordinates (illumination vector) 
% - bestscaling: matrix size (number of lightcurves) x 1 
% the difference between mean level of artificial lightcurve points for the
% selected scattering function and the  mean level of the observed point in 
% magnitudes, in other words the vertical shift needed to fit the points 
% onto the lightcurve 
% - datapoints:  matrix size (total number of lightcurve points) * 8 
% similar to phasedCurve, but only for rotational phases corresponding to
% the observed data points and only the magnitudes corresponding to the
% selected light scattering function are centered around 0
% - MikkoFile: matrix size (total number of lightcurve points) x 12
% the LCSFileVariable, but expanded by column 12 - rotational phase
% (same as plotphase in phasedCurve), and with the magnitude in column 10
% corrected by lightcurve-relevant bestscaling factor 
% - phaseangle: matrix size (number of lightcurves) x 1 - the phase angle 
% in degrees, the angle between positions of Earth and Sun as seen from the
% target
% - aspectAngle: matrix size (number of lightcurves) x 1 - the aspect
% angle, the angle between the pole orientation and the direction to Earth
% 
% * two plots with names like:
%   - 'pictures/',dataset,'_',num2str(lightcurveNo,'%02d') for the auto-scaled
% Y-axis
%   - 'pictures/'dataset'_'num2str(lightcurveNo,'%02d')'_fix' for y-axis
%   with upper scale limit -ylimits and lower limit ylimits; 
% 
% Output plots include:
% - data points (red) 
% - model lightcurve plotted with the selected scattering function (black line)

CBblue=[68 119 170]/255;
CBred=[238 102 119]/255;
CBgreen=[34 136 51]/255;
CBgrey=[187 187 187]/255;


lcinit=lightcurveNo;

if lightcurveNo > 0
    LCSfileVariable = LCSfileVariable(LCSfileVariable(:,9)==lightcurveNo,:);
    LCSfileVariable(:,9) = 1;
    lightcurveNo = 1;
else
    lightcurveNo = 1;
end
%convert lambda and beta to radians
lambda = lambda * pi /180;
beta = (90-beta)*pi /180;


% and phi
phi = phi*pi/180;


%convert P to days...!
P = P / 24;

%new variable for facet normals.  Not really necessary but is useful if
%this function is ever implemented as a script.  Will prevent overwriting
%orginal facet normals.
FNr = Fn;

%Original z axis.  Used for PAB determination.

%Set up the transformation matrices for lambda and beta.
Rlambda = [cos(lambda) -sin(lambda) 0;sin(lambda) cos(lambda) 0; 0 0 1];
Rlambda = Rlambda/ norm(Rlambda);


Rneglambda = [cos(lambda) sin(lambda) 0; -sin(lambda) cos(lambda) 0; 0 0 1];
Rneglambda = Rneglambda/ norm(Rneglambda);


Rbeta =  [cos(beta) 0 sin(beta);0 1 0;-sin(beta) 0 cos(beta) ];
Rbeta = Rbeta / norm(Rbeta);

Rnegbeta =  [cos(beta) 0 -sin(beta);0 1 0; sin(beta) 0 cos(beta) ];
Rnegbeta = Rnegbeta / norm(Rnegbeta);


%counting variables initialised
l=1;
m=0;

%Number of data points in the lightcurve variable
numberOfPoints = size(LCSfileVariable);

%Number of Lightcurves in the datafile.
numberOfLightcurves = LCSfileVariable(numberOfPoints(1),9);

%Difference between T0 and time... initialised to zero

%Facet area illuminated which can be seen by the observer
area_illuminated = []; 
%Lambertian brightness for each data point
SL_illumination=[];
%Lommel-Seeliger brightness for each datapoint
SLS_illumination = [];


%Output Return Matrix
results.phasedCurve =[];
correction = [ ];
results.angle2sun=[];
results.angle2earth=[];
results.Lambert=[];
results.LommelSeelinger=[];
results.Hapke=[];
results.observer=[];
results.sun=[];


lightcurvePoints = 1;

%Scaling between artifical and observed lightcurves
bestscaling = zeros(numberOfLightcurves,1);

phasedDatareturnMatrix = [];

%Phase angle of the observations
phaseangle = zeros(numberOfLightcurves,1);

%Aspect Angle
aspect = zeros(numberOfLightcurves,1);

times = LCSfileVariable(:,1);
%phases = (times -T0 )/P - floor((times-T0)/P);
diffs = times - T0;
alphas = 2 * pi /( P) .* (diffs) + 0.5*yorp*diffs.^2;
phases = alphas/(2*pi) - floor(alphas/(2*pi));

LCSfileVariable(:,12) = phases;

%operate from lightcurve number given in input arguments to whatever
%Lightcurve number you're interested in... this can be varied...
for lightcurveNo=lightcurveNo:numberOfLightcurves
    
    %Set up the sun and earth directional matrices.  
    
    sunDirection = mean(LCSfileVariable((LCSfileVariable(:,9)==lightcurveNo),3:5))*1e8;
    earthDirection = mean(LCSfileVariable((LCSfileVariable(:,9)==lightcurveNo),6:8))*1e8;
    
    sunDirN=norm(sunDirection);
    earthDirN=norm(earthDirection);
    
    sunDirection1=(Rnegbeta*Rneglambda*sunDirection')';
    earthDirection1=(Rnegbeta*Rneglambda*earthDirection')';
    
    %Start time for producing the artifical lightcurve.  The variable
    %lightcurvePoints is updated at the end of the loop with the numeber of
    %data points in the lightcurve...
    time0 = LCSfileVariable(lightcurvePoints,1);
    
    %Calculates the phase angle of this lightcurve observation.  Also
    %calculates the Hapke scattering functions
    phaseangle(lightcurveNo) = acos((dot(earthDirection,sunDirection))/(norm(earthDirection)*norm(sunDirection)));
    oppositionSurge = B0 / (1+ (tan(phaseangle(lightcurveNo)/2))/hwidth);
    PPF = (1-gF^2)/((1+2*gF*cos(phaseangle(lightcurveNo)) + gF^2)^1.5);
    BPPF = (1+oppositionSurge)*PPF;

    % Compute the lightcurve from 0 to 1 full phase
    time1 = time0 + P*101/100;
    
    %Sets up a matrix called return matrix
    artificalLightcurveData = 0;
    
    t = time0;
    i= 1;
    
    %Construct a full rotational lightcurve
    while t <= time1
        
        %Time difference between t and the initial time
        diff = t - T0;
                
        %Work out the plotting phases...
        wholecycles = floor((t - T0 )/(P));
        cycles = (t - T0 )/(P);
        phase = (cycles-wholecycles);
        %plotphase = phase;
        
        %Phase of the lightcurve
        alpha = 2 * pi /( P) * (diff) + 0.5*yorp*diff*diff;
        
        plotphase = alpha/(2*pi)-floor(alpha/(2*pi));
        
        
        %Calculate the rotation matrix due to rotation of object
        Rr = [cos(alpha+phi) -sin(alpha+phi) 0 ; sin(alpha+phi) cos(alpha+phi) 0 ;0 0 1];
   %     Rr = Rr / norm(Rr);
 
        Rnegr = [cos(alpha+phi) sin(alpha+phi) 0 ; -sin(alpha+phi) cos(alpha+phi) 0 ;0 0 1];
        %Rnegr = Rnegr / norm(Rnegr);

        
        %Recalculate facet normal orientations
        
        %FNr1 = (Rlambda * Rbeta *  Rr * FNr')';
        
        %Assign FNr1 to a new variable. Not strictly necessary but
        %potentially saves a headache
    
        %A = FNr1;
        A=FNr;
        
        %Create matrices which are the same size as A, but are populated
        %with the Earth and Sun direction vectors.
        
        %size(LCSfileVariable(LCSfileVariable(:,9)==lightcurveNo,:),1);
        
          sunDirection2=(Rnegr*sunDirection1')';
        earthDirection2=(Rnegr*earthDirection1')';
   
        if size(LCSfileVariable(LCSfileVariable(:,9)==lightcurveNo,:),1)>=2
            C = repmat(earthDirection2(1:3),size(A,1),1);
            D = repmat(sunDirection2(1:3),size(A,1),1);
        else
            C = earthDirection2;
            D = sunDirection2;
        end
        
        

        %Described in Section 6.1 of documentation
        %Find the angle between the earth direction and facet normals and
        %the sun direction and the facet normals.  Any angles which are
        %negative are assigned a value of zero.  Any angles which are
        %listed as nan are assigned as zeros.
        
        %E = cross(A,C)./normr(cross(A,C));
        
        E = dot(A,C,2);
        %angle2earth = cos(atan2(E(:,1),G));
        angle2earth=E/earthDirN;
        angle2earth(angle2earth<0)=0;
        angle2earth(isnan(angle2earth))=0;
        
        %E = cross(A,D)./normr(cross(A,D));
        
        G = dot(A,D,2);
%        angle2sun = cos(atan2(E(:,1),G));
        angle2sun=G/sunDirN;
        
        angle2sun(angle2sun<0)=0;
        angle2sun(isnan(angle2sun))=0;
        musmue = angle2sun.*angle2earth;

        
        
        %% Add shadowing
% denE and denS are taken before dividing by norms, as the un-normalized
% Earth and Sun directions are used
        
denE=1.0./E;
denS=1.0./G;




sizeFN=size(FNA);

    % i loops over facets (illuminated and visible)
    for ii = 1:sizeFN(1)
     if musmue(ii)==0
         continue
     end
% Loop over all facets above the i-th facet horizon,
% facing away from it but facing the Sun   
  
   rI=shadows.Dtab(:,ii).*denS;
   w = shadows.P0(:,:,ii)-rI*sunDirection2;
   
   wu=dot(w,shadows.u,2);
   wv=dot(w,shadows.v,2);
  sI= (shadows.uv.*wv-shadows.vv.*wu);      
  tI= (shadows.uv.*wu-shadows.uu.*wv);
  
% Determine whether the center of the facet in question is obstructed from 
% the Sun
   
   if sum(sI(shadows.logic(ii,:))>=0 & tI(shadows.logic(ii,:))>=0 & sI(shadows.logic(ii,:))+tI(shadows.logic(ii,:))<=1)>0 
       angle2sun(ii)=0;
       continue
   end
   
   rI=shadows.Dtab(:,ii).*denE;
   w = shadows.P0(:,:,ii)-rI*earthDirection2;
   wu=dot(w,shadows.u,2);
   wv=dot(w,shadows.v,2);
   sI= (shadows.uv.*wv-shadows.vv.*wu);      
   tI= (shadows.uv.*wu-shadows.uu.*wv);
       
% Determine whether the center of the facet in question is obstructed from
% the view
   
    if sum(sI(shadows.logic(ii,:))>=0 & tI(shadows.logic(ii,:))>=0 & sI(shadows.logic(ii,:))+tI(shadows.logic(ii,:))<=1)>0 
         angle2earth(ii)=0;
     end   
   
    end
    
       musmue = angle2sun.*angle2earth;


       %% Continue
        
        
           
        %More Hapke Terms
        H1 = (1 + 2*angle2sun)./(1+2*angle2sun*((1 - omega).^0.5));
        H2 = (1 + 2*angle2earth)./(1+2*angle2earth*((1 - omega).^0.5));
        
        %SLS parameter.
        musplusmue = angle2sun+angle2earth;
        OVER=musmue./musplusmue;  %SLS
        OVER(isnan(OVER))=0;

        %Total area illuminated
        area_illuminated = sum(FNA.*musmue);
        
        %Lommel-Selliger model
        SLS = sum((FNA.*OVER));
        
        %Hapke Scattering
        scatteringmodel = sum( ( (omega/(4*pi) ) .* FNA .* OVER .* ( BPPF + H1 .* H2 - 1) * cosd(rough) ) );
        scatteringmodel =-2.5*log10(scatteringmodel);

        SL_illumination= -2.5*log10((area_illuminated)) ;
        SLS_illumination=-2.5*log10(SLS);
        SLLS_illumination=-2.5*log10(SLS+omega*area_illuminated);
        
        
        %Data Matrix holding times, plotting phases, reflectances,
        %phases,lightcurve number and area illuminated.
 %       artificalLightcurveData(i,1:8) = [t-T0,plotphase,SL_illumination,SLS_illumination,scatteringmodel,phase,lightcurveNo,area_illuminated];
        artificalLightcurveData(i,1:8) = [t-T0,plotphase,SL_illumination,SLS_illumination,scatteringmodel,phase,lightcurveNo,SLLS_illumination];

        %Update new plotting time.
        t = t + P/100;
        
        %Update matrix indices
        i = i+1;
    end
    
    %These variables become important for scaling the artifical and
    %observed lightcurves
    stomp=i;
    data = LCSfileVariable(LCSfileVariable(:,9)==lightcurveNo,:);
    
    %How many datapoints in this lightcurve
    numberOfDataPoints = size(data(:,1));
    
    %Iterate over each data point in the lightcurve
    %Much of this loop is identical to the previous loop.  However, we're
    %iterating over the lightcurve data points now so that we can do a
    %straight compairsion to the artifical curves and teh observed data
    %points
    
     matrixForThePurposeOfScaling=[];

    
    
    for j = 1:numberOfDataPoints(1)
        
        t = data(j,1);
        diff = t - T0;
        
        wholecycles = floor((t - T0 )/(P));
        cycles = (t - T0 )/(P);
        phase = (cycles-wholecycles);
        %plotphase = phase;
        alpha = 2 * pi /( P) * (diff) + 0.5*yorp*diff*diff ;
        
        plotphase = alpha/(2*pi)-floor(alpha/(2*pi));
        
        Rr = [cos(alpha+phi) -sin(alpha+phi) 0 ; sin(alpha+phi) cos(alpha+phi) 0 ;0 0 1];
      %  Rr = Rr / norm(Rr);
        
        Rnegr = [cos(alpha+phi) sin(alpha+phi) 0 ; -sin(alpha+phi) cos(alpha+phi) 0 ;0 0 1];
       % Rnegr = Rnegr / norm(Rnegr);
        
        %FNr1 = (Rlambda * Rbeta *  Rr * FNr')';
        
        
        A=FNr;
        
        sunDirection2=(Rnegr*sunDirection1')';
        earthDirection2=(Rnegr*earthDirection1')';
        
        
        
        C = repmat(earthDirection2(1:3),size(A,1),1);
        D = repmat(sunDirection2(1:3),size(A,1),1);
  
        
        E = dot(A,C,2);
        %angle2earth = cos(atan2(E(:,1),G));
        angle2earth=E/earthDirN;
        angle2earth(angle2earth<0)=0;
        angle2earth(isnan(angle2earth))=0;
        
        %E = cross(A,D)./normr(cross(A,D));
        
        G = dot(A,D,2);
%        angle2sun = cos(atan2(E(:,1),G));
        angle2sun=G/sunDirN;
        
        angle2sun(angle2sun<0)=0;
        angle2sun(isnan(angle2sun))=0;


        musmue = angle2sun.*angle2earth;

        
          %% Add shadowing
% denE and denS are taken before dividing by norms, as the un-normalized
% Earth and Sun directions are used
        
denE=1.0./E;
denS=1.0./G;


sizeFN=size(FNA);

    % i loops over facets (illuminated and visible)
    for ii = 1:sizeFN(1)
     if musmue(ii)==0
         continue
     end
% Loop over all facets above the i-th facet horizon,
% facing away from it but facing the Sun   
  
   rI=shadows.Dtab(:,ii).*denS;
   w = shadows.P0(:,:,ii)-rI*sunDirection2;
   
   wu=dot(w,shadows.u,2);
   wv=dot(w,shadows.v,2);
  sI= (shadows.uv.*wv-shadows.vv.*wu);      
  tI= (shadows.uv.*wu-shadows.uu.*wv);
  
% Determine whether the center of the facet in question is obstructed from 
% the Sun
   
   if sum(sI(shadows.logic(ii,:))>=0 & tI(shadows.logic(ii,:))>=0 & sI(shadows.logic(ii,:))+tI(shadows.logic(ii,:))<=1)>0 
       angle2sun(ii)=0;
       continue
   end
   
   rI=shadows.Dtab(:,ii).*denE;
   w = shadows.P0(:,:,ii)-rI*earthDirection2;
   wu=dot(w,shadows.u,2);
   wv=dot(w,shadows.v,2);
   sI= (shadows.uv.*wv-shadows.vv.*wu);      
   tI= (shadows.uv.*wu-shadows.uu.*wv);
       
% Determine whether the center of the facet in question is obstructed from
% the view
   
    if sum(sI(shadows.logic(ii,:))>=0 & tI(shadows.logic(ii,:))>=0 & sI(shadows.logic(ii,:))+tI(shadows.logic(ii,:))<=1)>0 
         angle2earth(ii)=0;
     end   
   
    end
    
           musmue = angle2sun.*angle2earth;

       
         %% Continue
        
        
        H1 = (1 + 2*angle2sun)./(1+2*angle2sun*((1 - omega).^0.5));
        H2 = (1 + 2*angle2earth)./(1+2*angle2earth*((1 - omega).^0.5));
        
                
        
        musplusmue = angle2sun+angle2earth;
        OVER=musmue./musplusmue;  
        OVER(isnan(OVER))=0;

        area_illuminated = sum(FNA.*musmue);
        SLS = sum((FNA.*OVER));
        
        
        
        scatteringmodel = sum( ( (omega/(4*pi)) .*FNA.* OVER .* ( BPPF + (H1 .* H2) - 1) * cosd(rough) ) );
        scatteringmodel =-2.5*log10(scatteringmodel);

        SL_illumination= -2.5*log10((area_illuminated)) ;
        SLS_illumination=-2.5*log10(SLS);
        
        
         SLLS_illumination=-2.5*log10(SLS+omega*area_illuminated);
        
        artificalLightcurveData(i,1:8) = [t-T0,plotphase,SL_illumination,SLS_illumination,scatteringmodel,phase,lightcurveNo,SLLS_illumination];
        
%        artificalLightcurveData(i,1:8) = [t-T0,plotphase,SL_illumination,SLS_illumination,scatteringmodel,phase,lightcurveNo,area_illuminated];
        
        %This matrix holds the artifical curve produced at each datapoint.
        %It is used in scaling the lightcurves
%        matrixForThePurposeOfScaling(j,1:8) = [t-T0,plotphase,SL_illumination,SLS_illumination,scatteringmodel,phase,lightcurveNo,area_illuminated];
        matrixForThePurposeOfScaling(j,1:8) = [t-T0,plotphase,SL_illumination,SLS_illumination,scatteringmodel,phase,lightcurveNo,SLLS_illumination];
        
        i=i+1;
        
        
        results.angle2sun=[results.angle2sun angle2sun];
        results.angle2earth=[results.angle2earth  angle2earth];
        results.Lambert=[results.Lambert FNA.*musmue];
        results.LommelSeelinger=[results.LommelSeelinger  FNA.*OVER];       
        results.Hapke = [results.Hapke ( (omega/(4*pi)) .*FNA.* OVER .* ( BPPF + (H1 .* H2) - 1) * cosd(rough) )];
        results.observer=[results.observer; earthDirection2];
        results.sun=[results.sun; sunDirection2];
        
        
    end
    
  %Used for scaling purposes.  Along with stomp samples the start and end
    %of a lightcurve in terms of total ligthcurve points generated so far.
    chow = i;
    
    %Original Pole orientation of the object
    VPAB = [0 0 1];
    
    %Observed pole orientation of the object.
    poleorientation = Rlambda * Rbeta *  Rr *VPAB'
    
    %Aspect angle of the observation
    aspect(lightcurveNo) = acos((dot(earthDirection,poleorientation))/(norm(earthDirection)*norm(poleorientation)))*180/pi;
    
      
    %Now subtract the mean brightness of the model lightcurve
    %Subtract mean brightness from the observed data

    scales(lightcurveNo,1:5) = [...
       (max(artificalLightcurveData(artificalLightcurveData(:,7)==lightcurveNo,3))+min(artificalLightcurveData(artificalLightcurveData(:,7)==lightcurveNo,3)))/2,...
       (max(artificalLightcurveData(artificalLightcurveData(:,7)==lightcurveNo,4))+min(artificalLightcurveData(artificalLightcurveData(:,7)==lightcurveNo,4)))/2,...
       (max(artificalLightcurveData(artificalLightcurveData(:,7)==lightcurveNo,5))+min(artificalLightcurveData(artificalLightcurveData(:,7)==lightcurveNo,5)))/2,...
       mean(LCSfileVariable(LCSfileVariable(:,9)==lightcurveNo,10)),...
       (max(artificalLightcurveData(artificalLightcurveData(:,7)==lightcurveNo,8))+min(artificalLightcurveData(artificalLightcurveData(:,7)==lightcurveNo,8)))/2];
   
    %Now subtract the mean brightness of the model lightcurve
    % matrixForThePurposeOfScaling(matrixForThePurposeOfScaling(:,7)==lightcurveNo,4) =  matrixForThePurposeOfScaling(matrixForThePurposeOfScaling(:,7)==lightcurveNo,4)...
    %-mean(artificalLightcurveData(artificalLightcurveData(:,7)==lightcurveNo,4));
    
    
   ids = [0 0 1 2  3 0 0 5];
   
    
     matrixForThePurposeOfScaling(matrixForThePurposeOfScaling(:,7)==lightcurveNo,x) =  matrixForThePurposeOfScaling(matrixForThePurposeOfScaling(:,7)==lightcurveNo,x)...
         - scales(lightcurveNo,ids(x));
   
    
    phasedDatareturnMatrix = [phasedDatareturnMatrix;matrixForThePurposeOfScaling];
   
    %Subtract the mean brightness from each of the model lightcurves
    artificalLightcurveData(artificalLightcurveData(:,7)==lightcurveNo,3) = artificalLightcurveData(artificalLightcurveData(:,7)==lightcurveNo,3) ...
       - scales(lightcurveNo,1);
    artificalLightcurveData(artificalLightcurveData(:,7)==lightcurveNo,4) = artificalLightcurveData(artificalLightcurveData(:,7)==lightcurveNo,4) ...
       - scales(lightcurveNo,2);
    artificalLightcurveData(artificalLightcurveData(:,7)==lightcurveNo,5) = artificalLightcurveData(artificalLightcurveData(:,7)==lightcurveNo,5) ... 
       - scales(lightcurveNo,3);
   
    artificalLightcurveData(artificalLightcurveData(:,7)==lightcurveNo,8) = artificalLightcurveData(artificalLightcurveData(:,7)==lightcurveNo,8) ... 
       - scales(lightcurveNo,5);

    %Subtract mean brightness from the observed data
    LCSfileVariable(LCSfileVariable(:,9)==lightcurveNo,10) = LCSfileVariable(LCSfileVariable(:,9)==lightcurveNo,10) ... 
         - scales(lightcurveNo,4);
   
   
   
   %This allows you to enter your own scaling parameter. 
   if abs(bestscalingFixed) > 0
       bestsca = bestscalingFixed;
   else
        dataForScaling = LCSfileVariable(LCSfileVariable(:,9)==lightcurveNo,10);
        
        difAD=matrixForThePurposeOfScaling(matrixForThePurposeOfScaling(:,7)==lightcurveNo,x) - dataForScaling;
        tmpsca = mean(difAD);
        oldScalingChi= 1000000000;
        
        for sca=tmpsca-0.01:0.001:tmpsca+0.01
            dataForScaling = dataForScaling + sca;
            scalingChi = sum(( (difAD-sca )./LCSfileVariable(LCSfileVariable(:,9)==lightcurveNo,11)).^2);
            dataForScaling = dataForScaling - sca;
            if scalingChi < oldScalingChi
                oldScalingChi = scalingChi;
                bestsca = sca;
               
            end
        end

   end
   
   %Subtract the best or manual scaling factor from the observed data.
   LCSfileVariable(LCSfileVariable(:,9)==lightcurveNo,10) = LCSfileVariable(LCSfileVariable(:,9)==lightcurveNo,10) + bestsca;
   
   %Add the best scaling factor to the variable bestscaling.
   bestscaling(lightcurveNo) = bestsca;
   
   temp=artificalLightcurveData(stomp:stomp+j-1,:);
   
   %This lightcurve is added to the variable correction
   correction = [correction;temp];
   
   
   
   %Set up plot.
   plotfigure=figure('PaperSize',[20.984 20.984]);
   set(0,'defaulttextinterpreter','latex')
   
   %Order the lightcurve according to phase.
   artificalLightcurveData = sortrows(artificalLightcurveData,2);
   
   %Determine the maximum and minimum model lightcurve values.
   ymaxModel(2) = max(artificalLightcurveData(artificalLightcurveData(:,7)==lightcurveNo,4));
   yminModel(2) = min(artificalLightcurveData(artificalLightcurveData(:,7)==lightcurveNo,4));
   ymaxModel(1) = max(artificalLightcurveData(artificalLightcurveData(:,7)==lightcurveNo,3));
   yminModel(1) = min(artificalLightcurveData(artificalLightcurveData(:,7)==lightcurveNo,3));
   ymaxModel(3) = max(artificalLightcurveData(artificalLightcurveData(:,7)==lightcurveNo,5));
   yminModel(3) = min(artificalLightcurveData(artificalLightcurveData(:,7)==lightcurveNo,5));
   ymaxModel(5) = max(artificalLightcurveData(artificalLightcurveData(:,7)==lightcurveNo,8));
   yminModel(5) = min(artificalLightcurveData(artificalLightcurveData(:,7)==lightcurveNo,8));
   
   yminModel(4) = min(LCSfileVariable(LCSfileVariable(:,9)==lightcurveNo,10)-LCSfileVariable(LCSfileVariable(:,9)==lightcurveNo,11));
   ymaxModel(4) = max(LCSfileVariable(LCSfileVariable(:,9)==lightcurveNo,10)+LCSfileVariable(LCSfileVariable(:,9)==lightcurveNo,11));
   
   %Attempts to find the best yaxis scaling.  Does not work superbly well,
   %normally have to manually specify this.

   ymax=max([max(ymaxModel),-min(yminModel)])*1.25;
   
   ymin=-ymax;
   
   
   %Plot the model lightcurve

   plot(artificalLightcurveData((artificalLightcurveData(:,7)==lightcurveNo),2),(artificalLightcurveData((artificalLightcurveData(:,7)==lightcurveNo),x)),'Color','black');
   set(gca,'TickLabelInterpreter','latex')

   set(gca,'FontName','Helvetica','FontSize',16);
   %set(gca,'FontName','Helvetica')
   
    
   hold on
   

   %Plot the observed data
   
   errbar(LCSfileVariable(LCSfileVariable(:,9)==lightcurveNo,12),LCSfileVariable(LCSfileVariable(:,9)==lightcurveNo,10),LCSfileVariable(LCSfileVariable(:,9)==lightcurveNo,11),'Color', CBred);
   plot(LCSfileVariable(LCSfileVariable(:,9)==lightcurveNo,12),LCSfileVariable(LCSfileVariable(:,9)==lightcurveNo,10),'o','LineStyle','none','MarkerFaceColor',CBred, 'Color',CBred);


  % hold off
   
   %See, manually setting the yaxis limits
  
   ylim([real(ymin) real(ymax)])
   
%    ylim([-0.4 0.4])
%    ylim([-0.5 0.5])

   %How many datapoints have been used so far
   points = size((LCSfileVariable((LCSfileVariable(:,9)==lightcurveNo),1))) ;
   lightcurvePoints = lightcurvePoints + points(1) ;
   
   %Set labels
   xlabel('Rotational Phase','FontName','Helvetica','FontSize',20);
   ylabel('Relative Magnitude','FontName','Helvetica','FontSize',20);      
   startday = LCSfileVariable(LCSfileVariable(:,9)==lightcurveNo,1);
   
   dates=datefromjd(startday(1)); 

   
   %Various title parameters
   %      if lcinit > 0
   %title(strcat(num2str(lcinit),'---',num2str(years(1)),' --- ', num2str(startday(1))));
   %      else
   %title(strcat(num2str(lightcurveNo),'---',num2str(years(1)),' --- ', num2str(startday(1))));
   %      end
   % Full date title
   
      if lcinit > 0
%   title(strcat(num2str(lcinit),' \bullet ', dates ,' \bullet ', num2str(startday(1))));
  title(sprintf('%3d $\\bullet$ %14s $\\bullet$ %12.3f', lcinit, dates, startday(1)) );
          else
%   title(strcat(num2str(lightcurveNo),' $\bullet$ ', dates ,' $\bullet$ ', num2str(startday(1))));
  title(sprintf('%3d $\\bullet$ %14s $\\bullet$ %12.3f', lightcurveNo, dates, startday(1)) );
  end
        
         
         
   %title(strcat(num2str(lightcurveNo),'---',num2str(years(1))));
   %title(dataset,'FontName','Helvetica','FontSize',24)
   
  %Various annotations
%   annotation('textbox',[0.56 0.06 0.15 0.15],'String',{strcat('Aspect Angle = ',num2str(aspect(lightcurveNo),'%.1f'),'\circ')},'LineStyle','none','FontName','Helvetica','FontSize',12,'FitBoxToText','on');
%   annotation('textbox',[0.56 0.06 0.15 0.15],'String',{sprintf('Aspect Angle = %.1f\\circ',aspect(lightcurveNo))},'LineStyle','none','FontName','Helvetica','FontSize',12,'FitBoxToText','on');
   text(0.56,0.06,['Aspect Angle = $' num2str(aspect(lightcurveNo),'%.1f') '^{\circ}$'],'Interpreter','latex','FontName','Helvetica','FontSize',12,'Units','normalized')

   % xlabel('Measured Area ($m^{2}$)','Interpreter','latex')
   
   
   %annotation('textbox',[0.56 0.06 0.15 0.15],'String',{strcat('Phase Offset = ',num2str(phi,'%.1f'),'\circ')},'LineStyle','none','FontName','Helvetica','FontSize',20,'FitBoxToText','on');
%   annotation('textbox',[0.15 0.06 0.15 0.15],'String',{strcat('Phase Angle = ',num2str(phaseangle(lightcurveNo)*180/pi,'%.2f'),'\circ')},'LineStyle','none','FontName','Helvetica','FontSize',12,'FitBoxToText','on');
   text(0.06,0.06,['Phase Angle = $' num2str(phaseangle(lightcurveNo)*180/pi,'%.2f') '^\circ$'],'Interpreter','latex','FontName','Helvetica','FontSize',12,'Units','normalized')

   
%   annotation('textbox',[0.15 0.76 0.15 0.15],'String',{strcat('Model Peak-to-peak = ',num2str(ymaxModel(ids(x))-yminModel(ids(x)),'%.2f'),'^m')},'LineStyle','none','FontName','Helvetica','FontSize',12,'FitBoxToText','on');
   text(0.06,0.94,['Model Peak-to-peak = $' num2str(ymaxModel(ids(x))-yminModel(ids(x)),'%.2f') '^m$'],'Interpreter','latex','FontName','Helvetica','FontSize',12,'Units','normalized')

   %Reverse the direction of the yaxis
   set(gca,'YDir','reverse');
   
   
   %originaldata = [LCSfileVariable(LCSfileVariable(:,9)==lightcurveNo,12) LCSfileVariable(LCSfileVariable(:,9)==lightcurveNo,10) LCSfileVariable(LCSfileVariable(:,9)==lightcurveNo,11) LCSfileVariable(LCSfileVariable(:,9)==lightcurveNo,1)];
   
   %Various commands for printing or saving files.
   %filename = strcat('pictures/yorp_',num2str(lightcurveNo),'.eps')
   
     
      if lcinit > 0
          filename = strcat(dir,'pictures/',dataset,'_',num2str(lcinit,'%02d'))
      else
          filename = strcat(dir,'pictures/',dataset,'_',num2str(lightcurveNo,'%02d'))

      end
   
   
   %filename = strcat('pictures/all_PYORP-_',num2str(lightcurveNo),'.eps')
   
   %filename = strcat('pictures/NTT_lam',num2str(lambda/pi*180),'_bet',num2str(90-beta/pi*180),'.eps')
   
   format long;
   hold off
   %save('-ascii',strcat('67P_VLT_planning',num2str(lightcurveNo),'.dat'),'returnmatrix');
   %save('-ascii','-double',strcat('67P_steve_lcs_file_revised_',num2str(lightcurveNo),'.dat'),'originaldata');
   %legend('Lommel-Seeliger','Lambert-Lommel Seeliger combination','Observed','Location','SouthOutside')
   print(plotfigure,'-dpdf',filename);
   
 %  ylim([-1.25 1.25]) % JV6
   %ylim([-1 1]) % 68346
  % ylim([-0.25 0.25]) % 2102

  
    %if abs(ylims) > 0
    %      ylim([-ylims ylims])

    %   end
  
  
   %print(plotfigure,'-dpdf', strcat(filename, '_fix') );
   
   %close(plotfigure)
   
   
   
results.phasedCurve = [results.phasedCurve;artificalLightcurveData];
end

%Various output parameters.  Artifical lightcurve is returned etc.
%phasedDatareturnMatrix(1,:) = [];
%results.phasedCurve(1,:) = [];
%correction = sortrows(correction,1);
% results.correction=correction;
results.bestscaling = bestscaling;
results.datapoints = phasedDatareturnMatrix;
results.MikkoFile = LCSfileVariable;
results.phaseangle = phaseangle;
results.aspectAngle = aspect;
end