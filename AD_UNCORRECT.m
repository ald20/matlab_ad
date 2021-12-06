function results = AD_UNCORRECT(root_dir,date,target,obj_filename,V,F,FN,FNA,T0,lambda,beta,P,omega,gF,B0,hwidth,rough,lcfile,lcfilename,NUflag,magflag,lcid,ylims,x,plot_artificial_lcs,pics_file)

% This function generates artificial lightcurves for a given shape model
% At the epochs listed, and calculates the relative magnitude shift 
% Required for each epoch.


%% Read in lightcurve

directory=root_dir;
[ LCS_all, npoints ] = MikkoRead_AD(lcfilename,NUflag,magflag);

% How many lightcurves are read in?
lcmax=max(LCS_all(:,9));

%% Generate an example lightcurve 

% Define artificial LC filename based on which scattering model is used
if x==3
    scatt_name = '_lam';
end
if x==4
    scatt_name = '_LS';
end
if x==5
    scatt_name = '_hapke';
end

modstr = lcfile;
% With a non-convex shape
namecore=[target '_' lcfile '_' date scatt_name ];

Fn = FN;
FNA = FNA;
shadows = ShadowingGeometry(V,F,Fn);
LCSfileVariable = LCS_all;
yorp = 0;
dataset = namecore;
lightcurveNo = lcid;
phi = 0;
ylims = ylims;
dir = root_dir;

lcinit=lightcurveNo;
 
if lightcurveNo > 0
    LCSfileVariable = LCSfileVariable(LCSfileVariable(:,9)==lightcurveNo,:);
    npoints = npoints(lightcurveNo)
    LCSfileVariable(:,9) = 1;
    lightcurveNo = 1;
else
    lightcurveNo = 1;
end

%% Convert inputs from degrees to radians
%convert lambda, beta, phi to radians
lambda = lambda * pi /180;
beta = (90-beta)*pi /180;
phi = phi*pi/180;

%convert P to days
P = P / 24;

times = LCSfileVariable(:,1);
diffs = times - T0;
alphas = 2 * pi /( P) .* (diffs) + 0.5*yorp*diffs.^2;
phases = alphas/(2*pi) - floor(alphas/(2*pi));

CBblue=[68 119 170]/255;
CBred=[238 102 119]/255;
CBgreen=[34 136 51]/255;
CBgrey=[187 187 187]/255;


%new variable for facet normals.  Not really necessary but is useful if
%this function is ever implemented as a script.  Will prevent overwriting
%orginal facet normals.
FNr = Fn;

%% Set up the transformation matrices for lambda and beta.
Rlambda = [cos(lambda) -sin(lambda) 0;sin(lambda) cos(lambda) 0; 0 0 1];
Rlambda = Rlambda/ norm(Rlambda);

Rneglambda = [cos(lambda) sin(lambda) 0; -sin(lambda) cos(lambda) 0; 0 0 1];
Rneglambda = Rneglambda/ norm(Rneglambda);

Rbeta =  [cos(beta) 0 sin(beta);0 1 0;-sin(beta) 0 cos(beta) ];
Rbeta = Rbeta / norm(Rbeta);

Rnegbeta =  [cos(beta) 0 -sin(beta);0 1 0; sin(beta) 0 cos(beta) ];
Rnegbeta = Rnegbeta / norm(Rnegbeta);

%% Initialising various counts and vectors

%counting variables initialised
l=1;
m=0;

%Number of data points in the lightcurve variable
numberOfPoints = size(LCSfileVariable)

%Number of Lightcurves in the datafile.
numberOfLightcurves = LCSfileVariable(numberOfPoints(1),9);

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

phasedDatareturnMatrix = [];

pointOnLC={}

%Phase angle of the observations
phaseangle = zeros(numberOfLightcurves,1);

%Aspect Angle
aspect = zeros(numberOfLightcurves,1);

LCSfileVariable(:,12) = phases;

tic
counter = 0;

for lightcurveNo=lightcurveNo:numberOfLightcurves
    
    %Set up the sun and earth directional matrices.  
    
    if npoints(lightcurveNo,1) > 1
        sunDirection = mean(LCSfileVariable((LCSfileVariable(:,9)==lightcurveNo),3:5))*1e8;
        earthDirection = mean(LCSfileVariable((LCSfileVariable(:,9)==lightcurveNo),6:8))*1e8;

        sunDirN=norm(sunDirection);
        earthDirN=norm(earthDirection);

        sunDirection1=(Rnegbeta*Rneglambda*sunDirection')';
        earthDirection1=(Rnegbeta*Rneglambda*earthDirection')';

    else
        sunDirection = LCSfileVariable((LCSfileVariable(:,9)==lightcurveNo),3:5)*1e8;
        earthDirection = LCSfileVariable((LCSfileVariable(:,9)==lightcurveNo),6:8)*1e8;

        sunDirN=norm(sunDirection);
        earthDirN=norm(earthDirection);
        sunDirection1=(Rnegbeta*Rneglambda*sunDirection')';
        earthDirection1=(Rnegbeta*Rneglambda*earthDirection')';
    end
    
    %sunDirection = mean(LCSfileVariable((LCSfileVariable(:,9)==lightcurveNo),3:5))*1e8;
    %earthDirection = mean(LCSfileVariable((LCSfileVariable(:,9)==lightcurveNo),6:8))*1e8;
    
    %sunDirN=norm(sunDirection);
    %earthDirN=norm(earthDirection);
    
    %sunDirection1=(Rnegbeta*Rneglambda*sunDirection')';
    %earthDirection1=(Rnegbeta*Rneglambda*earthDirection')';
    
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
    
    % Define stepsize for generating artificial LC
    stepsize=50
    % Compute the lightcurve from 0 to 1 full phase
    time1 = time0 + P*((stepsize+1)/(stepsize));
    
    %Sets up a matrix called return matrix
    artificalLightcurveData = 0;
    
    t = time0;
    i= 1;
    
    %Construct a full rotational lightcurve
    while t <= time1
        
        counter = counter+1;
        percentage_complete = round((counter/(numberOfLightcurves*(stepsize+1)))*100,4)
        %Time difference between t and the initial time
        diff = t - T0;
                
        %Work out the plotting phases...
        wholecycles = floor((t - T0 )/(P));
        cycles = (t - T0 )/(P);
        phase = (cycles-wholecycles);
        %plotphase = phase;
        
        %Phase of the lightcurve
        alpha = 2 * pi /( P) * (diff);
        
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
   
        if size(LCSfileVariable(LCSfileVariable(:,9)==lightcurveNo,:),1)>=1
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

        %% ii loops over facets (illuminated and visible)
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


       %% Sum facets to obtain brightnesses depending on scattering law used
         
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
        %scatteringmodel=0;

        SL_illumination= -2.5*log10((area_illuminated)) ;
        SLS_illumination=-2.5*log10(SLS);
        SLLS_illumination=-2.5*log10(SLS+omega*area_illuminated);
        
        %Data Matrix holding times, plotting phases, reflectances,
        %phases,lightcurve number and area illuminated.
 %       artificalLightcurveData(i,1:8) = [t-T0,plotphase,SL_illumination,SLS_illumination,scatteringmodel,phase,lightcurveNo,area_illuminated];
        artificalLightcurveData(i,1:8) = [t-T0,plotphase,SL_illumination,SLS_illumination,scatteringmodel,phase,lightcurveNo,SLLS_illumination];

        %Update new plotting time.
        t = t + P/stepsize;
        
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
    %straight comparsion to the artifical curves and thh observed data
    %points
    
     matrixForThePurposeOfScaling=[];

    
    
    for j = 1:numberOfDataPoints(1) % N points in lightcurve in question
        
        t = data(j,1);
        diff = t - T0;
        wholecycles = floor((t - T0 )/(P));
        cycles = (t - T0 )/(P);
        phase = (cycles-wholecycles);
        alpha = 2 * pi /( P) * (diff) + 0.5*yorp*diff*diff;
        plotphase = alpha/(2*pi)-floor(alpha/(2*pi));
        Rr = [cos(alpha+phi) -sin(alpha+phi) 0 ; sin(alpha+phi) cos(alpha+phi) 0 ;0 0 1];
        Rnegr = [cos(alpha+phi) sin(alpha+phi) 0 ; -sin(alpha+phi) cos(alpha+phi) 0 ;0 0 1];
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
        scatteringmodel = -2.5*log10(scatteringmodel);
        SL_illumination = -2.5*log10((area_illuminated)) ;
        SLS_illumination =-2.5*log10(SLS);
        SLLS_illumination=-2.5*log10(SLS+omega*area_illuminated);
               
        artificalLightcurveData(i,1:8) = [t-T0,plotphase,SL_illumination,SLS_illumination,scatteringmodel,phase,lightcurveNo,SLLS_illumination];%area_illuminated];
        
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
    poleorientation = Rlambda * Rbeta *  Rr *VPAB';
    
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
   
   temp=artificalLightcurveData(stomp:stomp+j-1,:);
   %This lightcurve is added to the variable correction
   correction = [correction;temp];
   
   % make list of points on the lightcurve for generating observational LC
   pointOnLC{lightcurveNo} = matrixForThePurposeOfScaling(:,x)
   
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
       
       delta_m = ymaxModel(ids(x))-yminModel(ids(x));
   
   %% Plot the model lightcurve
   if plot_artificial_lcs==1
    
       %Set up plot.
       plotfigure=figure('PaperSize',[20.984 20.984]);
       set(0,'defaulttextinterpreter','latex')
       
       %Order the lightcurve according to phase.
       artificalLightcurveData = sortrows(artificalLightcurveData,2);
      
   
       
      plot(artificalLightcurveData((artificalLightcurveData(:,7)==lightcurveNo),2),(artificalLightcurveData((artificalLightcurveData(:,7)==lightcurveNo),x)), 'Color','black');

      set(gca,'TickLabelInterpreter','latex')

       set(gca,'FontName','Helvetica','FontSize',16);
       %set(gca,'FontName','Helvetica')

       hold on
       
      % Lambertian
      plot(artificalLightcurveData((artificalLightcurveData(:,7)==lightcurveNo),2),(artificalLightcurveData((artificalLightcurveData(:,7)==lightcurveNo),3)),'LineStyle',':','Color',CBblue);
      % Lommel-Seeliger
      plot(artificalLightcurveData((artificalLightcurveData(:,7)==lightcurveNo),2),(artificalLightcurveData((artificalLightcurveData(:,7)==lightcurveNo),4)),'LineStyle','-.','Color',CBred);
      % Combined LS, Lam
      plot(artificalLightcurveData((artificalLightcurveData(:,7)==lightcurveNo),2),(artificalLightcurveData((artificalLightcurveData(:,7)==lightcurveNo),8)),'LineStyle','-.','Color',CBgreen);

      
       %Plot the observed data

       %errbar(LCSfileVariable(LCSfileVariable(:,9)==lightcurveNo,12),LCSfileVariable(LCSfileVariable(:,9)==lightcurveNo,10),LCSfileVariable(LCSfileVariable(:,9)==lightcurveNo,11),'Color', [0.6350 0.0780 0.1840]);
       %%%%%%%%%%plot(LCSfileVariable(LCSfileVariable(:,9)==lightcurveNo,12),LCSfileVariable(LCSfileVariable(:,9)==lightcurveNo,10),'x','LineStyle','none','LineWidth',1.1, 'MarkerFaceColor',CBblue, 'Color',CBblue, 'Markersize', 12);

      % plot the individual brightnesses calculated for each distinct epoch
      % should lie on the LC
       plot(matrixForThePurposeOfScaling(matrixForThePurposeOfScaling(:,7)==lightcurveNo,2),matrixForThePurposeOfScaling(matrixForThePurposeOfScaling(:,7)==lightcurveNo,x),'o','LineStyle','none','MarkerFaceColor',CBred, 'Color',CBred);
       hold off

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
              title(sprintf('%3d $\\bullet$ %14s $\\bullet$ %12.3f', lcinit, dates, startday(1)) );
          else
              title(sprintf('%3d $\\bullet$ %14s $\\bullet$ %12.3f', lightcurveNo, dates, startday(1)) );
          end
          
          
        text(0.56,0.06,['Aspect Angle = $' num2str(aspect(lightcurveNo),'%.1f') '^{\circ}$'],'Interpreter','latex','FontName','Helvetica','FontSize',12,'Units','normalized')
        text(0.56,0.12,['Phase Angle = $' num2str(phaseangle(lightcurveNo)*180/pi,'%.2f') '^\circ$'],'Interpreter','latex','FontName','Helvetica','FontSize',12,'Units','normalized')
        
        %text(0.06,0.06,['Phase Angle = $' num2str(phaseangle(lightcurveNo)*180/pi,'%.2f') '^\circ$'],'Interpreter','latex','FontName','Helvetica','FontSize',12,'Units','normalized')
        text(0.06,0.94,['Model Peak-to-peak = $' num2str(delta_m,'%.2f') '^m$'],'Interpreter','latex','FontName','Helvetica','FontSize',12,'Units','normalized')

        %Reverse the direction of the yaxis
        set(gca,'YDir','reverse');

        if lcinit > 0
            filename = strcat(dir,pics_file,dataset,'_',num2str(lcinit,'%02d'))
        else
            filename = strcat(dir,pics_file,dataset,'_',num2str(lightcurveNo,'%02d'))

        end

       format long;
       hold off

       legend({'Hapke','Lambert', 'L-S', ['L-S+'+string(omega)+'L']}, 'Interpreter', 'latex', 'Fontsize', 10, 'Location', 'southwest')
       
       print(plotfigure,'-dpdf',filename);
       close(plotfigure)
   end
   
results.phasedCurve = [results.phasedCurve;artificalLightcurveData];
end

toc

%Various output parameters.  Artifical lightcurve is returned etc.
%phasedDatareturnMatrix(1,:) = [];
%results.phasedCurve(1,:) = [];
%correction = sortrows(correction,1);
% results.correction=correction;
results.pointOnLightCurve = pointOnLC'
results.datapoints = phasedDatareturnMatrix;
results.MikkoFile = LCSfileVariable;
results.phaseangle = phaseangle;
results.aspectAngle = aspect;
results.deltam = delta_m

end