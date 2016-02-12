function [cellarrayout] = RockPropertyPrediction(Zpred,Zwc,Pop,Sand1,Sand2,Sand3,Comp,GA)
%% Rock Physics Database and Rock Property Modelling Software
% Written by Matt Bolton December 2015
% Designed to incorporate all the rock properties of the sands, silts and mudrocks penetrated in the Tanzania offshore margin.
% Depth trends are then applied to the data after filtering for specific sand and caprock facies e.g. Age, Volume of Shale, Diagenesis, etc
% Used as inputs for Ratno's seismic forward modelling

%% FUNCTION INPUTS EXPLANATION
%Zpred = Input prospect depth below seabed (m)
%Zwc = Water depth for pressure calculation (m)
%Pop= Overpressure (psi)
%Sand = Reservoir shale content; 0 = Clean Sand, 1 = Shale. Recommend 0.05 for clean sand and 0.15 for silty sand
%Comp = Expected Sand diagenesis; 0=all, 1 = unconsolidated, 2 = consolidated, 3 = cemented
%GA = Expected reservoir age (Ma)

% Example Input
% RockPropertyPrediction(2950,1380,0,0.05,0.1,0.15,0,90)

%% KEY PARAMETERS
% Advanced parameters can be found on line 175 e.g. mineral properties, critical porosity, etc
%FitType=2; Not currently implemented
MaxZ=10000; % max depth (or pressure (psi)) to calculate to

% Data decimation
DEC=1; % Linear data decimation on LAS file read; 1=all data, 2=every other datapoint, etc
WinNo=1; % Non-linear, normalising data decimation - larger the number the greater the decimation of closely spaced data. 1=all data
OB_Scalar=1; % This can be used to increase non-linear decimation in higher sampled mudstone intervals

% Pressure filtering
Pco=1000; % Only include data in depth trend analysis which is at close to hydrostatic pressure (i.e. < Pco psi overpressued), 0=only hydrostatic - 10000=all data.

% Facies Classification
Sand=[Sand1 Sand2 Sand3];
Sfilt=0.05; % window around chosen Sandstone Vsh scenario, +/- Sand
COage=65; % Cut off Age where overburden and sandstones rock properties show different compaction trends. 65=K-T boundary
SndPorCO=0.08; % Net sand cut-off
SndVshCO=0.7; % Net sand cut-off
SndSwCO=0.8; % Net pay sand cut-off
SltVshCO=0.9; % Silt-Shale cut-off

% AVA modelling
MaxAngle=45; % Max angle to calcualte AVA over (degrees)
Iscalar=1520; % Model to Data scaling
Gscalar=1900; % Model to Data scaling

% Monte-Carlo Simulation
NI=5000; % Number of iteration for monte carlo simulation of rock properties
Dist=0.66; % Sample fraction error range - covariance to included in monte-carlo simulation, 1=full.


%% DATA READ

addpath(genpath('/apps/gsc/matlab-library/development/maps'));
filepath = {'/data/TZA/matlab/Depth_trend_analysis/LAS_Files/Input/'};
filename_string = 'las';

[files_in,nfiles] = directory_scan(filepath,filename_string); % files_in is a structure, .names cell array, .path cell array 
                                                              % directory_scan filters out the file names with the user supplied string and also returns the number of such files: nfiles that ...  
                                                              % ... can be used for index pre alocation in future steps 
%files_in.names = sort_nat(files_in.names);                    % Natural Sort file names 

[Well Zsb] = textread('Seabed.txt','%s %f','headerlines',1);

for i_file = 1:1:nfiles
    [well_data{i_file},well_header{i_file}] = read_las2_file([files_in.path{i_file},files_in.names{i_file}],0);
    fprintf('Read file %d\n',i_file);
      
    % Find NANS
    well_data{i_file}.logi =~ logical(sum(isnan(well_data{i_file}.curves),2));
    well_data{i_file}.logilength = sum(well_data{i_file}.logi,1);
    
    % Get rid of NANS
    well_data{i_file}.curves_wonans = well_data{i_file}.curves(well_data{i_file}.logi,:);
    
    % Collect number of data points in well
    well_no_datapoints(i_file,1) = size(well_data{i_file}.curves_wonans,1);
    
    % Add Seabed Depth to wells
    well_data{i_file}.Zsb = Zsb(strmatch(files_in.names{i_file},Well(:), 'exact'));
    
    % Error Check LAS file to Seabed depth match
    if well_data{i_file}.Zsb>0;
       continue
    else
        error('LAS file and Well seabed depth naming mismatch - please correct')
        end
    end



% Log Decimation

for i_file = 1:1:nfiles
    % Data decimation parameter
    well_data{i_file}.curves_dec = downsample(well_data{i_file}.curves_wonans,DEC);
     
    for i_curve = 1:size(well_data{i_file}.curves_dec,2);
        well_data{i_file}.(well_data{i_file}.curve_info{i_curve})=well_data{i_file}.curves_dec(:,i_curve);
    end
    
end

curve_info = well_data{1}.curve_info;
for i_curve = 1:size(well_data{i_file}.curves_dec,2);
   start_row = 1;
   for i_file = 1:1:nfiles
        %fprintf('Concatenating %s curves\n',well_data{1}.curve_info{i_curve})
        end_row = size(well_data{i_file}.(well_data{1}.curve_info{i_curve}),1)+start_row-1;
        curve_info{i_curve,4}(start_row:end_row,1) = well_data{i_file}.(well_data{1}.curve_info{i_curve});
        start_row = end_row + 1;
   end   
   eval(strcat(curve_info{i_curve,1},'= curve_info{',num2str(i_curve),',4};'))
end

%fprintf('Logs Read - If correct press any key to continue, else ctrl-c')
%well_data{1}.curve_info(:,1)
%pause
%PressureYN=input('Formation pressures available (Y) or assume hydrostatic (N)? Y/N:','s');
%AgeYN=input('Geological age data available (Y) or No (N)? Y/N:','s');
%fprintf('Continuing program......')
PressureYN='Y';
AgeYN='Y';

%% Initial Data filtering

% Geological Age filter
if AgeYN == 'Y'
    if GA>COage 
        GA_low=COage;
        GA_hi=200; % Base Jurassic
    else
        GA_low=0;
        GA_hi=COage;
    end
else
    Age(:)=0;
    COage=0;
    GA_low=0;
    GA_hi=0;
end


% Create diagensis log based on petrophyiscal cut offs and any other
% influencing factor.  Here Age is used, but this could also be T or P
% based.
for aa = 1:length(Vp); 
    if Comp == 0 & VSH(aa) < SndVshCO & PHIE(aa) >= SndPorCO
    Diag(aa,1)=0; % All Net sand points
        else if Comp>0 & VSH(aa) < SndVshCO & PHIE(aa) >= SndPorCO & Age(aa) <= COage
        Diag(aa,1)=1; % Unconsolidated Tertiary sands
            else if Comp>0 & VSH(aa) < SndVshCO & PHIE(aa) >= SndPorCO & Age(aa) >= COage
            Diag(aa,1)=2; % Consolidated Cretaceous sands
                else if Comp>0 & VSH(aa) < SndVshCO & PHIE(aa) < SndPorCO
                Diag(aa,1)=3; % Cemented sands
                    else
                    Diag(aa,1)=4;
                end
            end
        end
    end
end

% Overburden facies
RhoShl=DENS(VSH >= SltVshCO & PHIE < SndPorCO);

% Depth or pressure iteration
Z=0:1:MaxZ-1; % Z index with 1m sampling (requires all input depths are at this precision)

%% Advanced Parameters
% Rock Physics Modelling 
Por_min=0.0; % sandstone minimum porosity at infinite pressure
Por_min_ob=0.1; % mudrock minimum porosity at infinite pressure
%phiET=(-0.87*Sand)+0.87; % fraction of shale in pore space; 0 = shale matrix, 1 = all shale in pore space
phiET=0.25; % Fraction of Vsh which fills pore space
PhiC=0.5; % Clean Sand critical porosity
PhiCsh=0.7; % Clean Shale critical porosity
RhoMin=2; % Minimum shale density

% Rock property limits used to constrain the maximum possible properties (used in depth trend calculations)
RhoQ=2.65; % Quartz limits
VpQ=6038;
VsQ=4121;
K_sand=36.6; % Dry sandstone bulk moduli for fluid substitution
RhoI=2.8; % Illite rich shale limits
VpI=5264;
VsI=2719;
RhoS=2.5; % Smectite rich shale limits
VpS=3282;
VsS=2254;

sI=0.7; % Smecitite-Illite proportions for Young shales
sS=0.3;

RhoSI=(sI*RhoI)+(sS*RhoS); % Smectite-Illite shale limits
VpSI=(sI*VpI)+(sS*VpS);
VsSI=(sI*VsI)+(sS*VsS);

% Over pressure modelling coefficients
% See workflow in supporting documentation
% Mudstone parameters empirically derived
a_mud=1200;
b_mud=700;
c_mud=0.00143;

% Sandstone parameters
a_sst=1500;
b_sst=700;
c_sst=0.00005;

% Fluid properties
RhoSea=1.1; % Seawater density g/cc
HCcolumn=50; % m; this is required to determine an average gas buoyancy effect
RhoW=1.0;
RhoG=0.2;

% Temperature and Pressure
T_pred=5+(42*(Zpred/1000));
Pp=mean(MDT_Pform)*0.006895; % Average Pf for fluid property calculation (MPa)

% Temperature profile
Tp(1:MaxZ,1)=5;
for temp=Zwc:MaxZ;  
Tp(temp)=Tp(temp)+(42*((Z(temp)-Zwc)/1000));
end

%RhoW=mean(rhob) % Used as average overburden brine properties for prediction
RhoW=1;

% Average overburden profile
RhoOB=RhoShl; % Overburden density data
ZOB=Depth_TVDml(VSH >= SltVshCO & PHIE < SndPorCO);
%ini(MaxZ,1)=0;

% Exponential curve fit for shale trend
dRhoOB=ZOB;
GRhoOB=log((RhoOB-RhoMin)/(RhoMin-RhoSI));
mRhoOB=GRhoOB\dRhoOB;
betaRhoOB=real(-1/mRhoOB);

GrOB=exp(-betaRhoOB*ZOB);
d=RhoOB-RhoSI; % known data
mx=RhoSI; % max shale density
mRhoMudOB=ModelFitv2(d,GrOB,betaRhoOB,mx);
RhoL = mx -((mx-mRhoMudOB)*exp(-betaRhoOB*Z));

Pl(1:MaxZ)=0;
for pp=1:MaxZ;
  Pl(pp+1) = Pl(pp)+(1000*1E-6*9.81*(RhoL(pp)));  
end


% Integration of density trend (RhoL) to define average lithostatic


for bb=1:nfiles;
    %Pre-allocate the Lithstatic pressure array
    Plith = zeros(length(well_data{bb}.Depth_TVDml),1);
    %Calculate lithostatic pressure in overburden above first data point
    Plith(1) = Pl(round((well_data{bb}.Depth_TVDml(1)))); 
    %Calculate Overburden density values for each depth data point
    well_data{bb}.RhoL = mx-((mx-mRhoMudOB)*exp(-betaRhoOB*well_data{bb}.Depth_TVDml));   
    %Water Column Pressure
    wc(bb)=1000*1E-6*9.81*well_data{bb}.Zsb*RhoSea;

    for aa=1:length(well_data{bb}.Depth_TVDml)-1;
            %Create hydrostatic pressure points
            well_data{bb}.Phydro(aa,1) = 1000*1E-6*9.81*((RhoW*well_data{bb}.Depth_TVDml(aa))+(RhoSea*well_data{bb}.Zsb));
            % Calculate lithostatic pressure for each depth data point
            Plith(aa+1) = Plith(aa)+(1000*1E-6*9.81*(well_data{bb}.Depth_TVDml(aa+1)-well_data{bb}.Depth_TVDml(aa))*(well_data{bb}.RhoL(aa)));
            %Write to new cell array
            well_data{bb}.Plith(aa,1) = (Plith(aa)+wc(bb));
    end
    
    %Account for end of array
    well_data{bb}.Phydro(aa+1,1)=well_data{bb}.Phydro(aa,1);
    well_data{bb}.Plith(aa+1,1)= well_data{bb}.Plith(aa,1);
    
        if PressureYN == 'Y'
            %Calculate vertical effective stress for each data point
            well_data{bb}.Pves = (well_data{bb}.Plith-well_data{bb}.MDT_Pform);
            %Calculate over-pressure for each data point
            well_data{bb}.Pover = (well_data{bb}.MDT_Pform-well_data{bb}.Phydro);
        else
            %Calculate vertical effective stress for each data point
            well_data{bb}.Pves = (well_data{bb}.Plith-well_data{bb}.Phydro);  
            well_data{bb}.Pover = 0;
        end  
end

start_point = 1;
for i_file = 1:1:nfiles
        end_point = size(well_data{i_file}.curves_dec,1)+start_point-1;
        Pover(start_point:end_point,1) = 145.0377*well_data{i_file}.Pover;
        Pves(start_point:end_point,1) = 145.0377*well_data{i_file}.Pves;
        Plith(start_point:end_point,1) = 145.0377*well_data{i_file}.Plith;
        Phydro(start_point:end_point,1) = 145.0377*well_data{i_file}.Phydro;
        start_point = end_point + 1;
end

% Calculate pore pressure 
Phy(1:MaxZ)=1E-3*9.81*Z*RhoSea;
for hyd=Zwc:MaxZ
    Phy(hyd)=1E-3*(9.81*(Z(hyd)*RhoW)); %MPa
end

% Calculate pressures based on input parameters (Prospect modelling)
% Lithostatic
Plith_pred=145.0377*((Pl(Zpred))+((1000*1E-6*9.81)*(Zwc*RhoSea)));
% Hydrostatic
Phy_pred=145.0377*(1000*1E-6*9.81)*((Zpred*RhoW)+(Zwc*RhoSea));

% Prospect Pore pressure prediction
P_pred_sand=Phy_pred+Pop; %psi
Pp_pred_sand=0.006895*P_pred_sand; % MPa

% Buoyancy effect in gas leg: Pform_gas > Pform_wtr at same depth
Buoy=(RhoW-RhoG)*9.81*HCcolumn*0.145; % Gas buoyancy effect

% Prospect vertical effective stress
Pves_pred=Plith_pred-Phy_pred; %psi
PVES=Pves_pred*0.006895; %MPa


%% Overburden Trend Fitting

% Filter vertical effective stress points to get corresponding data for
% curve fitting to overburden
Pves_filt_slt_all=Pves(VSH >= SndVshCO & VSH < SltVshCO & PHIE < SndPorCO & Age>=GA_low & Age<=GA_hi & Pover < Pco);
Pves_filt_shl_all=Pves(VSH >= SltVshCO & PHIE < SndPorCO & Age>=GA_low & Age<=GA_hi & Pover < Pco);


% Filter Overburden facies
[Pves_filt_slt PhiT_filt_slt]=DataSamplingNormV2(Pves_filt_slt_all,NEUT(VSH >= SndVshCO & VSH < SltVshCO & PHIE < SndPorCO & Age>=GA_low & Age<=GA_hi & Pover < Pco),OB_Scalar*WinNo);
[Pves_filt_shl PhiT_filt_shl]=DataSamplingNormV2(Pves_filt_shl_all,NEUT(VSH >= SltVshCO & PHIE < SndPorCO & Age>=GA_low & Age<=GA_hi & Pover < Pco),OB_Scalar*WinNo);
[Pves_filt_shl RhoShl_filt]=DataSamplingNormV2(Pves_filt_shl_all,DENS(VSH >= SltVshCO & PHIE < SndPorCO & Age>=GA_low & Age<=GA_hi & Pover < Pco),OB_Scalar*WinNo);
[Pves_filt_shl VpShl_filt]=DataSamplingNormV2(Pves_filt_shl_all,Vp(VSH >= SltVshCO & PHIE < SndPorCO & Age>=GA_low & Age<=GA_hi & Pover < Pco),OB_Scalar*WinNo);
[Pves_filt_shl VsShl_filt]=DataSamplingNormV2(Pves_filt_shl_all,Vs_Fast(VSH >= SltVshCO & PHIE < SndPorCO & Age>=GA_low & Age<=GA_hi & Pover < Pco),OB_Scalar*WinNo);
[Pves_filt_slt RhoSlt_filt]=DataSamplingNormV2(Pves_filt_slt_all,DENS(VSH >= SndVshCO & VSH < SltVshCO & PHIE < SndPorCO & Age>=GA_low & Age<=GA_hi & Pover < Pco),OB_Scalar*WinNo);
[Pves_filt_slt VpSlt_filt]=DataSamplingNormV2(Pves_filt_slt_all,Vp(VSH >= SndVshCO & VSH < SltVshCO & PHIE < SndPorCO & Age>=GA_low & Age<=GA_hi & Pover < Pco),OB_Scalar*WinNo);
[Pves_filt_slt VsSlt_filt]=DataSamplingNormV2(Pves_filt_slt_all,Vs_Fast(VSH >= SndVshCO & VSH < SltVshCO & PHIE < SndPorCO & Age>=GA_low & Age<=GA_hi & Pover < Pco),OB_Scalar*WinNo);
[Pves_filt_slt Age_filt_slt]=DataSamplingNormV2(Pves_filt_slt_all,Age(VSH >= SndVshCO & VSH < SltVshCO & PHIE < SndPorCO & Age>=GA_low & Age<=GA_hi & Pover < Pco),OB_Scalar*WinNo);
[Pves_filt_shl Age_filt_shl]=DataSamplingNormV2(Pves_filt_shl_all,Age(VSH >= SltVshCO & PHIE < SndPorCO & Age>=GA_low & Age<=GA_hi & Pover < Pco),OB_Scalar*WinNo);
[Pves_filt_slt Vsh_filt_slt]=DataSamplingNormV2(Pves_filt_slt_all,VSH(VSH >= SndVshCO & VSH < SltVshCO & PHIE < SndPorCO & Age>=GA_low & Age<=GA_hi & Pover < Pco),OB_Scalar*WinNo);
[Pves_filt_shl Vsh_filt_shl]=DataSamplingNormV2(Pves_filt_shl_all,VSH(VSH >= SltVshCO & PHIE < SndPorCO & Age>=GA_low & Age<=GA_hi & Pover < Pco),OB_Scalar*WinNo);


% Caprock property fitting procedure:

% Caprock property trends
% Similar approach for caprock 
if GA >= COage
    RhoMud=RhoI;
    VpMud=VpI;
    VsMud=VsI;
else
    RhoMud=RhoSI;
    VpMud=VpSI;
    VsMud=VsSI;
end

dShPhiT=Pves_filt_shl;
GShPhiT=log((PhiT_filt_shl-Por_min_ob)/(Por_min_ob-PhiCsh));
mShPhiT=GShPhiT\dShPhiT;
betaShPhiT=real(-1/mShPhiT);
GSh=exp(-betaShPhiT*Pves_filt_shl); % Compaction gradient for caprock

Shl_por_pred = Por_min_ob-((Por_min_ob-PhiCsh)*exp(-betaShPhiT*(Pves_pred)));
Shl_por_trend = Por_min_ob-((Por_min_ob-PhiCsh)*exp(-betaShPhiT*Z));

d=RhoShl_filt-RhoMud; % known data
mx=RhoMud; % Quartz density
mRhoMud=ModelFitv2(d,GSh,betaShPhiT,mx);
RhoShl_pred = RhoMud-((RhoMud-mRhoMud)*exp(-betaShPhiT*Pves_pred));
RhoShl_trend = RhoMud-((RhoMud-mRhoMud)*exp(-betaShPhiT*Z));

d=VpShl_filt-VpMud; % known data
mx=VpMud; % Quartz density
mVpMud=ModelFitv2(d,GSh,betaShPhiT,mx);
VpShl_pred = VpMud-((VpMud-mVpMud)*exp(-betaShPhiT*Pves_pred));
VpShl_trend = VpMud-((VpMud-mVpMud)*exp(-betaShPhiT*Z));

d=VsShl_filt-VsMud; % known data
mx=VsMud; % Quartz density
mVsMud=ModelFitv2(d,GSh,betaShPhiT,mx);
VsShl_pred = VsMud-((VsMud-mVsMud)*exp(-betaShPhiT*Pves_pred));
VsShl_trend = VsMud-((VsMud-mVsMud)*exp(-betaShPhiT*Z));


if GA >= COage
    RhoMud=(((SndVshCO+SltVshCO)/2)*RhoI)+((1-((SndVshCO+SltVshCO)/2))*RhoQ);
    VpMud=(((SndVshCO+SltVshCO)/2)*VpI)+((1-((SndVshCO+SltVshCO)/2))*VpQ);
    VsMud=(((SndVshCO+SltVshCO)/2)*VsI)+((1-((SndVshCO+SltVshCO)/2))*VsQ);
else
    RhoMud=(((SndVshCO+SltVshCO)/2)*RhoSI)+((1-((SndVshCO+SltVshCO)/2))*RhoQ);
    VpMud=(((SndVshCO+SltVshCO)/2)*VpSI)+((1-((SndVshCO+SltVshCO)/2))*VpQ);
    VsMud=(((SndVshCO+SltVshCO)/2)*VsSI)+((1-((SndVshCO+SltVshCO)/2))*VsQ);
end

% Silt overburden rock properties
dSlPhiT=Pves_filt_slt;
GSlPhiT=log((PhiT_filt_slt-Por_min_ob)/(Por_min_ob-PhiCsh));
mSlPhiT=GSlPhiT\dSlPhiT;
betaSlPhiT=real(-1/mSlPhiT);
GSl=exp(-betaSlPhiT*Pves_filt_slt); % Compaction gradient for caprock

Slt_por_pred = Por_min_ob-((Por_min_ob-PhiCsh)*exp(-betaSlPhiT*(Pves_pred)));
Slt_por_trend = Por_min_ob-((Por_min_ob-PhiCsh)*exp(-betaSlPhiT*Z));

d=RhoSlt_filt-RhoMud; % known data
mx=RhoMud; % Quartz density
mRhoMud=ModelFitv2(d,GSl,betaSlPhiT,mx);
RhoSlt_pred = RhoMud-((RhoMud-mRhoMud)*exp(-betaSlPhiT*Pves_pred));
RhoSlt_trend = RhoMud-((RhoMud-mRhoMud)*exp(-betaSlPhiT*Z));

d=VpSlt_filt-VpMud; % known data
mx=VpMud; % Quartz density
mVpMud=ModelFitv2(d,GSl,betaSlPhiT,mx);
VpSlt_pred = VpMud-((VpMud-mVpMud)*exp(-betaSlPhiT*Pves_pred));
VpSlt_trend = VpMud-((VpMud-mVpMud)*exp(-betaSlPhiT*Z));

d=VsSlt_filt-VsMud; % known data
mx=VsMud; % Quartz density
mVsMud=ModelFitv2(d,GSl,betaSlPhiT,mx);
VsSlt_pred = VsMud-((VsMud-mVsMud)*exp(-betaSlPhiT*Pves_pred));
VsSlt_trend = VsMud-((VsMud-mVsMud)*exp(-betaSlPhiT*Z));

%% Sandstone Trend Fitting
for N=1:3

% Initial data filtering

% Volume of shale filter included +/- Sfilt of input value
if Sand(N)-Sfilt>0 % Limit low case to 0
    Sand_low=Sand(N)-Sfilt;
else
    Sand_low=0;
end
Sand_hi=Sand(N)+Sfilt;

Clay=Sand(N); % fraction; Volume of clay
RhoSd=((1-Clay)*RhoQ)+(Clay*RhoSI); % Sand limits
VpSd=((1-Clay)*VpQ)+(Clay*VpSI);
VsSd=((1-Clay)*VsQ)+(Clay*VsSI);
PhiC=(PhiC-(phiET*Clay)); % critical porosity can be reduced through shale in pore space and cementation


% Filter vertical effective stress points to get corresponding data for
% curve fitting to sand
Pves_filt_sand_gas_all=Pves(VSH <= Sand_hi & VSH >= Sand_low & Diag==Comp & SW <= SndSwCO & Pover < Pco);
Pves_filt_sand_wtr_all=Pves(VSH <= Sand_hi & VSH >= Sand_low & Diag==Comp & SW > SndSwCO & Pover < Pco);
Pves_filt_sand_all=Pves(VSH <= Sand_hi & VSH >= Sand_low & Diag==Comp & Pover < Pco);

% Filter sand points based on saturation, volume of shale and porosity
[Pves_filt_sand_gas VpGas_filt]=DataSamplingNormV2(Pves_filt_sand_gas_all,Vp(VSH <= Sand_hi & VSH >= Sand_low & Diag==Comp & SW <= SndSwCO & Pover < Pco),WinNo);
[Pves_filt_sand_wtr VpWtr_filt]=DataSamplingNormV2(Pves_filt_sand_wtr_all,Vp(VSH <= Sand_hi & VSH >= Sand_low & Diag==Comp & SW > SndSwCO & Pover < Pco),WinNo);
[Pves_filt_sand_gas VsGas_filt]=DataSamplingNormV2(Pves_filt_sand_gas_all,Vs_Fast(VSH <= Sand_hi & VSH >= Sand_low & Diag==Comp & SW <= SndSwCO & Pover < Pco),WinNo);
[Pves_filt_sand_wtr VsWtr_filt]=DataSamplingNormV2(Pves_filt_sand_wtr_all,Vs_Fast(VSH <= Sand_hi & VSH >= Sand_low & Diag==Comp & SW > SndSwCO & Pover < Pco),WinNo);
[Pves_filt_sand_gas RhoGas_filt]=DataSamplingNormV2(Pves_filt_sand_gas_all,DENS(VSH <= Sand_hi & VSH >= Sand_low & Diag==Comp & SW <= SndSwCO & Pover < Pco),WinNo);
[Pves_filt_sand_wtr RhoWtr_filt]=DataSamplingNormV2(Pves_filt_sand_wtr_all,DENS(VSH <= Sand_hi & VSH >= Sand_low & Diag==Comp & SW > SndSwCO & Pover < Pco), WinNo);
[Pves_filt_sand Por_filt]=DataSamplingNormV2(Pves_filt_sand_all,PHIT(VSH <= Sand_hi & VSH >= Sand_low & Diag==Comp & Pover < Pco),WinNo);
[Pves_filt_sand_wtr Por_filt_wtr]=DataSamplingNormV2(Pves_filt_sand_wtr_all,PHIT(VSH <= Sand_hi & VSH >= Sand_low & Diag==Comp & SW > SndSwCO & Pover < Pco),WinNo);
[Pves_filt_sand_gas Por_filt_gas]=DataSamplingNormV2(Pves_filt_sand_gas_all,PHIT(VSH <= Sand_hi & VSH >= Sand_low & Diag==Comp & SW <= SndSwCO & Pover < Pco),WinNo);
[Pves_filt_sand_gas Age_filt_gas]=DataSamplingNormV2(Pves_filt_sand_gas_all,Age(VSH <= Sand_hi & VSH >= Sand_low & Diag==Comp & SW <= SndSwCO & Pover < Pco),WinNo);
[Pves_filt_sand_wtr Age_filt_wtr]=DataSamplingNormV2(Pves_filt_sand_wtr_all,Age(VSH <= Sand_hi & VSH >= Sand_low & Diag==Comp & SW > SndSwCO & Pover < Pco),WinNo);
[Pves_filt_sand vsh_filt]=DataSamplingNormV2(Pves_filt_sand_all,VSH(VSH <= Sand_hi & VSH >= Sand_low & Diag==Comp & Pover < Pco),WinNo);


% Bulk Rock Properties for use in fluid substitution
[AI_Gas,VpVsR_Gas,Mu_Gas,K_Gas] = RockProps(VpGas_filt,VsGas_filt,RhoGas_filt);
[AI_Wtr,VpVsR_Wtr,Mu_Wtr,K_Wtr] = RockProps(VpWtr_filt,VsWtr_filt,RhoWtr_filt);


% CURVE FITTING (Vertical effective stress vs Rock Property)
% This version of the software uses a hybrid curve fitting procedure.
% In a later version it is intended to add multiple fit types.

% Sandstone fitting procedure:
% Limits to minimum and maximum porosities are manually input.
% Filtered data points included only
% A least squares solution to the exponential compaction gradient is then fitted to Sand Porosity v effective stress points
% This compaction gradient to porosity is then assumed for all sandstone 
% rock properties i.e. compaction trend is most sensitive to porosity
% decrease with increased effective stress.
% A singular initial constraint however is applied to the sandstone rock property
% trends - the maximum value at infinte effective stress.  Here the
% properties for quartz are used.

% Porosity Prediction
% Least squares solution to Por=Por_min-(Por_min-phic)*exp(-betaPor*Pves)
% solving for betaPor

dPor=Pves_filt_sand;
GPor=log((Por_filt-Por_min)/(Por_min-PhiC));
mPor=GPor\dPor;
betaPor=real(-1/mPor);

Porosity_predicted = Por_min-((Por_min-PhiC)*exp(-betaPor*(Pves_pred)));
Por_trend(:,N) = Por_min-((Por_min-PhiC)*exp(-betaPor*Z));


% Gradient matrix for all subsequent least squares depth trend equations
% Porosity depth trend used as a known compaction gradient
Ggas=exp(-betaPor*Pves_filt_sand_gas); % Compaction gradient for sands
Gwtr=exp(-betaPor*Pves_filt_sand_wtr);

% Gas Sand Rock Properties
% Rho Gas
d=RhoGas_filt-RhoSd; % known data
mx=RhoSd; % Quartz density
mRhoGas=ModelFitv2(d,Ggas,betaPor,mx);
RhoGas_pred = RhoSd-((RhoSd-mRhoGas)*exp(-betaPor*Pves_pred));
RhoGas_trend(:,N) = RhoSd-((RhoSd-mRhoGas)*exp(-betaPor*Z));


% Vp Gas
d=VpGas_filt-VpSd; % known data
mx=VpSd; % Quartz Vp
mVpGas=ModelFitv2(d,Ggas,betaPor,mx);
VpGas_pred = VpSd-((VpSd-mVpGas)*exp(-betaPor*Pves_pred));
VpGas_trend(:,N) = VpSd-((VpSd-mVpGas)*exp(-betaPor*Z));


% Vs Gas
d=VsGas_filt-VsSd; % known data
mx=VsSd; % Quartz Vs
mVsGas=ModelFitv2(d,Ggas,betaPor,mx);
VsGas_pred = VsSd-((VsSd-mVsGas)*exp(-betaPor*Pves_pred));
VsGas_trend(:,N) = VsSd-((VsSd-mVsGas)*exp(-betaPor*Z));


% Water Sand Rock Properties

% Rho Water
d=RhoWtr_filt-RhoSd; % known data
mx=RhoSd; % Quartz density
mRhoWtr=ModelFitv2(d,Gwtr,betaPor,mx);
RhoWtr_pred = RhoSd-((RhoSd-mRhoWtr)*exp(-betaPor*Pves_pred));
RhoWtr_trend(:,N) = RhoSd-((RhoSd-mRhoWtr)*exp(-betaPor*Z));


% Vp Water
d=VpWtr_filt-VpSd; % known data
mx=VpSd; % Quartz Vp
mVpWtr=ModelFitv2(d,Gwtr,betaPor,mx);
VpWtr_pred = VpSd-((VpSd-mVpWtr)*exp(-betaPor*Pves_pred));
VpWtr_trend(:,N) = VpSd-((VpSd-mVpWtr)*exp(-betaPor*Z));


% Vs Water
d=VsWtr_filt-VsSd; % known data
mx=VsSd; % Quartz Vs
mVsWtr=ModelFitv2(d,Gwtr,betaPor,mx);
VsWtr_pred = VsSd-((VsSd-mVsWtr)*exp(-betaPor*Pves_pred));
VsWtr_trend(:,N) = VsSd-((VsSd-mVsWtr)*exp(-betaPor*Z));

end

%% Outputs for Seismic forward modelling
cellarrayout{1,1} = 'Effective_Stress';
cellarrayout{1,2} = Z';
cellarrayout{2,1} = 'Rho_stiff_overburden';
cellarrayout{2,2} = RhoSlt_trend';
cellarrayout{3,1} = 'Vp_stiff_overburden';
cellarrayout{3,2} = VpSlt_trend';
cellarrayout{4,1} = 'Vs_stiff_overburden';
cellarrayout{4,2} = RhoSlt_trend';
cellarrayout{5,1} = 'Porosity_stiff_overburden';
cellarrayout{5,2} = Slt_por_trend';
cellarrayout{6,1} = 'Rho_soft_overburden';
cellarrayout{6,2} = RhoShl_trend';
cellarrayout{7,1} = 'Vp_soft_overburden';
cellarrayout{7,2} = VpShl_trend';
cellarrayout{8,1} = 'Vs_soft_overburden';
cellarrayout{8,2} = RhoShl_trend';
cellarrayout{9,1} = 'Porosity_soft_overburden';
cellarrayout{9,2} = Shl_por_trend';
cellarrayout{10,1} = 'Rho_gas_sand';
cellarrayout{10,2} = RhoGas_trend;
cellarrayout{11,1} = 'Vp_gas_sand';
cellarrayout{11,2} = VpGas_trend;
cellarrayout{12,1} = 'Vs_gas_sand';
cellarrayout{12,2} = VsGas_trend;
cellarrayout{13,1} = 'Rho_brine_sand';
cellarrayout{13,2} = RhoWtr_trend;
cellarrayout{14,1} = 'Vp_brine_sand';
cellarrayout{14,2} = VpWtr_trend;
cellarrayout{15,1} = 'Vs_brine_sand';
cellarrayout{15,2} = VsWtr_trend;
cellarrayout{16,1} = 'Porosity_sand';
cellarrayout{16,2} = Por_trend;
cellarrayout{19,1} = 'Sand_Facies';
cellarrayout{19,2} = Sand;

end



