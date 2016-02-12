function AVA_AllData_edit2_backus(Zpred,Zwc,Pop,Sand,Comp,GA,Shale,FitType,NI)
%% ------------------ Disclaimer  ------------------
% 
% BG Group plc or any of its respective subsidiaries, affiliates and 
% associated companies (or by any of their respective officers, employees 
% or agents) makes no representation or warranty, express or implied, in 
% respect to the quality, accuracy or usefulness of this repository. The code
% is this repository is supplied with the explicit understanding and 
% agreement of recipient that any action taken or expenditure made by 
% recipient based on its examination, evaluation, interpretation or use is 
% at its own risk and responsibility.
% 
% No representation or warranty, express or implied, is or will be made in 
% relation to the accuracy or completeness of the information in this 
% repository and no responsibility or liability is or will be accepted by 
% BG Group plc or any of its respective subsidiaries, affiliates and 
% associated companies (or by any of their respective officers, employees 
% or agents) in relation to it.
%% ------------------ License  ------------------ 
% GNU GENERAL PUBLIC LICENSE Version 3, 29 June 2007
%% github
% https://github.com/AnalysePrestackSeismic/
%% ------------------ FUNCTION DEFINITION ---------------------------------
%% Rock Physics Database and Depth Trend Modelling Software
% Written by Matt Bolton April 2015
% Designed to incorporate all the rock properties of the sands, silts and mudrocks penetrated in the Tanzania offshore margin.
% Depth trends are then applied to the data after filtering for specific sand and caprock facies e.g. Age, Volume of Shale, Diagenesis, etc
% Seismic interfaces for can then be modelled at any depth based on trends fitted. 
% Reflectivity with angle is estimated for the isotropic case using the full Zoeppritz solution.  An anisotropic case is also included based on the Shuey approximation.

%% FUNCTION INPUTS EXPLANATION
%Zpred = Input prospect depth below seabed (m)
%Zwc = Water depth for pressure calculation (m)
%Pop= Overpressure (psi)
%Sand = Reservoir shale content; 0 = Clean Sand, 1 = Shale. Recommend 0.05 for clean sand and 0.15 for silty sand
%Comp = Expected Sand diagenesis; 0=all, 1 = unconsolidated, 2 = consolidated, 3 = cemented
%GA = Expected reservoir age (Ma)
%Shale = Expected overburden; 0=Shale (soft) or 1=Slt (stiff)

% Example Input
% AVA_AllData_edit2(2950,1380,0,0.05,0,90,0,2,5000)

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

% Anisotropy
delta1=0.24; % Overburden 
delta2=0; % Reservoir
epsilon1=0.24; % Overburden
epsilon2=0; % Reservoir

% Depth trend intervals
Zmin=2000;
Zmax=8000;
Zit=(Zmax-Zmin)/5;
Zrange=[Zmin:Zit:Zmax]';

% Backus averaging
FL=20; %m
ds=0.1524; %m
F=ceil(FL/ds);

SeisF=input('Log Frequency (Y) or Seismic Frequency (N)? Y/N:','s');


%% DATA READ

addpath(genpath('/apps/gsc/matlab-library/development/maps'));
filepath = {'/data/TZA/matlab/Depth_trend_analysis/'};
filename_string = 'las';

[files_in,nfiles] = directory_scan(filepath,filename_string); % files_in is a structure, .names cell array, .path cell array 
                                                              % directory_scan filters out the file names with the user supplied string and also returns the number of such files: nfiles that ...  
                                                              % ... can be used for index pre alocation in future steps 
%files_in.names = sort_nat(files_in.names);                    % Natural Sort file names 

[Well Zsb] = textread('Seabed.txt','%s %f','headerlines',1);

for i_file = 1:1:nfiles
    [well_data{i_file},well_header{i_file}] = read_las2_file([files_in.path{i_file},files_in.names{i_file}],0);
    fprintf('Read file %d\n',i_file);
    
    if SeisF == 'N';
    % Arithmetic average
    [L W]=size(well_data{1}.curve_info); % determine number of curves for filtering
    [well_data{i_file}.curves] = ArithmNew_MB(well_data{i_file}.curves,F,L);
    end    
    
    % Find NANS
    well_data{i_file}.logi =~ logical(sum(isnan(well_data{i_file}.curves),2));
    well_data{i_file}.logilength = sum(well_data{i_file}.logi,1);
    
    % Get rid of NANS
    well_data{i_file}.curves_wonans = well_data{i_file}.curves(well_data{i_file}.logi,:);
    
    % Collect number of data points in well
    well_no_datapoints(i_file,1) = size(well_data{i_file}.curves_wonans,1);
    
    % Add Seabed Depth to wells
    well_data{i_file}.Zsb = Zsb(i_file);
end

%well_no_datapoints(:,2) = ceil(well_no_datapoints(:,1)./min(well_no_datapoints(:,1)));
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

well_data{1}.curve_info(:,1)
PressureYN=input('Formation pressures available (Y) or assume hydrostatic (N)? Y/N:','s');
AgeYN=input('Geological age data available (Y) or No (N)? Y/N:','s');
fprintf('Continuing program......')


%% Initial Data filtering

% Volume of shale filter included +/- Sfilt of input value
if Sand-Sfilt>0 % Limit low case to 0
    Sand_low=Sand-Sfilt;
else
    Sand_low=0;
end
Sand_hi=Sand+Sfilt;

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
Clay=Sand; % fraction; Volume of clay 
Feldspar=0; % fraction; Volume of feldpar
Calcite=0; % fraction; Volume of calcite
Por_min=0.0; % sandstone minimum porosity at infinite pressure
Por_min_ob=0.1; % mudrock minimum porosity at infinite pressure
%phiET=(-0.87*Sand)+0.87; % fraction of shale in pore space; 0 = shale matrix, 1 = all shale in pore space
phiET=0.25; % Fraction of Vsh which fills pore space
fcement=0.0; % cement fraction
PhiC=0.5; % Clean Sand critical porosity
PhiCsh=0.7; % Clean Shale critical porosity
PhiC=(PhiC-(phiET*Clay)); % critical porosity can be reduced through shale in pore space and cementation
shearfact=0.6; % Shear wave reduction factor
Coord=20; % Co-ordination number
RhoMin=2; % Minimum shale density

% Rock property limits used to constrain the maximum possible properties (used in depth trend calculations)
RhoQ=2.65; % Quartz limits
VpQ=6038;
VsQ=4121;
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

RhoSd=((1-Clay)*RhoQ)+(Clay*RhoSI); % Sand limits
VpSd=((1-Clay)*VpQ)+(Clay*VpSI);
VsSd=((1-Clay)*VsQ)+(Clay*VsSI);


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

% Batzle and Wang inputs
method=2; %1, use Gas Index in Oil other numbers, use GOR (L/L)
sal=23230; %ppm
og=15; % oil gravity (API no)
gg=0.57; % gas gravity (Specific gravity)
giib=1; % Gas Index in Brine
giio=1; % Gas Index in Oil
gor=0.1; % Gas Oil Ratio (L/L)
So=0; % Oil saturation
Sg=0.0; % Gas saturation
T_pred=5+(42*(Zpred/1000));
Pp=mean(MDT_Pform)*0.006895; % Average Pf for fluid property calculation (MPa)
%Tp=mean(T); % Average temperature of input data for fluid property calculation
Tp=90;

%% BACKUS AVERAGE
if SeisF == 'N'
   [Vp,Vs_Fast,DENS]=BackusNew_MB(Vp,Vs_Fast,DENS,F);
end

%% TEMPERATURE & PRESSURE CALCULATIONS
% All depth trends should be calculated wrt to vertical effective stress if
% basin is not hydrostatically pressured.

% EDITABLE CONTENT SEE LINE 56

% CALCULATIONS
% Batzle and wang calculations
[Kreuss,rhoeff,Kvoigt,vpb,rhob,Kb,vpo,rhoo,Ko,vpg,rhog,Kg,gor]=flprop(method,sal,og,gg,gor,giib,giio,Pp,Tp,So,Sg);

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

%dPhiOB=ZOB;
%GPhiOB=log((PhiTob-Por_min_ob)/(Por_min_ob-PhiCsh));
%mPhiOB=GPhiOB\dPhiOB;
%betaPhiOB=real(-1/mPhiOB);


%GrOB=exp(-betaPhiOB*ZOB);
%d=RhoOB-RhoI; % known data
%mx=RhoI; % max shale density
%mRhoMudOB=ModelFitv2(d,GrOB,betaPhiOB,mx);
%RhoL = RhoI-((RhoI-mRhoMudOB)*exp(-betaPhiOB*Z));

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
        %curve_info{i_curve,4}(start_row:end_row,1) = well_data{i_file}.(well_data{1}.curve_info{i_curve});
        Pover(start_point:end_point,1) = 145.0377*well_data{i_file}.Pover;
        Pves(start_point:end_point,1) = 145.0377*well_data{i_file}.Pves;
        Plith(start_point:end_point,1) = 145.0377*well_data{i_file}.Plith;
        Phydro(start_point:end_point,1) = 145.0377*well_data{i_file}.Phydro;
        start_point = end_point + 1;
end

% Calculate pressures based on input parameters (Prospect modelling)
% Lithostatic
Plith_pred=145.0377*((Pl(Zpred))+((1000*1E-6*9.81)*(Zwc*RhoSea)));
% Hydrostatic
Phy_pred=145.0377*(1000*1E-6*9.81)*((Zpred*RhoW)+(Zwc*RhoSea));

% Prospect Pore pressure prediction
P_pred_sand=Phy_pred+Pop %psi
Pp_pred_sand=0.006895*P_pred_sand; % MPa

% Prospect fluid properties predicted based on inputs
[Kreuss,rhoeff,Kvoigt,vpb,RhoW_pred,Kb,vpo,rhoo,Ko,vpg,RhoG_pred,Kg,gor]=flprop(method,sal,og,gg,gor,giib,giio,Pp_pred_sand,T_pred,So,Sg);

% Buoyancy effect in gas leg: Pform_gas > Pform_wtr at same depth
Buoy=(RhoW_pred-RhoG_pred)*9.81*HCcolumn*0.145; % Gas buoyancy effect

% Prospect vertical effective stress
Pves_pred=Plith_pred-Phy_pred; %psi
PVES=Pves_pred*0.006895; %MPa

%% DATA FILTERING
% Using prospect inputs the database is filtered to include only data which
% is appropriate for the modelling scenario.


% Filter vertical effective stress points to get corresponding data for
% curve fitting to sand
Pves_filt_sand_gas_all=Pves(VSH <= Sand_hi & VSH >= Sand_low & Diag==Comp & SW <= SndSwCO & Pover < Pco);
Pves_filt_sand_wtr_all=Pves(VSH <= Sand_hi & VSH >= Sand_low & Diag==Comp & SW > SndSwCO & Pover < Pco);
Pves_filt_sand_all=Pves(VSH <= Sand_hi & VSH >= Sand_low & Diag==Comp & Pover < Pco);

% Filter vertical effective stress points to get corresponding data for
% curve fitting to overburden
Pves_filt_slt_all=Pves(VSH >= SndVshCO & VSH < SltVshCO & PHIE < SndPorCO & Age>=GA_low & Age<=GA_hi & Pover < Pco);
Pves_filt_shl_all=Pves(VSH >= SltVshCO & PHIE < SndPorCO & Age>=GA_low & Age<=GA_hi & Pover < Pco);

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

if Shale == 0
    RhoOB_filt=RhoShl_filt;
    VpOB_filt=VpShl_filt;
    VsOB_filt=VsShl_filt;
    Pves_filt_OB=Pves_filt_shl;
    PhiT_filt_OB=PhiT_filt_shl;
    
    if GA >= COage
            RhoMud=RhoI;
            VpMud=VpI;
            VsMud=VsI;
        else
            RhoMud=RhoSI;
            VpMud=VpSI;
            VsMud=VsSI;
        end
   
else
    
    Pves_filt_OB=Pves_filt_slt;
    RhoOB_filt=RhoSlt_filt;
    VpOB_filt=VpSlt_filt;
    VsOB_filt=VsSlt_filt;
    PhiT_filt_OB=PhiT_filt_slt;
    % Voight average for silt end member rock properties using Slt and Shl cut-off
    if GA >= COage
        RhoMud=(((SndVshCO+SltVshCO)/2)*RhoI)+((1-((SndVshCO+SltVshCO)/2))*RhoQ);
        VpMud=(((SndVshCO+SltVshCO)/2)*VpI)+((1-((SndVshCO+SltVshCO)/2))*VpQ);
        VsMud=(((SndVshCO+SltVshCO)/2)*VsI)+((1-((SndVshCO+SltVshCO)/2))*VsQ);
    else
        RhoMud=(((SndVshCO+SltVshCO)/2)*RhoSI)+((1-((SndVshCO+SltVshCO)/2))*RhoQ);
        VpMud=(((SndVshCO+SltVshCO)/2)*VpSI)+((1-((SndVshCO+SltVshCO)/2))*VpQ);
        VsMud=(((SndVshCO+SltVshCO)/2)*VsSI)+((1-((SndVshCO+SltVshCO)/2))*VsQ);
    end
end



%% Bulk Rock Properties

[AI_Gas,VpVsR_Gas,Mu_Gas,K_Gas] = RockProps(VpGas_filt,VsGas_filt,RhoGas_filt);
[AI_Wtr,VpVsR_Wtr,Mu_Wtr,K_Wtr] = RockProps(VpWtr_filt,VsWtr_filt,RhoWtr_filt);

%% CURVE FITTING (Vertical effective stress vs Rock Property)
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
Por_trend = Por_min-((Por_min-PhiC)*exp(-betaPor*Z));

% % Similar approach for caprock however as PHIT not always present, RHOB
% % used as a proxy
dPhiTob=Pves_filt_OB;
GPhiTob=log((PhiT_filt_OB-Por_min_ob)/(Por_min_ob-PhiCsh));
mPhiTob=GPhiTob\dPhiTob;
betaPhiTob=real(-1/mPhiTob);

OB_por_pred = Por_min_ob-((Por_min_ob-PhiCsh)*exp(-betaPhiTob*(Pves_pred)));
OB_por_trend = Por_min_ob-((Por_min_ob-PhiCsh)*exp(-betaPhiTob*Z));


% Gradient matrix for all subsequent least squares depth trend equations
% Porosity depth trend used as a known compaction gradient
Gr=exp(-betaPhiTob*Pves_filt_OB); % Compaction gradient for caprock
Ggas=exp(-betaPor*Pves_filt_sand_gas); % Compaction gradient for sands
Gwtr=exp(-betaPor*Pves_filt_sand_wtr);

% Gas Sand Rock Properties
% Rho Gas
d=RhoGas_filt-RhoSd; % known data
mx=RhoSd; % Quartz density
mRhoGas=ModelFitv2(d,Ggas,betaPor,mx);
RhoGas_pred = RhoSd-((RhoSd-mRhoGas)*exp(-betaPor*Pves_pred));
RhoGas_trend = RhoSd-((RhoSd-mRhoGas)*exp(-betaPor*Z));
RhoGas_dt = RhoSd-((RhoSd-mRhoGas)*exp(-betaPor*Zrange));

% Vp Gas
d=VpGas_filt-VpSd; % known data
mx=VpSd; % Quartz Vp
mVpGas=ModelFitv2(d,Ggas,betaPor,mx);
VpGas_pred = VpSd-((VpSd-mVpGas)*exp(-betaPor*Pves_pred));
VpGas_trend = VpSd-((VpSd-mVpGas)*exp(-betaPor*Z));
VpGas_dt = VpSd-((VpSd-mVpGas)*exp(-betaPor*Zrange));

% Vs Gas
d=VsGas_filt-VsSd; % known data
mx=VsSd; % Quartz Vs
mVsGas=ModelFitv2(d,Ggas,betaPor,mx);
VsGas_pred = VsSd-((VsSd-mVsGas)*exp(-betaPor*Pves_pred));
VsGas_trend = VsSd-((VsSd-mVsGas)*exp(-betaPor*Z));
VsGas_dt = VsSd-((VsSd-mVsGas)*exp(-betaPor*Zrange));

% Water Sand Rock Properties

% Rho Water
d=RhoWtr_filt-RhoSd; % known data
mx=RhoSd; % Quartz density
mRhoWtr=ModelFitv2(d,Gwtr,betaPor,mx);
RhoWtr_pred = RhoSd-((RhoSd-mRhoWtr)*exp(-betaPor*Pves_pred));
RhoWtr_trend = RhoSd-((RhoSd-mRhoWtr)*exp(-betaPor*Z));
RhoWtr_dt = RhoSd-((RhoSd-mRhoWtr)*exp(-betaPor*Zrange));

% Vp Water
d=VpWtr_filt-VpSd; % known data
mx=VpSd; % Quartz Vp
mVpWtr=ModelFitv2(d,Gwtr,betaPor,mx);
VpWtr_pred = VpSd-((VpSd-mVpWtr)*exp(-betaPor*Pves_pred));
VpWtr_trend = VpSd-((VpSd-mVpWtr)*exp(-betaPor*Z));
VpWtr_dt = VpSd-((VpSd-mVpWtr)*exp(-betaPor*Zrange));

% Vs Water
d=VsWtr_filt-VsSd; % known data
mx=VsSd; % Quartz Vs
mVsWtr=ModelFitv2(d,Gwtr,betaPor,mx);
VsWtr_pred = VsSd-((VsSd-mVsWtr)*exp(-betaPor*Pves_pred));
VsWtr_trend = VsSd-((VsSd-mVsWtr)*exp(-betaPor*Z));
VsWtr_dt = VsSd-((VsSd-mVsWtr)*exp(-betaPor*Zrange));

% Overpressure Modelling
% POROSITY
[Porosity_predicted]=Overpressure(a_sst,b_sst,c_sst,Pves_pred,Pop,0,Porosity_predicted,PhiC,Por_min)

% GAS SAND
[VpGas_pred]=Overpressure(a_sst,b_sst,c_sst,Pves_pred,Pop,Buoy,VpGas_pred,mVsGas,VpSd);
[VsGas_pred]=Overpressure(a_sst,b_sst,c_sst,Pves_pred,Pop,Buoy,VsGas_pred,mVsGas,VsSd);
[RhoGas_pred]=Overpressure(a_sst,b_sst,c_sst,Pves_pred,Pop,Buoy,RhoGas_pred,mRhoGas,RhoSd);

[VpGas_dt]=Overpressure(a_sst,b_sst,c_sst,Pves_pred,Pop,Buoy,VpGas_dt,mVsGas,VpSd);
[VsGas_dt]=Overpressure(a_sst,b_sst,c_sst,Pves_pred,Pop,Buoy,VsGas_dt,mVsGas,VsSd);
[RhoGas_dt]=Overpressure(a_sst,b_sst,c_sst,Pves_pred,Pop,Buoy,RhoGas_dt,mRhoGas,RhoSd);

% BRINE SAND
[VpWtr_pred]=Overpressure(a_sst,b_sst,c_sst,Pves_pred,Pop,0,VpWtr_pred,mVsWtr,VpSd);
[VsWtr_pred]=Overpressure(a_sst,b_sst,c_sst,Pves_pred,Pop,0,VsWtr_pred,mVsWtr,VsSd);
[RhoWtr_pred]=Overpressure(a_sst,b_sst,c_sst,Pves_pred,Pop,0,RhoWtr_pred,mRhoWtr,RhoSd);

[VpWtr_dt]=Overpressure(a_sst,b_sst,c_sst,Pves_pred,Pop,0,VpWtr_dt,mVsWtr,VpSd);
[VsWtr_dt]=Overpressure(a_sst,b_sst,c_sst,Pves_pred,Pop,0,VsWtr_dt,mVsWtr,VsSd);
[RhoWtr_dt]=Overpressure(a_sst,b_sst,c_sst,Pves_pred,Pop,0,RhoWtr_dt,mRhoWtr,RhoSd);
% Caprock property fitting procedure:

% Filtered overburden rock properties

d=RhoOB_filt-RhoMud; % known data
mx=RhoMud; % Quartz density
mRhoMud=ModelFitv2(d,Gr,betaPhiTob,mx);
RhoOB_pred = RhoMud-((RhoMud-mRhoMud)*exp(-betaPhiTob*Pves_pred));
RhoOB_trend = RhoMud-((RhoMud-mRhoMud)*exp(-betaPhiTob*Z));
RhoOB_dt = RhoMud-((RhoMud-mRhoMud)*exp(-betaPhiTob*Zrange));

d=VpOB_filt-VpMud; % known data
mx=VpMud; % Quartz density
mVpMud=ModelFitv2(d,Gr,betaPhiTob,mx);
VpOB_pred = VpMud-((VpMud-mVpMud)*exp(-betaPhiTob*Pves_pred));
VpOB_trend = VpMud-((VpMud-mVpMud)*exp(-betaPhiTob*Z));
VpOB_dt = VpMud-((VpMud-mVpMud)*exp(-betaPhiTob*Zrange));

d=VsOB_filt-VsMud; % known data
mx=VsMud; % Quartz density
mVsMud=ModelFitv2(d,Gr,betaPhiTob,mx);
VsOB_pred = VsMud-((VsMud-mVsMud)*exp(-betaPhiTob*Pves_pred));
VsOB_trend = VsMud-((VsMud-mVsMud)*exp(-betaPhiTob*Z));
VsOB_dt = VsMud-((VsMud-mVsMud)*exp(-betaPhiTob*Zrange));

% Caprock overpressure modelling
% TZA data at Jodari suggests density measurements not effected by
% overpressure.  This corroborates published research and datasets.
% Therefore for the caprock we assume Rho_pred(Pves)=Rho_pred(Pves-Pover)

[VpOB_pred]=Overpressure(a_mud,b_mud,c_mud,Pves_pred,Pop,0,VpOB_pred,mVpMud,VpMud);
[VsOB_pred]=Overpressure(a_mud,b_mud,c_mud,Pves_pred,Pop,0,VsOB_pred,mVsMud,VsMud);

[VpOB_dt]=Overpressure(a_mud,b_mud,c_mud,Pves_pred,Pop,0,VpOB_dt,mVpMud,VpMud);
[VsOB_dt]=Overpressure(a_mud,b_mud,c_mud,Pves_pred,Pop,0,VsOB_dt,mVsMud,VsMud);

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
RhoShl_dt = RhoMud-((RhoMud-mRhoMud)*exp(-betaShPhiT*Zrange));

d=VpShl_filt-VpMud; % known data
mx=VpMud; % Quartz density
mVpMud=ModelFitv2(d,GSh,betaShPhiT,mx);
VpShl_pred = VpMud-((VpMud-mVpMud)*exp(-betaShPhiT*Pves_pred));
VpShl_trend = VpMud-((VpMud-mVpMud)*exp(-betaShPhiT*Z));
VpShl_dt = VpMud-((VpMud-mVpMud)*exp(-betaShPhiT*Zrange));

d=VsShl_filt-VsMud; % known data
mx=VsMud; % Quartz density
mVsMud=ModelFitv2(d,GSh,betaShPhiT,mx);
VsShl_pred = VsMud-((VsMud-mVsMud)*exp(-betaShPhiT*Pves_pred));
VsShl_trend = VsMud-((VsMud-mVsMud)*exp(-betaShPhiT*Z));
VsShl_dt = VpMud-((VsMud-mVsMud)*exp(-betaShPhiT*Zrange));

% Overpressure modelling
[VpShl_pred]=Overpressure(a_mud,b_mud,c_mud,Pves_pred,Pop,0,VpShl_pred,mVpMud,VpMud);
[VsShl_pred]=Overpressure(a_mud,b_mud,c_mud,Pves_pred,Pop,0,VsShl_pred,mVsMud,VsMud);

[VpShl_dt]=Overpressure(a_mud,b_mud,c_mud,Pves_pred,Pop,0,VpShl_dt,mVpMud,VpMud);
[VsShl_dt]=Overpressure(a_mud,b_mud,c_mud,Pves_pred,Pop,0,VsShl_dt,mVsMud,VsMud);

if GA >= COage
    RhoMud=(((SndVshCO+SltVshCO)/2)*RhoI)+((1-((SndVshCO+SltVshCO)/2))*RhoSd);
    VpMud=(((SndVshCO+SltVshCO)/2)*VpI)+((1-((SndVshCO+SltVshCO)/2))*VpSd);
    VsMud=(((SndVshCO+SltVshCO)/2)*VsI)+((1-((SndVshCO+SltVshCO)/2))*VsSd);
else
    RhoMud=(((SndVshCO+SltVshCO)/2)*RhoSI)+((1-((SndVshCO+SltVshCO)/2))*RhoSd);
    VpMud=(((SndVshCO+SltVshCO)/2)*VpSI)+((1-((SndVshCO+SltVshCO)/2))*VpSd);
    VsMud=(((SndVshCO+SltVshCO)/2)*VsSI)+((1-((SndVshCO+SltVshCO)/2))*VsSd);
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
RhoSlt_dt = RhoMud-((RhoMud-mRhoMud)*exp(-betaSlPhiT*Zrange));

d=VpSlt_filt-VpMud; % known data
mx=VpMud; % Quartz density
mVpMud=ModelFitv2(d,GSl,betaSlPhiT,mx);
VpSlt_pred = VpMud-((VpMud-mVpMud)*exp(-betaSlPhiT*Pves_pred));
VpSlt_trend = VpMud-((VpMud-mVpMud)*exp(-betaSlPhiT*Z));
VpSlt_dt = VpMud-((VpMud-mVpMud)*exp(-betaSlPhiT*Zrange));

d=VsSlt_filt-VsMud; % known data
mx=VsMud; % Quartz density
mVsMud=ModelFitv2(d,GSl,betaSlPhiT,mx);
VsSlt_pred = VsMud-((VsMud-mVsMud)*exp(-betaSlPhiT*Pves_pred));
VsSlt_trend = VsMud-((VsMud-mVsMud)*exp(-betaSlPhiT*Z));
VsSlt_dt = VsMud-((VsMud-mVsMud)*exp(-betaSlPhiT*Zrange));

% Overpressure modelling
[VpSlt_pred]=Overpressure(a_mud,b_mud,c_mud,Pves_pred,Pop,0,VpSlt_pred,mVpMud,VpMud);
[VsSlt_pred]=Overpressure(a_mud,b_mud,c_mud,Pves_pred,Pop,0,VsSlt_pred,mVsMud,VsMud);

[VpSlt_dt]=Overpressure(a_mud,b_mud,c_mud,Pves_pred,Pop,0,VpSlt_dt,mVpMud,VpMud);
[VsSlt_dt]=Overpressure(a_mud,b_mud,c_mud,Pves_pred,Pop,0,VsSlt_dt,mVsMud,VsMud);

%% PREDICTED ROCK PROPERTIES AND AVA RESPONSE

% Predicted Bulk rock properties 
[AI_Gas_pred,VpVsR_Gas_pred,Mu_Gas_pred,K_Gas_pred] = RockProps(VpGas_pred,VsGas_pred,RhoGas_pred);
[AI_Wtr_pred,VpVsR_Wtr_pred,Mu_Wtr_pred,K_Wtr_pred] = RockProps(VpWtr_pred,VsWtr_pred,RhoWtr_pred);
[AI_Shl_pred,VpVsR_Shl_pred,Mu_Shl_pred,K_Shl_pred] = RockProps(VpShl_pred,VsShl_pred,RhoShl_pred);
[AI_Slt_pred,VpVsR_Slt_pred,Mu_Slt_pred,K_Slt_pred] = RockProps(VpSlt_pred,VsSlt_pred,RhoSlt_pred);

% Angle dependent reflectivity based on depth trend fitting
[RPPpred_gas,RPSpred_gas,TPPpred_gas,TPSPpred_gas]=Zoeppritz(MaxAngle,VpOB_pred,VsOB_pred,RhoOB_pred,VpGas_pred,VsGas_pred,RhoGas_pred);
[RPPpred_wtr,RPSpred_wtr,TPPpred_wtr,TPSPpred_wtr]=Zoeppritz(MaxAngle,VpOB_pred,VsOB_pred,RhoOB_pred,VpWtr_pred,VsWtr_pred,RhoWtr_pred);
[RPPpred_gwc,RPSpred_gwc,TPPpred_gwc,TPSPpred_gwc]=Zoeppritz(MaxAngle,VpGas_pred,VsGas_pred,RhoGas_pred,VpWtr_pred,VsWtr_pred,RhoWtr_pred);
[RPPpred_shlslt,RPSpred_shlslt,TPPpred_shlslt,TPSPpred_shlslt]=Zoeppritz(MaxAngle,VpShl_pred,VsShl_pred,RhoShl_pred,VpSlt_pred,VsSlt_pred,RhoSlt_pred);
[RPPpred_sltshl,RPSpred_sltshl,TPPpred_sltshl,TPSPpred_sltshl]=Zoeppritz(MaxAngle,VpSlt_pred,VsSlt_pred,RhoSlt_pred,VpShl_pred,VsShl_pred,RhoShl_pred);

% Calculation of approximation to Zoeppritz equations
[A_Pred_gas_iso,B_Pred_gas_iso,C_Pred_gas_iso,A_Pred_gas_aniso,B_Pred_gas_aniso,C_Pred_gas_aniso,RPPshuey_Pred_gas_iso,RPPshuey_Pred_gas_aniso] = Shuey(MaxAngle,VpOB_pred,VsOB_pred,RhoOB_pred,VpGas_pred,VsGas_pred,RhoGas_pred,delta1,delta2,epsilon1,epsilon2);
[A_Pred_wtr_iso,B_Pred_wtr_iso,C_Pred_wtr_iso,A_Pred_wtr_aniso,B_Pred_wtr_aniso,C_Pred_wtr_aniso,RPPshuey_Pred_wtr_iso,RPPshuey_Pred_wtr_aniso] = Shuey(MaxAngle,VpOB_pred,VsOB_pred,RhoOB_pred,VpWtr_pred,VsWtr_pred,RhoWtr_pred,delta1,delta2,epsilon1,epsilon2);
[A_Pred_gwc_iso,B_Pred_gwc_iso,C_Pred_gwc_iso,A_Pred_gwc_aniso,B_Pred_gwc_aniso,C_Pred_gwc_aniso,RPPshuey_Pred_gwc_iso,RPPshuey_Pred_gwc_aniso] = Shuey(MaxAngle,VpGas_pred,VsGas_pred,RhoGas_pred,VpWtr_pred,VsWtr_pred,RhoWtr_pred,0,0,0,0);
[A_Pred_shlslt_iso,B_Pred_shlslt_iso,C_Pred_shlslt_iso,A_Pred_shlslt_aniso,B_Pred_shlslt_aniso,C_Pred_shlslt_aniso,RPPshuey_Pred_shlslt_iso,RPPshuey_Pred_shlslt_aniso] = Shuey(MaxAngle,VpShl_pred,VsShl_pred,RhoShl_pred,VpSlt_pred,VsSlt_pred,RhoSlt_pred,delta1,delta2,epsilon1,epsilon2);
[A_Pred_sltshl_iso,B_Pred_sltshl_iso,C_Pred_sltshl_iso,A_Pred_sltshl_aniso,B_Pred_sltshl_aniso,C_Pred_sltshl_aniso,RPPshuey_Pred_sltshl_iso,RPPshuey_Pred_sltshl_aniso] = Shuey(MaxAngle,VpSlt_pred,VsSlt_pred,RhoSlt_pred,VpShl_pred,VsShl_pred,RhoShl_pred,delta1,delta2,epsilon1,epsilon2);

% Full Stack Response
% Mean of Zoeppritz reflectivity with incidence angle over full angle range
FS_Gas(1:MaxAngle+1,1)=mean(RPPpred_gas); % Gas
FS_Wtr(1:MaxAngle+1,1)=mean(RPPpred_wtr); % Brine

% Estimated depth trends
[A_dt_gas_iso,B_dt_gas_iso,C_dt_gas_iso,A_dt_gas_aniso,B_dt_gas_aniso,C_dt_gas_aniso,RPPshuey_dt_gas_iso,RPPshuey_dt_gas_aniso] = Shuey(MaxAngle,VpOB_dt,VsOB_dt,RhoOB_dt,VpGas_dt,VsGas_dt,RhoGas_dt,delta1,delta2,epsilon1,epsilon2);
[A_dt_wtr_iso,B_dt_wtr_iso,C_dt_wtr_iso,A_dt_wtr_aniso,B_dt_wtr_aniso,C_dt_wtr_aniso,RPPshuey_dt_wtr_iso,RPPshuey_dt_wtr_aniso] = Shuey(MaxAngle,VpOB_dt,VsOB_dt,RhoOB_dt,VpWtr_dt,VsWtr_dt,RhoWtr_dt,delta1,delta2,epsilon1,epsilon2);



%% ROCK PHYSICS MODELS
%EDITABLE CONTENT SEE LINE 39
kf=1e9*Kb;
rhof=1000*RhoW_pred;
rhos=1000*RhoShl_pred;
Ks=K_Shl_pred*1000;
Mus=Mu_Shl_pred*1000;

% Hertz Mindlin Soft Sediment model
[Phi,VpDry,VsDry,RHOBDry,VpSat,VsSat,RHOBSat]=softsediments(Clay,Feldspar,Calcite,PVES,PhiC,Coord,kf,rhof,Ks,Mus,rhos,shearfact);
AISat=VpSat.*(RHOBSat/1000);
VpVs_R_Sat=VpSat./VsSat;

% Cemented Sand Model
[Phic,Rhodc,Vpdc,Vsdc,Rhosc,Vpsc,Vssc,Phis,Rhods,Vpds,Vsds,Rhoss,Vpss,Vsss]=cementedsand(Clay,Feldspar,Calcite,PhiC,fcement,kf,rhof,shearfact);
AIsc=Vpsc.*(Rhosc/1000);
VpVs_R_sc=Vpsc./Vssc;
AIss=Vpss.*(Rhoss/1000);
VpVs_R_ss=Vpss./Vsss;


%% MONTE CARLO SIMULATION

% Calculate residuals between model and filtered data points for Vp, Vs and Rho. Residuals are used as the uncertainty distribution as this removes the compaction effect.
% Generate a normal random distribution of rock properties based on the covariance between the residuals
[VpOB_r VsOB_r RhoOB_r VpGas_r VsGas_r RhoGas_r Por_r_ob Por_r_gas] = RockPropDist(NI,Dist,MaxZ,Pves_filt_OB,Pves_filt_sand_gas,VpOB_pred,VsOB_pred,RhoOB_pred,VpGas_pred,VsGas_pred,RhoGas_pred,VpOB_filt,VsOB_filt,RhoOB_filt,VpGas_filt,VsGas_filt,RhoGas_filt,VpOB_trend,VsOB_trend,RhoOB_trend,VpGas_trend,VsGas_trend,RhoGas_trend,PhiT_filt_OB,OB_por_pred,OB_por_trend,Por_filt_gas,Porosity_predicted,Por_trend);
[VpOB_r VsOB_r RhoOB_r VpWtr_r VsWtr_r RhoWtr_r Por_r_ob Por_r_wtr] = RockPropDist(NI,Dist,MaxZ,Pves_filt_OB,Pves_filt_sand_wtr,VpOB_pred,VsOB_pred,RhoOB_pred,VpWtr_pred,VsWtr_pred,RhoWtr_pred,VpOB_filt,VsOB_filt,RhoOB_filt,VpWtr_filt,VsWtr_filt,RhoWtr_filt,VpOB_trend,VsOB_trend,RhoOB_trend,VpWtr_trend,VsWtr_trend,RhoWtr_trend,PhiT_filt_OB,OB_por_pred,OB_por_trend,Por_filt_wtr,Porosity_predicted,Por_trend);
[VpSlt_r VsSlt_r RhoSlt_r VpShl_r VsShl_r RhoShl_r Por_r_ob1 Por_r_ob2] = RockPropDist(NI,Dist,MaxZ,Pves_filt_slt,Pves_filt_shl,VpSlt_pred,VsSlt_pred,RhoSlt_pred,VpShl_pred,VsShl_pred,RhoShl_pred,VpSlt_filt,VsSlt_filt,RhoSlt_filt,VpShl_filt,VsShl_filt,RhoShl_filt,VpSlt_trend,VsSlt_trend,RhoSlt_trend,VpShl_trend,VsShl_trend,RhoShl_trend,PhiT_filt_slt,Slt_por_pred,Slt_por_trend,PhiT_filt_shl,Shl_por_pred,Shl_por_trend);

% Shuey A, B and C
[A_Rand_gas_iso,B_Rand_gas_iso,C_Rand_gas_iso,A_Rand_gas_aniso,B_Rand_gas_aniso,C_Rand_gas_aniso,RPPshuey_Rand_gas_iso,RPPshuey_Rand_gas_aniso] = Shuey(MaxAngle,VpOB_r,VsOB_r,RhoOB_r,VpGas_r,VsGas_r,RhoGas_r,delta1,delta2,epsilon1,epsilon2);
[A_Rand_wtr_iso,B_Rand_wtr_iso,C_Rand_wtr_iso,A_Rand_wtr_aniso,B_Rand_wtr_aniso,C_Rand_wtr_aniso,RPPshuey_Rand_wtr_iso,RPPshuey_Rand_wtr_aniso] = Shuey(MaxAngle,VpOB_r,VsOB_r,RhoOB_r,VpWtr_r,VsWtr_r,RhoWtr_r,delta1,delta2,epsilon1,epsilon2);
[A_Rand_shlslt_iso,B_Rand_shlslt_iso,C_Rand_shlslt_iso,A_Rand_shlslt_aniso,B_Rand_shlslt_aniso,C_Rand_shlslt_aniso,RPPshuey_Rand_shlslt_iso,RPPshuey_Rand_shlslt_aniso] = Shuey(MaxAngle,VpShl_r,VsShl_r,RhoShl_r,VpSlt_r,VsSlt_r,RhoSlt_r,delta1,delta2,epsilon1,epsilon2);

A_Rand_bg_iso=vertcat(-A_Rand_shlslt_iso',A_Rand_shlslt_iso');
B_Rand_bg_iso=vertcat(-B_Rand_shlslt_iso',B_Rand_shlslt_iso');

% Rock Property Distribution
[AI_Gas_r,VpVsR_Gas_r,Mu_Gas_r,K_Gas_r] = RockProps(VpGas_r,VsGas_r,RhoGas_r);
[AI_Wtr_r,VpVsR_Wtr_r,Mu_Wtr_r,K_Wtr_r] = RockProps(VpWtr_r,VsWtr_r,RhoWtr_r);
[AI_Shl_r,VpVsR_Shl_r,Mu_Shl_r,K_Shl_r] = RockProps(VpShl_r,VsShl_r,RhoShl_r);
[AI_Slt_r,VpVsR_Slt_r,Mu_Slt_r,K_Slt_r] = RockProps(VpSlt_r,VsSlt_r,RhoSlt_r);
Lith1(1:length(AI_Gas_r),1)=1;
Lith2(1:length(AI_Gas_r),1)=2;
Lith3(1:length(AI_Gas_r),1)=3;
Lith4(1:length(AI_Gas_r),1)=4;
Lith=vertcat(Lith1,Lith2,Lith3,Lith4);
AI_r=vertcat(AI_Wtr_r,AI_Shl_r,AI_Gas_r,AI_Slt_r);
VpVsR_r=vertcat(VpVsR_Wtr_r,VpVsR_Shl_r,VpVsR_Gas_r,VpVsR_Slt_r);

% Chi - Minimum Energy Angle
% Minimum energy angle calculated based on monte carlo shale-silt interface
I0=max(A_Rand_shlslt_iso); % find the max intercept value
G0=min(B_Rand_shlslt_iso); % find the min gradient value
Chi_min_r=atan((-I0*Iscalar)/(G0*Gscalar)); % use cluster extreme to fit trend line through data points
Chi_min=180*Chi_min_r/pi

% Chi - Fluid Angle
% Minimum energy angle calculated based on monte carlo caprock-brine sand interface
I1=max(A_Rand_wtr_iso);
G1=min(B_Rand_wtr_iso);
Chi_minl_r=atan((-I1*Iscalar)/(G1*Gscalar));
Chi_minl=180*Chi_minl_r/pi

% Maximum energy angle calcualted based on predicted GWC interface
Chi_max_r=atan((-A_Pred_gwc_iso*Iscalar)/(B_Pred_gwc_iso*Gscalar));
Chi_max=180*Chi_max_r/pi

% Predicted Extended Elastic Reflectivity calculation based on calculated Minimum Energy Chi angle
MinEERgas=A_Pred_gas_iso*Iscalar*cos(Chi_min_r)+B_Pred_gas_iso*Gscalar*sin(Chi_min_r);
MinEERwtr=A_Pred_wtr_iso*Iscalar*cos(Chi_min_r)+B_Pred_wtr_iso*Gscalar*sin(Chi_min_r);
MinEERgwc=A_Pred_gwc_iso*Iscalar*cos(Chi_min_r)+B_Pred_gwc_iso*Gscalar*sin(Chi_min_r);

% Monte Carlo Extended Elastic Reflectivity calculation based on calculated Minimum Energy Chi angle
MinEERgas_rnd=A_Rand_gas_iso.*Iscalar.*cos(Chi_min_r)+B_Rand_gas_iso.*Gscalar.*sin(Chi_min_r);
MinEERwtr_rnd=A_Rand_wtr_iso.*Iscalar.*cos(Chi_min_r)+B_Rand_wtr_iso.*Gscalar.*sin(Chi_min_r);
MinEERbg_rnd=A_Rand_bg_iso.*Iscalar.*cos(Chi_min_r)+B_Rand_bg_iso.*Gscalar.*sin(Chi_min_r);
 
% Monte Carlo Extended Elastic Reflectivity calculation based on calculated Fluid Chi angle
MinlEERgas_rnd=A_Rand_gas_iso.*Iscalar.*cos(Chi_minl_r)+B_Rand_gas_iso.*Gscalar.*sin(Chi_minl_r);
MinlEERwtr_rnd=A_Rand_wtr_iso.*Iscalar.*cos(Chi_minl_r)+B_Rand_wtr_iso.*Gscalar.*sin(Chi_minl_r);
MinlEERbg_rnd=A_Rand_bg_iso.*Iscalar.*cos(Chi_minl_r)+B_Rand_bg_iso.*Gscalar.*sin(Chi_minl_r);

% Fluid Seperation
FluidSep=mean(MinEERgas_rnd)-mean(MinEERwtr_rnd);
GasSep=mean(MinEERgas_rnd)-mean(MinEERbg_rnd);
WtrSep=mean(MinEERwtr_rnd)-mean(MinEERbg_rnd);
FluidDiscrim=FluidSep/WtrSep

Vp_gas_min=min(VpGas_r)
Vp_gas_mid=mean(VpGas_r)
Vp_gas_max=max(VpGas_r)

%% PLOTTING INPUTS

% AXES LIMITS
% I-G 
MinI=-750;
MaxI=750;
MinG=-1500;
MaxG=1500;
Step=250;

% Reflectivity vs Angle
MinR=-0.3;
MaxR=0.3;

% Rock properties vs Vertical Effective Stress
MinVes=0;
MaxVes=MaxZ;
MinRho=1.7;
MaxRho=2.7;
MinVp=1500;
MaxVp=5500;
MinVs=500;
MaxVs=3500;

% Background trend
BGmin=[MinI -MinI/tan(Chi_min_r); MaxI -MaxI/tan(Chi_min_r)];
BGminl=[MinI -MinI/tan(Chi_minl_r); MaxI -MaxI/tan(Chi_minl_r)];
BGmax=[MinI -MinI/tan(Chi_max_r); MaxI -MaxI/tan(Chi_max_r)];

YAxis=[0 MaxG*Gscalar; 0 MinG*Gscalar];
XAxis=[MinI*Iscalar 0; MaxI*Iscalar 0];

% Plot positions
Width=750; %double monitor 1000
LengthSh=500; %double monitor 725
LengthLg=1000; %double monitor 1600
%Width=1000; %double monitor 1000
%LengthSh=725; %double monitor 725
%LengthLg=1600; %double monitor 1600


%% FIGURES
%% FIGURE 1 I-G crossplot showing depth trends for modelled facies
figure1 = figure;
set(figure1, 'Position', [0 LengthSh+(0.1*LengthSh) Width LengthSh])
% Create subplot
subplot1 = subplot(1,1,1,'Parent',figure1,'XTick',[MinI:Step:MaxI],'YTick',[MinG:Step:MaxG]);
box(subplot1,'on');
hold(subplot1,'all');
xlim(subplot1,[MinI MaxI]);
ylim(subplot1,[MinG MaxG]);
 
% Create plot
scatter(A_dt_gas_iso*Iscalar,B_dt_gas_iso*Gscalar,200,Zrange,'o','filled','Parent',subplot1,'LineWidth',4,'MarkerEdgeColor','red');
scatter(A_dt_wtr_iso*Iscalar,B_dt_wtr_iso*Gscalar,200,Zrange,'o','filled','Parent',subplot1,'LineWidth',4,'MarkerEdgeColor','blue');
 
% Plot Axes
plot(XAxis(:,1),XAxis(:,2),'Parent',subplot1,'LineWidth',2,'color','black');
plot(YAxis(:,1),YAxis(:,2),'Parent',subplot1,'LineWidth',2,'color','black');

% Create title
title('Intercept v Gradient (Coloured by Effective Stress)');
 

%% FIGURE 2 Reservoir rock properties v Effective Stress coloured by saturation
figure2 = figure;
set(figure2, 'Position', [2*Width 0 0.75*Width LengthLg])
% Create subplot

% Create subplot
subplot1 = subplot(1,3,1,'Parent',figure2,'Ydir','reverse','XTick',[MinRho:0.25:MaxRho],'YTick',[MinVes:500:MaxVes]);
box(subplot1,'on');
hold(subplot1,'all');
xlim(subplot1,[MinRho MaxRho]);
ylim(subplot1,[MinVes MaxVes]);

% Create plot
plot(RhoGas_filt,Pves_filt_sand_gas,'.','MarkerSize',5,'Parent',subplot1,'color','red');
plot(RhoWtr_filt,Pves_filt_sand_wtr,'.','MarkerSize',5,'Parent',subplot1,'color','blue');
plot(RhoGas_trend,Z,'--','Parent',subplot1,'color','red','LineWidth',2);
plot(RhoWtr_trend,Z,'--','Parent',subplot1,'color','blue','LineWidth',2);
plot(RhoGas_pred,Pves_pred-Pop-Buoy,'o','MarkerSize',10,'Parent',subplot1,'MarkerFaceColor','black','LineWidth',3,'MarkerEdgeColor','red');
plot(RhoWtr_pred,Pves_pred-Pop,'o','MarkerSize',10,'Parent',subplot1,'MarkerFaceColor','black','LineWidth',3,'MarkerEdgeColor','blue');

% Create title
title('Reservoir density vs Effective Stress');

% Create subplot
subplot2 = subplot(1,3,2,'Parent',figure2,'Ydir','reverse','XTick',[MinVp:1000:MaxVp],'YTick',[MinVes:500:MaxVes]);
box(subplot2,'on');
hold(subplot2,'all');
xlim(subplot2,[MinVp MaxVp]);
ylim(subplot2,[MinVes MaxVes]);

% Create plot
plot(VpGas_filt,Pves_filt_sand_gas,'.','MarkerSize',5,'Parent',subplot2,'color','red');
plot(VpWtr_filt,Pves_filt_sand_wtr,'.','MarkerSize',5,'Parent',subplot2,'color','blue');
plot(VpGas_trend,Z,'--','Parent',subplot2,'color','red','LineWidth',2);
plot(VpWtr_trend,Z,'--','Parent',subplot2,'color','blue','LineWidth',2);
plot(VpGas_pred,Pves_pred-Pop-Buoy,'o','MarkerSize',10,'Parent',subplot2,'MarkerFaceColor','black','LineWidth',3,'MarkerEdgeColor','red');
plot(VpWtr_pred,Pves_pred-Pop,'o','MarkerSize',10,'Parent',subplot2,'MarkerFaceColor','black','LineWidth',3,'MarkerEdgeColor','blue');

% Create title
title('Vp reservoir v Effective Stress');

% Create subplot
subplot3 = subplot(1,3,3,'Parent',figure2,'Ydir','reverse','XTick',[MinVs:1000:MaxVs],'YTick',[MinVes:500:MaxVes]);
box(subplot3,'on');
hold(subplot3,'all');
xlim(subplot3,[MinVs MaxVs]);
ylim(subplot3,[MinVes MaxVes]);

% Create plot
plot(VsGas_filt,Pves_filt_sand_gas,'.','MarkerSize',5,'Parent',subplot3,'color','red');
plot(VsWtr_filt,Pves_filt_sand_wtr,'.','MarkerSize',5,'Parent',subplot3,'color','blue');
plot(VsGas_trend,Z,'--','Parent',subplot3,'color','red','LineWidth',2);
plot(VsWtr_trend,Z,'--','Parent',subplot3,'color','blue','LineWidth',2);
plot(VsGas_pred,Pves_pred-Pop-Buoy,'o','MarkerSize',10,'Parent',subplot3,'MarkerFaceColor','black','LineWidth',3,'MarkerEdgeColor','red');
plot(VsWtr_pred,Pves_pred-Pop,'o','MarkerSize',10,'Parent',subplot3,'MarkerFaceColor','black','LineWidth',3,'MarkerEdgeColor','blue');

% Create title
title('Vs reservoir v Effective Stress');
   

%% FIGURE 3 Reservoir rock properties v Effective Stress coloured by age
figure3 = figure;
set(figure3, 'Position', [2.75*Width 0 0.75*Width LengthLg])
% Create subplot
subplot1 = subplot(1,3,1,'Parent',figure3,'Ydir','reverse','XTick',[MinRho:0.25:MaxRho],'YTick',[MinVes:500:MaxVes]);
box(subplot1,'on');
hold(subplot1,'all');
xlim(subplot1,[MinRho MaxRho]);
ylim(subplot1,[MinVes MaxVes]);

% Create plot
scatter(RhoGas_filt,Pves_filt_sand_gas,20,Age_filt_gas,'o','filled','Parent',subplot1);
plot(RhoGas_trend,Z,'--','Parent',subplot1,'color','red','LineWidth',2);
plot(RhoGas_pred,Pves_pred-Pop-Buoy,'o','MarkerSize',10,'Parent',subplot1,'MarkerFaceColor','black','LineWidth',4,'MarkerEdgeColor','red');


% Create title
title('Reservoir density vs Effective Stress');

% Create subplot
subplot2 = subplot(1,3,2,'Parent',figure3,'Ydir','reverse','XTick',[MinVp:1000:MaxVp],'YTick',[MinVes:500:MaxVes]);
box(subplot2,'on');
hold(subplot2,'all');
xlim(subplot2,[MinVp MaxVp]);
ylim(subplot2,[MinVes MaxVes]);

% Create plot
%scatter(VpGas_filt,Pves_filt_sand_gas,100,'.','Parent',subplot2);
scatter(VpGas_filt,Pves_filt_sand_gas,20,Age_filt_gas,'o','filled','Parent',subplot2);
plot(VpGas_trend,Z,'--','Parent',subplot2,'color','red','LineWidth',2);
plot(VpGas_pred,Pves_pred-Pop-Buoy,'o','MarkerSize',10,'Parent',subplot2,'MarkerFaceColor','black','LineWidth',4,'MarkerEdgeColor','red');


% Create title
title('Vp reservoir v Effective Stress');

% Create subplot
subplot3 = subplot(1,3,3,'Parent',figure3,'Ydir','reverse','XTick',[MinVs:1000:MaxVs],'YTick',[MinVes:500:MaxVes]);
box(subplot3,'on');
hold(subplot3,'all');
xlim(subplot3,[MinVs MaxVs]);
ylim(subplot3,[MinVes MaxVes]);

% Create plot
%scatter(VsGas_filt,Pves_filt_sand_gas,100,'.','Parent',subplot3);
scatter(VsGas_filt,Pves_filt_sand_gas,20,Age_filt_gas,'o','filled','Parent',subplot3);
plot(VsGas_trend,Z,'--','Parent',subplot3,'color','red','LineWidth',2);
plot(VsGas_pred,Pves_pred-Pop-Buoy,'o','MarkerSize',10,'Parent',subplot3,'MarkerFaceColor','black','LineWidth',4,'MarkerEdgeColor','red');


% Create title
title('Vs reservoir v Effective Stress');

%% FIGURE 4 Overburden Porosity v Effective Stress
figure4 = figure;
set(figure4, 'Position', [3.5*Width 0 0.75*Width LengthLg])

% Create subplot
subplot1 = subplot(1,3,1,'Parent',figure4,'Ydir','reverse','XTick',[2.2:0.2:2.8],'YTick',[MinVes:500:MaxVes]);
box(subplot1,'on');
hold(subplot1,'all');
xlim(subplot1,[2.2 2.8]);
ylim(subplot1,[MinVes MaxVes]);

% Create plot
scatter(RhoShl_filt,Pves_filt_shl,100,Vsh_filt_shl,'.','Parent',subplot1);
scatter(RhoSlt_filt,Pves_filt_slt,100,Vsh_filt_slt,'.','Parent',subplot1);
plot(RhoOB_trend,Z,'--','Parent',subplot1);
plot(RhoOB_pred,Pves_pred-Pop,'o','MarkerSize',10,'Parent',subplot1,'MarkerFaceColor','black','LineWidth',4,'MarkerEdgeColor','green');

% Create title
title('Overburden density v Effective Stress');


% Create subplot
subplot2 = subplot(1,3,2,'Parent',figure4,'Ydir','reverse','XTick',[MinVp:1000:MaxVp],'YTick',[MinVes:500:MaxVes]);
box(subplot2,'on');
hold(subplot2,'all');
xlim(subplot2,[MinVp MaxVp]);
ylim(subplot2,[MinVes MaxVes]);

% Create plot
scatter(VpShl_filt,Pves_filt_shl,100,Vsh_filt_shl,'.','Parent',subplot2);
scatter(VpSlt_filt,Pves_filt_slt,100,Vsh_filt_slt,'.','Parent',subplot2);
plot(VpOB_trend,Z,'--','Parent',subplot2);
plot(VpOB_pred,Pves_pred-Pop,'o','MarkerSize',10,'Parent',subplot2,'MarkerFaceColor','black','LineWidth',4,'MarkerEdgeColor','green');

% Create title
title('Vp overburden v Effective Stress');

% Create subplot
subplot3 = subplot(1,3,3,'Parent',figure4,'Ydir','reverse','XTick',[MinVs:1000:MaxVs],'YTick',[MinVes:500:MaxVes]);
box(subplot3,'on');
hold(subplot3,'all');
xlim(subplot3,[MinVs MaxVs]);
ylim(subplot3,[MinVes MaxVes]);

% Create plot
scatter(VsShl_filt,Pves_filt_shl,100,Vsh_filt_shl,'.','Parent',subplot3);
scatter(VsSlt_filt,Pves_filt_slt,100,Vsh_filt_slt,'.','Parent',subplot3);
plot(VsOB_trend,Z,'--','Parent',subplot3);
plot(VsOB_pred,Pves_pred-Pop,'o','MarkerSize',10,'Parent',subplot3,'MarkerFaceColor','black','LineWidth',4,'MarkerEdgeColor','green');


% Create title
title('Vs overburden v Effective Stress');



%% FIGURE 5 Reservoir Porosity v Effective Stress
figure5 = figure;
set(figure5, 'Position', [4.25*Width 0 0.75*Width LengthLg])

% Create subplot
subplot1 = subplot(1,2,1,'Parent',figure5,'Ydir','reverse','XTick',[0:0.05:0.5],'YTick',[MinVes:500:MaxVes]);
box(subplot1,'on');
hold(subplot1,'all');
xlim(subplot1,[0 0.5]);
ylim(subplot1,[MinVes MaxVes]);

% Create plot
scatter(Por_filt,Pves_filt_sand,20,vsh_filt,'o','filled','Parent',subplot1);
plot(Por_trend,Z,'Parent',subplot1);
plot(Porosity_predicted,Pves_pred-Pop,'o','MarkerSize',10,'Parent',subplot1,'MarkerFaceColor','red','LineWidth',4,'MarkerEdgeColor','black');


% Create title
title('Reservoir Porosity v Effective Stress');

% Create subplot
subplot2 = subplot(1,2,2,'Parent',figure5,'Ydir','reverse','XTick',[0:0.05:0.7],'YTick',[MinVes:500:MaxVes]);
box(subplot2,'on');
hold(subplot2,'all');
xlim(subplot2,[0 0.7]);
ylim(subplot2,[MinVes MaxVes]);

% Create plot
scatter(PhiT_filt_shl,Pves_filt_shl,20,Age_filt_shl,'o','filled','Parent',subplot2);
scatter(PhiT_filt_slt,Pves_filt_slt,20,Age_filt_slt,'o','filled','Parent',subplot2);
plot(OB_por_trend,Z,'Parent',subplot2);
plot(OB_por_pred,Pves_pred-Pop,'o','MarkerSize',10,'Parent',subplot2,'MarkerFaceColor','red','LineWidth',4,'MarkerEdgeColor','black');


% Create title
title('Overburden Porosity v Effective Stress');

%% FIGURE 6 I-G interface models
figure6 = figure;
set(figure6, 'Position', [0 0 2*Width LengthSh])
% Create subplot
subplot1 = subplot(1,2,1,'Parent',figure6,'XTick',[MinI:Step:MaxI],'YTick',[MinG:Step:MaxG]);
box(subplot1,'on');
hold(subplot1,'all');
xlim(subplot1,[MinI MaxI]);
ylim(subplot1,[MinG MaxG]);

% Plot predicted I and G values based on user input
plot(A_Rand_bg_iso*Iscalar,B_Rand_bg_iso*Gscalar,'o','MarkerSize',2,'Parent',subplot1,'MarkerFaceColor','green','LineWidth',1,'MarkerEdgeColor','green');
plot(A_Rand_gas_iso*Iscalar,B_Rand_gas_iso*Gscalar,'o','MarkerSize',2,'Parent',subplot1,'MarkerFaceColor','red','LineWidth',1,'MarkerEdgeColor','red');
plot(A_Rand_wtr_iso*Iscalar,B_Rand_wtr_iso*Gscalar,'o','MarkerSize',2,'Parent',subplot1,'MarkerFaceColor','blue','LineWidth',1,'MarkerEdgeColor','blue');
plot(A_Pred_wtr_iso*Iscalar,B_Pred_wtr_iso*Gscalar,'o','MarkerSize',10,'Parent',subplot1,'MarkerFaceColor','blue','LineWidth',4,'MarkerEdgeColor','black');
plot(A_Pred_gwc_iso*Iscalar,B_Pred_gwc_iso*Gscalar,'o','MarkerSize',10,'Parent',subplot1,'MarkerFaceColor','magenta','LineWidth',4,'MarkerEdgeColor','black');
plot(A_Pred_gas_iso*Iscalar,B_Pred_gas_iso*Gscalar,'o','MarkerSize',10,'Parent',subplot1,'MarkerFaceColor','red','LineWidth',4,'MarkerEdgeColor','black');

% Plot Background trend
plot(BGmin(:,1),BGmin(:,2),'--','Parent',subplot1,'LineWidth',3,'color','green');
plot(BGmax(:,1),BGmax(:,2),'--','Parent',subplot1,'LineWidth',3,'color','magenta');
plot(BGminl(:,1),BGminl(:,2),'--','Parent',subplot1,'LineWidth',3,'color','black');

% Plot Axes
plot(XAxis(:,1),XAxis(:,2),'Parent',subplot1,'LineWidth',2,'color','black');
plot(YAxis(:,1),YAxis(:,2),'Parent',subplot1,'LineWidth',2,'color','black');

% Create title
title('Intercept v Gradient');

% Create subplot
subplot2 = subplot(1,2,2,'Parent',figure6,'YTick',[MinR:0.1:MaxR]);
box(subplot2,'on');
hold(subplot2,'all');
xlim(subplot2,[0 MaxAngle]);
ylim(subplot2,[MinR,MaxR]);

% Create plot
plot(0:MaxAngle,RPPpred_gas,'color','red','LineWidth',2,'Parent',subplot2);
plot(0:MaxAngle,RPPpred_wtr,'color','blue','LineWidth',2,'Parent',subplot2);
plot(0:MaxAngle,RPPpred_shlslt,'color','green','LineWidth',2,'Parent',subplot2);
plot(0:MaxAngle,RPPpred_sltshl,'color','cyan','LineWidth',2,'Parent',subplot2);
plot(0:MaxAngle,RPPpred_gwc,'color','magenta','LineWidth',2,'Parent',subplot2);
plot(0:MaxAngle,FS_Gas,'color','red','LineWidth',2,'LineStyle',':','Parent',subplot2);
plot(0:MaxAngle,FS_Wtr,'color','blue','LineWidth',2,'LineStyle',':','Parent',subplot2);
plot(0:MaxAngle,RPPshuey_Pred_gas_aniso,'color','red','LineWidth',2,'LineStyle','--','Parent',subplot2);
plot(0:MaxAngle,RPPshuey_Pred_wtr_aniso,'color','blue','LineWidth',2,'LineStyle','--','Parent',subplot2);
plot(0:MaxAngle,RPPshuey_Pred_shlslt_aniso,'color','green','LineWidth',2,'LineStyle','--','Parent',subplot2);
plot(0:MaxAngle,RPPshuey_Pred_sltshl_aniso,'color','cyan','LineWidth',2,'LineStyle','--','Parent',subplot2);

% Create title
title('Reflectivity vs Angle');

%% Figure 7 Rock physics crossplots
figure7 = figure;
set(figure7, 'Position', [Width LengthSh+(0.1*LengthSh) Width LengthSh])

% Create subplot
subplot1 = subplot(1,2,1,'Parent',figure7,'XTick',[0:0.05:0.4],'YTick',[MinVp:500:MaxVp]);
box(subplot1,'on');
hold(subplot1,'all');
xlim(subplot1,[0 0.4]);
ylim(subplot1,[MinVp MaxVp]);

% Create plot
scatter(Por_filt_wtr,VpWtr_filt,100,Pves_filt_sand_wtr,'.','Parent',subplot1);
plot(Porosity_predicted,VpWtr_pred,'o','MarkerSize',10,'Parent',subplot1,'MarkerFaceColor','black','LineWidth',4,'MarkerEdgeColor','magenta');
plot(Phic,Vpsc,'Parent',subplot1,'color','red');
plot(Phis,Vpss,'Parent',subplot1,'color','red');
plot(Phi,VpSat,'Parent',subplot1,'color','blue');

% Create title
title('Porosity v Vp (m/s)');

%plot(A_Shl_Wtr_iso,B_Shl_Wtr_iso,'.','MarkerSize',30,'Parent',subplot1,'color','blue');
%plot(A_Pred_gas_iso,B_Pred_gas_iso,'o','MarkerSize',10,'Parent',subplot1,'MarkerFaceColor','red','LineWidth',3,'MarkerEdgeColor','black');
%plot(A_Pred_wtr_iso,B_Pred_wtr_iso,'o','MarkerSize',10,'Parent',subplot1,'MarkerFaceColor','blue','LineWidth',3,'MarkerEdgeColor','black');
  
% Create subplot
subplot2 = subplot(1,2,2,'Parent',figure7,'XTick',[4000:1000:12000],'YTick',[1:0.5:4.5]);
box(subplot2,'on');
hold(subplot2,'all');
xlim(subplot2,[4000 12000]);
ylim(subplot2,[1 4.5]);

% Create plot
scatter(AI_Wtr,VpVsR_Wtr,100,Pves_filt_sand_wtr,'.','Parent',subplot2);
plot(AI_Wtr_pred,VpVsR_Wtr_pred,'o','MarkerSize',10,'Parent',subplot2,'MarkerFaceColor','black','LineWidth',4,'MarkerEdgeColor','magenta');
plot(AIss,VpVs_R_ss,'Parent',subplot2,'color','red');
plot(AISat,VpVs_R_Sat,'Parent',subplot2,'color','blue');


% Create title
title('AI v Vp/Vs Ratio');


%% Figure 8 EER projection histograms
figure8 = figure;
set(figure8, 'Position', [Width LengthSh+(0.1*LengthSh) Width LengthSh])

% Create subplot
subplot1 = subplot(2,1,1,'Parent',figure8);
box(subplot1,'on');
hold(subplot1,'all');

hist(MinEERgas_rnd,30);
hold on
hist(MinEERwtr_rnd,30);
hold on
hist(MinEERbg_rnd,30);
h = findobj(gca,'Type','patch');
set(h(1),'FaceColor','green','EdgeColor','green','facealpha',0.5,'Parent',subplot1)
set(h(2),'FaceColor','blue','EdgeColor','blue','facealpha',0.5,'Parent',subplot1)
set(h(3),'FaceColor','red','EdgeColor','red','facealpha',1,'Parent',subplot1)

% Create title
title('Extended Elastic Reflectivity: Minimum Energy Projection');

% Create subplot
subplot2 = subplot(2,1,2,'Parent',figure8);
box(subplot2,'on');
hold(subplot2,'all');

hist(MinlEERgas_rnd,30);
hold on
hist(MinlEERwtr_rnd,30);
hold on
hist(MinlEERbg_rnd,30);
h = findobj(gca,'Type','patch');
set(h(1),'FaceColor','green','EdgeColor','green','facealpha',0.5,'Parent',subplot2)
set(h(2),'FaceColor','blue','EdgeColor','blue','facealpha',0.5,'Parent',subplot2)
set(h(3),'FaceColor','red','EdgeColor','red','facealpha',1,'Parent',subplot2)

% Create title
title('Extended Elastic Reflectivity: Fluid Projection');


%% Figure 9 AI v Vp/Vs
figure9 = figure;
set(figure9, 'Position', [0 LengthSh+(0.1*LengthSh) Width LengthSh])

% Create plot
scatterhist(AI_r,VpVsR_r,'Group',Lith,'LineStyle',{'-','-','-','-'},'LineWidth',[2,2,2,2]);

% Create title
title('AI v Vp/Vs Ratio');

%% Figure 10 I-G Crossplot for gas and brine facies coloured by Porosity
figure10 = figure;
set(figure10, 'Position', [Width LengthSh+(0.1*LengthSh) Width LengthSh])

% Create subplot
subplot1 = subplot(1,2,1,'Parent',figure10,'XTick',[MinI:Step:MaxI],'YTick',[MinG:Step:MaxG]);
box(subplot1,'on');
hold(subplot1,'all');
xlim(subplot1,[MinI MaxI]);
ylim(subplot1,[MinG MaxG]);

% Plot Axes
plot(XAxis(:,1),XAxis(:,2),'Parent',subplot1,'LineWidth',2,'color','black');
plot(YAxis(:,1),YAxis(:,2),'Parent',subplot1,'LineWidth',2,'color','black');

% Plot Data
scatter(A_Rand_gas_iso*Iscalar,B_Rand_gas_iso*Gscalar,100,Por_r_gas,'.','Parent',subplot1)

% Create subplot
subplot2 = subplot(1,2,2,'Parent',figure10,'XTick',[MinI:Step:MaxI],'YTick',[MinG:Step:MaxG]);
box(subplot2,'on');
hold(subplot2,'all');
xlim(subplot2,[MinI MaxI]);
ylim(subplot2,[MinG,MaxG]);

% Plot Axes
plot(XAxis(:,1),XAxis(:,2),'Parent',subplot2,'LineWidth',2,'color','black');
plot(YAxis(:,1),YAxis(:,2),'Parent',subplot2,'LineWidth',2,'color','black');

% Plot Data
scatter(A_Rand_wtr_iso*Iscalar,B_Rand_wtr_iso*Gscalar,100,Por_r_wtr,'.','Parent',subplot2)

end
function downlogs = downsamplelog(logs,sampint)
    downlogs = logs(1:sampint:end,:);
end


