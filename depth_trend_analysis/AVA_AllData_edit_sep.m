%% Rock Physics Database and Depth Trend Modelling Software
% Written by Matt Bolton April 2015
% Designed to incorporate all the rock properties of the sands, silts and mudrocks penetrated in the Tanzania offshore margin.
% Depth trends are then applied to the data after filtering for specific sand and caprock facies e.g. Age, Volume of Shale, Diagenesis, etc
% Seismic interfaces for can then be modelled at any depth based on trends fitted. 
% Reflectivity with angle is estimated for the isotropic case using the full Zoeppritz solution.  An anisotropic case is also included based on the Shuey approximation.

%% DATA READ
function AVA_AllData_edit()
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

well_no_datapoints(:,2) = ceil(well_no_datapoints(:,1)./min(well_no_datapoints(:,1)));
    
for i_file = 1:1:nfiles
    % Data decimation parameter
    %well_data{i_file}.curves_dec = downsample(well_data{i_file}.curves_wonans,well_no_datapoints(i_file,2));
    well_data{i_file}.curves_dec = downsample(well_data{i_file}.curves_wonans,well_no_datapoints(i_file,2));
     
    for i_curve = 1:size(well_data{i_file}.curves_dec,2);
        well_data{i_file}.(well_data{i_file}.curve_info{i_curve})=well_data{i_file}.curves_dec(:,i_curve);
    end
    
end

curve_info = well_data{1}.curve_info;
for i_curve = 1:size(well_data{i_file}.curves_dec,2);
   start_row = 1;
   for i_file = 1:1:nfiles
        fprintf('Concatenating %s curves\n',well_data{1}.curve_info{i_curve})
        end_row = size(well_data{i_file}.(well_data{1}.curve_info{i_curve}),1)+start_row-1;
        curve_info{i_curve,4}(start_row:end_row,1) = well_data{i_file}.(well_data{1}.curve_info{i_curve});
        start_row = end_row + 1;
   end   
   eval(strcat(curve_info{i_curve,1},'= curve_info{',num2str(i_curve),',4};'))
end
  

%Read in log files containing petrophysical interpretation and Vp, Vs and Rho measurements
% [WellA,las_hdr]=read_las2_file('Well_A_DepthTrends.las',0);
% [WellB,las_hdr]=read_las2_file('Well_B_DepthTrends.las',0);
% [WellC,las_hdr]=read_las2_file('Well_C_DepthTrends.las',0);
% [WellI,las_hdr]=read_las2_file(logiA=~logical(sum(isnan(WellA.curves),2));
% [WellJ,las_hdr]=read_las2_file('Well_J_DepthTrends.las',0);
% [WellK,las_hdr]=read_las2_file('Well_K_DepthTrends.las',0);
% [WellX,las_hdr]=read_las2_file('Well_X_DepthTrends.las',0);


% % Data decimation parameter
% N=10;
% WellA.curves=downsample(WellA.curves,N);
% WellB.curves=downsample(WellB.curves,N);
% WellC.curves=downsample(WellC.curves,N);
% WellI.curves=downsample(WellI.curves,N);
% WellJ.curves=downsample(WellJ.clogiA=~logical(sum(isnan(WellA.curves),2));
% WellK.curves=downsample(WellK.curves,N);
% WellX.curves=downsample(WellX.curves,N);

% Determine log names
% lognamesA=WellA.curve_info(:,1);
% lognamesB=WellB.curve_info(:,1);
% lognamesC=WellC.curve_info(:,1);
% lognamesI=WellI.curve_info(:,1);
% lognamesJ=WellJ.curve_info(:,1);
% lognamesK=WellK.curve_info(:,1);
% lognamesX=WellX.curve_info(:,1);

% Find samples where data exists for all logs
% logiA=~logical(sum(isnan(WellA.curves),2));
% lengthA=sum(logiA,1);
% logiB=~logical(sum(isnan(WellB.curves),2));
% lengthB=sum(logiB,1);
% logiC=~logical(sum(isnan(WellC.curves),2));
% lengthC=sum(logiC,1);
% logiI=~logical(sum(isnan(WellI.curves),2));
% lengthI=sum(logiI,1);
% logiJ=~logical(sum(isnan(WellJ.curves),2));
% lengthJ=sum(logiJ,1);
% logiK=~logical(sum(isnan(WellK.curves),2));
% lengthK=sum(logiK,1);
% logiX=~logical(sum(isnan(WellX.curves),2));
% lengthX=sum(logiX,1);

% logi=vertcat(logiA,logiB,logiC,logiI,logiJ,logiK,logiX);
% 
% for ii=1:size(WellA.curves,2);
% WellA.(lognamesA{ii})=WellA.curves(:,ii);
% WellB.(lognamesB{ii})=WellB.curves(:,ii);
% WellC.(lognamesC{ii})=WellC.curves(:,ii);
% WellI.(lognamesI{ii})=WellI.curves(:,ii);
% WellJ.(lognamesJ{ii})=WellJ.curves(:,ii);
% WellK.(lognamesK{ii})=WellK.curves(:,ii);
% WellX.(lognamesX{ii})=WellX.curves(:,ii);
% end

% Vp=vertcat(WellA.Vp,WellB.Vp,WellC.Vp,WellI.Vp,WellJ.Vp,WellK.Vp,WellX.Vp);
% Vs=vertcat(WellA.Vs_Fast,WellB.Vs_Fast,WellC.Vs_Fast,WellI.Vs_Fast,WellJ.Vs_Fast,WellK.Vs_Fast,WellX.Vs_Fast);
% Rho=vertcat(WellA.DENS,WellB.DENS,WellC.DENS,WellI.DENS,WellJ.DENS,WellK.DENS,WellX.DENS);
% TVDml=vertcat(WellA.Depth_TVDml,WellB.Depth_TVDml,WellC.Depth_TVDml,WellI.Depth_TVDml,WellJ.Depth_TVDml,WellK.Depth_TVDml,WellX.Depth_TVDml);
% P=vertcat(WellA.MDT_Pform,WellB.MDT_Pform,WellC.MDT_Pform,WellI.MDT_Pform,WellJ.MDT_Pform,WellK.MDT_Pform,WellX.MDT_Pform);
% PhiE=vertcat(WellA.PHIE,WellB.PHIE,WellC.PHIE,WellI.PHIE,WellJ.PHIE,WellK.PHIE,WellX.PHIE);
% PhiT=vertcat(WellA.PHIT,WellB.PHIT,WellC.PHIT,WellI.PHIT,WellJ.PHIT,WellK.PHIT,WellX.PHIT);
% Sw=vertcat(WellA.SW,WellB.SW,WellC.SW,WellI.SW,WellJ.SW,WellK.SW,WellX.SW);
% Vsh=vertcat(WellA.VSH,WellB.VSH,WellC.VSH,WellI.VSH,WellJ.VSH,WellK.VSH,WellX.VSH);
% Age=vertcat(WellA.Age,WellB.Age,WellC.Age,WellI.Age,WellJ.Age,WellK.Age,WellX.Age);
% GR=vertcat(WellA.GR,WellB.GR,WellC.GR,WellI.GR,WellJ.GR,WellK.GR,WellX.GR);
% NEU=vertcat(WellA.NEUT,WellB.NEUT,WellC.NEUT,WellI.NEUT,WellJ.NEUT,WellK.NEUT,WellX.NEUT);

% Vs_Fast=Vs_Fast(logi);
% DENS=DENS(logi);
% Depth_TVDml=Depth_TVDml(logi)
% MDT_Pform=MDT_Pform(logi)./0.006895;
% PHIE=PHIE(logi);
% PHIT=PHIT(logi);
% SW=SW(logi);
% VSH=VSH(logi);
% % Age=Age(logi);
% % GR=GR(logi);
% NEUT=NEUT(logi);

Zml=round(Depth_TVDml)

%% PROSPECT INPUTS
% Geologically constrained depth trend analysis
% Key criteria for selecting the most appropriate sands and caprock for the prospect being modelled

%Zpred=input('Enter prospect TVDml(m)'); % Input prospect depth
%Zwc=input('Enter prospect water depth (m)'); % Water depth for pressure
%Pop=input('Overpressure expected in sand (psi); 0 = normal'); % Overpressure
%Sand=input('Enter Reservoir shale content; 0 = Clean Sand, 1 = Shale'); % Reservoir facies
%Diag=input('Enter expected Sand diagenesis; 0=all, 1 = unconsolidated, 2 = consolidated, 3 = cemented'); % Diagenesis
%GA=input('Enter expected reservoir age (Ma)'); % Geological age
%Shale=input('Enter expected overburden; 0=Shl or 1=Slt'); % Overburden type
%FitType=input('Enter depth trend fitting method; 1=Data driven,2=Constrained, 3=Rock Physics Model;
%NI=input('Enter number of monte carlo simulations');

% Test Inputs
Zpred=4300;
Zwc=1700;
Pop=1000;
Sand=0.05;
Comp=0;
GA=130;
Shale=1;
FitType=2;
NI=2500;

%% BASIC INPUTS
%Zsb=1153;
Pco=1000; % Only include data in depth trend analysis which is at close to hydrostatic pressure (i.e. < Pco psi overpressued)

% Facies Classification
SndPorCO=0.08; % Net sand cut-off
SndVshCO=0.7; % Net sand cut-off
SndSwCO=0.5; % Net pay sand cut-off
SltVshCO=0.9; % Silt-Shale cut-off

% Create diagensis log based on petrophyiscal cut offs and any other
% influencing factor.  Here Age is used, but this could also be T or P
% based.
for aa = 1:length(Vp); 
    if Comp == 0 & VSH(aa) < SndVshCO & PHIE(aa) >= SndPorCO
    Diag(aa,1)=0; % All Net sand points
        else if Comp>0 & VSH(aa) < SndVshCO & PHIE(aa) >= SndPorCO & Age(aa) <= 66
        Diag(aa,1)=1; % Unconsolidated Tertiary sands
            else if Comp>0 & VSH(aa) < SndVshCO & PHIE(aa) >= SndPorCO & Age(aa) > 66
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

% Create facies logs based on cut-offs
% Sandstone facies
VpGas=Vp(Diag==Comp & SW <= SndSwCO);
VpWtr=Vp(Diag==Comp & SW > SndSwCO);
VsGas=Vs_Fast(Diag==Comp & SW <= SndSwCO);
VsWtr=Vs_Fast(Diag==Comp & SW > SndSwCO);
RhoGas=DENS(Diag==Comp & SW <= SndSwCO);
RhoWtr=DENS(Diag==Comp & SW > SndSwCO);

% Overburden facies
VpSlt=Vp(VSH >= SndVshCO & VSH < SltVshCO & PHIE < SndPorCO);
VpShl=Vp(VSH >= SltVshCO & PHIE < SndPorCO);
VsSlt=Vs_Fast(VSH >= SndVshCO & VSH < SltVshCO & PHIE < SndPorCO);
VsShl=Vs_Fast(VSH >= SltVshCO & PHIE < SndPorCO);
RhoSlt=DENS(VSH >= SndVshCO & VSH < SltVshCO & PHIE < SndPorCO);
RhoShl=DENS(VSH >= SltVshCO & PHIE < SndPorCO);
Pob=MDT_Pform(VSH >= SltVshCO & PHIE < SndPorCO);

% Max incidence angle to calculate AVA over
MaxAngle=45; %degrees

% Depth or pressure iteration
MaxZ=10000; % max depth (or pressure (psi)) to calculate to
Z=0:1:MaxZ-1; % Z index with 1m sampling (requires all input depths are at this precision)

% Rock Physics Modelling
Clay=Sand; % fraction; Volume of clay 
Feldspar=0; % fraction; Volume of feldpar
Calcite=0; % fraction; Volume of calcite
Por_min=0.0; % sandstone minimum porosity at infinite pressure
Por_min_ob=0.1; % mudrock minimum porosity at infinite pressure
%phiET=(-0.87*Sand)+0.87; % fraction of shale in pore space; 0 = shale matrix, 1 = all shale in pore space
phiET=0.25; % Fraction of Vsh which fills pore space
fcement=0.03; % cement fraction
PhiC=0.45; % Clean Sand critical porosity
PhiCsh=0.55; % Clean Shale critical porosity
PhiC=(PhiC-(phiET*Clay)); % critical porosity can be reduced through shale in pore space and cementation
shearfact=0.6; % Shear wave reduction factor
Coord=20; % Co-ordination number
RhoMin=2; % Minimum shale density

% Model to Data scaling (This needs to be improved using a method such as covariance matrix matching between the modelled and the seismic data.
Iscalar=1520;
Gscalar=1800;

% Rock property limits used to constrain the maximum possible properties (used in depth trend calculations)
RhoQ=2.65; % Sand limits
VpQ=6038;
VsQ=4121;
RhoI=2.8; % Illite rich shale limits
VpI=5264;
VsI=2719;
RhoS=2.5; % Smectite rich shale limits
VpS=3282;
VsS=2254;
RhoOBmax=(RhoI+RhoS)/2; % Average shale limit

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
HCcolumn=200; % m; this is required to determine an average gas buoyancy effect

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

% Anisotropy
delta1=0.12;
delta2=0;
epsilon1=0.2;
epsilon2=0;

%% TEMPERATURE & PRESSURE CALCULATIONS
% All depth trends should be calculated wrt to vertical effective stress if
% basin is not hydrostatically pressured.

% EDITABLE CONTENT SEE LINE 56

% CALCULATIONS
% Batzle and wang calculations
[Kreuss,rhoeff,Kvoigt,vpb,rhob,Kb,vpo,rhoo,Ko,vpg,rhog,Kg,gor]=flprop(method,sal,og,gg,gor,giib,giio,Pp,Tp,So,Sg);

RhoWtr=mean(rhob) % Used as average overburden brine properties for prediction
RhoW=1;

% Average overburden profile
RhoOB=RhoShl; % Overburden density data
VpOB=VpShl; % Overburden Vp data
VsOB=VsShl; % Overburden Vs data
ZOB=Zml(VSH >= SltVshCO & PHIE < SndPorCO);
PhiTob=NEUT(VSH >= SltVshCO & PHIE < SndPorCO);
ini(MaxZ,1)=0;

% Exponential curve fit for shale trend
dRhoOB=ZOB;
GRhoOB=log((RhoOB-RhoMin)/(RhoMin-RhoI));
mRhoOB=GRhoOB\dRhoOB;
betaRhoOB=real(-1/mRhoOB);

GrOB=exp(-betaRhoOB*ZOB);
d=RhoOB-RhoI; % known data
mx=RhoI; % max shale density
mRhoMudOB=ModelFitv2(d,GrOB,betaRhoOB,mx);
RhoL = mx -((mx-mRhoMudOB)*exp(-betaRhoOB*Z));

%dPhiOB=ZOB;
%GPhiOB=log((PhiTob-Por_min_ob)/(Por_min_ob-PhiCsh));
%mPhiOB=GPhiOB\dPhiOB;
%betaPhiOB=real(-1/mPhiOB);


%GrOB=exp(-betaPhiOB*ZOB);
%d=RhoOB-RhoI; % known data
%mx=RhoI; % max shale density
%mRhoMudOB=ModelFitv2(d,GrOB,betaPhiOB,mx);
%RhoL = RhoI-((RhoI-mRhoMudOB)*exp(-betaPhiOB*Z));

% Integration of density trend (RhoL) to define average lithostatic trend
for z = 1:MaxZ-2
    Plith(z+1)=1000*(Z(z+2)-Z(z+1))*(RhoL(z+1));
    ini(z+1)=ini(z)+Plith(z+1);
end

ini(MaxZ)=ini(MaxZ-1); % account for array indexing

%Convert from Pa to PSI: PSI = 0.000145Pa and add water column.
for bb=1:nfiles;
    for aa=1:length(Z);
        %Create Hydrostatic curve
        Phydro(aa,bb)=1000*1E-6*9.81*((RhoW*Z(aa))+(RhoSea*Zsb(bb)));
        %Create Lithostatic curve
        Plith_f(aa,bb)=(1E-6*9.81)*(ini(aa,1)+(Zsb(bb)*RhoSea*1000));
    end
end

% Calculate Vertical Effective Stress for each point in the database using lithostatic curve (Plith) and Formation Pressure (P)
tic
for w = 1:nfiles;
    for y = 1:length(well_data{w}.curves_dec);
        for v = 1:MaxZ       
            if Z(v) == round(well_data{w}.Depth_TVDml(y))
            Pv = Plith_f(v,w)-well_data{w}.MDT_Pform(y);
            Pves(y,w)=Pv;
            Pover(y,w)=well_data{w}.MDT_Pform(y)-Phydro(v,w);
            end
        end
    end
end
toc
%for i_con = 1:length(Pves)
%   start_row = 1;
%   for i_file = 1:1:nfiles
%        end_row = well_no_datapoints{i_file,1}/well_no_datapoints{i_file,2}+start_row-1;
%        Pover(i_con) = Pover(i_con,I_file).(well_data{1}.curve_info{i_curve});
%        start_row = end_row + 1;
%   end   
%   eval(strcat(curve_info{i_curve,1},'= curve_info{',num2str(i_curve),',4};'))
%end

P1=Pover(1:lengthA,1);
P2=Pover(lengthA+1:lengthA+lengthB,2);
P3=Pover(lengthA+lengthB+1:lengthA+lengthB+lengthC,3);
P4=Pover(lengthA+lengthB+lengthC+1:lengthA+lengthB+lengthC+lengthI,4);
P5=Pover(lengthA+lengthB+lengthC+lengthI+1:lengthA+lengthB+lengthC+lengthI+lengthJ,5);
P6=Pover(lengthA+lengthB+lengthC+lengthI+lengthJ+1:lengthA+lengthB+lengthC+lengthI+lengthJ+lengthK,6);
P7=Pover(lengthA+lengthB+lengthC+lengthI+lengthJ+lengthK+1:lengthA+lengthB+lengthC+lengthI+lengthJ+lengthK+lengthX,7);

Pover=vertcat(P1,P2,P3,P4,P5,P6,P7);

P1=Pves(1:lengthA,1);
P2=Pves(lengthA+1:lengthA+lengthB,2);
P3=Pves(lengthA+lengthB+1:lengthA+lengthB+lengthC,3);
P4=Pves(lengthA+lengthB+lengthC+1:lengthA+lengthB+lengthC+lengthI,4);
P5=Pves(lengthA+lengthB+lengthC+lengthI+1:lengthA+lengthB+lengthC+lengthI+lengthJ,5);
P6=Pves(lengthA+lengthB+lengthC+lengthI+lengthJ+1:lengthA+lengthB+lengthC+lengthI+lengthJ+lengthK,6);
P7=Pves(lengthA+lengthB+lengthC+lengthI+lengthJ+lengthK+1:lengthA+lengthB+lengthC+lengthI+lengthJ+lengthK+lengthX,7);

Pves=vertcat(P1,P2,P3,P4,P5,P6,P7);

% Calculate pressures based on input parameters (Prospect modelling)
% Lithostatic
Plith_pred=(0.000145*9.81)*(ini(:,1)+(Zwc*RhoSea*1000));
% Hydrostatic
Phy_pred=1000*0.000145*RhoW*9.81*((Zpred*RhoW)+(Zwc*RhoSea));

% Prospect Pore pressure prediction
P_pred_sand=Phy_pred+Pop; %psi
Pp_pred_sand=0.006895*P_pred_sand; % MPa

% Prospect fluid properties predicted based on inputs
[Kreuss,rhoeff,Kvoigt,vpb,RhoW_pred,Kb,vpo,rhoo,Ko,vpg,RhoG_pred,Kg,gor]=flprop(method,sal,og,gg,gor,giib,giio,Pp_pred_sand,T_pred,So,Sg);

% Buoyancy effect in gas leg: Pform_gas > Pform_wtr at same depth
Buoy=(RhoW_pred-RhoG_pred)*9.81*HCcolumn*0.145; % Gas buoyancy effect

% Prospect vertical effective stress
Pves_pred=Plith_pred(Zpred)-Phy_pred; %psi
PVES=Pves_pred*0.006895; %MPa

%% DATA FILTERING
% Using prospect inputs the database is filtered to include only data which
% is appropriate for the modelling scenario.

% Volume of shale filter included +/- 5% of input value
if Sand-0.05>0 % Limit low case to 0
    Sand_low=Sand-0.05;
else
    Sand_low=0;
end
Sand_hi=Sand+0.05;

% Geological Age filter
COage=65; % Cut off Age which overburden rock properties show different compaction trends
if GA>COage 
    GA_low=COage;
    GA_hi=200; % Base Jurassic
else
    GA_low=0;
    GA_hi=COage;
end

% Filter sand points based on saturation, volume of shale and porosity
VpGas_filt=Vp(VSH <= Sand_hi & VSH >= Sand_low & Diag==Comp & SW <= SndSwCO & Pover < Pco);
VpWtr_filt=Vp(VSH <= Sand_hi & VSH >= Sand_low & Diag==Comp & SW > SndSwCO & Pover < Pco);
VsGas_filt=Vs_Fast(VSH <= Sand_hi & VSH >= Sand_low & Diag==Comp & SW <= SndSwCO & Pover < Pco);
VsWtr_filt=Vs_Fast(VSH <= Sand_hi & VSH >= Sand_low & Diag==Comp & SW > SndSwCO & Pover < Pco);
RhoGas_filt=DENS(VSH <= Sand_hi & VSH >= Sand_low & Diag==Comp & SW <= SndSwCO & Pover < Pco);
RhoWtr_filt=DENS(VSH <= Sand_hi & VSH >= Sand_low & Diag==Comp & SW > SndSwCO & Pover < Pco);
Por_filt=PHIT(VSH <= Sand_hi & VSH >= Sand_low & Diag==Comp & Pover < Pco);
Por_filt_wtr=PHIT(VSH <= Sand_hi & VSH >= Sand_low & Diag==Comp & SW > SndSwCO & Pover < Pco);
Por_filt_gas=PHIT(VSH <= Sand_hi & VSH >= Sand_low & Diag==Comp & SW <= SndSwCO & Pover < Pco);
Age_filt_gas=Age(VSH <= Sand_hi & VSH >= Sand_low & Diag==Comp & SW <= SndSwCO & Pover < Pco);
Age_filt_wtr=Age(VSH <= Sand_hi & VSH >= Sand_low & Diag==Comp & SW > SndSwCO & Pover < Pco);
%Facies_filt_gas=Facies(Facies==Diag & vsh <= Sand_hi & vsh >= Sand_low);
vsh_filt=VSH(VSH <= Sand_hi & VSH >= Sand_low & Diag==Comp & Pover < Pco);

% Filter vertical effective stress points to get corresponding data for
% curve fitting to sand
Pves_filt_sand_gas=Pves(VSH <= Sand_hi & VSH >= Sand_low & Diag==Comp & SW <= SndSwCO & Pover < Pco);
Pves_filt_sand_wtr=Pves(VSH <= Sand_hi & VSH >= Sand_low & Diag==Comp & SW > SndSwCO & Pover < Pco);
Pves_filt_sand=Pves(VSH <= Sand_hi & VSH >= Sand_low & Diag==Comp & Pover < Pco);

% Filter Overburden facies
Pves_filt_slt=Pves(VSH >= SndVshCO & VSH < SltVshCO & PHIE < SndPorCO & Age>GA_low & Age<GA_hi & Pover < Pco);
Pves_filt_shl=Pves(VSH >= SltVshCO & PHIE < SndPorCO & Age>GA_low & Age<GA_hi & Pover < Pco);
PhiT_filt_slt=NEUT(VSH >= SndVshCO & VSH < SltVshCO & PHIE < SndPorCO & Age>GA_low & Age<GA_hi & Pover < Pco);
PhiT_filt_shl=NEUT(VSH >= SltVshCO & PHIE < SndPorCO & Age>GA_low & Age<GA_hi & Pover < Pco);
RhoShl_filt=DENS(VSH >= SltVshCO & PHIE < SndPorCO & Age>GA_low & Age<GA_hi & Pover < Pco);
VpShl_filt=Vp(VSH >= SltVshCO & PHIE < SndPorCO & Age>GA_low & Age<GA_hi & Pover < Pco);
VsShl_filt=Vs_Fast(VSH >= SltVshCO & PHIE < SndPorCO & Age>GA_low & Age<GA_hi & Pover < Pco);
RhoSlt_filt=DENS(VSH >= SndVshCO & VSH < SltVshCO & PHIE < SndPorCO & Age>GA_low & Age<GA_hi & Pover < Pco);
VpSlt_filt=Vp(VSH >= SndVshCO & VSH < SltVshCO & PHIE < SndPorCO & Age>GA_low & Age<GA_hi & Pover < Pco);
VsSlt_filt=Vs_Fast(VSH >= SndVshCO & VSH < SltVshCO & PHIE < SndPorCO & Age>GA_low & Age<GA_hi & Pover < Pco);
Age_filt_slt=Age(VSH >= SndVshCO & VSH < SltVshCO & PHIE < SndPorCO & Age>GA_low & Age<GA_hi & Pover < Pco);
Age_filt_shl=Age(VSH >= SltVshCO & PHIE < SndPorCO & Age>GA_low & Age<GA_hi & Pover < Pco);
GR_filt_slt=GR(VSH >= SndVshCO & VSH < SltVshCO & PHIE < SndPorCO & Age>GA_low & Age<GA_hi & Pover < Pco);
GR_filt_shl=GR(VSH >= SltVshCO & PHIE < SndPorCO & Age>GA_low & Age<GA_hi & Pover < Pco);

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
            RhoMud=RhoS;
            VpMud=VpS;
            VsMud=VsS;
        end
   
else
    RhoOB_filt=DENS(VSH >= SndVshCO & VSH < SltVshCO & PHIE < SndPorCO & Age>GA_low & Age<GA_hi & Pover < Pco);
    VpOB_filt=Vp(VSH >= SndVshCO & VSH < SltVshCO & PHIE < SndPorCO & Age>GA_low & Age<GA_hi & Pover < Pco);
    VsOB_filt=Vs_Fast(VSH >= SndVshCO & VSH < SltVshCO & PHIE < SndPorCO & Age>GA_low & Age<GA_hi & Pover < Pco);
    Pves_filt_OB=Pves_filt_slt;
    PhiT_filt_OB=PhiT_filt_slt;
    % Voight average for silt end member rock properties using Slt and Shl cut-off
    if GA >= COage
        RhoMud=(((SndVshCO+SltVshCO)/2)*RhoI)+((1-((SndVshCO+SltVshCO)/2))*RhoQ);
        VpMud=(((SndVshCO+SltVshCO)/2)*VpI)+((1-((SndVshCO+SltVshCO)/2))*VpQ);
        VsMud=(((SndVshCO+SltVshCO)/2)*VsI)+((1-((SndVshCO+SltVshCO)/2))*VsQ);
    else
        RhoMud=(((SndVshCO+SltVshCO)/2)*RhoS)+((1-((SndVshCO+SltVshCO)/2))*RhoQ);
        VpMud=(((SndVshCO+SltVshCO)/2)*VpS)+((1-((SndVshCO+SltVshCO)/2))*VpQ);
        VsMud=(((SndVshCO+SltVshCO)/2)*VsS)+((1-((SndVshCO+SltVshCO)/2))*VsQ);
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

% Similar approach for caprock however as PHIT not always present, RHOB
% used as a proxy
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
d=RhoGas_filt-RhoQ; % known data
mx=RhoQ; % Quartz density
mRhoGas=ModelFitv2(d,Ggas,betaPor,mx);
RhoGas_pred = RhoQ-((RhoQ-mRhoGas)*exp(-betaPor*Pves_pred));
RhoGas_trend = RhoQ-((RhoQ-mRhoGas)*exp(-betaPor*Z));

% Vp Gas
d=VpGas_filt-VpQ; % known data
mx=VpQ; % Quartz Vp
mVpGas=ModelFitv2(d,Ggas,betaPor,mx);
VpGas_pred = VpQ-((VpQ-mVpGas)*exp(-betaPor*Pves_pred));
VpGas_trend = VpQ-((VpQ-mVpGas)*exp(-betaPor*Z));

% Vs Gas
d=VsGas_filt-VsQ; % known data
mx=VsQ; % Quartz Vs
mVsGas=ModelFitv2(d,Ggas,betaPor,mx);
VsGas_pred = VsQ-((VsQ-mVsGas)*exp(-betaPor*Pves_pred));
VsGas_trend = VsQ-((VsQ-mVsGas)*exp(-betaPor*Z));

% Water Sand Rock Properties

% Rho Water
d=RhoWtr_filt-RhoQ; % known data
mx=RhoQ; % Quartz density
mRhoWtr=ModelFitv2(d,Gwtr,betaPor,mx);
RhoWtr_pred = RhoQ-((RhoQ-mRhoWtr)*exp(-betaPor*Pves_pred));
RhoWtr_trend = RhoQ-((RhoQ-mRhoWtr)*exp(-betaPor*Z));

% Vp Water
d=VpWtr_filt-VpQ; % known data
mx=VpQ; % Quartz Vp
mVpWtr=ModelFitv2(d,Gwtr,betaPor,mx);
VpWtr_pred = VpQ-((VpQ-mVpWtr)*exp(-betaPor*Pves_pred));
VpWtr_trend = VpQ-((VpQ-mVpWtr)*exp(-betaPor*Z));

% Vs Water
d=VsWtr_filt-VsQ; % known data
mx=VsQ; % Quartz Vs
mVsWtr=ModelFitv2(d,Gwtr,betaPor,mx);
VsWtr_pred = VsQ-((VsQ-mVsWtr)*exp(-betaPor*Pves_pred));
VsWtr_trend = VsQ-((VsQ-mVsWtr)*exp(-betaPor*Z));

% Overpressure Modelling
% POROSITY
[Porosity_predicted]=Overpressure(a_sst,b_sst,c_sst,Pves_pred,Pop,0,Porosity_predicted,PhiC,Por_min)

% GAS SAND
[VpGas_pred]=Overpressure(a_sst,b_sst,c_sst,Pves_pred,Pop,Buoy,VpGas_pred,mVsGas,VpQ);
[VsGas_pred]=Overpressure(a_sst,b_sst,c_sst,Pves_pred,Pop,Buoy,VsGas_pred,mVsGas,VsQ);
[RhoGas_pred]=Overpressure(a_sst,b_sst,c_sst,Pves_pred,Pop,Buoy,RhoGas_pred,mRhoGas,RhoQ);

% BRINE SAND
[VpWtr_pred]=Overpressure(a_sst,b_sst,c_sst,Pves_pred,Pop,0,VpWtr_pred,mVsWtr,VpQ);
[VsWtr_pred]=Overpressure(a_sst,b_sst,c_sst,Pves_pred,Pop,0,VsWtr_pred,mVsWtr,VsQ);
[RhoWtr_pred]=Overpressure(a_sst,b_sst,c_sst,Pves_pred,Pop,0,RhoWtr_pred,mRhoWtr,RhoQ);

% Caprock property fitting procedure:

% Filtered overburden rock properties

d=RhoOB_filt-RhoMud; % known data
mx=RhoMud; % Quartz density
mRhoMud=ModelFitv2(d,Gr,betaPhiTob,mx);
RhoOB_pred = RhoMud-((RhoMud-mRhoMud)*exp(-betaPhiTob*Pves_pred));
RhoOB_trend = RhoMud-((RhoMud-mRhoMud)*exp(-betaPhiTob*Z));

d=VpOB_filt-VpMud; % known data
mx=VpMud; % Quartz density
mVpMud=ModelFitv2(d,Gr,betaPhiTob,mx);
VpOB_pred = VpMud-((VpMud-mVpMud)*exp(-betaPhiTob*Pves_pred));
VpOB_trend = VpMud-((VpMud-mVpMud)*exp(-betaPhiTob*Z));

d=VsOB_filt-VsMud; % known data
mx=VsMud; % Quartz density
mVsMud=ModelFitv2(d,Gr,betaPhiTob,mx);
VsOB_pred = VsMud-((VsMud-mVsMud)*exp(-betaPhiTob*Pves_pred));
VsOB_trend = VsMud-((VsMud-mVsMud)*exp(-betaPhiTob*Z));

% Caprock overpressure modelling
% TZA data at Jodari suggests density measurements not effected by
% overpressure.  This corroborates published research and datasets.
% Therefore for the caprock we assume Rho_pred(Pves)=Rho_pred(Pves-Pover)

[VpOB_pred]=Overpressure(a_mud,b_mud,c_mud,Pves_pred,Pop,0,VpOB_pred,mVpMud,VpMud);
[VsOB_pred]=Overpressure(a_mud,b_mud,c_mud,Pves_pred,Pop,0,VsOB_pred,mVsMud,VsMud);


% Caprock property trends
% Similar approach for caprock 
if GA >= COage
    RhoMud=RhoI;
    VpMud=VpI;
    VsMud=VsI;
else
    RhoMud=RhoS;
    VpMud=VpS;
    VsMud=VsS;
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

% Overpressure modelling
[VpShl_pred]=Overpressure(a_mud,b_mud,c_mud,Pves_pred,Pop,0,VpShl_pred,mVpMud,VpMud);
[VsShl_pred]=Overpressure(a_mud,b_mud,c_mud,Pves_pred,Pop,0,VsShl_pred,mVsMud,VsMud);

if GA >= COage
    RhoMud=(((SndVshCO+SltVshCO)/2)*RhoI)+((1-((SndVshCO+SltVshCO)/2))*RhoQ);
    VpMud=(((SndVshCO+SltVshCO)/2)*VpI)+((1-((SndVshCO+SltVshCO)/2))*VpQ);
    VsMud=(((SndVshCO+SltVshCO)/2)*VsI)+((1-((SndVshCO+SltVshCO)/2))*VsQ);
else
    RhoMud=(((SndVshCO+SltVshCO)/2)*RhoS)+((1-((SndVshCO+SltVshCO)/2))*RhoQ);
    VpMud=(((SndVshCO+SltVshCO)/2)*VpS)+((1-((SndVshCO+SltVshCO)/2))*VpQ);
    VsMud=(((SndVshCO+SltVshCO)/2)*VsS)+((1-((SndVshCO+SltVshCO)/2))*VsQ);
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

% Overpressure modelling
[VpSlt_pred]=Overpressure(a_mud,b_mud,c_mud,Pves_pred,Pop,0,VpSlt_pred,mVpMud,VpMud);
[VsSlt_pred]=Overpressure(a_mud,b_mud,c_mud,Pves_pred,Pop,0,VsSlt_pred,mVsMud,VsMud);

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
[VpOB_r VsOB_r RhoOB_r VpGas_r VsGas_r RhoGas_r Por_r_ob Por_r_gas] = RockPropDist(NI,MaxZ,Pves_filt_OB,Pves_filt_sand_gas,VpOB_pred,VsOB_pred,RhoOB_pred,VpGas_pred,VsGas_pred,RhoGas_pred,VpOB_filt,VsOB_filt,RhoOB_filt,VpGas_filt,VsGas_filt,RhoGas_filt,VpOB_trend,VsOB_trend,RhoOB_trend,VpGas_trend,VsGas_trend,RhoGas_trend,PhiT_filt_OB,OB_por_pred,OB_por_trend,Por_filt_gas,Porosity_predicted,Por_trend);
[VpOB_r VsOB_r RhoOB_r VpWtr_r VsWtr_r RhoWtr_r Por_r_ob Por_r_wtr] = RockPropDist(NI,MaxZ,Pves_filt_OB,Pves_filt_sand_wtr,VpOB_pred,VsOB_pred,RhoOB_pred,VpWtr_pred,VsWtr_pred,RhoWtr_pred,VpOB_filt,VsOB_filt,RhoOB_filt,VpWtr_filt,VsWtr_filt,RhoWtr_filt,VpOB_trend,VsOB_trend,RhoOB_trend,VpWtr_trend,VsWtr_trend,RhoWtr_trend,PhiT_filt_OB,OB_por_pred,OB_por_trend,Por_filt_wtr,Porosity_predicted,Por_trend);
[VpSlt_r VsSlt_r RhoSlt_r VpShl_r VsShl_r RhoShl_r Por_r_ob1 Por_r_ob2] = RockPropDist(NI,MaxZ,Pves_filt_slt,Pves_filt_shl,VpSlt_pred,VsSlt_pred,RhoSlt_pred,VpShl_pred,VsShl_pred,RhoShl_pred,VpSlt_filt,VsSlt_filt,RhoSlt_filt,VpShl_filt,VsShl_filt,RhoShl_filt,VpSlt_trend,VsSlt_trend,RhoSlt_trend,VpShl_trend,VsShl_trend,RhoShl_trend,PhiT_filt_slt,Slt_por_pred,Slt_por_trend,PhiT_filt_shl,Shl_por_pred,Shl_por_trend);

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
MinRho=2.0;
MaxRho=3.0;
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
Width=1000; %double monitor 1000
LengthSh=725; %double monitor 725
LengthLg=1600; %double monitor 1600

%% FIGURES

% FIGURE 1
figure1 = figure;
set(figure1, 'Position', [0 LengthSh+(0.1*LengthSh) Width LengthSh])
% Create subplot
subplot1 = subplot(2,2,1,'Parent',figure1,'XTick',[MinI:Step:MaxI],'YTick',[MinG:Step:MaxG]);
box(subplot1,'on');
hold(subplot1,'all');
xlim(subplot1,[MinI MaxI]);
ylim(subplot1,[MinG MaxG]);
 
% Create plot
%plot(A_Shl_Gas_iso*Iscalar,B_Shl_Gas_iso*Gscalar,'.','MarkerSize',30,'Parent',subplot1,'color','red');
%plot(A_Shl_Wtr_iso*Iscalar,B_Shl_Wtr_iso*Gscalar,'.','MarkerSize',30,'Parent',subplot1,'color','blue');
plot(A_Pred_gas_iso*Iscalar,B_Pred_gas_iso*Gscalar,'o','MarkerSize',10,'Parent',subplot1,'MarkerFaceColor','red','LineWidth',3,'MarkerEdgeColor','black');
plot(A_Pred_wtr_iso*Iscalar,B_Pred_wtr_iso*Gscalar,'o','MarkerSize',10,'Parent',subplot1,'MarkerFaceColor','blue','LineWidth',3,'MarkerEdgeColor','black');
 
 
% Plot Axes
plot(XAxis(:,1),XAxis(:,2),'Parent',subplot1,'LineWidth',2,'color','black');
plot(YAxis(:,1),YAxis(:,2),'Parent',subplot1,'LineWidth',2,'color','black');
 
% Create title
title('Intercept v Gradient (coloured by saturation)');
 
% Create subplot
subplot2 = subplot(2,2,2,'Parent',figure1,'XTick',[MinI:Step:MaxI],'YTick',[MinG:Step:MaxG]);
box(subplot2,'on');
hold(subplot2,'all');
xlim(subplot2,[MinI MaxI]);
ylim(subplot2,[MinG MaxG]);

% Create plot
%scatter(A_Shl_Gas_iso*Iscalar,B_Shl_Gas_iso*Gscalar,1000,por,'.','Parent',subplot2);
scatter(A_Pred_gas_iso*Iscalar,B_Pred_gas_iso*Gscalar,300,Porosity_predicted,'s','filled','Parent',subplot2);


% Plot Axes
plot(XAxis(:,1),XAxis(:,2),'Parent',subplot2,'LineWidth',2,'color','black');
plot(YAxis(:,1),YAxis(:,2),'Parent',subplot2,'LineWidth',2,'color','black');
 
% Create title
title('Intercept v Gradient (coloured by porosity)');
 
% Create subplot
subplot3 = subplot(2,2,3,'Parent',figure1,'XTick',[MinI:Step:MaxI],'YTick',[MinG:Step:MaxG]);
box(subplot3,'on');
hold(subplot3,'all');
xlim(subplot3,[MinI MaxI]);
ylim(subplot3,[MinG MaxG]);
 
% Create plot
%scatter(A_Shl_Gas_iso*Iscalar,B_Shl_Gas_iso*Gscalar,1000,Pves,'.','Parent',subplot3);
scatter(A_Pred_gas_iso*Iscalar,B_Pred_gas_iso*Gscalar,300,Pves_pred,'s','filled','Parent',subplot3);


% Plot Axes
plot(XAxis(:,1),XAxis(:,2),'Parent',subplot3,'LineWidth',2,'color','black');
plot(YAxis(:,1),YAxis(:,2),'Parent',subplot3,'LineWidth',2,'color','black');

% Create title
title('Intercept v Gradient (coloured by vertical effective stress)');

% Create subplot
subplot4 = subplot(2,2,4,'Parent',figure1,'XTick',[MinI:Step:MaxI],'YTick',[MinG:Step:MaxG]);
box(subplot4,'on');
hold(subplot4,'all');
xlim(subplot4,[MinI MaxI]);
ylim(subplot4,[MinG MaxG]);


% Create plot
%scatter(A_Shl_Gas_iso*Iscalar,B_Shl_Gas_iso*Gscalar,1000,vsh,'.','Parent',subplot4);
%scatter(A_Pred_gas_iso*Iscalar,B_Pred_gas_iso*Gscalar,300,Sand,'s','filled','Parent',subplot4);


% Plot Axes
plot(XAxis(:,1),XAxis(:,2),'Parent',subplot4,'LineWidth',2,'color','black');
plot(YAxis(:,1),YAxis(:,2),'Parent',subplot4,'LineWidth',2,'color','black');


% Create title
title('Intercept v Gradient (coloured by volume of shale)');

% FIGURE 2
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
%plot(RhoGas,Pves,'.','MarkerSize',20,'color','red','Parent',subplot1);
%plot(RhoWtr,Pves,'.','MarkerSize',20,'color','blue','Parent',subplot1);
plot(RhoGas_filt,Pves_filt_sand_gas,'o','MarkerSize',5,'Parent',subplot1,'MarkerFaceColor','red');
plot(RhoWtr_filt,Pves_filt_sand_wtr,'o','MarkerSize',5,'Parent',subplot1,'MarkerFaceColor','blue');
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
%plot(VpGas,Pves,'.','MarkerSize',20,'color','red','Parent',subplot2);
%plot(VpWtr,Pves,'.','MarkerSize',20,'color','blue','Parent',subplot2);
plot(VpGas_filt,Pves_filt_sand_gas,'o','MarkerSize',5,'Parent',subplot2,'MarkerFaceColor','red');
plot(VpWtr_filt,Pves_filt_sand_wtr,'o','MarkerSize',5,'Parent',subplot2,'MarkerFaceColor','blue');
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
%plot(VsGas,Pves,'.','MarkerSize',20,'color','red','Parent',subplot3);
%plot(VsWtr,Pves,'.','MarkerSize',20,'color','blue','Parent',subplot3);
plot(VsGas_filt,Pves_filt_sand_gas,'o','MarkerSize',5,'Parent',subplot3,'MarkerFaceColor','red');
plot(VsWtr_filt,Pves_filt_sand_wtr,'o','MarkerSize',5,'Parent',subplot3,'MarkerFaceColor','blue');
plot(VsGas_trend,Z,'--','Parent',subplot3,'color','red','LineWidth',2);
plot(VsWtr_trend,Z,'--','Parent',subplot3,'color','blue','LineWidth',2);
plot(VsGas_pred,Pves_pred-Pop-Buoy,'o','MarkerSize',10,'Parent',subplot3,'MarkerFaceColor','black','LineWidth',3,'MarkerEdgeColor','red');
plot(VsWtr_pred,Pves_pred-Pop,'o','MarkerSize',10,'Parent',subplot3,'MarkerFaceColor','black','LineWidth',3,'MarkerEdgeColor','blue');

% Create title
title('Vs reservoir v Effective Stress');
   

% FIGURE 3
figure3 = figure;
set(figure3, 'Position', [2.75*Width 0 0.75*Width LengthLg])
% Create subplot
subplot1 = subplot(1,3,1,'Parent',figure3,'Ydir','reverse','XTick',[MinRho:0.25:MaxRho],'YTick',[MinVes:500:MaxVes]);
box(subplot1,'on');
hold(subplot1,'all');
xlim(subplot1,[MinRho MaxRho]);
ylim(subplot1,[MinVes MaxVes]);

% Create plot
%scatter(RhoGas_filt,Pves_filt_sand_gas,1000,'.','Parent',subplot1);
scatter(RhoGas_filt,Pves_filt_sand_gas,50,Age_filt_gas,'o','filled','Parent',subplot1,'LineWidth',1,'MarkerEdgeColor','red');
%scatter(RhoWtr_filt,Pves_filt_sand_wtr,100,Age_filt_wtr,'o','filled','Parent',subplot1,'LineWidth',3,'MarkerEdgeColor','red');
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
scatter(VpGas_filt,Pves_filt_sand_gas,50,Age_filt_gas,'o','filled','Parent',subplot2,'LineWidth',1,'MarkerEdgeColor','red');
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
scatter(VsGas_filt,Pves_filt_sand_gas,50,Age_filt_gas,'o','filled','Parent',subplot3,'LineWidth',1,'MarkerEdgeColor','red');
plot(VsGas_trend,Z,'--','Parent',subplot3,'color','red','LineWidth',2);
plot(VsGas_pred,Pves_pred-Pop-Buoy,'o','MarkerSize',10,'Parent',subplot3,'MarkerFaceColor','black','LineWidth',4,'MarkerEdgeColor','red');


% Create title
title('Vs reservoir v Effective Stress');

% FIGURE 4
figure4 = figure;
set(figure4, 'Position', [3.5*Width 0 0.75*Width LengthLg])

% Create subplot
subplot1 = subplot(1,3,1,'Parent',figure4,'Ydir','reverse','XTick',[MinRho:0.25:MaxRho],'YTick',[MinVes:500:MaxVes]);
box(subplot1,'on');
hold(subplot1,'all');
xlim(subplot1,[MinRho MaxRho]);
ylim(subplot1,[MinVes MaxVes]);

% Create plot
scatter(RhoShl_filt,Pves_filt_shl,100,Age_filt_shl,'.','Parent',subplot1);
scatter(RhoSlt_filt,Pves_filt_slt,100,Age_filt_slt,'.','Parent',subplot1);
%scatter(RhoOB_filt,Pves_filt_OB,100,Age_filt,'o','filled','Parent',subplot1,'LineWidth',3,'MarkerEdgeColor','black');
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
scatter(VpShl_filt,Pves_filt_shl,100,Age_filt_shl,'.','Parent',subplot2);
scatter(VpSlt_filt,Pves_filt_slt,100,Age_filt_slt,'.','Parent',subplot2);
%scatter(VpOB_filt,Pves_filt_OB,100,Age_filt,'o','filled','Parent',subplot2,'LineWidth',3,'MarkerEdgeColor','black');
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

% Create plot/apps/gsc/matlab-library/development/maps
scatter(VsShl_filt,Pves_filt_shl,100,Age_filt_shl,'.','Parent',subplot3);
scatter(VsSlt_filt,Pves_filt_slt,100,Age_filt_slt,'.','Parent',subplot3);
%scatter(VsOB_filt,Pves_filt_OB,100,Age_filt,'o','filled','Parent',subplot3,'LineWidth',3,'MarkerEdgeColor','black');
plot(VsOB_trend,Z,'--','Parent',subplot3);
plot(VsOB_pred,Pves_pred-Pop,'o','MarkerSize',10,'Parent',subplot3,'MarkerFaceColor','black','LineWidth',4,'MarkerEdgeColor','green');


% Create title
title('Vs overburden v Effective Stress');



% FIGURE 5
figure5 = figure;
set(figure5, 'Position', [4.25*Width 0 0.75*Width LengthLg])

% Create subplot
subplot1 = subplot(1,2,1,'Parent',figure5,'Ydir','reverse','XTick',[0:0.05:0.4],'YTick',[MinVes:500:MaxVes]);
box(subplot1,'on');
hold(subplot1,'all');
xlim(subplot1,[0 0.4]);
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
%scatter(Por_filt,Pves_filt_sand,100,Facies_filt,'o','filled','Parent',subplot2,'LineWidth',3,'MarkerEdgeColor','black');
plot(OB_por_trend,Z,'Parent',subplot2);
plot(OB_por_pred,Pves_pred-Pop,'o','MarkerSize',10,'Parent',subplot2,'MarkerFaceColor','red','LineWidth',4,'MarkerEdgeColor','black');


% Create title
title('Overburden Porosity v Effective Stress');

% FIGURE 6
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
plot(A_Rand_wtr_iso*Iscalar,B_Rand_wtr_iso*Gscalar,'o','M/apps/gsc/matlab-library/development/mapsarkerSize',2,'Parent',subplot1,'MarkerFaceColor','blue','LineWidth',1,'MarkerEdgeColor','blue');
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
%plot(0:MaxAngle,RPPshg,'color','black','LineWidth',1,'Parent',subplot2);
%plot(0:MaxAngle,RPPsig,'color','black','LineWidth',1,'Parent',subplot2);
%plot(0:MaxAngle,RPPsig(:,INT),'color','red','LineWidth',4,'Parent',subplot2);
%plot(0:MaxAngle,RPPshg(:,INT),'color','red','LineWidth',4,'Parent',subplot2);
%plot(0:MaxAngle,RPPsiw(:,INT),'color','blue','LineWidth',4,'Parent',subplot2);
%plot(0:MaxAngle,RPPshw(:,INT),'color','blue','LineWidth',4,'Parent',subplot2);
plot(0:MaxAngle,RPPpred_gas,'color','red','LineWidth',2,'Parent',subplot2);
plot(0:MaxAngle,RPPpred_wtr,'color','blue','LineWidth',2,'Parent',subplot2);
plot(0:MaxAngle,RPPpred_shlslt,'color','green','LineWidth',2,'Parent',subplot2);
plot(0:MaxAngle,RPPpred_sltshl,'color','cyan','LineWidth',2,'Parent',subplot2);
plot(0:MaxAngle,RPPpred_gwc,'color','magenta','LineWidth',2,'Parent',subplot2);
plot(0:MaxAngle,FS_Gas,'color','red','LineWidth',2,'LineStyle',':','Parent',subplot2);
plot(0:MaxAngle,FS_Wtr,'color','blue','LineWidth',2,'LineStyle',':','Parent',subplot2);
%plot(0:MaxAngle,RPPshuey_Slt_Gas_aniso(:,INT),'color','red','LineWidth',4,'LineStyle','--','Parent',subplot2);
%plot(0:MaxAngle,RPPshuey_Shl_Gas_aniso(:,INT),'color','red','LineWidth',4,'LineStyle','--','Parent',subplot2);
%plot(0:MaxAngle,RPPshuey_Slt_Wtr_aniso(:,INT),'color','blue','LineWidth',4,'LineStyle','--','Parent',subplot2);
%plot(0:MaxAngle,RPPshuey_Shl_Wtr_aniso(:,INT),'color','blue','LineWidth',4,'LineStyle','--','Parent',subplot2);
plot(0:MaxAngle,RPPshuey_Pred_gas_aniso,'color','red','LineWidth',2,'LineStyle','--','Parent',subplot2);
plot(0:MaxAngle,RPPshuey_Pred_wtr_aniso,'color','blue','LineWidth',2,'LineStyle','--','Parent',subplot2);
plot(0:MaxAngle,RPPshuey_Pred_shlslt_aniso,'color','green','LineWidth',2,'LineStyle','--','Parent',subplot2);
plot(0:MaxAngle,RPPshuey_Pred_sltshl_aniso,'color','cyan','LineWidth',2,'LineStyle','--','Parent',subplot2);
%plot(0:MaxAngle,RPP,'color','yellow','LineWidth',4,'Parent',subplot2);

% Create title
title('Reflectivity vs Angle');

% Figure 7
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


% Figure 8
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
well_data{i_file}.curves_dec


% Create Figure
figure9 = figure;
set(figure9, 'Position', [0 LengthSh+(0.1*LengthSh) Width LengthSh])

% Create plot
scatterhist(AI_r,VpVsR_r,'Group',Lith,'LineStyle',{'-','-','-','-'},'LineWidth',[2,2,2,2]);

% Create title
title('AI v Vp/Vs Ratio');

% Figure 10
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


