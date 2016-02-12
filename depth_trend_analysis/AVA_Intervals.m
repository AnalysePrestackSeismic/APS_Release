function AVA_Intervals(Zpred,Zwc,Pop,Sand,Diag,GA,Shale,FitType,NI)
%% Rock Physics Database and Depth Trend Modelling Software
% Written by Matt Bolton April 2015
% Designed to incorporate all the rock properties of the sands, silts and mudrocks penetrated in the Tanzania offshore margin.
% Depth trends are then applied to the data after filtering for specific sand and caprock facies e.g. Age, Volume of Shale, Diagenesis, etc
% Seismic interfaces for can then be modelled at any depth based on trends fitted. 
% Reflectivity with angle is estimated for the isotropic case using the full Zoeppritz solution.  An anisotropic case is also included based on the Shuey approximation.

%% DATA READ
%Read in txt file containing Vp, Vs and Rho averages for reservoir and overburden intervals
%fileID = fopen('RockPhysicsAverageSets.txt')
[well, int, Zsb, Zml, por, vsh, VpGas, VsGas, RhoGas, VpWtr, VsWtr, RhoWtr, VpShl, VsShl, RhoShl, VpSlt, VsSlt, RhoSlt, P, T, Facies, Age, phit] = textread('RockPhysicsAverageSets2.txt','%s %s %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f ','headerlines',1);
[SeisI SeisG] = textread('BackgroundI-G.txt','%f %f');
%fileID = fopen('RockPhysicsAverageSets.txt');
%C = textscan(fileID,'%s %s %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f','headerLines',1);

%% PROSPECT INPUTS
% Geologically constrained depth trend analysis
% Key criteria for selecting the most appropriate sands and caprock for the prospect being modelled

% EXAMPLE PROSPECT ENTRY:
% AVA_Intervals(Zpred,Zwc,Pop,Sand,Diag,GA,Shale,FitType,NI)
% AVA_Intervals(1770,1500,0,0.05,1,60,1,2,5000)
%Zpred=input('Enter prospect TVDml(m)'); % Input prospect depth
%Zwc=input('Enter prospect water depth (m)'); % Water depth for pressure
%Pop=input('Overpressure expected in sand (psi); 0 = normal'); % Overpressure
%Sand=input('Enter Reservoir shale content; 0 = Clean Sand, 1 = Shale'); % Reservoir facies
%Diag=input('Enter expected Sand diagenesis; 1 = unconsolidated, 2 = consolidated, 3 = cemented'); % Diagenesis
%GA=input('Enter expected reservoir age (Ma)'); % Geological age
%Shale=input('Enter expected overburden; 0=Shl or 1=Slt'); % Overburden type
%FitType=input('Enter depth trend fitting method; 1=Data driven,2=Constrained, 3=Rock Physics Model;
%NI=input('Enter number of monte carlo simulations');

%% BASIC INPUTS
% Max incidence angle to calculate AVA over
MaxAngle=45; %degrees

% Depth or pressure iteration
MaxZ=10000; % max depth (or pressure (psi)) to calculate to
Z=0:1:MaxZ-1; % Z index with 1m sampling (requires all input depths are at this precision)
Pco=500; % Only include data in depth trend analysis which is at close to hydrostatic pressure (i.e. < Pco psi overpressued)

% Rock Physics Modelling
Clay=Sand; % fraction; Volume of clay 
Feldspar=0; % fraction; Volume of feldpar
Calcite=0; % fraction; Volume of calcite
Por_min=0.0; % minimum porosity at infinite pressure
%phiET=(-0.87*Sand)+0.87; % fraction of shale in pore space; 0 = shale matrix, 1 = all shale in pore space
phiET=0.25;
fcement=0.03; % cement fraction
PhiC=0.4; % Clean Sand critical porosity
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
Pp=P*0.006895; % convert Pore pressure from psi to Mpa

%% ROCK PHYSICS DATABASE CALCULATIONS
% Average Vp, Vs and Rho for each reservoir interval penetrated with
% consistent fluid substitution to 100% water and 95% gas.

% Shuey 3 terms approximation
% Equations to solve for Zero offset reflectivity R0, AVA gradient G and
% AVA curvature C.
% Modified set of equations to take velocity anisotropy into account

% Thomsen Anisotropy Parameter input
delta1=0.2; % Overburden
delta2=0; % Reservoir
epsilon1=0.1; % Overburden
epsilon2=0; % Reservoir
Ddelta=delta2-delta1; % form for use in anisotropic equation
Deps=epsilon2-epsilon1; % form for use in anisotropic equation

% Zoeppritz solution
% Overburden-Sand interface reflectivity
[RPPshg,RPSshg,TPPshg,TPSshg]=Zoeppritz(MaxAngle,VpShl,VsShl,RhoShl,VpGas,VsGas,RhoGas);
[RPPsig,RPSsig,TPPsig,TPSsig]=Zoeppritz(MaxAngle,VpSlt,VsSlt,RhoSlt,VpGas,VsGas,RhoGas);
[RPPshw,RPSshw,TPPshw,TPSshw]=Zoeppritz(MaxAngle,VpShl,VsShl,RhoShl,VpWtr,VsWtr,RhoWtr);
[RPPsiw,RPSsiw,TPPsiw,TPSsiw]=Zoeppritz(MaxAngle,VpSlt,VsSlt,RhoSlt,VpWtr,VsWtr,RhoWtr);

% Gas-Water contact interface reflectivity
[RPPgwc,RPSgwc,TPPgwc,TPSgwc]=Zoeppritz(MaxAngle,VpGas,VsGas,RhoGas,VpWtr,VsWtr,RhoWtr);

% Background shale-silt interface reflectivity
[RPPshsl,RPSshsl,TPPshsl,TPSshsl]=Zoeppritz(MaxAngle,VpShl,VsShl,RhoShl,VpSlt,VsSlt,RhoSlt);
RPPslsh=-RPPshsl;

% Shuey 3 term approximation
% Overburden-Sand interface Shuey coefficients
[A_Slt_Gas_iso,B_Slt_Gas_iso,C_Slt_Gas_iso,A_Slt_Gas_aniso,B_Slt_Gas_aniso,C_Slt_Gas_aniso,RPPshuey_Slt_Gas_iso,RPPshuey_Slt_Gas_aniso] = Shuey(MaxAngle,VpSlt,VsSlt,RhoSlt,VpGas,VsGas,RhoGas,delta1,delta2,epsilon1,epsilon2);
[A_Slt_Wtr_iso,B_Slt_Wtr_iso,C_Slt_Wtr_iso,A_Slt_Wtr_aniso,B_Slt_Wtr_aniso,C_Slt_Wtr_aniso,RPPshuey_Slt_Wtr_iso,RPPshuey_Slt_Wtr_aniso] = Shuey(MaxAngle,VpSlt,VsSlt,RhoSlt,VpWtr,VsWtr,RhoWtr,delta1,delta2,epsilon1,epsilon2);
[A_Shl_Gas_iso,B_Shl_Gas_iso,C_Shl_Gas_iso,A_Shl_Gas_aniso,B_Shl_Gas_aniso,C_Shl_Gas_aniso,RPPshuey_Shl_Gas_iso,RPPshuey_Shl_Gas_aniso] = Shuey(MaxAngle,VpShl,VsShl,RhoShl,VpGas,VsGas,RhoGas,delta1,delta2,epsilon1,epsilon2);
[A_Shl_Wtr_iso,B_Shl_Wtr_iso,C_Shl_Wtr_iso,A_Shl_Wtr_aniso,B_Shl_Wtr_aniso,C_Shl_Wtr_aniso,RPPshuey_Shl_Wtr_iso,RPPshuey_Shl_Wtr_aniso] = Shuey(MaxAngle,VpShl,VsShl,RhoShl,VpWtr,VsWtr,RhoWtr,delta1,delta2,epsilon1,epsilon2);

% Gas-Water contact interface Shuey coefficients (isotropic case assumed)
[A_Gas_Wtr_iso,B_Gas_Wtr_iso,C_Gas_Wtr_iso,A_Gas_Wtr_aniso,B_Gas_Wtr_aniso,C_Gas_Wtr_aniso,RPPshuey_Gas_Wtr_iso,RPPshuey_Gas_Wtr_aniso] = Shuey(MaxAngle,VpGas,VsGas,RhoGas,VpWtr,VsWtr,RhoWtr,0,0,0,0);

% Background shale-silt interface Shuey coefficients
[A_Shl_Slt_iso,B_Shl_Slt_iso,C_Shl_Slt_iso,A_Shl_Slt_aniso,B_Shl_Slt_aniso,C_Shl_Slt_aniso,RPPshuey_Shl_Slt_iso,RPPshuey_Shl_Slt_aniso] = Shuey(MaxAngle,VpShl,VsShl,RhoShl,VpSlt,VsSlt,RhoSlt,delta1,delta2,epsilon1,epsilon2);

% silt-shale interface is the negative form of shale-silt
A_Slt_Shl_iso = -A_Shl_Slt_iso;
B_Slt_Shl_iso = -B_Shl_Slt_iso;
C_Slt_Shl_iso = -C_Shl_Slt_iso;

A_Slt_Shl_aniso = -A_Shl_Slt_aniso;
B_Slt_Shl_aniso = -B_Shl_Slt_aniso;
C_Slt_Shl_aniso = -C_Shl_Slt_aniso;

% Bulk rock properties 
[AI_Gas,VpVsR_Gas,Mu_Gas,K_Gas] = RockProps(VpGas,VsGas,RhoGas);
[AI_Wtr,VpVsR_Wtr,Mu_Wtr,K_Wtr] = RockProps(VpWtr,VsWtr,RhoWtr);
[AI_Shl,VpVsR_Shl,Mu_Shl,K_Shl] = RockProps(VpShl,VsShl,RhoShl);
[AI_Slt,VpVsR_Slt,Mu_Slt,K_Slt] = RockProps(VpSlt,VsSlt,RhoSlt);

%% TEMPERATURE & PRESSURE CALCULATIONS
% All depth trends should be calculated wrt to vertical effective stress if
% basin is not hydrostatically pressured.

% EDITABLE CONTENT SEE LINE 56

% CALCULATIONS
% Batzle and wang calculations
[Kreuss,rhoeff,Kvoigt,vpb,rhob,Kb,vpo,rhoo,Ko,vpg,rhog,Kg,gor]=flprop(method,sal,og,gg,gor,giib,giio,Pp,T,So,Sg);

RhoW=1;
%RhoW=mean(rhob); % Used as average overburden brine properties for prediction
% Lithostatic Trend
% Overburden Weights - must add up to one.  Vary depending on expected overburden lithology
Shl=1;
Slt=0;
Snd=0;

% Average overburden profile
RhoOB=(Shl*RhoShl)+(Slt*RhoSlt); % Overburden density data
VpOB=(Shl*VpShl)+(Slt*VpSlt); % Overburden Vp data
VsOB=(Shl*VsShl)+(Slt*VsSlt); % Overburden Vs data

ZOB=Zml;
ini(MaxZ,1)=0;

% Exponential curve fit to average overburden density points
%beta = [1000 1000]; % initial model
%x=RhoOB;
%betaRhoOB = nlinfit(ZOB,x,@ModelFit,beta);
%RhoL = (betaRhoOB(1)*log(Z)-betaRhoOB(2));

dRhoOB=ZOB;
GRhoOB=log((RhoOB-RhoMin)/(RhoMin-RhoI));
mRhoOB=GRhoOB\dRhoOB;
betaRhoOB=real(-1/mRhoOB);

GrOB=exp(-betaRhoOB*ZOB);
d=RhoOB-RhoI; % known data
mx=RhoI; % max shale density
mRhoMudOB=ModelFitv2(d,GrOB,betaRhoOB,mx);
RhoL = RhoI-((RhoI-mRhoMudOB)*exp(-betaRhoOB*Z));

% Integration of density trend (RhoL) to define average lithostatic trend
for z = 1:MaxZ-2
    Plith(z+1)=1000*(Z(z+2)-Z(z+1))*(RhoL(z+1));
    ini(z+1)=ini(z)+Plith(z+1);
end

ini(MaxZ)=ini(MaxZ-1); % account for array indexing

%Convert from Pa to PSI: PSI = 0.000145Pa and add water column.
for aa=1:length(Zml);
    %Create Hydrostatic curve
    Phydro(aa,1)=1000*0.000145*9.81*((RhoW*Zml(aa))+(RhoSea*Zsb(aa)));
    for bb=1:length(Z);
    Plith_f(aa,bb)=(0.000145*9.81)*(ini(bb,1)+(Zsb(aa)*RhoSea*1000));

end
end

% Calculate Vertical Effective Stress for each point in the database using lithostatic curve (Plith) and Formation Pressure (P)
for y = 1:length(Zml);
    for v = 1:MaxZ       
        if Z(v) == Zml(y)
        Pv = Plith_f(y,v)-P(y);
        Pves(y,1)=Pv;       
        end
    end
end

% Overpressure calculated for each point in the database
Pover=P-Phydro;

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
COage=66; % Cut off Age which overburden rock properties show different compaction trends
if GA>COage % 66Ma Base Tertiary
    GA_low=COage;
    GA_hi=200; % Base Jurassic
else
    GA_low=0;
    GA_hi=COage;
end

% Filter sand points based on volume of shale and diagenesis
RhoGas_filt=RhoGas(Facies==Diag & vsh <= Sand_hi & vsh >= Sand_low & Pover < Pco);
VpGas_filt=VpGas(Facies==Diag & vsh <= Sand_hi & vsh >= Sand_low & Pover < Pco);
VsGas_filt=VsGas(Facies==Diag & vsh <= Sand_hi & vsh >= Sand_low & Pover < Pco);
RhoWtr_filt=RhoWtr(Facies==Diag & vsh <= Sand_hi & vsh >= Sand_low & Pover < Pco);
VpWtr_filt=VpWtr(Facies==Diag & vsh <= Sand_hi & vsh >= Sand_low & Pover < Pco);
VsWtr_filt=VsWtr(Facies==Diag & vsh <= Sand_hi & vsh >= Sand_low & Pover < Pco);
Por_filt=por(Facies==Diag & vsh <= Sand_hi & vsh >= Sand_low & Pover < Pco);
Facies_filt=Facies(Facies==Diag & vsh <= Sand_hi & vsh >= Sand_low & Pover < Pco);
vsh_filt=vsh(Facies==Diag & vsh <= Sand_hi & vsh >= Sand_low & Pover < Pco);

% Filter shale points based on facies and geological age
if Shale == 0
    RhoOB_filt=RhoShl(Age>GA_low & Age<GA_hi & Pover < Pco);
    VpOB_filt=VpShl(Age>GA_low & Age<GA_hi & Pover < Pco);
    VsOB_filt=VsShl(Age>GA_low & Age<GA_hi & Pover < Pco);
    % Smectite-Illite constraints
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
    RhoOB_filt=RhoSlt(Age>GA_low & Age<GA_hi & Pover < Pco);
    VpOB_filt=VpSlt(Age>GA_low & Age<GA_hi & Pover < Pco);
    VsOB_filt=VsSlt(Age>GA_low & Age<GA_hi & Pover < Pco);
    
    if GA >= COage
        % Smectite-Illite constraints
        % In silt overburden case mudstone minerals are Voight mixed with quartz
            RhoMud=(0.8*RhoI)+(0.2*RhoQ);
            VpMud=(0.8*VpI)+(0.2*VpQ);
            VsMud=(0.8*VsI)+(0.2*VsQ);
    else
            RhoMud=(0.8*RhoS)+(0.2*RhoQ);
            VpMud=(0.8*VpS)+(0.2*VpQ);
            VsMud=(0.8*VsS)+(0.2*VsQ);
    end
end


% Filter Overburden facies on geological age
RhoShl_filt=RhoShl(Age>GA_low & Age<GA_hi & Pover < Pco);
VpShl_filt=VpShl(Age>GA_low & Age<GA_hi & Pover < Pco);
VsShl_filt=VsShl(Age>GA_low & Age<GA_hi & Pover < Pco);
RhoSlt_filt=RhoSlt(Age>GA_low & Age<GA_hi & Pover < Pco);
VpSlt_filt=VpSlt(Age>GA_low & Age<GA_hi & Pover < Pco);
VsSlt_filt=VsSlt(Age>GA_low & Age<GA_hi & Pover < Pco);
Age_filt=Age(Age>GA_low & Age<GA_hi & Pover < Pco);

% Filter vertical effective stress points to get corresponding data for
% curve fitting to overburden
Pves_filt_OB=Pves(Age>GA_low & Age<GA_hi & Pover < Pco);

% Filter vertical effective stress points to get corresponding data for
% curve fitting to sand
Pves_filt_sand=Pves(Facies==Diag & vsh <= Sand_hi & vsh >= Sand_low & Pover < Pco);


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
% Least squares solution to Por=Por_min-(Por_min-PhiC)*exp(-betaPor*Pves)
% solving for betaPor

dPor=Pves_filt_sand;
GPor=log((Por_filt-Por_min)/(Por_min-PhiC));
mPor=GPor\dPor;
betaPor=real(-1/mPor);

Porosity_predicted = Por_min-((Por_min-PhiC)*exp(-betaPor*(Pves_pred)));
Por_trend = Por_min-((Por_min-PhiC)*exp(-betaPor*Z));

% Similar approach for caprock however as PHIT not always present, RHOB
% used as a proxy
dRho=Pves_filt_OB;
GRho=log((RhoOB_filt-RhoMin)/(RhoMin-RhoMud));
mRho=GRho\dRho;
betaRho=real(-1/mRho);

% Gradient matrix for all subsequent least squares depth trend equations
% Porosity depth trend used as a known compaction gradient
Gr=exp(-betaRho*Pves_filt_OB); % Compaction gradient for caprock
G=exp(-betaPor*Pves_filt_sand); % Compaction gradient for sands

% Gas Sand Rock Properties

% Rho Gas
d=RhoGas_filt-RhoQ; % known data
mx=RhoQ; % Quartz density
mRhoGas=ModelFitv2(d,G,betaPor,mx);
RhoGas_pred = RhoQ-((RhoQ-mRhoGas)*exp(-betaPor*Pves_pred));
RhoGas_trend = RhoQ-((RhoQ-mRhoGas)*exp(-betaPor*Z));

% Vp Gas
d=VpGas_filt-VpQ; % known data
mx=VpQ; % Quartz Vp
mVpGas=ModelFitv2(d,G,betaPor,mx);
VpGas_pred = VpQ-((VpQ-mVpGas)*exp(-betaPor*Pves_pred));
VpGas_trend = VpQ-((VpQ-mVpGas)*exp(-betaPor*Z));

% Vs Gas
d=VsGas_filt-VsQ; % known data
mx=VsQ; % Quartz Vs
mVsGas=ModelFitv2(d,G,betaPor,mx);
VsGas_pred = VsQ-((VsQ-mVsGas)*exp(-betaPor*Pves_pred));
VsGas_trend = VsQ-((VsQ-mVsGas)*exp(-betaPor*Z));

% Water Sand Rock Properties

% Rho Water
d=RhoWtr_filt-RhoQ; % known data
mx=RhoQ; % Quartz density
mRhoWtr=ModelFitv2(d,G,betaPor,mx);
RhoWtr_pred = RhoQ-((RhoQ-mRhoWtr)*exp(-betaPor*Pves_pred));
RhoWtr_trend = RhoQ-((RhoQ-mRhoWtr)*exp(-betaPor*Z));

% Vp Water
d=VpWtr_filt-VpQ; % known data
mx=VpQ; % Quartz Vp
mVpWtr=ModelFitv2(d,G,betaPor,mx);
VpWtr_pred = VpQ-((VpQ-mVpWtr)*exp(-betaPor*Pves_pred));
VpWtr_trend = VpQ-((VpQ-mVpWtr)*exp(-betaPor*Z));

% Vs Water
d=VsWtr_filt-VsQ; % known data
mx=VsQ; % Quartz Vs
mVsWtr=ModelFitv2(d,G,betaPor,mx);
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
mRhoMud=ModelFitv2(d,Gr,betaRho,mx);
RhoOB_pred = RhoMud-((RhoMud-mRhoMud)*exp(-betaRho*Pves_pred));
RhoOB_trend = RhoMud-((RhoMud-mRhoMud)*exp(-betaRho*Z));

d=VpOB_filt-VpMud; % known data
mx=VpMud; % Quartz density
mVpMud=ModelFitv2(d,Gr,betaRho,mx);
VpOB_pred = VpMud-((VpMud-mVpMud)*exp(-betaRho*Pves_pred));
VpOB_trend = VpMud-((VpMud-mVpMud)*exp(-betaRho*Z));

d=VsOB_filt-VsMud; % known data
mx=VsMud; % Quartz density
mVsMud=ModelFitv2(d,Gr,betaRho,mx);
VsOB_pred = VsMud-((VsMud-mVsMud)*exp(-betaRho*Pves_pred));
VsOB_trend = VsMud-((VsMud-mVsMud)*exp(-betaRho*Z));

% Caprock overpressure modelling
% TZA data at Jodari suggests density measurements not effected by
% overpressure.  This corroborates published research and datasets.
% Therefore for the caprock we assume Rho_pred(Pves)=Rho_pred(Pves-Pover)

[VpOB_pred]=Overpressure(a_mud,b_mud,c_mud,Pves_pred,Pop,0,VpOB_pred,mVpMud,VpMud);
[VsOB_pred]=Overpressure(a_mud,b_mud,c_mud,Pves_pred,Pop,0,VsOB_pred,mVsMud,VsMud);


% Caprock property trends

% Similar approach for caprock however as PHIT not always present, RHOB
% used as a proxy

if GA >= COage
            RhoMud=RhoI;
            VpMud=VpI;
            VsMud=VsI;
        else
            RhoMud=RhoS;
            VpMud=VpS;
            VsMud=VsS;
end

dShRho=Pves_filt_OB;
GShRho=log((RhoShl_filt-RhoMin)/(RhoMin-RhoMud));
mShRho=GShRho\dShRho;
betaShRho=real(-1/mShRho);
GSh=exp(-betaShRho*Pves_filt_OB); % Compaction gradient for caprock

d=RhoShl_filt-RhoMud; % known data
mx=RhoMud; % Quartz density
mRhoMud=ModelFitv2(d,GSh,betaShRho,mx);
RhoShl_pred = RhoMud-((RhoMud-mRhoMud)*exp(-betaShRho*Pves_pred));
RhoShl_trend = RhoMud-((RhoMud-mRhoMud)*exp(-betaShRho*Z));

d=VpShl_filt-VpMud; % known data
mx=VpMud; % Quartz density
mVpMud=ModelFitv2(d,GSh,betaShRho,mx);
VpShl_pred = VpMud-((VpMud-mVpMud)*exp(-betaShRho*Pves_pred));
VpShl_trend = VpMud-((VpMud-mVpMud)*exp(-betaShRho*Z));

d=VsShl_filt-VsMud; % known data
mx=VsMud; % Quartz density
mVsMud=ModelFitv2(d,GSh,betaShRho,mx);
VsShl_pred = VsMud-((VsMud-mVsMud)*exp(-betaShRho*Pves_pred));
VsShl_trend = VsMud-((VsMud-mVsMud)*exp(-betaShRho*Z));

% Overpressure modelling
[VpShl_pred]=Overpressure(a_mud,b_mud,c_mud,Pves_pred,Pop,0,VpShl_pred,mVpMud,VpMud);
[VsShl_pred]=Overpressure(a_mud,b_mud,c_mud,Pves_pred,Pop,0,VsShl_pred,mVsMud,VsMud);

    if GA >= COage
        % Smectite-Illite constraints
        % In silt overburden case mudstone minerals are Voight mixed with quartz
            RhoMud=(0.8*RhoI)+(0.2*RhoQ);
            VpMud=(0.8*VpI)+(0.2*VpQ);
            VsMud=(0.8*VsI)+(0.2*VsQ);
    else
            RhoMud=(0.8*RhoS)+(0.2*RhoQ);
            VpMud=(0.8*VpS)+(0.2*VpQ);
            VsMud=(0.8*VsS)+(0.2*VsQ);
    end


% Silt overburden rock properties
dSlRho=Pves_filt_OB;
GSlRho=log((RhoSlt_filt-RhoMin)/(RhoMin-RhoMud));
mSlRho=GSlRho\dSlRho;
betaSlRho=real(-1/mSlRho);
GSl=exp(-betaSlRho*Pves_filt_OB); % Compaction gradient for caprock

d=RhoSlt_filt-RhoMud; % known data
mx=RhoMud; % Quartz density
mRhoMud=ModelFitv2(d,GSl,betaSlRho,mx);
RhoSlt_pred = RhoMud-((RhoMud-mRhoMud)*exp(-betaSlRho*Pves_pred));
RhoSlt_trend = RhoMud-((RhoMud-mRhoMud)*exp(-betaSlRho*Z));

d=VpSlt_filt-VpMud; % known data
mx=VpMud; % Quartz density
mVpMud=ModelFitv2(d,GSl,betaSlRho,mx);
VpSlt_pred = VpMud-((VpMud-mVpMud)*exp(-betaSlRho*Pves_pred));
VpSlt_trend = VpMud-((VpMud-mVpMud)*exp(-betaSlRho*Z));

d=VsSlt_filt-VsMud; % known data
mx=VsMud; % Quartz density
mVsMud=ModelFitv2(d,GSl,betaSlRho,mx);
VsSlt_pred = VsMud-((VsMud-mVsMud)*exp(-betaSlRho*Pves_pred));
VsSlt_trend = VsMud-((VsMud-mVsMud)*exp(-betaSlRho*Z));

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
[PhiC,Rhodc,Vpdc,Vsdc,Rhosc,Vpsc,Vssc,Phis,Rhods,Vpds,Vsds,Rhoss,Vpss,Vsss]=cementedsand(Clay,Feldspar,Calcite,PhiC,fcement,kf,rhof,shearfact);
AIsc=Vpsc.*(Rhosc/1000);
VpVs_R_sc=Vpsc./Vssc;
AIss=Vpss.*(Rhoss/1000);
VpVs_R_ss=Vpss./Vsss;


%% MONTE CARLO SIMULATION

% Calculate residuals between model and filtered data points for Vp, Vs and Rho. Residuals are used as the uncertainty distribution as this removes the compaction effect.
% Generate a normal random distribution of rock properties based on the covariance between the residuals
[VpOB_r VsOB_r RhoOB_r VpGas_r VsGas_r RhoGas_r Por_r_ob Por_r_gas] = RockPropDist(NI,1,MaxZ,Pves_filt_OB,Pves_filt_sand,VpOB_pred,VsOB_pred,RhoOB_pred,VpGas_pred,VsGas_pred,RhoGas_pred,VpOB_filt,VsOB_filt,RhoOB_filt,VpGas_filt,VsGas_filt,RhoGas_filt,VpOB_trend,VsOB_trend,RhoOB_trend,VpGas_trend,VsGas_trend,RhoGas_trend,ones(length(Pves_filt_OB),1),1,ones(length(Z),1),Por_filt,Porosity_predicted,Por_trend);
[VpOB_r VsOB_r RhoOB_r VpWtr_r VsWtr_r RhoWtr_r Por_r_ob Por_r_wtr] = RockPropDist(NI,1,MaxZ,Pves_filt_OB,Pves_filt_sand,VpOB_pred,VsOB_pred,RhoOB_pred,VpWtr_pred,VsWtr_pred,RhoWtr_pred,VpOB_filt,VsOB_filt,RhoOB_filt,VpWtr_filt,VsWtr_filt,RhoWtr_filt,VpOB_trend,VsOB_trend,RhoOB_trend,VpWtr_trend,VsWtr_trend,RhoWtr_trend,ones(length(Pves_filt_OB),1),1,ones(length(Z),1),Por_filt,Porosity_predicted,Por_trend);
[VpSlt_r VsSlt_r RhoSlt_r VpShl_r VsShl_r RhoShl_r Por_r_ob1 Por_r_ob2] = RockPropDist(NI,1,MaxZ,Pves_filt_OB,Pves_filt_sand,VpSlt_pred,VsSlt_pred,RhoSlt_pred,VpShl_pred,VsShl_pred,RhoShl_pred,VpSlt_filt,VsSlt_filt,RhoSlt_filt,VpShl_filt,VsShl_filt,RhoShl_filt,VpSlt_trend,VsSlt_trend,RhoSlt_trend,VpShl_trend,VsShl_trend,RhoShl_trend,ones(length(Pves_filt_OB),1),1,ones(length(Z),1),ones(length(Pves_filt_OB),1),1,ones(length(Z),1));

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
MinRho=1.5;
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
plot(A_Shl_Gas_iso*Iscalar,B_Shl_Gas_iso*Gscalar,'.','MarkerSize',30,'Parent',subplot1,'color','red');
plot(A_Shl_Wtr_iso*Iscalar,B_Shl_Wtr_iso*Gscalar,'.','MarkerSize',30,'Parent',subplot1,'color','blue');
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
scatter(A_Shl_Gas_iso*Iscalar,B_Shl_Gas_iso*Gscalar,1000,por,'.','Parent',subplot2);
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
scatter(A_Shl_Gas_iso*Iscalar,B_Shl_Gas_iso*Gscalar,1000,Pves,'.','Parent',subplot3);
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
scatter(A_Shl_Gas_iso*Iscalar,B_Shl_Gas_iso*Gscalar,1000,vsh,'.','Parent',subplot4);
scatter(A_Pred_gas_iso*Iscalar,B_Pred_gas_iso*Gscalar,300,Sand,'s','filled','Parent',subplot4);


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
plot(RhoGas,Pves,'.','MarkerSize',20,'color','red','Parent',subplot1);
plot(RhoWtr,Pves,'.','MarkerSize',20,'color','blue','Parent',subplot1);
plot(RhoGas_filt,Pves_filt_sand,'o','MarkerSize',8,'Parent',subplot1,'MarkerFaceColor','red','LineWidth',3,'MarkerEdgeColor','black');
plot(RhoWtr_filt,Pves_filt_sand,'o','MarkerSize',8,'Parent',subplot1,'MarkerFaceColor','blue','LineWidth',3,'MarkerEdgeColor','black');
plot(RhoGas_trend,Z,'--','Parent',subplot1,'color','red','LineWidth',2);
plot(RhoWtr_trend,Z,'--','Parent',subplot1,'color','blue','LineWidth',2);
plot(RhoGas_pred,Pves_pred-Pop-Buoy,'o','MarkerSize',10,'Parent',subplot1,'MarkerFaceColor','black','LineWidth',4,'MarkerEdgeColor','red');
plot(RhoWtr_pred,Pves_pred-Pop,'o','MarkerSize',10,'Parent',subplot1,'MarkerFaceColor','black','LineWidth',4,'MarkerEdgeColor','blue');

% Create title
title('Reservoir density vs Effective Stress');

% Create subplot
subplot2 = subplot(1,3,2,'Parent',figure2,'Ydir','reverse','XTick',[MinVp:1000:MaxVp],'YTick',[MinVes:500:MaxVes]);
box(subplot2,'on');
hold(subplot2,'all');
xlim(subplot2,[MinVp MaxVp]);
ylim(subplot2,[MinVes MaxVes]);

% Create plot
plot(VpGas,Pves,'.','MarkerSize',20,'color','red','Parent',subplot2);
plot(VpWtr,Pves,'.','MarkerSize',20,'color','blue','Parent',subplot2);
plot(VpGas_filt,Pves_filt_sand,'o','MarkerSize',8,'Parent',subplot2,'MarkerFaceColor','red','LineWidth',3,'MarkerEdgeColor','black');
plot(VpWtr_filt,Pves_filt_sand,'o','MarkerSize',8,'Parent',subplot2,'MarkerFaceColor','blue','LineWidth',3,'MarkerEdgeColor','black');
plot(VpGas_trend,Z,'--','Parent',subplot2,'color','red','LineWidth',2);
plot(VpWtr_trend,Z,'--','Parent',subplot2,'color','blue','LineWidth',2);
plot(VpGas_pred,Pves_pred-Pop-Buoy,'o','MarkerSize',10,'Parent',subplot2,'MarkerFaceColor','black','LineWidth',4,'MarkerEdgeColor','red');
plot(VpWtr_pred,Pves_pred-Pop,'o','MarkerSize',10,'Parent',subplot2,'MarkerFaceColor','black','LineWidth',4,'MarkerEdgeColor','blue');

% Create title
title('Vp reservoir v Effective Stress');

% Create subplot
subplot3 = subplot(1,3,3,'Parent',figure2,'Ydir','reverse','XTick',[MinVs:1000:MaxVs],'YTick',[MinVes:500:MaxVes]);
box(subplot3,'on');
hold(subplot3,'all');
xlim(subplot3,[MinVs MaxVs]);
ylim(subplot3,[MinVes MaxVes]);

% Create plot
plot(VsGas,Pves,'.','MarkerSize',20,'color','red','Parent',subplot3);
plot(VsWtr,Pves,'.','MarkerSize',20,'color','blue','Parent',subplot3);
plot(VsGas_filt,Pves_filt_sand,'o','MarkerSize',8,'Parent',subplot3,'MarkerFaceColor','red','LineWidth',3,'MarkerEdgeColor','black');
plot(VsWtr_filt,Pves_filt_sand,'o','MarkerSize',8,'Parent',subplot3,'MarkerFaceColor','blue','LineWidth',3,'MarkerEdgeColor','black');
plot(VsGas_trend,Z,'--','Parent',subplot3,'color','red','LineWidth',2);
plot(VsWtr_trend,Z,'--','Parent',subplot3,'color','blue','LineWidth',2);
plot(VsGas_pred,Pves_pred-Pop-Buoy,'o','MarkerSize',10,'Parent',subplot3,'MarkerFaceColor','black','LineWidth',4,'MarkerEdgeColor','red');
plot(VsWtr_pred,Pves_pred-Pop,'o','MarkerSize',10,'Parent',subplot3,'MarkerFaceColor','black','LineWidth',4,'MarkerEdgeColor','blue');

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
scatter(RhoGas,Pves,1000,Facies,'.','Parent',subplot1);
scatter(RhoGas_filt,Pves_filt_sand,100,Facies_filt,'o','filled','Parent',subplot1,'LineWidth',3,'MarkerEdgeColor','black');
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
scatter(VpGas,Pves,100,Facies,'o','filled','Parent',subplot2);
scatter(VpGas_filt,Pves_filt_sand,100,Facies_filt,'o','filled','Parent',subplot2,'LineWidth',3,'MarkerEdgeColor','black');
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
scatter(VsGas,Pves,1000,Facies,'.','Parent',subplot3);
scatter(VsGas_filt,Pves_filt_sand,100,Facies_filt,'o','filled','Parent',subplot3,'LineWidth',3,'MarkerEdgeColor','black');
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
scatter(RhoShl,Pves,1000,Age,'.','Parent',subplot1);
scatter(RhoSlt,Pves,1000,Age,'.','Parent',subplot1);
scatter(RhoOB_filt,Pves_filt_OB,100,Age_filt,'o','filled','Parent',subplot1,'LineWidth',3,'MarkerEdgeColor','black');
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
scatter(VpShl,Pves,1000,Age,'.','Parent',subplot2);
scatter(VpSlt,Pves,1000,Age,'.','Parent',subplot2);
scatter(VpOB_filt,Pves_filt_OB,100,Age_filt,'o','filled','Parent',subplot2,'LineWidth',3,'MarkerEdgeColor','black');
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
scatter(VsShl,Pves,1000,Age,'.','Parent',subplot3);
scatter(VsSlt,Pves,1000,Age,'.','Parent',subplot3);
scatter(VsOB_filt,Pves_filt_OB,100,Age_filt,'o','filled','Parent',subplot3,'LineWidth',3,'MarkerEdgeColor','black');
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
scatter(por,Pves,1000,vsh,'.','Parent',subplot1);
scatter(Por_filt,Pves_filt_sand,100,vsh_filt,'o','filled','Parent',subplot1,'LineWidth',3,'MarkerEdgeColor','black');
plot(Por_trend,Z,'Parent',subplot1);
plot(Porosity_predicted,Pves_pred-Pop,'o','MarkerSize',10,'Parent',subplot1,'MarkerFaceColor','red','LineWidth',4,'MarkerEdgeColor','black');


% Create title
title('Reservoir Porosity v Effective Stress');

% Create subplot
subplot2 = subplot(1,2,2,'Parent',figure5,'Ydir','reverse','XTick',[0:0.05:0.4],'YTick',[MinVes:500:MaxVes]);
box(subplot2,'on');
hold(subplot2,'all');
xlim(subplot2,[0 0.4]);
ylim(subplot2,[MinVes MaxVes]);

% Create plot
scatter(por,Pves,1000,Facies,'.','Parent',subplot2);
scatter(Por_filt,Pves_filt_sand,100,Facies_filt,'o','filled','Parent',subplot2,'LineWidth',3,'MarkerEdgeColor','black');
plot(Por_trend,Z,'Parent',subplot2);
plot(Porosity_predicted,Pves_pred-Pop,'o','MarkerSize',10,'Parent',subplot2,'MarkerFaceColor','red','LineWidth',4,'MarkerEdgeColor','black');


% Create title
title('Reservoir Porosity v Effective Stress');

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
scatter(Por_filt,VpWtr_filt,1000,Pves_filt_sand,'.','Parent',subplot1);
plot(Porosity_predicted,VpWtr_pred,'o','MarkerSize',10,'Parent',subplot1,'MarkerFaceColor','black','LineWidth',4,'MarkerEdgeColor','magenta');
plot(PhiC,Vpsc,'Parent',subplot1,'color','red');
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
scatter(AI_Wtr,VpVsR_Wtr,1000,Pves,'.','Parent',subplot2);
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



% Create Figure
figure9 = figure;
set(figure9, 'Position', [0 LengthSh+(0.1*LengthSh) Width LengthSh])

% Create plot
%scatterhist(AI_r,VpVsR_r,'Group',Lith,'LineStyle',{'-','-','-','-'},'LineWidth',[2,2,2,2]);

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

