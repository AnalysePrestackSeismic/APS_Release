%% AnalysePrestackSeismic
% A library of MATLAB code for analysing post-migration pre-stack seismic 
% data.

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

%% ------------------ SCRIPT DEFINITION ---------------------------------
% Defining global settings used by APS functions
% This file is loaded by all APS functions

%% Column numbers define output format of .mat_lite file. 
% Should make this a global format definition.  
pkey_loc = 1;                               % column numbers needs to be implemented Primary Key
skey_loc = 2;                               % Secondary Key
byte_loc = 3;                               % Byte location
skey_max_loc = 4;                           % Secondary Key Maximum
skey_inc_loc = 5;                           % Secondary Key Increment  
tkey_loc = 6;                               % Tertiary Key
tkey_max_loc = 7;                           % Tertiary Key Maximum
tkey_inc_loc = 8;                           % Tertiary Key Increment 