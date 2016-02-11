function [] = add_horizon_to_job_meta(job_meta_path,varargin)
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
%% Funtion to add a horizon to job meta file. For use in limiting the IG inversion when becomes unstable at depth.
% INPUTS:
%   job_meta_path = path of job meta file
%
%   if you have a 3D horizon (for 3D surveys only: varargin =
%   '1','horizon_path',(Note:horizon should be a text of the  format: IL XL Z as columns
%   if you want a constant horizon:    varargin{1} ='0','horizon',varargin{2} = 'horizon value in TWT or depth' 
%
%
% OUPUT:
%   void
%
% WRITES TO DISK:
%   saves path to horizon in job meta file
%   if constant horizon is given creates a horizon file in outputdir/user_horizon.

%%
    

  job_meta = load(job_meta_path);   % load job meta file
  flag_wb= str2double(varargin{1}); % flag says whether wb grid present or to create a contant grid
  horizon=varargin{2};
  
  if isfield(job_meta,horizon)
      fprintf( 'Horizon already exists in job meta file or invalid horizon name\n');
      return;
  
  end
  % if water bottom 3D grid file is given
  if flag_wb==1
      job_meta.(varargin{2}) = varargin{3};
      fprintf( strcat('Horizon ',horizon,' added to meta file\n'));
  
  elseif flag_wb ==0
      horz_value= str2double(varargin{3});
      il=job_meta.pkey_min:job_meta.pkey_inc:job_meta.pkey_max;
      xl=job_meta.skey_min:job_meta.skey_inc:job_meta.skey_max;
      length_horz=length(il)*length(xl);
      horz=zeros(length_horz,3);
      k=1;
      
      % create a conatnt horizon grid at all inline and xline locations
      for i = 1: length(il)
          for j =1:length(xl)
              horz(k,1)=il(i);
              horz(k,2)=xl(j);
              horz(k,3)=horz_value;
              k=k+1;
          end
      end
      
      
      % create directory if not already present
      if exist(strcat(job_meta.output_dir,'user_horizon/'),'dir') == 0
          mkdir(strcat(job_meta.output_dir,'user_horizon/'))
      end
      % same path to job meta file
      job_meta.(horizon) = strcat(job_meta.output_dir,'user_horizon/',horizon,'.txt');
      % save wb gridfile in the right place
      eval(strcat('save(job_meta.',horizon,',''horz'',''-ascii'',''-tabs'');'));
      fprintf( 'Constant horizon added to meta file');
  else
      fprintf( 'Needs more inputs');
  end
  
  
  save(job_meta_path,'-struct','job_meta','-v7.3');  
  
end

