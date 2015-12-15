function [] = add_water_bottom_to_job_meta(job_meta_path,varargin)
%% Funtion to add water bottom to job meta file. Use when water bottom picking is likely to go wrong. e.g, very shallow water

% INPUTS:
%   job_meta_path = path of job meta file

%   if you have a 3D water bottom (for 3D surveys only:           varargin = '1','wb_path' file should be of the  format: IL XL Z
%   if you want a constant waterbottom:    varargin{1} ='0',varargin{2} = 'wbvalue in TWT or depth' 
%   if you have a 2D waterbottom for 2D line surveys: varargin ='2','CDP_no column no.','z column number'. 


% OUPUT:
%   void

% WRITES TO DISK:
%   saves wb_path in job meta file
%   if constant water bottom is given creates an wb file in outputdir/wb_use.
%   if its a 2D wb , it manipulates the file and writes another ascii file to
%   disk which has its adress in job meta file. Donot run this twice on a
%   2D line.

%%
    

  job_meta = load(job_meta_path);   % load job meta file
  flag_wb= str2double(varargin{1}); % flag says whether wb grid present or to create a contant grid
  
  % if water bottom 3D grid file is given
  if flag_wb==1
      wb_path = varargin{2};
      job_meta.wb_path = wb_path;
      fprintf( 'Water bottom added to meta file');
  
  % if water bottom 2D grid file is given
  elseif flag_wb==2
      
      wb_path = varargin{2};
      wb=dlmread(wb_path);% load water bottom file
      wb2=zeros(size(wb,1),3);
      wb2(:,2)=wb(:,str2double(varargin{3}));
      wb2(:,3)=wb(:,str2double(varargin{4}))*1000;
      wb2(:,1)=ones(size(wb,1),1)*job_meta.block_keys(1,1);
      
      if exist(strcat(job_meta.output_dir,'2D_wb_user/'),'dir') == 0
          mkdir(strcat(job_meta.output_dir,'2D_wb_user/'))
      end
      
      wb_path2=strcat(job_meta.output_dir,'2D_wb_user/2D_wb_user.txt');
      dlmwrite(wb_path2, wb2);
      job_meta.wb_path = wb_path2;
      % if a constant water bottom grid needs to be created
  elseif flag_wb ==0
      wb_value= str2double(varargin{2});
      il=job_meta.pkey_min:job_meta.pkey_inc:job_meta.pkey_max;
      xl=job_meta.skey_min:job_meta.skey_inc:job_meta.skey_max;
      length_wb=length(il)*length(xl);
      wb=zeros(length_wb,3);
      k=1;
      
      % create a conatnt wb grid at all inline and xline locations
      for i = 1: length(il)
          for j =1:length(xl)
              wb(k,1)=il(i);
              wb(k,2)=xl(j);
              wb(k,3)=wb_value;
              k=k+1;
          end
      end
      
      
      % create directory if not already present
      if exist(strcat(job_meta.output_dir,'wb_user/'),'dir') == 0
          mkdir(strcat(job_meta.output_dir,'wb_user/'))
      end
      % same path to job meta file
      job_meta.wb_path =strcat(job_meta.output_dir,'wb_user/wb_user.txt') ;
      % save wb gridfile in the right place
      save (job_meta.wb_path,'wb','-ascii', '-tabs');
      fprintf( 'Constant water bottom added to meta file');
  else
      fprintf( 'Needs more inputs');
  end
  
  
  save(job_meta_path,'-struct','job_meta','-v7.3');  
  
end

