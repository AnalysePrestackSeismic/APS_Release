function [ block_maxzout] = make_horizon_mask( job_meta_path,horizon_path,horizon_shift, flag_plot)
%------------------------------------------------------------------
% MAKE_HORIZON_MASK = Creates a mask of horizon decimated at block
% resolution

% Arguments:Pass all arguments a strings under single quotes
% job_meta_path: job_mat_path = path of metadata .mat file.
% horizon_path = path of the horizon (Horizon must be in IL, XL, ZValue
% Format)
% horizon_shift = Amount of shift in horizon in ms or m
% flag_plot = '1' if you want ot plot mask

%   Ouput: block_maxzout: Array of max_zout for live_blocks in number of
%   samples
%%   
%-------------LOAD AND INTIALIZE-----------------------------------------
job_meta = load(job_meta_path);                                         % Load Job Meta file
horizon = dlmread(horizon_path);                                        % Import ASCII Horizon
MAX_horizon = max(horizon(:,3));                                        % Deepest Value of horizon
block_maxzout = MAX_horizon * ones (size(job_meta.liveblocks,1),1);     % Initializing out put array indexed to live blocks to deepest horizon value
n = size(horizon,1);
horizon_shift = str2double (horizon_shift);
flag_plot = str2double (flag_plot);

if(horizon_shift~=0)
    horizon (:,3) = horizon(:,3)+horizon_shift;                         % Shift the horizon z values by user defined value
end

%%
%--------------LOOP TO FIND MAXZOUT FOR BLOCKS----------------------------
nliveblocks = size( job_meta.liveblocks,1);                     % Number of live blocks
dec = max(1,floor(n/40000));                                    % Decimation for search (minimum decimation is st to 1)
%Loop to find out shallowest z for each block
for row = 1:dec:n                                               % Scan though the horizon points
    %fprintf ('+,\n');
    for i_live_block = 1:nliveblocks                            % Loop though the blocks
        i_block = job_meta.liveblocks(i_live_block);            % Reference the block numbers for live blocks
        il = horizon (row,1);                                   % Inline Value for point
        xl = horizon (row,2);                                   % Inline Value for point
        if( (job_meta.block_keys(i_block,1)<il && (job_meta.block_keys(i_block,2)>il) && (job_meta.block_keys(i_block,3)<xl) && (job_meta.block_keys(i_block,4)>xl))) % check if the point is in this block
            %fprintf ('-,');
            z = horizon (row,3);                                % z value at the horizon point
            if(z < block_maxzout(i_live_block))
                block_maxzout(i_live_block) = z;                % update max zout of the block if the current point is shallower
            end                        
        end
    end
end
    
block_maxzout = floor((block_maxzout *1000)/job_meta.s_rate);   % Convert maxzout array into number of samples dividing by sampling rate
%%
% ---------------PLOT THE MASK--------------------------
if(flag_plot == 1)
    cjxdata = zeros(4,size(job_meta.block_keys,1));             % preallocation
    cjydata = zeros(4,size(job_meta.block_keys,1));             % preallocation
    cdata = NaN*ones (size(job_meta.block_keys,1),1);           % preallocation

    % Defining the blocks by the bounding rectangle
    for i_block = 1:1:size(job_meta.block_keys,1)
    cjxdata(:,i_block) = [job_meta.block_keys(i_block,1); job_meta.block_keys(i_block,2); job_meta.block_keys(i_block,2); job_meta.block_keys(i_block,1)];  % block inlines
    cjydata(:,i_block) = [job_meta.block_keys(i_block,3); job_meta.block_keys(i_block,3); job_meta.block_keys(i_block,4); job_meta.block_keys(i_block,4)];  % block crosslines
    end

    for lpi=1:nliveblocks
    i_block = job_meta.liveblocks(lpi);                         % Reference the block numbers for live blocks
    cdata(i_block,1) = floor(block_maxzout(lpi));               % Rreferencing colour data to z values
    end
    close                                                       % Close any previously open figure
    figure                                                      % Initiate a graphical object
    p = patch(cjxdata,cjydata,'b');                             % Plot the framework of blocks 
    colormap (jet);
    set(p,'EdgeColor','b','FaceColor','flat','FaceVertexCData',cdata,'CDataMapping','scaled'); % colour the blocks by the maxzout values
end
%---------------------------------------------------------
end

