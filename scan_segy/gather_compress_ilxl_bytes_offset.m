function compress_ilxl_bytes = gather_compress_ilxl_bytes_offset(trace_ilxl_bytes,blocktr)
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
% gather_compress_ilxl_bytes: function to compress scanned SEGY trace
% headers by finding repeating patterns. Only for pre-stack datasets that 
% are not angle gathers.
%   Arguments:
%       trace_ilxl_bytes = scan of trace headers created by segy_make_structure
%       blocktr = number of input traces  
%   
%   Outputs:
%       compress_ilxl_bytes = compressed version of trace_ilxl_bytes
%
%   Writes to Disk:
%       nothing
pkey_loc = 1; % column numbers needs to be implemented
skey_loc = 2;
byte_loc = 3;
tkey_loc = 4;

trace_ilxl_bytes(1:end-1,5) = diff(trace_ilxl_bytes(:,4));

start_idx = 1;
count = 0;
tcount = 1;
row_i = 1;
pkey_prev = -995837;
skey_prev = -9999437;

if blocktr > 1
    for row_i = start_idx:blocktr
        pkey = trace_ilxl_bytes(row_i,pkey_loc);
        skey = trace_ilxl_bytes(row_i,skey_loc);
        tkey = trace_ilxl_bytes(row_i,tkey_loc);
        tbyte = trace_ilxl_bytes(row_i,byte_loc);
        tkey_inc = trace_ilxl_bytes(row_i,5);
        
        if pkey == pkey_prev % same inline
            if skey == skey_prev
                tcount = tcount + 1;
                compress_offset_bytes(count,5:6) = [tcount 1];
                %count = count + 1;
                %tkey_inc_prev = tkey_inc;
                skey_prev = skey;
            else
                tcount = 1;
                count = count + 1;
                compress_offset_bytes(count,:) = [ pkey skey tbyte tcount tcount 1];
                %count = count + 1;
                %tkey_inc_prev = tkey_inc;
                skey_prev = skey;
            end
            
%             if tkey_inc == tkey_inc_prev
%                 tkey_inc_prev = tkey_inc;
%                 skey_prev = skey;
%             elseif skey == skey_prev
%                 compress_offset_bytes(count,5:6) = [tkey tkey_inc_prev];
%                 count = count + 1;
%                 tkey_inc_prev = tkey_inc;
%                 skey_prev = skey;
%             else
%                 compress_offset_bytes(count,:) = [ pkey skey tbyte tkey tkey 1];
%                 tkey_inc_prev = tkey_inc;
%                 skey_prev = skey;
%             end
            
        else % pkey ~= pkey_prev
            count = count + 1;
            tcount = 1;
            compress_offset_bytes(count,:) = [ pkey skey tbyte tcount tcount 1];
            pkey_prev = pkey;
            skey_prev = skey;
            %tkey_inc_prev = tkey_inc;
        end
        
    end
    
    % Now compress inlines, crosslines with same offset range
    blocktr = size(compress_offset_bytes,1);
    
    if blocktr > 1
        start_idx = 1;
        count = 0;
        row_i = 1;
        
        %blocktr = size(trace_ilxl_bytes,1);
        pkey_prev = -995837;
        skey_prev = -9999437;
        skey_inc = -999971;
        cur_inc = -27389995;
        tkey_min_prev = -999971;
        tkey_max_prev = 999971;
        tkey_inc_prev = 999971;
        tkey_inc_prev = -27389995;
        
        %compress_ilxl_bytes = zeros([blocktr,5],'int64');
        for row_i = start_idx:blocktr
            pkey = compress_offset_bytes(row_i,pkey_loc);
            skey = compress_offset_bytes(row_i,skey_loc);
            tbyte = compress_offset_bytes(row_i,byte_loc);
            tkey_min = compress_offset_bytes(row_i,4);
            tkey_max = compress_offset_bytes(row_i,5);
            tkey_inc = compress_offset_bytes(row_i,6);
            
            if pkey == pkey_prev && tkey_min == tkey_min_prev && tkey_max == tkey_max_prev && tkey_inc == tkey_inc_prev
                cur_inc = skey - skey_prev;
                if cur_inc ~= skey_inc
                    if cur_inc == 0
                        count = count + 1;
                        compress_ilxl_bytes(count,:) = [ pkey skey tbyte skey 1 tkey_min tkey_max tkey_inc];
                        skey_inc = -999971;
                        skey_prev = skey;
                        tkey_min_prev = tkey_min;
                        tkey_max_prev = tkey_max;
                        tkey_inc_prev = tkey_inc;
                    else
                        if skey_inc == -999971 % cur_inc is first time or after duplicate trace
                            compress_ilxl_bytes(count,4:5) = [ skey cur_inc ];
                            skey_inc = cur_inc;
                            skey_prev = skey;
                            tkey_min_prev = tkey_min;
                            tkey_max_prev = tkey_max;
                            tkey_inc_prev = tkey_inc;
                        else % cur_inc is not 0
                            count = count + 1;
                            %compress_ilxl_bytes(count,:) = [ pkey skey tbyte skey cur_inc tkey_min tkey_max tkey_inc];
                            compress_ilxl_bytes(count,:) = [ pkey skey tbyte skey 1 tkey_min tkey_max tkey_inc];
                            skey_inc = cur_inc;
                            skey_prev = skey;
                            tkey_min_prev = tkey_min;
                            tkey_max_prev = tkey_max;
                            tkey_inc_prev = tkey_inc;
                        end
                    end
                else % cur_inc == skey_inc
                    %compress_ilxl_bytes(count,4) = skey;
                    compress_ilxl_bytes(count,4:5) = [ skey cur_inc ];
                    skey_prev = skey;
                    tkey_min_prev = tkey_min;
                    tkey_max_prev = tkey_max;
                    tkey_inc_prev = tkey_inc;
                end
                
                
            else % pkey ~= pkey_prev
                count = count + 1;
                compress_ilxl_bytes(count,:) = [ pkey skey tbyte skey 1 tkey_min tkey_max tkey_inc];
                pkey_prev = pkey;
                skey_prev = skey;
                tkey_min_prev = tkey_min;
                tkey_max_prev = tkey_max;
                tkey_inc_prev = tkey_inc;
                skey_inc = -999971;
            end
            
        end
    else
        count = count + 1;
        compress_ilxl_bytes(count,:) = [ compress_offset_bytes(row_i,pkey_loc) compress_offset_bytes(row_i,skey_loc) compress_offset_bytes(row_i,byte_loc) compress_offset_bytes(row_i,skey_loc) 1 compress_offset_bytes(row_i,4) compress_offset_bytes(row_i,5) compress_offset_bytes(row_i,6)];
    end
    
else % for blocktr = 1
    count = count + 1;
    compress_ilxl_bytes(count,:) = [ trace_ilxl_bytes(row_i,pkey_loc) trace_ilxl_bytes(row_i,skey_loc) trace_ilxl_bytes(row_i,byte_loc) trace_ilxl_bytes(row_i,skey_loc) 1 trace_ilxl_bytes(row_i,4) trace_ilxl_bytes(row_i,4) 1];
end
    
end