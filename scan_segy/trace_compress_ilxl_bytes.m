function compress_ilxl_bytes  = trace_compress_ilxl_bytes(trace_ilxl_bytes,blocktr)
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
% trace_compress_ilxl_bytes: function to compress scanned SEGY trace
% headers by finding repeating patterns. Only for post-stack datasets.
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

start_idx = 1;
count = 0;
row_i = 1;

%blocktr = size(trace_ilxl_bytes,1);
pkey_prev = -995837;
skey_prev = -9999437;
skey_inc = -999971;
cur_inc = -27389995;

%compress_ilxl_bytes = zeros([blocktr,5],'int64');

if blocktr > 1
    for row_i = start_idx:blocktr
        pkey = trace_ilxl_bytes(row_i,pkey_loc);
        skey = trace_ilxl_bytes(row_i,skey_loc);
        tbyte = trace_ilxl_bytes(row_i,byte_loc);
        
        if pkey == pkey_prev
            cur_inc = skey - skey_prev;
            if cur_inc ~= skey_inc
                if cur_inc == 0
                    count = count + 1;
                    compress_ilxl_bytes(count,:) = [ pkey skey tbyte skey 1 ];
                    skey_inc = -999971;
                    skey_prev = skey;
                else
                    if skey_inc == -999971 % cur_inc is first time or after duplicate trace
                        compress_ilxl_bytes(count,4:5) = [ skey cur_inc ];
                        skey_inc = cur_inc;
                        skey_prev = skey;
                    else % cur_inc is not 0
                        count = count + 1;
                        compress_ilxl_bytes(count,:) = [ pkey skey tbyte skey cur_inc ];
                        skey_inc = cur_inc;
                        skey_prev = skey;
                    end
                end
            else % cur_inc == skey_inc
                compress_ilxl_bytes(count,4) = skey;
                skey_prev = skey;
            end
            
            
        else % pkey ~= pkey_prev
            count = count + 1;
            compress_ilxl_bytes(count,:) = [ pkey skey tbyte skey 1 ];
            pkey_prev = pkey;
            skey_prev = skey;
            skey_inc = -999971;
        end
        
    end
    
else % for blocktr = 1
    count = count + 1;
    compress_ilxl_bytes(count,:) = [ trace_ilxl_bytes(row_i,pkey_loc) trace_ilxl_bytes(row_i,skey_loc) trace_ilxl_bytes(row_i,byte_loc) trace_ilxl_bytes(row_i,skey_loc) 1 ];
end
% now truncate the array
%compress_ilxl_bytes(count+1:end,:) = [];

end
