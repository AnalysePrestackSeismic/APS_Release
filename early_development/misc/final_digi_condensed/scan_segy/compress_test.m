function [compress_ilxl_bytes] = compress_test()
%COMPRESS_TEST Summary of this function goes here
%   Detailed explanation goes here


trace_ilxl_bytes_orig = [0 3 80; 0 3 81; 0 3 82; 1 7 83; 1 8 84; 1 9 85; 1 9 86; 2 4 87; 2 6 88; 2 7 89; 2 8 90; 2 8 91; 2 9 92; 3 1 93; 3 2 94; 3 3 95; 3 6 96; 3 8 97; 3 9 98; 4 4 99 ];
trace_ilxl_bytes = repmat(trace_ilxl_bytes_orig,200,1);
    pkey_loc = 1; % column numbers needs to be implemented
    skey_loc = 2;
    byte_loc = 3;
    skey_max_loc = 4;
    skey_inc_loc = 5;
    tkey_loc = 6;
    tkey_max_loc = 7;
    tkey_inc_loc = 8;

    start_idx = 1;
    count = 1;
    constant = 10000;
    testc = 1;
    

    
    if testc == 0;
        
        row_i = 0;
        blocktr = size(trace_ilxl_bytes,1);
        compress_ilxl_bytes = zeros(blocktr,5,'int64');
        
        while start_idx < blocktr
            cdp_1 = (trace_ilxl_bytes(start_idx,1)+constant)*(trace_ilxl_bytes(start_idx,2));
            if start_idx < blocktr-1
                cdp_2 = (trace_ilxl_bytes(start_idx+1,1)+constant)*(trace_ilxl_bytes(start_idx+1,2));
            end
            if start_idx < blocktr-2
                cdp_3 = (trace_ilxl_bytes(start_idx+2,1)+constant)*(trace_ilxl_bytes(start_idx+2,2));
            end
            if (cdp_2 - cdp_1) == (cdp_3 - cdp_2)
                row_i = row_i+1;
            elseif trace_ilxl_bytes(start_idx,1) == trace_ilxl_bytes(start_idx+1,1)
                
                compress_ilxl_bytes(count,1) = trace_ilxl_bytes(start_idx-row_i,1);
                compress_ilxl_bytes(count,2) = trace_ilxl_bytes(start_idx-row_i,2);
                compress_ilxl_bytes(count,3) = trace_ilxl_bytes(start_idx-row_i,3);
                compress_ilxl_bytes(count,4) = trace_ilxl_bytes(start_idx+1,2);
                
                xl_inc = trace_ilxl_bytes(start_idx-row_i+1,2)-trace_ilxl_bytes(start_idx-row_i,2);
                compress_ilxl_bytes(count,5) = xl_inc;
                
                row_i = 0;
                count = count + 1;
            end
            start_idx = start_idx + 1;
        end
        compress_ilxl_bytes(end,4) = trace_ilxl_bytes(end,2);
        
    else
        
        start_idx = 1;
        count = 0;
        row_i = 1;
        blocktr = size(trace_ilxl_bytes,1);
        compress_ilxl_bytes = zeros([(blocktr*3),5],'int64');
        
        pkey_prev = -995837;
        skey_prev = -9999437;
        skey_inc = -999971;
        cur_inc = -27389995;
        
        if blocktr > 1
            for row_i = start_idx:blocktr
                pkey = trace_ilxl_bytes(row_i,1);
                skey = trace_ilxl_bytes(row_i,2);
                tbyte = trace_ilxl_bytes(row_i,3);
                
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
            compress_ilxl_bytes(count,:) = [ trace_ilxl_bytes(row_i,1) trace_ilxl_bytes(row_i,2) trace_ilxl_bytes(row_i,1) trace_ilxl_bytes(row_i,2) 0 ];
            
        end
        compress_ilxl_bytes(count+1:end,:) = [];
        
    end
end

