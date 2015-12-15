function [traces,ilxl_read,offset_read] = segy_read_test()
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

%seismicpath = '/bgdata/msc_tg_2013/data/gathers_final_model_il1822_xl2100-3100_0-6000.segy';
%seismic.n_samples = int32(601);
%seismic.n_samples_str = '601*float32=>float32';

%seismicpath = '/data/CNS/segy/CGG_Cornerstone_TomoML_re-pro/PSDM_Final_Stack_Depth/Full_Offset_Stack_Depth';

seismicpath = '/segy/TZA/cj_segy_ibmfl.sgy';
seismic.n_samples = int32(1601);
seismic.n_samples_str = '1601*float32=>float32';

start_byte = 3600;
n_traces_to_read = int32(200);

seismic.file_type = 5;
seismic.pkey = 73;
seismic.skey = 77;
seismic.tkey = 37;
ll.bytes = double(single(seismic.n_samples)*4*single(n_traces_to_read));

    fid = fopen(char(seismicpath),'r','b');
    fseek(fid,start_byte,'bof');

    
%    if seismic.file_type == 1
        % Convert traces from IBM32FP read as UINT32 into IEEE64FP (doubles) - need to make it singles
        %traces_tmp = fread(fid,[60+seismic.n_samples,n_traces_to_read],strcat(num2str(seismic.n_samples),'*uint32=>uint32'));
        
        traces_tmp_orig2 = fread(fid,[60+seismic.n_samples,n_traces_to_read],'*uint32');
        traces_tmp = traces_tmp_orig2;
        tic
        
        
        %ilxl_read = traces_tmp(48:49,:)'; % what happens if the inline and crossline are not in this location  
        trchead = traces_tmp(1:60,:);
        [trace_header bytes_to_samples] = interpretbe(reshape(typecast(trchead(:),'uint16'),120,[]));
        
        % Inline / crossline compression step
        ilxl_read(:,1) = int32(trace_header(bytes_to_samples == seismic.pkey,:))';
        ilxl_read(:,2) = int32(trace_header(bytes_to_samples == seismic.skey,:))';
        offset_read = int32(trace_header(bytes_to_samples == seismic.tkey,:))';
        
        %offset_read = traces_tmp(10,:)';
        
        traces = single((1-2*double(bitget(traces_tmp(61:end,:),32))).*16.^ ...
        (double(bitshift(bitand(traces_tmp(61:end,:),uint32(hex2dec('7f000000'))),-24))-64).* ...
        (double(bitand(traces_tmp(61:end,:),uint32(hex2dec('00ffffff'))))/2^24));
    
    
    toc
fprintf('type 1 %d MB of data read (and some written) at %d MB/sec in %-10.5f seconds\n',round((ll.bytes-3600)/(1024*1024)),round(((ll.bytes-3600)/(1024*1024))/toc),toc);

    
%    elseif seismic.file_type == 2 
        %disp('This seismic file type is not currently supported. Please speak to Charles Jones.');

        
%         
%   % grab sign bit 32
% sgn = double(bitget(in,32));
% 
% % grab characteristic (exponent) bits 25-31
% %
% % Bitand isolates bits 25-31, dividing by 2^24 shifts them down to the
% % correct position
% 
% expo = double(bitand(in,2130706432))/2^24;
% 
% % grab fractional (mantissa) bits 1-24 with bitand 
% frac=double(bitand(in,16777215));  
% 
% % remove bias from exponent
% expo = expo - 64;
% 
% %normalize fractional
% frac=frac/2^24;
% 
% % put together for output as a float
% out = (1-2*sgn) .* 16.^expo .* frac;
% 
% 
% <<(x, n) % Left bit shift operator.
% >>(x, n) % Right bit shift operator, preserving the sign of x.
% >>>(x, n) % Unsigned right bit shift operator
% |(x, y) %  Bitwise or
% 
%     Unsigned right bit shift operator.
% reinterpret(Float32, uint32(fr >>> 9 | exp << 23 | sgn))

% 
% %special case: NaNs, Infinities and Zeros 
% temp = (expo==63) & (frac~=0);
% if any(temp) && nanflag
%     out(temp)= NaN;
%     warning('EGLTools:ibm2f:InputNaN','IBM2F: NaNs detected in IBM input; to override run IBM2F with nanflag=0');
% end
% 
% %special case: Inf and -Inf
% pos_inf_index = (expo==63) & (frac==0) & (sgn==0);
% out(pos_inf_index) = Inf;
% 
% neg_inf_index = (expo==63) & (frac==0) & (sgn==1);
% out(neg_inf_index) = -Inf;
% 
% %special case: zero
% zero_index = (expo==-64) & (frac==0);
% out(zero_index) = 0;
% 
% %error checks
% err = (frac==0) & (expo~=-64 | sgn~=0);   % Invalid zero input formats
% if any(err)
%    warning('EGLTools:ibm2f:IvalidZeroFormat','IBM2F: format for ZERO value invalid; is input 32-bit IBM FP?')
% end
% 
% err = (frac~=0) & (frac<1/16 | frac>=1);
% if any(err)
%    warning('EGLTools:ibm2f:InvalidFracValue','WARNING: IBM2F: fractional outside valid range; is input 32-bit IBM FP?')
% end   
%       
%         
%   ieeeOfPieces(fr, exp, sgn2) = reinterpret(Float32, uint32(fr >>> 9 | exp << 23 | sgn))
%   local fr::Uint32 = ntoh(reinterpret(Uint32, ibm))
%   local sgn::Uint32 = fr & 0x8000000
%   fr <<= 1 # shift sign out
%   local exp::Int32 = int32(fr >>> 25) # save exponent
%   fr <<= 7 # shift exponent out
%   
%   
%   
%   if fr == 0 # short-circuit for zero
%     ieeeOfPieces(0, 0, sgn)
%   else
%     # adjust exponent from base 16 offset 64 radix point before first digit to base 2 offset 127 radix point after first digit
%     # (exp - 64) * 4 + 127 - 1
%     exp = (exp - 64) * 4 + 127 - 1
%     while (fr < 0x80000000) != 0 # (re)normalize, 3 times max for normalized input
%       exp -= 1
%       fr <<= 1
%     end
%  
%     if exp <= 0
%       # complete underflow, return properly signed zero
%       # OR partial underflow, return denormalized number
%       fr = exp < -24 ? zero(Uint32) : (fr >>> -exp)
%       ieeeOfPieces(fr, 0, sgn)
%     elseif exp >= 255 then # overflow - return infinity
%       ieeeOfPieces(0, 255, sgn)
%     else # just a plain old number - remove the assumed high bit
%       ieeeOfPieces(fr << 1, exp, sgn)
%     end
%   end 
% end    

% <<(x, n) % Left bit shift operator.
% >>(x, n) % Right bit shift operator, preserving the sign of x.
% >>>(x, n) % Unsigned right bit shift operator
% |(x, y) %  Bitwise or
% &(x, y) %  Bitwise and

%traces(400:400,1)
%dec2bin(typecast(traces(400:400,1),'uint32'))
%   local fr::Uint32 = ntoh(reinterpret(Uint32, ibm))
%fseek(fid,start_byte,'bof');
%traces_tmp = fread(fid,[60+seismic.n_samples,n_traces_to_read],'*uint32');
traces_tmp = traces_tmp_orig2;
tic

%   local sgn::Uint32 = fr & 0x8000000
%sgn = uint32(bitget(traces_tmp(61:end,:),32));
%dec2bin(traces_tmp(460:460,1))
sgn = bitand(traces_tmp(61:end,:),hex2dec('80000000'),'uint32');
%dec2bin(sgn(400:400,1))

%origtraces_tmp = traces_tmp(61:end,:);
%   fr <<= 1 # shift sign out
traces_tmp(61:end,:) = bitshift(traces_tmp(61:end,:),1,'uint32');
%dec2bin(traces_tmp(460:460,1))

%   local exp::Int32 = int32(fr >>> 25) # save exponent
expo = int32(bitshift(traces_tmp(61:end,:),-25,'uint32'));

%dec2bin(expo(400:400,1))

expo = bitshift(expo,2,'int32') - 130;
%expo = (expo - 64) .* 4 + 127 - 1;
%dec2bin(expo(400:400,1))

%   fr <<= 7 # shift exponent out
%traces_tmp(61:end,:) = bitshift(traces_tmp(61:end,:),7,'uint32');
frac = bitshift(traces_tmp(61:end,:),7,'uint32');
%dec2bin(frac(400:400,1))

%reinterpret(Float32, uint32(fr >>> 9 | exp << 23 | sgn))
%traces = typecast(bitor(bitshift(traces_tmp(61:end,:),-9,'uint32'),bitshift(expo,23,'uint32'),sgn,'uint32'),'single');

% this needs some etra reshares to make vectors
%traces = typecast(bitor(bitor(bitshift(traces_tmp(61:end,:),-9,'uint32'),bitshift(expo,23,'uint32'),'uint32'),sgn,'uint32'),'single'); 

%       while (fr < 0x80000000) != 0 # (re)normalize, 3 times max for normalized input
%while
logfrac = frac < hex2dec('80000000');
%       exp -= 1
expo(logfrac) = expo(logfrac) - 1;
%       fr <<= 1
frac(logfrac) = bitshift(frac(logfrac),1,'uint32');
%     end
%dec2bin(frac(400:400,1))

%tracesout = arrayfun(exptest,expo,frac,sgn,'UniformOutput',true);
%tracesout = zeros(seismic.n_samples,n_traces_to_read,'single');
tracesout = reshape(typecast(bitor(bitor(bitshift(bitshift(frac(:),1,'uint32'),-9,'uint32'),typecast(bitshift(expo(:),23,'int32'),'uint32'),'uint32'),sgn(:),'uint32'),'single'),seismic.n_samples,n_traces_to_read);

% for t = 1:n_traces_to_read
%     for n = 1:seismic.n_samples
%         %tracesout(n, t) = exptest(expo(n,t),frac(n,t),sgn(n,t));
%         
%         if expo(n,t) <= 0
%             %       # complete underflow, return properly signed zero
%             %       # OR partial underflow, return denormalized number
%             %       fr = exp < -24 ? zero(Uint32) : (fr >>> -exp)
%             %       ieeeOfPieces(fr, 0, sgn)
%             if expo(n,t) < -24
%                 frac(n,t) = uint32(0);
%             else
%                 frac(n,t) = bitshift(frac(n,t),(-1 * expo(n,t)),'uint32');
%             end
%             out = typecast(bitor(bitor(bitshift(frac(n,t),-9,'uint32'),0),sgn(n,t),'uint32'),'single');
%             
%         elseif expo(n,t) >= 255 % overflow - return infinity
%             %       ieeeOfPieces(0, 255, sgn)
%             out = typecast(bitor(bitor(0,bitshift(255,23,'uint32'),'uint32'),sgn(n,t),'uint32'),'single');
%             
%         else % just a plain old number - remove the assumed high bit
%             %       ieeeOfPieces(fr << 1, exp, sgn)
%             %out = typecast(bitor(bitor(bitshift(bitshift(fr,1,'uint32'),-9,'uint32'),typecast(bitshift(exp,23,'int32'),'uint32'),'uint32'),sgnb,'uint32'),'single');
%             out = bitor(bitshift(bitshift(frac(n,t),1,'uint32'),-9,'uint32'),typecast(bitshift(expo(n,t),23,'int32'),'uint32'));
%             out = typecast(bitor(out,sgn(n,t),'uint32'),'single');
%         end
%         tracesout(n, t) = out;
%     end
% end

%tracesout = exptest(expo,frac,frac);
%     if exp <= 0
%       # complete underflow, return properly signed zero
%       # OR partial underflow, return denormalized number
%       fr = exp < -24 ? zero(Uint32) : (fr >>> -exp)
%       ieeeOfPieces(fr, 0, sgn)
% logexpo = expo <= 0;
% bsxfun
% frac(logexpo) =  

%     elseif exp >= 255 then # overflow - return infinity
%       ieeeOfPieces(0, 255, sgn)
%     else # just a plain old number - remove the assumed high bit
%       ieeeOfPieces(fr << 1, exp, sgn)



% frac = bitshift(frac,-9,'uint32');
% dec2bin(frac(400:400,1))
% 
% dec2bin(expo(400:400,1))
% expo = bitshift(expo,23,'int32');
% dec2bin(frac(400:400,1))
% dec2bin(expo(400:400,1))
% dec2bin(sgn(400:400,1))
% dec2bin(typecast(traces(400:400,1),'uint32'))
% 
% 
% traces = typecast(bitor(bitor(bitshift(frac(:),-8,'uint32'),typecast(bitshift(expo(:),23,'int32'),'uint32'),'uint32'),sgn(:),'uint32'),'single');
% 
% %     end
% %   end 
% 
% 
% 
% %ieeeOfPieces(fr, exp, sgn2) = reinterpret(Float32, uint32(fr >>> 9 | exp << 23 | sgn))
% traces = typecast(bitor(bitor(bitshift(frac(:),-9,'uint32'),bitshift(expo(:),23,'int32'),'uint32'),sgn(:),'uint32'),'single'); 

    toc
fprintf('type 2 %d MB of data read (and some written) at %d MB/sec in %-10.5f seconds\n',round((ll.bytes-3600)/(1024*1024)),round(((ll.bytes-3600)/(1024*1024))/toc),toc);

%   ieeeOfPieces(fr, exp, sgn2) = reinterpret(Float32, uint32(fr >>> 9 | exp << 23 | sgn))
%   function convert(::Type{Float32}, ibm::IbmFloat32)
%   local fr = ntoh(reinterpret(Uint32, ibm))
%   local sgn::Uint32 = fr & 0x8000000 # save sign
%   fr <<= 1 # shift sign out
%   local exp::Int32 = int32(fr >>> 25) # save exponent
%   fr <<= 7 # shift exponent ou
%   if fr == 0 # short-circuit for zero
%     ieeeOfPieces(0, 0, sgn)
%   else
%     # adjust exponent from base 16 offset 64 radix point before first digit to base 2 offset 127 radix point after first digit
%     # (exp - 64) * 4 + 127 - 1
%     exp = (exp - 64) * 4 + 127 - 1
%     while (fr < 0x80000000) != 0 # (re)normalize, 3 times max for normalized input
%       exp -= 1
%       fr <<= 1
%     end
%  
%     if exp <= 0
%       # complete underflow, return properly signed zero
%       # OR partial underflow, return denormalized number
%       fr = exp < -24 ? zero(Uint32) : (fr >>> -exp)
%       ieeeOfPieces(fr, 0, sgn)
%     elseif exp >= 255 then # overflow - return infinity
%       ieeeOfPieces(0, 255, sgn)
%     else # just a plain old number - remove the assumed high bit
%       ieeeOfPieces(fr << 1, exp, sgn)
%     end
%   end 
% end    



% Only integer shifting and masking are used.
% *************************************************************************
% Credits: CWP: Brian Sumner,  c.1985
% *************************************************************************/
% {
%     register int fconv, fmant, i, t;
% 
%     for (i = 0;i < n; ++i) {
% 
% 	fconv = from[i];
% 
% 	/* if little endian, i.e. endian=0 do this */
% 	if (endian == 0) fconv = (fconv << 24) | ((fconv >> 24) & 0xff) |
% 		((fconv & 0xff00) << 8) | ((fconv & 0xff0000) >> 8);
% 
% 	if (fconv) {
% 	    fmant = 0x00ffffff & fconv;
% 	    if (fmant == 0)
% 		fconv = 0;
% 	    else {
% 	        t = (int) ((0x7f000000 & fconv) >> 22) - 130;
% 	        while (!(fmant & 0x00800000)) { --t; fmant <<= 1; }
% 	        if (t > 254) fconv = (0x80000000 & fconv) | 0x7f7fffff;
% 	        else if (t <= 0) fconv = 0;
% 	        else fconv =   (0x80000000 & fconv) | (t << 23)
% 			 | (0x007fffff & fmant);
% 	    }
% 	}
% 	to[i] = fconv;
%     }
%     return;
% }


%    elseif seismic.file_type == 5
        % Traces are IEEE32FP (singles)   
        %traces = fread(fid,[60+seismic.n_samples,n_traces_to_read],strcat(num2str(seismic.n_samples),'*float32=>float32'));
        %cjtypestr = strcat(seismic.n_samples_str,'*float32=>float32');
        fseek(fid,start_byte,'bof');
        traces = fread(fid,[60+seismic.n_samples,n_traces_to_read],seismic.n_samples_str);
        tic
        %trace_headers = typecast(single(reshape(traces(1:60,:),1,60*n_traces_to_read)),'int32');  
        %trace_headers = reshape(trace_headers,60,n_traces_to_read);
        
        trchead = typecast(single(reshape(traces(1:60,:),1,60*n_traces_to_read)),'uint16');  
        %trchead = reshape(trchead,120,n_traces_to_read);
        
        [trace_header bytes_to_samples] = interpretbe(reshape(trchead,120,[]));
        
        % Inline / crossline compression step
        ilxl_read(:,1) = int32(trace_header(bytes_to_samples == seismic.pkey,:))';
        ilxl_read(:,2) = int32(trace_header(bytes_to_samples == seismic.skey,:))';
        offset_read = int32(trace_header(bytes_to_samples == seismic.tkey,:))';
        
        %ilxl_read = trace_headers(48:49,:)';        
        %offset_read = trace_headers(10,:)';        
        traces = traces(61:end,:);      
            toc
fprintf('type 3 %d MB of data read (and some written) at %d MB/sec in %-10.5f seconds\n',round((ll.bytes-3600)/(1024*1024)),round(((ll.bytes-3600)/(1024*1024))/toc),toc);

%     else
%         disp('This seismic file type is not currently supported. Please speak to Charles Jones.');
%     end

    fclose(fid);  
end

function [out] = exptest(exp,fr,sgnb)
if exp <= 0
    %       # complete underflow, return properly signed zero
    %       # OR partial underflow, return denormalized number
    %       fr = exp < -24 ? zero(Uint32) : (fr >>> -exp)
    %       ieeeOfPieces(fr, 0, sgn)
    if exp < -24
        fr = uint32(0);
    else
        fr = bitshift(fr,(-1 * exp),'uint32');
    end
    out = typecast(bitor(bitor(bitshift(fr,-9,'uint32'),0),sgnb,'uint32'),'single');
    
elseif exp >= 255 % overflow - return infinity
    %       ieeeOfPieces(0, 255, sgn)
    out = typecast(bitor(bitor(0,bitshift(255,23,'uint32'),'uint32'),sgnb,'uint32'),'single');
    
else % just a plain old number - remove the assumed high bit
    %       ieeeOfPieces(fr << 1, exp, sgn)
    %out = typecast(bitor(bitor(bitshift(bitshift(fr,1,'uint32'),-9,'uint32'),typecast(bitshift(exp,23,'int32'),'uint32'),'uint32'),sgnb,'uint32'),'single');
    out = bitor(bitshift(bitshift(fr,1,'uint32'),-9,'uint32'),typecast(bitshift(exp,23,'int32'),'uint32'));
    out = typecast(bitor(out,sgnb,'uint32'),'single');
end
end

function [trace_header bytes_to_samples] = interpretbe(tmptrheader)
byte_type = [ ...
    2*ones(7,1); ones(4,1);
    2*ones(8,1); ones(2,1);
    2*ones(4,1); ones(46,1);
    2*ones(5,1); ones(2,1);
    2*ones(1,1); ones(5,1);
    2*ones(1,1); ones(1,1);
    2*ones(1,1); ones(2,1);
    2*ones(1,1); 2*ones(1,1)];

ntr = size(tmptrheader,2);
trace_header = zeros(91,ntr);
bytes_to_samples = zeros(91,1);

count =1;
for ii = 1:91
    bytes_to_samples(ii,1) = 2*count-1;
    if byte_type(ii) == 1
        trace_header(ii,:) = double(tmptrheader(count,:));
        count = count+1;
    elseif byte_type(ii) == 2
        trace_header(ii,:) = double(tmptrheader(count+1,:))*2^16 + double(tmptrheader(count,:)); % note this is big-endian and different to one in segy_make_structure
        count = count+2;
    end
end

trace_header(21,:) = trace_header(21,:)-2^16;

end
