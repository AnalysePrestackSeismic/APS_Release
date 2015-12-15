function [filter,output] = match_wavelets(wref,wout,NF,mu) 
  % wref wout
%mu:  prewhitening as percentage of the zero lag autocorrelation

 NW = length(wref);              % length of the input
 
 if floor(NF/2) == NF/2;
     NF = NF + 1;
 end
  
 % Compute auto-correlation of reference wavelet
 % Cref = convmtx(wref,NF);
 % autoref3 =  Cref'*Cref;   
 autoref = xcorr(wref,wref,floor(NF/2));
 autoref = spdiags(repmat(autoref,1,NF),[-floor(NF/2):1:floor(NF/2)]);
  
 % add pre-whitening
 autoref = autoref+((autoref(1,1)*mu/100)*eye(NF));
 
 % Compute cross-correlation of reference wavelet with desired wavelet
 crossref = xcorr(wref,wout,floor(NF/2));
 %pad = NW+NF-1-length(wout);
 %crossref = Cref'*[wout; zeros(pad,1)];
 
 filter = autoref\crossref;   
 
 output = conv(wref,filter,'same'); 
 
 figure
 subplot(2,2,1); plot(wref)
 subplot(2,2,2); plot(wout)
 subplot(2,2,3); plot(filter)
 subplot(2,2,4); plot(wref); hold all; plot(wout); plot(output);
 hold off
                              



