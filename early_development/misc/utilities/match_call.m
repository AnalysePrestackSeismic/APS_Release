function wfilter = match_call(reference,shapee,prewhite)

% USE:      type: filter = match(reference,shapee,prewhite,'file_out');
%
% INPUTS:   reference = wavelet to match to
%           shapee = wavelet to shape
%           prewhite = prewhitenning as dB down from 0 power (should be a 
%                           negative number)
%           file_out = 'filename for output filter' (must be entered within
%                           quotation marks)
%
% OUTPUTS:  filter = time domain least-squares matching filter to match
%                       shapee to reference. Filter is scaled such that the
%                       RMS of shapee is the same after filtering as
%                       before.

warning off

% get number of samples of inputs
nt = 1000;
%nt=nt*4;

% prewhiten shapee
fa_reference = abs(fft(reference-mean(reference),nt));
fp_reference = 20*log10(fa_reference/max(fa_reference));
bandwidth = (fp_reference > prewhite);

fa_shapee = abs(fft(shapee,nt));
fp_shapee = 20*log10(fa_shapee/max(fa_shapee));
edges = logical(abs([0;diff(bandwidth)]));
edges_idx = cumsum(edges);
low = (edges_idx == 0);
mid = (edges_idx == 2);
high = (edges_idx == 4);
fp_shapee(low) = fp_shapee(and(high,edges));
fp_shapee(mid) = fp_shapee(and(mid,edges));
fp_shapee(high) = fp_shapee(and(high,edges));
faw_shapee = (10.^(fp_shapee/20))*max(fa_shapee);
nt = length(reference);
wshapee = circshift(ifft(faw_shapee,nt,'symmetric'),floor(nt/2));

shapee = circshift(ifft(fa_shapee,nt,'symmetric'),floor(nt/2));

white = zeros(nt,1);
white(floor((nt/2)+1),1) = max(wshapee)-max(shapee);
wshapee = shapee+white;

% make convolution matrix
% wG = convmtx(wshapee,nt);   % with pre-whitenning
% G = convmtx(shapee,nt);     % without pre-whitenning
clip = (((2*nt)-1)-nt)/2;
if round(clip)==clip    
%     wG = wG(clip+1:end-clip,:);
%     G = G(clip+1:end-clip,:);    
% code for when user has no signal processing toolbox (avoids convmtx)
    for i = 1:nt
        wG(i,:) = circshift(wshapee,i+clip);
        G(i,:) = circshift(shapee,i+clip);
    end
    A = spdiags(wG,(-clip:1:clip));
    wG = spdiags(A,(-clip:1:clip),nt,nt);
    B = spdiags(G,(-clip:1:clip));
    G = spdiags(B,(-clip:1:clip),nt,nt);
else 
    clip = round(clip);
%     wG = wG(clip+1:end-clip+1,:);   % with pre-whitenning
%     G = G(clip+1:end-clip+1,:);     % without pre-whitenning
% code for when user has no signal processing toolbox (avoids convmtx)    
    for i = 1:nt
        wG(i,:) = circshift(wshapee,i+clip-1);
        G(i,:) = circshift(shapee,i+clip-1);
    end
    A = spdiags(wG,(-clip:1:clip));
    wG = spdiags(A,(-clip:1:clip),nt,nt);
    B = spdiags(G,(-clip:1:clip));
    G = spdiags(B,(-clip:1:clip),nt,nt);
end

reference = double(reference);

% calculate least-squares matching filter
[wfilter flag] = lsqr(wG,reference,1e-4,1000);  % with pre-whitenning
[filter flag] = lsqr(G,reference,1e-4,1000);    % without pre-whitenning

% ensure RMS of shapee remains the same after filtering
wscalar = norm(shapee)/norm(G*wfilter);     % with pre-whitenning
wfilter = wfilter*wscalar;                  % with pre-whitenning
scalar = norm(shapee)/norm(G*filter);       % without pre-whitenning
filter = filter*scalar;                     % without pre-whitenning

end