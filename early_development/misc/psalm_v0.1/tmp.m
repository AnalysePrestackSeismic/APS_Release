clear all

nr = 100;
nc = 200;

winr = 10;
winc = 20;

A = zeros(nr+winr,nc+winc);
B = zeros(nr+winr,nc+winc);
C = zeros(nr+winr,nc+winc);
M = zeros(nr,nc);

tmp = rand(nr,nc)-0.5;

A(winr/2:nr-1+winr/2,winc/2:nc-1+winc/2) = tmp;

tmp = 2*tmp+0.25*rand(nr,nc)-0.125;

B(winr/2:nr-1+winr/2,winc/2:nc-1+winc/2) = tmp;

tic

tmp = reshape((1:1:nr*nc)',nr,nc);

C(winr/2:nr-1+winr/2,winc/2:nc-1+winc/2) = tmp;

for ii=1:winr
    for jj=1:winc
        D{ii,jj} = reshape(cell2mat(reshape(mat2cell(A(ii:nr+ii-1,jj:nc+jj-1),winr*ones(1,nr/winr),winc*ones(1,nc/winc)),1,[])),winr*winc,[]);
        E{ii,jj} = reshape(cell2mat(reshape(mat2cell(B(ii:nr+ii-1,jj:nc+jj-1),winr*ones(1,nr/winr),winc*ones(1,nc/winc)),1,[])),winr*winc,[]);
        F{ii,jj} = reshape(cell2mat(reshape(mat2cell(C(ii:nr+ii-1,jj:nc+jj-1),winr*ones(1,nr/winr),winc*ones(1,nc/winc)),1,[])),winr*winc,[]);
    end
end

G = cell2mat(reshape(D,1,[]));
H = cell2mat(reshape(E,1,[]));
I = cell2mat(reshape(F,1,[]));

I(I==0) = nr*nc+1;

[~, J] = min(I);
K = I(J(1,1),:);

L = sum(bsxfun(@minus,G,mean(G)).*bsxfun(@minus,H,mean(H)))./sum(bsxfun(@minus,G,mean(G)).^2);

M(K) = L;

toc