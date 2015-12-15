% number of rows and columns
nr = 100;
nc = 200;

% window size in rows and columns
winr = 10;
winc = 20;

% pre-allocate memory
A = zeros(nr+winr,nc+winc);
B = zeros(nr+winr,nc+winc);
C = zeros(nr+winr,nc+winc);
D = cell(winr,winc);
E = cell(winr,winc);
F = cell(winr,winc);
M = zeros(nr,nc);

% make some random numbers
tmp = rand(nr,nc)-0.5;

% put the random numbers into A -> data1
A(winr/2:nr-1+winr/2,winc/2:nc-1+winc/2) = tmp;

% make some more random numbers, correlated to the first lot
tmp = 2*tmp+0.25*rand(nr,nc)-0.125;

% put the random numbers into B -> data2
B(winr/2:nr-1+winr/2,winc/2:nc-1+winc/2) = tmp;

% start timing
tic

% make the linear indices for our data
tmp = reshape((1:1:nr*nc)',nr,nc);

% put the linear indices in a matrix of the same size as A and B
C(winr/2:nr-1+winr/2,winc/2:nc-1+winc/2) = tmp;

% re-sort our data into columns; each column is the window of data we need
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

% do the vectorized calculation
L = sum(bsxfun(@minus,G,mean(G)).*bsxfun(@minus,H,mean(H)))./sum(bsxfun(@minus,G,mean(G)).^2);

% get the linear indices, which have been re-sorted in exactly the same way as the data
I(I==0) = nr*nc+1;
[~, J] = min(I);
K = I(J(1,1),:);

% put the result in the correct location
M(K) = L;

toc