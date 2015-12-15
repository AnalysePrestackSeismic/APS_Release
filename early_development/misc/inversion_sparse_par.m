%tic

matlabpool
clear all
load workspace.mat
warning off all

angles = [10;20;30]*(pi/180);
nangles = length(angles);

% Import the file
% rawData1 = importdata('n.txt');
% name = 'ndata';
% newData1.(genvarname(name)) = rawData1;
% vars = fieldnames(newData1);
% for i = 1:length(vars)
%     assignin('base', vars{i}, newData1.(vars{i}));
% end
% rawData2 = importdata('m.txt');
% name = 'mdata';
% newData2.(genvarname(name)) = rawData2;
% vars = fieldnames(newData2);
% for i = 1:length(vars)
%     assignin('base', vars{i}, newData2.(vars{i}));
% end
% rawData3 = importdata('f.txt');
% name = 'fdata';
% newData3.(genvarname(name)) = rawData3;
% vars = fieldnames(newData3);
% for i = 1:length(vars)
%     assignin('base', vars{i}, newData3.(vars{i}));
% end

% inline = ndata(:,1);
% xline = ndata(:,2);
% ndata = ndata(:,3:end)';
% mdata = mdata(:,3:end)';
% fdata = fdata(:,3:end)';

% [nt nrec] = size(ndata);

nt = 751;
nrec = 1687;
nil = 768;
block = 60;

zpmod = ones(nt+1,nrec*block)*6892.6;
zsmod = ones(nt+1,nrec*block)*3376.41;
rhomod = ones(nt+1,nrec*block)*2.41;

tstart = 0; 
tend = 3000;
tsample = 4;
sstart = (tstart/tsample)+1;
send = (tend/tsample)+1;
no_rho = 1;

k = 1.37602;
m = 0.18682;
kc = -4.02995;
mc = -0.796873;
gamma = 0.4964;

c1 = 1+(tan(angles).*tan(angles));
c2 = (-8)*(gamma^2)*(tan(angles).*tan(angles));
c3 = ((-0.5)*(tan(angles).*tan(angles)))+(2*(gamma^2)*(sin(angles).*sin(angles)));
c1 = (0.5.*c1)+(0.5*k.*c2)+(m.*c3);
c2 = 0.5.*c2;

w1 = nearw;
w2 = midw;
w3 = farw;

nw = length(w1);

Diff = sparse(1:nt,1:nt,-1*ones(1,nt),nt,nt+1) + sparse(1:nt,2:nt+1,1*ones(1,nt),nt,nt+1);

Rp = Diff*log(zsinv(sstart:send+1));
Rs = Diff*log(zpinv(sstart:send+1));
Rd = Diff*log(rhoinv(sstart:send+1));

Rn=(c1(1,1)*Rp)+(c2(1,1)*Rs)+(c3(1,1)*Rd);
Rm=(c1(2,1)*Rp)+(c2(2,1)*Rs)+(c3(2,1)*Rd);
Rf=(c1(3,1)*Rp)+(c2(3,1)*Rs)+(c3(3,1)*Rd);

% nscale = sum(abs(ndata(:,821)))/(sum(abs(Rn))*sum(abs(w1)));
% mscale = sum(abs(mdata(:,821)))/(sum(abs(Rm))*sum(abs(w2)));
% fscale = sum(abs(fdata(:,821)))/(sum(abs(Rd))*sum(abs(w3)));

nscale = 152/(sum(abs(Rn))*sum(abs(w1)));
mscale = 154/(sum(abs(Rm))*sum(abs(w2)));
fscale = 135/(sum(abs(Rd))*sum(abs(w3)));

fudge = 1;

w1 = w1*nscale*fudge;
w2 = w2*mscale*fudge;
w3 = w3*fscale*fudge;

W1 = convmtx(w1,nt);
W2 = convmtx(w2,nt);
W3 = convmtx(w3,nt);

[nrw ncw] = size(W1);
cropw = ceil((nrw-ncw)/2);

W1 = sparse(W1(cropw:nrw-cropw,:));
W2 = sparse(W2(cropw:nrw-cropw,:));
W3 = sparse(W3(cropw:nrw-cropw,:));

if no_rho == 1
    c3=c3*0;
end

%build operator
G = [(c1(1,1)*W1*Diff) (c2(1,1)*W1*Diff) (c3(1,1)*W1*Diff); (c1(2,1)*W2*Diff) (c2(2,1)*W2*Diff) (c3(2,1)*W2*Diff); (c1(3,1)*W3*Diff) (c2(3,1)*W3*Diff) (c3(3,1)*W3*Diff)];

itermax = 30;
mu = 0;
tol = 0.025;

Lpmod = log(zpmod);
dLsmod = log(zsmod)-(k*Lpmod)-kc;
dLdmod = log(rhomod)-(m*Lpmod)-mc;

Model = [Lpmod(sstart:send+1,:); dLsmod(sstart:send+1,:); dLdmod(sstart:send+1,:)];

nfid = fopen('near_10_full.dat');
mfid = fopen('mid_20_full.dat');
ffid = fopen('far_30_full.dat');
zpfid = fopen('zp_full.dat','w');
zsfid = fopen('zs_full.dat','w');
rhofid = fopen('rho_full.dat','w');

tic

for il = 1:ceil(nil/block)
    ndata = fread(nfid, [nt,nrec*block], 'float32');
    mdata = fread(mfid, [nt,nrec*block], 'float32');
    fdata = fread(ffid, [nt,nrec*block], 'float32');

    Data = [ndata(sstart:send,:); mdata(sstart:send,:); fdata(sstart:send,:)];
    [row col] = size(Data);
    fprintf('processing inlines %d to %d of %d\n',(il*block)-block+1,min(il*block,nil),nil);
    
    zp_section = zeros(nt+1,col);
    zs_section = zeros(nt+1,col);
    rho_section = zeros(nt+1,col);

    parfor xl = 1:col
        Data(:,xl) = (1-mu)*Data(:,xl) + mu*G*Model(:,xl);
        [solution,flag,relres,itexit,resvec] = lsqr(G,Data(:,xl),tol,itermax,[],[],Model(:,xl));

        zp_section(:,xl) = exp(solution(1:nt+1,end));
        zs_section(:,xl) = exp(solution(nt+2:2*(nt+1),end)+(k*solution(1:nt+1,end))+kc);
        rho_section(:,xl) = exp(solution((2*(nt+1))+1:end,end)+(m*solution(1:nt+1,end))+mc);

        %fprintf('iteration %d of %d completed\n',xl,nrec);
    end
    
    fwrite(zpfid, zp_section, 'float32');
    fwrite(zsfid, zs_section, 'float32');
    fwrite(rhofid, rho_section, 'float32');
end

toc

fclose all;
matlabpool close
