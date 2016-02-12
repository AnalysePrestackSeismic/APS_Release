function [ZNorm PropNorm] = DataSamplingNormV2(Z,Prop,WinNo)
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
%DATASAMPLINGNORM Summary of this function goes here
%   Function which outputs a uniformly distributed data sample

DataIn=[Z Prop];
[Y,I]=sort(DataIn(:,1));
DataSort=DataIn(I,:);

N=length(Z(:,1)); % Determine number of datapoints
%Zmin=min(Z)-1;
%Zmax=max(Z);
%Zint=(Zmax-Zmin)/WinNo;
Nint=floor(N/WinNo);
Nstart=1;
Nend=0;
for ii=1:WinNo
    NumPoints(ii,1)=ii;
    Nend=Nint+Nstart-1;
    NumPoints(ii,2)=min(DataSort(Nstart:Nend,1)); % Find minimum Z in Window N 
    NumPoints(ii,3)=max(DataSort(Nstart:Nend,1)); % Find maxium Z in Window N
    Nstart=Nend+1;
end

NumPoints(:,4)=NumPoints(:,3)-NumPoints(:,2); % Find Z range in Window N
NumPoints(:,5)=round(max(NumPoints(:,4))./NumPoints(:,4)); % Calculate decimation based on Z range over length N samples (Smaller Z range = denser sampling)
NumPoints(:,6)=ceil(Nint./NumPoints(:,5));

start_n=1;
for jj=1:WinNo
    end_n=start_n+NumPoints(jj,6)-1;
    ZNorm(start_n:end_n,1)=DataSort(((jj-1)*Nint)+1:NumPoints(jj,5):jj*Nint,1);
    PropNorm(start_n:end_n,1)=DataSort(((jj-1)*Nint)+1:NumPoints(jj,5):jj*Nint,2);
    start_n=end_n+1;
end

end


