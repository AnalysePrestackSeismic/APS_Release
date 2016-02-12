function [ downlogs ] = downsample(logs,sampint)
%DOWNSAMPLE Summary of this function goes here
%   Detailed explanation goes here
 downlogs = logs(1:sampint:end,:);
end

