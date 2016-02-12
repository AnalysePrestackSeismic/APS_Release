function [m] = ModelFitv2(d,G,betaPor,max)
%MODELFIT Summary of this function goes here
%   Detailed explanation goes here

m1=G\d;
m=m1+max;



