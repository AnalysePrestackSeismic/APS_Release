function perltest( )
%PERTEST Summary of this function goes here
%   Detailed explanation goes here
%
%[result, status] = perl(/apps/gsc/scripts/perl_parfor,arg1,arg2);
%
[result, status] = perl('/apps/gsc/scripts/perl_matlab_test.pl','/apps/gsc/scripts/threadtest.sh','1','2','3');
testa = result;
testb = status;
%
end

