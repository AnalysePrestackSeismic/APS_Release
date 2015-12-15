
function p = cdf1(name,x,varargin)
%CDF Cumulative distribution function for a specified distribution.
%   Y = CDF(NAME,X,A) returns an array of values of the cumulative
%   distribution function for the one-parameter probability distribution
%   specified by NAME with parameter values A, evaluated at the values in X.
%
%   Y = CDF(NAME,X,A,B) or Y = CDF(NAME,X,A,B,C) returns values of the
%   cumulative distribution function for a two- or three-parameter probability
%   distribution with parameter values A, B (and C).
%
%   The size of Y is the common size of the input arguments.  A scalar input
%   functions as a constant matrix of the same size as the other inputs.  Each
%   element of Y contains the cumulative distribution evaluated at the
%   corresponding elements of the inputs.
%
%   NAME can be:
%
%      'beta'  or 'Beta',
%      'bino'  or 'Binomial',
%      'chi2'  or 'Chisquare',
%      'exp'   or 'Exponential',
%      'ev'    or 'Extreme Value',
%      'f'     or 'F',
%      'gam'   or 'Gamma',
%      'gev'   or 'Generalized Extreme Value',
%      'gp'    or 'Generalized Pareto',
%      'geo'   or 'Geometric',
%      'hyge'  or 'Hypergeometric',
%      'logn'  or 'Lognormal',
%      'nbin'  or 'Negative Binomial',
%      'ncf'   or 'Noncentral F',
%      'nct'   or 'Noncentral t',
%      'ncx2'  or 'Noncentral Chi-square',
%      'norm'  or 'Normal',
%      'poiss' or 'Poisson',
%      'rayl'  or 'Rayleigh',
%      't'     or 'T',
%      'unif'  or 'Uniform',
%      'unid'  or 'Discrete Uniform',
%      'wbl'   or 'Weibull'.
%
%   CDF calls many specialized routines that do the calculations.
%
%   See also ICDF, MLE, PDF, RANDOM.

%   Copyright 1993-2006 The MathWorks, Inc.
%   $Revision: 1.1.8.2 $  $Date: 2010/10/08 17:22:39 $

if nargin<2, error(message('stats:cdf:TooFewInputs')); end

if ~ischar(name)
   error('stats:cdf:BadDistribution',...
         'First argument must be distribution name');
end

if nargin<5
    a3=0;
else
    a3 = varargin{3};
end
if nargin<4
    a2=0;
else
    a2 = varargin{2};
end
if nargin<3
    a1=0;
else
    a1 = varargin{1};
end

if     strcmpi(name,'beta'),
    p = betacdf(x,a1,a2);
elseif strcmpi(name,'bino') || strcmpi(name,'Binomial'),
    p = binocdf(x,a1,a2);
elseif strcmpi(name,'chi2') || strcmpi(name,'Chisquare'),
    p = chi2cdf(x,a1);
elseif strcmpi(name,'exp') || strcmpi(name,'Exponential'),
    p = expcdf(x,a1);
elseif strcmpi(name,'ev') || strcmpi(name,'Extreme Value'),
    p = evcdf(x,a1,a2);
elseif strcmpi(name,'f'),
    p = fcdf(x,a1,a2);
elseif strcmpi(name,'gam') || strcmpi(name,'Gamma'),
    p = gamcdf(x,a1,a2);
elseif strcmpi(name,'gev') || strcmpi(name,'Generalized Extreme Value'),
    p = gevcdf(x,a1,a2,a3);
elseif strcmpi(name,'gp') || strcmpi(name,'Generalized Pareto'),
    p = gpcdf(x,a1,a2,a3);
elseif strcmpi(name,'geo') || strcmpi(name,'Geometric'),
    p = geocdf(x,a1);
elseif strcmpi(name,'hyge') || strcmpi(name,'Hypergeometric'),
    p = hygecdf(x,a1,a2,a3);
elseif strcmpi(name,'logn') || strcmpi(name,'Lognormal'),
    p = logncdf(x,a1,a2);
elseif strcmpi(name,'nbin') || strcmpi(name,'Negative Binomial'),
    p = nbincdf(x,a1,a2);
elseif strcmpi(name,'ncf') || strcmpi(name,'Noncentral F'),
    p = ncfcdf(x,a1,a2,a3);
elseif strcmpi(name,'nct') || strcmpi(name,'Noncentral T'),
    p = nctcdf(x,a1,a2);
elseif strcmpi(name,'ncx2') || strcmpi(name,'Noncentral Chi-square'),
    p = ncx2cdf(x,a1,a2);
elseif strcmpi(name,'norm') || strcmpi(name,'Normal'),
    p = normcdf(x,a1,a2);
elseif strcmpi(name,'poiss') || strcmpi(name,'Poisson'),
    p = poisscdf(x,a1);
elseif strcmpi(name,'rayl') || strcmpi(name,'Rayleigh'),
    p = raylcdf(x,a1);
elseif strcmpi(name,'t'),
    p = tcdf(x,a1);
elseif strcmpi(name,'unid') || strcmpi(name,'Discrete Uniform'),
    p = unidcdf(x,a1);
elseif strcmpi(name,'unif')  || strcmpi(name,'Uniform'),
    p = unifcdf(x,a1,a2);
elseif strcmpi(name,'weib') || strcmpi(name,'Weibull') || strcmpi(name,'wbl')
    if strcmpi(name,'weib') || strcmpi(name,'Weibull')
        warning(message('stats:cdf:ChangedParameters'));
    end
    p = wblcdf(x,a1,a2);
else
    spec = dfgetdistributions(name);
    if isempty(spec)
       error('stats:cdf:BadDistribution',...
             'Unrecognized distribution name: ''%s''.',name);
    elseif length(spec)>1
       error('stats:cdf:BadDistribution',...
             'Ambiguous distribution name: ''%s''.',name);
    end

    p = feval(spec.cdffunc,x,varargin{:});
end