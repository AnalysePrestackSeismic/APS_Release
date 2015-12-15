function [h, theta, rho] = hough(varargin)
%HOUGH Hough transform.
%   HOUGH implements the Standard Hough Transform. HOUGH is designed
%   to detect lines. It uses the parametric representation of a line:
%
%                       rho = x*cos(theta) + y*sin(theta).
%
%   The variable rho is the distance from the origin to the line along a
%   vector perpendicular to the line.  Theta is the angle between
%   the x-axis and this vector.
%
%   The Standard Hough Transform (SHT) is a parameter space matrix whose
%   rows and columns correspond to rho and theta values respectively. The
%   elements in the SHT represent accumulator cells. Initially, each cell
%   is set to zero. Then, for every nonbackground point in the image, rho
%   is calculated for every theta. Rho is rounded off to the nearest
%   allowed row in SHT. That accumulator cell is incremented. At the end of
%   this procedure, a value of Q in SHT(r,c) means that Q points in the XY
%   plane lie on the line specified by the theta(c) and rho(r). Peak values
%   in the SHT represent potential lines in the input image.
%
%   [H, THETA, RHO] = HOUGH(BW) computes the SHT of the binary image BW.
%   THETA (in degrees) and RHO are the arrays of rho and theta values over 
%   which the Hough transform matrix, H, was generated.
%
%   [H, THETA, RHO] = HOUGH(BW,PARAM1,VAL1,PARAM2,VAL2) sets various
%   parameters.  Parameter names can be abbreviated, and case does not
%   matter. Each string parameter is followed by a value as indicated
%   below:
%
%   'RhoResolution'   Real scalar between 0 and norm(size(BW)), exclusive.
%                     'RhoResolution' specifies the spacing of the Hough
%                     transform bins along the rho axis.
%
%                     Default: 1.
%
%   'Theta'           Vector of Hough transform theta values. Each element
%                     of the vector specifies the theta value for the
%                     corresponding column of the output matrix H. Theta
%                     values must be within the range [-90, 90) degrees.
%
%                     Default: -90:89
%
%   Notes
%   -----
%   The Hough transform matrix, H, is NRHO-by-NTHETA.
%
%   NRHO = (2*ceil(D/RhoResolution)) + 1, where
%   D = sqrt((numRowsInBW - 1)^2 + (numColsInBW - 1)^2).
%
%   RHO values range from -DIAGONAL to DIAGONAL where
%   DIAGONAL = RhoResolution*ceil(D/RhoResolution).
%
%   THETA values are within the range [-90, 90) degrees.
%
%   Class Support
%   -------------
%   BW can be logical or numeric and it must be real, 2-D, and nonsparse.
%
%   Example 1
%   ---------
%   Compute and display the Hough transform of the gantrycrane.png image
%
%      RGB = imread('gantrycrane.png');
%      I  = rgb2gray(RGB); % convert to intensity
%      BW = edge(I,'canny'); % extract edges
%      [H,T,R] = hough(BW,'RhoResolution',0.5,'Theta',-90:0.5:89.5);
%
%      % display the original image
%      subplot(2,1,1);
%      imshow(RGB);
%      title('gantrycrane.png');
%
%      % display the hough matrix
%      subplot(2,1,2);
%      imshow(imadjust(mat2gray(H)),'XData',T,'YData',R,...
%             'InitialMagnification','fit');
%      title('Hough transform of gantrycrane.png');
%      xlabel('\theta'), ylabel('\rho');
%      axis on, axis normal, hold on;
%      colormap(hot);
%
%   Example 2
%   ---------
%   Compute and display the Hough transform for a limited range of theta
%   values
%
%      RGB = imread('gantrycrane.png');
%      I  = rgb2gray(RGB); % convert to intensity
%      BW = edge(I,'canny'); % extract edges
%      [H,T,R] = hough(BW,'Theta',44:0.5:46);
% 
%      figure;
%      subplot(2,1,1);
%      imshow(RGB);
%      title('gantrycrane.png');
% 
%      subplot(2,1,2); 
%      imshow(imadjust(mat2gray(H)),'XData',T,'YData',R,...
%             'InitialMagnification','fit');
%      title('Limited theta range Hough transform of gantrycrane.png');
%      xlabel('\theta'), ylabel('\rho');
%      axis on, axis normal;
%      colormap(hot)
%
%   See also HOUGHPEAKS and HOUGHLINES.

%   Copyright 1993-2011 The MathWorks, Inc.
%   $Revision: 1.1.8.11.2.1 $  $Date: 2011/12/19 06:30:27 $

%   References:
%   Rafael C. Gonzalez, Richard E. Woods, Steven L. Eddins, "Digital
%   Image Processing Using MATLAB", 2nd ed., Gatesmark Publishing, 2009

%   The function HOUGH changed in Image Processing Toolbox version 6.4
%   (R2009b) to include the 'Theta' parameter. Previous versions allowed
%   more limited control of THETA via the 'ThetaResolution' parameter.
%   'ThetaResolution' may be removed in a future release. The help text
%   below can help you understand existing code that uses
%   'ThetaResolution'. Consider updating your code to use 'Theta' instead.
%
%       'ThetaResolution'   Real scalar between 0 and 90, exclusive.
%                           'ThetaResolution' specifies the spacing 
%                           (in degrees) of the Hough transform bins along 
%                           the theta axis.
% 
%                           Default: 1.
% 
%       NTHETA = 2*ceil(90/ThetaResolution). Theta angle values are in the
%       range [-90, 90) degrees. If 90/ThetaResolution is not an integer,
%       the actual angle spacing will be 90/ceil(90/ThetaResolution).

[bw, rho, theta] = parseInputs(varargin{:});

h = houghmex_(bw,rho,theta*pi/180);

%% Parse Input Parameters
function [bw, rho, theta] = parseInputs(varargin)

%narginchk(1,5);

bw = varargin{1};
validateattributes(bw, {'numeric','logical'},...
    {'real', '2d', 'nonsparse', 'nonempty'}, ...
    mfilename, 'BW', 1);

if ~islogical(bw)
    bw = bw~=0;
end

% Set the defaults
thetaResolution = 1;
rhoResolution = 1;
theta = [];

% Process parameter-value pairs
validStrings = {'RhoResolution', 'Theta', 'ThetaResolution'};
if nargin > 1
    inputStrings = {};
    idx = 2;
    while idx <= nargin
        input = varargin{idx};
        inputStr = validatestring(input, validStrings, mfilename, 'PARAM', idx);
        
        idx = idx+1; %advance index to point to the VAL portion of the input
        
        validateValue(idx, nargin, mfilename, inputStr);
              
        switch inputStr
            case 'RhoResolution'
                rhoResolution = varargin{idx};
                validateRhoResolution(bw, rhoResolution, mfilename, inputStr, idx);
                
            case 'Theta'
                theta = varargin{idx};
                validateTheta(theta, mfilename, inputStr, idx);                            
                
            case 'ThetaResolution'
                thetaResolution = varargin{idx};
                validateThetaResolution(thetaResolution, mfilename, inputStr, idx);
                
            otherwise
                %should never get here
                error(message('images:hough:internalError'))
        end      
        
        idx=idx+1;
        inputStrings = [inputStrings inputStr]; %#ok<AGROW>       
    end
    
    validateThetaStrings(inputStrings, mfilename);
end

% Compute theta and rho
[M,N] = size(bw);

if (isempty(theta))
    theta = linspace(-90, 0, ceil(90/thetaResolution) + 1);
    theta = [theta -fliplr(theta(2:end - 1))];
end

D = sqrt((M - 1)^2 + (N - 1)^2);
q = ceil(D/rhoResolution);
nrho = 2*q + 1;
rho = linspace(-q*rhoResolution, q*rhoResolution, nrho);


%% Validate existence of a value
function validateValue(idx, nargin, ~, inputStr)

if idx > nargin
    error(message('images:hough:missingParamValue', inputStr))
end


%% Validate 'Rho Resolution' parameter
function validateRhoResolution(BW, rhoResolution, mfilename, inputStr, idx)

validateattributes(rhoResolution, {'double'}, {'real', 'scalar', ...
    'finite','positive'}, mfilename, inputStr, idx);

if (rhoResolution >= norm(size(BW)))
    error(message('images:hough:invalidRhoRes', inputStr))
end


%% Validate 'Theta' parameter
function validateTheta(theta, mfilename, inputStr, idx)

validateattributes(theta, {'double'}, {'nonempty', 'real',...
    'vector','finite'}, mfilename, inputStr, idx);

if (min(theta) < -90)
    error(message('images:hough:invalidThetaMin', inputStr))
end

if (max(theta) >= 90)
    error(message('images:hough:invalidThetaMax', inputStr))
end

%% Validate 'Theta Resolution' parameter
function validateThetaResolution(thetaResolution, mfilename, inputStr, idx)

validateattributes(thetaResolution, {'double'}, {'real', 'scalar', ...
    'finite','positive'}, mfilename, inputStr, idx);

if (thetaResolution >= 90)
    error(message('images:hough:invalidThetaRes', inputStr))
end

%% Validate Mutual Exclusivity of 'Theta' and 'Theta Resolution' parameters
function validateThetaStrings(inputStrings, ~)

if (length(strmatch('Theta', inputStrings)) > 1)
    error(message('images:hough:multipleThetaParams'))
end
