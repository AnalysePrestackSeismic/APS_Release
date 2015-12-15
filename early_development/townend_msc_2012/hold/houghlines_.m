function lines = houghlines(varargin)
%HOUGHLINES Extract line segments based on Hough transform.
%   LINES = HOUGHLINES(BW, THETA, RHO, PEAKS) extracts line segments
%   in the image BW associated with particular bins in a Hough 
%   transform.  THETA and RHO are vectors returned by function HOUGH.
%   Matrix PEAKS, which is returned by function HOUGHPEAKS,
%   contains the row and column coordinates of the Hough transform 
%   bins to use in searching for line segments. HOUGHLINES returns
%   LINES structure array whose length equals the number of merged
%   line segments found. Each element of the structure array has
%   these fields: 
%
%      point1  End-point of the line segment; two-element vector
%      point2  End-point of the line segment; two-element vector
%      theta   Angle (in degrees) of the Hough transform bin
%      rho     Rho-axis position of the Hough transform bin
%
%   The end-point vectors contain [X, Y] coordinates.
%
%   LINES = HOUGHLINES(...,PARAM1,VAL1,PARAM2,VAL2) sets various
%   parameters. Parameter names can be abbreviated, and case 
%   does not matter. Each string parameter is followed by a value 
%   as indicated below:
%
%   'FillGap'   Positive real scalar.
%               When HOUGHLINES finds two line segments associated
%               with the same Hough transform bin that are separated
%               by less than 'FillGap' distance, HOUGHLINES merges
%               them into a single line segment.
%
%               Default: 20
%
%   'MinLength' Positive real scalar.
%               Merged line segments shorter than 'MinLength'
%               are discarded.
%
%               Default: 40
%
%   Class Support
%   -------------
%   BW can be logical or numeric and it must be real, 2-D, and nonsparse.
%
%   Example
%   -------
%   Search for line segments corresponding to five peaks in the Hough 
%   transform of the rotated circuit.tif image. Additionally, highlight
%   the longest segment.
%
%      I  = imread('circuit.tif');
%      rotI = imrotate(I,33,'crop');
%      BW = edge(rotI,'canny');
%      [H,T,R] = hough(BW);
%      imshow(H,[],'XData',T,'YData',R,'InitialMagnification','fit');
%      xlabel('\theta'), ylabel('\rho');
%      axis on, axis normal, hold on;
%      P  = houghpeaks(H,5,'threshold',ceil(0.3*max(H(:))));
%      x = T(P(:,2)); 
%      y = R(P(:,1));
%      plot(x,y,'s','color','white');
%
%      % Find lines and plot them
%      lines = houghlines(BW,T,R,P,'FillGap',5,'MinLength',7);
%      figure, imshow(rotI), hold on
%      max_len = 0;
%      for k = 1:length(lines)
%        xy = [lines(k).point1; lines(k).point2];
%        plot(xy(:,1),xy(:,2),'LineWidth',2,'Color','green');
%
%        % plot beginnings and ends of lines
%        plot(xy(1,1),xy(1,2),'x','LineWidth',2,'Color','yellow');
%        plot(xy(2,1),xy(2,2),'x','LineWidth',2,'Color','red');
%
%        % determine the endpoints of the longest line segment 
%        len = norm(lines(k).point1 - lines(k).point2);
%        if ( len > max_len)
%          max_len = len;
%          xy_long = xy;
%        end
%      end
%
%      % highlight the longest line segment
%      plot(xy_long(:,1),xy_long(:,2),'LineWidth',2,'Color','cyan');
%
%   See also HOUGH and HOUGHPEAKS.

%   Copyright 1993-2010 The MathWorks, Inc.
%   $Revision: 1.1.8.10 $  $Date: 2011/08/09 17:50:04 $

%   References:
%   Rafael C. Gonzalez, Richard E. Woods, Steven L. Eddins, "Digital
%   Image Processing Using MATLAB", Prentice Hall, 2003

[nonzeropix,theta,rho,peaks,fillgap,minlength] = parseInputs(varargin{:});

minlength_sq = minlength^2;
fillgap_sq = fillgap^2;
numlines = 0; 
lines = struct;

for k = 1:size(peaks,1)

  % Get all pixels associated with Hough transform cell.
  [r, c] = houghpixels(nonzeropix, theta, rho, peaks(k,:));
  if isempty(r) 
    continue 
  end
  
  % Compute distance^2 between the point pairs
  xy = [c r]; % x,y pairs in coordinate system with the origin at (1,1)
  diff_xy_sq = diff(xy,1,1).^2;
  dist_sq = sum(diff_xy_sq,2);
  
  % Find the gaps larger than the threshold.
  fillgap_idx = find(dist_sq > fillgap_sq);
  idx = [0; fillgap_idx; size(xy,1)];
  for p = 1:length(idx) - 1
    p1 = xy(idx(p) + 1,:); % offset by 1 to convert to 1 based index
    p2 = xy(idx(p + 1),:); % set the end (don't offset by one this time)

    linelength_sq = sum((p2-p1).^2);
    if linelength_sq >= minlength_sq
      numlines = numlines + 1;
      lines(numlines).point1 = p1;
      lines(numlines).point2 = p2;
      lines(numlines).theta = theta(peaks(k,2));
      lines(numlines).rho = rho(peaks(k,1));
    end
  end
end

%-----------------------------------------------------------------------------
function [r, c] = houghpixels(nonzeropix, theta, rho, peak)
%HOUGHPIXELS Compute image pixels belonging to Hough transform bin.
%   [R, C] = HOUGHPIXELS(NONZEROPIX, THETA, RHO, PEAK) computes the
%   row-column indices (R, C) for nonzero pixels NONZEROPIX that map
%   to a particular Hough transform bin, PEAK which is a two element
%   vector [RBIN CBIN].  RBIN and CBIN are scalars indicating the 
%   row-column bin location in the Hough transform matrix returned by
%   function HOUGH.  THETA and RHO are the second and third output 
%   arguments from the HOUGH function.

x = nonzeropix(:,1);
y = nonzeropix(:,2);

theta_c = theta(peak(2)) * pi / 180;
rho_xy = x*cos(theta_c) + y*sin(theta_c);
nrho = length(rho);
slope = (nrho - 1)/(rho(end) - rho(1));
rho_bin_index = round(slope*(rho_xy - rho(1)) + 1);

idx = find(rho_bin_index == peak(1));

r = y(idx) + 1; 
c = x(idx) + 1;
[r,c] = reSortHoughPixels(r, c);

%--------------------------------------------------------------------------
function [r_new, c_new] = reSortHoughPixels(r, c)
% make sure that r an c are in the order along the line segment

if isempty(r)
    r_new = r;
    c_new = c;
    return;
end

r_range = max(r) - min(r);
c_range = max(c) - min(c);
if r_range > c_range
    % Sort first on r, then on c
    sorting_order = [1 2];
else
    % Sort first on c, then on r
    sorting_order = [2 1];
end

[rc_new] = sortrows([r c], sorting_order);
r_new = rc_new(:,1);
c_new = rc_new(:,2);

%-----------------------------------------------------------------------------
function [nonzeropix,theta,rho,peaks,fillgap,minlength] = ...
    parseInputs(varargin)

narginchk(1,8);

idx = 1;
bw = varargin{idx};
validateattributes(bw, {'numeric','logical'},...
              {'real', '2d', 'nonsparse', 'nonempty'}, ...
              mfilename, 'BW', idx);

idx = idx+1;
theta = varargin{idx};
validateattributes(theta, {'double'}, {'real','vector','finite',...
                    'nonsparse','nonempty'}, ...
              mfilename, 'THETA', idx);

idx = idx+1;
rho = varargin{idx};
validateattributes(rho, {'double'}, {'real','vector','finite',...
                    'nonsparse','nonempty'}, ...
              mfilename, 'RHO', idx);

idx = idx+1;
peaks = varargin{idx};
validateattributes(peaks, {'double'}, {'real','2d','nonsparse','integer'}, ...
              mfilename, 'PEAKS', idx);

if size(peaks,2) ~= 2
  error(message('images:houghlines:invalidPEAKS'))
end

% Set the defaults
fillgap = 20; 
minlength = 40; 

% Process parameter-value pairs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
validStrings = {'FillGap','MinLength'};
idx = idx+1;

if nargin > idx-1 % we have parameter/value pairs
  done = false;

  while ~done
    input = varargin{idx};
    inputStr = validatestring(input, validStrings,mfilename,'PARAM',idx);
    
    idx = idx+1; %advance index to point to the VAL portion of the input 
    
    if idx > nargin
      error(message('images:houghlines:valForhoughlinesMissing', inputStr))
    end
    
    switch inputStr
      
     case 'FillGap'
      fillgap = varargin{idx};
      validateattributes(fillgap, {'double'}, {'finite','real', 'scalar', ...
                          'positive'}, mfilename, inputStr, idx);
     
     case 'MinLength'
      minlength = varargin{idx};
      validateattributes(minlength, {'double'}, {'finite','real', 'scalar', ...
                          'positive'}, mfilename, inputStr, idx);
      
     otherwise
      %should never get here
      error(message('images:houghlines:internalError'))
    end
    
    if idx >= nargin
      done = true;
    end
    
    idx=idx+1;
  end
end

% Compute the required parameters
[y, x] = find(bw);
nonzeropix = [x, y] - 1;
