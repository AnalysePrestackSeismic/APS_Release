function out = segy_types(varargin)

% Index of file type
% Description of file type

    types = ...
    {
    1 'Full stack';
    2 'Angle stack';
    3 'RMS velocity';
    4 'Post migration gathers'; % could have sub types
    5 'Frequency gathers'; % could have sub types
    6 'Inline dip|Crossline dip|Azimuth gathers'; % other dips can be calculated from this
    7 'Time reference';
    8 'Other attribute';
    };

    if nargin
        out = types(cell2mat(varargin),2);
    else
        fprintf('\nThe following files types are permitted:\n')
        for ii=1:1:length(types)
            fprintf('[%d] %s\n',types{ii,1},types{ii,2});
        end
    end
end