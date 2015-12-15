classdef dTectBinaryRead < handle
    
   properties (SetAccess = protected)
      FileName
      FilePath
      NumberInlines
      NumberCrosslines
      NumberSamples
      Angle
      FileSize
      FileFormat
      MachineFormat = 'native';
   end % public properties
   
   % Property data is private to the class
    properties (SetAccess = private, GetAccess = private)
        FileID;                   % file id for reading
    end % private properties
   
   methods
      function obj = dTectBinaryRead(FileName)
            switch nargin
                case 0 % create object but don't link to a file
                    return; % return default object
                case 1 % filename and path provided
                    obj.FileName = FileName;
                    obj.FileFormat = 'float32';
                otherwise
                    error('SeismicFileReader:TooManyInputs',...
                        'Input should be a FileName,NumberInlines,NumberCrosslines,NumberSamples,Angle or empty.')
            end
            
        openFile(obj,FileName)
        end % constructor 
        
        function d = readTraceData(obj,fmt,varargin)
%READTRACEDATA reads the trace numeric data
%   H = READTRACEDATA(OBJ,TRACES,OFFSET) reads in the file header
%   information for OBJ (the SEISMICFILEREADER object) according to the
%   number of traces in the scalar or vector TRACES. OFFSET is
%   the byte offset from the beginning of file to begin the read from.
%
%   Example:
%   % read in part of the trace header for a SEG2 formatted file.
%   s = SeismicFileReader('WFLT0001')
%
%   % read in the file header information to get the pointers to the trace
%   % data
%   fmt = {'uint16',{'FileDescriptor','RevisionNumber',...
%                     'TracePointerSize','NumberOfTraces'};
%          'uint8', {6, 'Terminators'} };
%   fileHeader = readFileHeader(s,fmt)
%
%   % get the number of traces and the trace locations in the file
%   fmt = {'uint32', {fileHeader.NumberOfTraces, 'TracePointers'}};
%   offset = 32; % bytes from beginning of file
%   locInFile = readFileHeader(s,fmt,offset);
%
%   % read in the header for trace 1
%   fmt = {'uint16',{'TraceDescriptor','TraceBlockSize'};...
%          'uint32',{'TraceDataSize','NumberOfSamples'};...
%          'uint8', {'TraceDataFormatCode'} };
%   offset = locInFile.TracePointers(1);      
%   traceHeader = readTraceHeader(s,fmt,offset)
%
%   % read in the data for trace 1
%   fmt = {'int32', traceHeader.NumberOfSamples}; % 32-bit integer 
%   offset = 1348; % trace 1 starts at byte 1348
%   traceData = readTraceData(s,fmt,offset)
%
%   See also SEISMICFILEREADER, READFILEHEADER, READTRACEHEADER

            % trace data should be a single numeric array
            if ~isOpen(obj)
                % open file if not already
                openFile(obj)
            end
            
            % fmt is assumed to be a 1x2 or 2x1 vector
            if ~isvector(fmt) || length(fmt) ~= 2
                error('SeismicFileReader:InvalideTraceFormat',...
                    'Format specification must be a vecor of length 2')
            end
            if nargin > 2
                offset = varargin{1};
                if offset >= 0
                    fseek(obj.FileID,offset,'bof');
                end % otherwise do nothing, use current position.
            end
            d = readFile(obj);
        end % readTraceData

   end% methods 
   
   methods (Access = private)
   
   function openFile(obj,fileName)
        %OPENFILE open the file for reading
        %   OPENFILE(OBJ,FILENAME) opens the file in FILENAME and returns
        %   the file id to OBJ.FILEID.  An error is returned if the file
        %   cannot be opened.
            if nargin == 1
                fileName = fullfile(obj.FilePath,obj.FileName);
            end
            
            % If filename is provided, open it for reading if exists
            if ~isempty(fileName)
                assert(exist(fileName,'file') == 2, ...
                    'SeismicFileReader:FileNotFound',...
                    'Could not locate file: %s', fileName)
                
                obj.FileID = fopen(fileName,'r');
                
                assert(obj.FileID > 0, ...
                    'SeismicFileReader:CouldNotOpenFile',...
                    'Could not open file: %s', fileName);
            end
            
            % Set object properties from file, use which to capture path if
            % file is on the path
            
            ff = which(fileName);
            if isempty(ff)
                ff = fileName;
            end
               [obj.FilePath, obj.FileName, fileExt] = fileparts(ff);
               obj.FileName = [obj.FileName fileExt];
               s = dir(ff);
               obj.FileSize = s.bytes;
        end % openFile
   
    function st = isOpen(obj)
        %ISOPEN returns true or false if the file is open.
        %   ISOPEN(OBJ) returns TRUE if the file in OBJ.FILENAME is open,
        %   otherwise it returns FALSE.
        
            % Check if the file ID is valid
            try
                ferror(obj.FileID);
                st = true;
            catch e
                if sum(strcmpi(e.identifier,{'MATLAB:FileIO:InvalidFid',...
                        'MATLAB:badfid_mx'}))
                    st = false;
                else
                    throw(e)
                end
            end 
        end % isOpen
        
    function f = readFile(obj,sizeA,precision)
%READFILE reads format similar to fread
%   F = READFILE(OBJ,SIZEA,PRECISION) reads in data from the file linked to
%   OBJ using SIZEA and PRECISION specifiers (see FREAD for details).
%
%   See also FREAD
            if ~isOpen(obj)
                % open file if not already
                openFile(obj)
            end
            skip = 0;
            f = fread(obj.FileID,sizeA,precision,skip,obj.MachineFormat);
    end
    
    end % private methods
        
end % class dTectBinaryRead