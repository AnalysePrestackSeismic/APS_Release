function seis = segy_read_files(scan_files,index_files,files_in,il_byte,xl_byte,varargin)


% for extra bytes
if (size(varargin,2) > 0)
        if (scan_files == 1)
                for i = 1:length(index_files)                   
                    seis{i} = segy_make_index(cell2mat(files_in.path(index_files(i))),cell2mat(files_in.names(index_files(i))), ...
                        il_byte,xl_byte,cell2mat(varargin));                
                end
            elseif (scan_files == 0)
             
            seis{1} = segy_make_index(cell2mat(files_in.path(index_files(1))),cell2mat(files_in.names(index_files(1))), ...
                il_byte,xl_byte,cell2mat(varargin));   

            for i = 2:length(index_files)
                seis{i} = seis{1};
                seis{i}.filepaths = strcat(files_in.path(index_files(i)),'/',cell2mat(files_in.names(index_files(i))));
            end
        end
        
else
        if (scan_files == 1)
            for i = 1:length(index_files)
                seis{i} = segy_make_index(cell2mat(files_in.path(index_files(i))),cell2mat(files_in.names(index_files(i))),il_byte,xl_byte);
            end
        elseif (scan_files == 0)
               seis{1} = segy_make_index(cell2mat(files_in.path(index_files(1))),cell2mat(files_in.names(index_files(1))),il_byte,xl_byte);   

            for i = 2:length(index_files)
                seis{i} = seis{1};
                seis{i}.filepaths = strcat(files_in.path(index_files(i)),'/',cell2mat(files_in.names(index_files(i))));
            end
        end   
end

% Add option to save scanned files as matlab .mat


% Get more information about the type of segy for error checking
for i = 1:length(index_files)
    fprintf('\nFor file %s enter: \n',seis{i}.filepaths)
    seis{i}.type = input('Enter file type [1 - Full stack, 2 - Angle stack, 3 - Velocity, 4 - Attribute, 5 - Background model]: ');
        if  (seis{i}.type == 2)
            seis{i}.angle = input('Angle in degrees: ');
        end
        if  (seis{i}.type == 4)
            seis{i}.attribtype = input('Enter attribute description: ');
        end
end

end