function [byteArray,padded] = get_object(credentials_path,bucketName,key,range)
%% Inputs
% original_type, string specifying for example 'single'
% should we store original size of object in the bucket?

%% Import JAVA Classes
import java.io.*;
import java.util.UUID;
import org.apache.commons.io.*;

import com.amazonaws.AmazonClientException;
import com.amazonaws.AmazonServiceException;
import com.amazonaws.auth.*;
import com.amazonaws.regions.Region;
import com.amazonaws.regions.Regions;
import com.amazonaws.services.s3.AmazonS3;
import com.amazonaws.services.s3.AmazonS3Client;
import com.amazonaws.services.s3.model.*;
import com.amazonaws.util.*;

%% Use JAVA Class to create AWS credentials class
[aws_id,aws_key] = get_credentials(credentials_path);
awscred = BasicAWSCredentials(aws_id,aws_key);

s3 = AmazonS3Client(awscred);
s3.setEndpoint('s3-eu-west-1.amazonaws.com');

%% Get Object
if size(range,2) == 2
    if range(1) < 0 || range(2) < 0
        fprintf('Range cannot have negative values.\n');
        byteArray = NaN;
        padded.start = NaN;
        padded.end = NaN;
    elseif range(1)/4 ~= floor(range(1)/4) || range(2)/4 ~= floor(range(2)/4)
        fprintf('Ranges are not divisible by 4\n');
        byteArray = NaN;
        padded.start = NaN;
        padded.end = NaN;
    else
        % Get object
        object_request = GetObjectRequest(bucketName,key);
        if range(1) > 0 || range(2) > 0
            range_request = range(2)-range(1)+1;
            
            % check if request range is multiple of 4 elements
            if range_request/4 == floor(range_request/4);
                object_request.setRange(range(1),range(2));
            else
                fprintf('Error the range request is not correct for conversion from int8 to single.\n');
                floor_range = rem(range_request,4);
                range(2) = range(2)-floor_range;
                fprintf('Flooring to nearest element %i.\n',range(2));
                object_request.setRange(range(1),range(2));
            end
        end
        
        object = s3.getObject(object_request);
        % Interpret the meta data
        temp_meta = char(object.getObjectMetadata().getUserMetadata());
        temp_meta = strsplit(temp_meta,{',','{','}'});
        row_i = 1;
        for i_field = 1:1:size(temp_meta,2)
            temp_cell = strsplit(temp_meta{i_field},'=');
            if size(temp_cell,2) == 2;
                meta_array{row_i,1} = strtrim(temp_cell{1}); %using { because these are strings
                meta_array{row_i,2} = strtrim(temp_cell{2});
                row_i = row_i+1;
            end
        end
        
        % Use the meta data
        meta.orig_type = cell2mat(meta_array(find(~cellfun('isempty',strfind(meta_array,'orig-type'))),2));
        meta.orig_rows = str2double(cell2mat(meta_array(find(~cellfun('isempty',strfind(meta_array,'orig-rows'))),2)));
        meta.orig_cols = str2double(cell2mat(meta_array(find(~cellfun('isempty',strfind(meta_array,'orig-cols'))),2)));
        meta.type_length = str2double(cell2mat(meta_array(find(~cellfun('isempty',strfind(meta_array,'type-length'))),2)));
        
        % Get Object
        byteArray = IOUtils.toByteArray(object.getObjectContent());
        
        % Convert object back to original type and size
        byteArray = typecast(byteArray,meta.orig_type);
        
        if size(byteArray,1) == meta.orig_rows*meta.orig_cols
            byteArray = reshape(byteArray,meta.orig_rows,meta.orig_cols);
            padded.start = 0;
            padded.end = 0;
        else
            fprintf('Restricted range request. Padding with NaNs to whole number of rows.\n')
            % Limit length of ByteArray to complete rows, since this will
            % mean complete traces
            start_pad = NaN(range(1)/4,1,'single');
            byteArray = [start_pad; byteArray];
            end_pad = NaN(ceil(size(byteArray,1)/meta.orig_rows)*meta.orig_rows-size(byteArray,1),1,'single');
            byteArray = [byteArray; end_pad];
            byteArray = reshape(byteArray,meta.orig_rows,[]);
            padded.start = size(start_pad,1);
            padded.end = size(end_pad,1);
        end
        fprintf('Object downloaded.\n');
    end
elseif size(range,2) ~= 2
    fprintf('Range should be 1x2 matrix, start byte | end byte.\n');
    byteArray = NaN;
    padded.start = NaN;
    padded.end = NaN;
end

end
