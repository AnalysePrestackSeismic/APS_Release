function [byteArray] = get_object(credentials_path,bucketName,key,range)
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
    if range(1) == 0 && range(1) == 0 % Download all of object
        % Get object
        object = s3.getObject(GetObjectRequest(bucketName,key));
    elseif range(1) > 0 || range(2) > 0
        % Get object
        object_request = GetObjectRequest(bucketName,key);
        object_request.setRange(range(1),range(2));
        object = s3.getObject(object_request);
    elseif range(1) < 0 || range(2) < 0
        fprintf('Range cannot have negative values.\n');
    end    
elseif size(range,2) ~= 2
    fprintf('Range should be 1x2 matrix, start byte | end byte.\n');
end
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

% Get meta data
orig_type = cell2mat(meta_array(find(~cellfun('isempty',strfind(meta_array,'orig-type'))),2));
orig_rows = str2double(cell2mat(meta_array(find(~cellfun('isempty',strfind(meta_array,'orig-rows'))),2)));
orig_cols = str2double(cell2mat(meta_array(find(~cellfun('isempty',strfind(meta_array,'orig-cols'))),2)));

% Get Object
byteArray = IOUtils.toByteArray(object.getObjectContent());

% Convert object back to original type and size
byteArray = typecast(byteArray,orig_type);
byteArray = reshape(byteArray,orig_rows,orig_cols);

fprintf('Object downloaded.\n');

end

