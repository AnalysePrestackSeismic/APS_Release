function [byteArray] = get_object(credentials_path,bucketName,key,original_type,range)
%% Inputs
% original_type, string specifying for example 'single'
% should we store original size of object in the bucket?

%% Import JAVA Classes
import java.io.*;
import java.util.UUID;

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
        object = s3.getObject(GetObjectRequest(bucketName,key));
    elseif range(1) > 0 || range(2) > 0
        object_request = GetObjectRequest(bucketName,key);
        object_request.setRange(range(1),range(2));
        object = s3.getObject(object_request);
    elseif range(1) < 0 || range(2) < 0
        fprintf('Range cannot have negative values.\n');
    end    
elseif size(range,2) ~= 2
    fprintf('Range should be 1x2 matrix, start byte | end byte.\n');
end

byteArray = IOUtils.toByteArray(object.getObjectContent());
% Convert object back to original type
byteArray = typecast(byteArray,original_type);
fprintf('Object downloaded.\n');

end

