function [object_meta] = create_object(credentials_path,bucketName,key,object_data)
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

%% Create Object
object_meta.orig_type = 'single';
%object_data = typecast(object_data,object_meta.orig_type); % convert to single so that we always know the original type
object_data = single(object_data);
[object_meta.orig_rows,object_meta.orig_cols] = size(object_data); % could store these alongside or in array?
% convert to signed int8 for use with ByteArray*
object_data = typecast(object_data(:),'int8');
object_meta.typecast_length = length(object_data);

% Set object meta data
meta = ObjectMetadata();
meta.setContentLength(object_meta.typecast_length);
meta.addUserMetadata('orig-type',object_meta.orig_type);
meta.addUserMetadata('orig-rows',num2str(object_meta.orig_rows));
meta.addUserMetadata('orig-cols',num2str(object_meta.orig_cols));
meta.addUserMetadata('type-length',num2str(object_meta.typecast_length));

object_data_ByteArray = ByteArrayInputStream(object_data);
s3.putObject(PutObjectRequest(bucketName, key, object_data_ByteArray, meta));

fprintf('Object created.\n');

end

