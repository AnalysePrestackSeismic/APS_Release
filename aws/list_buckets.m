function [bucket_array] = list_buckets(credentials_path)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
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

%% List buckets
temp_bucketlist = char(s3.listBuckets()); % list buckets and convert to string
temp_bucketlist = strsplit(temp_bucketlist,'S3Bucket');
n_buckets = size(temp_bucketlist,2)-1;

for i_bucket=2:1:n_buckets+1
    temp_bucket = strsplit(temp_bucketlist{i_bucket},{',','[',']'});
    for i_field = 1:1:size(temp_bucket,2)
        temp_cell = strsplit(temp_bucket{i_field},'=');
        if size(temp_cell,2) == 2;
            bucket_array{i_bucket-1}(i_field-1,1) = strtrim(temp_cell(1));
            bucket_array{i_bucket-1}(i_field-1,2) = strtrim(temp_cell(2));
        end
    end
    fprintf('Bucket: %s\n',char(bucket_array{i_bucket-1}(1,2)));
end
fprintf('Total buckets = %i\n',n_buckets);

end

