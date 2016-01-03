clear all

import java.io.*;
% import java.io.BufferedReader;
% import java.io.File;
% import java.io.FileOutputStream;
% import java.io.IOException;
% import java.io.InputStream;
% import java.io.InputStreamReader;
% import java.io.OutputStreamWriter;
% import java.io.Writer;
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
% import com.amazonaws.services.s3.model.GetObjectRequest;
% import com.amazonaws.services.s3.model.ListObjectsRequest;
% import com.amazonaws.services.s3.model.ObjectListing;
% import com.amazonaws.services.s3.model.PutObjectRequest;
% import com.amazonaws.services.s3.model.S3Object;
% import com.amazonaws.services.s3.model.S3ObjectSummary;
import com.amazonaws.util.*;

%% Load AWS Credentials
fileID = fopen('/bgdata/git/APS_Release/early_development/aws/creden.txt');
delimiter = '=';
formatSpec = '%s%[^\n\r]';
creden = textscan(fileID, formatSpec, 'Delimiter', delimiter,  'ReturnOnError', false);
fclose(fileID);

%% Use JAVA Class to create AWS credentials class
awscred = BasicAWSCredentials(aws_access_key_id,);

s3 = AmazonS3Client(awscred);
s3.setEndpoint('s3-eu-west-1.amazonaws.com');

%% Create Bucket
% start buckets with random survey id/block_id/
bucketName = ['survey',num2str(randi(100000))];
block_id = '1';
key = [block_id,'/block/','seismic_gathers_block',block_id,'.segy']; % essentially the filename of the segy
s3.createBucket(bucketName);

%% List buckets
temp_bucketlist = char(s3.listBuckets()); % list buckets and convert to string
temp_bucketlist = strsplit(temp_bucketlist,'S3Bucket');
n_buckets = size(temp_bucketlist,2)-1;

for i_bucket=2:1:n_buckets+1
    temp_bucket = strsplit(temp_bucketlist{i_bucket},{',','[',']'});
    for i_field = 1:1:size(temp_bucket,2)
        temp_cell = strsplit(temp_bucket{i_field},'=');
        if size(temp_cell,2) == 2;
            bucket_array{i_bucket-1}(i_field-1,1) = temp_cell(1);
            bucket_array{i_bucket-1}(i_field-1,2) = temp_cell(2);
        end
    end
    fprintf('Bucket: %s\n',char(bucket_array{i_bucket-1}(1,2)));
end
fprintf('Total buckets = %i\n',n_buckets);
clearvars temp*

%% Read some seismic data
load('matlab_gather_test.mat');
gather_scan = repmat(gather_scan,10,40);
gather_scan = single(gather_scan);
figure(1)
subplot(2,1,1); imagesc(gather_scan);
[orig_rows,orig_cols] = size(gather_scan); % could store these alongside or in array?
% convert to int8
gather_scan = typecast(gather_scan(:),'int8');
length_type = length(gather_scan);

y = ByteArrayInputStream(gather_scan);
%y = IOUtils.toInputStream(gather_scan);
meta = ObjectMetadata();
meta.setContentLength(length_type);
s3.putObject(PutObjectRequest(bucketName, key, y, meta));

% %% Upload to Bucket
% %y = 'james';
% y = IOUtils.toInputStream(char(y));
% meta = ObjectMetadata();
% meta.setContentLength(length(y));
% s3.putObject(PutObjectRequest(bucketName, key, y, meta));

%% Download from bucket
object = s3.getObject(GetObjectRequest(bucketName,key));
byteArray = IOUtils.toByteArray(object.getObjectContent());
% str2double(char(byteArray)); % i get j back
y_down = reshape(typecast(byteArray,'single'),orig_rows,orig_cols);
figure(1)
subplot(2,1,2); imagesc(y_down);

% AmazonS3 s3Client = new AmazonS3Client(new ProfileCredentialsProvider());        
% 
% GetObjectRequest rangeObjectRequest = new GetObjectRequest(
% 		bucketName, key);
% rangeObjectRequest.setRange(0, 10); // retrieve 1st 11 bytes.
% S3Object objectPortion = s3Client.getObject(rangeObjectRequest);
% 
% InputStream objectData = objectPortion.getObjectContent();
% // Process the objectData stream.
% objectData.close();

%% Delete object then bucket
% s3.deleteObject(DeleteObjectRequest(bucketName, key));	
% s3.deleteBucket(bucketName);