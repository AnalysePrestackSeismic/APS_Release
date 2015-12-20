% javaaddpath('/bgdata/git/AWS-MLEP/lib/aws-java-sdk-1.4.7.jar')
% javaaddpath('/bgdata/git/AWS-MLEP/lib/httpclient-4.1.1.jar')
% javaaddpath('/bgdata/git/AWS-MLEP/lib/httpcore-4.1.jar')
% javaaddpath('/bgdata/git/commons-io-2.4/commons-io-2.4.jar')
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

%% 

awscred = BasicAWSCredentials('','');

s3 = AmazonS3Client(awscred);
s3.setEndpoint('s3-eu-west-1.amazonaws.com');

%% Create Bucket

bucketName = ['my-first-bucket-matlab-',num2str(randi(100000))];
key = 'MyObjectKeyMatlab/y_up.char';
s3.createBucket(bucketName);
s3.listBuckets()
%% Upload to Bucket
y = 'james';
y = IOUtils.toInputStream(char(y));
meta = ObjectMetadata();
meta.setContentLength(length(y));
s3.putObject(PutObjectRequest(bucketName, key, y, meta));

%% Download from bucket
%s3.listBuckets
%bucketName = 'download-matlab-91339';
%key = 'test/aws_s3_seismic.segy';
object = s3.getObject(GetObjectRequest(bucketName,key));

byteArray = IOUtils.toByteArray(object.getObjectContent());
char(byteArray); % i get j back
%% Delete object then bucket
s3.deleteObject(DeleteObjectRequest(bucketName, key));	
s3.deleteBucket(bucketName);