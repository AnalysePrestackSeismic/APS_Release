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