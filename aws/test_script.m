%% Test script
load('test_data.mat')

data_up = vol_3d_orig(:,:,10);

figure(1)
subplot(2,1,1); imagesc(data_up);

%%
credentials_path = '/bgdata/git/APS_Release/early_development/aws/creden.txt';
bucketName = '011968-jamess-test';
[bucketName] = create_bucket(credentials_path,bucketName);

%%
key = 'selvage/1-vol_orig.segy';
[obj_meta] = create_object(credentials_path,bucketName,key,data_up);

%%
[object_out,pad] = get_object(credentials_path,bucketName,key,[0,0]);

figure(1)
subplot(2,1,2); imagesc(object_out);

%%
[object_out,pad] = get_object(credentials_path,bucketName,key,[280,180000]);

figure(2)
imagesc(object_out);