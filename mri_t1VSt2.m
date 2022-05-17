fldr_MRI = '/home/jingjing/Documents/Data/MRI_images/';
dicom_t1 = [fldr_MRI '/002-3D_BrainVIEW_T1W/DICOM/IM_0001'];
dicom_t2 = [fldr_MRI '/003-3D_Brain_VIEW_T2/DICOM/IM_0001']
% dicom_MRI = [fldr_MRI '/004-s3DI_MC_HR/DICOM/IM_0001']
info_t1 = dicominfo(dicom_t1);
dicom_t1 = dicomread(dicom_t1);
img_t1 = squeeze(dicom_t1);
info_t2 = dicominfo(dicom_t2);
dicom_t2 = dicomread(dicom_t2);
img_t2 = squeeze(dicom_t2);

[nx1 ny1 nz1] = size(img_t1);
[nx2 ny2 nz2] = size(img_t2);

volumeViewer(img_t1)


%%
V = squeeze(VTopSeg);
