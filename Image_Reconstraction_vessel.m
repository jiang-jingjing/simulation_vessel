% - add PATHS (temporary)
addpath(genpath('/home/jingjing/Documents/code2022/PioneerImageReconstruction/'))
% define nirfast paths
% Nirfaster: GPU-fascilitated model
% Nirfast8: old CPU version
pathNirfaster = '/home/jingjing/Documents/SoftwarePackages/Nirfaster';
pathNirfast8 = '//home/jingjing/Documents/SoftwarePackages/nirfast8';
% folder for generated results
fldr_out = '/home/jingjing/Documents/code2022/IR_vessel/output/'
%% load the mesh
% addpath(genpath(pathNirfast8))
fn_mesh = [fldr_out 'mesh_HeadTop_fine.mat'];
% mesh=load(fn_mesh);
% mesh = mesh.mesh;

% %%  mesh from Djazia
% addpath(genpath(pathNirfaster))
% fldr_mesh_Djazia='/home/jingjing/Documents/code2022/IR_vessel/fromDjazia/mesh_new/mesh';
% mesh_pre = load_mesh([fldr_mesh_Djazia '/img_1_nirfast_mesh'])
% nirfast2vtk([fldr_mesh_Djazia '/img_1_nirfast_mesh'], mesh_pre)
%% create a mesh 
addpath(genpath(pathNirfast8))
% load source and detectors
load('detSrcPosFEM.mat')
% voxel.sx=0.39;
% voxel.sy=0.55;
% voxel.sz=0.39;
voxel.sx=0.48828125; %for dimension 512
voxel.sy=0.5; %for dimension 360
voxel.sz=0.348827778; %136
gis_args.medfilter=1;
fldr_mask = '/home/jingjing/Documents/code2022/IR_vessel/tiff_0517/';
% fldr_mask = '/home/jingjing/Documents/code2022/IR_vessel/tiff/';
fn_mask = 'img_1';
type_mask = '.tif';
path_mask = {fldr_mask, fn_mask, type_mask}
pos.center = modelCenter;
pos.src.coordinates = PosSrcFEM;
pos.det.coordinates = PosDetFEM;
mesh = create_mesh_HeadTop(pos, ...
    path_mask,...
    voxel,...  
    2, 1  , 1 , 2);
% mesh.region = mesh.region-1;
save(fn_mesh,'mesh')
%% visualize mesh with sources and detectors
h_mesh = figure()
% % plot boundary nodes
% plot3(mesh.nodes(logical(mesh.bndvtx),1),...
%     mesh.nodes(logical(mesh.bndvtx),2),...
%     mesh.nodes(logical(mesh.bndvtx),3),'c.');
% hold on
% plot sources
plot3(mesh.source.coord(:,1),mesh.source.coord(:,2),...
    mesh.source.coord(:,3),'rx','LineWidth',2,'MarkerSize',8);
hold on
%plot detectors
plot3(mesh.meas.coord(:,1), mesh.meas.coord(:,2),mesh.meas.coord(:,3),...
    'bo','LineWidth',2,'MarkerSize',8);

%% changing the Optical Pro
val_p12.mua= 0.001;% 1 and 2 are the scalp and the skull
val_p12.mus=1;
val_p12.ri=1.4;
val_p3.mua= 0.002;% 3 is the CSF
val_p3.mus=0.3;
val_p3.ri=1.4;
val_p45.mua= 0.01;% 4  are the GM and the WM
val_p45.mus=1;
val_p45.ri=1.37;
val_p7.mua= 0.235;% 5 is the veins 
val_p7.mus=0.512;
val_p7.ri=1.4;

mesh_1 = mesh;
mesh_1 = set_mesh(mesh_1,1, val_p12);
mesh_1 = set_mesh(mesh_1,2, val_p12);
mesh_1 = set_mesh(mesh_1,3, val_p3); 
mesh_1 = set_mesh(mesh_1,4, val_p45);
% mesh_1 = set_mesh(mesh_1,5, val_p45);

% mesh_1 = set_mesh(mesh_1,6, val_p45);% not verified ventricle


mesh_1 = set_mesh(mesh_1,5, val_p12); % assume the same with scalp

%mesh_data = femdata_stnd_tr(mesh_1,7e-9,0.05e-9);%ploting the temprol data
fn_mesh_base =  [fldr_out, 'mesh_base'];
addpath(genpath(pathNirfaster))
nirfast2vtk(fn_mesh_base, mesh_1)
sources2vtk([fn_mesh_base '_source'] ...
    ,mesh_1)
detectors2vtk([fn_mesh_base '_detectors'],mesh_1)



% %% check mask
% gis_args.medfilter=1;
% fldr_mask = '/home/jingjing/Documents/code2022/IR_vessel/tiff/';
% fn_mask = 'img_1';
% type_mask = '.tif';
% [mask param] = GetImageStack([fldr_mask '/' fn_mask type_mask],...
%     gis_args);
% volumeViewer(mask)
% 
% %%
% figure,
% for ii = 35:500
% imagesc(squeeze(mask(:,:,ii)))
% title(num2str(ii))
% % pause(0.3)
% end


%%
% % gis_args.medfilter=1;
% % fldr_mask = '/home/jingjing/Documents/code2022/IR_vessel/fromDjazia/DICOMTOP'
% % type_mask = '.dcm';
% % [mask param] = GetImageStack([fldr_mask '/' fn_mask type_mask],...
% %     gis_args);
%% case 1: contrast : vessel
mesh_2=mesh_1;
mesh_2 = set_mesh(mesh_2,5, val_p7);
fn_mesh_contrast = [fldr_out '/mesh_contrast']
nirfast2vtk(fn_mesh_contrast, mesh_2)
%mesh_1
 %add path for GPU-version Nirfaster
addpath(genpath(pathNirfaster))
% Frequency domain
OPTIONS = solver_options;
OPTIONS.tolerance = 1e-12;
frequencies = 100*1e6;
 
% val_test.mua = 0.01;
% val_test.mus = 1;
% val_test.ri =1.37;
% mesh_test = set_mesh(mesh_1,0, val_test);
% for iregion = 1:7
% mesh_test = set_mesh(mesh_test,iregion, val_test);
% end
tic; [dataSIM_fd_base] = femdata_stnd_FD(mesh_1,...
    frequencies(1), OPTIONS); 
toc; 

tic; [dataSIM_fd_contrast] = femdata_stnd_FD(mesh_2,...
    frequencies(1), OPTIONS); 
toc; 
figure,
subplot(121)
semilogy(dataSIM_fd_base.amplitude)
hold on
semilogy(dataSIM_fd_contrast.amplitude)
title('amplitude')
subplot(122)
plot(dataSIM_fd_base.phase)
hold on
plot(dataSIM_fd_contrast.phase)
title('phase')
legend('base','vessel')
%% image reconstruction
xc = mean(mesh_1.nodes(:,1));
yc = mean(mesh_1.nodes(:,2));
zmax = max(mesh_1.nodes(:,3));
z0 = zmax -5;

mesh_coarse = [200 165 50] /3;
n_iter = 15;
lambda =  1e-1;
% STND
tic
[meshRec_fd_SIM, pj_error] = reconstruct_stnd_FD(...
            mesh_1,frequencies(1), ...
            dataSIM_fd_contrast, 'mua',...
            [],mesh_coarse,n_iter, lambda);
toc
 
coords4plot = [xc yc  z0];
plot_truthVSrecon(mesh_2, ...
        meshRec_fd_SIM,  coords4plot); 
nirfast2vtk( './mesh_SIM_rec', meshRec_fd_SIM)

%% visualize the results in paraview
% generation of the base line forward results
