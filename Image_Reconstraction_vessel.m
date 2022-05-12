% - add PATHS (temporary)
addpath(genpath('/home/jingjing/Documents/code2022/PioneerImageReconstruction/'))
% define nirfast paths
% Nirfaster: GPU-fascilitated model
% Nirfast8: old CPU version
pathNirfaster = '/home/jingjing/Documents/SoftwarePackages/Nirfaster';
pathNirfast8 = '//home/jingjing/Documents/SoftwarePackages/nirfast8';

%% load the mesh
% addpath(genpath(pathNirfast8))
fn_mesh = 'mesh_HeadTop_fine.mat';
% mesh=load(fn_mesh);
% mesh = mesh.mesh;

% % mesh from Djazia
% mesh_old = load_mesh('/home/jingjing/Documents/code2022/IR_vessel/fromDjazia/Mesh/img_1_nirfast_mesh')
nirfast2vtk('mesh_old', mesh_old)
%% create a mesh 
addpath(genpath(pathNirfast8))
% load source and detectors
load('detSrcPosFEM.mat')
fldr_mask = '../fromDjazia/Test1/bmp/';
fn_mask = 'img_1';
type_mask = '.bmp';
path_mask = {fldr_mask, fn_mask, type_mask}
pos.center = modelCenter;
pos.src.coordinates = PosSrcFEM;
pos.det.coordinates = PosDetFEM;
mesh = create_mesh_HeadTop(pos, path_mask,  2, 1.5 , 1.5, 2);
mesh.region = mesh.region-1;
save(fn_mesh,'mesh')
%% visualize mesh with sources and detectors
h_mesh = figure()
% plot boundary nodes
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
val_p45.mua= 0.01;% 4 and 5 are the GM and the WM
val_p45.mus=1;
val_p45.ri=1.37;
val_p7.mua= 0.235;% 7 is the veins 
val_p7.mus=0.512;
val_p7.ri=1.4;

mesh_1 = mesh;
mesh_1 = set_mesh(mesh_1,1, val_p12);
mesh_1 = set_mesh(mesh_1,2, val_p12);
mesh_1 = set_mesh(mesh_1,3, val_p3);
mesh_1 = set_mesh(mesh_1,4, val_p45);
mesh_1 = set_mesh(mesh_1,5, val_p45);
mesh_1 = set_mesh(mesh_1,7, val_p7);
%mesh_data = femdata_stnd_tr(mesh_1,7e-9,0.05e-9);%ploting the temprol data
fn_mesh_anom = 'mesh_anom'
addpath(genpath(pathNirfaster))
nirfast2vtk(fn_mesh_anom, mesh_1)
sources2vtk([fn_mesh_anom '_source'] ...
    ,mesh_1)
detectors2vtk([fn_mesh_anom '_detectors'],mesh_1)




%% check bmps
gis_args.medfilter=1;
[mask param] = GetImageStack([fldr_mask '/' fn_mask type_mask],...
    gis_args);
volumeViewer(mask)






%% case 1: base
mesh_2=mesh;
mesh_2 = set_mesh(mesh_2,1, val_p12);
mesh_2 = set_mesh(mesh_2,2, val_p12);
mesh_2 = set_mesh(mesh_2,3, val_p3);
mesh_2 = set_mesh(mesh_2,4, val_p45);
mesh_2 = set_mesh(mesh_2,5, val_p45);
mesh_2 = set_mesh(mesh_2,7, val_p12);


%% case 2: contrast
%mesh_1
 %add path for GPU-version Nirfaster
addpath(genpath(pathNirfaster))
% Frequency domain
OPTIONS = solver_options;
OPTIONS.tolerance = 1e-12;
frequencies = 100*1e6;
 
val_test.mua = 0.01;
val_test.mus = 1;
val_test.ri =1.37;
mesh_test = set_mesh(mesh_1,0, val_test);
for iregion = 1:7
mesh_test = set_mesh(mesh_test,iregion, val_test);
end
 
tic; [dataSIM_fd_contrast] = femdata_stnd_FD(mesh_test,...
    frequencies(1), OPTIONS); 
toc; 
plot_data(dataSIM_fd_contrast)
%% image reconstruction


%% visualize the results in paraview
% generation of the base line forward results
