function mesh = create_mesh_HeadTop(pos, path_mask, facet_distance, ...
    facet_size, cell_size, cell_radius_edge)
gis_args.medfilter=1;
fldr_mask = path_mask{1};
fn_mask = path_mask{2};
type_mask = path_mask{3};
[mask param] = GetImageStack([fldr_mask '/' fn_mask type_mask],...
    gis_args);

param.facet_angle = (25.0);
param.facet_distance = facet_distance;
% param.facet_size = (1);
param.facet_size = facet_size;


param.medfilter = 0;
param.pad = 0;
% param.cell_size = (0.9);
param.cell_size = cell_size;

% param.cell_radius_edge = (3.0);
param.cell_radius_edge = cell_radius_edge;
param.special_subdomain_label = (0);
param.special_subdomain_size  = (0);

sx=(1);
sy=(1);
sz=(1);

param.PixelDimensions(1) = sx;
param.PixelDimensions(2) = sy;
param.PixelDimensions(3) = sz;
param.PixelSpacing(1) = sx;
param.PixelSpacing(2) = sy;
param.SliceThickness  = sz;

outfn = [fldr_mask '/' fn_mask '.ele'];
param.tmppath = fileparts(outfn);
if isempty(param.tmppath), param.tmppath = getuserdir(); end
param.delmedit = 0;

[e p] = RunCGALMeshGenerator(mask,param);
if size(e,2) > 4, mat = e(:,5); else, mat = ones(size(e,1),1); end

q1 = simpqual(p, e, 'Min_Sin_Dihedral');
genmesh.ele = e;
genmesh.node = p;
genmesh.ele(:,5) = mat; genmesh.nnpe = 4; genmesh.dim = 3;
fn_mesh = [fldr_mask '/' fn_mask '_mesh'];
solidmesh2nirfast(genmesh,fn_mesh,'stnd');
fprintf('done.\n');

mesh = load_mesh(fn_mesh);
clear e f1 f2 genmesh info mask mat e p sx sy sz cr1 cr2 genmesh outfn param
clear save_fn
%--------------------%

mesh_tmp = load_mesh(fn_mesh);
% mesh_tmp.link =link_all;
xc = mean(mesh_tmp.nodes(:,1));
yc = mean(mesh_tmp.nodes(:,2));
zmax = max(mesh_tmp.nodes(:,3));
x_shift = xc-pos.center(1);
y_shift = yc-pos.center(2) + 3;
PosSrcNIRFAST = pos.src.coordinates;
PosSrcNIRFAST(:,1) = PosSrcNIRFAST(:,1)+x_shift;
PosSrcNIRFAST(:,2) = PosSrcNIRFAST(:,2)+y_shift;
PosSrcNIRFAST(:,3) = zmax+1;
PosDetNIRFAST = pos.det.coordinates;
PosDetNIRFAST(:,1) = PosDetNIRFAST(:,1)+x_shift; 
PosDetNIRFAST(:,2) = PosDetNIRFAST(:,2)+y_shift;  
PosDetNIRFAST(:,3) = zmax+1;
mesh_tmp.source.coord =PosSrcNIRFAST;
mesh_tmp.source.num = (1:size(PosSrcNIRFAST,1))';
mesh_tmp.source.fwhm = zeros(size(PosSrcNIRFAST,1),1);
mesh_tmp.source.fixed =0;
mesh_tmp.source.distributed =0;
mesh_tmp.meas.coord = PosDetNIRFAST;
mesh_tmp.meas.num = (1:size(PosDetNIRFAST,1))';
mesh_tmp.meas.fixed =0;
mesh_tmp = minband_opt(mesh_tmp);
save_mesh(mesh_tmp,fn_mesh);
clear mesh_tmp
mesh = load_mesh(fn_mesh);
num_det = size(PosDetNIRFAST, 1)
num_src = size(PosSrcNIRFAST, 1)
 
link_all = [];
for is = 1: num_src
    for id = 1: num_det
        link_tmp = [is id 1];
        link_all = [link_all;link_tmp];
    end
end
mesh.link = link_all;
save_mesh(mesh, fn_mesh);
end