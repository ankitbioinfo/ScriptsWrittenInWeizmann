% this script generates individuals and grid scale cells features and
% produces 3D maps and spatial profiles

% This script run with or without registration (Alignment_matrix.dat) file 
% If you provide 3x3 rotation matrix inside each data folder such that Z axis 
% is aligend to proximal-distal axis 
% then spatial profile will be plotted along proximal-distal axis. 
% If you do not provide rotation matrix for registration then it will take
% 3x3 identity matrix and spatial profile might plot in wrong way. 

% add path
root_dir = '';
%initialize path 
dtwtpath=[];
ptwtpath=[];
dtmtpath=[];
ptmtpath=[];
duwtpath=[];
humerus_wt_path=[];
humerus_mt_path=[];

% Subfolders
addpath(fullfile(root_dir, 'data'));
addpath(fullfile(root_dir, 'utils'));


presets



humerus_wt_path={'data/Nuclei_and_Cells_DU_S51_m2_wt/'};          
%humerus_wt_path={'data/Nuclei_and_Cells_PH_S93_m1_het_E185/'};     
% humerus_mt_path={'data\Nuclei_and_Cells_PH_S93_m2_mut_E185\'};
%duwtpath={'data\Nuclei_and_Cells_DU_S84_m2_wt\'};


%opt.path=[dtwtpath,ptwtpath,dtmtpath,ptmtpath,duwtpath];
% opt.path=[humerus_wt_path,humerus_mt_path];
opt.path = [humerus_wt_path];

%temp=strsplit(dtwtpath{1},'Nuc');
mainpath='data\';  % directory for saving statistical results

% The number of zeros should be equivalent to number of sample. 
opt.flip_x_axis = {0;0};

% dimensions of the grid
opt.delta_x = 15;
opt.delta_y = 15;
opt.delta_z = 15;
opt.minimum_no_object_in_grid=1;


% number of bins along P-D axis (normalized averaged individual spatial
% profiles)
opt.ProximalDistal_binsize=50;

% number of bins for Grid average 
opt.grid_bin=6;


% 3D view
opt.az = 63;
opt.el = -7;

% maps for the nuclei - value: 0 or 1 (if they should be included)
opt.nuclei = 1;

% maps for the cells - value: 0 or 1 (if they should be included)
opt.cells = 0;

% maps for crossed features (alignment, shift etc..)
opt.crossed = 0;

% save figs
opt.save_figs = 1;
opt.matlabfig= 1;

opt.plot_3D_map_14_features=1;
opt.plot_3D_map_log_features=1;
opt.plot_PC_orientations=1;
opt.plot_PC_cross_section_orientationX=1;
opt.plot_PC_cross_section_orientationY=1;
opt.plot_PC_cross_section_orientationZ=1;
generate_3D_maps(opt);


temp=strsplit(opt.path{1},'Nuc');
mainpath=temp{1}; % directory for saving statistical results
wildTypesBones={dtwtpath, ptwtpath, duwtpath, humerus_wt_path};
mutantBones={dtmtpath, ptmtpath, humerus_mt_path};
wt_bone_name={'dt-wt','pt-wt','du-wt','ph-wt'};
mt_bone_name={'dt-mt','pt-mt','ph-mt'};


outputpath=[mainpath,'StatisticalTestOutputs\'];
opt.allometric.limit=[8,32,40,opt.ProximalDistal_binsize];

opt.AllStatisticalTest=0;
opt.onlyAllometricTest=0;
%perform_statistical_test_on_avg_profile(wildTypesBones, mutantBones,wt_bone_name, mt_bone_name,outputpath,opt);








