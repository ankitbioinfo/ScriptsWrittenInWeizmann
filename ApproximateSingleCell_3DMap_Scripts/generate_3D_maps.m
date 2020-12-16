function [] = generate_3D_maps(opt)
%GENERATE_3D_MAPS Generate 3D maps of nuclei and cells features
%   the input opt contains information about the type of maps that should
%   be generated.
%   the output of the function is a set of maps

gp = length(opt.path);
for i = 1:gp
    % compute features on the grid
    if opt.nuclei
        if exist([opt.path{i},'nuclei_grid.mat'], 'file') == 0
            % compute and save nuclei features on the 3D grid
            nuclei_features(opt,i,'nuclei');
        end
    end
    
    if opt.cells
        if exist([opt.path{i},'cells_grid.mat'], 'file') == 0
            % compute and save cell features on the 3D grid
            cells_features(opt,i,'cell');
        end
    end
    
    if opt.crossed
        if exist([opt.path{i},'crossed_grid.mat'], 'file') == 0
            % compute and save crossed features on the 3D grid
            crossed_features(opt,i);
        end
    end
end

% compute visualizations
% compute min and max
[visualization_nuclei,visualization_cells,visualization] = compute_min_max(opt);

% draw maps
for i = 1:gp
       disp(opt.path{i})
       
       % Each gp RZ should be negative z-axis and HZ should be positive
       % axis. In this bonetype=0 
       % If RZ is postive z-axis and HZ is negative z-axis side then growth
       % need to flipped. 
       % Till now Distal Tibia is opposite If any other also then 
       % need to flipped 
       s=strfind(opt.path{i},'_DT_'); % find whether it is DT or PT, If DT then bonetype=1 else it is 0 
       if length(s)>0; bonetype=1; else bonetype=0; end
       
       s=strfind(opt.path{i},'_PH_');
       if length(s)>0; is_bone_humerus=1; else is_bone_humerus=0; end
       
    
    if opt.nuclei
        % load nuclei features on the 3D grid
        res_saved = load([opt.path{i},'nuclei_grid.mat']);
        load([opt.path{i},'all_cells_nuclei.mat']);

        opt.save_folder = [opt.path{i},'figures/nuclei/3D_Maps/'];
        draw_maps_nuclei_cells(all_cells_nuclei,res_saved.nuclei_grid,res_saved.grid,visualization_nuclei,opt,i,'nuclei');
        
        
        opt.save_folder = [opt.path{i},'figures/nuclei/spatial_profiles/'];
        draw_spatial_profiles_new(res_saved.nuclei_grid,bonetype, res_saved.grid,opt);
        
        opt.save_folder = [opt.path{i},'figures/nuclei/normalized_avg_individual_spatial_profiles/'];
        draw_normalized_averaged_individual_spatial_profiles(all_cells_nuclei,bonetype,is_bone_humerus,opt,opt.path{i},'nuclei')

    end
    
    if opt.cells
        % load cell features on the 3D grid
        res_saved = load([opt.path{i},'cells_grid.mat']);
        load([opt.path{i},'all_cells.mat']);
        
        opt.save_folder = [opt.path{i},'figures/cells/3D_Maps/'];
        draw_maps_nuclei_cells(all_cells,res_saved.cells_grid,res_saved.grid,visualization_cells,opt,i,'cell');
        
        opt.save_folder = [opt.path{i},'figures/cells/spatial_profiles/'];
        draw_spatial_profiles_new(res_saved.cells_grid,bonetype,res_saved.grid,opt);

        opt.save_folder = [opt.path{i},'figures/cells/normalized_avg_individual_spatial_profiles/'];
        draw_normalized_averaged_individual_spatial_profiles(all_cells,bonetype,is_bone_humerus,opt,opt.path{i},'cell')

    end
    
    if opt.crossed
        % load crossed features on the 3D grid
        res_saved = load([opt.path{i},'crossed_grid.mat']);
        
        opt.save_folder = [opt.path{i},'figures/crossed/3D_Maps/'];
        draw_maps_crossed(res_saved.crossed_grid,res_saved.grid,visualization,opt);
        
        opt.save_folder = [opt.path{i},'figures/crossed/spatial_profiles/'];
        draw_spatial_profiles_new(res_saved.crossed_grid,bonetype,res_saved.grid,opt);
      
        
        opt.save_folder = [opt.path{i},'figures/crossed/normalized_avg_individual_spatial_profiles/'];
        load([opt.path{i},'all_crossed.mat']);     
        draw_normalized_averaged_individual_spatial_profiles(all_crossed,bonetype,is_bone_humerus,opt,opt.path{i},'crossed')


    end
end

end

