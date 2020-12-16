function [] = crossed_features(opt,gp)
%CROSSED_FEATURES Computes individual crossed features and their aggregated
%values on a 3D grid

if exist([opt.path{gp},'all_crossed.mat'], 'file') == 0
    
    % load coordinates
    [numbers,txt,raw] = xlsread([opt.path{gp},'Tile_coordinates.xlsx']);
    
    if exist([opt.path{gp},'Alignment_matrix.dat'],'file')==0
        registration=eye(3);
    else
	M = load([opt.path{gp},'Alignment_matrix.dat']);
	registration=M(:,[1:3]);
    end
    disp('Registration Matrix')
    registration 
    
    
    coordinates = zeros(size(txt,1)-3,5);
    for i = 4:size(txt,1),
        temp =  char(txt(i,1));
        res = strsplit(temp,'_POS');
        coordinates(i-3,1) = str2num(char(res(2)));
        coordinates(i-3,2:5) = numbers(i-3,:);
    end
    
    disp('Computing individual crossed features')
    
    % compute individual cell features
    all_crossed = [];
    
    % loop over all the positions
    for position = coordinates(:,1)'
        disp(num2str(position));
        
        % load individual cell features
        load([opt.path{gp},'c_n_pos',sprintf('%01d',position),' (Characteristics).mat']);
        
        % here need to check if nuclear cell ratio is good
        indtemp = find(G.inter.volume_ratio>1);
        if length(indtemp)>0,
            
            all_crossed_temp = zeros(length(indtemp),7);
            all_crossed_temp(:,1) = position*ones(length(indtemp),1);
            all_crossed_temp(:,2:4) = G.cel.centroids(indtemp,:);
            
            % organizing individual cell position in the global coordinate
            % system
            indtemp1 = find(coordinates(:,1) == position);
            
            if opt.flip_x_axis{gp},
                all_crossed_temp(:,2) = - all_crossed_temp(:,2) - coordinates(indtemp1,3)*ones(size(all_crossed_temp,1),1);
            else
                all_crossed_temp(:,2) = all_crossed_temp(:,2)+coordinates(indtemp1,3)*ones(size(all_crossed_temp,1),1);
            end
            all_crossed_temp(:,3) = all_crossed_temp(:,3)-coordinates(indtemp1,2)*ones(size(all_crossed_temp,1),1);
            all_crossed_temp(:,4) = all_crossed_temp(:,4)+coordinates(indtemp1,4)*ones(size(all_crossed_temp,1),1);
            
            % proportion of corresponding cells compared to number of
            % nuclei only
            indtemp1 = find(G.nuc.index);
            all_crossed_temp(:,5) = (length(indtemp)/length(indtemp1))*ones(length(indtemp),1);
            
            % n/c ratio
            all_crossed_temp(:,6) = 1./G.inter.volume_ratio(indtemp);
            
            % Centroid shift
            %all_crossed_temp(:,7) = sqrt(G.inter.centroid_shift_physical_coords(indtemp,1).^2+G.inter.centroid_shift_physical_coords(indtemp,2).^2+G.inter.centroid_shift_physical_coords(indtemp,3).^2);
            Paul_value=sqrt(G.inter.centroid_shift_physical_coords(indtemp,1).^2+G.inter.centroid_shift_physical_coords(indtemp,2).^2+G.inter.centroid_shift_physical_coords(indtemp,3).^2);
            
            % method 1 
%             volume= G.cel.volume(1:size(G.inter.volume_ratio,1),:);
%             eq_sphere_radii=(3*volume/(4*pi)).^(1/3);
%             all_crossed_temp(:,7) = 100*value./eq_sphere_radii(indtemp);
            
            % method 2 
            %ellipsoid_radii=G.cel.ellipsoid_radii(1:size(G.inter.volume_ratio,1),:);
            %maximum_ellipsoid_radii=(prod(ellipsoid_radii,2)).^(1/3);
            %all_crossed_temp(:,7) = 100*value./maximum_ellipsoid_radii(indtemp);

            % method 3 
            syms 'qt';
            for indexOfInter=1:size(G.inter.volume_ratio,1)
%              ellipsoid equation Ax^2 + By^2 + Cz^2 + 2Dxy + 2Exz + 2Fyz + 2Gx + 2Hy + 2Iz + J = 0
%              [X,Y,Z,1]  [A D E G;   [X
%                          D B F H;    Y         ( 1 x 4 ) (4 x 4) (4 x 1)
%                          E G C I;    Z                    = (1 x 1) 
%                          G H I J]    1]; 
%
%                     A = [ v(1) v(4) v(5) v(7); ...
%                           v(4) v(2) v(6) v(8); ...
%                           v(5) v(6) v(3) v(9); ...
%                           v(7) v(8) v(9) v(10) ];
%                    A=ep(1);B=ep(2);C=ep(3);D=ep(4);E=ep(5);F=ep(6);G=ep(7);H=ep(8);I=ep(9);J=ep(10); 
                    
                    ep=G.cel.ellipsoid_v{indexOfInter};
                    nucCentroids =G.nuc.centroids(indexOfInter,:);
                    celCentroids =G.cel.centroids(indexOfInter,:);
                    l=[celCentroids-nucCentroids];

                    qx=celCentroids(1)+qt*l(1); qy=celCentroids(2)+qt*l(2); qz=celCentroids(3)+qt*l(3);
                    eqn=ep(1)*qx^2 + ep(2)*qy^2 + ep(3)*qz^2 + 2*ep(4)*qx*qy + 2*ep(5)*qx*qz + 2*ep(6)*qy*qz + 2*ep(7)*qx + 2*ep(8)*qy + 2*ep(9)*qz + ep(10) == 0;
                    solt = vpasolve(eqn,qt);
                    
                    for i=1:length(solt)
                        var=double(solt(i));
                        P(i,:)=[celCentroids(1)+var*l(1),celCentroids(2)+var*l(2),celCentroids(3)+var*l(3)];
                    end
                    
                    if imag(solt(1))==0
                        normalizedDist(indexOfInter,1)=pdist(P);
                    else
                        normalizedDist(indexOfInter,1)=max(G.cel.ellipsoid_radii(indexOfInter,:)); % some extreme case not able to find the solution. 
                        %So dividing by largest radii of ellipsoid
                    end
            end
            
%             size(Paul_value)
%             size(normalizedDist)
            all_crossed_temp(:,7) = 100*Paul_value./normalizedDist(indtemp);
            clear normalizedDist
          
           
            

            [NPCA_coeff, NPCA_latent] = calc_principal_components(N.masks, N.spacing, N.surfaces, registration);

            % PC nuclei
            vals = cell(1,3);
            for pc = 1 : 3
                %tmp = cell2mat(cellfun(@(x) x(:,pc)', G.nuc.PCA_coeff(indtemp), 'UniformOutput', false));
                tmp = cell2mat(cellfun(@(x) x(:,pc)', NPCA_coeff(indtemp), 'UniformOutput', false));
                vals{pc} = [vals{pc}; tmp];
            end
            
            all_crossed_temp(:,8:10) = vals{1}(:,1:3);
            all_crossed_temp(:,11:13) = vals{2}(:,1:3);
            all_crossed_temp(:,14:16) = vals{3}(:,1:3);
            %all_crossed_temp(:,17:19) = G.nuc.PCA_latent(indtemp,:);
            all_crossed_temp(:,17:19) = NPCA_latent(indtemp,:);
            clear NPCA_coeff
            clear NPCA_latent            

            
            [CPCA_coeff, CPCA_latent] = calc_principal_components(C.masks, C.spacing, C.surfaces, registration);

            % PC Cells
            vals = cell(1,3);
            for pc = 1 : 3
                %tmp = cell2mat(cellfun(@(x) x(:,pc)', G.cel.PCA_coeff(indtemp), 'UniformOutput', false));
                tmp = cell2mat(cellfun(@(x) x(:,pc)', CPCA_coeff(indtemp), 'UniformOutput', false));
                vals{pc} = [vals{pc}; tmp];
            end
            
            all_crossed_temp(:,20:22) = vals{1}(:,1:3);
            all_crossed_temp(:,23:25) = vals{2}(:,1:3);
            all_crossed_temp(:,26:28) = vals{3}(:,1:3);
            all_crossed_temp(:,29:31) = CPCA_latent(indtemp,:);
            %all_crossed_temp(:,29:31) = G.cel.PCA_latent(indtemp,:);
            clear CPCA_coeff
            clear CPCA_latent  
            
            if opt.flip_x_axis{gp},
                all_crossed_temp(:,8) = - all_crossed_temp(:,8);
                all_crossed_temp(:,11) = - all_crossed_temp(:,11);
                all_crossed_temp(:,14) = - all_crossed_temp(:,14);
                all_crossed_temp(:,20) = - all_crossed_temp(:,20);
                all_crossed_temp(:,23) = - all_crossed_temp(:,23);
                all_crossed_temp(:,26) = - all_crossed_temp(:,26);
            end
            
            all_crossed = [all_crossed;all_crossed_temp];
        end
    end
    
    DTPT=all_crossed(:,[2,3,4]);
    all_crossed(:,[2,3,4])=DTPT*registration;



    disp('Done');
    
    save([opt.path{gp},'all_crossed.mat'],'all_crossed');
else
    load([opt.path{gp},'all_crossed.mat']);
end

% compute crossed features on 3D grid
characteristics_on_grid_crossed(all_crossed,opt,gp);

end
