function [coeff, latent] = calc_principal_components(masks, spacings, surfaces,rotation)

% Note: in 'coeff', the first column (i.e. coeff(:,1)) is the orientation vector of PC1.

coeff = cell(length(masks),1);
latent = nan(length(masks),3);

for i = 1 : length(masks)
    [x,y,z] = ind2sub(size(masks{i}), find(masks{i}));
    xyz = [x,y,z] .* repmat(spacings(i,:), length(x), 1);
    xyz = xyz - repmat(mean(xyz,1), size(xyz,1), 1);
    try
        [coeff{i,1}, ~, latent(i,:)] = pca(xyz*rotation);
    catch
        keyboard;
    end
    
    if 0
        figure;
%         scatter3(xyz(1:100:end,1), xyz(1:100:end,2), xyz(1:100:end,3), '.');
        s = surfaces(i);
        s.vertices = s.vertices - repmat(mean(s.vertices,1), size(s.vertices,1), 1);
        show_surface(s);
        hold on;
        quiver3(0, 0, 0, coeff{i}(1,1), coeff{i}(2,1), coeff{i}(3,1), 10);
        quiver3(0, 0, 0, coeff{i}(1,2), coeff{i}(2,2), coeff{i}(3,2), 10);
        quiver3(0, 0, 0, coeff{i}(1,3), coeff{i}(2,3), coeff{i}(3,3), 10);
        hold off;
        cameratoolbar;
        axis image;
        axis vis3d;
        axis off;
        close;
    end
end
