

clear all 


 fid = fopen(['beforeFusion_zone80.dat'],'rt');
        tline = fgetl(fid);
        count=1;
        while ischar(tline)
            line= split(tline,' ');            
             if length(line)>3
                for j=1:length(line)-1
                    LCC{count}(j)=str2num(line{j});
                end
            end     
            
            tline = fgetl(fid);
            count=count+1;
        end
        fclose(fid);
        

size(LCC)


a=load('centroid_and_surface_cells.mat');

fid=fopen('quaternionCells.dat','w');

for i=1:length(a.nuc)
    V=a.nuc{i};
    [PC,~,latent(i,:)]=pca(V); pc1(i,:)=PC(:,1)';pc2(i,:)=PC(:,2)';pc3(i,:)=PC(:,3)';
    %[vec,val]=eig(cov(V));
    warning off; 
    [center(i,:), radii(i,:), evecs{i}, v{i}, chi2(i,:)] = ellipsoid_fit_new(V);
    quat(i,:)=rotm2quat(evecs{i});
    %ankit{i}=evecs{i};
    fprintf(fid,'%0.15f\t%0.15f\t%0.15f\t%0.15f\n',quat(i,[2,3,4,1]));
end
fclose(fid);

fid=fopen('quaternionCluster.dat','w');

for i=1:length(LCC)
    lccid=LCC{i};
    vals=pc1(lccid,:);
    vals=[vals;-vals];
    [coeff,~,~] = pca(vals);  ankit(i,1)=mean(latent(lccid,1));
    q1=rotm2quat(coeff);
    
    vals=pc2(lccid,:);
    vals=[vals;-vals];
    [coeff,~,~] = pca(vals);  ankit(i,2)=mean(latent(lccid,2));
    q2=rotm2quat(coeff);
    
    vals=pc3(lccid,:);
    vals=[vals;-vals];
    [coeff,~,~] = pca(vals);  ankit(i,3)=mean(latent(lccid,3));
    q3=rotm2quat(coeff);
    fprintf(fid,'%0.15f\t%0.15f\t%0.15f\t%0.15f\t%0.15f\t%0.15f\t%0.15f\t%0.15f\t%0.15f\t%0.15f\t%0.15f\t%0.15f\n',[q1([2,3,4,1]),q2([2,3,4,1]),q3([2,3,4,1])]         );
end
    
fclose(fid);

dlmwrite('quaternionLatent.dat',ankit,'\t');
