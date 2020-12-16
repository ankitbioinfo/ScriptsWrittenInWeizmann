
clear all 

a=load('Tibia/dataSave/SaveVariablesForPlot_medium.mat');
b=load('NullModel/SaveVariablesForPlot20_medium.mat');
n=size(a.dataAvgMeanFeatures,1); % # of features 
m=size(a.dataAvgMeanFeatures,2); % # of xbin
% 
% 
% savename={'clusterVolume' ,'clusterSurfaceArea','clusterSize','clusterSphericity',  ...
%             'clusterPC1','clusterPC2','clusterPC3','clusterPC2_by_PC1','clusterPC3_by_PC1','clusterPC3_by_PC2','radiusOfGyration',...
%             'volumefraction','LorderPar1','LorderPar2','LorderPar3','GorderPar1','GorderPar2','GorderPar3','plane_of_division_growthAxis',...
%             'plane_of_division_to_PC1', 'plane_of_division_to_PC2', 'plane_of_division_to_PC3','highest_mode', 'smallest_mode','clusterRadius'};
% 
% 
%         
%         
%         



        
for i=13   %1:n
    for j=1:m
        v1=a.dataAvgMeanFeatures{i,j};
        ankit(j)=length(v1);
        v2=b.realization_AvgMeanFeatures{i}(j,1);
        ankur(j)=length(v2);
    end
end
        
radius=load('Tibia/dataSave/ParameterForRandomModel_basedon_nuc_medium.dat');
radius=radius(1:481,:);



random=load('NullModel/MakeListColumnarStructurePrediction/DT_E185_nuclei/centroid_and_surface_cells.mat');
real=load('Tibia/MakeListColumnarStructurePrediction/centroid_and_surface_cells.mat');
x_real=load('Tibia/dataSave/AllFeaturesSave_medium.mat');
x_random=load('NullModel/Realization_medium/AllFeaturesSave3.mat');


%index=find(radius>0.7 & radius<1.0);
index=1:length(radius);

% d1=x_real.LOP(index,:);
% d2=x_random.LOP(index,:);


d1=x_real.clusterPC1(index,:);
d2=x_random.clusterPC1(index,:);

[sa,sb1]=sort(abs(d1(:,1)-d2(:,1)),'descend');

sb=sb1(1);
id=index(sb);
oop1=d1(sb,1);
oop2=d2(sb,1);

d1_rg=x_real.rg;
d2_rg=x_random.rg;

d1_phi=x_real.VolumeFraction;
d2_phi=x_random.VolumeFraction;


dir2=strcat('Visual_clusters_random_iteration_','/');
if ~exist([dir2],'dir')
    mkdir([dir2]);
end



for id=1:length(radius)
    h=figure;
    set(gcf, 'PaperSize', [10 3]); %7
    set(gcf, 'PaperPosition', [0 0 10 3]);

    height=radius(id,1);
    id1=x_real.LCC{id};
    id2=x_random.LCC{id};

    subplot(1,2,1)
    for i=1:length(id1)
        v=real.nuc{id1(i)};
        plot3(v(:,1),v(:,2),v(:,3),'.');
        hold on 
    end
    %title(['real(local OOP 1): ', sprintf('%0.2f',oop1)])
    %title(['real(cluster PC 1): ', sprintf('%0.2f',oop1)])
    title(['real:',' R_g = ',sprintf('%0.2f',d1_rg(id)),', \phi = ',sprintf('%0.2f',d1_phi(id))]) 
    xlabel('X'); ylabel('Y'); zlabel('Z')


    subplot(1,2,2)
    for i=1:length(id2)
        v=random.nuc{id2(i)};
        plot3(v(:,1),v(:,2),v(:,3),'.');
        hold on 
    end
    %title(['random(local OOP 1): ', sprintf('%0.2f',oop2)])
    %title(['random(cluster PC 1): ', sprintf('%0.2f',oop2)])
    title(['random:',' R_g = ',sprintf('%0.2f',d2_rg(id)),', \phi = ',sprintf('%0.2f',d2_phi(id))]) 
    xlabel('X'); ylabel('Y'); zlabel('Z')

    saveas(h,[dir2,'/Cluster_',num2str(id),'_X=',sprintf('%0.2f',height),'.png'])
    close all 
end




