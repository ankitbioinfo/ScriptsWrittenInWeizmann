clear all 
 
dirname='MakeListColumnarStructurePrediction/';
load([dirname,'centroid_and_surface_cells.mat']);

load('MakeListColumnarStructurePrediction/alledges.mat');

% %nnetwork_edges_large: [5911?2 double]
%     network_edges_medium: [3589?2 double]
%      network_edges_small: [1494?2 double]

%mycutoffEdges=network_edges_large; [LCC,LCC1]=readClusterFile([dirname,'ClusterLargerCutoff.dat']); outputcutoff='larger';
mycutoffEdges=network_edges_medium; [LCC,LCC1]=readClusterFile([dirname,'ClusterMediumCutoff.dat']); outputcutoff='medium'; 
%mycutoffEdges=network_edges_small; [LCC,LCC1]=readClusterFile([dirname,'ClusterSmallerCutoff.dat']); outputcutoff='smaller'; 


 



bone_mu=mean(centroid);
c=centroid-bone_mu;
bonetype=1;
if (bonetype==3)|(bonetype==1)
       c=[-c(:,1),c(:,2:3)];
       RotateFactor=[-1 0 0;0 1 0 ;0 0 1];
else 
       c=[c(:,1:2),c(:,3)];
end
RZ=min(c(:,1))
HZ=max(c(:,1))

RZ=-562;
HZ=677;
[RZ, HZ]

plot3(c(:,1),c(:,2),c(:,3),'.');
hold on 
axis image
%clear centroid

GVec=averagePC(nuc)  

[vec,val]=eig(cov(c));
mu=mean(c);d = sqrt(diag(val));
factor=3; quiver3(mu(1),mu(2),mu(3),factor*vec(1,1),factor*vec(2,1),factor*vec(3,1),d(1),'r','LineWidth',10);
factor=3; quiver3(mu(1),mu(2),mu(3),factor*vec(1,2),factor*vec(2,2),factor*vec(3,2),d(2),'k','LineWidth',5);
factor=3; quiver3(mu(1),mu(2),mu(3),factor*vec(1,3),factor*vec(2,3),factor*vec(3,3),d(3),'m','LineWidth',2);
 
factor=-3; quiver3(mu(1),mu(2),mu(3),factor*vec(1,1),factor*vec(2,1),factor*vec(3,1),d(1),'r','LineWidth',10);
factor=-3; quiver3(mu(1),mu(2),mu(3),factor*vec(1,2),factor*vec(2,2),factor*vec(3,2),d(2),'k','LineWidth',5);
factor=-3; quiver3(mu(1),mu(2),mu(3),factor*vec(1,3),factor*vec(2,3),factor*vec(3,3),d(3),'m','LineWidth',2);

BoneGrowthAxis=vec(:,3)
%[0.8709, -0.1407, 0.4708]';


root_dir = '';
addpath(fullfile(root_dir, 'utils'));


dir2=strcat('dataSave','/');
if ~exist([dir2],'dir')
    mkdir([dir2]);
end



 close all 
 
% for ci=330
%     id=LCC{ci};
%     
%     if length(id)>3
%     
%             mycell=[];
%             for j=1:length(id)
%                 mycell=[mycell;nuc{id(j)}];
%             end
% 
%             K = convhull(mycell); bnd_pnts = mycell(K,:); 
% 
%             pos=centroid(id,:);
%             h=figure;
%             [vornb,vorvx,T,Aaug,Baug] = polybnd_voronoi(pos,bnd_pnts); 
%             for i = 1:size(vorvx,2)
%                 col(i,:)= rand(1,3);
%             end
% 
%              for i = 1:size(pos,1)
%                     K = convhulln(vorvx{i});
%                     trisurf(K,vorvx{i}(:,1),vorvx{i}(:,2),vorvx{i}(:,3),'FaceColor',col(i,:),'FaceAlpha',1,'EdgeAlpha',0)
%                     hold on;
%              end
%              %scatter3(pos(:,1),pos(:,2),pos(:,3),'Marker','o','MarkerFaceColor',[0 .75 .75], 'MarkerEdgeColor','k');
%              axis('equal')
%              saveas(h,['polybound/cluster',num2str(ci),'.fig'])
%              saveas(h,['polybound/cluster',num2str(ci),'.png'])
%              close all 
%     end
% end

% plot3(vorvx{1}(:,1),vorvx{1}(:,2),vorvx{1}(:,3),'r.'); hold on 
% plot3(vorvx{2}(:,1),vorvx{2}(:,2),vorvx{2}(:,3),'b.'); hold on
% plot3(vorvx{3}(:,1),vorvx{3}(:,2),vorvx{3}(:,3),'k.'); hold on
% plot3(vorvx{4}(:,1),vorvx{4}(:,2),vorvx{4}(:,3),'g.'); hold on
% plot3(pos(:,1),pos(:,2),pos(:,3),'ko','markersize',30)

% 
% for ci=330
%         id=LCC{ci};
%         N=c(id,:);     
%         % calculating the triangulation and the volume of each triangle:
%         TRI = delaunay(N(:,1), N(:,2), N(:,3));
%         clear neighbor
%         for i = 1 : size(N,1)
%             temp=[];
%             for j=1:size(TRI,1)
%                 for k=1:size(TRI,2)
%                     if TRI(j,k)==i
%                         temp=[temp,TRI(j,:)];
%                     end
%                 end
%             end
%             neighbor=setdiff(unique(temp),i);
%             for j=1:length(neighbor)
%                 centroidVector= N([i,neighbor(j)],:);
%                 normCentroidVector=centroidVector/norm(centroidVector);
%                 vec=find_perp(normcentroidVector); 
%             end 
%         end
% end

       



%cluster properties calculate 
 for j=1:length(LCC)
     ClusterShape=[];
     mergeVol=0;
     for i=1:length(LCC{j})
            ClusterShape=[  ClusterShape;nuc{LCC{j}(i)}];
            mergeVol=mergeVol+celvolume(LCC{j}(i));
     end
     combinedCellVolumeInCluster(j,1)=mergeVol;
     [K,V]=convhull(ClusterShape);
     volume(j,1)=V;
     surface_area(j,1)=CalculateSurfaceArea(K,ClusterShape);
     clusterSize(j,1)=length(LCC{j});
     [clusterPCvector,~,latent]=pca(ClusterShape);
     
     [ccenter,clusterRadius(j,1)] = sphereFit(ClusterShape);
     clusterCenter(j,:)=ccenter-mu;
         
     
     AngleBetweenClusterPC1AndBone_PD(j,1)=180/pi*oangle(clusterPCvector(:,1),BoneGrowthAxis);
     
     clusterPC1(j,1)=latent(1);
     clusterPC2(j,1)=latent(2);
     clusterPC3(j,1)=latent(3);
     
     mainVec=averagePC(nuc(LCC{j}));
     
     clear LNS 
     clear GNS
     clear theta 
     clear biaxial
     
     for i=1:length(LCC{j})
          id=LCC{j}(i);
          [PC,~,~]= pca(nuc{id});
          %oangle(PC(:,1),mainVec(:,1))*180/pi
          LNS(i,1)= 3*(cos(  oangle(PC(:,1),mainVec(:,1))   )^2) -1;
          LNS(i,2)= 3*(cos(  oangle(PC(:,2),mainVec(:,2))  )^2) -1;
          LNS(i,3)= 3*(cos(  oangle(PC(:,3),mainVec(:,3))  )^2)-1;
          
          GNS(i,1)= 3*(cos(  oangle(PC(:,1),GVec(:,1))  )^2) -1;
          GNS(i,2)= 3*(cos(  oangle(PC(:,2),GVec(:,2))  )^2) -1;
          GNS(i,3)= 3*(cos(  oangle(PC(:,3),GVec(:,3))  )^2)-1;     
          
          
          biaxial(i,:)=UniaxialAndBiaxialOOP( PC, mainVec);
          
     end
     
     LOP(j,1)=0.5*(mean(LNS(:,1))); 
     LOP(j,2)=0.5*(mean(LNS(:,2))); 
     LOP(j,3)=0.5*(mean(LNS(:,3))); 
     
     GOP(j,1)=0.5*(mean(GNS(:,1))); 
     GOP(j,2)=0.5*(mean(GNS(:,2))); 
     GOP(j,3)=0.5*(mean(GNS(:,3))); 
     
     LocalBiaxial(j,:)=mean(biaxial,1);
     rg(j,1)=radiusOfGyration(c,LCC{j});
     edges=search_edges(mycutoffEdges,LCC{j});

     
        unique(edges(:));
        no_of_nodes=length(unique(LCC{j}));
     
        [coordNumber(j,1), diameter(j,1),stepsize]=coordination_number(c,edges);
        % this is basically R_g of random walk 
        averageStepSize(j,1) =   sqrt( 1/6*(no_of_nodes*(stepsize^2)));
        
       
        e=normal_modes(c,LCC{j});
        % [length(LCC{j}),nnz(e)]
        highest_mode(j,1)=max(e);
        smallest_mode(j,1)=min(e);
        
        id=LCC{j};
        theta1=[];
        theta2=[];
        theta3=[];
        theta4=[];
        for i=1:size(edges,1)
            flag=[0,0];
            for k=1:length(id)
                if id(k)==edges(i,1)
                    flag(1)=1;
                end
                if id(k)==edges(i,2)
                    flag(2)=1;
                end
            end
            if sum(flag)==2
                 P1= c(edges(i,1),:);  P2= c(edges(i,2),:);
                 P= nuc{edges(i,1)}(1,:);
                 R = cross(P1-P2, P1-P);
                 S = cross(R, P1-P2);
                 alphap=1; betap=1;
                 plane= mean([P1;P2])' + alphap*R' + betap*S';
                 %normCentroidVector=centroidVector/norm(centroidVector);
                 %vec=find_perp(normCentroidVector); vec=vec/norm(vec); [normCentroidVector', vec']
                 %[sum((P1-P2).*S), sum((P1-P2).*R)]
                 
                 
                
                 norm_plane=norm(plane);
                 
                 [PC1,~,~]= pca(nuc{edges(i,1)});  [PC2,~,~]= pca(nuc{edges(i,2)});
                 
                 theta(1)=asin( abs(sum(plane.*PC1(:,1))/norm_plane))*180/pi;
                 theta(2)=asin( abs(sum(plane.*PC2(:,1))/norm_plane))*180/pi;
                 theta1 =[theta1;mean(theta)];
                 
                 theta(1)=asin( abs(sum(plane.*PC1(:,2))/norm_plane))*180/pi;
                 theta(2)=asin( abs(sum(plane.*PC2(:,2))/norm_plane))*180/pi;
                 theta2 =[theta2;mean(theta)];
           
                 theta(1)=asin( abs(sum(plane.*PC1(:,3))/norm_plane))*180/pi;
                 theta(2)=asin( abs(sum(plane.*PC2(:,3))/norm_plane))*180/pi;
                 theta3 = [theta3;mean(theta)];
                 
                 
                 theta4=[theta4;asin( abs(sum(plane.*BoneGrowthAxis)/norm_plane))*180/pi];
                 

            end
        end
        
            
        
        plane_of_division_to_cell(j,:)=[mean(theta1),mean(theta2),mean(theta3)];
        plane_of_division_to_bone(j,1)=mean(theta4);
   
 end
 
 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 for j=1:length(LCC1)
     ClusterShape=[];
     mergeVol=0;
     for i=1:length(LCC1{j})
            ClusterShape=[  ClusterShape;nuc{LCC1{j}(i)}];
            mergeVol=mergeVol+celvolume(LCC1{j}(i));
     end
     warning off;
     [ccenter,clusterRadius1(j,1)] = sphereFit(ClusterShape);
     clusterCenter1(j,:)=(ccenter-bone_mu)*RotateFactor;
     [ccenter2(j,:),tempR] = sphereFit(c(LCC1{j},:));
     if tempR>0
        temp_cR(j,1)=min([tempR,clusterRadius1(j,1)]);
     else
        temp_cR(j,1)=clusterRadius1(j,1);
     end
     
 end
 

binsize=51;zone=linspace(RZ,HZ,binsize);Xbin=linspace(0,1,binsize-1);
for i=1:length(zone)-1
    for j=1:length(LCC1)
         cluster= clusterCenter1(j,:);
             if ((cluster(1) >zone(i)) &  (cluster(1) <=zone(i+1)))
                  clusterPositionAndRadius(j,:) = [Xbin(i),  clusterRadius1(j)];
             end
    end
end
dlmwrite([dir2,'ParameterForRandomModel_basedon_nuc_',outputcutoff,'.dat'], clusterPositionAndRadius,'\t');

% for i=1:length(zone)-1
%     for j=1:length(LCC1)
%          cluster=  ccenter2(j,:);
%              if ((cluster(1) >zone(i)) &  (cluster(1) <=zone(i+1)))
%                   clusterPositionAndRadius(j,:) = [Xbin(i),  temp_cR(j)];
%              end
%     end
% end
% dlmwrite('ParameterForRandomModel_basedon_centroid.dat', clusterPositionAndRadius,'\t');
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 

sphericity= (pi^(1/3))*((6*volume).^(2/3)) ./ surface_area;
clusterPC2_by_PC1=clusterPC2./clusterPC1;
clusterPC3_by_PC1=clusterPC3./clusterPC1;
clusterPC3_by_PC2=clusterPC3./clusterPC2;

VolumeFraction=  combinedCellVolumeInCluster./volume;


save([dir2,'AllFeaturesSave_',outputcutoff,'.mat'],'volume','surface_area','clusterSize','sphericity','clusterPC1', 'clusterPC2', 'clusterPC3','clusterPC2_by_PC1','clusterPC3_by_PC1','clusterPC3_by_PC2',...
    'rg','VolumeFraction','LOP','GOP','plane_of_division_to_bone','plane_of_division_to_cell','highest_mode', 'diameter',...
    'coordNumber','smallest_mode','clusterRadius', 'averageStepSize','LCC','LocalBiaxial','AngleBetweenClusterPC1AndBone_PD');



Features={volume,surface_area,clusterSize,sphericity,clusterPC1, clusterPC2, clusterPC3,clusterPC2_by_PC1,clusterPC3_by_PC1,clusterPC3_by_PC2,...
    rg,VolumeFraction,LOP(:,1),LOP(:,2),LOP(:,3),GOP(:,1),GOP(:,2),GOP(:,3), plane_of_division_to_bone,plane_of_division_to_cell(:,1),...
    plane_of_division_to_cell(:,2),plane_of_division_to_cell(:,3), highest_mode, smallest_mode,clusterRadius, ...,
    LocalBiaxial(:,1), LocalBiaxial(:,2), LocalBiaxial(:,3), LocalBiaxial(:,4), AngleBetweenClusterPC1AndBone_PD,...
    coordNumber,diameter,averageStepSize};
% savename={'CheegerInequality','clusterVolume' ,'clusterSurfaceArea','clusterSize','clusterSphericity',  ...
%             'clusterPC1','clusterPC2','clusterPC3','clusterPC2_by_PC1','clusterPC3_by_PC1','clusterPC3_by_PC2',...
%             'conductance','radiusOfGyration','volumefraction','orderParameter1','orderParameter2','orderParameter3','GOP1','GOP2','GOP3',...
%             'InteractionEnergy1','InteractionEnergy2','InteractionEnergy3'};
% savename{43}='Cluster_Radius';        
% ylabelname={'Cheeger''s Inequality bounds of graphlets','<cluster volume>', '<cluster surface area>','<cluster size>', '<cluster sphericity>',...
%   '<cluster PC1>','<cluster PC2>','<cluster PC3>','<cluster PC2/PC1>','<cluster PC3/PC1>','<cluster PC3/PC2>',...
%  '<conductance>','<r_g>', '<volume fraction>', '<orientational OP1>','<orientational OP2>','<orientational OP3>',...
%  '<orientational OP1>','<orientational OP2>','<orientational OP3>','<E_1>','<E_2>','<E_3>' };
% ylabelname{43}='<cluster radius>';



binsize=21;
zone=linspace(RZ,HZ,binsize);
Xbin=linspace(0,1,binsize-1);
Ybin=zeros(binsize-1,1); 

radius=clusterPositionAndRadius;

for i=1:length(zone)-1
    cellsIdInXInterval=[];
%     for j=1:length(LCC)
%          cluster=mean(c(LCC{j},:));
%              if ((cluster(1) >zone(i)) &  (cluster(1) <=zone(i+1)))
%                  Ybin(i)=Ybin(i)+1;
%                  cellsIdInXInterval=[cellsIdInXInterval,j];
%              end
%     end
    
      for j=1:length(LCC)
           if ((radius(j,1) >Xbin(i)) &  (radius(j,1) <=Xbin(i+1)))
                 Ybin(i)=Ybin(i)+1;
                 cellsIdInXInterval=[cellsIdInXInterval,j];
           end
       end
    
    
    
      for k=1:length(Features)          
                AvgMeanFeatures{k}(i,1)=nanmean(Features{k}(  cellsIdInXInterval,1));
                AvgStdFeatures{k}(i,1)=nanstd(Features{k}(  cellsIdInXInterval,1));
                dataAvgMeanFeatures{k,i}=Features{k}(  cellsIdInXInterval,1);
      end
end


save([dir2,'SaveVariablesForPlot_',outputcutoff,'.mat'],   'AvgMeanFeatures','AvgStdFeatures','Xbin','Ybin','dataAvgMeanFeatures');





function [avgdeg,diameter,averageStepSize]=coordination_number(centroids,edges)

             id=unique(edges(:));
             for j=1:length(id)
                 new2old(j)=id(j);
                 old2new(id(j))=j;
             end
               bondDist=[];
               for i=1:size(edges,1)
                 newedges(i,1)=old2new(edges(i,1));
                 newedges(i,2)=old2new(edges(i,2));
                 bondDist=[bondDist;pdist(centroids(edges(i,:),:))];
               end
               
            averageStepSize=mean(bondDist);


            G=graph(newedges(:,1),newedges(:,2));
            deg=G.degree();
            avgdeg=mean(deg);
            %[length(nodes),length(deg)]
            
            d=distances(G);
            diameter=max(d(:));
end




function e=normal_modes(centroids,nodes)

%          id=unique(edges(:));
%          for j=1:length(id)
%              new2old(j)=id(j);
%              old2new(id(j))=j;
%          end
%          N=length(id);
%          M=zeros(3*N,3*N);
%          for i=1:size(edges,1)
%              r1=centroids(edges(i,1),:);
%              r2=centroids(edges(i,2),:);
%              
%              p1=old2new(edges(i,1));
%              p2=old2new(edges(i,2));
%              
%              q1=3*p1-2; q2=3*p2-2;
%              
%              dist=pdist([r1;r2])^2;
%              
%              M(q1:q1+2,q2:q2+2)=-(1/dist)*   ( (r1-r2)'*(r1-r2) );
%              M(q2:q2+2,q1:q1+2)=-(1/dist)*   ( (r2-r1)'*(r2-r1) );
%          end
%     
         


         id=unique(nodes);
         for j=1:length(id)
             new2old(j)=id(j);
             old2new(id(j))=j;
         end
         N=length(id);
         M=zeros(3*N,3*N);
         for i=1:N
             r1=centroids(id(i),:);
             for j=i+1:N
                     r2=centroids(id(j),:);
             
                     p1=old2new(id(i));
                     p2=old2new(id(j));
             
                     q1=3*p1-2; q2=3*p2-2;

                     dist=pdist([r1;r2])^2;

                     M(q1:q1+2,q2:q2+2)=-(1/dist)*   ( (r1-r2)'*(r1-r2) );
                     M(q2:q2+2,q1:q1+2)=-(1/dist)*   ( (r2-r1)'*(r2-r1) );
             end
         end
         
         
         
         
         %sum(diag(M))
         for i=1:size(M,1)
             M(i,i)=sum(M(i,:));
         end
         
         e=eig(M);
         
end
                



function biaxial=UniaxialAndBiaxialOOP( MoleculePC, LaboratoryPC)
            
         n=MoleculePC(:,1); l=MoleculePC(:,2); m=MoleculePC(:,3); 
         w=LaboratoryPC(:,1);  u=LaboratoryPC(:,2);  v=LaboratoryPC(:,3); 
         
         S = 0.5*  (  (3*(dot(n,w)/(norm(n)*norm(w)) )^2) -1) ; 

         P = 3/2*  (  (dot(n,v)/(norm(n)*norm(v)) )^2 -  (dot(n,u)/(norm(n)*norm(u)) )^2        ) ; 
         
         D = 3/2*  (  (dot(l,w)/(norm(l)*norm(w)) )^2 -  (dot(m,w)/(norm(m)*norm(w)) )^2        ) ;
         
         C=  3/2*  (  (dot(l,v)/(norm(l)*norm(v)) )^2 -  (dot(l,u)/(norm(l)*norm(u)) )^2    +  ...
                      (dot(m,u)/(norm(m)*norm(u)) )^2 -  (dot(m,v)/(norm(m)*norm(v)) )^2                 ) ;
                  
                  
         biaxial=[S,P,D,C];         
end




function rg=radiusOfGyration(centroid,nodes)
        c=centroid(nodes,:);
        N= size(c,1);

        
        center=mean(c,1);
        for i=1:N
            r(i,:)=sum((c(i,:) - center).^2);
        end
        rg=sqrt(sum(r)/N);

%           rg=[]; 
%           for i=1:N
%               for j=1:N
%                   rg=[rg;norm((c(i,:) - c(j,:)))^2  ];
%               end
%           end
%           
%           rg= sqrt(  1/(2*N*N) * sum(rg));
          
          
        
end
            


 function GVec=averagePC(nuc)        
            %length(nuc)
            for i=1:length(nuc)
                [pc,~,~]=pca(nuc{i});
                val1(i,:)=pc(:,1)';
                val2(i,:)=pc(:,2)';
                val3(i,:)=pc(:,3)';
            end

            pc1=pca( [val1;-val1]);
            pc2=pca( [val2;-val2]);
            pc3=pca( [val3;-val3]);

            GVec(:,1)=pc1(:,1);
            GVec(:,2)=pc2(:,1);
            GVec(:,3)=pc3(:,1);
 end




function Area=CalculateSurfaceArea(K,v)
          Area=0;
          for surfind=1:size(K,1)
              pointA=v(K(surfind,1),:);
              pointB=v(K(surfind,2),:);
              pointC=v(K(surfind,3),:);
              TriangleVertex=[pointA; pointB; pointC];
              Ts=pdist(TriangleVertex,'euclidean');
              p=sum(Ts)/2;
              Area=Area+sqrt(p*(p-Ts(1))*(p-Ts(2))*(p-Ts(3)));
          end
          
end


function graphEnergy=calculateEnergy(node,cent)
    
    A=zeros(length(node));
    D=zeros(length(node));
    centroid=cent(node,:);
    for i=1:size(centroid,1)
        for j=i+1:size(centroid,1)
            dist=pdist(centroid([i,j],:));
            A(i,j)= 1/dist;
            A(j,i)= 1/dist;
        end
    end


    for i=1:size(A,1)
        D(i,i)=sum(A(i,:));
    end

    normalizedAdjacency=D^(-1/2)*A*(D^(-1/2));
    L=D-A;
    normalizedLaplacian=D^(-1/2)*L*(D^(-1/2));


    E1=eig(normalizedAdjacency);
    E2=eig(normalizedLaplacian);

    %E2([1:3,end])

    energyAdjacency=sum(abs(E1));

    %totalWeight=sum(sum(A));
    totalWeight=sum(E2);

    energyLaplacian=sum( abs( E2 - (totalWeight/size(A,1))));

    secondSmallestLaplacianEigenvalue=E2(2);

    graphEnergy=[ secondSmallestLaplacianEigenvalue,energyAdjacency,energyLaplacian,min(E1),max(E1)];


end


 function  angle= oangle(u,v)
      angle=atan2(norm(cross(u,v)),dot(u,v));
      if angle>(pi/2)
          angle=pi-angle;
      end
 end
 
 
 function [LCC,LCC1]=readClusterFile(name) 
        fid = fopen(name,'rt');
        tline = fgetl(fid);
        count=1;
         while ischar(tline)
            line= split(tline,' ');
             if length(line)>3
                for j=1:length(line)-1
                    LCC{count}(j)=str2num(line{j});
                end
             end    
            
             if length(line)>2
                for j=1:length(line)-1
                    LCC1{count}(j)=str2num(line{j});
                end
             end    
             
            
            tline = fgetl(fid);
            count=count+1;
        end
        fclose(fid);
 end

 
                         
function myedges=search_edges(edges,vertices) 
    [~,ia]=unique(edges,'rows');
    edges=edges(ia,:);
	count=1;
    myedges=[];
	for j=1: size(edges,1)
		flag1=0;
		flag2=0;
		for i=1:length(vertices)
			if (edges(j,1)==vertices(i))
				flag1=1;
            end
			if (edges(j,2)==vertices(i))
				flag2=1;
            end
            if (flag1+flag2)==2
                myedges(count,:)=edges(j,:);
                count=count+1;
            end
        end
    end
    
    [~,ia]=unique(myedges,'rows');
    myedges=myedges(ia,:);
    
 end
