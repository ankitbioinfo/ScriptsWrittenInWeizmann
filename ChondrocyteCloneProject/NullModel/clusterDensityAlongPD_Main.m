clear all 
 
dirname='MakeListColumnarStructurePrediction/DT_E185_nuclei/';
load([dirname,'centroid_and_surface_cells.mat']);
mu=mean(centroid);
c=centroid-mu;
bonetype=1;
if (bonetype==3)|(bonetype==1)
       c=[-c(:,1),c(:,2:3)];
else
       c=[c(:,1:2),c(:,3)];
end
RZ=min(c(:,1))
HZ=max(c(:,1))
RZ=-562;
HZ=677;
[RZ, HZ]

plot3(c(:,1),c(:,2),c(:,3),'.');
axis image
clear centroid

 





for realization=1:100
    if mod(realization,10)==0
        realization
    end
    %LCC=readClusterFile(strcat('Realization_DT_E185_nuclei/Random_cluster_',num2str(realization),'.dat'));
    [LCC,allid]=readClusterFile(strcat('Realization_medium/Random_cluster_',num2str(realization),'.dat'));
    [realization,length(LCC),length(allid)]

    GVec=averagePC(nuc(allid)); 
    [Xbin,Ybin_iter, AvgMeanFeatures_iter,AvgStdFeatures_iter]=FindAllFeatures(LCC,mu,nuc,c,celvolume,GVec,RZ,HZ,realization);
    realization_Ybin(:,realization)=Ybin_iter;    
    for i=1:length(AvgMeanFeatures_iter)
        realization_AvgMeanFeatures{i}(:,realization)=AvgMeanFeatures_iter{i};
    end
    
end

Ybin=nanmean(realization_Ybin,2);
YbinStd=nanstd(realization_Ybin')';
for i=1:length(realization_AvgMeanFeatures)
    AvgMeanFeatures{i}=nanmean(realization_AvgMeanFeatures{i},2);
    AvgStdFeatures{i}=nanstd(realization_AvgMeanFeatures{i}')';
end

save('SaveVariablesForPlot20_medium.mat', 'AvgMeanFeatures','AvgStdFeatures','Xbin','Ybin','YbinStd','realization_Ybin','realization_AvgMeanFeatures');



       
function [Xbin,Ybin, AvgMeanFeatures,AvgStdFeatures]=FindAllFeatures(LCC,mu,nuc,c,celvolume,GVec,RZ,HZ,realization)

    BoneGrowthAxis=[0.8709, -0.1407, 0.4708]';

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
     
                   
        edges=load(strcat('Realization_medium/Iteration_',num2str(realization),'/Edges_',num2str(j),'.dat'));
        no_of_nodes=length(unique(LCC{j}));
        [coordNumber(j,1), diameter(j,1),stepsize]=coordination_number(c,edges);
        % this is basically R_g of random walk 
        averageStepSize(j,1) =   sqrt( 1/6*(no_of_nodes*(stepsize^2)));
       
        e=normal_modes(c,LCC{j});
        %[length(LCC{j}),nnz(e)]
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

sphericity= (pi^(1/3))*((6*volume).^(2/3)) ./ surface_area;
clusterPC2_by_PC1=clusterPC2./clusterPC1;
clusterPC3_by_PC1=clusterPC3./clusterPC1;
clusterPC3_by_PC2=clusterPC3./clusterPC2;

VolumeFraction=  combinedCellVolumeInCluster./volume;


save(['Realization_medium/AllFeaturesSave',num2str(realization),'.mat'],'volume','surface_area','clusterSize','sphericity','clusterPC1', 'clusterPC2', 'clusterPC3','clusterPC2_by_PC1','clusterPC3_by_PC1','clusterPC3_by_PC2',...
    'rg','VolumeFraction','LOP','GOP','plane_of_division_to_bone','plane_of_division_to_cell','highest_mode', 'diameter',...
    'coordNumber','smallest_mode','clusterRadius', 'averageStepSize','LCC','LocalBiaxial','AngleBetweenClusterPC1AndBone_PD');






Features={volume,surface_area,clusterSize,sphericity,clusterPC1, clusterPC2, clusterPC3,clusterPC2_by_PC1,clusterPC3_by_PC1,clusterPC3_by_PC2,...
    rg,VolumeFraction,LOP(:,1),LOP(:,2),LOP(:,3),GOP(:,1),GOP(:,2),GOP(:,3),plane_of_division_to_bone,plane_of_division_to_cell(:,1),...
    plane_of_division_to_cell(:,2),plane_of_division_to_cell(:,3),highest_mode, smallest_mode,clusterRadius,...
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


radius=load('../Tibia/dataSave/ParameterForRandomModel_basedon_nuc_medium.dat');

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
      end
end

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
               
                



function rg=radiusOfGyration(centroid,nodes)
        c=centroid(nodes,:);
        
        center=mean(c,1);
        N= size(c,1);
        
        for i=1:N
            r(i,:)=sum((c(i,:) - center).^2);
        end
        rg=sqrt(sum(r)/N);
        
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
 
 
 function [LCC,allCellId]=readClusterFile(name) 
        fid = fopen(name,'rt');
        tline = fgetl(fid);
        count=1;
        allCellId=[];
         while ischar(tline)
            line= split(tline);
             if length(line)>3
                for j=1:length(line)-1
                    LCC{count}(j)=str2num(line{j});
                    allCellId=[allCellId;LCC{count}(j)];
                end
            end     
            
            tline = fgetl(fid);
            count=count+1;
        end
        fclose(fid);
 end
