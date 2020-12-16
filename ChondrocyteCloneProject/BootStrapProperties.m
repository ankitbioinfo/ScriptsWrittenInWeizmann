clear all 
 
 
% 1-5  CheegerInequality,volume,surface_area,clusterSize,sphericity,
% 6-10 clusterPC1, clusterPC2, clusterPC3,clusterPC2_by_PC1,clusterPC3_by_PC1,
% 11-14 clusterPC3_by_PC2,cycle,diameter,Res3};


% 1: index 
% 2: diameter 
% 3: transitivity 
% 4: clustering coefficient 1
% 5: clustering coefficient 2
% 6: average degree 1
% 7: average degree 2
% 8: efficiency global 
% 9: efficiency local 
% 10: wiener 1 
% 11: wiener 2 
% 12: estrada_index
% 13: subgraph_centrality
% 14: node_connectivity 
% 15: edge_connectivity 
% 16: current_flow_betweenness_centrality_node_1 
% 17: current_flow_betweenness_centrality_node_2
% 18: current_flow_betweenness_centrality_edge_1 
% 19: current_flow_betweenness_centrality_edge_2
% 20: resistance_distance1
% 21: resistance_distance2
% 22: information_centrality1
% 23: information_centrality2
% 24: average degree connectivity 1 
% 25: average degree connectiivty 2 
% 26: robustness 





savename={'clusterVolume' ,'clusterSurfaceArea','clusterSize','clusterSphericity',  ...
            'clusterPC1','clusterPC2','clusterPC3','clusterPC2_by_PC1','clusterPC3_by_PC1','clusterPC3_by_PC2','radiusOfGyration',...
            'volumefraction','LorderPar1','LorderPar2','LorderPar3','GorderPar1','GorderPar2','GorderPar3','plane_of_division_growthAxis',...
            'plane_of_division_to_PC1', 'plane_of_division_to_PC2', 'plane_of_division_to_PC3','highest_mode', 'smallest_mode','clusterRadius',...
          'biaxial_S','biaxial_P','biaxial_D','biaxial_C',  'AngleBetweenClusterPC1AndBone_PD'  };
        
ylabelname={'<cluster volume>', '<cluster surface area>','<cluster size>', '<cluster sphericity>',...
  '<cluster PC1>','<cluster PC2>','<cluster PC3>','<cluster PC2/PC1>','<cluster PC3/PC1>','<cluster PC3/PC2>', '<R_g>',...
   '<volume fraction>', '<orientational OP 1>','<orientational OP 2>','<orientational OP 3>',...
   '<orientational OP 1>','<orientational OP 2>','<orientational OP 3>',  '<angle(plane, P-D axis)>','<angle(plane,PC1)>',...
            '<angle(plane, PC2)>', '<angle(plane,PC3)>', '<largest eigenvalue of Hessian>', '<smallest eigenvalue of Hessian>',...
             '<cluster radius>', '< OOP S>', '< OOP P>', '< OOP D>', '< OOP C>', '<\alpha>'};



dir2=strcat('ClusterFeatures','/');
if ~exist([dir2],'dir')
    mkdir([dir2]);
end


dataFile={'Tibia/','NullModel/'};
mycolor={'ro-','bs-'};

d1=load([dataFile{1},'dataSave/SaveVariablesForPlot_medium.mat']);
d2=load([dataFile{2},'SaveVariablesForPlot20_medium.mat']);

d={d1,d2};

if exist(['SaveBootstrapPvalue.mat'], 'file') == 0  
        bootstrapN=200000;
        for i=1:size(d1.dataAvgMeanFeatures,1) % 19 properites 
            for j=1:size(d1.dataAvgMeanFeatures,2) % xbin data 
                    data1=d1.dataAvgMeanFeatures{i,j};   % individual point 
                    data2=(d2.realization_AvgMeanFeatures{i}(j,:))';
                    %[size(data1), size(data2)]
                    [h,p]=ttest2(data1,data2);
                    ttestpvaluetest{i}(j,1)=p;

                    ObservedDifferences= abs(nanmean(data1) - nanmean(data2));

                    pooldata=[data2; mean(data1)];
                     randomizedDifferences=zeros(bootstrapN,1);
                    for repeat=1:bootstrapN
                        sampling2=randi(length(pooldata),100,1);
                        sampling1=randi(length(pooldata),1,1);
                        randomizedDifferences(repeat,1)= abs(nanmean(pooldata(sampling1)) - nanmean(pooldata(sampling2)));
                    end

                    good=sum( randomizedDifferences >=  ObservedDifferences);
                    empiricalPvalue= good/ bootstrapN;
                    pvaluetest{i}(j,1)=empiricalPvalue;

            end
           % [ttestpvaluetest{i}, pvaluetest{i}]
        end
        save('SaveBootstrapPvalue.mat', 'pvaluetest','ttestpvaluetest')
else
       load(['SaveBootstrapPvalue.mat']);
end


%pvaluetest=ttestpvaluetest;




h=figure;
set(gcf, 'PaperSize', [5 3]); %7
set(gcf, 'PaperPosition', [0 0 5 3]);
for k=1:length(dataFile)
    if k==1
    q(k)=plot(d1.Xbin,d1.Ybin,mycolor{k},'linewidth',1);  
    hold on    
    else
         q(k)=errorbar(d2.Xbin,d2.Ybin,d2.YbinStd,mycolor{k},'linewidth',1); 
    end
    
        
    
end

legname={'cre cluster','background cluster'};

legend(q,legname,'location','northeast')
xlabel('Long axis of bone')
ylabel('cluster density');
saveas(h,[dir2,'clusterDensity.png'])
close all 




for k=1:length(savename)
    h=figure;
    set(gcf, 'PaperSize', [5 3]); %7
    set(gcf, 'PaperPosition', [0 0 5 3]);
    for i=1:length(dataFile)
             data=d{i};
           
                 
                   index=1:length(data.Xbin);
                   
                   q(i)=errorbar(data.Xbin(index),data.AvgMeanFeatures{k}(index),data.AvgStdFeatures{k}(index),mycolor{i},'linewidth',1);
          
                 
                   
            
            
            hold on 
    end
    
   %legend(q,legname,'location','northeast')  %'northeast'
    
    
                   for j=1:length(pvaluetest{k})
                       if pvaluetest{k}(j,1)<0.05
                           plot(d{1}.Xbin(j),d{1}.AvgMeanFeatures{k}(j),'ks','linewidth',2,'markersize',15,'markerfacecolor','none');
                       end
                   end
    
               
    
    if (k>=13)&(k<=15)
        title('local','fontweight','normal')
        legend(q,legname,'location','south')  %'northeast'
        ylim([0,1.01])
    end
    
     if (k>=16)&(k<=18)
        title('global','fontweight','normal')
        legend(q,legname,'location','south')
        ylim([-0.69,1])
     end
     
     if (k>=19)&(k<=22)
	if k==22
        	legend(q,legname,'location','northwest')  %'northeast'
 	else
		legend(q,legname,'location','northeast')  %'northeast'
	end
         ylim([10,80])
     end
    
    
    
    
                  
    
   
    
    xlim([0,1])
    xlabel('Long axis of bone')
    ylabel(ylabelname{k});
    saveas(h,[dir2,savename{k},'.png'])
    close all 
end



% All graph properties 


















% 
% 
%  h=figure;
% XL=0.08;XR=0.02;XGap=0.05;Row=2;
% YT=0.06;YB=0.12;YGap=0.12;Col=4;
% Width=(1-XL-XR-((Col-1)*XGap))/Col;
% Height=(1-YT-YB-((Row-1)*YGap))/Row;
% YPos=1-YT-Height; 
% 
% set(gcf, 'PaperSize', [9 4]); %7
% set(gcf, 'PaperPosition', [0 0 9 4]);
% %mycolor={'b.-','g.-','c.-','r.-','k.-','m.-','b*-','r*-'};
% %mycolor={'r*-','r*-','r*-','r*-','r*-','r*-','r*-','r*-'};
% ylabelname={'3 tris','4 clique', '2 star',  '4 chordcycle', '4 tailed tris', '4 cycle', '3 star','4 path'};
% 
% for i=1:Row
%     XPos=XL;
%     for j=1:Col
%         chro=j+(i-1)*Col;
%             marray=[XPos,YPos,Width,Height];
%             subplot('Position',marray);
%             for k=1:length(dataFile)
%             
%                      load([dataFile{k},'SaveVariablesForPlot.mat'])
% 
%                      %p(k)=plot(Xbin,AvgMeanGraphlet{k},mycolor{k});
%                     plot(Xbin,AvgMeanGraphlet{chro}/max(AvgMeanGraphlet{chro}),mycolor{k},'linewidth',1);
%                     hold on 
%             %xlim([-0.05,1.05])
%             end
%             if i==2
%                 xlabel('Long axis of bone');
%             end
%             if j==1
%               ylabel('frequency/ maximum');
%             end
%            
%             title(ylabelname{chro},'fontweight','normal')
% 
% 	XPos=XPos+Width+XGap;
%     end
%     YPos=YPos-YGap-Height;
% end
% 
% %legend(p,ylabelname,'location','north');
% saveas(h,[dir2,'GraphletFrequency','.png'])
% close all





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

