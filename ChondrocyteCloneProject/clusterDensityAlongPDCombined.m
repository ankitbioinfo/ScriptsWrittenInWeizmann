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
          'biaxial_S','biaxial_P','biaxial_D','biaxial_C',  'AngleBetweenClusterPC1AndBone_PD' ,'coordination_number','diameter' };
        
ylabelname={'<cluster volume>', '<cluster surface area>','<cluster size>', '<cluster sphericity>',...
  '<cluster PC1>','<cluster PC2>','<cluster PC3>','<cluster PC2/PC1>','<cluster PC3/PC1>','<cluster PC3/PC2>', '<R_g>',...
   '<volume fraction>', '<Local OOP 1>','<Local OOP 2>','Local <OOP 3>',...
   '<Global OOP 1>','<Global OOP 2>','<Global OOP 3>',  '<angle(plane, P-D axis)>','<angle(plane,PC1)>',...
            '<angle(plane, PC2)>', '<angle(plane,PC3)>', '<largest eigenvalue of Hessian>', '<smallest eigenvalue of Hessian>',...
             '<cluster radius>', '< OOP S>', '< OOP P>', '< OOP D>', '< OOP C>', '<\alpha>','<coordination number>','<diameter>'};

ylabelname={'<volume>', '<surface area>','<size>', '<sphericity>',...
  '<PC1 coeff>','<PC2 coeff>','<PC3 coeff>','<PC2/PC1 coeff>','<PC3/PC1 coeff>','<PC3/PC2 coeff>', '<R_g>',...
   '<\phi>', '<Local OOP 1>','<Local OOP 2>','Local <OOP 3>',...
   '<Global OOP 1>','<Global OOP 2>','<Global OOP 3>',  '<angle(plane, P-D axis)>','<angle(plane,PC1)>',...
            '<angle(plane, PC2)>', '<angle(plane,PC3)>', '<largest eigenvalue of Hessian>', '<smallest eigenvalue of Hessian>',...
             '<cluster radius>', '< OOP S>', '< OOP P>', '< OOP D>', '< OOP C>', '<\alpha>','<coordination number>','<diameter>'};
         

dir2=strcat('ClusterFeatures','/');
if ~exist([dir2],'dir')
    mkdir([dir2]);
end


dataFile={'Tibia/','NullModel/'};
mycolor={'ro-','bs-'};

d1=load([dataFile{1},'dataSave/SaveVariablesForPlot_medium.mat']);
d2=load([dataFile{2},'SaveVariablesForPlot20_medium.mat']);

d={d1,d2};


% for i=1:size(d1.dataAvgMeanFeatures,1) % 19 properites 
%     for j=1:size(d1.dataAvgMeanFeatures,2) % xbin data 
%             data1=d1.dataAvgMeanFeatures{i,j};   % individual point 
%             data2=(d2.realization_AvgMeanFeatures{i}(j,:))';
%             %[size(data1), size(data2)]
%             [h,p]=ttest2(data1,data2);
%             pvaluetest{i}(j,1)=p;
%             
%     end
% end

radius=load('Tibia/dataSave/ParameterForRandomModel_basedon_nuc_medium.dat');



% coordination number 
h=figure;
XL=0.09;XR=0.02;XGap=0.06;Row=2;
YT=0.06;YB=0.1;YGap=0.09;Col=2;
Width=(1-XL-XR-((Col-1)*XGap))/Col;
Height=(1-YT-YB-((Row-1)*YGap))/Row;
YPos=1-YT-Height; 
set(gcf, 'PaperSize', [6 4]); %7
set(gcf, 'PaperPosition', [0 0 6 4]);
zones={'RZ','PZ','PHZ','HZ'};
for i=1:Row
    XPos=XL;
    for j=1:Col
        chro=j+(i-1)*Col;
        marray=[XPos,YPos,Width,Height];
        subplot('Position',marray);
        FunctionOfCellSize(d1.dataAvgMeanFeatures, d1.Xbin,1,3,31,1,chro);
        %title(zones{chro},'fontweight','normal')
        if j==1
        ylabel(ylabelname{31})
        end
        if i==Row
        xlabel('cluster size (# of cells)')
        else
            set(gca,'xticklabel',[])
        end
        %ylim([5,28])
        xlim([0,40])
        XPos=XPos+Width+XGap;
    end
    YPos=YPos-YGap-Height;
end
saveas(h,[dir2,'CoordinationNumber.png'])
close all 






%Radius of gyration from Random walk model 
h=figure;
XL=0.09;XR=0.02;XGap=0.06;Row=2;
YT=0.06;YB=0.1;YGap=0.09;Col=2;
Width=(1-XL-XR-((Col-1)*XGap))/Col;
Height=(1-YT-YB-((Row-1)*YGap))/Row;
YPos=1-YT-Height; 
set(gcf, 'PaperSize', [6 4]); %7
set(gcf, 'PaperPosition', [0 0 6 4]);
zones={'RZ','PZ','PHZ','HZ'};
for i=1:Row
    XPos=XL;
    for j=1:Col
        chro=j+(i-1)*Col;
        marray=[XPos,YPos,Width,Height];
        subplot('Position',marray);

        p(1)=FunctionOfCellSize(d1.dataAvgMeanFeatures, d1.Xbin,1,3,11,0,chro);
        hold on 
        p(2)=FunctionOfCellSize(d1.dataAvgMeanFeatures, d1.Xbin,1,3,33,0,chro);
        if chro==1
        legend(p,'cre cluster', 'random walk theory','location','northwest')
        end
        title(zones{chro},'fontweight','normal')
        if j==1
        ylabel(ylabelname{11})
        end
        if i==Row
        xlabel('cluster size (# of cells)')
        else
            set(gca,'xticklabel',[])
        end
        ylim([5,28])
        xlim([0,40])
        XPos=XPos+Width+XGap;
    end
    YPos=YPos-YGap-Height;
end
saveas(h,[dir2,'RandomWalkModel.png'])
close all 


% Property is function of cluster size 
h=figure;
XL=0.05;XR=0.02;XGap=0.06;Row=4;
YT=0.06;YB=0.1;YGap=0.09;Col=5;
Width=(1-XL-XR-((Col-1)*XGap))/Col;
Height=(1-YT-YB-((Row-1)*YGap))/Row;
YPos=1-YT-Height; 
set(gcf, 'PaperSize', [12 7]); %7
set(gcf, 'PaperPosition', [0 0 12 7]);
Prop=[1,2,4,5:18,30,31,32];
for i=1:Row
    XPos=XL;
    for j=1:Col
        chro=j+(i-1)*Col;
        marray=[XPos,YPos,Width,Height];
        
        if chro<=length(Prop)
        subplot('Position',marray);

        FunctionOfCellSize(d1.dataAvgMeanFeatures, d1.Xbin,1,3,Prop(chro),1,2);
        ylabel(ylabelname{Prop(chro)})
        
        if i==Row 
            xlabel('cluster size (# of cells)')
        end

        end
        XPos=XPos+Width+XGap;
    end
    YPos=YPos-YGap-Height;
end
saveas(h,[dir2,'VariousMetricFunctionOfClusterSize.png'])
close all 








% coordination number 
% real_size= mydata{3,j};
% real_rg=mydata{11,j};
% avg_deg=mydata{31,j};
% phi=mydata{12,j};
 h=figure; 
CoordinatinNumberCalculation(d1.dataAvgMeanFeatures, d1.Xbin,1,12,31)
xlabel('Volume fraction')
ylabel('Coordination number');
saveas(h,[dir2,'phi_vs_avgCoordinatinNumber.png'])
close all 


 




 
 h=figure;
 plot(d1.Xbin, d1.AvgMeanFeatures{31}, 'r.-');
 hold on 
 errorbar(d1.Xbin,d1.AvgMeanFeatures{31},d1.AvgStdFeatures{31},'r','linewidth',1);
 ylabel('Coordination number');     
 xlabel('Long axis of bone')
 saveas(h,[dir2,'SpatialProfile_AvgCoordinatinNumber.png'])
 close all 

 
 
 

 
 


% fractal like pattern 
 h=figure;
XL=0.09;XR=0.02;XGap=0.05;Row=1;
YT=0.06;YB=0.15;YGap=0.12;Col=2;
Width=(1-XL-XR-((Col-1)*XGap))/Col;
Height=(1-YT-YB-((Row-1)*YGap))/Row;
YPos=1-YT-Height; 
set(gcf, 'PaperSize', [7 3]); %7
set(gcf, 'PaperPosition', [0 0 7 3]);
for i=1:Row
    XPos=XL;
    for j=1:Col
        chro=j+(i-1)*Col;
        marray=[XPos,YPos,Width,Height];
        subplot('Position',marray);
        
        if chro==1
               fractalCalculation( d1.dataAvgMeanFeatures, d1.Xbin , 1 )
                title('cre clusters') 
               ylabel('log(R_g)'); 
        end
        
        if chro==2
               data=[];
               Xbin=[];
               for realization=1:50
                   d3=load([dataFile{2},'Realization_medium/AllFeaturesSave',num2str(realization),'.mat']);
                   data=[data; [d3.clusterSize d3.rg]];  
                   Xbin=[Xbin; radius(1:length(d3.rg),1)];
               end
               
               
               
                fractalCalculation( data, Xbin, 2  )
                title('background clusters') 
        end

        
        ylim([1.3,3.3])
        xlim([1,4])
        
        XPos=XPos+Width+XGap;
    end
    YPos=YPos-YGap-Height;
end
   
saveas(h,[dir2,'fractalDimension.png'])
close all 




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


function fractalCalculation( mydata, Xbin,flag)

%    data1=d1.dataAvgMeanFeatures{i,j};   % individual point 
%    data2=(d2.realization_AvgMeanFeatures{i}(j,:))';
%    
%     
        
    g=fittype(@(m,c,x) (m*x+c));
    RZ=[];
    PZ=[];
    PHZ=[];
    HZ=[];
    for j=1:length(Xbin)
        
        if flag==1
            real_size= mydata{3,j};
            real_rg=mydata{11,j};
        end    
        
        if flag==2
            real_size=mydata(j,1);
            real_rg=mydata(j,2);
        end
            
        if Xbin(j)<=0.2
            RZ=[RZ; [log(real_size),log(real_rg)]];
        elseif (Xbin(j)>0.2)&(Xbin(j)<=0.5)
            PZ=[PZ; [log(real_size),log(real_rg)]];
        elseif (Xbin(j)>0.5)&(Xbin(j)<=0.7)
            PHZ=[PHZ; [log(real_size),log(real_rg)]];
        elseif (Xbin(j)>0.7)&(Xbin(j)<=1)
            HZ=[HZ; [log(real_size),log(real_rg)]];
        end 

    end
    
    zone={RZ,PZ,PHZ,HZ};
    zonecolor={'r.','b.','g.','k.'};
    legname={'RZ','PZ','PHZ','HZ'};
    
    for i=1:length(zone)
        data=zone{i};
        
        
        t=min(data(:,1)): 0.1: max(data(:,1))+0.000001;
        %t=unique(data(:,1));
        
        X=[];
        Y=[];
        Z=[];
        for j=1:length(t)-1
            index=find((data(:,1)>= t(j) ) & (  data(:,1)<t(j+1)));
            %index=find(data(:,1)==t(j));
            X=[X; t(j)];
            Y=[Y;  nanmean(data(index,2))];
            Z=[Z;  nanstd(data(index,2))];
        end
            
        
      
        %plot(data(:,1),data(:,2),zonecolor{i});
            
         
        index = find(~isnan(Y));
        x=X(index); y=Y(index); 
        plot(x,y,zonecolor{i});
        hold on 
        
        %[x,y]
        
        min_y1=y-Z(index); max_y1=y+Z(index);  index=1:length(x);
        stacky2=(min_y1);stacky1=(max_y1);
        fillxx=x(index([1:end end:-1:1 1]));fillyy=[stacky1(index); stacky2(index(end:-1:1)); stacky1(1)];
        h=fill(fillxx,fillyy,zonecolor{i}(1),'EdgeColor','none','facealpha',0.2);
        
        
       
        [f,~]=fit(x,y,g,'startpoint',[1,0]);
        p(i)=plot(x,f(x),[zonecolor{i}(1),'-'],'linewidth',0.5);

        Df=1/f.m; 
        lname{i}=strcat(legname{i},' :', 'D_f =',sprintf('%0.2f',Df));
    
    end

        legend(p, lname,'location','southeast')
        xlabel('log(cluster size)');
        

   
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



