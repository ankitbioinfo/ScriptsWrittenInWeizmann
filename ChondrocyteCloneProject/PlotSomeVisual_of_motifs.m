 
clear all 
 
a=load('Tibia/dataSave/SaveVariablesForPlot_medium.mat');
b=load('NullModel/SaveVariablesForPlot20_medium.mat');
n=size(a.dataAvgMeanFeatures,1); % # of features 
m=size(a.dataAvgMeanFeatures,2); % # of xbin





random=load('NullModel/MakeListColumnarStructurePrediction/DT_E185_nuclei/centroid_and_surface_cells.mat');
real=load('Tibia/MakeListColumnarStructurePrediction/centroid_and_surface_cells.mat');
x_real=load('Tibia/dataSave/AllFeaturesSave_medium.mat');
x_random=load('NullModel/Realization_medium/AllFeaturesSave3.mat');


radius=load('Tibia/dataSave/ParameterForRandomModel_basedon_nuc_medium.dat');
radius=radius(1:481,:);


[d1,~]=readClusterFile('Real_node_name/name1.dat');
[d2,~]=readClusterFile('Random_node_name/Random_1/name1.dat');

radius=radius(1:length(d2),:);

index=find(radius>0 & radius<1);
for i=1:length(index)
   
    x(i,1)=sum(d1{index(i)}==-1);
    y(i,1)=sum(d2{index(i)}==-1);
end

[sa,sb1]=sort(abs(x-y),'descend');

sb=sb1(13);
id=index(sb);
oop1=x(sb,1);
oop2=y(sb,1);
[sa(13), oop1,oop2,length(id1),length(id2)]



id1=d1{id};
id2=d2{id};


nucid1=x_real.LCC{id};
nucid2=x_random.LCC{id};


factor=length(id1)/oop1; 

count=1;
for j=1:factor:length(id1)
    subplot(3,4,count)
    
    for i=1:length(nucid1)
        flag=0;
        first=nucid1(i);
        for k=j:j+factor-2
            second=id1(k);
            if first==second
                   flag=1;
            end
        end
        
        v=real.nuc{first};
        if flag==0      
                    plot3(v(:,1),v(:,2),v(:,3),'r.','markersize',0.1);
        else
                    plot3(v(:,1),v(:,2),v(:,3),'b.','markersize',0.1);
        end
        hold on 
    end
    
   
    view(-6,66)
    count=count+1;
     if count>13
        break
    end
    
end

figure


count=1;
for j=1:factor:length(id2)
    subplot(3,3,count)
    
    for i=1:length(nucid2)
        flag=0;
        first=nucid2(i);
        for k=j:j+factor-2
            second=id2(k);
            if first==second
                   flag=1;
            end
        end
        
        v=random.nuc{first};
        if flag==0      
                    plot3(v(:,1),v(:,2),v(:,3),'r.','markersize',0.1);
        else
                    plot3(v(:,1),v(:,2),v(:,3),'b.','markersize',0.1);
        end
        hold on 
    end
    
   view(-91,15)
    
    count=count+1;
     if count>9
        break
    end
    
end











%     
% subplot(1,2,2)
% for i=1:length(id2)
%     v=random.nuc{id2(i)};
%     plot3(v(:,1),v(:,2),v(:,3),'.');
%     hold on 
% end
% %title(['random(local OOP 1): ', sprintf('%0.2f',oop2)])
% title(['random(cluster PC 1): ', sprintf('%0.2f',oop2)])
% 
% 



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
