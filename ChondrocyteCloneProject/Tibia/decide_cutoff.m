

clear 

load('MakeListColumnarStructurePrediction/centroid_and_surface_cells.mat');

bone_mu=mean(centroid);
c=centroid-bone_mu;
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


binsize=11;
zone=linspace(RZ,HZ,binsize);
[distances,neighbor] = calculate_nuclei_density(c);

% cell2nuclei_median_ratio=[1.2800
%     1.2654
%     1.2545
%     1.2894
%     1.3246
%     1.4792
%     1.6774
%     1.7599
%     1.7362
%     1.6931];



network_edges_small=[];
network_edges_medium=[];
network_edges_large=[];



for i=1:length(zone)-1
      index=find((c(:,1) >zone(i)) &  (c(:,1) <=zone(i+1)));
      GP_zone_cell_diameters(i,:)=2*mean(fitellipsoid(index,:));
      
      cutoff_large=1.25*GP_zone_cell_diameters(i,1);
      cutoff_medium=GP_zone_cell_diameters(i,1);
      cutoff_small=0.75*GP_zone_cell_diameters(i,1);
      
      ankit(i,:)=[cutoff_small,cutoff_medium,cutoff_large];
      
      
      network_edges_small=[network_edges_small; find_edges_for_a_cutoff(neighbor,distances,index,cutoff_small)];
      network_edges_medium=[network_edges_medium; find_edges_for_a_cutoff(neighbor,distances,index,cutoff_medium)];
      network_edges_large=[network_edges_large; find_edges_for_a_cutoff(neighbor,distances,index,cutoff_large)];
end


dlmwrite('MakeListColumnarStructurePrediction/cells_cutoff.dat',ankit,'\t');




[~,ia]=unique(network_edges_small,'rows');  network_edges_small=network_edges_small(ia,:);
[~,ia]=unique(network_edges_medium,'rows');  network_edges_medium=network_edges_medium(ia,:);
[~,ia]=unique(network_edges_large,'rows');  network_edges_large=network_edges_large(ia,:);

[size(network_edges_small,1),size(network_edges_small,1), size(network_edges_large,1)]

save('MakeListColumnarStructurePrediction/alledges.mat',   'network_edges_small',  'network_edges_medium','network_edges_large');


LCC=LargestConnectedComponents(network_edges_small);  
fid=fopen(['MakeListColumnarStructurePrediction/ClusterSmallerCutoff.dat'],'w');
for i=1:length(LCC)
        for j=1:length(LCC{i})
                fprintf(fid,'%d ',LCC{i}(j));
        end
        fprintf(fid,'\n');
end
         

LCC=LargestConnectedComponents(network_edges_medium);  
fid=fopen(['MakeListColumnarStructurePrediction/ClusterMediumCutoff.dat'],'w');
for i=1:length(LCC)
        for j=1:length(LCC{i})
                fprintf(fid,'%d ',LCC{i}(j));
        end
        fprintf(fid,'\n');
end
       


LCC=LargestConnectedComponents(network_edges_large);  
fid=fopen(['MakeListColumnarStructurePrediction/ClusterLargerCutoff.dat'],'w');
for i=1:length(LCC)
        for j=1:length(LCC{i})
                fprintf(fid,'%d ',LCC{i}(j));
        end
        fprintf(fid,'\n');
end
       




function myedges=find_edges_for_a_cutoff(neighbor,distance,mylist,cutoff)

    myedges=[];
    count=1;
    for i=1:length(mylist)
        node=mylist(i);
        neigh=neighbor{node};
        for j=1:length(neigh)
            if distance{node}(j,1)<=cutoff
                myedges(count,:)=sort([node,neigh(j)]);
                count=count+1;
            end
        end
    end
end
                
        




function LCC=LargestConnectedComponents(edges)
        cellIds=unique(edges(:));
        for j=1:length(cellIds)
            old2new(cellIds(j),1)=j;
            new2old(j,1)=cellIds(j);    
        end
        [length(old2new),length(new2old),length(cellIds)];

        for i=1:size(edges,1)
            for j=1:2 
                newedgename(i,j)= old2new(edges(i,j));
            end
        end
        G=graph(newedgename(:,1),newedgename(:,2));
        bins=conncomp(G);
        % number of connected components 
        nocomp=unique(bins);
        %disp(['# of connected components  ', num2str(length(nocomp))]);
        for i=1:length(nocomp)
            numberOfObjectsInConnectedComponents(i)=sum(nocomp(i)==bins);
        end
        
        
        [sa,sb]=sort(numberOfObjectsInConnectedComponents,'descend');
        
        index=1;
        for i=1:length(sa)
            if sa(i)>1
                LCCIds=find(bins==nocomp(sb(i)));
                LCC{index}=new2old(LCCIds);
                index=index+1;
            end
        end
end








function [dist,neighborList] = calculate_nuclei_density(N)
% if the user didn't define 'exclude_boundary_points' we set it to false:
if ~exist('exclude_boundary_points', 'var')
    exclude_boundary_points = false;
end
% calculating the triangulation and the volume of each triangle:
TRI = delaunay(N(:,1), N(:,2), N(:,3));

clear neighbor
clear dist
for i = 1 : size(N,1)
    temp=[];
    for j=1:size(TRI,1)
        for k=1:size(TRI,2)
            if TRI(j,k)==i
                temp=[temp,TRI(j,:)];
            end
        end
    end

   
    ids=setdiff(unique(temp),i);
    neighborList{i}=ids;
  
    %neighbor(i,1)=length(neighborList{i});
      
    for k=1:length(ids)
        dist{i}(k,1)=pdist(N([i,ids(k)],:)); 
    end
%    
%     [sa,sb]=sort(dist);
%     neighbor(i,1)=min(dist);
%     
%     edges1(i,:) = [sort([i,ids(sb(1))]) sa(1)  ];  
%     edges2{i}=[ids(sb(2:end))',sa(2:end)'];
end



end
