
function goodcellindex=RemoveBadCells(centroid, volume,surfaces)
    tic
    [~,firstgood]=unique(centroid,'rows');
    %[length(centroid),length(firstgood)]
% spacing=[1, 1, 1];
% delta_spacing=2;
% if the user didn't define 'exclude_boundary_points' we set it to false:
if ~exist('exclude_boundary_points', 'var')
    exclude_boundary_points = false;
end
% calculating the triangulation and the volume of each triangle:
N=centroid(firstgood,:);

TRI = delaunay(N(:,1), N(:,2), N(:,3));
clear neighbor
edges=[];
count=1;
for i = 1 : size(N,1)
    temp=[];
    for j=1:size(TRI,1)
        for k=1:size(TRI,2)
            if TRI(j,k)==i
                temp=[temp,TRI(j,:)];
            end
        end
    end

   
    neighborList=setdiff(unique(temp),i);
    %neighbor(i,1)=length(neighborList{i});
    ids= neighborList;
    
    for j=1:length(ids)
        a=min(ids(j),i);
        b=max(ids(j),i);
        edges(count,:)=[firstgood(a),firstgood(b)];
        count=count+1;
    end
 
end

%size(edges)
[~,ia]=unique(edges,'rows');
edges=edges(ia,:);
%size(edges)

 badcell=[];
 for j=1:size(edges,1)
     
     dist=pdist(centroid(edges(j,:),:));
     if dist<8     
             cell1=surfaces{edges(j,1)}; actualVol(1)=volume( edges(j,1)   );      
             cell2=surfaces{edges(j,2)}; actualVol(2)=volume( edges(j,2)   );  
             combined=[cell1;cell2];
             
             [~,combV]=convhull(combined);
             [~,cv(1)]=convhull(cell1);
             [~,cv(2)]=convhull(cell2);

             %volumeFraction= combV/sum(actualVol);
             volumeFraction = combV/sum(cv); 
             
             roundness = actualVol./cv; 
             
             
             if volumeFraction<1
                        %[sa,sb]=min(actualVol);
                        [sa,sb]=min(roundness);
                        if sb==1
                            badid=edges(j,1);
                        end
                        if sb==2
                            badid=edges(j,2);
                        end
                        badcell=[badcell;badid];
             end
     end
 end    

  goodcellindex=setdiff(firstgood,badcell);
  toc
  %length(goodcellindex)  
  
end
