clear

bonetype=100;
% a1=load('N/all_cells_nuclei.mat');
% b1=load(['mt/all_cells.mat']);


% profilesize=51;        
% myinterval=linspace(0,1,profilesize);  
% 
%  nuc=a1.all_cells_nuclei;
%  %cel=b1.all_cells;
%  
%  %celcent(:,1:3)=[cel(:,5),cel(:,6),cel(:,7)];
%  nuccent(:,1:3)=[nuc(:,5),nuc(:,6),nuc(:,7)];

 a1=load('Nuclei_and_Cells_S126_m3_PT_DV/AlignedXYZ.mat');
 %a1=load('Nuclei_and_Cells_S17_m2_pt_thresh05/AlignedXYZ.mat');
 
 nuccent=a1.alphaShapePT*a1.vec;
 nuccent=nuccent-mean(nuccent);
 
 if bonetype==1
     nuccent(:,3)=-nuccent(:,3);
 end

 
 %celcent=celcent-mean(celcent);
 
 % plot3(nuccent(:,1),nuccent(:,2),nuccent(:,3),'b.'); hold on 

% tileid=unique(nuc(:,1));
 
 shp = alphaShape(nuccent);
 [tetrahedron,V] = alphaTriangulation(shp);
trep = triangulation(tetrahedron, V);
[tri, V] = freeBoundary(trep);
 
  plot3(V(:,1),V(:,2),V(:,3),'r.'); hold on 
   xlabel('X'); ylabel('Y'); zlabel('Z');
   
   cutoff=50;
  
  data=nuccent( find(nuccent(:,3)<cutoff),:) ; 
  binsize=100;
  index=find(V(:,3)<cutoff);
  
  gridx= (min(V(index,1))-binsize/2):binsize:(max(V(index,1))+binsize/2); 
  gridy= (min(V(index,2))-binsize/2):binsize:(max(V(index,2))+binsize/2);
 

count=1;
for i=1:length(gridx)-1
    for j=1:length(gridy)-1
            index= find( (data(:,1)>gridx(i)) & (data(:,1)<gridx(i+1)) & (data(:,2)>gridy(j)) & (data(:,2)<gridy(j+1)) );
            dist=[];
            for k=1:length(index)
                dist(k)=pdist( [data(index(k),:); [ gridx(i),gridy(j),2000   ]   ]);
            end
            if length(dist)>0
                [sa,sb]=min(dist);
                boundaryPoints(count,:)= data(index(sb),:);
                count=count+1;
            end
    end
end
    
hold on 
plot3(boundaryPoints(:,1),boundaryPoints(:,2),boundaryPoints(:,3),'bo','markersize',20,'markerfacecolor','b');

DM=[boundaryPoints(:,1),boundaryPoints(:,2),ones(size(boundaryPoints,1),1)];
B=DM\boundaryPoints(:,3);
[X,Y] = meshgrid(linspace(min(boundaryPoints(:,1)),max(boundaryPoints(:,1)),10), linspace(min(boundaryPoints(:,2)),max(boundaryPoints(:,2)),10));
Z = B(1)*X + B(2)*Y + B(3)*ones(size(X));
  
%meshc(X, Y, Z)
  
P1=[X(:),Y(:),Z(:)] ;
B0=mean(P1);
P = P1([1,100,50],:);
V1 = P(2,:)-P(1,:);
V2 = P(3,:)-P(1,:);
% add one point on the plane defined by the previous three points
P = [P; P(1,:)+V1+V2];

% orthogonal vector
W = cross(V2,V1);
% unit vector
U = W/sqrt(sum(W.^2,2));
% project B by the distance k
k = 1500;
B = B0+k*U;
% plot
quiver3(B0(1),B0(2),B0(3),k*U(1),k*U(2),k*U(3),'AutoScale','off');
plot3(B(1),B(2),B(3),'bs');


index=1;
for k=-3000:3000
   centerAxisLine(index,:)  = B0+k*0.5*U;
   index=index+1;
end

dlmwrite('centerAxisLine.dat',centerAxisLine,'\t');
dlmwrite('centroid.dat',nuccent,'\t');

plot3(centerAxisLine(:,1),centerAxisLine(:,2),centerAxisLine(:,3),'g.-')

