
clear all 
real=load('MakeListColumnarStructurePrediction/centroid_and_surface_cells.mat');

random=load('../NullModel/MakeListColumnarStructurePrediction/DT_E185_nuclei/centroid_and_surface_cells.mat','centroid');

dirname='MakeListColumnarStructurePrediction/';
alledges=load([dirname,'alledges.mat']);

c=real.centroid- mean(real.centroid);

bonetype=1;
if (bonetype==3)|(bonetype==1)
       c=[-c(:,1),c(:,2:3)];
else 
       c=[c(:,1:2),c(:,3)];
end



shp=alphaShape(c,100)

plot(shp,'edgecolor','none','facealpha',0.1);
hold on 

text(550,0,0,'Hypertropic zone')
text(-400,70,0,'Resting zone')

    [vec,val]=eig(cov(random.centroid));
    ovec=vec;
    oval=val;
    mu=mean(c);
    d = sqrt(diag(val));
%    factor=2; quiver3(mu(1),mu(2),mu(3),factor*vec(1,1),factor*vec(2,1),factor*vec(3,1),d(1),'b','LineWidth',2);
%     factor=2; quiver3(mu(1),mu(2),mu(3),factor*vec(1,2),factor*vec(2,2),factor*vec(3,2),d(2),'g','LineWidth',3);
     factor=1.9; quiver3(mu(1),mu(2),mu(3),factor*vec(1,3),factor*vec(2,3),factor*vec(3,3),d(3),'m','LineWidth',5);

%    factor=-2; quiver3(mu(1),mu(2),mu(3),factor*vec(1,1),factor*vec(2,1),factor*vec(3,1),d(1),'b','LineWidth',2);
%     factor=-2; quiver3(mu(1),mu(2),mu(3),factor*vec(1,2),factor*vec(2,2),factor*vec(3,2),d(2),'g','LineWidth',3);
    factor=-2.6; quiver3(mu(1),mu(2),mu(3),factor*vec(1,3),factor*vec(2,3),factor*vec(3,3),d(3),'m','LineWidth',5);






[LCC,LCC1]=readClusterFile([dirname,'ClusterMediumCutoff.dat']);

id1=LCC{295};

plot3(c(id1,1),c(id1,2),c(id1,3),'b.')

view([90,90])

set(findobj(gcf,'type','axes'),'Visible','off')



print(gcf,'bisectorplane1.png','-dpng','-r300');

figure 




for i=1:length(id1)
    v=real.nuc{id1(i)};
    maxdim(i,:)=max(v);
    mindim(i,:)=min(v);
    plot3(v(:,1),v(:,2),v(:,3),'b.');
    hold on 
    
    alignment(v)
    
end



%      network_edges_large: [5911?2 double]
%     network_edges_medium: [3589?2 double]
%      network_edges_small: [1494?2 double]
myedge=search_edges(alledges.network_edges_medium,id1);

for i=1:length(myedge)
    c=real.centroid(myedge(i,:),:);
    d(i,1)=pdist(c);
    %plot3(c(:,1),c(:,2),c(:,3),'k-')
    hold on 
    
    
     P1= c(1,:);  P2= c(2,:);
                 P= real.nuc{myedge(i,1)}(1,:);
                 R = cross(P1-P2, P1-P);
                 S = cross(R, P1-P2);
                 R=R/norm(R);
                 S=S/norm(S);
%                  alphap=1; betap=1;
%                  plane= mean([P1;P2])' + alphap*R' + betap*S';
                 %normCentroidVector=centroidVector/norm(centroidVector);
                 %vec=find_perp(normCentroidVector); vec=vec/norm(vec); [normCentroidVector', vec']
                 %[sum((P1-P2).*S), sum((P1-P2).*R)]
                 
                 
                
%                  norm_plane=norm(plane);
                 
                
                 
                 theta=linspace(0,2*pi,1000);
                 mid=mean([P1;P2]);
                 r=4; 
                 for k=1:length(theta)
                     qt=mid' + r* cos(theta(k))* R' +  r*sin(theta(k))* S'; 
                     q(k,:)=qt;
                     
                 end
                 
                 K=convhull(q(:,1),q(:,2),q(:,3));
                 
                 %plot3(q(:,1),q(:,2),q(:,3),'k.-');
                 trisurf(K,q(:,1),q(:,2),q(:,3),'FaceColor','black','edgecolor','none','facealpha',0.6);
                     
    
    
    
end


maxy=max(maxdim);
miny=min(mindim);

axis([miny(1),maxy(1),miny(2),maxy(2),miny(3),maxy(3)])

set(findobj(gcf,'type','axes'),'Visible','off')

view([-94,42])

print(gcf,'bisectorplane2.png','-dpng','-r300');



function alignment(C)
    [vec,val]=eig(cov(C));
    ovec=vec;
    oval=val;
    mu=mean(C);
    d = sqrt(diag(val));
    factor=2; quiver3(mu(1),mu(2),mu(3),factor*vec(1,1),factor*vec(2,1),factor*vec(3,1),d(1),'r','LineWidth',2);
    factor=2; quiver3(mu(1),mu(2),mu(3),factor*vec(1,2),factor*vec(2,2),factor*vec(3,2),d(2),'g','LineWidth',3);
    factor=2; quiver3(mu(1),mu(2),mu(3),factor*vec(1,3),factor*vec(2,3),factor*vec(3,3),d(3),'m','LineWidth',5);

    factor=-2; quiver3(mu(1),mu(2),mu(3),factor*vec(1,1),factor*vec(2,1),factor*vec(3,1),d(1),'r','LineWidth',2);
    factor=-2; quiver3(mu(1),mu(2),mu(3),factor*vec(1,2),factor*vec(2,2),factor*vec(3,2),d(2),'g','LineWidth',3);
    factor=-2; quiver3(mu(1),mu(2),mu(3),factor*vec(1,3),factor*vec(2,3),factor*vec(3,3),d(3),'m','LineWidth',5);

end







                
function myedges=search_edges(edges,vertices)
%     edges
%     vertices 
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