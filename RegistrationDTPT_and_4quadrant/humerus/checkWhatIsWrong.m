
clear


a1=load('wt/all_cells_nuclei.mat');
b1=load(['wt/all_cells.mat']);


profilesize=51;        
myinterval=linspace(0,1,profilesize);  

 nuc=a1.all_cells_nuclei;
 %cel=b1.all_cells;
 
 %celcent(:,1:3)=[cel(:,5),cel(:,6),cel(:,7)];
 nuccent(:,1:3)=[nuc(:,5),nuc(:,6),nuc(:,7)];
 
 nuccent=nuccent-mean(nuccent);
 %celcent=celcent-mean(celcent);
 
%   plot3(nuccent(:,1),nuccent(:,2),nuccent(:,3),'b.'); hold on 

 tileid=unique(nuc(:,1));
 
 shp = alphaShape(nuccent);
 [tetrahedron,V] = alphaTriangulation(shp);
trep = triangulation(tetrahedron, V);
[tri, V] = freeBoundary(trep);
%  
  plot3(V(:,1),V(:,2),V(:,3),'r.'); hold on 
  
  for i=1:length(tileid)
      index=find(nuc(:,1)==tileid(i));      
      ori=mean(nuccent(index,:));
      text(ori(1),ori(2),ori(3),num2str(tileid(i)), 'fontsize',20)
  end
 
  removetile=[28,29,46,30,45,47,59,31,44,48,58,43,49,50,56,57,58,55]; %wt 

%  removetile=[28,29,30,31,41,27,32,17,40,42]; mt 
  for i=1:length(removetile)
        index=find(nuc(:,1)==removetile(i));
        plot3(nuccent(index,1),nuccent(index,2),nuccent(index,3),'b.'); 
  end
  
  
  figure
 working=setdiff(tileid,removetile)';
   for i=1:length(working)
        index=find(nuc(:,1)==working(i));
        plot3(nuccent(index,1),nuccent(index,2),nuccent(index,3),'r.'); 
        hold on 
   end
  
   for i=1:length(tileid)
      index=find(nuc(:,1)==tileid(i));      
      ori=mean(nuccent(index,:));
      text(ori(1),ori(2),ori(3),num2str(tileid(i)), 'fontsize',20)
  end
   
 
 
%  [vec,val]=eig(cov(nuccent));
%  
% mu=mean(V);
% d = sqrt(diag(val));
% factor=3; quiver3(mu(1),mu(2),mu(3),factor*vec(1,1),factor*vec(2,1),factor*vec(3,1),d(1),'r','LineWidth',10);
% factor=3; quiver3(mu(1),mu(2),mu(3),factor*vec(1,2),factor*vec(2,2),factor*vec(3,2),d(2),'k','LineWidth',5);
% factor=3; quiver3(mu(1),mu(2),mu(3),factor*vec(1,3),factor*vec(2,3),factor*vec(3,3),d(3),'m','LineWidth',2);
%  
% factor=-3; quiver3(mu(1),mu(2),mu(3),factor*vec(1,1),factor*vec(2,1),factor*vec(3,1),d(1),'r','LineWidth',10);
% factor=-3; quiver3(mu(1),mu(2),mu(3),factor*vec(1,2),factor*vec(2,2),factor*vec(3,2),d(2),'k','LineWidth',5);
% factor=-3; quiver3(mu(1),mu(2),mu(3),factor*vec(1,3),factor*vec(2,3),factor*vec(3,3),d(3),'m','LineWidth',2);
% 
%  [bone_curvature,~]=movingMeanAverage(nuccent);
%  plot3(bone_curvature(:,1),bone_curvature(:,2),bone_curvature(:,3),'g-','linewidth',10)
% 
%  view(-90,0)
%  xlabel('X')
%  ylabel('Y')
%  zlabel('Z')
 
 
 
%  
%  ankur=V*vec(:,[3,2,1]);
%  
%  figure
%  plot3(ankur(:,1),ankur(:,2),ankur(:,3),'b.')
 




            
    
%    %USE THIS COMMAND FOR SHIFTING THE CENTROIDS 
% 
%  nV=centroid_shift_of_bone(bone_curvature,V);
%      
%         
%  figure
%  plot3(V(:,1),V(:,2),V(:,3),'b.');
%  hold on 
%  plot3(bone_curvature(:,1),bone_curvature(:,2),bone_curvature(:,3),'r-','linewidth',20)
%  axis image
%  
%  figure
%  plot3(nV(:,1),nV(:,2),nV(:,3),'b.')
%  hold on 
%  plot3(bone_curvature(:,1),bone_curvature(:,2),bone_curvature(:,3),'r-','linewidth',20)
%  axis image

 
 
%  [y1,x1,ny2,x2]= makeProfile(nuccent(:,3),nuc(:,2),profilesize); 
%  [y1,x1,cy2,x2]= makeProfile(celcent(:,3),cel(:,2),profilesize); 
%  ankit=TilemakeProfile(nuccent(:,3),nuc(:,1));
%  
%  nucavg=ny2;
%  celavg=cy2;
%  
%  plot(myinterval,celavg,'b.-')
%  
%  temp=[];
%  for i=1:length(ankit)
%      temp=[temp;ankit{i}];
%      ankit{i}'
%  end
%  length(unique(temp))
%  
%  
 
 
 
function [data,centroidZ,Xinterval,youtput]= makeProfile(centroidZ,data,profilesize)
                %centroidZ=centroid(:,3);
                minl=min(centroidZ); maxl=max(centroidZ);
                Yinterval=linspace(minl-0.00001,maxl+0.00001,profilesize);
                Xinterval=zeros(1,length(Yinterval));
                binsize=Yinterval(2)-Yinterval(1);
                horizontalCount=cell(1,length(Yinterval));
                for k=1:length(Yinterval)
                    horizontalCount{k}=[];
                    for tt=1:length(centroidZ)
                        if (centroidZ(tt) > Yinterval(k)) & (centroidZ(tt) < Yinterval(k+1))
                            horizontalCount{k}=[horizontalCount{k},tt];
                        end
                    end
                    Xinterval(k)=mean(data(horizontalCount{k}));
                    ankit(k)=length(horizontalCount{k});
                end
                youtput=Yinterval+binsize/2;
                ankit 
end





 
function [Xinterval]= TilemakeProfile(centroidZ,data,tileid)
                %centroidZ=centroid(:,3);
                profilesize=11;
                minl=min(centroidZ); maxl=max(centroidZ);
                Yinterval=linspace(minl-0.00001,maxl+0.00001,profilesize);
                %Xinterval=zeros(1,length(Yinterval));
                binsize=Yinterval(2)-Yinterval(1);
                horizontalCount=cell(1,length(Yinterval));
                for k=1:length(Yinterval)
                    horizontalCount{k}=[];
                    for tt=1:length(centroidZ)
                        if (centroidZ(tt) > Yinterval(k)) & (centroidZ(tt) < Yinterval(k+1))
                            horizontalCount{k}=[horizontalCount{k},tt];
                        end
                    end
                    Xinterval{k}=unique(data(horizontalCount{k}));
                end
                youtput=Yinterval+binsize/2;
end



        

function nv=centroid_shift_of_bone(bone_curvature,V)
	x=bone_curvature(:,3);
	nV=zeros(size(V)); 
	for i=1:size(V,1)
        for j=1:length(x)
             dist(j)=abs(V(i,3)-x(j));
        end
        [sa,sb]=min(dist);
        localcenter=bone_curvature(sb,[1:2]);    
        nv(i,[1,2])=V(i,[1,2])-localcenter;
        nv(i,3)=V(i,3);
    end
end 




function [CurvatureAxisLine,dd]=movingMeanAverage(nuccent)
    zmin=min(nuccent(:,3)); zmax=max(nuccent(:,3));
    interval=zmin-0.001:1:zmax+0.001;
    count=1;
     for i=1:length(interval)
            left=interval(i)-100;
            right=interval(i)+100;
            index=find((nuccent(:,3)>=left) & (nuccent(:,3)<right));
            if length(index)>=100
                CurvatureAxisLine(count,:)=mean(nuccent(index,:));
                count=count+1;
            end
     end
     
%        d=CurvatureAxisLine;
%        CS = cat(1,0,cumsum(sqrt(sum(diff(d,[],1).^2,2))));
%        dd = interp1(CS, d, unique([CS(:)' linspace(0,CS(end),100)]),'pchip');
         dd=[];
     
     
end