
clear


a1=load('mt/all_cells_nuclei.mat');
b1=load(['mt/all_cells.mat']);
bone_curvature=load('mt/centerAxisLine.dat');



profilesize=51;        
myinterval=linspace(0,1,profilesize);  

 nuc=a1.all_cells_nuclei;
 cel=b1.all_cells;
 
 celcent(:,1:3)=[cel(:,5),cel(:,6),cel(:,7)];
 nuccent(:,1:3)=[nuc(:,5),nuc(:,6),nuc(:,7)];
 
 nuccent=nuccent-mean(nuccent);
 celcent=celcent-mean(celcent);
 

 
 [vec,val]=eig(cov(bone_curvature));
 nuccent=nuccent*vec;  
 celcent=celcent*vec;
    
%USE THIS COMMAND FOR SHIFTING THE CENTROIDS 
% nuccent1=centroid_shift_of_bone(bone_curvature,nuccent);
% celcent1=centroid_shift_of_bone(bone_curvature,celcent);

     
 theta=angleCompute(vec(:,1),[0;0;1])

        
 figure
 %plot3(V(:,1),V(:,2),V(:,3),'b.');
 plot3(nuccent(:,1),nuccent(:,2),nuccent(:,3),'b.'); 
 hold on 
 plot3(bone_curvature(:,1),bone_curvature(:,2),bone_curvature(:,3),'r-','linewidth',10)
 axis image
 
%  figure
%  plot3(nuccent1(:,1),nuccent1(:,2),nuccent1(:,3),'b.'); 
%  hold on 
%  plot3(bone_curvature(:,1),bone_curvature(:,2),bone_curvature(:,3),'r-','linewidth',10)
%  axis image

 
 figure
 [y1,x1,ny2,x2]= makeProfile(nuccent(:,3),nuc(:,2),profilesize); 
 [y1,x1,cy2,x2]= makeProfile(celcent(:,3),cel(:,2),profilesize); 
 ankit=TilemakeProfile(nuccent(:,3),nuc(:,1));
 
 
 
 nucavg=ny2;
 celavg=cy2;
 
 plot(myinterval,celavg,'b.-'); hold on 
 plot(myinterval,nucavg,'r.-'); 
 
 temp=[];
 for i=1:length(ankit)
     temp=[temp;ankit{i}];
     ankit{i}';
 end
 
 
 
 
 
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

