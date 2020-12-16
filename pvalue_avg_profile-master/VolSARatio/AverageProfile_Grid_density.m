
clear all 

dt_path_wt={ '../data/Nuclei_and_Cells_DT_S18_m6_wt/', '../data/Nuclei_and_Cells_DT_S17_m2_wt/',...
             '../data/Nuclei_and_Cells_DT_S84_m3_wt/', '../data/Nuclei_and_Cells_DT_S51_m2_wt/',...
             '../data/Nuclei_and_Cells_DT_S84_m4_wt/'};

pt_path_wt = {  '../data/Nuclei_and_Cells_PT_S18_m6_wt/','../data/Nuclei_and_Cells_PT_S17_m2_wt/',...  
                '../data/Nuclei_and_Cells_PT_S84_m3_wt/','../data/Nuclei_and_Cells_PT_S51_m2_wt/',...
                '../data/Nuclei_and_Cells_PT_S84_m4_wt/'};

dt_path_mut= {'../data/Nuclei_and_Cells_DT_S17_m1_mut/', '../data/Nuclei_and_Cells_DT_S18_m2_mut/' ,...
              '../data/Nuclei_and_Cells_DT_S84_m1_mut/', '../data/Nuclei_and_Cells_DT_S84_m5_mut/'};
    
pt_path_mut = {'../data/Nuclei_and_Cells_PT_S17_m1_mut/', '../data/Nuclei_and_Cells_PT_S18_m2_mut/',...
               '../data/Nuclei_and_Cells_PT_S84_m1_mut/', '../data/Nuclei_and_Cells_PT_S84_m5_mut/', };
               
du_path_wt={'../data/Nuclei_and_Cells_DU_S51_m2_wt/','../data/Nuclei_and_Cells_DU_S84_m2_wt/','../data/Nuclei_and_Cells_DU_S84_m3_wt/'};

       
allpath={dt_path_wt; pt_path_wt; dt_path_mut; pt_path_mut; du_path_wt}; 
%allpath={du_path_wt};
  
% the columns contains the individual nuclei features
% 1 - stack id
% 2 - volume
% 3 - surface area
% 4 - sphericity
% 5-7 - centroid x,y,z coordinates
% 8-10 - PC1 x,y,z orientation
% 11-13 - PC2 x,y,z orientation
% 14-16 - PC3 x,y,z orientation
% 17-19 - PC1,PC2,PC3 latent coefficient
% 20 - Delaunay density

mycolor={'r.','b.','g.','m.','k.'};
nucallcolor={'r--','b--','g--','m--','k--'};

fcelallcolor={'ro-','bo-','go-','mo-','ko-'};
fnucallcolor={'ro--','bo--','go--','mo--','ko--'};

profilesize=51;        
myinterval=linspace(0,1,profilesize);       
h1=figure();
%count=1;

BigData=[];
for bonetype=1:5
        path=allpath{bonetype};
        clear celavg
        clear nucavg
        
for fi=1:length(path)
        a1=load(['../',path{fi},'all_cells_nuclei.mat']);
        b1=load(['../',path{fi},'all_cells.mat']);
        nuc=a1.all_cells_nuclei;
        cel=b1.all_cells;
        clear celcent
        clear nuccent
       
        if (bonetype==3)|(bonetype==1)
            celcent(:,1:3)=[cel(:,5),cel(:,6),-cel(:,7)];
            nuccent(:,1:3)=[nuc(:,5),nuc(:,6),-nuc(:,7)];
        else
            celcent(:,1:3)=[cel(:,5),cel(:,6),cel(:,7)];
            nuccent(:,1:3)=[nuc(:,5),nuc(:,6),nuc(:,7)];
        end
        
        nuccent=nuccent-mean(nuccent);
        celcent=celcent-mean(celcent);
%         
%         b1= load(['../',path{fi},'cells_grid.mat']);
%         a1= load(['../',path{fi},'nuclei_grid.mat']);
%         ankit{bonetype,fi}=b1.grid.z_bins;
%         nucgrid=a1.
%         celgrid=b1.

            
       
        
       
        
        resting=min(nuccent(:,3));
        hypertrophic=max(nuccent(:,3));
        
        [y1,x1,ny2,x2]= makeProfile(nuccent(:,3),ones(size(nuccent,1),1),profilesize,nuccent(:,1),nuccent(:,2)); 
        [y1,x1,cy2,x2]= makeProfile(celcent(:,3),ones(size(celcent,1),1),profilesize,celcent(:,1),celcent(:,2)); 
        %sum(isnan(data{5}.cel))
        %ankit(count,:)=[mean(nuc(1:200,3)),mean(cel(1:200,3))];
        %count=count+1;
        celavg(:,fi)=cy2;
        nucavg(:,fi)=ny2;
        
end


data{bonetype}.cel=celavg;
data{bonetype}.nuc=nucavg;
end

xlabel('X');
ylabel('Y');
zlabel('Z');

% doubletPC=pca(BigData);
% GC=mean(BigData);
% 
% factor=2000;quiver3(GC(1,1),GC(1,2),GC(1,3),factor*doubletPC(1,1),factor*doubletPC(2,1),factor*doubletPC(3,1),1,'k:','LineWidth',2);
% factor=-2000;quiver3(GC(1,1),GC(1,2),GC(1,3),factor*doubletPC(1,1),factor*doubletPC(2,1),factor*doubletPC(3,1),1,'k:','LineWidth',2);
% 
% % 
% 
% l= [-389.5,-297.7,-193.9;389.5,297.7,193.9];
% plot3(l(:,1),l(:,2),l(:,3),'b.-');

save('avg_Bone_GridDensity.mat','data')
%save('avgBoneData_SA2vol50.mat','data')

% axis image 
% view(90,0);
% xlabel('X');
% ylabel('y');
% zlabel('z');
% saveas(h1,'all.png');


function [data,centroidZ,Xinterval,youtput]= makeProfile(centroidZ,data,profilesize,x,y)
                %centroidZ=centroid(:,3);
                minl=min(centroidZ); maxl=max(centroidZ);
                Yinterval=linspace(minl-0.00000001,maxl+0.00000001,profilesize);
                Xinterval=zeros(1,length(Yinterval));
                binsize=Yinterval(2)-Yinterval(1);
                horizontalCount=cell(1,length(Yinterval));
                for k=1:length(Yinterval)
                    horizontalCount=[];
                    xscale=[];
                    yscale=[];
                    for tt=1:length(centroidZ)
                        if (centroidZ(tt) > Yinterval(k)) & (centroidZ(tt) <= Yinterval(k+1))
                            horizontalCount=[horizontalCount,tt];
                            xscale=[xscale,x(tt)];
                            yscale=[yscale,y(tt)];
                        end
                    end
                    xbin=max(xscale)-min(xscale);
                    ybin=max(yscale)-min(yscale);
                    grid_volume_bin=binsize*xbin*ybin;
                    if length(grid_volume_bin)>0
                    Xinterval(k)=sum(data(horizontalCount))/grid_volume_bin;
                    end
                end
                youtput=Yinterval+binsize/2;
end


function b=binarySum(a)
          for i=1:size(a,1)
              if ((a(i,1)==1) & (a(i,2)==1))
                  b(i)=1;
              elseif ((a(i,1)==0) & (a(i,2)==1))
                  b(i)=2;
              elseif ((a(i,1)==0) & (a(i,2)==0))    
                  b(i)=3;
              elseif ((a(i,1)==1) & (a(i,2)==0))
                  b(i)=4;
              end
          end
end
       