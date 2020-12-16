
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

       

dt_path_wt={ '../data/Nuclei_and_Cells_DT_S18_m6_wt/', '../data/Nuclei_and_Cells_DT_S84_m3_wt/',...
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

profilesize=50;        
myinterval=linspace(0,1,profilesize);       
h1=figure();
%count=1;
for bonetype=1:5
        path=allpath{bonetype};
        clear celavg
        clear nucavg
        
for fi=1:length(path)
        a1=load([path{fi},'all_cells_nuclei.mat']);
        b1=load([path{fi},'all_cells.mat']);
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
            
    
        plot3(nuccent(:,1),nuccent(:,2),nuccent(:,3),mycolor{fi})
        resting=min(nuccent(:,3));
        hypertrophic=max(nuccent(:,3));
        
        [y1,x1,ny2,x2]= makeProfile(nuccent(:,3),nuc(:,20),profilesize); 
        [y1,x1,cy2,x2]= makeProfile(celcent(:,3),cel(:,20),profilesize); 
        %ankit(count,:)=[mean(nuc(1:200,3)),mean(cel(1:200,3))];
        %count=count+1;
        celavg(:,fi)=cy2;
        nucavg(:,fi)=ny2;
        hold on 
        
end


data{bonetype}.cel=celavg;
data{bonetype}.nuc=nucavg;
end
save('avgBoneDataDelaunay.mat','data')

axis image 
view(90,0);
xlabel('X');
ylabel('y');
zlabel('z');
saveas(h1,'all.png');


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
                end
                youtput=Yinterval+binsize/2;
end
