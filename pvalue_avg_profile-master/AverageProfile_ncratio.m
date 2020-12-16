
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
  
% each row of all_crossed is nucleus
% the columns contains the individual nuclei features
% 1 - stack id
% 2-4 - centroid x,y,z coordinates
% 5 - proportion corresponding
% 6 - n/c ratio
% 7 - centroid shift
% 8-10 - PC1 x,y,z orientation - nucleus
% 11-13 - PC1 x,y,z orientation - nucleus
% 14-16 - PC1 x,y,z orientation - nucleus
% 17-19 - PC1,PC2,PC3 latent coefficient - nucleus
% 20-22 - PC1 x,y,z orientation - cells
% 23-25 - PC1 x,y,z orientation - cells
% 26-28 - PC1 x,y,z orientation - cells
% 29-31 - PC1,PC2,PC3 latent coefficient - cells

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
        clear crossavg
        
for fi=1:length(path)
        a1=load([path{fi},'all_crossed.mat']);
        cross=a1.all_crossed;
       
        clear crosscent
       
        if (bonetype==3)|(bonetype==1)
            crosscent(:,1:3)=[cross(:,2),cross(:,3),-cross(:,4)];
        else
            crosscent(:,1:3)=[cross(:,2),cross(:,3),cross(:,4)];
        end
        
        crosscent=crosscent-mean(crosscent);
            
    
        plot3(crosscent(:,1),crosscent(:,2),crosscent(:,3),mycolor{fi})
%         resting=min(nuccent(:,3));
%         hypertrophic=max(nuccent(:,3));
        
        [y1,x1,cy2,x2]= makeProfile(crosscent(:,3),cross(:,6),profilesize); 
        %ankit(count,:)=[mean(nuc(1:200,3)),mean(cel(1:200,3))];
        %count=count+1;
        crossavg(:,fi)=cy2;
        hold on 
        
end


data{bonetype}.cross=crossavg;
end
save('avgBoneData_ncratio.mat','data')

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
