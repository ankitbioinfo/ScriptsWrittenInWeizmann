
clear all 

% dt_path_wt={  '../data/Nuclei_and_Cells_DT_S17_m2_wt/',...
%              '../data/Nuclei_and_Cells_DT_S84_m3_wt/', '../data/Nuclei_and_Cells_DT_S51_m2_wt/'};
% 
% pt_path_wt = { '../data/Nuclei_and_Cells_PT_S17_m2_wt/',...  
%                 '../data/Nuclei_and_Cells_PT_S84_m3_wt/','../data/Nuclei_and_Cells_PT_S51_m2_wt/',...
%                 '../data/Nuclei_and_Cells_PT_S84_m4_wt/'};
% 
% dt_path_mut= {'../data/Nuclei_and_Cells_DT_S17_m1_mut/', '../data/Nuclei_and_Cells_DT_S18_m2_mut/' ,...
%               '../data/Nuclei_and_Cells_DT_S84_m1_mut/', };
%     
% pt_path_mut = {'../data/Nuclei_and_Cells_PT_S17_m1_mut/', '../data/Nuclei_and_Cells_PT_S18_m2_mut/',...
%                '../data/Nuclei_and_Cells_PT_S84_m1_mut/',  };


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
        
%         
%         npc=pca(nuccent(:,[1,2]));
%         cpc=pca(celcent(:,[1,2]));
%         ankit{bonetype,fi}=[npc,cpc];
            
    
        plot3(nuccent(:,1),nuccent(:,2),nuccent(:,3),mycolor{fi})
        resting=min(nuccent(:,3));
        hypertrophic=max(nuccent(:,3));
        
               
        
        nucVec=findMainVector(nuc(:,8:10));
        celVec=findMainVector(cel(:,8:10));
        
        
        
        clear phi
        for i=1:size(nuc,1)
            phi(i)=angleCompute(nuc(i,8:10)', nucVec' );
            %phi(i)=angleCompute(nuc(i,8:10)', [1;0;0] );
        end
        
        clear theta 
        for i=1:size(cel,1)
            theta(i)=angleCompute(cel(i,8:10)', celVec' );
            %theta(i)=angleCompute(cel(i,8:10)', [1;0;0] );
        end
        
        
        [y1,x1,ny2,x2]= makeProfile(nuccent(:,3),phi,profilesize); 
        [y1,x1,cy2,x2]= makeProfile(celcent(:,3),theta,profilesize); 
        %ankit(count,:)=[mean(nuc(1:200,3)),mean(cel(1:200,3))];
        %count=count+1;
        celavg(:,fi)=cy2;
        nucavg(:,fi)=ny2;
        hold on 
        
end


data{bonetype}.cel=celavg;
data{bonetype}.nuc=nucavg;
end
save('local_avgBonePC1.mat','data')

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



%angleCompute=@(u,v) atan2(norm(cross(u,v)),dot(u,v));
function value=angleCompute(u,v) 
         value=atan2(norm(cross(u,v)),dot(u,v));
         value=180/pi*value;
         
         if value>90
             value=180-value;
         end
         
end




function colorVec=findMainVector(vec)   
         vals = [vec; -vec];
         [coeff,score,latent] = pca(vals);
         indtemp2 = find(dot(vals(:,:)',repmat(coeff(:,1),[1,size(vals(:,:),1)]))>0);
         R = sqrt((sum(vals(indtemp2,1))).^2+(sum(vals(indtemp2,2))).^2+(sum(vals(indtemp2,3))).^2);
         colorVec = [sum(vals(indtemp2,1))/R, sum(vals(indtemp2,2))/R,sum(vals(indtemp2,3))/R];      
end

