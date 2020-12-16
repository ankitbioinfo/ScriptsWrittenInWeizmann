
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
% 
% 
%         
% for fi=5
%       nuccent=[];
%    for bonetype=[1,2]
%         path=allpath{bonetype}{fi};
%         a1=load(['../',path,'all_cells_nuclei.mat']);
%         nuc=a1.all_cells_nuclei;       
%         nuccent=[nuccent;nuc(:,5:7)];
%    end
%    globalmean(fi,:)=mean(nuccent);
% end
% 
% plot3(nuccent(:,1),nuccent(:,2),nuccent(:,3),'g.')
% 
% 



        
for fi=4
      h=figure; 
   for bonetype=[3,4]
        path=allpath{bonetype}{fi};
        a1=load(['../',path,'all_cells_nuclei.mat']);
        nuc=a1.all_cells_nuclei;
        clear nuccent
        clear cent        
        
        if (bonetype==3)|(bonetype==1)
            nuccent(:,1:3)=[nuc(:,5),nuc(:,6),nuc(:,7)];
            %cent(:,1:3)=[nuc(:,5),nuc(:,6),-nuc(:,7)];
        else
            nuccent(:,1:3)=[nuc(:,5),nuc(:,6),nuc(:,7)];
            %cent(:,1:3)=[nuc(:,5),nuc(:,6),nuc(:,7)];
        end
        
        C=mean(nuccent); zgap=max(nuccent(:,3))-min(nuccent(:,3));
        nuccent1=nuccent-C;
        
%         figure
%         plot3(nuccent1(:,1),nuccent1(:,2),nuccent1(:,3),'g.')
%         figure 
%         
        
        
%         cent1=cent-mean(cent);
%         
%         ankit=[mean(cent(:,3)) C(:,3)]
%         ankur=[min(cent1(:,3)) max(cent1(:,3)) min(nuccent1(:,3)) max(nuccent1(:,3))]

        
        
            
        s=strsplit(path,'Nuclei_and_Cells_');
        outputpath=strcat('../../data/ColumnRelated/only_PZ_column2/Columnar_Structure_Prediction_8_100/',s{2});
        a=load([outputpath,'Deviation_PCs.dat']);
        Z=a(:,[9,12]); zmin=min(Z(:)); zmax=max(Z(:));
        
        [zmin,zmax,zgap];
        if (bonetype==3)|(bonetype==1)
        zmin=-zmin+C(1,3);  zmax=-zmax+C(1,3);
        else
        zmin=zmin+C(1,3);  zmax=zmax+C(1,3);
        end
        [zmin,zmax,zgap];
        
        
        xneg=find(nuccent1(:,1) < 0); yneg=find(nuccent1(:,2) < 0);
        Pair=ones(length(nuccent),2);   Pair(xneg,1)=0;  Pair(yneg,2)=0;  
        signPair=binarySum(Pair);
        
        ind=find(signPair==1);
        plot3(nuccent(ind,1),nuccent(ind,2),nuccent(ind,3),'r.'); hold on 
        ind=find(signPair==2);
        plot3(nuccent(ind,1),nuccent(ind,2),nuccent(ind,3),'b.')
        
        ind=find(signPair==3);
        plot3(nuccent(ind,1),nuccent(ind,2),nuccent(ind,3),'y.')
        ind=find(signPair==4);
        plot3(nuccent(ind,1),nuccent(ind,2),nuccent(ind,3),'g.')
        
        xmin=min(nuccent(:,1)); xmax=max(nuccent(:,1));
        ymin=min(nuccent(:,2)); ymax=max(nuccent(:,2));

        P=[xmin,ymin,zmin; xmax,ymin,zmin; xmax,ymax,zmin;xmin,ymax,zmin; xmin,ymin,zmin];
        plot3(P(:,1),P(:,2),P(:,3),'k-');
        P=[xmin,ymin,zmax; xmax,ymin,zmax; xmax,ymax,zmax;xmin,ymax,zmax; xmin,ymin,zmax];
        plot3(P(:,1),P(:,2),P(:,3),'k-');
        
      
        
   end
        saveas(h, ['GP_Quadrant/GP_quadrant_',s{2}(4:strlength(s{2})-1),'.fig'])
        %saveas(h, ['GP_Quadrant/GP_quadrant_',s{2}(4:strlength(s{2})-1),'.png'])

        %close all   
end

xlabel('X');
ylabel('Y');
zlabel('Z');










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
                        if (centroidZ(tt) > Yinterval(k)) & (centroidZ(tt) <= Yinterval(k+1))
                            horizontalCount{k}=[horizontalCount{k},tt];
                        end
                    end
                    Xinterval(k)=sum(data(horizontalCount{k}));
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
       