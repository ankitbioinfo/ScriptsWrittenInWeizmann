function [] = draw_normalized_averaged_individual_spatial_profiles(cel,bonetype,is_bone_humerus,opt,manualRotationFilePath, cellornucleiorcrossed)



%names = fieldnames(features);
profilesize= opt.ProximalDistal_binsize+1;        
myinterval=linspace(0,1,profilesize-1);   

if strmatch(cellornucleiorcrossed,'crossed')
    
    % 1 - stack id
    % 2-4 - centroid x,y,z coordinates
    % 5 - proportion corresponding
    % 6 - n/c ratio
    % 7 - centroid shift
    % 8-10 - PC1 x,y,z orientation - nucleus
    % 11-13 - PC2 x,y,z orientation - nucleus
    % 14-16 - PC3 x,y,z orientation - nucleus
    % 17-19 - PC1,PC2,PC3 latent coefficient - nucleus
    % 20-22 - PC1 x,y,z orientation - cells
    % 23-25 - PC2 x,y,z orientation - cells
    % 26-28 - PC3 x,y,z orientation - cells
    % 29-31 - PC1,PC2,PC3 latent coefficient - cells
    

    
    
      for f=[1,5,6,7,8,9:11]
          if f==1
               data=ones(size(cel,1),1);
               featureName='number density';
               unit='';
          elseif f==5
              data=cel(:,f);       
              featureName='proportion corresponding';
              unit='';
          elseif f==6
              data=cel(:,f);       
              featureName='nc ratio';
              unit='';
          elseif f==7
              data=cel(:,f); featureName='centroid shift'; unit='\mum';
          elseif f==8
              data=log(cel(:,7)); featureName='log(centroid shift)'; unit='';
          elseif f==9    
              data=abs(dot(cel(:,8:10),cel(:,20:22),2));
              featureName='PC1 alignment';unit='';
          elseif f==10    
              data=abs(dot(cel(:,11:13),cel(:,23:25),2));
              featureName='PC2 alignment';unit='';
          elseif f==11    
              data=abs(dot(cel(:,14:16),cel(:,26:28),2));
              featureName='PC3 alignment';unit='';
              
          end
      
     
        if (bonetype==1)
            celcent(:,1:3)=[cel(:,2),cel(:,3),-cel(:,4)];
        else
            celcent(:,1:3)=[cel(:,2),cel(:,3),cel(:,4)];
        end
        
        celcent=celcent-repmat(mean(celcent),size(celcent,1),1  )  ;
        
        % for humerus bone need to rotate the bone to correct the spatial profile
        if is_bone_humerus
            %disp(['humerus', cellornucleiorcrossed])
            manualCenterAxisLine=load([manualRotationFilePath,'centerAxisLine.dat']);
            [vec,val]=eig(cov( manualCenterAxisLine));
            celcent=celcent*vec;
        end
        
        
        
        
        if f==1
            [ydata,xdata]= makeProfileDensity(celcent(:,3),data,profilesize); 
        else
            [ydata,xdata]= makeProfile(celcent(:,3),data,profilesize);
        end
        
        celavg=ydata(1:end-1);
        
        h=figure;
        plot(myinterval,celavg,'b.-'); 
        xlabel('Bone long axis RZ-HZ','fontsize',10)
        ylabel(['average ',featureName,' ', unit],'fontsize',10)
        set(gca,'fontsize',10)

        
        title(['Bin average along P-D axis: ',featureName],'Interpreter','None','fontweight','normal','fontsize',9);
        if opt.save_figs,
            if ~exist([opt.save_folder],'dir')
               mkdir([opt.save_folder]);
            end
            dlmwrite([opt.save_folder,'norm_avg_ind_spat_prof_',strrep(featureName,' ','_'),'.dat'],[myinterval',xdata(1:end-1)',ydata(1:end-1)'],'\t')
            saveas(h,[opt.save_folder,strrep(featureName,' ','_'),'_norm_avg_ind_spat_prof']);
            saveas(h,[opt.save_folder,strrep(featureName,' ','_'),'_norm_avg_ind_spat_prof','.png']);
        end
        close all
      end

else
    
    
% each row of all_cells_nuclei is nucleus
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
    

        for f = [1,2,3,4,5,6,8,11,14,17:20,21,22:24] 
                if f==1
                   data=ones(size(cel,1),1);
                   featureName='number density'; saveName=featureName;
                   unit='';
                elseif f==2 
                   data=log(cel(:,f));       
                   featureName='log(volume)'; saveName=featureName;
                   unit='';
                elseif f==3 
                   data=log(cel(:,f));    
                   featureName='log(surface area)'; saveName=featureName;      
                   unit='';
                elseif f==4 
                   data=cel(:,f);    
                   featureName='sphericity';  saveName=featureName;  
                   unit='';
                elseif f==5
                   data=cel(:,2);       
                   featureName='volume'; saveName=featureName;
                   unit='\mum^3';
                elseif f==6
                   data=cel(:,3);    
                   featureName='surface area'; saveName=featureName;      
                   unit='\mum^2';   
                elseif f==20    
                   data=75^3*cel(:,f);    %delaunay density multipled with 75 micron cube that means how many 
                   % cell/nuclei present in 75^3 cube 
                   featureName='delaunay density';  saveName=featureName; 
                   unit='[per 75 \mum^3]';
                elseif f==17 
                   data=cel(:,f);    
                   featureName='PC1 coefficient';  unit='';  saveName=featureName;       
                elseif f==18 
                   data=cel(:,f);    
                   featureName='PC2 coefficient';  unit=''; saveName=featureName;
                elseif f==19
                   data=cel(:,f);    
                   featureName='PC3 coefficient';  unit=''; saveName=featureName;
                   
                elseif f==22
                   data=cel(:,18)./cel(:,17);    
                   featureName='PC2/PC1 coefficient';  unit='';  saveName='Ratio_PC2_PC1_coefficients'; 
                   
                elseif f==23
                   data=cel(:,19)./cel(:,17);    
                   featureName='PC3/PC1 coefficient';  unit='';  saveName='Ratio_PC3_PC1_coefficients';
                   
                elseif f==24
                   data=cel(:,19)./cel(:,18);    
                   featureName='PC3/PC2 coefficient';  unit='';  saveName='Ratio_PC3_PC2_coefficients';   
                   
                elseif f==21
                   data=(cel(:,2).^(2/3))./cel(:,3);
                   featureName='[volume^{2/3} / surface area]'; unit=''; saveName='Allometry_vol_sa';
                   
                elseif f==8
                   vec=findMainVector(cel(:,8:10));
                   data=angleCompute(cel(:,8:10), vec' );                
                   featureName='deviation in PC1 orientation';  unit='\circ';  saveName=featureName;
                   
                elseif f==11
                   vec=findMainVector(cel(:,11:13));
                   data=angleCompute(cel(:,11:13), vec' );                
                   featureName='deviation in PC2 orientation';  unit='\circ';  saveName=featureName;
                   
                elseif f==14
                   vec=findMainVector(cel(:,14:16));
                   data=angleCompute(cel(:,14:16), vec' );                
                   featureName='deviation in PC3 orientation';  unit='\circ';   saveName=featureName;
                   
                end       
        
       
        if (bonetype==1)
            celcent(:,1:3)=[cel(:,5),cel(:,6),-cel(:,7)];
        else
            celcent(:,1:3)=[cel(:,5),cel(:,6),cel(:,7)];
        end
        
        celcent=celcent-repmat(mean(celcent),size(celcent,1),1  )  ;
        
        % for humerus bone need to rotate the bone to correct the spatial profile
        if is_bone_humerus
            %disp(['humerus', cellornucleiorcrossed])
            manualCenterAxisLine=load([manualRotationFilePath,'centerAxisLine.dat']);
            [vec,val]=eig(cov( manualCenterAxisLine));
            celcent=celcent*vec;
        end
        
        
        
        if f==1
            [ydata,xdata]= makeProfileDensity(celcent(:,3),data,profilesize); 
            CreateTableForTilePositionAcrossSpatialProfile(celcent(:,3),cel(:,1),profilesize,opt);
        else
            [ydata,xdata]= makeProfile(celcent(:,3),data,profilesize);
        end

        celavg=ydata(1:end-1);
        
        h=figure;
        plot(myinterval,celavg,'b.-'); 
        xlabel('Bone long axis RZ-HZ','fontsize',10)
        ylabel(['average ',featureName,' ', unit],'fontsize',10)
        set(gca,'fontsize',10)
        
        
        title(['Bin average along P-D axis: ',featureName],'Interpreter','None','fontweight','normal','fontsize',9);
        if opt.save_figs,
            if ~exist([opt.save_folder],'dir')
               mkdir([opt.save_folder]);
            end
            dlmwrite([opt.save_folder,'norm_avg_ind_spat_prof_',strrep(saveName,' ','_'),'.dat'],[myinterval',xdata(1:end-1)',ydata(1:end-1)'],'\t')
            saveas(h,[opt.save_folder,strrep(saveName,' ','_'),'_norm_avg_ind_spat_prof']);
            saveas(h,[opt.save_folder,strrep(saveName,' ','_'),'_norm_avg_ind_spat_prof','.png']);
        end
        close all
        end
end

end



function [Xinterval,youtput]= makeProfile(centroidZ,data,profilesize)
                minl=min(centroidZ); maxl=max(centroidZ);
                Yinterval=linspace(minl,maxl,profilesize);
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
                    Xinterval(k)=mean(data(horizontalCount{k}));
                end
                youtput=Yinterval+binsize/2;
end



function [Xinterval,youtput]= makeProfileDensity(centroidZ,data,profilesize)
                minl=min(centroidZ); maxl=max(centroidZ);
                Yinterval=linspace(minl,maxl,profilesize);
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


function CreateTableForTilePositionAcrossSpatialProfile(centroidZ,data,profilesize,opt)
                minl=min(centroidZ); maxl=max(centroidZ);
                profilesize=11;
                Interval_for_table=linspace(0,1,profilesize);
                Yinterval=linspace(minl,maxl,profilesize);
                binsize=Yinterval(2)-Yinterval(1);
                horizontalCount=cell(1,length(Yinterval));
                for k=1:(length(Yinterval)-1)
                    horizontalCount{k}=[];
                    for tt=1:length(centroidZ)
                        if (centroidZ(tt) > Yinterval(k)) & (centroidZ(tt) <= Yinterval(k+1))
                            horizontalCount{k}=[horizontalCount{k},tt];
                        end
                    end
                    Xinterval{k}=unique(data(horizontalCount{k}));
                    columnsize(k,1)=length(Xinterval{k});
                end
                A=nan(max(columnsize),length(columnsize));
                for k=1:(length(Yinterval)-1)
                    A((1:columnsize(k)),k) =  Xinterval{k};
                    var_name{k}=strcat('Xint_', num2str(10*Interval_for_table(k)),'_', num2str(10*Interval_for_table(k+1)));
                end
                
                
           
                if ~exist([opt.save_folder],'dir')
                   mkdir([opt.save_folder]);
                end
          
            
             outputfilename = [opt.save_folder,'tiles_spatial_profile_name_divide_by_10.xlsx'];
             data=array2table(A);
             data.Properties.VariableNames = var_name;    
             writetable(data,outputfilename,'Sheet','Tiles spatial profile');%,'WriteVariableNames',false);
           
                
end



function angleDeviation=angleCompute(uarray,v) 
         for i=1:size(uarray,1)
             u=uarray(i,:)';
             value=atan2(norm(cross(u,v)),dot(u,v));
             value=180/pi*value;
             if value>90
                 value=180-value;
             end
             angleDeviation(i,1)=value;
         end
end


function colorVec=findMainVector(vec)   
         vals = [vec; -vec];
         [coeff,score,latent] = pca(vals);
         indtemp2 = find(dot(vals(:,:)',repmat(coeff(:,1),[1,size(vals(:,:),1)]))>0);
         R = sqrt((sum(vals(indtemp2,1))).^2+(sum(vals(indtemp2,2))).^2+(sum(vals(indtemp2,3))).^2);
         colorVec = [sum(vals(indtemp2,1))/R, sum(vals(indtemp2,2))/R,sum(vals(indtemp2,3))/R];      
end
