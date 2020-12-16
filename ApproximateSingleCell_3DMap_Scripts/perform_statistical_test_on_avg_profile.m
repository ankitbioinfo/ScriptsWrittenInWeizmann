function  perform_statistical_test_on_avg_profile(wtpath,mtpath,wtBoneName,mtBoneName,savepath,opt)


for i=1:length(wtpath)
    if length(wtpath{i})~=0
        wtbonedata{i}=Read_GrowthPlate(wtpath{i});
    else
        wtbonedata{i}=[];
    end
end


for i=1:length(mtpath)
    if length(mtpath{i})~=0
        mtbonedata{i}=Read_GrowthPlate(mtpath{i});
    else
        mtbonedata{i}=[];
    end
end



     
%%% WT -MT test 
path=[savepath, 'WT_MT/'];
for i=1:length(wtbonedata)
    for j=1:length(mtbonedata)
            bone1=wtbonedata{i};
            bone2=mtbonedata{j};
            WTpairs=strcat(wtBoneName{i},'-',mtBoneName{j});
            legendTitle={ wtBoneName{i}, mtBoneName{j} };            
            if ((length(bone1)==3)&(length(bone2)==3))
                perform_statistical_test(bone1,bone2,path,WTpairs,legendTitle,opt);
            end
    end
end


% WT -WT test 
path=[savepath, 'WT_WT/'];
for i=1:length(wtbonedata)
    for j=i+1:length(wtbonedata)
            bone1=wtbonedata{i};
            bone2=wtbonedata{j};
            
            WTpairs=strcat(wtBoneName{i},'-',wtBoneName{j});
            legendTitle={ wtBoneName{i}, wtBoneName{j} };            

            if ((length(bone1)==3)&(length(bone2)==3))
                perform_statistical_test(bone1,bone2,path,WTpairs,legendTitle,opt);
            end
    end
end

 

%%%% MT -MT test 
path=[savepath,'MT_MT/'];
for i=1:length(mtbonedata)
    for j=i+1:length(mtbonedata)
            bone1=mtbonedata{i};
            bone2=mtbonedata{j};
            
            WTpairs=strcat(mtBoneName{i},'-',mtBoneName{j});
            legendTitle={ mtBoneName{i}, mtBoneName{j} };            

            if ((length(bone1)==3)&(length(bone2)==3))
                perform_statistical_test(bone1,bone2,path,WTpairs,legendTitle,opt);
            end
    end
end



%%%% MT test 
path=[savepath,'MT/'];
for i=1:length(mtbonedata)
    bone1=mtbonedata{i};
    WTpairs=strcat(mtBoneName{i});
    legendTitle={mtBoneName{i}};            
    if (length(bone1)==3)
        perform_individual_statistical_test(bone1,path,WTpairs,legendTitle,opt);
    end
end


%%%% WT test 
path=[savepath,'WT/'];
for i=1:length(wtbonedata)
    bone1=wtbonedata{i};
    WTpairs=strcat(wtBoneName{i});
    legendTitle={wtBoneName{i}};            
    if (length(bone1)==3)
        perform_individual_statistical_test(bone1,path,WTpairs,legendTitle,opt);
    end
end



end



function  perform_individual_statistical_test(dtwt,savepath,growthplateTitle,legendTitle,opt)
    ylabelCelorNuc{17}='[volume^{2/3} / surface area]';      
    units{17}=''; 
    titlename={'Cell', 'Nuclei','Crossed'};    
    Feature.Legend=legendTitle;
    
     %Volume^2/3 over surface area statistics for individual gp  
     if opt.onlyAllometricTest==1
     for object=1:2
        for feature=17
            clear wtFeature
            clear mtFeature 
            for sample=1:length(dtwt{object})
                 wtFeature(:,sample)=dtwt{object}{sample}{feature}(:,3);
                 X(:,1)=dtwt{object}{sample}{feature}(:,1);
            end

            Feature.TitleName=strcat(growthplateTitle,' - ',titlename{object});
            Feature.Ylabel= ylabelCelorNuc{feature};
            Feature.SaveName=[growthplateTitle,'_','allometry_vol_power_2_over_3_divide_by_sa'];
            Feature.save_folder= strcat(savepath, titlename{object},'/');
            Feature.Unit=units{feature};
            statistical_test_for_vol_over_sa_one_gp(wtFeature,X,Feature,opt);
        end
     end
     end


end




function  perform_statistical_test(dtwt,dtmt,savepath,growthplateTitle,legendTitle,opt)
    % First index: [Cell,nuclei,cross], 
    % Second index: Sample number 
    % Third Index: Features 
    % run over feature 11 [cell,nuclei]; 4 [ crossed]
    ylabelCelorNuc={'number density','delaunay density','ln(volume)','ln(surface area)','sphericity',...
          'PC1 coefficient', 'PC2 coefficient',  'PC3 coefficient', 'deviation in PC1 orientation',... 
          'deviation in PC2 orientation', 'deviation in PC3 orientation','volume','surface area', ...
         'PC2/PC1 ratio','PC3/PC1 ratio','PC3/PC2 ratio', '[volume^{2/3} / surface area]'};
      
    units={'', '[per 75 \mum^3]','','','','','','','\circ','\circ','\circ','\mum^3','\mum^2','','','',''}; 
    titlename={'Cell', 'Nuclei','Crossed'};
    
    Feature.Legend=legendTitle; 
    
   
    if opt.AllStatisticalTest==1
    for object=1:2
        for feature=1:16
            clear wtFeature
            clear mtFeature 
            for sample=1:length(dtwt{object})
                 wtFeature(:,sample)=dtwt{object}{sample}{feature}(:,3);
                 X(:,1)=dtwt{object}{sample}{feature}(:,1);
            end

            for sample=1:length(dtmt{object})
                 mtFeature(:,sample)=dtmt{object}{sample}{feature}(:,3);
            end
            Feature.TitleName=strcat(growthplateTitle,' - ',titlename{object});
            Feature.Ylabel= ylabelCelorNuc{feature};
            Feature.SaveName=[growthplateTitle,'_',strrep(ylabelCelorNuc{feature},' ','_')];
            Feature.SaveName=strrep(Feature.SaveName,'/','_');
            Feature.save_folder= strcat(savepath, titlename{object},'/');
            Feature.Unit=units{feature};
            statistical_test_function(wtFeature,mtFeature,X,Feature);
        end
    end
    end
    
     % Volume^2/3 over surface area statistics between WT and MT 
     if opt.onlyAllometricTest==1
     for object=1:2
        for feature=17
            clear wtFeature
            clear mtFeature 
            for sample=1:length(dtwt{object})
                 wtFeature(:,sample)=dtwt{object}{sample}{feature}(:,3);
                 X(:,1)=dtwt{object}{sample}{feature}(:,1);
            end

            for sample=1:length(dtmt{object})
                 mtFeature(:,sample)=dtmt{object}{sample}{feature}(:,3);
            end
            Feature.TitleName=strcat(growthplateTitle,' - ',titlename{object});
            Feature.Ylabel= ylabelCelorNuc{feature};
            Feature.SaveName=[growthplateTitle,'_','allometry_vol_power_2_over_3_divide_by_sa'];
            Feature.save_folder= strcat(savepath, titlename{object},'/');
            Feature.Unit=units{feature};
            statistical_test_for_vol_over_sa_two_gp(wtFeature,mtFeature,X,Feature,opt);
        end
     end
     end
    
    
    
     ylabelCross={'centroid shift','log(centroid shift)','nc ratio','number density',...
      'proportion corresponding', 'PC1 alignment', 'PC2 alignment', 'PC3 alignment'};
     units={'\mum','', '','','','','',''}; 
     if opt.AllStatisticalTest==1
     for object=3
        for feature=1:8
            clear wtFeature
            clear mtFeature 
            for sample=1:length(dtwt{object})
                 wtFeature(:,sample)=dtwt{object}{sample}{feature}(:,3);
                 X(:,1)=dtwt{object}{sample}{feature}(:,1);
            end

            for sample=1:length(dtmt{object})
                 mtFeature(:,sample)=dtmt{object}{sample}{feature}(:,3);
            end
            Feature.TitleName=strcat(growthplateTitle,' - ',titlename{object});
            Feature.Ylabel= ylabelCross{feature};
            Feature.SaveName=[growthplateTitle,'_',strrep(ylabelCross{feature},' ','_')];
            Feature.save_folder= strcat(savepath, titlename{object},'/');
            Feature.Unit=units{feature};
            statistical_test_function(wtFeature,mtFeature,X,Feature);
        end
     end
     end
    
    
    
    
end





function data=Read_GrowthPlate(gp)
for fi=1:length(gp)
    celpath=[gp{fi},'figures/cells/normalized_avg_individual_spatial_profiles/'];
    nucpath=[gp{fi},'figures/nuclei/normalized_avg_individual_spatial_profiles/'];
    crosspath=[gp{fi},'figures/crossed/normalized_avg_individual_spatial_profiles/'];
    celdata{fi}=Read_CellorNucleiAverageProfiles(celpath);
    nucdata{fi}=Read_CellorNucleiAverageProfiles(nucpath);
    crossdata{fi}=Read_CrossedAverageProfiles(crosspath); 
end
    data={celdata,nucdata,crossdata};
end





function feature=Read_CellorNucleiAverageProfiles(path)
        feature{1}=load([path,'norm_avg_ind_spat_prof_number_density.dat']);
        feature{2}=load([path,'norm_avg_ind_spat_prof_delaunay_density.dat']);
        feature{3}=load([path,'norm_avg_ind_spat_prof_log(volume).dat']);
        feature{4}=load([path,'norm_avg_ind_spat_prof_log(surface_area).dat']);
        feature{5}=load([path,'norm_avg_ind_spat_prof_sphericity.dat']);
        feature{6}=load([path,'norm_avg_ind_spat_prof_PC1_coefficient.dat']);
        feature{7}=load([path,'norm_avg_ind_spat_prof_PC2_coefficient.dat']);
        feature{8}=load([path,'norm_avg_ind_spat_prof_PC3_coefficient.dat']);
        feature{9}=load([path,'norm_avg_ind_spat_prof_Deviation_in_PC1_orientation.dat']);
        feature{10}=load([path,'norm_avg_ind_spat_prof_Deviation_in_PC2_orientation.dat']);
        feature{11}=load([path,'norm_avg_ind_spat_prof_Deviation_in_PC3_orientation.dat']);
        feature{12}=load([path,'norm_avg_ind_spat_prof_volume.dat']);
        feature{13}=load([path,'norm_avg_ind_spat_prof_surface_area.dat']);
        feature{14}=load([path,'norm_avg_ind_spat_prof_Ratio_PC2_PC1_coefficients.dat']);
        feature{15}=load([path,'norm_avg_ind_spat_prof_Ratio_PC3_PC1_coefficients.dat']);
        feature{16}=load([path,'norm_avg_ind_spat_prof_Ratio_PC3_PC2_coefficients.dat']);
        feature{17}=load([path,'norm_avg_ind_spat_prof_Allometry_vol_sa.dat']);
end


function feature=Read_CrossedAverageProfiles(path)
        feature{1}=load([path,'norm_avg_ind_spat_prof_centroid_shift.dat']);
        feature{2}=load([path,'norm_avg_ind_spat_prof_log(centroid_shift).dat']);
        feature{3}=load([path,'norm_avg_ind_spat_prof_nc_ratio.dat']);
        feature{4}=load([path,'norm_avg_ind_spat_prof_number_density.dat']);
        feature{5}=load([path,'norm_avg_ind_spat_prof_proportion_corresponding.dat']);
        feature{6}=load([path,'norm_avg_ind_spat_prof_PC1_alignment.dat']);
        feature{7}=load([path,'norm_avg_ind_spat_prof_PC2_alignment.dat']);
        feature{8}=load([path,'norm_avg_ind_spat_prof_PC3_alignment.dat']);
end












