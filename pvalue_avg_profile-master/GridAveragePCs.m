
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
               
du_path_wt={'../DUdata/DU_S51_m2_wt/','../DUdata/DU_S84_m2_wt/','../DUdata/DU_S84_m3_wt/'};

       
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
        
        N=load([path{fi},'nuclei_grid.mat']);
        C=load([path{fi},'cells_grid.mat']);
        
        v=N.nuclei_grid.object_number.vals(:,:,:,1);
        t=v(:);  temp=find(t>10);
        ng=sum(t(temp))/length(temp);
        
        v=C.cells_grid.object_number.vals(:,:,:,1);
        t=v(:);  temp=find(t>10);
        cg=sum(t(temp))/length(temp);
       
        celavg(fi)=cg;
        nucavg(fi)=ng;
        
        
end


data{bonetype}=[celavg;nucavg;nucavg./celavg];
end



