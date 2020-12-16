

clear all 

d1=load('dataSave/SaveVariablesForPlot_medium.mat');



% %Features={volume,surface_area,clusterSize,sphericity,clusterPC1, clusterPC2, clusterPC3,clusterPC2_by_PC1,clusterPC3_by_PC1,clusterPC3_by_PC2,...
%     rg,VolumeFraction,LOP(:,1),LOP(:,2),LOP(:,3),GOP(:,1),GOP(:,2),GOP(:,3), plane_of_division_to_bone,plane_of_division_to_cell(:,1),...
%     plane_of_division_to_cell(:,2),plane_of_division_to_cell(:,3), highest_mode, smallest_mode,clusterRadius, ...,
%     

%26-30  LocalBiaxial(:,1), LocalBiaxial(:,2), LocalBiaxial(:,3), LocalBiaxial(:,4), AngleBetweenClusterPC1AndBone_PD};
 





h=figure();
XL=0.07;XR=0.01;XGap=0.07;Row=2;
YT=0.05;YB=0.09;YGap=0.07;Col=2;
Width=(1-XL-XR-((Col-1)*XGap))/Col;
Height=(1-YT-YB-((Row-1)*YGap))/Row;
YPos=1-YT-Height; 

set(gcf, 'PaperSize', [9 6]); %7
set(gcf, 'PaperPosition', [0 0 9 6]);

ylabelname={'<S>', '<P>', '<D>', '<C>'};

for i=1:Row
    XPos=XL;
    for j=1:Col
        chro=j+(i-1)*Col;
        if chro<=6
            marray=[XPos,YPos,Width,Height];
            subplot('Position',marray);
            
            
            if chro==1 
                plot_OOP(d1,26); title(['cutoff: c'], 'fontweight','normal','fontsize',14)
            end  
            
            if chro==2
                plot_OOP(d1,27); title(['cutoff: c'], 'fontweight','normal','fontsize',14)
            end   
            
            if chro==3 
                plot_OOP(d1,28); title(['cutoff: c'], 'fontweight','normal','fontsize',14)
            end   
            
            if chro==4
                plot_OOP(d1,29); title(['cutoff: c'], 'fontweight','normal','fontsize',14)
            end  
            
           
%                 
%                 
%             if chro>3
%                 ylim([-0.55,1.05])
%             else
%                 ylim([0,1.1])
%                 set(gca,'xticklabel',[])
%             end
%             
           
            
                
           
            ylabel(ylabelname{chro},'fontsize',14)
            if i==Row
                    xlabel('Long axis of bone','fontsize',14)
            else
                    set(gca,'xticklabel',[]);
            end
        
      end
	XPos=XPos+Width+XGap;
    end
    YPos=YPos-YGap-Height;
end


saveas(h,'OOP_local_global.png'); close all 
              
                
                
function plot_OOP(d,vector)  


mycolor={'r.-','b.-','g.-'};
facealpha=[0.3,0.2,0.15];

count=1;
for k=vector
    a=d.AvgMeanFeatures{k};
    X=d.Xbin;
    p(count)=plot(X,a,mycolor{count}); 
    hold on 
    
    
     index = find(~isnan(d.AvgMeanFeatures{k}));
                    x=(d.Xbin(index))';
                    min_y1=d.AvgMeanFeatures{k}(index)-d.AvgStdFeatures{k}(index); max_y1=d.AvgMeanFeatures{k}(index)+d.AvgStdFeatures{k}(index);  
                    index=1:length(x);
                    stacky2=(min_y1);stacky1=(max_y1);
                    fillxx=x(index([1:end end:-1:1 1]));fillyy=[stacky1(index); stacky2(index(end:-1:1)); stacky1(1)];
                    h=fill(fillxx,fillyy,mycolor{count}(1),'EdgeColor','none','facealpha',facealpha(count));
           
    
    
    
    
    
    count=count+1;
end



%legend(p,'PC1','PC2','PC3','location','best')
end


