

clear all 

d3=load('dataSave/SaveVariablesForPlot_larger.mat');
d2=load('dataSave/SaveVariablesForPlot_medium.mat');
d1=load('dataSave/SaveVariablesForPlot_smaller.mat');


 



h=figure();
XL=0.06;XR=0.01;XGap=0.012;Row=2;
YT=0.05;YB=0.09;YGap=0.07;Col=3;
Width=(1-XL-XR-((Col-1)*XGap))/Col;
Height=(1-YT-YB-((Row-1)*YGap))/Row;
YPos=1-YT-Height; 

set(gcf, 'PaperSize', [12 7]); %7
set(gcf, 'PaperPosition', [0 0 12 7]);

for i=1:Row
    XPos=XL;
    for j=1:Col
        chro=j+(i-1)*Col;
        if chro<=6
            marray=[XPos,YPos,Width,Height];
            subplot('Position',marray);
            
            
            
            
            
            
            h=figure;
            set(gcf, 'PaperSize', [5 3]); %7
            set(gcf, 'PaperPosition', [0 0 5 3]);
            
            
            if chro==1 
                plot_OOP(d1,13:15); title(['OOP: local, cutoff: 0.75*c'], 'fontweight','normal','fontsize',14)
            end  
            
            if chro==2
                plot_OOP(d2,13:15); title(['OOP: local, cutoff: c'], 'fontweight','normal','fontsize',14)
            end   
            
            if chro==3 
                plot_OOP(d3,13:15); title(['OOP: local, cutoff: 1.25*c'], 'fontweight','normal','fontsize',14)
            end   
            
            if chro==4
                plot_OOP(d1,16:18); title(['OOP: global, cutoff: 0.75*c'], 'fontweight','normal','fontsize',14)
            end  
            
            if chro==5
                plot_OOP(d2,16:18); title(['OOP: global, cutoff: c'], 'fontweight','normal','fontsize',14)
            end   
            
            if chro==6 
                plot_OOP(d3,16:18); title(['OOP: global, cutoff: 1.25*c'], 'fontweight','normal','fontsize',14)
            end   
                
                
            if chro>3
                ylim([-0.55,1.05])
            else
                ylim([0,1.1])
                set(gca,'xticklabel',[])
            end
            
           
            
                
            %if j==1
            ylabel('<OOP>','fontsize',14)
            %else
                set(gca,'yticklabel',[])
            %end
            %if i==Row
                    xlabel('Long axis of bone','fontsize',14)
            %end
            
            
            
            saveas(h,['OOP_',num2str(chro),'.png'])
            close all 
            
            
            
        
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



legend(p,'PC1','PC2','PC3','location','best')
end


