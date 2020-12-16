

clear all 

d=load('SaveVariablesForPlot20_medium.mat');



h=figure;
set(gcf, 'PaperSize', [5 3]); %7
set(gcf, 'PaperPosition', [0 0 5 3]);

mycolor={'r.-','b.-','g.-'};
facealpha=[0.3,0.2,0.15];

count=1;
h=figure;
for k=13:15
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
title('local')

saveas(h,'local_med.png'); close all 


