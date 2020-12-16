

clear all 

a=load('Feature16_L_M_S.dat');
b=load('Feature17_L_M_S.dat');
c=load('Feature18_L_M_S.dat');


h=figure;

p(1)=plot(a(:,2),'r.-'); hold on 
p(2)=plot(b(:,2),'b.-');
p(3)=plot(c(:,2),'g.-');


legend(p,'PC1','PC2','PC3','location','best')

saveas(h,'global_med.png'); close all 


h=figure;

p(1)=plot(a(:,1),'r.-'); hold on 
p(2)=plot(b(:,1),'b.-');
p(3)=plot(c(:,1),'g.-');


legend(p,'PC1','PC2','PC3','location','best')

saveas(h,'global_larger.png'); close all 


h=figure;

p(1)=plot(a(:,3),'r.-'); hold on 
p(2)=plot(b(:,3),'b.-');
p(3)=plot(c(:,3),'g.-');


legend(p,'PC1','PC2','PC3','location','best')

saveas(h,'global_smaller.png'); close all 