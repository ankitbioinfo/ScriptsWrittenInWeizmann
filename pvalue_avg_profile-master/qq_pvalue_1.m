clear

% a1=load('pvalueMain/pvalue_Volume.dat');
% a2=load('pvalueMain/pvalue_SA.dat');
% 
% a3=load('pvalueMain/pvalue_Sphericity.dat');
% a4=load('pvalueMain/pvalue_delaunay.dat');
% 
% a5=load('pvalueMain/pvalue_ncRatio.dat');
% 
% 
% b=[a1;a2;a3;a4;a5];


%KStest ln(volume)
% b=[3.8007e-39
% 2.9177e-197
% 6.0733e-141
% 2.0620e-44
% 4.8279e-315
% 4.8861e-263
% 3.6888e-13
% 3.8310e-19
% 1.2601e-157
% 7.5292e-144 
% 1.7359e-11
% 4.8405e-49
% 1.9002e-78
% 8.4983e-31
% 4.0037e-68
% 2.4699e-105
% 1.7198e-16
% 3.9146e-07
% 6.1140e-72
% 1.3736e-95];

%T test ln(volume)
b=[2.9603e-24
2.8886e-153
2.1466e-111
4.6751e-24
1.2881e-247
2.6471e-198
8.3114e-05
1.8630e-05
1.2900e-152
1.9387e-143 
1.1309e-12
1.8469e-76
1.6336e-133
4.7373e-28
3.7763e-126
1.1752e-192
9.6558e-12
1.2915e-08
3.6177e-91
4.7889e-128];




fig=figure();
%data=-log10(b(5:end,:));
data=-log10(b);

pd = makedist('Uniform');
h=qqplot(data,pd);
set(h, 'Color', 'w','MarkerEdgeColor', 'b')
ylabel('Observed -log_{10}(p-values)','fontsize',11)
xlabel('Expected -log_{10}(p-values)','fontsize',11)

axObjs = fig.Children;
dataObjs = axObjs.Children;
x = dataObjs(1).XData;
y = dataObjs(1).YData;
z = dataObjs(1).ZData;

hold on 
index=1:5;
par=polyfit(x(index),y(index),1);
myx=linspace(0,2,100);
yfit=polyval(par,myx);
plot(myx,yfit,'r-','linewidth',2);     
%ht=text(5.5,7,strcat('Y=',sprintf('%0.2f',par(1)),'*x +',sprintf('%0.2f',par(2))),'color','b');
       
axis([0,1,0,80])
set(gca,'fontsize',10)
plot([0,1],[12,12],'k--');
title('QQ plot of p-value versus Uniform Distribution','fontweight','normal','fontsize',11)

% saveas(gcf, 'QQplot_pvalues')
% print('QQplot_pvalues','-dpng','-r300');

saveas(gcf, 'QQ_pvalue_Ttest')
print('QQ_pvalues_Ttest','-dpng','-r300');




