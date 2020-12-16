
clear all 
%https://blog.minitab.com/blog/adventures-in-statistics-2/how-to-compare-regression-lines-between-different-models
% red=load('redData.dat');
% blue=load('blueData.dat');

red(:,1)=[60,65,70,75,80,85];
red(:,2)=[1.76,2.78,3.93,5.06,6.55,7.24];
blue(:,1)=red(:,1);
blue(:,2)=[2.98,3.68,5.54,6.8,8.55,9.55];


plot(red(:,1),red(:,2),'r.')
hold on 
plot(blue(:,1),blue(:,2),'b.')

out1=errorAndPvalue(red(:,1),red(:,2));
  q=plot(red(:,1),out1.yfit,'r-','linewidth',0.6);
out2=errorAndPvalue(blue(:,1),blue(:,2));
  q=plot(blue(:,1),out2.yfit,'b-','linewidth',0.6);

tvalue=abs(out1.p(1) - out2.p(1)) / sqrt( out1.SE_slope^2 + out2.SE_slope^2)

 mypvalue=2*tcdf(tvalue,length(blue)+length(red)-4, 'upper')



function  output=errorAndPvalue(x,y)
               g=fittype(@(m,c,x) (m*x+c));              
               [p,s]=polyfit(x,y,1);
               [yfit,d] = polyval(p,x,s); 

               yresid = y-yfit;
               SSresid = sum(yresid.^2);
               SStotal = sum( (y-mean(y)).^2); % same thing   (length(y)-1) * var(y)
               nrsq = 1 - SSresid/SStotal;
               nrsq_adj = 1 - SSresid/SStotal * (length(y)-1)/(length(y)-length(p));
                 
               [f,stat]=fit(x,y,g,'startpoint',[-0.1,1]);
               
               % Two Sided 2*tcdf(2.29,99,'upper') 
               % One Sided   tcdf(2.29,99,'upper')
               %https://handbook-5-1.cochrane.org/chapter_7/7_7_3_2_obtaining_standard_deviations_from_standard_errors_and.htm
               %https://stattrek.com/regression/slope-test.aspx   
            
               SE = sqrt(SSresid/(length(y)-2)) / sqrt( sum((x-mean(x)).^2) ); % Standard Error 
               SE1 = sqrt(SSresid/(length(y)-2)) * sqrt( sum(x.^2) /( length(y)*sum((x-mean(x)).^2)) ); % Standard Error 
               mypvalue=2*tcdf(abs(f.m/SE),stat.dfe, 'upper');
               [f.m/SE,stat.dfe,mypvalue];
               output.pvalue=mypvalue;
               output.yfit=yfit;
               output.Rsq1=nrsq;
               output.SE_slope=SE;
               output.SE_intercept=SE1;
               output.p=p;
end
