function []=statistical_test_function(cel,myinterval,Feature,opt)

titlename=Feature.TitleName;
ylabelfeature=Feature.Ylabel;
FigureLegend=Feature.Legend;
unit=Feature.Unit;

limit1=opt.allometric.limit(1);
limit2=opt.allometric.limit(2);
limit3=opt.allometric.limit(3);
limit4=opt.allometric.limit(4);
          
%celallcolor={'r-','b-','g-','m-','k-'};
mycolor={'r','b','g','m','c'};
lt={'.','*','x','+','s'};

fcelallcolor={'ro','ro','ro','ro','ro','ro'};
fnucallcolor={'b^','b^','b^','b^','b^','b^'};

%tname={'wt', 'mut'};     

h1=figure();
XL=0.07;XR=0.03;XGap=0.07;Row=1;
YT=0.06;YB=0.11;YGap=0.15;Col=1;
Width=(1-XL-XR-((Col-1)*XGap))/Col;
Height=(1-YT-YB-((Row-1)*YGap))/Row;
YPos=1-YT-Height; 

set(gcf, 'PaperSize', [13 7]);
set(gcf, 'PaperPosition', [0 0 13 7]);


for i=1:Row
    XPos=XL;
    for j=1:Col
        chro=j+(i-1)*Col;
      
            marray=[XPos,YPos,Width,Height];
            subplot('Position',marray);
            testwithstd=chro;
           


        celavg=nanmean(cel,2);
        celstd=(nanstd(cel'))';
        
        p(1)=plot(myinterval,celavg,strcat(fcelallcolor{1},'-'),'linewidth',0.2);
        hold on 
        
     
        index = 1:length(myinterval); min_y1=celavg-celstd; max_y1=celavg+celstd;        
        x=myinterval;stacky2=(min_y1)';stacky1=(max_y1)';
        fillxx=x(index([1:end end:-1:1 1]));fillyy=[stacky1(index) stacky2(index(end:-1:1)) stacky1(1)];
        h=fill(fillxx,fillyy,'r','EdgeColor','none','facealpha',0.3);
        

        globalmin=min([min(min_y1)]);
        globalmax=max([max(max_y1)]);
      
 

            set(gca,'fontsize',11);            
            gap=[globalmin,globalmax];
            axis([0,1, globalmin, 1.2*globalmax]);
            
%              limit1=8;limit2=32; limit3=40;
               b(1)=myinterval(limit1);
               b(2)=myinterval(limit2);
               b(3)=myinterval(limit3);
             
               
             plot([b(1),b(1)],[0,110],'k:','linewidth',0.1);
             plot([b(2),b(2)],[0,110],'k:','linewidth',0.1);
             plot([b(3),b(3)],[0,110],'k:','linewidth',0.1);
               
               index=(1:limit1)'; 
               out1=dofitting(index,myinterval',celavg,0.202,0.198,0.194,'k','g-');           
               %text(myinterval(index(1)),0.210, ['pvalue = ', sprintf('%0.4f',mypvalue)],'fontsize',5)

               index=(limit1:limit2)';
               out1=dofitting(index,myinterval',celavg,0.188,0.184,0.18,'k','g-');
               %text(myinterval(index(1)),0.21, ['pvalue = ', sprintf('%0.4f',mypvalue)],'fontsize',5)
               index=(limit2:limit3)';
               out1=dofitting(index,myinterval',celavg,0.202,0.198,0.194,'k','g-');
               %text(myinterval(index(1)),0.21, ['pvalue = ', sprintf('%0.4f',mypvalue)],'fontsize',5
               index=(limit3:limit4)';
               out1=dofitting(index,myinterval',celavg,0.188,0.184,0.18,'k','g-');
               %text(myinterval(index(1)),0.21, ['pvalue = ', sprintf('%0.4f',mypvalue)],'fontsize',5);
               p(2)=out1.plot;
               
                
            if j==1
                ylabel(['average ',ylabelfeature,' ',unit]);
            end
          
            xlabel('Bone long axis RZ-HZ')
            n=length(FigureLegend);
            for lname=1:n
                FigureLegend{n+lname}=  strcat(FigureLegend{lname},':fit');
            end
            legend(p,FigureLegend,'location','north','fontsize',7);
            legend 'boxoff'; 
            %legend(q,legendarray,'location','northwest','fontsize',7);

        
    
   
	XPos=XPos+Width+XGap;
    end
    YPos=YPos-YGap-Height;
end


            if ~exist([Feature.save_folder],'dir')
               mkdir([Feature.save_folder]);
            end
            saveas(h1,[Feature.save_folder, Feature.SaveName]);
            saveas(h1,[Feature.save_folder, Feature.SaveName,'.png']);

            close all 


end



function   boundary=statisticalTest(celavg,nucavg)
              n=size(celavg,1);
              significant=[];
              for ii=1:n-4
                  ind2=ii:ii+4;
                  [h,pP]=ttest2(celavg(ind2),nucavg(ind2));
                  if pP<0.05
                      significant=[significant,ii];
                  end
              end    
              
              if length(significant)>1
                  boundary=[significant(1)];
              else
                  boundary=[];
              end
              for ii=2:length(significant)
                  if significant(ii)-significant(ii-1)~=1
                      boundary=[boundary,significant(ii-1),significant(ii)];
                  end
              end
              if length(significant)>1
                boundary=[boundary,significant(end)+4];
              end
end


function output=dofitting(index,myinterval,celavg,pos1,pos2,pos3,mycolor,plotcolor)
               x=(myinterval(index))';
               y=celavg(index);
               output=errorAndPvalue(x,y); 
               q=plot(x,output.yfit,plotcolor,'linewidth',0.6);
               text(x(1),pos1, ['Y = ',sprintf('%0.2f',output.p(1)),'X + ',sprintf('%0.2f',output.p(2))],'fontsize',5,'color',mycolor)
               text(x(1),pos2, ['R^2 = ', sprintf('%0.2f',output.Rsq1)],'fontsize',5,'color',mycolor)
               text(x(1),pos3, ['pvalue = ', sprintf('%0.2e',output.pvalue)],'fontsize',5,'color',mycolor)
               [output.Rsq1,output.Rsq2];
               output.plot=q;
end



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
               mypvalue=2*tcdf(abs(f.m/SE),stat.dfe, 'upper');
               [f.m/SE,stat.dfe,mypvalue];
               output.pvalue=mypvalue;
               output.yfit=yfit;
               output.Rsq1=nrsq;
               output.Rsq2=stat.rsquare;
               output.SE=SE;
               output.p=p;
end