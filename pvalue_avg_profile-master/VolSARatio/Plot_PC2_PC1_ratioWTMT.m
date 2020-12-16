
clear all 

a21=load('../CreateMatFiles/avgBoneData_RatioPC2_PC1.mat');
a31=load('../CreateMatFiles/avgBoneData_RatioPC3_PC1.mat');
a32=load('../CreateMatFiles/avgBoneData_RatioPC3_PC2.mat');

% cell   <V> = 832.4551,  <A> = 547 
% nucleu <V> =  343.2248, <A> = 261 

           
%celallcolor={'r-','b-','g-','m-','k-'};
mycolor={'k','k','k','k','k'};
lt={'.','*','x','+','s'};


% fcelallcolor={'ro-','bo-','go-','mo-','co-'};
% fnucallcolor={'r^--','b^--','g^--','m^--','c^--'};

fcelallcolor={'r-','r-','r-','r-','r-'};
fnucallcolor={'b-','b-','b-','b-','b-'};

% tname={'DT:S18 m6 wt','DT:S17 m2 wt','DT:S84 m3 wt','DT:S51 m2 wt','DT:S84 m4 wt',...
%        'PT:S18 m6 wt','PT:S17 m2 wt','PT:S84 m3 wt','PT:S51 m2 wt','PT:S84 m4 wt',...
%        'DT:S17 m1 mut','DT:S18 m2 mut','DT:S84 m1 mut','DT: S84 m5 mut',...
%        'PT:S17 m1 mut','PT:S18 m2 mut','PT:S84 m1 mut','PT: S84 m5 mut',...      
%         'DU:S51 m2 wt', 'DU:S84 m2 wt', 'DU:S84 m3 wt'};
tname={'DT', 'PT','DT','PT','DU wt'};

g=fittype(@(alpha,beta,constant,A0,x) (alpha/beta + constant/A0*exp(-beta*x))    );
logistic=fittype(@(a,b,c,d,x) (d+b./(c+exp(-a*x))));


ylabels={'PC2/PC1 ratio','PC3/PC1 ratio','PC3/PC2 ratio'};
   


lname={'wt', 'mt','wt-fit','mt-fit'};
testwithstd=2;
profilesize=51;        
myinterval=linspace(0,1,profilesize-1);  
%myinterval1=myinterval(1:50);


h1=figure();
XL=0.07;XR=0.01;XGap=0.08;Row=3;
YT=0.05;YB=0.09;YGap=0.07;Col=2;
Width=(1-XL-XR-((Col-1)*XGap))/Col;
Height=(1-YT-YB-((Row-1)*YGap))/Row;
YPos=1-YT-Height; 

set(gcf, 'PaperSize', [10 12]); %7
set(gcf, 'PaperPosition', [0 0 10 12]);

for i=1:Row
    XPos=XL;
    for j=1:Col
        chro=j+(i-1)*Col;
        if chro<=6
            marray=[XPos,YPos,Width,Height];
            subplot('Position',marray);
            
            
            if chro==1 
                     bonetype=3;
                     cel=a21.data{1}.cel;
                     nuc=a21.data{3}.cel;
            elseif chro==2
                     bonetype=4;
                     cel=a21.data{2}.cel;
                     nuc=a21.data{4}.cel;
            elseif chro==3
                     bonetype=5;
                     cel=a31.data{1}.cel;
                     nuc=a31.data{3}.cel;
            elseif chro==4
                     cel=a31.data{2}.cel;
                     nuc=a31.data{4}.cel;
            elseif chro==5
                     cel=a32.data{1}.cel;
                     nuc=a32.data{3}.cel;
            elseif chro==6
                     cel=a32.data{2}.cel;
                     nuc=a32.data{4}.cel;
            end
           
         

            celavg=nanmean(cel,2);
            nucavg=nanmean(nuc,2);
            celstd=(nanstd(cel'))';
            nucstd=(nanstd(nuc'))';
           
            p(1)=plot(myinterval,celavg(1:profilesize-1),fcelallcolor{bonetype},'linewidth',0.5);% 'r'
            hold on 
            p(2)=plot(myinterval,nucavg(1:profilesize-1),fnucallcolor{bonetype},'linewidth',0.2);% 'b'

          

            
            
            index = 1:profilesize-1; min_y1=celavg-celstd; max_y1=celavg+celstd;        
            x=myinterval;stacky2=(min_y1)';stacky1=(max_y1)';
            fillxx=x(index([1:end end:-1:1 1]));fillyy=[stacky1(index) stacky2(index(end:-1:1)) stacky1(1)];
            h=fill(fillxx,fillyy,'r','EdgeColor','none','facealpha',0.3);
            
            
            min_y2=nucavg-nucstd; max_y2=nucavg+nucstd;  
            stacky2=(min_y2)';stacky1=(max_y2)';
            fillxx=x(index([1:end end:-1:1 1]));fillyy=[stacky1(index) stacky2(index(end:-1:1)) stacky1(1)];
            h=fill(fillxx,fillyy,'b','EdgeColor','none','facealpha',0.2);

         
       % legendarray2={'WT','MT'};
        set(gca,'fontsize',11)
%         legend(p,legendarray2,'location','southeast','fontsize',7);
%         legend 'boxoff';
        %legend(q,legendarray,'location','northwest','fontsize',7);
        ylabel(ylabels{i},'fontsize',14); 
        xlabel('Bone long axis','fontsize',14)
        if i==1
        title(tname{bonetype},'fontweight','normal','fontsize',14)
        end
        
        if i==1        
             gap=[0.25,0.7];
        elseif i==2
             gap=[0.05,0.4];
        else
             gap=[0.1,0.7];
        end
        
        axis([0,1.02,gap(1),gap(2)])

            
             if testwithstd==20
              boundary=statisticalTest(celavg,nucavg);
              value1=0.95*gap(2);
              howmanybound=0;
              for ii=1:2:length(boundary)
                  %disp('bone');
                  ind2=boundary(ii):boundary(ii+1);
                  if length(ind2)>3
                          [h,pR]=ttest2(celavg(ind2),nucavg(ind2)); %stat= mes(celavg(ind2),nucavg(ind2),'hedgesg' );
                       if pR<0.05 
                          if mod(howmanybound,2)==0
                              value=value1;
                          else
                              value=value1-0.0*diff(gap);
                          end
                          if pR<0.001
                              text(ind2(1)/50,value,strcat('p=',sprintf('%0.1e',pR)));
                          else
                              text(ind2(1)/50,value,strcat('p=',sprintf('%0.3f',pR)));
                          end
                          howmanybound=howmanybound+1;
              %text(ind2(1)/50,value-0.09*diff(gap),strcat('g=',sprintf('%0.3f',abs(stat.hedgesg))));
              plot([ind2(1)/50,ind2(1)/50],[0,110],'k:','linewidth',0.1);
              plot([ind2(end)/50,ind2(end)/50],[0,110],'k:','linewidth',0.1);
              %[ind2(1)/50,ind2(end)/50]
              %rectangle('Position',[ind2(1)/50,0, (ind2(end)-ind2(1))/50, 11 ],'FaceColor',[0 .5 .5],'EdgeColor','none','LineWidth',0.1,'facealpha',0.02);
              
              x=[ind2(1)/50,ind2(end)/50];index=[1,2];stacky1=[0,0];stacky2=[110,110];
              fillxx=x(index([1:end end:-1:1 1]));fillyy=[stacky1(index) stacky2(index(end:-1:1)) stacky1(1)];
              h=fill(fillxx,fillyy,'k','EdgeColor','none','facealpha',0.05);
                       end
                  end
              end
             end
             
             
               limit1=9;limit2=20; limit3=32; limit4=40; 
               b(1)=myinterval(limit1);
               b(2)=myinterval(limit2);
               b(3)=myinterval(limit3);
               b(4)=myinterval(limit4);
             
               
             plot([b(1),b(1)],[0,110],'k:','linewidth',0.1);
             plot([b(2),b(2)],[0,110],'k:','linewidth',0.1);
             plot([b(3),b(3)],[0,110],'k:','linewidth',0.1);
             plot([b(4),b(4)],[0,110],'k:','linewidth',0.1);
               
               index=1:limit1;
               out1=dofitting(index,myinterval,celavg,0.202,0.198,0.194,'g','g-');
               out2=dofitting(index,myinterval,nucavg,0.132,0.128,0.124,'k','k-');
               p(3)=out1.plot;  p(4)=out2.plot; 
               legend(p,lname,'location','southwest','fontsize',7);
               %legend 'boxoff';
               tvalue=abs(out1.p(1) - out2.p(1)) / sqrt( out1.SE^2 + out2.SE^2);
               mypvalue=2*tcdf(tvalue,length(celavg)+length(nucavg)-4, 'upper');
               text(myinterval(index(1)),0.95*gap(2), ['pvalue = ', sprintf('%0.4f',mypvalue)],'fontsize',6)


               index=1+limit1:limit2;
               out1=dofitting(index,myinterval,celavg,0.188,0.184,0.18,'g','g-');
               out2=dofitting(index,myinterval,nucavg,0.118,0.114,0.11,'k','k-');
               tvalue=abs(out1.p(1) - out2.p(1)) / sqrt( out1.SE^2 + out2.SE^2);
               mypvalue=2*tcdf(tvalue,length(celavg)+length(nucavg)-4, 'upper');
               text(myinterval(index(1)),0.92*gap(2), ['pvalue = ', sprintf('%0.4f',mypvalue)],'fontsize',6)

               
               index=limit2+1:limit3;
               out1=dofitting(index,myinterval,celavg,0.202,0.198,0.194,'g','g-');
               out2=dofitting(index,myinterval,nucavg,0.132,0.128,0.124,'k','k-');
               tvalue=abs(out1.p(1) - out2.p(1)) / sqrt( out1.SE^2 + out2.SE^2);
               mypvalue=2*tcdf(tvalue,length(celavg)+length(nucavg)-4, 'upper');
               text(myinterval(index(1)),0.95*gap(2), ['pvalue = ', sprintf('%0.4f',mypvalue)],'fontsize',6)

               
               
               index=limit3+1:limit4;
               out1=dofitting(index,myinterval,celavg,0.188,0.184,0.18,'g','g-');
               out2=dofitting(index,myinterval,nucavg,0.118,0.114,0.11,'k','k-');
               tvalue=abs(out1.p(1) - out2.p(1)) / sqrt( out1.SE^2 + out2.SE^2);
               mypvalue=2*tcdf(tvalue,length(celavg)+length(nucavg)-4, 'upper');
               text(myinterval(index(1)),0.92*gap(2), ['pvalue = ', sprintf('%0.4f',mypvalue)],'fontsize',6)


               index=limit4+1:50;
               out1=dofitting(index,myinterval,celavg,0.188,0.184,0.18,'g','g-');
               out2=dofitting(index,myinterval,nucavg,0.118,0.114,0.11,'k','k-');
               tvalue=abs(out1.p(1) - out2.p(1)) / sqrt( out1.SE^2 + out2.SE^2);
               mypvalue=2*tcdf(tvalue,length(celavg)+length(nucavg)-4, 'upper');
               text(myinterval(index(1)),0.95*gap(2), ['pvalue = ', sprintf('%0.4f',mypvalue)],'fontsize',6)

             
        
        
      end
	XPos=XPos+Width+XGap;
    end
    YPos=YPos-YGap-Height;
end





saveas(h1,['PCs_RatioMT']);
saveas(h1,['PCs_RatioMT.png']);



function   boundary=statisticalTest(celavg,nucavg)
              significant=[];
              for ii=1:46
                  ind2=ii:ii+4;
                  [h,pP]=ttest2(celavg(ind2),nucavg(ind2));
                  if pP<0.05
                      significant=[significant,ii];
                  end
              end    
              
              boundary=[significant(1)];
              for ii=2:length(significant)
                  if significant(ii)-significant(ii-1)~=1
                      boundary=[boundary,significant(ii-1),significant(ii)];
                  end
              end
              boundary=[boundary,significant(end)+4];
end

function output=dofitting(index,myinterval,celavg,pos1,pos2,pos3,mycolor,plotcolor)
               x=(myinterval(index))';
               y=celavg(index);
               output=errorAndPvalue(x,y); 
               q=plot(x,output.yfit,plotcolor,'linewidth',0.6);
%                text(x(1),pos1, ['Y = ',sprintf('%0.2f',output.p(1)),'X + ',sprintf('%0.2f',output.p(2))],'fontsize',5,'color',mycolor)
%                text(x(1),pos2, ['R^2 = ', sprintf('%0.2f',output.Rsq1)],'fontsize',5,'color',mycolor)
%                text(x(1),pos3, ['pvalue = ', sprintf('%0.2e',output.pvalue)],'fontsize',5,'color',mycolor)
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