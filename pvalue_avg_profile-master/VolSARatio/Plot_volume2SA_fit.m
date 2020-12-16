
clear all 

load('../CreateMatFiles/avg_individual_vol2sa100.mat')
profilesize=101;        


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
tname={'DT wt', 'PT wt','DT mut','PT mut','DU wt'};

g=fittype(@(alpha,beta,constant,A0,x) (alpha/beta + constant/A0*exp(-beta*x))    );
logistic=fittype(@(a,b,c,d,x) (d+b./(c+exp(-a*x))));

g=fittype(@(m,c,x) (m*x+c));




scaling=1;
myinterval=linspace(0,scaling,profilesize-1);  

% 1 for std deviaton test and 2 for mean test 
testwithstd=2;


h1=figure();
XL=0.05;XR=0.02;XGap=0.05;Row=1;
YT=0.06;YB=0.12;YGap=0.15;Col=3;
Width=(1-XL-XR-((Col-1)*XGap))/Col;
Height=(1-YT-YB-((Row-1)*YGap))/Row;
YPos=1-YT-Height; 

set(gcf, 'PaperSize', [13 4]); %7
set(gcf, 'PaperPosition', [0 0 13 4]);

for i=1:Row
    XPos=XL;
    for j=1:Col
        chro=j+(i-1)*Col;
        if chro<6
            marray=[XPos,YPos,Width,Height];
            subplot('Position',marray);
            
            
            if chro==1 
                     bonetype=chro;
            elseif chro==2
                     bonetype=2;
            elseif chro==3
                     bonetype=5;
            end
                
       
            cel=data{bonetype}.cel;
            nuc=data{bonetype}.nuc;

           
       
            %         for j=1:size(cel,2)
            %                 q(count)=plot(myinterval,cel(:,j),strcat(mycolor{bonetype},lt{j},'-'),'linewidth',0.5);
            %                 legendarray{count}=strcat('cel:',tname{count});
            %                 hold on 
            %                 q(length(tname)+count)=plot(myinterval,nuc(:,j),strcat(mycolor{bonetype},lt{j},'--'),'linewidth',0.5);
            %                 legendarray{length(tname)+count}=strcat('nuc:',tname{count});
            %                 count=count+1;
            %         end

            if chro==3
                %celavg=nanmean(cel(:,[1,3]),2);
                celavg=nanmean(cel,2);
            else            
            celavg=nanmean(cel,2);
            end
           % nucavg=nanmean(nuc,2);
            celstd=(nanstd(cel'))';
            %nucstd=(nanstd(nuc'))';
            %dlmwrite(['ankit',num2str(chro),'.dat'],[myinterval', celavg(1:100),celstd(1:100)],'\t')
            dlmwrite(['ankit',num2str(chro),'.dat'],[myinterval', cel(1:100,:)],'\t')

            
          
        


            p(1)=plot(myinterval,celavg(1:profilesize-1),fcelallcolor{bonetype},'linewidth',0.5);% 'r'
            hold on 
            %p(2)=plot(myinterval,nucavg,fnucallcolor{bonetype},'linewidth',0.2);% 'b'

%             d=celavg>nucavg;  limit=round(profilesize/2);   d1=max(find(d(1:limit))) ; d2=limit+min(find(d((limit+1):(profilesize-1))));
%             plot([myinterval(d1),myinterval(d1)],[0,3],'k:','linewidth',0.2)
%             hold on 
%             plot([myinterval(d2),myinterval(d2)],[0,3],'k:','linewidth',0.2)
            
          
            index = 1:profilesize-1; min_y1=celavg-celstd; max_y1=celavg+celstd;        
            x=myinterval;stacky2=(min_y1)';stacky1=(max_y1)';
            fillxx=x(index([1:end end:-1:1 1]));fillyy=[stacky1(index) stacky2(index(end:-1:1)) stacky1(1)];
            h=fill(fillxx,fillyy,'r','EdgeColor','none','facealpha',0.3);
            
            
%             min_y2=nucavg-nucstd; max_y2=nucavg+nucstd;  
%             stacky2=(min_y2)';stacky1=(max_y2)';
%             fillxx=x(index([1:end end:-1:1 1]));fillyy=[stacky1(index) stacky2(index(end:-1:1)) stacky1(1)];
%             h=fill(fillxx,fillyy,'b','EdgeColor','none','facealpha',0.2);
% 

        


            index=1:15;
            %[f,~]=fit(myinterval(index)',celavg(index),g,'startpoint',[5,10,1,500]);
            %plot(myinterval(index),f(myinterval(index)),strcat('r','-'),'linewidth',2)
            legendarray2{1}=strcat('C');
           % ctext=strcat('\color{red}','\alpha =',sprintf('%0.1f',f.alpha),', \beta =',sprintf('%0.1f',f.beta), ', \alpha/\beta =',sprintf('%0.1f',f.alpha/f.beta) );

            index=1:25;
            %[f,~]=fit(myinterval(index)',nucavg(index),g,'startpoint',[5,10,1,200]);
           % plot(myinterval(index),f(myinterval(index)),strcat('b','-'),'linewidth',2)
            legendarray2{2}=strcat('N');
            %ntext=strcat('\color{blue}','\alpha =',sprintf('%0.1f',f.alpha),', \beta =',sprintf('%0.1f',f.beta),', \alpha/\beta =',sprintf('%0.1f',f.alpha/f.beta) );


        set(gca,'fontsize',11)
%         legend(p,legendarray2,'location','southeast','fontsize',7);
%         legend 'boxoff';
        %legend(q,legendarray,'location','northwest','fontsize',7);
        ylabel('Vol/SA','fontsize',14); 
        xlabel('Bone long axis','fontsize',14)
        title(tname{bonetype},'fontweight','normal','fontsize',14)
        axis([0,1*scaling,1,2.9])
        
        gap=[1,2.9];
        
             if chro==1
             windowLen=11;  % 6    11 
             elseif chro==2
                 windowLen=4;  %3   4 
             else 
                 windowLen=7;  %5   7
             end
             for ii=1:profilesize-1-windowLen
                  ind2=ii:ii+ windowLen;
                  x=(myinterval(ind2))';
                  y=celavg(ind2);
                 
                  p=polyfit(x,y,1);
                  yfit = polyval(p,x);
                  slope(ii,1)=p(1);
                  boundary(ii,:)=[min(x),max(x)];
                  yresid = y-yfit;
                  SSresid = sum(yresid.^2);
                  SStotal = (length(y)-1) * var(y);
                  rsq(ii,1) = 1 - SSresid/SStotal; 
                  rsq_adj(ii,1) = 1 - SSresid/SStotal * (length(y)-1)/(length(y)-length(p));
                 
%                   if  rsq(ii,1)<0.9
%                     %if ((0.3<min(x)) & (0.5>max(x)))
%                     if p(1)>0
%                         plot(x,yfit,'k-','linewidth',3);
%                     else
%                         %plot(x,yfit,'b-','linewidth',3);
%                     end 
%                   end
             end
             ankit{chro}=rsq;
           
             index=[];
             region=[];
             for k=1:length(slope)
                 if (slope(k)<0) & (myinterval(k)<0.7*scaling)
                     region=[region; boundary(k,:)];
                     index=[index,k:k+windowLen];
                 end
             end
             b(1)=min(region(:,1));
             b(2)=max(region(:,2));
             
             
              plot([b(1),b(1)],[0,110],'k:','linewidth',0.1);
              plot([b(2),b(2)],[0,110],'k:','linewidth',0.1);
              
             
               index=unique(index)
               x=(myinterval(index))';y=celavg(index);
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
               [f.m/SE,stat.dfe,mypvalue]

               
               plot(x,yfit,'b-','linewidth',0.5);
               text(x(1),2.8, ['Y = ',sprintf('%0.2f',p(1)),'X + ',sprintf('%0.2f',p(2))])
               text(x(1),2.7, ['R^2 = ', sprintf('%0.2f',nrsq)])
               text(x(1),2.6, ['pvalue = ', sprintf('%0.2e',mypvalue)])

               
               ax1=gca;
               ax2 = axes('Position',[XPos+0.3*x(1),YPos+0.45, 0.25*(x(end)-x(1)), .15]);
                box on;
                plot(myinterval,celavg(1:profilesize-1),fcelallcolor{bonetype},'linewidth',0.5);% 'r'
                hold on 
                h=fill(fillxx,fillyy,'r','EdgeColor','none','facealpha',0.3);
                plot(x,yfit,'b-','linewidth',0.5);
                errorbar(x,yfit,d,'color','b','LineStyle', 'None');
                axis([x(1),x(end),1.25,1.5])
                set(gca,'xtick',[])
                %hold(ax1, 'on');
               
              
              
             if testwithstd==50  
              boundary=statisticalTest(celavg,yfit);
              value1=2.85;
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
              
        
        
        
        
        
        
        
        
        
      
        
        
        
      end
	XPos=XPos+Width+XGap;
    end
    YPos=YPos-YGap-Height;
end



% beta=1;
% alpha=1.2;
% constant=1;
% V0=3;
% t=linspace(0.5,1,100);
% sa_by_v= beta/alpha + constant/V0*exp(-alpha*t);
% plot(t,0.4*sa_by_v,'m.-')



if testwithstd==1
saveas(h1,['Volume2SurfaceArea_test_stddev']);
saveas(h1,['Volume2SurfaceArea_test_stddev.png']);
end

if testwithstd==2
saveas(h1,['Volume2SurfaceArea_test_cell']);
saveas(h1,['Volume2SurfaceArea_test_cell.png']);
end


% saveas(h1,['Volume2SurfaceArea']);
% saveas(h1,['Volume2SurfaceArea.png']);


function   boundary=statisticalTest(celavg,linearfit)
              significant=[];
              for ii=1:46
                  ind2=ii:ii+4;
                  [h,pP]=ttest2(celavg(ind2),linearfit{ii});
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
