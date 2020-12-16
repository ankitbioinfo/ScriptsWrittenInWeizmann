
clear all 

pc1=load('avg_BoneData_lambda1.mat');
pc2=load('avg_BoneData_lambda2.mat');
pc3=load('avg_BoneData_lambda3.mat');



           
%celallcolor={'r-','b-','g-','m-','k-'};
mycolor={'r','b','g','m','c'};
lt={'.','*','x','+','s'};





fcelallcolor={'ro','ro','ro','ro','ro','ro'};
fnucallcolor={'b^','b^','b^','b^','b^','b^'};

            tname{1}={'wt', 'mut'};
            tname{2}={'wt', 'mut'};
            tname{3}={'wt', 'mut'};
            tname{4}={'wt' ,'mut'};
            tname{5}={'wt' ,'mut'};
            tname{6}={'wt' ,'mut'};
            
            
titlename={'PT Cell','DT Cell'};



profilesize=51;        
myinterval=linspace(0,1,profilesize-1);  
%myinterval=myinterval(1:49);



h1=figure();
XL=0.07;XR=0.03;XGap=0.07;Row=3;
YT=0.06;YB=0.08;YGap=0.07;Col=2;
Width=(1-XL-XR-((Col-1)*XGap))/Col;
Height=(1-YT-YB-((Row-1)*YGap))/Row;
YPos=1-YT-Height; 

set(gcf, 'PaperSize', [8 10]);
set(gcf, 'PaperPosition', [0 0 8 10]);

% 1 for std deviaton test and 2 for mean test 
testwithstd=2;



% 
% Bad sample  
% DT S84 m5 mut  3,4
% DT S84 m4 wt   1,5 
% PT S84 m5 mut  4,4
% DU S84 m3 wt   5,3

for i=1:Row
    XPos=XL;
    for j=1:Col
        chro=j+(i-1)*Col;
        if chro<=6
            marray=[XPos,YPos,Width,Height];
            subplot('Position',marray);
            bonetype=chro;
             
            
            if chro==1
                    cel=pc1.data{2}.cel;
                    nuc=pc1.data{4}.cel; 
            end            

            if chro==2                
                    cel=pc1.data{1}.cel; 
                    nuc=pc1.data{3}.cel;            
            end
            
            if chro==3       
                    cel=pc2.data{2}.cel;
                    nuc=pc2.data{4}.cel;                 
            end
            
            if chro==4              
                    cel=pc2.data{1}.cel; 
                    nuc=pc2.data{3}.cel;            
            end
            
            
            if chro==5
                    cel=pc3.data{2}.cel;
                    nuc=pc3.data{4}.cel;   
            end
            
            if chro==6                
                    cel=pc3.data{1}.cel; 
                    nuc=pc3.data{3}.cel;            
            end
        


        celavg=mean(cel(1:profilesize-1,:),2);
        nucavg=mean(nuc(1:profilesize-1,:),2);
        
        celstd=(std(cel(1:profilesize-1,:)'))';
        nucstd=(std(nuc(1:profilesize-1,:)'))';
        
        p(1)=plot(myinterval,celavg,strcat(fcelallcolor{bonetype},'-'),'linewidth',0.2);
        legendarray2{1}=tname{bonetype};
        hold on 
      
        p(2)=plot(myinterval,nucavg,strcat(fnucallcolor{bonetype},'--'),'linewidth',0.2);
        legendarray2{2}=tname{bonetype};
        
        index = 1:profilesize-1; min_y1=celavg-celstd; max_y1=celavg+celstd;        
        x=myinterval;stacky2=(min_y1)';stacky1=(max_y1)';
        fillxx=x(index([1:end end:-1:1 1]));fillyy=[stacky1(index) stacky2(index(end:-1:1)) stacky1(1)];
        h=fill(fillxx,fillyy,'r','EdgeColor','none','facealpha',0.3);

        min_y2=nucavg-nucstd; max_y2=nucavg+nucstd;  
        stacky2=(min_y2)';stacky1=(max_y2)';
        fillxx=x(index([1:end end:-1:1 1]));fillyy=[stacky1(index) stacky2(index(end:-1:1)) stacky1(1)];
        h=fill(fillxx,fillyy,'b','EdgeColor','none','facealpha',0.2);

              
    
        
%         index=1:15;
%         [f,~]=fit(myinterval(index)',nucavg(index),g,'startpoint',[0.5,1,0.3,1000]);
%         plot(myinterval(index),f(myinterval(index)),strcat(mycolor{bonetype},'--'))
%         
%         [f,~]=fit(myinterval(index)',celavg(index),g,'startpoint',[0.1,1,0.3,1000]);
%         plot(myinterval(index),f(myinterval(index)),strcat(mycolor{bonetype},'-'))
        
%         index=21:39;
%         [f,~]=fit(myinterval(index)',nucavg(index),g,'startpoint',[10,5,100,1000])
%         plot(myinterval(index),f(myinterval(index)),'b.-')

            set(gca,'fontsize',11);
            if i==1
                gap=[0,93];
            elseif i==2
                gap=[0,40];
            else
                gap=[0,20];          
            end
            axis([0,1,gap(1), gap(2)]);
            if i==1
            title(titlename{bonetype},'fontweight','normal','fontsize',14)
            end
            

            if testwithstd==1
              boundary=statisticalTest(celstd,nucstd);
              value1=90;
              howmanybound=0;
              for ii=1:2:length(boundary)
                  %disp('bone');
                  ind2=boundary(ii):boundary(ii+1);
                  if length(ind2)>3
                          [h,pR]=ttest2(celstd(ind2),nucstd(ind2)); %stat= mes(celavg(ind2),nucavg(ind2),'hedgesg' );
                       if pR<0.05 
                          if mod(howmanybound,2)==0
                              value=value1;
                          else
                              value=value1-0.09*diff(gap);
                          end
                          if pR<0.001
                              text(ind2(1)/50,value,strcat('p=',sprintf('%0.1e',pR)));
                          else
                              text(ind2(1)/50,value,strcat('p=',sprintf('%0.3f',pR)));
                          end
                          howmanybound=howmanybound+1;
              %text(ind2(1)/50,value-0.09*diff(gap),strcat('g=',sprintf('%0.3f',abs(stat.hedgesg))));
              plot([ind2(1)/50,ind2(1)/50],[0,110],'k:','linewidth',0.5);
              plot([ind2(end)/50,ind2(end)/50],[0,110],'k:','linewidth',0.5);
              %[ind2(1)/50,ind2(end)/50]
              %rectangle('Position',[ind2(1)/50,0, (ind2(end)-ind2(1))/50, 11 ],'FaceColor',[0 .5 .5],'EdgeColor','none','LineWidth',0.1,'facealpha',0.02);
              
              x=[ind2(1)/50,ind2(end)/50];index=[1,2];stacky1=[0,0];stacky2=[110,110];
              fillxx=x(index([1:end end:-1:1 1]));fillyy=[stacky1(index) stacky2(index(end:-1:1)) stacky1(1)];
              h=fill(fillxx,fillyy,'k','EdgeColor','none','facealpha',0.05);
                       end
                  end
              end
            end
              
              
             if testwithstd==2  
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
                              value=value1-0.09*diff(gap);
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
              
            
          
            if i==1
                ylabel('PC1 coefficient');
            elseif i==2
                ylabel('PC2 coefficient');
            else
                ylabel('PC3 coefficient');
            end
            xlabel('Bone long axis')
            legend(p,tname{bonetype},'location','west','fontsize',7);
            legend 'boxoff'; 
            %legend(q,legendarray,'location','northwest','fontsize',7);

        
    
      end
	XPos=XPos+Width+XGap;
    end
    YPos=YPos-YGap-Height;
end


if testwithstd==1
saveas(h1,['Cells_PCs_test_stddev']);
saveas(h1,['Cells_PCs_test_stddev.png']);
end

if testwithstd==2
saveas(h1,['Cells_PC_lambda_WT_MT_test_mean']);
saveas(h1,['Cells_PC_lambda_WT_MT_test_mean.png']);
end


% saveas(h1,['IndividualProfileSurfaceArea']);
% saveas(h1,['IndividualProfileSurfaceArea.png']);

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
