
clear all 

load('avgBoneDataDelaunay.mat')

%celallcolor={'r-','b-','g-','m-','k-'};
mycolor={'r','b','g','m','c'};
lt={'.','*','x','+','s'};


% fcelallcolor={'ro-','bo-','go-','mo-','co-'};
% fnucallcolor={'r^--','b^--','g^--','m^--','c^--'};

fcelallcolor={'ro','bo','go','mo','co'};
fnucallcolor={'r^','b^','g^','m^','c^'};


fcelallcolor={'ro','ro','ro','ro','ro'};
fnucallcolor={'b^','b^','b^','b^','b^'};

            tname{1}={'wt', 'mut'};
            tname{2}={'wt', 'mut'};
            tname{3}={'wt', 'mut'};
            tname{4}={'wt' ,'mut'};
            %tname{3}={'Cell','Nuclei'};
            titlename={'DT Nuclei','PT Nuclei'};



profilesize=50;        
myinterval=linspace(0,1,profilesize);  
myinterval=myinterval(1:49);


h1=figure();
XL=0.07;XR=0.01;XGap=0.07;Row=1;
YT=0.07;YB=0.12;YGap=0.15;Col=2;
Width=(1-XL-XR-((Col-1)*XGap))/Col;
Height=(1-YT-YB-((Row-1)*YGap))/Row;
YPos=1-YT-Height; 

set(gcf, 'PaperSize', [9 4]);
set(gcf, 'PaperPosition', [0 0 9 4]);

for i=1:Row
    XPos=XL;
    for j=1:Col
        chro=j+(i-1)*Col;
        if chro<6
            marray=[XPos,YPos,Width,Height];
            subplot('Position',marray);
            bonetype=chro;
            
          

         
            
            if chro==1            
                    cel=data{1}.nuc;
                    nuc=data{3}.nuc;               
            end
            if chro==2  
                    cel=data{2}.nuc;
                    nuc=data{4}.nuc;
            end
        
%             if chro==3 
%                   cel=data{5}.cel;
%                   nuc=data{5}.nuc; 
%             end
            
        if chro==1
            celavg=mean(cel(1:profilesize-1,:),2); celstd=(std(cel(1:profilesize-1,:)'))';
            nucavg2=mean(nuc(1:profilesize-1,[1,3]),2);     nucstd2=(std(cel(1:profilesize-1,[1,3])'))';
            nucavg=mean(nuc(1:profilesize-1,:),2); nucstd=(std(cel(1:profilesize-1,:)'))';       
            nucavg(2:4)=nucavg2(2:4);  nucstd(2:4)=nucstd2(2:4);
        else
            celavg=mean(cel(1:profilesize-1,:),2);
            nucavg=mean(nuc(1:profilesize-1,:),2);
            celstd=(std(cel(1:profilesize-1,:)'))';
            nucstd=(std(nuc(1:profilesize-1,:)'))';
        end
        
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
     
       
    

            set(gca,'fontsize',11);
            gap=[0,1.4*10^-3];
            axis([0,1,gap(1), gap(2)]);

           
            title(titlename{bonetype},'fontweight','normal','fontsize',14)
            
            
            
              significant=[];
              for ii=1:45
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
              
              value=0.88*gap(2);
              for ii=1:2:length(boundary)
                  ind2=boundary(ii):boundary(ii+1);
              [h,pR]=ttest2(celavg(ind2),nucavg(ind2)); stat= mes(celavg(ind2),nucavg(ind2),'hedgesg' );
              if pR<0.001
                  text(ind2(1)/50,value+0.07*diff(gap),strcat('p=',sprintf('%0.1e',pR)));
              else
                  text(ind2(1)/50,value+0.07*diff(gap),strcat('p=',sprintf('%0.3f',pR)));
              end

              text(ind2(1)/50,value,strcat('g=',sprintf('%0.3f',abs(stat.hedgesg))));
              plot([ind2(1)/50,ind2(1)/50],[0,11],'k:','linewidth',0.5);
              plot([ind2(end)/50,ind2(end)/50],[0,11],'k:','linewidth',0.5);
              %[ind2(1)/50,ind2(end)/50]
              %rectangle('Position',[ind2(1)/50,0, (ind2(end)-ind2(1))/50, 11 ],'FaceColor',[0 .5 .5],'EdgeColor','none','LineWidth',0.1,'facealpha',0.02);
              
              x=[ind2(1)/50,ind2(end)/50];index=[1,2];stacky1=[0,0];stacky2=[11,11];
              fillxx=x(index([1:end end:-1:1 1]));fillyy=[stacky1(index) stacky2(index(end:-1:1)) stacky1(1)];
              h=fill(fillxx,fillyy,'k','EdgeColor','none','facealpha',0.02);
     
              
              end
            
              
            
            
          

            ylabel('Delaunay density');
            xlabel('Bone long axis')
            legend(p,tname{bonetype},'location','southwest','fontsize',7);
            legend 'boxoff'; 
            %legend(q,legendarray,'location','northwest','fontsize',7);

        
    
      end
	XPos=XPos+Width+XGap;
    end
    YPos=YPos-YGap-Height;
end





saveas(h1,['AverageProfileDelaunay']);
saveas(h1,['AverageProfileDelaunay.png']);

% saveas(h1,['IndividualProfileSurfaceArea']);
% saveas(h1,['IndividualProfileSurfaceArea.png']);



