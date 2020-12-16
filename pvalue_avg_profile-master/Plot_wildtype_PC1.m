
clear all 

load('local_avgBonePC1.mat')

%celallcolor={'r-','b-','g-','m-','k-'};
mycolor={'r','b','g','m','c'};
lt={'.','*','x','+','s'};


% fcelallcolor={'ro-','bo-','go-','mo-','co-'};
% fnucallcolor={'r^--','b^--','g^--','m^--','c^--'};

fcelallcolor={'ro','bo','go','mo','co'};
fnucallcolor={'r^','b^','g^','m^','c^'};


fcelallcolor={'ro','ro','ro','ro','ro','ro'};
fnucallcolor={'b^','b^','b^','b^','b^','b^'};

            tname{1}={'DT wt', 'PT wt'};
            tname{2}={'DT wt', 'DU wt'};
            tname{3}={'PT wt', 'DU wt'};
            tname{4}={'DT wt', 'PT wt'};
            tname{5}={'DT wt', 'DU wt'};
            tname{6}={'PT wt', 'DU wt'};
            titlename={'Cell', 'Cell','Cell','Nuclei','Nuclei','Nuclei'};



profilesize=50;        
myinterval=linspace(0,1,profilesize);  
myinterval=myinterval(1:49);


h1=figure();
XL=0.07;XR=0.03;XGap=0.07;Row=2;
YT=0.06;YB=0.11;YGap=0.15;Col=3;
Width=(1-XL-XR-((Col-1)*XGap))/Col;
Height=(1-YT-YB-((Row-1)*YGap))/Row;
YPos=1-YT-Height; 

set(gcf, 'PaperSize', [13 7]);
set(gcf, 'PaperPosition', [0 0 13 7]);

for i=1:Row
    XPos=XL;
    for j=1:Col
        chro=j+(i-1)*Col;
        if chro<=6
            marray=[XPos,YPos,Width,Height];
            subplot('Position',marray);
            bonetype=chro;
            
          

            if chro==1                
                    cel=data{1}.cel;
                    nuc=data{2}.cel;               
            end
            if chro==2   
                    cel=data{1}.cel;
                    nuc=data{5}.cel;
            end
            
            if chro==3               
                    cel=data{2}.cel;
                    nuc=data{5}.cel;               
            end
            if chro==4   
                    cel=data{1}.nuc;
                    nuc=data{2}.nuc;
            end
        
            if chro==5
                  cel=data{1}.nuc;
                  nuc=data{5}.nuc; 
            end
            
            if chro==6               
                    cel=data{2}.nuc;
                    nuc=data{5}.nuc;               
            end

      
            
            
        if (chro==5) | (chro==6)
             celavg=mean(cel(1:profilesize-1,:),2);          celstd=(std(cel(1:profilesize-1,:)'))';
             nucavg=mean(nuc(1:profilesize-1,:),2);          nucstd=(std(nuc(1:profilesize-1,:)'))';    
             nucavg2=mean(nuc(1:profilesize-1,[1,3]),2);     nucstd2=(std(nuc(1:profilesize-1,[1,3])'))';
             nucavg(2:3)=nucavg2(2:3);  nucstd(2:3)=nucstd2(2:3);
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
        
        ankit(chro,:)=fillyy;

        min_y2=nucavg-nucstd; max_y2=nucavg+nucstd;  
        stacky2=(min_y2)';stacky1=(max_y2)';
        fillxx=x(index([1:end end:-1:1 1]));fillyy=[stacky1(index) stacky2(index(end:-1:1)) stacky1(1)];
        h=fill(fillxx,fillyy,'b','EdgeColor','none','facealpha',0.2);

              
        ankur(chro,:)=fillyy;
    
        
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
          
             
          
              gap=[0,93];
              axis([0,1,0,93])
                 title(titlename{bonetype},'fontweight','normal','fontsize',14)
                 
            
                  significant=[];
              for ii=1:45
                  ind2=ii:ii+4;
                  [h,pP]=ttest2(celavg(ind2),nucavg(ind2));
                  ankit(ii)=pP;
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
              
              value=90;
              for ii=1:2:length(boundary)
                  %disp('bone');
                  ind2=boundary(ii):boundary(ii+1);
                  if length(ind2)>3
              [h,pR]=ttest2(celavg(ind2),nucavg(ind2)); stat= mes(celavg(ind2),nucavg(ind2),'hedgesg' );
              
              if length(ind2)>30
                  pos=mean(ind2)/50;
              else
                  pos=ind2(1)/50;
              end
              
              if pR<0.001
                  text(pos,value,strcat('p=',sprintf('%0.1e',pR)));
              else
                  text(pos,value,strcat('p=',sprintf('%0.3f',pR)));
              end

              text(pos,value-0.09*diff(gap),strcat('g=',sprintf('%0.3f',abs(stat.hedgesg))));
              plot([ind2(1)/50,ind2(1)/50],[0,110],'k:','linewidth',0.5);
              plot([ind2(end)/50,ind2(end)/50],[0,110],'k:','linewidth',0.5);
              %[ind2(1)/50,ind2(end)/50]
              %rectangle('Position',[ind2(1)/50,0, (ind2(end)-ind2(1))/50, 11 ],'FaceColor',[0 .5 .5],'EdgeColor','none','LineWidth',0.1,'facealpha',0.02);
              
              x=[ind2(1)/50,ind2(end)/50];index=[1,2];stacky1=[0,0];stacky2=[110,110];
              fillxx=x(index([1:end end:-1:1 1]));fillyy=[stacky1(index) stacky2(index(end:-1:1)) stacky1(1)];
              h=fill(fillxx,fillyy,'k','EdgeColor','none','facealpha',0.08);
     
                  end
              end
              

            ylabel('PC1 orientation');
            xlabel('Bone long axis')
            legend(p,tname{bonetype},'location','southwest','fontsize',7);
            legend 'boxoff'; 
            %legend(q,legendarray,'location','northwest','fontsize',7);

        
    
      end
	XPos=XPos+Width+XGap;
    end
    YPos=YPos-YGap-Height;
end




saveas(h1,['wildType_PC1']);
saveas(h1,['wildType_PC1.png']);

% saveas(h1,['IndividualProfileSurfaceArea']);
% saveas(h1,['IndividualProfileSurfaceArea.png']);



