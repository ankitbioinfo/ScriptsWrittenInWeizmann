
clear all 

%load('../CreateMatFiles/avg_Bone_Data_Volume.mat')
load('../CreateMatFiles/avg_Bone_Data_Volume2Sa.mat')



%celallcolor={'r-','b-','g-','m-','k-'};
mycolor={'r','b','g','m','c'};
lt={'.','*','x','+','s'};


% fcelallcolor={'ro-','bo-','go-','mo-','co-'};
% fnucallcolor={'r^--','b^--','g^--','m^--','c^--'};

fcelallcolor={'ro','bo','go','mo','co'};
fnucallcolor={'r^','b^','g^','m^','c^'};


fcelallcolor={'r','r','r','r','r','ro'};
fnucallcolor={'b','b','b','b','b','b^'};

            tname{1}={'DT wt', 'PT wt'};
            tname{2}={'DT wt', 'DU wt'};
            tname{3}={'PT wt', 'DU wt'};
            tname{4}={'DT wt', 'PT wt'};
            tname{5}={'DT wt', 'DU wt'};
            tname{6}={'PT wt', 'DU wt'};
            titlename={'Cell', 'Cell','Cell','Nuclei','Nuclei','Nuclei'};

tname={'DT wt', 'PT wt','DT mut','PT mut','DU wt'};

profilesize=51;        
myinterval=linspace(0,1,profilesize);  
myinterval=myinterval(1:50);


h1=figure();
XL=0.07;XR=0.01;XGap=0.07;Row=1;
YT=0.06;YB=0.11;YGap=0.15;Col=3;
Width=(1-XL-XR-((Col-1)*XGap))/Col;
Height=(1-YT-YB-((Row-1)*YGap))/Row;
YPos=1-YT-Height; 

set(gcf, 'PaperSize', [13 4]);
set(gcf, 'PaperPosition', [0 0 13 4]);

for i=1:Row
    XPos=XL;
    for j=1:Col
        chro=j+(i-1)*Col;
        if chro<=6
            marray=[XPos,YPos,Width,Height];
            subplot('Position',marray);
            bonetype=chro;
            
          

           if chro==1 
                     bonetype=chro;
            elseif chro==2
                     bonetype=3;
            elseif chro==3
                     bonetype=4;
           end           
            %        celV: {51?5 cell}
            %        nucV: {51?5 cell}
            %        celS: {51?5 cell}
            %        nucS: {51?5 cell}
            %     celV_Sa: {51?5 cell}
            %     nucV_Sa: {51?5 cell}
                
            clear nuc
            clear cel 
            for k=1:size(data{bonetype}.celV,2)
                  for index=1:50                  
                      cel(index,k)=mean(data{bonetype}.celS{index,k});
                      nuc(index,k)=mean(data{bonetype}.nucS{index,k});                      
                  end
            end     
           
           

        celavg=nanmean(cel,2);
        nucavg=nanmean(nuc,2);        
        celstd=(nanstd(cel'))';
        nucstd=(nanstd(nuc'))';
              
        
        p(1)=plot(myinterval,celavg,strcat(fcelallcolor{bonetype},'-'),'linewidth',1.5);
        legendarray2{1}=strcat('C');
        hold on 
      
        p(2)=plot(myinterval,nucavg,strcat(fnucallcolor{bonetype},'-'),'linewidth',1.5);
        legendarray2{2}=strcat('N');
        
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
%                 if chro>3
%                   gap=[4.8,7.5];
%                   axis([0,1,gap(1),gap(2)]);
%                 elseif chro<=3
%                    gap=[4.8,9];
%                    axis([0,1,gap(1),gap(2)])    
%                 else
%                    axis([0,1,5.2,9])           
%                 end
                 title(tname{bonetype},'fontweight','normal','fontsize',14)
                 
%                  if chro~=13
%                      ind=1:20;
%               [h,pR]=ttest2(celavg(ind),nucavg(ind)); stat= mes(celavg(ind),nucavg(ind),'hedgesg' );
%               text(0.01,5.0+0.1*diff(gap),strcat('p=',sprintf('%0.3f',pR)));
%               text(0.01,5.0,strcat('g=',sprintf('%0.3f',abs(stat.hedgesg))));
%               
%               ind2=21:42;
%               [h,pP]=ttest2(celavg(ind2),nucavg(ind2)); stat= mes(celavg(ind2),nucavg(ind2),'hedgesg' );
%               text(ind2(1)/50,5.0+0.1*diff(gap),strcat('p=',sprintf('%0.3f',pP)));
%               text(ind2(1)/50,5.0,strcat('g=',sprintf('%0.3f',abs(stat.hedgesg))));
%               
%               ind3=43:49;
%               [h,pH]=ttest2(celavg(ind3),nucavg(ind3));stat= mes(celavg(ind3),nucavg(ind3),'hedgesg' );
%               text(ind3(1)/50,5.0+0.1*diff(gap),strcat('p=',sprintf('%0.3f',pH)));
%               text(ind3(1)/50,5.0,strcat('g=',sprintf('%0.3f',abs(stat.hedgesg))));
% 
%               plot([ind2(1)/50,ind2(1)/50],[0,11],'k--');
%               plot([ind3(1)/50,ind3(1)/50],[0,11],'k--');
%                  end
              
              

            ylabel('Surface Area');
            xlabel('Bone long axis')
            legend(p,legendarray2,'location','northwest','fontsize',7);
            legend 'boxoff'; 
            %legend(q,legendarray,'location','northwest','fontsize',7);

        
    
      end
	XPos=XPos+Width+XGap;
    end
    YPos=YPos-YGap-Height;
end




%saveas(h1,['wildType_SurfArea']);
%saveas(h1,['wildType_SurfArea.png']);

% saveas(h1,['IndividualProfileSurfaceArea']);
% saveas(h1,['IndividualProfileSurfaceArea.png']);



