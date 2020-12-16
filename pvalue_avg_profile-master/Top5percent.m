
clear all 

%load('avg_Bone_Data_Volume.mat')
%load('50avg/avgBoneDataVolume.mat')
load('avg_Bone_Data_Volume_80.mat')

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




profilesize=81;        
myinterval=linspace(0,1,profilesize);  
myinterval=myinterval(1:80);

myindex=62:72;

for chro=1:3
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
           

        if (chro==1)
            celavg=mean(cel(1:profilesize-1,:),2);   celstd=(std(cel(1:profilesize-1,:)'))';
            nucavg=mean(nuc(1:profilesize-1,:),2);   nucstd=(std(nuc(1:profilesize-1,:)'))';
%             celavg2=mean(cel(1:profilesize-1,[2:5]),2); celstd2=(std(cel(1:profilesize-1,[2:5])'))';
%             celavg(99)=celavg2(99); celstd(99)=celstd2(99);
            
            Hypertrophic_DT=[celavg,celstd];
          
        elseif (chro==3)
            celavg=mean(cel(1:profilesize-1,:),2);   celstd=(std(cel(1:profilesize-1,:)'))';
            nucavg=mean(nuc(1:profilesize-1,:),2);   nucstd=(std(nuc(1:profilesize-1,:)'))';
%             nucavg2=mean(nuc(1:profilesize-1,[2:3]),2); nucstd2=(std(nuc(1:profilesize-1,[2:3])'))';
%             nucavg(77)=nucavg2(77); nucstd(77)=nucstd2(77); 
            nucavg2=mean(nuc(1:profilesize-1,[2:3]),2); nucstd2=(std(nuc(1:profilesize-1,[2:3])'))';
            nucavg(62)=nucavg2(62); nucstd(62)=nucstd2(62);  
            Hypertrophic_PT=[celavg,celstd];
            Hypertrophic_DU=[nucavg,nucstd];
        end
end





h1=figure();
XL=0.15;XR=0.01;XGap=0.07;Row=1;
YT=0.06;YB=0.11;YGap=0.15;Col=1;
Width=(1-XL-XR-((Col-1)*XGap))/Col;
Height=(1-YT-YB-((Row-1)*YGap))/Row;
YPos=1-YT-Height; 

set(gcf, 'PaperSize', [5 4]);
set(gcf, 'PaperPosition', [0 0 5 4]);

topPercent=10;
tname={'DT wt', 'PT wt','DU wt'};

for i=1:Row
    XPos=XL;
    for j=1:Col
        chro=j+(i-1)*Col;
        if chro<=6
            marray=[XPos,YPos,Width,Height];
            subplot('Position',marray);
%             
%             [sa,sb]=sort(Hypertrophic_DT(:,1),'descend');
%              DT=sa(1:topPercent); DTind=sort(sb(1:topPercent)) 
%              
%             [sa,sb]=sort(Hypertrophic_PT(:,1),'descend');
%              PT=sa(1:topPercent); PTind=sort(sb(1:topPercent))
%         
%             [sa,sb]=sort(Hypertrophic_DU(:,1),'descend');
%              DU=sa(1:topPercent); DUind=sort(sb(1:topPercent))

            DTind=62:72; PTind=62:72; DUind=62:72;
            DT=Hypertrophic_DT(DTind,1);  DU=Hypertrophic_DU(DUind,1);
            PT=Hypertrophic_PT(PTind,1);

        
            X=[mean(DT),mean(PT),mean(DU)]; 
            %Y=[mean(Hypertrophic_DT(DTind,2)),mean(Hypertrophic_PT(PTind,2)),mean(Hypertrophic_DU(DUind,2))];
            Y=[std(DT), std(PT), std(DU)];
            hb=bar(X);
            hold on 
            pause(2); %pause allows the figure to be created
            for ib = 1:numel(hb)
                xData = hb(ib).XData+hb(ib).XOffset;
                errorbar(xData,X,Y,'k.')
            end           
    
        
            [h,pR]=ttest2(DT,PT);     text(0.6,8000,strcat('DT-PT(p=',sprintf('%0.4f',pR),')'  ));     
            [h,pR]=ttest2(DT,DU);     text(0.6,7500,strcat('DT-DU(p=',sprintf('%0.4f',pR),')'  ));   
            [h,pR]=ttest2(PT,DU);     text(0.6,7000,strcat('PT-DU(p=',sprintf('%0.4f',pR),')'  ));   
            
             ylim([0,8200])
%           
%               
%                 if chro>3
%                   gap=[0,2500];
%                   axis([0,1,gap(1),gap(2)]);
%                 elseif chro<=3
%                    gap=[0,9000];
%                    axis([0,1,gap(1),gap(2)])    
%                 else
%                         
%                 end
%                  title(titlename{bonetype},'fontweight','normal','fontsize',14)
%                  
%                 
              

            ylabel('Hypertrophic volume');
            set(gca,'xticklabel',tname);
            set(gca,'fontsize',11);
%             legend(p,tname{bonetype},'location','southeast','fontsize',7);
%             legend 'boxoff'; 
            %legend(q,legendarray,'location','northwest','fontsize',7);

        
    
      end
	XPos=XPos+Width+XGap;
    end
    YPos=YPos-YGap-Height;
end




saveas(h1,['Hypertrophic_Volume_Top_10_percent']);
saveas(h1,['Hypertrophic_Volume_Top_10_percent.png']);

% saveas(h1,['IndividualProfileSurfaceArea']);
% saveas(h1,['IndividualProfileSurfaceArea.png']);


