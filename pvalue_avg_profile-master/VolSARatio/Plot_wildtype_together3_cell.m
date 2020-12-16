
clear all 

%load('../CreateMatFiles/avg_Bone_Data_Volume.mat')
load('../CreateMatFiles/avg_Bone_Data_Volume2Sa.mat');
d=load('avg_Bone_GridDensity.mat');




fcelallcolor={'r','g','b','b','b'};
fnucallcolor={'r','g','b','b','b'};

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
XL=0.15;XR=0.02;XGap=0.07;Row=3;
YT=0.03;YB=0.05;YGap=0.03;Col=1;
Width=(1-XL-XR-((Col-1)*XGap))/Col;
Height=(1-YT-YB-((Row-1)*YGap))/Row;
YPos=1-YT-Height; 

set(gcf, 'PaperSize', [6 10]);
set(gcf, 'PaperPosition', [0 0 6 10]);

for i=1:Row
    XPos=XL;
    for j=1:Col
        chro=j+(i-1)*Col;
        if chro<=6
            marray=[XPos,YPos,Width,Height];
            subplot('Position',marray);
              
            
            legendarray2{1}=tname{1};
           legendarray2{2}=tname{2};
           legendarray2{3}=tname{5};
                
            count=1;
            for bonetype=[1,2,5]
                 clear nuc
                 clear cel 
                 if chro==1          
                    for k=1:size(data{bonetype}.celV,2)
                          for index=1:50                  
                              cel(index,k)=mean(data{bonetype}.celV{index,k});
                              nuc(index,k)=mean(data{bonetype}.nucV{index,k});                      
                          end
                    end
                elseif chro==2                 
                    for k=1:size(data{bonetype}.celV,2)
                          for index=1:50                  
                              cel(index,k)=mean(data{bonetype}.celS{index,k});
                              nuc(index,k)=mean(data{bonetype}.nucS{index,k});                      
                          end
                    end     
                 else
                        cel=d.data{bonetype}.cel(1:50,:);
                        nuc=d.data{bonetype}.nuc(1:50,:);
                 end
           

      

         
        celavg=nanmean(cel,2);
        nucavg=nanmean(nuc,2);        
        celstd=(nanstd(cel'))';
        nucstd=(nanstd(nuc'))';
        if (chro==3) &(bonetype==1)
                    tcel=cel(:,[2:5]); tcelstd=(nanstd(tcel'))';
                    celstd(50,:)=tcelstd(50,:);
                    celavg(50,:)=mean(tcel(50,:));
        end
              
        
        p(count)=plot(myinterval,celavg,strcat(fcelallcolor{bonetype},'-'),'linewidth',1.5); hold on
        
        index = 1:profilesize-1; min_y1=celavg-celstd; max_y1=celavg+celstd;        
        x=myinterval;stacky2=(min_y1)';stacky1=(max_y1)';
        fillxx=x(index([1:end end:-1:1 1]));fillyy=[stacky1(index) stacky2(index(end:-1:1)) stacky1(1)];
        h=fill(fillxx,fillyy,fcelallcolor{bonetype},'EdgeColor','none','facealpha',0.2);
        
        
        
%         p(count)=plot(myinterval,nucavg,strcat(fnucallcolor{bonetype},'-'),'linewidth',1.5); hold on 
%         index = 1:profilesize-1; x=myinterval;
%         min_y2=nucavg-nucstd; max_y2=nucavg+nucstd;  
%         stacky2=(min_y2)';stacky1=(max_y2)';
%         fillxx=x(index([1:end end:-1:1 1]));fillyy=[stacky1(index) stacky2(index(end:-1:1)) stacky1(1)];
%         h=fill(fillxx,fillyy,fnucallcolor{bonetype},'EdgeColor','none','facealpha',0.2);

              
    
        

            set(gca,'fontsize',11);               
            %     title(tname{bonetype},'fontweight','normal','fontsize',14)
                 count=count+1;
            end         

              
            if chro==1
                ylabel('Volume [\mum^3]');
                axis([0,1,0,7200])  
            end
            if chro==2
                ylabel('Surface Area [\mum^2]');
                axis([0,1,0,2700])  
            end
            if chro==3
                ylabel('Density [per \mum^3]');
                axis([0,1,0,3.5*10^-4]) 
            end
            if chro==3
                xlabel('Bone long axis')
            else
                set(gca,'xticklabel',[])
            end
            if chro==3
                 legend(p,legendarray2,'location','northeast','fontsize',7);
            else
                 legend(p,legendarray2,'location','northwest','fontsize',7);
            end
            legend 'boxoff'; 

        
    
      end
	XPos=XPos+Width+XGap;
    end
    YPos=YPos-YGap-Height;
end




saveas(h1,['wildType_Together_cell']);
saveas(h1,['wildType_Together_cell.png']);

% saveas(h1,['IndividualProfileSurfaceArea']);
% saveas(h1,['IndividualProfileSurfaceArea.png']);



