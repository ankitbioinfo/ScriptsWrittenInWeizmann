
clear all 

load('avgBoneDataSurfaceArea.mat')
           
%celallcolor={'r-','b-','g-','m-','k-'};
mycolor={'r','b','g','m','c'};
lt={'.','*','x','+','s'};


% fcelallcolor={'ro-','bo-','go-','mo-','co-'};
% fnucallcolor={'r^--','b^--','g^--','m^--','c^--'};

fcelallcolor={'ro','bo','go','mo','co'};
fnucallcolor={'r^','b^','g^','m^','c^'};

fcelallcolor={'ro','ro','ro','ro','ro'};
fnucallcolor={'b^','b^','b^','b^','b^'};

tname{1}={'S18 m6','S17 m2','S84 m3','S51 m2','S84 m4'};
tname{2}={'S18 m6','S17 m2','S84 m3','S51 m2','S84 m4'};
tname{3}={'S17 m1','S18 m2','S84 m1','S84 m5'};
tname{4}={'S17 m1','S18 m2','S84 m1','S84 m5'};      
tname{5}={'S51 m2','S84 m2','S84 m3'};

titlename={'DT wt', 'PT wt','DT mut','PT mut','DU wt'};



profilesize=50;        
myinterval=linspace(0,1,profilesize);  
globalindex=1:49;
myinterval=myinterval(globalindex);


h1=figure();
XL=0.07;XR=0.01;XGap=0.07;Row=1;
YT=0.06;YB=0.11;YGap=0.15;Col=1;
Width=(1-XL-XR-((Col-1)*XGap))/Col;
Height=(1-YT-YB-((Row-1)*YGap))/Row;
YPos=1-YT-Height; 

set(gcf, 'PaperSize', [13 7]);
set(gcf, 'PaperPosition', [0 0 13 7]);

clear q
clear legendarray 
count=1;

for i=1:Row
    XPos=XL;
    for j=1:Col
        chro=j+(i-1)*Col;
        if chro<6
            marray=[XPos,YPos,Width,Height];
            subplot('Position',marray);
            bonetype=chro;

            
            for bonetype=1:5

            cel=data{bonetype}.cel;
            nuc=data{bonetype}.nuc;
       
          
            for j=1:size(cel,2)
                    q(count)=plot(myinterval,cel(globalindex,j),strcat(mycolor{bonetype},lt{j},'-'),'linewidth',0.5);
                    legendarray{count}=strcat('cel:',titlename{bonetype}, ':',tname{bonetype}{j});
                    hold on 
                    q(21+count)=plot(myinterval,nuc(globalindex,j),strcat(mycolor{bonetype},lt{j},'--'),'linewidth',0.5);
                    legendarray{21+count}=strcat('nuc:',titlename{bonetype}, ':',tname{bonetype}{j});
                    count=count+1;
            end

            end
%         celavg=mean(cel(1:profilesize-1,:),2);
%         nucavg=mean(nuc(1:profilesize-1,:),2);
%         
%         celstd=(std(cel(1:profilesize-1,:)'))';
%         nucstd=(std(cel(1:profilesize-1,:)'))';
%         
%         p(1)=plot(myinterval,celavg,strcat(fcelallcolor{bonetype},'-'),'linewidth',0.2);
%         legendarray2{1}=strcat('C');%,tname{bonetype});
%         hold on 
%       
%         p(2)=plot(myinterval,nucavg,strcat(fnucallcolor{bonetype},'--'),'linewidth',0.2);
%         legendarray2{2}=strcat('N');%,tname{bonetype});
%         
%         
%         index = 1:profilesize-1; min_y1=celavg-celstd; max_y1=celavg+celstd;        
%         x=myinterval;stacky2=(min_y1)';stacky1=(max_y1)';
%         fillxx=x(index([1:end end:-1:1 1]));fillyy=[stacky1(index) stacky2(index(end:-1:1)) stacky1(1)];
%         h=fill(fillxx,fillyy,'r','EdgeColor','none','facealpha',0.2);
% 
%         min_y2=nucavg-nucstd; max_y2=nucavg+nucstd;  
%         stacky2=(min_y2)';stacky1=(max_y2)';
%         fillxx=x(index([1:end end:-1:1 1]));fillyy=[stacky1(index) stacky2(index(end:-1:1)) stacky1(1)];
%         h=fill(fillxx,fillyy,'b','EdgeColor','none','facealpha',0.2);
% 
%         
        
        
%         index=1:15;
%         [f,~]=fit(myinterval(index)',nucavg(index),g,'startpoint',[0.5,1,0.3,1000]);
%         plot(myinterval(index),f(myinterval(index)),strcat(mycolor{bonetype},'--'))
%         
%         [f,~]=fit(myinterval(index)',celavg(index),g,'startpoint',[0.1,1,0.3,1000]);
%         plot(myinterval(index),f(myinterval(index)),strcat(mycolor{bonetype},'-'))
        
%         index=21:39;
%         [f,~]=fit(myinterval(index)',nucavg(index),g,'startpoint',[10,5,100,1000])
%         plot(myinterval(index),f(myinterval(index)),'b.-')


%             legend(p,legendarray2,'location','northwest','fontsize',7);
%             legend 'boxoff'; 
            %legend(q,legendarray,'location','northwest','fontsize',7);

        
    
      end
	XPos=XPos+Width+XGap;
    end
    YPos=YPos-YGap-Height;
end



            set(gca,'fontsize',9);
            axis([0,1,5.2,9])
            %title(strcat(titlename{bonetype},', Cell(-), Nuclei(--)'),'fontweight','normal','fontsize',14)
          
            hl=legend(q,legendarray,'location','northwest','fontsize',5);
            legend 'boxoff'; 
            
            %hl.NumColumns=2;
            %set(legend, 'NumColumns' ,3)
          

            ylabel('ln(Surface area)');
            xlabel('Bone long axis')


saveas(h1,['Individual_Profile_SurfaceArea']);
saveas(h1,['Individual_Profile_SurfaceArea.png']);



