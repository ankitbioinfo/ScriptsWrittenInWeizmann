
clear all 

load('avgBoneDataVolume.mat')
           
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

%g=fittype(@(alpha,beta,constant,V0,x) (beta/alpha + constant/V0*exp(-alpha*x))    );



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
           
        
        
            
            

             set(gca,'fontsize',11);
            axis([0,1,5.2,9])
            %title(strcat(titlename{bonetype},', Cell(-), Nuclei(--)'),'fontweight','normal','fontsize',14)
          
            legend(q,legendarray,'location','northwest','fontsize',5);
            legend 'boxoff';          

            ylabel('ln(Volume)');
            xlabel('Bone long axis')
%             legend(p,legendarray2,'location','northwest','fontsize',7);
%             legend 'boxoff'; 
            %legend(q,legendarray,'location','northwest','fontsize',7);

        
    
      end
	XPos=XPos+Width+XGap;
    end
    YPos=YPos-YGap-Height;
end






saveas(h1,['Individual_Profile_Volume']);
saveas(h1,['Individual_Profile_Volume.png']);



