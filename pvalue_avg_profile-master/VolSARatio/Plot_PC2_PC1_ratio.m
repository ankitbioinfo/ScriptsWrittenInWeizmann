
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
tname={'DT wt', 'PT wt','DT mut','PT mut','DU wt'};

g=fittype(@(alpha,beta,constant,A0,x) (alpha/beta + constant/A0*exp(-beta*x))    );
logistic=fittype(@(a,b,c,d,x) (d+b./(c+exp(-a*x))));





profilesize=51;        
myinterval=linspace(0,1,profilesize);  
myinterval1=myinterval(1:50);


h1=figure();
XL=0.07;XR=0.01;XGap=0.07;Row=1;
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
                
        
             cel21=a21.data{bonetype}.cel;
             cel31=a31.data{bonetype}.cel;
             cel32=a32.data{bonetype}.cel;
            
            
%            cel21=a21.data{bonetype}.nuc;
%            cel31=a31.data{bonetype}.nuc;
%            cel32=a32.data{bonetype}.nuc;
            
            
            
            %nuc=data{bonetype}.nuc;

          
       
            %         for j=1:size(cel,2)
            %                 q(count)=plot(myinterval,cel(:,j),strcat(mycolor{bonetype},lt{j},'-'),'linewidth',0.5);
            %                 legendarray{count}=strcat('cel:',tname{count});
            %                 hold on 
            %                 q(length(tname)+count)=plot(myinterval,nuc(:,j),strcat(mycolor{bonetype},lt{j},'--'),'linewidth',0.5);
            %                 legendarray{length(tname)+count}=strcat('nuc:',tname{count});
            %                 count=count+1;
            %         end

            celavg=nanmean(cel21,2);
            nucavg=nanmean(cel31,2);
            objavg=nanmean(cel32,2);
            celstd=(nanstd(cel21'))';
            nucstd=(nanstd(cel31'))';
            objstd=(nanstd(cel32'))';
           
            p(1)=plot(myinterval,celavg,fcelallcolor{bonetype},'linewidth',0.5);% 'r'
            hold on 
            p(2)=plot(myinterval,nucavg,fnucallcolor{bonetype},'linewidth',0.2);% 'b'
            p(3)=plot(myinterval,objavg,'k-','linewidth',0.2);% 'b'

            legendarray2{1}=strcat('PC2/PC1');
            legendarray2{2}=strcat('PC3/PC1');
            legendarray2{3}=strcat('PC3/PC2');

            

%             d=celavg>nucavg;  limit=round(profilesize/2);   d1=max(find(d(1:limit))) ; d2=limit+min(find(d((limit+1):(profilesize-1))));
%             plot([myinterval(d1),myinterval(d1)],[0,3],'k:','linewidth',0.2)
%             hold on 
%             plot([myinterval(d2),myinterval(d2)],[0,3],'k:','linewidth',0.2)
            
            x=myinterval1;
            index = 1:profilesize-1; min_y1=celavg-celstd; max_y1=celavg+celstd;        
            x=myinterval;stacky2=(min_y1)';stacky1=(max_y1)';
            fillxx=x(index([1:end end:-1:1 1]));fillyy=[stacky1(index) stacky2(index(end:-1:1)) stacky1(1)];
            %h=fill(fillxx,fillyy,'r','EdgeColor','none','facealpha',0.3);
            
            
            min_y2=nucavg-nucstd; max_y2=nucavg+nucstd;  
            stacky2=(min_y2)';stacky1=(max_y2)';
            fillxx=x(index([1:end end:-1:1 1]));fillyy=[stacky1(index) stacky2(index(end:-1:1)) stacky1(1)];
            %h=fill(fillxx,fillyy,'b','EdgeColor','none','facealpha',0.2);

            min_y2=objavg-objstd; max_y2=objavg+objstd;  
            stacky2=(min_y2)';stacky1=(max_y2)';
            fillxx=x(index([1:end end:-1:1 1]));fillyy=[stacky1(index) stacky2(index(end:-1:1)) stacky1(1)];
            %h=fill(fillxx,fillyy,'k','EdgeColor','none','facealpha',0.2);


        


            index=1:15;
            [f,~]=fit(myinterval(index)',celavg(index),g,'startpoint',[5,10,1,500]);
            %plot(myinterval(index),f(myinterval(index)),strcat('r','-'),'linewidth',2)
           % ctext=strcat('\color{red}','\alpha =',sprintf('%0.1f',f.alpha),', \beta =',sprintf('%0.1f',f.beta), ', \alpha/\beta =',sprintf('%0.1f',f.alpha/f.beta) );

            index=1:25;
            [f,~]=fit(myinterval(index)',nucavg(index),g,'startpoint',[5,10,1,200]);
           % plot(myinterval(index),f(myinterval(index)),strcat('b','-'),'linewidth',2)
            %ntext=strcat('\color{blue}','\alpha =',sprintf('%0.1f',f.alpha),', \beta =',sprintf('%0.1f',f.beta),', \alpha/\beta =',sprintf('%0.1f',f.alpha/f.beta) );

%         index=8:39;
%         if bonetype==4
%              [f,~]=fit(myinterval(index)',celavg(index),logistic,'startpoint',[1,1,1,0.6]);
%         else 
%              [f,~]=fit(myinterval(index)',celavg(index),logistic,'startpoint',[1,1,1,0.6]);
%         end
%         plot(myinterval(index),f(myinterval(index)),'r.-');
%         ctext=strcat('C:','a =',sprintf('%0.1f',f.a) ,',b=',sprintf('%0.1e', f.b ),',c=',sprintf('%0.1e', f.c ),',d=',sprintf('%0.2f', f.d) );      
%         
%         [f,~]=fit(myinterval(index)',nucavg(index),logistic,'startpoint',[1,1,1,0.6]);
%         plot(myinterval(index),f(myinterval(index)),'r.--');
%         ntext=strcat('N:','a =',sprintf('%0.1f',f.a) ,',b=',sprintf('%0.1e', f.b ),',c=',sprintf('%0.1e', f.c ),',d=',sprintf('%0.2f', f.d) );
%        
      
        %text(0.02,2.2,ntext,'fontsize',9)
        %text(0.02,2.3,ctext,'fontsize',9)

        set(gca,'fontsize',11)
        legend(p,legendarray2,'location','southeast','fontsize',7);
        legend 'boxoff';
        %legend(q,legendarray,'location','northwest','fontsize',7);
        ylabel('Ratio of PC coefficients','fontsize',14); 
        xlabel('Bone long axis','fontsize',14)
        title(tname{bonetype},'fontweight','normal','fontsize',14)
        %axis([0,1,1,2.9])
        
        
        
%         ax2 = axes('Position',[XPos+0.025,YPos+0.2,0.44*Width,0.33*Height]);
%         box on;
%         plot(myinterval,celavg,strcat(fcelallcolor{bonetype},'-'),'linewidth',0.5); hold on 
%         plot(myinterval,nucavg,strcat(fnucallcolor{bonetype},'-'),'linewidth',0.2);% 'b'
%         axis([0,0.66,1.2,1.6])
%         set(gca,'fontsize',6)
%         
        
        
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


%saveas(h1,['PCs_Ratio']);
%saveas(h1,['PCs_Ratio.png']);



