clear all 
close all 
 
dirname='./Tibia/MakeListColumnarStructurePrediction/';
load([dirname,'centroid_and_surface_cells.mat']);
n=479;

real=load('RealSubgraphIsomorphism.dat');
random=load('Random_Mean_SubgraphIsomorphism.dat');


[LCC,LCC1]=readClusterFile([dirname,'ClusterMediumCutoff.dat']);


c=centroid-mean(centroid);
bonetype=1;
if (bonetype==3)|(bonetype==1)
       c=[-c(:,1),c(:,2:3)];
else
       c=[c(:,1:2),c(:,3)];
end

HZ=min(c(:,1))
RZ=max(c(:,1))
disp('Resting And HZ')
RZ=-562;
HZ=677;
[RZ, HZ]


%pause

binsize=16;
zone=linspace(RZ,HZ,binsize);
Xbin=linspace(0,1,binsize-1);
Ybin=zeros(binsize-1,1); 

radius=load('Tibia/dataSave/ParameterForRandomModel_basedon_nuc_medium.dat');

dir2=strcat('MotifSpatialProfile','/');
load([dir2,'SaveBootstrapPvalueMotifs.mat']);



n=length( AvgMeanReal);

halfbinwidth=(Xbin(2)-Xbin(1))*0.5;


for k=1:n
                   for i=1:length(Xbin)
                         if ftest{k}(i)<0.05
                            Highest(i,k)=AvgMeanReal{k}(i);
                            HighestVar(i,k)=AvgStdReal{k}(i);
                            
                            randomHighest(i,k)=AvgMeanRandom{k}(i);
                            randomHighestVar(i,k)=AvgStdRandom{k}(i);
                        end
                   end

   
end



[motifType,highFreqMotif,stdFreqMotif,random_highFreqMotif,random_stdFreqMotif,x]=getHighest(Highest,HighestVar,randomHighest,randomHighestVar,Xbin);
motifType

h=figure;
    set(gcf, 'PaperSize', [5 3]); %7
    set(gcf, 'PaperPosition', [0 0 5 3]);

    p(1)=plot(x,highFreqMotif,'ro-','linewidth',1,'markersize',4,'markerfacecolor','r'); hold on
    p(2)=plot(x,random_highFreqMotif,'bo-','linewidth',0.5,'markersize',2,'markerfacecolor','b'); hold on
    hold on 
    
    
                
%                     min_y1=highFreqMotif-stdFreqMotif; max_y1=highFreqMotif+stdFreqMotif;  
%                     index=1:length(x);
%                     stacky2=(min_y1)';stacky1=(max_y1)';
%                     fillxx=x(index([1:end end:-1:1 1]));fillyy=[stacky1(index); stacky2(index(end:-1:1)); stacky1(1)];
%                     h=fill(fillxx,fillyy,'r','EdgeColor','none','facealpha',0.3);
                    
                    
                    min_y1=random_highFreqMotif-random_stdFreqMotif; max_y1=random_highFreqMotif+random_stdFreqMotif;  
                    index=1:length(x);
                    stacky2=(min_y1)';stacky1=(max_y1)';
                    fillxx=x(index([1:end end:-1:1 1]));fillyy=[stacky1(index); stacky2(index(end:-1:1)); stacky1(1)];
                    h=fill(fillxx,fillyy,'b','EdgeColor','none','facealpha',0.2);
                    
                    
                
    ax = gca;
    ax.YAxis.Exponent = -3;                
                    
    legend(p,{'cre','background'},'location','northeast')
    
    xlim([0,1])                
    xlabel('Long axis of bone')
    ylabel('Largest expectation value');                
                    
    title('Largest expected motifs in real network','fontweight','normal')
    %saveas(h,'HighestMotifType.png')
    print(gcf,'Largest_expected_real.png','-dpng','-r300');
    






function [LCC,LCC1]=readClusterFile(name) 
        fid = fopen(name,'rt');
        tline = fgetl(fid);
        count=1;
         while ischar(tline)
            line= split(tline,' ');
             if length(line)>3
                for j=1:length(line)-1
                    LCC{count}(j)=str2num(line{j});
                end
             end    
            
             if length(line)>2
                for j=1:length(line)-1
                    LCC1{count}(j)=str2num(line{j});
                end
             end    
             
            
            tline = fgetl(fid);
            count=count+1;
        end
        fclose(fid);
end


function [motifType,real_highFreqMotif,real_stdFreqMotif,random_highFreqMotif,random_stdFreqMotif,x]=getHighest(real_Highest,real_HighestVar,random_Highest,random_HighestVar,Xbin)
    for i=1:size(real_Highest,1)
        [sa,sb]=max(real_Highest(i,:));
        real_highFreqMotif(i)=sa;
        random_highFreqMotif(i)=random_Highest(i,sb);
        x(i,1)=Xbin(i);
        real_stdFreqMotif(i)=real_HighestVar(i,sb);
        random_stdFreqMotif(i)=random_HighestVar(i,sb);
        motifType(i)=sb;
    end
end