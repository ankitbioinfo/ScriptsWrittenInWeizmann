clear all 
 
dirname='./Tibia/MakeListColumnarStructurePrediction/';
load([dirname,'centroid_and_surface_cells.mat']);


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

binsize=8;
zone=linspace(RZ,HZ,binsize);
Xbin=linspace(0,1,binsize-1);
Ybin=zeros(binsize-1,1); 

radius=load('Tibia/dataSave/ParameterForRandomModel_basedon_nuc_medium.dat');


for i=1:length(zone)-1
    cellsIdInXInterval=[];
%     for j=1:length(LCC)
%          cluster=mean(centroid(LCC{j},:));
%              if ((cluster(1) >zone(i)) &  (cluster(1) <=zone(i+1)))
%                  Ybin(i)=Ybin(i)+1;
%                  cellsIdInXInterval=[cellsIdInXInterval,j];
%              end
%     end
    
      for j=1:length(LCC)
           if ((radius(j,1) >Xbin(i)) &  (radius(j,1) <=Xbin(i+1)))
                 Ybin(i)=Ybin(i)+1;
                 cellsIdInXInterval=[cellsIdInXInterval,j];
           end
       end
           
      
      
      for k=2:size(real,2)
          a=real(cellsIdInXInterval,k)./random(cellsIdInXInterval,k);
          index=find(a~=Inf);
          AvgMeanReal{k-1}(i,1)= nanmean(a(index));
          AvgStdReal{k-1}(i,1)= nanstd( a(index) );
          
          
         
         stats = bootstrp(10000, @(x) [mean(x) std(x)], a(index));
          
          
          BootStrapMean{k-1}(i,1)= nanmean(stats(:,1));
          %BootStrapStd{k-1}(i,1)=  nanmean(stats(:,2));
          
          
          
%           variance = AvgStdReal{k-1}(i,1)^2;
%           [h,p]=vartest(a(index), 0);
%           ftest{k-1}(i,1)=p;
                    
      
      end

end


% AvgMeanReal{1}
% pause

dir2=strcat('MotifSpatialProfile','/');
if ~exist([dir2],'dir')
    mkdir([dir2]);
end


n=length( AvgMeanReal)

halfbinwidth=(Xbin(2)-Xbin(1))*0.5;


for k=1:n
    h=figure;
    set(gcf, 'PaperSize', [3 1.8]); %7
    set(gcf, 'PaperPosition', [0 0 3 1.8]);
    
           
                  
                   index=1:length(Xbin);
                    q(1)=plot(Xbin(index),AvgMeanReal{k}(index),'ro-','linewidth',1,'markersize',4,'markerfacecolor','r'); hold on 
   
                   plot([0,1],[1,1],'k-');
                   
                   for i=1:length(index)
                        zscore(i)=(AvgMeanReal{k}(i)-1)/AvgStdReal{k}(i);       
                        %zscore(i)=(AvgMeanReal{k}(i)-AvgMeanRandom{k}(i))/sqrt( (AvgStdReal{k}(i)^2)/Ybin(i) + (AvgStdRandom{k}(i)^2)/Ybin(i));
                         if ftest{k}(i)<0.05
                            plot(Xbin(i),AvgMeanReal{k}(i),'ks','linewidth',1,'markersize',15,'markerfacecolor','none'); 
                        end
                   end
%                    
                
%                    q(1)=errorbar(Xbin(index),AvgMeanReal{k}(index),zeros(size(Xbin(index))),AvgStdReal{k}(index),'ro-','linewidth',1,'markersize',10); hold on 
%                    q(2)=errorbar(Xbin(index),AvgMeanRandom{k}(index),zeros(size(Xbin(index))), AvgStdRandom{k}(index),'b-','linewidth',1,'markersize',10);
%    
                
                  
                   
                  
                    index = find(~isnan(AvgMeanReal{k})); 
                    if length(index)>0
                    min_y1=AvgMeanReal{k}(index)-AvgStdReal{k}(index); max_y1=AvgMeanReal{k}(index)+AvgStdReal{k}(index);  
                    x=Xbin(index);
                    index=1:length(x);
                    stacky2=(min_y1)';stacky1=(max_y1)';
                    fillxx=x(index([1:end end:-1:1 1]));fillyy=[stacky1(index) stacky2(index(end:-1:1)) stacky1(1)];
                    h=fill(fillxx,fillyy,'r','EdgeColor','none','facealpha',0.3);
                    
                    end
                    
                    hold on 
                   
                   p_two = 2*normcdf(-abs(zscore));                                
                   %ankit(k)=sum(p_two<0.05)
                   
%                    for i=1:length(p_two)
%                        if p_two(i)<=0.05
%                                 disp(' I am here')
%                                 plot([Xbin(i)-halfbinwidth Xbin(i)+halfbinwidth],[AvgMeanReal{k}(i)-AvgStdReal{k}(i)  AvgMeanReal{k}(i)+AvgStdReal{k}(i)],'b','linewidth',7);
%                          
%                        end
%                    end
%                    
                  
    

       %  legend(q,{'real','random'},'location','northwest')
       
    maxy=max(AvgMeanReal{k}+AvgStdReal{k}) ;
    miny=min(AvgMeanReal{k}-AvgStdReal{k}) ;
    
    set(gca,'fontsize',10)
    %xlabel('Long axis of bone','fontsize',10)
    %ylabel('<frequency>','fontsize',10);
     xlim([0,1])
     set(gca,'xtick',[0,0.5,1])
     
%      if miny~=maxy
%      set(gca,'ytick',[ 0, ceil(maxy/2)])
%      ylim([0,ceil(maxy)])
%      else
%          set(gca,'ytick',[0, 1]);
%          ylim([-1,1])
%      end
    % set(h, 'PaperPositionMode', 'auto');
    saveas(gcf,[dir2,'motif',num2str(k),'.png'])
    close all 
end



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
