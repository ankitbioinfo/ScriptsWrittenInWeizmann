clear all 
 
dirname='./Tibia/MakeListColumnarStructurePrediction/';
load([dirname,'centroid_and_surface_cells.mat'],'centroid');


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

binsize=16; %16
zone=linspace(RZ,HZ,binsize);
Xbin=linspace(0,1,binsize-1);
Ybin=zeros(binsize-1,1);



[LCC,LCC1]=readClusterFile([dirname,'ClusterMediumCutoff.dat']);


real=load('RealSubgraphIsomorphism.dat');
[real,realWpb]=Freq2Pb(real);


for i=1:100
    temp=load(['Random_node_name/Random_',num2str(i),'/Frequency.dat']);
    [randomExp,randomWpb]=Freq2Pb(temp);
    random{i}=randomExp;
end


dir2=strcat('MotifSpatialProfile','/');
if ~exist([dir2],'dir')
    mkdir([dir2]);
end
 

radius=load('Tibia/dataSave/ParameterForRandomModel_basedon_nuc_medium.dat');

 bootstrapN=200000;
 
if exist([dir2,'SaveBootstrapPvalueMotifs.mat'], 'file') == 0   
    
    disp('I am here')

for i=1:length(zone)-1
    cellsIdInXInterval=[];
      for j=1:length(LCC)
           if ((radius(j,1) >Xbin(i)) &  (radius(j,1) <=Xbin(i+1)))
                 Ybin(i)=Ybin(i)+1;
                 cellsIdInXInterval=[cellsIdInXInterval,j];
           end
      end
           
      
      for k=1:size(real,2) 
          clear data 
          for randset=1:length(random)
              f1=real(cellsIdInXInterval,k);
              f2=random{randset}(cellsIdInXInterval,k);
%               a=f1./f2;
%               index=find(a~=Inf);              
              data2(randset,1)=nansum(f2);
          end
          
          value=nansum(f1);
          ObservedDifferences= abs(value - nanmean(data2));
          
          pooldata=[data2; value];
          randomizedDifferences=zeros(bootstrapN,1);
          for repeat=1:bootstrapN
              sampling2=randi(length(pooldata),100,1);
              sampling1=randi(length(pooldata),1,1);
              randomizedDifferences(repeat,1)= abs(nanmean(pooldata(sampling1)) - nanmean(pooldata(sampling2)));
          end

          good=sum( randomizedDifferences >=  ObservedDifferences);
          empiricalPvalue= good/ bootstrapN;
          ftest{k}(i,1)=empiricalPvalue;
          
          
       
%           stats = bootstrp(10000, @(x) [mean(x) std(x)], data);
%           BootstrapMean{k-1}(i,1)=mean(stats(:,1));
%           BootStrapStd{k-1}(i,1)= mean(stats(:,2));
%           if isnan(nanmean(data))
%               data
%           end
          
          AvgMeanReal{k}(i,1)=value;
          AvgStdReal{k}(i,1)=nanstd(f1);
          AvgMeanRandom{k}(i,1)= nanmean( data2);
          AvgStdRandom{k}(i,1)= nanstd( data2);
         
          
%           goodIndex=find(~isnan(data));%           
%           shifteddata= data(goodIndex) -mean(data(goodIndex));%           
%           [h,p]=ttest(shifteddata);
%           ftest{k}(i,1)=p; 
      end      
%           variance = AvgStdReal{k-1}(i,1)^2;
%           [h,p]=vartest(a(index), 0);
%           ftest{k-1}(i,1)=p;  

end

save([dir2,'SaveBootstrapPvalueMotifs.mat'], 'ftest','AvgMeanReal','AvgStdReal','AvgMeanRandom','AvgStdRandom');

else
    disp('I am not here');
    load([dir2,'SaveBootstrapPvalueMotifs.mat']);
    end







n=length( AvgMeanReal)

halfbinwidth=(Xbin(2)-Xbin(1))*0.5;


for k=1:n
    h=figure;
    set(gcf, 'PaperSize', [3 1.8]); %7
    set(gcf, 'PaperPosition', [0 0 3 1.8]);
    
           
                  
                   index=1:length(Xbin);
                    q(1)=plot(Xbin(index),AvgMeanReal{k}(index),'ro-','linewidth',1,'markersize',4,'markerfacecolor','r'); hold on 
                    q(2)=plot(Xbin(index),AvgMeanRandom{k}(index),'bo-','linewidth',0.5,'markersize',2,'markerfacecolor','b');
                   %plot([0,1],[1,1],'k-');
                   
                   for i=1:length(index)
                        %zscore(i)=(AvgMeanReal{k}(i)-1)/AvgStdReal{k}(i);       
                        %zscore(i)=(AvgMeanReal{k}(i)-AvgMeanRandom{k}(i))/sqrt( (AvgStdReal{k}(i)^2)/Ybin(i) + (AvgStdRandom{k}(i)^2)/Ybin(i));
                         if (ftest{k}(i)<0.05)&(ftest{k}(i)>0)
                            plot(Xbin(i),AvgMeanReal{k}(i),'ks','linewidth',1,'markersize',15,'markerfacecolor','none'); 
                        end
                   end
                   
                   
                   
                   
                   
                   
%                    
                
%                    q(1)=errorbar(Xbin(index),AvgMeanReal{k}(index),zeros(size(Xbin(index))),AvgStdReal{k}(index),'ro-','linewidth',1,'markersize',10); hold on 
%                    q(2)=errorbar(Xbin(index),AvgMeanRandom{k}(index),zeros(size(Xbin(index))), AvgStdRandom{k}(index),'b-','linewidth',1,'markersize',10);
%    
                
                  
                   
                  
%                     index = find(~isnan(AvgMeanReal{k})); 
%                     if length(index)>0
%                     min_y1=AvgMeanReal{k}(index)-AvgStdReal{k}(index); max_y1=AvgMeanReal{k}(index)+AvgStdReal{k}(index);  
%                     x=Xbin(index);
%                     index=1:length(x);
%                     stacky2=(min_y1)';stacky1=(max_y1)';
%                     fillxx=x(index([1:end end:-1:1 1]));fillyy=[stacky1(index) stacky2(index(end:-1:1)) stacky1(1)];
%                     h=fill(fillxx,fillyy,'r','EdgeColor','none','facealpha',0.3);
%                     end

                     index1 = find(~isnan(AvgMeanReal{k})); index2= find(~isnan(AvgMeanRandom{k}));
                    index = intersect(index1,index2);

                    min_y1=AvgMeanRandom{k}(index)-AvgStdRandom{k}(index); max_y1=AvgMeanRandom{k}(index)+AvgStdRandom{k}(index);  
                    x=Xbin(index);
                    index=1:length(x);
                    stacky2=(min_y1)';stacky1=(max_y1)';
                    fillxx=x(index([1:end end:-1:1 1]));fillyy=[stacky1(index) stacky2(index(end:-1:1)) stacky1(1)];
                    h=fill(fillxx,fillyy,'b','EdgeColor','none','facealpha',0.2);



                    
                    hold on 
                   
                   %p_two = 2*normcdf(-abs(zscore));                                
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
       
      %maxy=max(AvgMeanRandom{k}+AvgStdRandom{k}) ;  
      maxy=max([max(max_y1),max(AvgMeanReal{k})]);
      miny=min([min(min_y1), min(AvgMeanReal{k})]) ;
    
      
    ankit(k,:)=[miny,maxy];  
    
    ax = gca;
    ax.YAxis.Exponent = -3;
    
    yticks([0 .010 .020 .030])
    ytickformat('%01.0f')
    %ylim([miny,maxy])
    ylim([-0.0045,0.0398])
    
    
    set(gca,'fontsize',10)
    %xlabel('Long axis of bone','fontsize',10)
    %ylabel('<frequency>','fontsize',10);
     xlim([0,1])
     set(gca,'xtick',[0,0.5,1])
     
%      if miny~=maxy
%      set(gca,'ytick',[ 0, ceil(maxy/2)])
     
%      else
%          set(gca,'ytick',[0, 1]);
%          ylim([-1,1])
%      end
    % set(h, 'PaperPositionMode', 'auto');
    saveas(gcf,[dir2,'motif',num2str(k),'.png'])
    close all 
end




function [b,weightedPb]=Freq2Pb(a)
    % a is input frequencies of 65 types of motifs 
    total=sum(a);
    for i=1:size(a,1)
        b(i,:)=a(i,:)./total;
    end 


    totalsum=sum(total);
    for i=1:size(a,2)
        total(i)=total(i)/totalsum;
    end

    weightedPb=total; 
    
    for i=1:size(a,2)
         b(:,i)=weightedPb(i)*b(:,i);
    end 

    
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
