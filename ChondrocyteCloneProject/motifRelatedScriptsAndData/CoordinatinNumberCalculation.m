function CoordinatinNumberCalculation( mydata, Xbin,flag,P1,P2)

%    data1=d1.dataAvgMeanFeatures{i,j};   % individual point 
%    data2=(d2.realization_AvgMeanFeatures{i}(j,:))';
%    
%     
        
    g=fittype(@(m,c,x) (m*x+c));
    RZ=[];
    PZ=[];
    PHZ=[];
    HZ=[];
    for j=1:length(Xbin)
        
        if flag==1
            real_size= mydata{P1,j};
            real_rg=mydata{P2,j};
          
        end    
        
        if flag==2
            real_size=mydata(j,1);
            real_rg=mydata(j,2);
        end
            
        if Xbin(j)<=0.2
            RZ=[RZ; [real_size,real_rg]];
        elseif (Xbin(j)>0.2)&(Xbin(j)<=0.5)
            PZ=[PZ; [real_size,real_rg]];
        elseif (Xbin(j)>0.5)&(Xbin(j)<=0.7)
            PHZ=[PHZ; [real_size,real_rg]];
        elseif (Xbin(j)>0.7)&(Xbin(j)<=1)
            HZ=[HZ; [real_size,real_rg]];
        end 

    end
    
    zone={RZ,PZ,PHZ,HZ};
    zonecolor={'r.','b.','g.','k.'};
    legname={'RZ','PZ','PHZ','HZ'};
    
    for i=1:length(zone)
        data=zone{i};
      
        plot(data(:,1),data(:,2),zonecolor{i});
        hold on 

        index=find(~isnan(data(:,2)));
        x=data(index,1);
        [f,~]=fit(x,data(index,2),g,'startpoint',[1,0]);
        p(i)=plot(x,f(x),[zonecolor{i}(1),'-']);

        Df=f.m; 
        lname{i}=strcat(legname{i},' :', 'slope =',sprintf('%0.2f',Df));
    
    end

        legend(p, lname,'location','northeast')
        %xlabel('log(cluster size)');
        

   
end