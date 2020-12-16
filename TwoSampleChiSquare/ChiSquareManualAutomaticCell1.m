

clear all 
path = './data/Nuclei_and_Cells_DT_S17_m2_wt/';

a=load([path,'all_cells.mat']);
cel=a.all_cells;
%DT_wt_S17 = RZ 16-25; PZ 11-15; PHZ  8-10; HZ 1-7

for i=1:length(cel)
	for RZ=16:25
		if cel(i,1)==RZ			
		   indRZ(i)=1;
		end
	end 	
	for PZ=11:15
		if cel(i,1)==PZ			
		   indPZ(i)=1;
		end
	end 	
	
	for HZ=1:10
		if cel(i,1)==HZ			
		   indHZ(i)=1;
		end
	end 	
end

indRZ=find(indRZ==1);
indPZ=find(indPZ==1);
indHZ=find(indHZ==1);	

[length(indRZ), length(indPZ), length(indHZ)]
size(cel)

% the columns contains the individual nuclei features
% 1 - stack id
% 2 - volume
% 3 - surface area
% 4 - sphericity
% 5-7 - centroid x,y,z coordinates
% 8-10 - PC1 x,y,z orientation
% 11-13 - PC2 x,y,z orientation
% 14-16 - PC3 x,y,z orientation
% 17-19 - PC1,PC2,PC3 latent coefficient
% 20 - Delaunay density

manual='./data/validating_segmentation_manualconfocal_vs_automaticlightsheet/cells/';

T= readtable([manual,'pz_cells_S17_m4_seg.xls'], 'FileType', 'spreadsheet', 'Sheet','Nuclei');
manPZ=[T.surface_area T.volume];	
T= readtable([manual,'rz_cells_S17_m4_seg.xls'], 'FileType', 'spreadsheet', 'Sheet','Nuclei');
manRZ=[T.surface_area T.volume];
T= readtable([manual,'hz_cells_S17_m4_seg.xls'], 'FileType', 'spreadsheet', 'Sheet','Nuclei');
mHZ=[T.surface_area T.volume];
T= readtable([manual,'phz_cells_S17_m4_seg.xls'], 'FileType', 'spreadsheet', 'Sheet','Nuclei');
mPHZ=[T.surface_area T.volume];
manHZ=[mPHZ;mHZ];
data={manRZ, manPZ, manHZ};
% 
% disp('volume')
% x={indRZ, indPZ, indPHZ, indHZ}; 
% for i=1:4 
%      index=x{i};  allmanual=data{i};    
%      zonevol(i,:)=[min(cel(index,2)) max(cel(index,2)) min(allmanual(:,2)) max(allmanual(:,2))];
% end
% 
% disp('surface area')
% for i=1:4
%      index=x{i};  allmanual=data{i}; 
%      zsa(i,:)=[min(cel(index,3)) max(cel(index,3)) min(allmanual(:,1)) max(allmanual(:,1))];
% end





h2=figure();

set(gcf, 'PaperSize', [10 4]);
set(gcf, 'PaperPosition', [0 0 10 4]);

for figplot=1:2
subplot(1,2,figplot) 

if figplot==1
    index=2; ind=2;
    gap=200:100:4000; binsize=length(gap);
    interval=linspace(200,4000,binsize+1);bingap=interval(3)-interval(2);
else
    index=1; ind=3;
    gap=200:100:2000; binsize=length(gap);
    interval=linspace(200,2000,binsize+1);bingap=interval(3)-interval(2);
end

poolsize=20;

clear cellauto 
clear cellman

for iter=1:10000
t=randi(length(indRZ),poolsize,1); lRZ=indRZ(t);
t=randi(length(indPZ),poolsize,1); lPZ=indPZ(t);
t=randi(length(indHZ),poolsize,1); lHZ=indHZ(t);

y1=manRZ(randi(length(manRZ),poolsize,1),index); 
y2=manPZ(randi(length(manPZ),poolsize,1),index);
y4=manHZ(randi(length(manHZ),poolsize,1),index);

x1=(cel(lRZ,ind)); 
x2=(cel(lPZ,ind)); 
x4=(cel(lHZ,ind)); 

manual=[x1;x2;x4];  auto=[y1; y2; y4];
%auto=[manRZ(:,index); manPZ(:,index);manPHZ(:,index); manHZ(:,index)];

cellauto(iter,:)=repeat1000times(interval,binsize,auto); 
cellman(iter,:)=repeat1000times(interval,binsize,manual); 

end
%[min(ankit(:,1)), min(ankur(:,1)), max(ankit(:,2)), max(ankur(:,2))]


freq_auto=mean(cellauto);
freq_man=mean(cellman);

            
if figplot==1
    ind=1:binsize;  
    barWidth=bingap*0.006;
else
    ind=1:binsize;
    barWidth=bingap*0.006;
end
p1=bar(interval(ind),freq_auto(ind), 'BarWidth', barWidth); hold on 
p2=bar(interval(ind)-0.1*bingap,freq_man(ind), 'BarWidth', barWidth);

ind=find((freq_man>0)&(freq_auto>0));

K1= sqrt( sum(freq_man(ind))/sum(freq_auto(ind))   ); dof=length(ind)-1;
K2= sqrt( sum(freq_auto(ind))/sum(freq_man(ind))   ); 

chi2=sum(((K1*freq_auto(ind) - K2*freq_man(ind)).^2)./(freq_auto(ind)+freq_man(ind)));

pvalue=1-chi2cdf(chi2,dof);

set(p2,'FaceColor','red','facealpha',0.2);
set(p1,'FaceColor','blue','facealpha',0.2);

set(gca,'fontsize',7)
legend([p1,p2],'Automatic','Manual','location','northeast')

   if figplot==1
       xlabel('Volume (\mum^3)');
       ylabel('<Cell Number>')
   else
       xlabel('Surface area (\mum^2)');
       ylabel('<Cell Number>')
   end
    set(gca,'XTickLabelRotation',45)
   title(strcat('\chi^2=',sprintf('%0.2f',chi2),  ', p-value : ', sprintf('%0.2f',pvalue)));

end

saveas(h2,['Cell_segmentation_accuracy']);
saveas(h2,['Cell_segmentation_accuracy.png']);


function freq_man=repeat1000times(interval,binsize,manual)

freq_man=zeros(binsize,1);
for i=1:length(manual)
    for j=2:length(interval)
        if manual(i)< interval(j)
            freq_man(j-1)=freq_man(j-1)+1;
            break 
        end 
    end 
end 
end
