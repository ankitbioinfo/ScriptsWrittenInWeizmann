
% map=colormap(parula(10));
% 
% map=[0,0,0;    1,0.4,0.4; 0.0582    0.4677    0.8589;   0.5 0.5 0; ...
%      1 1 0.8784; 1,1,0; 0,1,0;...
%       1 0.6 0.73;...
%       0,1,1; 0.29 0 0.51;0.25 0.25 0];  %101
%   
  
% map= flipud(redbluecmap(10));
 map=flipud(colormap(jet(7)));
 %map= flipud(redgreencmap(6,'Interpolation','linear' ));

%theta=[0,5,10,15,20,25,30,45,60,75,90];
theta=[0,5,10,15,25,40,60,90];


h=figure();
set(gcf, 'PaperSize', [2 5]);
set(gcf, 'PaperPosition', [0 0 2 5]);


%tname={'<5','5-10','10-15','15-20','20-25','25-30','30-45','45-60','60-75','75-90'};
tname={'<5','5-10','10-15','15-25','25-40','40-60','60-90'};

for i=2:length(theta)
plot(0.1,-i*0.7,'marker','o','markersize',20, 'markerfacecolor',map(i-1,:),'color',map(i-1,:))
    hold on;
    text(0.3,-i*0.7,tname{i-1},'fontsize',18);
end

axis([0,1,-6.2,-0.8])
title('Orientation colormap')

set(findobj(gcf,'type','axes'),'visible','off');

saveas(h,['OrientationColormap']);
saveas(h,['OrientationColormap','.png']);