% 

a=load('data/Nuclei_and_Cells_DT_S17_m1_mut/c_n_pos1 (Characteristics).mat');

C=a.C.surfaces(1).vertices; 
C=C-mean(C);
plot3(C(:,1),C(:,2),C(:,3),'r.')

[vec,val]=eig(cov(C));
ovec=vec;
oval=val;
mu=mean(C);
d = sqrt(diag(val));
hold on;
factor=3; quiver3(mu(1),mu(2),mu(3),factor*vec(1,1),factor*vec(2,1),factor*vec(3,1),d(1),'r','LineWidth',10);
factor=3; quiver3(mu(1),mu(2),mu(3),factor*vec(1,2),factor*vec(2,2),factor*vec(3,2),d(2),'k','LineWidth',5);
factor=3; quiver3(mu(1),mu(2),mu(3),factor*vec(1,3),factor*vec(2,3),factor*vec(3,3),d(3),'m','LineWidth',2);
 
factor=-3; quiver3(mu(1),mu(2),mu(3),factor*vec(1,1),factor*vec(2,1),factor*vec(3,1),d(1),'r','LineWidth',10);
factor=-3; quiver3(mu(1),mu(2),mu(3),factor*vec(1,2),factor*vec(2,2),factor*vec(3,2),d(2),'k','LineWidth',5);
factor=-3; quiver3(mu(1),mu(2),mu(3),factor*vec(1,3),factor*vec(2,3),factor*vec(3,3),d(3),'m','LineWidth',2);

hold off;

xlabel('x');ylabel('y');zlabel('z');



r= [ 0.7123   -0.3094   -0.6300
    0.4597    0.8840    0.0855
    0.5304   -0.3505    0.7719];


figure
C=C*r;

plot3(C(:,1),C(:,2),C(:,3),'r.')

[vec,val]=eig(cov(C));
svec=vec;
sval=val;
mu=mean(C);
d = sqrt(diag(val));
hold on;
quiver3(mu(1),mu(2),mu(3),vec(1,2),vec(2,2),vec(3,2),d(2),'k','LineWidth',5);
quiver3(mu(1),mu(2),mu(3),vec(1,1),vec(2,1),vec(3,1),d(1),'r','LineWidth',5);
quiver3(mu(1),mu(2),mu(3),vec(1,3),vec(2,3),vec(3,3),d(3),'m','LineWidth',5);
hold off;

xlabel('x');ylabel('y');zlabel('z');

figure 
plot3(C(:,1),C(:,2),C(:,3),'r.')
vec=ovec*r;
hold on;
quiver3(mu(1),mu(2),mu(3),vec(1,2),vec(2,2),vec(3,2),d(2),'k','LineWidth',5);
quiver3(mu(1),mu(2),mu(3),vec(1,1),vec(2,1),vec(3,1),d(1),'r','LineWidth',5);
quiver3(mu(1),mu(2),mu(3),vec(1,3),vec(2,3),vec(3,3),d(3),'m','LineWidth',5);
hold off;

xlabel('x');ylabel('y');zlabel('z');


%subplot(1,2,2)

% 
% clear ;
% s  = [2 2] ;
% set = randn(200,1);
% x = normrnd(s(1).*set,1)+3;
% x = zscore(x); % Standardize
% y = normrnd(s(1).*set,1)+2;
% y= zscore(y);%Standardize
% x_0 = mean(x);
% y_0 = mean (y) ;
% c = linspace(1,100,length(x)); % color
% 
% scatter(x,y,100,c,'filled')
% xlabel('1st Feature : x')
% ylabel('2nd Feature : y')
% title('2D_dataset')
% 
% grid on
% % gettign the covariance matrix 
% covariance = cov([x,y]);
% % getting the eigenvalues and the  eigenwert 
% [eigen_vector, eigen_values] = eig(covariance); 
% eigen_vector_1 = eigen_vector(:,1);
% eigen_vector_2 = eigen_vector(:,2);
% d = sqrt(diag(eigen_values));
% 
% % ploting the eigenvectors ! 
% hold on 
% % x_0 = repmat(x_0,size(eigen_vector_2,1),1);
% % y_0 = repmat(y_0,size(eigen_vector_1,1),1);
% 
% quiver(x_0,y_0,eigen_vector(1,2),eigen_vector(2,2),d(2),'k','LineWidth',5);
% quiver(x_0,y_0,eigen_vector(1,1),eigen_vector(2,1),d(1),'r','LineWidth',5);
% hold off;
% 











