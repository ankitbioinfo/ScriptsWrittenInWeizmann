clear all 


% pt_data_path={'Nuclei_and_Cells_PT_S17_m2_wt/','Nuclei_and_Cells_PT_S84_m3_wt/','Nuclei_and_Cells_PT_S18_m2_mut/',...
%            'Nuclei_and_Cells_PT_S17_m1_mut/', 'Nuclei_and_Cells_PT_S84_m1_mut/', 'Nuclei_and_Cells_PT_S84_m5_mut/' ,...
%            'Nuclei_and_Cells_PT_S18_m6_wt/', 'Nuclei_and_Cells_PT_S51_m2_wt/', 'Nuclei_and_Cells_PT_S84_m4_wt/'  };
% 
%        
% dt_data_path={'Nuclei_and_Cells_DT_S17_m2_wt/','Nuclei_and_Cells_DT_S84_m3_wt/','Nuclei_and_Cells_DT_S18_m2_mut/',...
%            'Nuclei_and_Cells_DT_S17_m1_mut/', 'Nuclei_and_Cells_DT_S84_m1_mut/', 'Nuclei_and_Cells_DT_S84_m5_mut/' ,...
%            'Nuclei_and_Cells_DT_S18_m6_wt/', 'Nuclei_and_Cells_DT_S51_m2_wt/', 'Nuclei_and_Cells_DT_S84_m4_wt/'  };



dt_data_path={  'data\Nuclei_and_Cells_PH_S93_m1_het_E18.5\', 'data\Nuclei_and_Cells_PH_S93_m2_mut_E18.5\'};
legendname={'m1 het', 'm2 mut','S18 m2 mt', 'S17 m1 mt','S84 m1 mut','S84 m5 mt','S18 m6 wt','S51 m2 wt','S84 m4 wt'};
mycolor={'r*','b*','g*','m*','c*','k*','y*'};



bigvec{1}=[1.1798e-17	-1.0053e-17	1
-0.76114	0.64859	1.55e-17
0.64859	0.76114	0];

bigvec{2}=[-2.299e-17	6.3439e-18	1
-0.96397	0.266	-2.3849e-17
-0.266	-0.96397	0];



h1=figure;
count=1;
for i=[1,2];%put all the bone here one by one to check alignment with sample 1 
    data=load([dt_data_path{i},'AlignedXYZ.mat']);
    %vec=data.vec;
    vec=bigvec{i};
    V=[data.alphaShapeAll];
    
    
    %The angle of rotation are chosen manually to visualize the bone that
    %they bulge are aligned together 
    if i==1
         %Rotation should done only for required sample 
         vec=vec*RotationMatrix_X(pi)*RotationMatrix_Y(pi);
    end

    if i==2 
         vec=vec*RotationMatrix_Z(pi)*RotationMatrix_Y(pi);
    end
    
    if i==3
         vec=vec*RotationMatrix_Z(180/180*pi)*RotationMatrix_X(0/180*pi);
    end
    
    if i==4
         vec=vec*RotationMatrix_Z(0/180*pi)*RotationMatrix_X(180/180*pi);
    end
    
    if i==5
          vec=vec*RotationMatrix_Z(180/180*pi)*RotationMatrix_X(pi);
    end
    
    if i==6
         vec=vec*RotationMatrix_Z(0/180*pi)*RotationMatrix_X(0/180*pi);
    end
    
    if i==7
         vec=vec*RotationMatrix_Z(180/180*pi)*RotationMatrix_X(0/180*pi);
    end
    
    if i==8
         vec=vec*RotationMatrix_Z(0/180*pi)*RotationMatrix_X(0/180*pi);
    end
    
    if i==9
         vec=vec*RotationMatrix_Z(0/180*pi)*RotationMatrix_X(0/180*pi);
    end
    
    
    V=V*vec;
    p(count)=plot3(V(:,1),V(:,2),V(:,3),mycolor{count},'markersize',2);
    hold on ;
    legendarray{count}=legendname{i};
    count=count+1;
    
  
    s=strsplit(dt_data_path{i},'Nuclei_and_Cells_');
    dlmwrite(strcat(s{2}(1:strlength(s{2})-1),'Alignment_matrix.dat'),vec,'\t');
end
    
legend(p,legendarray,'location','northeast');
   
    
xlabel('x')
ylabel('y')
zlabel('z')
axis image
view(42,11)



%saveas(h1,['PT_registered_bone_Together_']);
%saveas(h1,['PT_registered_bone_Together_.png']);


    
    
% factor=2.5;
% quiver3(mu(1),mu(2),mu(3),factor*vec(1,2),factor*vec(2,2),factor*vec(3,2),d(2),'k','LineWidth',3);
% quiver3(mu(1),mu(2),mu(3),-factor*vec(1,2),-factor*vec(2,2),-3*vec(3,2),d(2),'k','LineWidth',3);
% 
% quiver3(mu(1),mu(2),mu(3),factor*vec(1,1),factor*vec(2,1),3*vec(3,1),d(1),'b','LineWidth',3);
% quiver3(mu(1),mu(2),mu(3),-factor*vec(1,1),-factor*vec(2,1),-3*vec(3,1),d(1),'b','LineWidth',3);
% 
% quiver3(mu(1),mu(2),mu(3),factor*vec(1,3),factor*vec(2,3),factor*vec(3,3),d(3),'m','LineWidth',3);
% quiver3(mu(1),mu(2),mu(3),-factor*vec(1,3),-factor*vec(2,3),-factor*vec(3,3),d(3),'m','LineWidth',3);
% hold off;
%         

% 
% 
% [vec,val]=eig(cov(V2.V));
% V2=V2.V;%*vec;  V2=V2*RotationMatrix(pi);
% d = sqrt(diag(val));
% hold on;
% [vec,val]=eig(cov(V2));
% %plot(shp)
% mu=mean(V2);
% plot3(V2(:,1),V2(:,2),V2(:,3),'g*','markersize',2)
% % factor=0.5;
% quiver3(mu(1),mu(2),mu(3),factor*vec(1,2),factor*vec(2,2),factor*vec(3,2),d(2),'k','LineWidth',7);
% quiver3(mu(1),mu(2),mu(3),-factor*vec(1,2),-factor*vec(2,2),-3*vec(3,2),d(2),'k','LineWidth',7);
% 
% quiver3(mu(1),mu(2),mu(3),factor*vec(1,1),factor*vec(2,1),3*vec(3,1),d(1),'b','LineWidth',7);
% quiver3(mu(1),mu(2),mu(3),-factor*vec(1,1),-factor*vec(2,1),-3*vec(3,1),d(1),'b','LineWidth',7);
% 
% quiver3(mu(1),mu(2),mu(3),factor*vec(1,3),factor*vec(2,3),factor*vec(3,3),d(3),'m','LineWidth',7);
% quiver3(mu(1),mu(2),mu(3),-factor*vec(1,3),-factor*vec(2,3),-factor*vec(3,3),d(3),'m','LineWidth',7);
% hold off;
%         
% 
% 
% [vec,val]=eig(cov(V3.V));
% V3=V3.V;%*vec;
% d = sqrt(diag(val));
% hold on;
% [vec,val]=eig(cov(V3));
% %plot(shp)
% mu=mean(V3);
% plot3(V3(:,1),V3(:,2),V3(:,3),'b*','markersize',2)
% 





function R=RotationMatrix_Z(theta)
         R=[cos(theta) -sin(theta) 0;
            sin(theta) cos(theta) 0;
            0 0 1];
end



function R=RotationMatrix_X(theta)
         R=[1 0 0;
            0 cos(theta) -sin(theta);
            0 sin(theta) cos(theta)];
end

function R=RotationMatrix_Y(theta)
         R=[cos(theta) 0 sin(theta);
             0  1  0;
            -sin(theta) 0 cos(theta)];
end


