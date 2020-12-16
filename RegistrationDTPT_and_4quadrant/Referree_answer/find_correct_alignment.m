

clear all 
close all 


a1=load('Nuclei_and_Cells_S126_m3_PT_DV/AlignedXYZ.mat');
a2=load('Nuclei_and_Cells_S126_m3_PT_DV/centerAxisLine.dat');
a3=load('Nuclei_and_Cells_S126_m3_PT_DV/centroid.dat');


b1=load('Nuclei_and_Cells_S17_m2_pt_thresh05/AlignedXYZ.mat');
b2=load('Nuclei_and_Cells_S17_m2_pt_thresh05/centerAxisLine.dat');
b3=load('Nuclei_and_Cells_S17_m2_pt_thresh05/centroid.dat');



[vec1,val1]= eig(cov([a3;a2])); 
[vec2,val2]= eig(cov([b3;b2]));



a3=a3*vec1;
b3=b3*vec2;

figure 
plot3(a3(:,1),a3(:,2),a3(:,3),'b.')
hold on 
plot3(b3(:,1),b3(:,2),b3(:,3),'r.')


figure 

v=a1.alphaShapePT;
vector1=a1.vec*vec1;
v=v*vector1;

plot3(v(:,1),v(:,2),v(:,3),'b.')
hold on 


v=b1.alphaShapePT;
vector2=b1.vec*vec2*RotationMatrix_Z(pi)*RotationMatrix_X(pi);
v=v*vector2;


plot3(v(:,1),v(:,2),v(:,3),'r.')




% d=[1 1 1];
% d*RotationMatrix_Y(pi)*RotationMatrix_X(pi)
% d*RotationMatrix_X(pi)
% d*RotationMatrix_Y(pi)





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

