clear variables, close all
% three random points in 3d
P = rand(3,3);
% vectors
V1 = P(2,:)-P(1,:);
V2 = P(3,:)-P(1,:);
% add one point on the plane defined by the previous three points
P = [P; P(1,:)+V1+V2];
% plot points
figure, hold on, axis equal
plot3(P(:,1),P(:,2),P(:,3),'o');
quiver3(P(1,1),P(1,2),P(1,3),V1(1),V1(2),V1(3),'AutoScale','off')
quiver3(P(1,1),P(1,2),P(1,3),V2(1),V2(2),V2(3),'AutoScale','off')
% plane center
B0 = mean(P,1);
% plot
plot3(B0(1),B0(2),B0(3),'*');
% orthogonal vector
W = cross(V2,V1);
% unit vector
U = W/sqrt(sum(W.^2,2));
% project B by the distance k
k = 2;
B = B0+k*U;
% plot
quiver3(B0(1),B0(2),B0(3),k*U(1),k*U(2),k*U(3),'AutoScale','off');
plot3(B(1),B(2),B(3),'s');