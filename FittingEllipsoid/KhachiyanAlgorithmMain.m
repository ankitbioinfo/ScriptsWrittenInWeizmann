function output=KhachiyanAlgorithmMain(P)

N=size(P,1);
d=size(P,2);
Q=  [P'; ones(1,N)] ;
QT = Q';

tolerance=0.01;
err=1+tolerance;

U=zeros(N);
for i=1:N
    U(i,i)=1/N;
end

while (err> tolerance)
     V=Q*U*QT;
     Mtemp=  QT*inv(V)*Q;
     M=diag(Mtemp);
     [maximum, argmax] = max(M);
     step_size = (maximum - d - 1.0) / ((d + 1.0) * (maximum - 1.0));
     newU = (1.0 - step_size) * U;
     newU(argmax,argmax)=newU(argmax,argmax) + step_size;
     err = norm(diag(U)-diag(newU));
     U=newU;
end

center=P'*diag(U);
output.center=center;

A= inv(P'*U*P  -  center*center')/d ; 

[U,s,rotation]=svd(A);

radii= 1./sqrt(diag(s));
output.radii=radii;

newaxes = diag(radii);  
%rotate accordingly
newaxes = newaxes*rotation';
% plot axes

for i=1:length(newaxes) 
	p=newaxes(i,:);
	X3 = linspace(-p(1), p(1), 100) + center(1);
	Y3 = linspace(-p(2), p(2), 100) + center(2);
        Z3 = linspace(-p(3), p(3), 100) + center(3);
	[sa,sb]=min(radii);
	[sa1,sb1]=max(radii);
	if sb==i
		%plot3(X3,Y3,Z3,'k')
		Minor_axis=[X3',Y3',Z3'];
	elseif sb1==i
		%plot3(X3,Y3,Z3,'b')
		Major_axis=[X3',Y3',Z3'];
	else
		%plot3(X3,Y3,Z3,'y')
		Medium_axis=[X3',Y3',Z3'];
	end
	%hold on 
end

output.major=Major_axis;
output.minor=Minor_axis;
output.medium=Medium_axis;

u = linspace(0.0, 2.0 * pi, 100);
v = linspace(0.0, pi, 100);
        
% cartesian coordinates that correspond to the spherical angles:
x = radii(1) * cos(u)'*sin(v);       %np.outer(np.cos(u), np.sin(v))
y = radii(2) * sin(u)'*sin(v);			    %np.outer(np.sin(u), np.sin(v))
z = radii(3) * ones(1,length(u))'*cos(v);                     %np.outer(np.ones_like(u), np.cos(v))



% rotate accordingly
for i=1:length(x)
     for j=1:length(x)
          r = [x(i,j),y(i,j),z(i,j)]*rotation' + center';  
	  x(i,j)=r(1);
	  y(i,j)=r(2);
	  z(i,j)=r(3);  
	  %np.dot([x[i,j],y[i,j],z[i,j]], rotation) + center
     end
end

output.x=x;
output.y=y;
output.z=z;

%[x, y, z] = ellipsoid(0,0,0,radii(1),radii(2),radii(3),40);
%tt=[x(:) y(:) z(:)]*inv(U);
%nx=reshape(tt(:,1),size(x,1),size(x,2));
%ny=reshape(tt(:,2),size(x,1),size(x,2));
%nz=reshape(tt(:,3),size(x,1),size(x,2));
%hSurface=surf(nx+center(1), ny+center(2), nz+center(3), 'FaceColor','r','EdgeColor','none','FaceAlpha',0.6);


%hSurface=surf(x,y,z,'FaceColor','r','EdgeColor','none','FaceAlpha',0.6);
%hold on 


%plot3(P(:,1),P(:,2),P(:,3),'k.-','markersize',3,'linewidth',0.25)

%print(gcf,'-djpeg','-r400',strcat('view3d',num2str(500)));
%hold off
%close all

