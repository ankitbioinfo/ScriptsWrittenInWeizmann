# FittingEllipsoid
Fit ellipsoid to 3d data points using Khachiyan Algorithm 

This code will find a minimum volume of ellipsoid which contains set of all the 3d data points. 


 output=KhachiyanAlgorithmMain(data);
 plot3(data(:,1),data(:,2),data(:,3),'b.')
 hold on 
 hsurface=surf(output.x,output.y,output.z,'FaceColor','b','EdgeColor','none','FaceAlpha',0.1);
 radii=output.radii;
 center=output.center;
 
% plot principal axis  
plot3(output.minor(:,1),output.minor(:,2),output.minor(:,3),'k','linewidth',0.2)
plot3(output.major(:,1),output.major(:,2),output.major(:,3),'k','linewidth',5)
plot3(output.medium(:,1),output.medium(:,2),output.medium(:,3),'k','linewidth',2)
 
 
major_axis_vector=output.major(1,:)-output.major(end,:);
medium_axis_vector=output.medium(1,:)-output.medium(end,:);
minor_axis_vector=output.minor(1,:)-output.minor(end,:);


if you have any doubt related to using the code please ask me here ankitbioinfo[at]gmail.com
