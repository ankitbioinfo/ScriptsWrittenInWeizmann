

n=100;
count=1;
for i=1:n
    for j=i+1:n
        for k=j+1:n
            m=sort([i,j,k],'descend');
            if m(1)/m(2)>1.5 
            radii(count,:)=
            count=count+1;
        end
    end
end

pc21=radii(:,2)./radii(:,1);
pc31=radii(:,3)./radii(:,1);
pc32=radii(:,3)./radii(:,2);

result=[radii, pc21,pc31,pc32];

[yr,xr]=histnorm(result(:,4),50);plot(xr,yr,'r.-');
hold on 
[yr,xr]=histnorm(result(:,5),50);
plot(xr,yr,'b.-')
[yr,xr]=histnorm(result(:,6),50);
plot(xr,yr,'k.-')
axis([0,1,0,5])