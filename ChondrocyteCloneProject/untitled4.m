
clear all 

a=load('RealSubgraphIsomorphism.dat');

total=sum(a);

for i=1:size(a,1)
    b(i,:)=a(i,:)./total;
end 
    

noc=[

totalsum=sum(total);

for i=1:size(a,2)
    total(i)=total(i)/totalsum;
end

weightedPb=total; 

for i=1:size(a,2)
    b(:,i)=weightedPb(i)*b(:,i);
end 




dlmwrite('RealSubgraphPb.dat',b,'\t')
