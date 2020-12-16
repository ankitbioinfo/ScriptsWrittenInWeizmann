

syms 'qt';

a=3;
b=2;

p=[1,1];
p2=[-1,1];
l=[p2-p];

eqn= (p(1)+t*l(1))^2/a^2  +  (p(2)+t*l(2))^2/b^2 ==1;

solt = vpasolve(eqn,t);

for i=1:length(solt)
    var=double(solt(i));
    P(i,:)=[p(1)+var*l(1),p(2)+var*l(2)];
end

P