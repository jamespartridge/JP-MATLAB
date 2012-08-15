clear all
clc

gamma_G=0.90;    ...prob. pays off if good
gamma_B=0.50;    ...prob. pays off if bad
y=1.5;           ...output
D=1;             ...investment size
r=1.01;          ...risk free rate
p1=0.5;          ...Pr(A|H)=Pr(C|L)=p1+(p2+p3)pi
p2=0.3;          ...Pr(B)=p2(1-pi)
p3=0.2;          ...Pr(C|H)=Pr(A|L)=p3(1-pi)

%parameter vector
par=[gamma_G;   %1
     gamma_B;   %2
     y;         %3
     D;         %4
     r;         %5
     p1;        %6
     p2;        %7
     p3];       %8

w=0.7;
l=0.7;          ...fraction of firms that are "good"
z=1e-3;

pi=0:0.1:1;
lp=length(pi);

for i=1:lp
    
    x(1,i)=par(5)*(w*l*(par(6)+(par(7)+par(8))*pi(i))+(1-w)*(1-l)*par(8)*(1-pi(i)))/(w*l*(par(6)+(par(7)+par(8))*pi(i))*par(1)+(1-w)*(1-l)*par(8)*(1-pi(i))*par(2));
    x(2,i)=par(5)*(w*l+(1-w)*(1-l))/(w*l*par(1)+(1-w)*(1-l)*par(2));
    x(3,i)=par(5)*(w*l*par(8)*(1-pi(i))+(1-w)*(1-l)*(par(6)+(par(7)+par(8))*pi(i)))/(w*l*par(8)*(1-pi(i))*par(1)+(1-w)*(1-l)*(par(6)+(par(7)+par(8))*pi(i))*par(2));
    y(1,i)=par(5)*((1-w)*l*(par(6)+(par(7)+par(8))*pi(i))+w*(1-l)*par(8)*(1-pi(i)))/((1-w)*l*(par(6)+(par(7)+par(8))*pi(i))*par(1)+w*(1-l)*par(8)*(1-pi(i))*par(2));
    y(2,i)=par(5)*((1-w)*l+w*(1-l))/((1-w)*l*par(1)+w*(1-l)*par(2));
    y(3,i)=par(5)*((1-w)*l*par(8)*(1-pi(i))+w*(1-l)*(par(6)+(par(7)+par(8))*pi(i)))/((1-w)*l*par(8)*(1-pi(i))*par(1)+w*(1-l)*(par(6)+(par(7)+par(8))*pi(i))*par(2));
    for j=1:3
        RETh(j,i)=par(3)-par(4)*x(j,i);
        RETl(j,i)=par(3)-par(4)*y(j,i);
    end
    B(1,i)=((par(5)*(w*l+(1-w)*(1-l)))^(-1))*(...
         (w*l*par(1)*(par(6)+(par(7)+par(8))*pi(i))+(1-w)*(1-l)*par(2)*par(8)*(1-pi(i)))*RETh(1,i)*(RETh(1,i)>0)...
        +(w*l*par(1)+(1-w)*(1-l)*par(2))*par(7)*(1-pi(i))*RETh(2,i)*(RETh(2,i)>0)...
        +(w*l*par(1)*par(8)*(1-pi(i))+(1-w)*(1-l)*par(2)*(par(6)+(par(7)+par(8))*pi(i)))*RETh(3,i)*(RETh(3,i)>0));
    
    B(2,i)=((par(5)*((1-w)*l+w*(1-l)))^(-1))*(...
         ((1-w)*l*par(1)*(par(6)+(par(7)+par(8))*pi(i))+w*(1-l)*par(2)*par(8)*(1-pi(i)))*RETl(1,i)*(RETl(1,i)>0)...
        +((1-w)*l*par(1)+w*(1-l)*par(2))*par(7)*(1-pi(i))*RETl(2,i)*(RETl(2,i)>0)...
        +((1-w)*l*par(1)*par(8)*(1-pi(i))+w*(1-l)*par(2)*(par(6)+(par(7)+par(8))*pi(i)))*RETl(3,i)*(RETl(3,i)>0));
    
    M(i)=z*(pi(i)/(1-pi(i)));
    
end

for i=1:lp
    [pi(i),RETh(1,i),RETh(2,i),RETh(3,i),B(1,i)]
end
for i=1:lp
    [pi(i),RETl(1,i),RETl(2,i),RETl(3,i),B(2,i)]
end

figure(1)
subplot(2,2,1)
plot(pi,RETh(1,:),pi,RETh(2,:),'--',pi,RETh(3,:),':')
title('EQ returns for H')
xlabel('\pi')
legend('A','B','C','Location','SouthWest')

subplot(2,2,2)
plot(pi,B(1,:))
xlabel('\pi')
title('Exp. return for H')

subplot(2,2,3)
plot(pi,RETl(1,:),pi,RETl(2,:),pi,RETl(3,:),':')
title('EQ returns for L')
xlabel('\pi')
legend('A','B','C','Location','SouthWest')

subplot(2,2,4)
plot(pi,B(2,:))
xlabel('\pi')
title('Exp. return for L')

