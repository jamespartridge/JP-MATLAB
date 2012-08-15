clear all

pi=0:0.001:1;
lp=length(pi);

gamma_G=0.90;    ...prob. pays off if good
gamma_B=0.50;    ...prob. pays off if bad
y=1.5;           ...output
D=1;             ...investment size
r=1.01;          ...risk free rate
p1=0.5;          ...Pr(A|H)=Pr(C|L)=p1+(p2+p3)pi
p2=0.3;          ...Pr(B)=p2(1-pi)
p3=0.2;          ...Pr(C|H)=Pr(A|L)=p3(1-pi)

z=1/1000;
alf=1.5;
delta=2;

l=0.6;
w=0.7;

%parameter vector
par=[gamma_G;   %1
     gamma_B;   %2
     y;         %3
     D;         %4
     r;         %5
     p1;        %6
     p2;        %7
     p3;        %8
     z;         %9
     alf;       %10
     delta;     %11
     l];        %12

for i=1:lp
    %High signal interest rates
    x(1,i)=par(5)*(w*l*(par(6)+(par(7)+par(8))*pi(i))+(1-w)*(1-l)*par(8)*(1-pi(i)))/(w*l*(par(6)+(par(7)+par(8))*pi(i))*par(1)+(1-w)*(1-l)*par(8)*(1-pi(i))*par(2));
    x(2,i)=par(5)*(w*l+(1-w)*(1-l))/(w*l*par(1)+(1-w)*(1-l)*par(2));
    x(3,i)=par(5)*(w*l*par(8)*(1-pi(i))+(1-w)*(1-l)*(par(6)+(par(7)+par(8))*pi(i)))/(w*l*par(8)*(1-pi(i))*par(1)+(1-w)*(1-l)*(par(6)+(par(7)+par(8))*pi(i))*par(2));
    %Low signal interest rates
    y(1,i)=par(5)*((1-w)*l*(par(6)+(par(7)+par(8))*pi(i))+w*(1-l)*par(8)*(1-pi(i)))/((1-w)*l*(par(6)+(par(7)+par(8))*pi(i))*par(1)+w*(1-l)*par(8)*(1-pi(i))*par(2));
    y(2,i)=par(5)*((1-w)*l+w*(1-l))/((1-w)*l*par(1)+w*(1-l)*par(2));
    y(3,i)=par(5)*((1-w)*l*par(8)*(1-pi(i))+w*(1-l)*(par(6)+(par(7)+par(8))*pi(i)))/((1-w)*l*par(8)*(1-pi(i))*par(1)+w*(1-l)*(par(6)+(par(7)+par(8))*pi(i))*par(2));
end


plot(pi,x(1,:),pi,x(2,:),pi,x(3,:),pi,y(1,:),pi,y(2,:),pi,y(3,:))
legend('Rah','Rbh','Rch','Ral','Rbl','Rcl')