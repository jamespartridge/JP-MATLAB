gamma_G=0.90;    ...prob. pays off if good
gamma_B=0.30;    ...prob. pays off if bad
y=1.5;           ...output
D=1;             ...investment size
r=1.01;          ...risk free rate
p1=0.5;          ...Pr(A|H)=Pr(C|L)=p1+(p2+p3)pi
p2=0.3;          ...Pr(B)=p2(1-pi)
p3=0.2;          ...Pr(C|H)=Pr(A|L)=p3(1-pi)

% z=1/13;
% alf=1;
% delta=1.2;

l=0.1;

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

wz=[0.6,0.8];

pi=0:0.01:1;
lp=length(pi);

G=zeros(3,lp);
B=zeros(3,lp);
for i=1:lp
   G(1,i)=par(6)+(par(7)+par(8))*pi(i);
   G(2,i)=par(7)*(1-pi(i));
   G(3,i)=par(8)*(1-pi(i));
   B(1,i)=par(8)*(1-pi(i));
   B(2,i)=par(7)*(1-pi(i));
   B(3,i)=par(6)+(par(7)+par(8))*pi(i);
end

for j=1:2
w=wz(j);
x=zeros(3,lp);
y=zeros(3,lp);
piRh=zeros(3,lp);
piRl=zeros(3,lp);
for i=1:lp
    %High signal interest rates
    x(1,i)=par(5)*(w*par(12)*G(1,i)+(1-w)*(1-par(12))*B(1,i))/(w*par(12)*G(1,i)*par(1)+(1-w)*(1-par(12))*B(1,i)*par(2));
    x(2,i)=par(5)*(w*par(12)*G(2,i)+(1-w)*(1-par(12))*B(2,i))/(w*par(12)*G(2,i)*par(1)+(1-w)*(1-par(12))*B(2,i)*par(2));
    x(3,i)=par(5)*(w*par(12)*G(3,i)+(1-w)*(1-par(12))*B(3,i))/(w*par(12)*G(3,i)*par(1)+(1-w)*(1-par(12))*B(3,i)*par(2));
    %Low signal interest rates
    y(1,i)=par(5)*((1-w)*par(12)*G(1,i)+w*(1-par(12))*B(1,i))/((1-w)*par(12)*G(1,i)*par(1)+w*(1-par(12))*B(1,i)*par(2));
    y(2,i)=par(5)*((1-w)*par(12)*G(2,i)+w*(1-par(12))*B(2,i))/((1-w)*par(12)*G(2,i)*par(1)+w*(1-par(12))*B(2,i)*par(2));
    y(3,i)=par(5)*((1-w)*par(12)*G(3,i)+w*(1-par(12))*B(3,i))/((1-w)*par(12)*G(3,i)*par(1)+w*(1-par(12))*B(3,i)*par(2));

piRh(1,i)=(w*par(12)*par(6)*(1-x(1,i)*par(1)/par(5))+(1-w)*(1-par(12))*par(8)*(1-x(1,i)*par(2)/par(5)))/...
            (w*par(12)*(par(7)+par(8))*(x(1,i)*par(1)/par(5)-1)+(1-w)*(1-par(12))*(-par(8))*(x(1,i)*par(2)/par(5)-1));
        
piRh(3,i)=(w*par(12)*par(8)*(1-x(3,i)*par(1)/par(5))+(1-w)*(1-par(12))*par(6)*(1-x(3,i)*par(2)/par(5)))/...
            (w*par(12)*(-par(8))*(x(3,i)*par(1)/par(5)-1)+(1-w)*(1-par(12))*(par(7)+par(8))*(x(3,i)*par(2)/par(5)-1));
            
piRl(1,i)=((1-w)*par(12)*par(6)*(1-y(1,i)*par(1)/par(5))+w*(1-par(12))*par(8)*(1-y(1,i)*par(2)/par(5)))/...
            ((1-w)*par(12)*(par(7)+par(8))*(y(1,i)*par(1)/par(5)-1)+w*(1-par(12))*(-par(8))*(y(1,i)*par(2)/par(5)-1));
        
piRl(3,i)=((1-w)*par(12)*par(8)*(1-y(3,i)*par(1)/par(5))+w*(1-par(12))*par(6)*(1-y(3,i)*par(2)/par(5)))/...
            ((1-w)*par(12)*(-par(8))*(y(3,i)*par(1)/par(5)-1)+w*(1-par(12))*(par(7)+par(8))*(y(3,i)*par(2)/par(5)-1));
end

PiRh(j,1,:)=piRh(1,:);
PiRh(j,3,:)=piRh(3,:);
PiRl(j,1,:)=piRl(1,:);
PiRl(j,3,:)=piRl(3,:);

X(j,1,:)=x(1,:);
X(j,3,:)=x(3,:);
Y(j,1,:)=y(1,:);
Y(j,3,:)=y(3,:);

end

piH1(1,:)=PiRh(1,1,:);
piH1(3,:)=PiRh(1,3,:);
piH2(1,:)=PiRh(2,1,:);
piH2(3,:)=PiRh(2,3,:);
piL1(1,:)=PiRl(1,1,:);
piL1(3,:)=PiRl(1,3,:);
piL2(1,:)=PiRl(2,1,:);
piL2(3,:)=PiRl(2,3,:);

X1(1,:)=X(1,1,:);
X1(3,:)=X(1,3,:);
X2(1,:)=X(2,1,:);
X2(3,:)=X(2,3,:);
Y1(1,:)=Y(1,1,:);
Y1(3,:)=Y(1,3,:);
Y2(1,:)=Y(2,1,:);
Y2(3,:)=Y(2,3,:);


figure(1)
plot(X1(1,:),piH1(1,:),'b',X2(1,:),piH2(1,:),'--b',X1(3,:),piH1(3,:),'r',X2(3,:),piH2(3,:),'--r')
title('High')
legend('A 0.6','A 0.8','C 0.6','C 0.8')


figure(2)
plot(Y1(1,:),piL1(1,:),'b',Y2(1,:),piL2(1,:),'--b',Y1(3,:),piL1(3,:),'r',Y2(3,:),piL2(3,:),'--r')
title('Low')
legend('A 0.6','A 0.8','C 0.6','C 0.8')
