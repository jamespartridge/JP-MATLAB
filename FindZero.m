tic
clear all

gG=1.00;    ...prob. pays off if good
gB=0.30;    ...prob. pays off if bad
y=3;           ...output
D=1;             ...investment size
r=1.01;          ...risk free rate
g1=0.5;          ...Pr(A|G)=g1+(g2+g3)pi
g2=0.2;          ...Pr(B|G)=g2(1-pi)
g3=0.3;          ...Pr(C|G)=g3(1-pi)
b1=0.5;          ...Pr(A|L)=b1+(b2+b3)pi
b2=0.2;          ...Pr(B|L)=b2(1-pi)
b3=0.3;          ...Pr(C|L)=b3(1-pi)
alf=5;           ...c(pi)=1/alpha * pi^alpha
l=0.6;

par=[gG;
    gB;
    y;
    D;
    r;
    g1;
    g2;
    g3;
    b1;
    b2;
    b3;
    alf;
    l];

%omega grid
w=0.5:0.001:1;
% w=0.5:0.05:1;
lw=length(w);

%initialize matrix sizes
piH=zeros(1,lw);
piL=zeros(1,lw);
Rh=zeros(3,lw);
Rl=zeros(3,lw);

%calculate EQ values at each omega
for k=1:lw
    [piH(k),piL(k),Rh(1,k),Rh(2,k),Rh(3,k),...
        Rl(1,k),Rl(2,k),Rl(3,k)]=FP(w(k),par);
end

RETh(1,:)=y-D*Rh(1,:);
RETh(2,:)=y-D*Rh(2,:);
RETh(3,:)=y-D*Rh(3,:);
RETl(1,:)=y-D*Rl(1,:);
RETl(2,:)=y-D*Rl(2,:);
RETl(3,:)=y-D*Rl(3,:);

MBHbyw=zeros(1,lw);
MBLbyw=zeros(1,lw);
for i=1:lw
aMBHw1(i)=(RETh(1,i)>0)*...
         (y-D*Rh(1,i))*l*(1-l)*(gG*(g2+g3)-gB*(-b3))/((w(i)*l+(1-w(i))*(1-l))^2);
bMBHw1(i)=(RETh(2,i)>0)*...
         (y-D*Rh(2,i))*l*(1-l)*(gG*(-g2)-gB*(-b2))/((w(i)*l+(1-w(i))*(1-l))^2);
cMBHw1(i)=(RETh(3,i)>0)*...
         (y-D*Rh(3,i))*l*(1-l)*(gG*(-g3)-gB*(b2+b3))/((w(i)*l+(1-w(i))*(1-l))^2);
aMBHw2(i)=(RETh(1,i)>0)*...
          ((w(i)*l*gG*(g2+g3)+(1-w(i))*(1-l)*gB*(-b3))/(w(i)*l+(1-w(i))*(1-l)))*(-D*r*l*(1-l)*(g1+(g2+g3)*piH(i))*b3*(1-piH(i))*(gB-gG)/((w(i)*l*(g1+(g2+g3)*piH(i))*gG+(1-w(i))*(1-l)*b3*(1-piH(i))*gB)^2));
bMBHw2(i)=(RETh(2,i)>0)*...
          ((w(i)*l*gG*(-g2)+(1-w(i))*(1-l)*gB*(-b2))/(w(i)*l+(1-w(i))*(1-l)))*(-D*r*l*(1-l)*(gB-gG)/((w(i)*l*gG+(1-w(i))*(1-l)*gB)^2));
cMBHw2(i)=(RETh(3,i)>0)*...
          ((w(i)*l*gG*(-g3)+(1-w(i))*(1-l)*gB*(b2+b3))/(w(i)*l+(1-w(i))*(1-l)))*(-D*r*l*(1-l)*g3*(1-piH(i))*(b1+(b2+b3)*piH(i))*(gB-gG)/((w(i)*l*b3*(1-piH(i))*gG+(1-w(i))*(1-l)*(b1+(b2+b3)*piH(i))*gB)^2));

aMBHbyw(i)=aMBHw1(i)+aMBHw2(i);
bMBHbyw(i)=bMBHw1(i)+bMBHw2(i);
cMBHbyw(i)=cMBHw1(i)+cMBHw2(i);

MBHbyw(i)=aMBHbyw(i)+bMBHbyw(i)+cMBHbyw(i);

%rating terms of marginal benefit (ys cancel)
aMBH(i)=(RETh(1,i)>0)*...
         (y-D*Rh(1,i))*(w(i)*l*gG*(g2+g3)+(1-w(i))*(1-l)*gB*(-b3))/(w(i)*l+(1-w(i))*(1-l));
bMBH(i)=(RETh(2,i)>0)*...
         (y-D*Rh(2,i))*(w(i)*l*gG*(-g2)+(1-w(i))*(1-l)*gB*(-b2))/(w(i)*l+(1-w(i))*(1-l));
cMBH(i)=(RETh(3,i)>0)*...
         (y-D*Rh(3,i))*(w(i)*l*gG*(-g3)+(1-w(i))*(1-l)*gB*(b2+b3))/(w(i)*l+(1-w(i))*(1-l));


end



% MBH=aMBH+bMBH+cMBH;
% piHbyw=(1/(alf-1))*(MBH.^((2-alf)/(alf-1))).*(MBHbyw);
% figure(1)
% plot(w,aMBH,w,bMBH,w,cMBH,w,MBH.^(-2/3))
% legend('a','b','c','total')
% title('MBH')
figure(2)
plot(w,piH,w,piL)
title('EQ pi')
% figure(3)
%     w,aMBHw1,w,bMBHw1,w,cMBHw1,...
%     w,aMBHw2,w,bMBHw2,w,cMBHw2,...
% plot(...
%     w,aMBHw1+bMBHw1+cMBHw1,...
%     w,aMBHw2+bMBHw2+cMBHw2,...
%     w,MBHbyw)
% line([0.5 1],[0 0],'color','k')
% title('d MBH / d w')
% legend('total1','total2','total','Location','best')

%EQ signal probabilities
probH=w*l+(1-w)*(1-l); 
probL=(1-w)*l+w*(1-l);

%EQ rating probabilities conditional on signal, type
probAgivGH=par(6)+(par(7)+par(8))*piH;
probAgivGL=par(6)+(par(7)+par(8))*piL;
probBgivGH=par(7)*(1-piH);
probBgivGL=par(7)*(1-piL);
probCgivGH=par(8)*(1-piH);
probCgivGL=par(8)*(1-piL);
probAgivBH=par(8)*(1-piH);
probAgivBL=par(8)*(1-piL);
probBgivBH=par(7)*(1-piH);
probBgivBL=par(7)*(1-piL);
probCgivBH=par(6)+(par(7)+par(8))*piH;
probCgivBL=par(6)+(par(7)+par(8))*piL;

%EQ rating probabilities conditional on signal
probAgivH=probAgivGH*l+probAgivBH*(1-l);
probAgivL=probAgivGL*l+probAgivBL*(1-l);
probBgivH=probBgivGH*l+probBgivBH*(1-l);
probBgivL=probBgivGL*l+probBgivBL*(1-l);
probCgivH=probCgivGH*l+probCgivBH*(1-l);
probCgivL=probCgivGL*l+probCgivBL*(1-l);

%EQ rating probabilities
numA=probAgivH.*probH+probAgivL.*probL;
numB=probBgivH.*probH+probBgivL.*probL;
numC=probCgivH.*probH+probCgivL.*probL;

figure(5)
plot(w,numA,w,numB,w,numC)
legend('A','B','C')

figure(6)
plot(w,Rl(1,:)-Rh(1,:),'--g',w,Rl(2,:)-Rh(2,:),'-b')
title('Interest rate spreads, within rating')
xlabel('\omega')
legend('A','B','Location','best')


hc=min(RETh(3,:))
lc=min(RETl(3,:))