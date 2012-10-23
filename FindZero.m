tic
clear all

gG=1.00;    ...prob. pays off if good
gB=0.30;    ...prob. pays off if bad
y=3;           ...output
D=1;             ...investment size
r=1.01;          ...risk free rate
g1=0.4;          ...Pr(A|G)=g1+(g2+g3)pi
g2=0.55;          ...Pr(B|G)=g2(1-pi)
g3=0.05;          ...Pr(C|G)=g3(1-pi)
b2=g2;          ...Pr(B|L)=b2(1-pi)
b3=g3;          ...Pr(C|L)=b1+b3(1-pi)
b1=1-b2-b3;          ...Pr(A|L)=(b2+b3)pi
alf=3;           ...c(pi)=1/alpha * pi^alpha
l=0.34;

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
w=0.5:0.0001:1;
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

%sometimes interest rate is 0/0 at omega=1
Rl(isnan(Rl(:,lw)))=r/gB;
Rh(isnan(Rh(:,lw)))=r/gG;

RETh(1,:)=y-D*Rh(1,:);
RETh(2,:)=y-D*Rh(2,:);
RETh(3,:)=y-D*Rh(3,:);
RETl(1,:)=y-D*Rl(1,:);
RETl(2,:)=y-D*Rl(2,:);
RETl(3,:)=y-D*Rl(3,:);

MBHbyw=zeros(1,lw);
MBLbyw=zeros(1,lw);
for i=1:lw
% %first deriv of MBH by omega
% %high by rating, first term
% aMBHw1(i)=(RETh(1,i)>0)*...
%          (y-D*Rh(1,i))*l*(1-l)*(gG*(g2+g3)-gB*(b2+b3))/((w(i)*l+(1-w(i))*(1-l))^2);
% bMBHw1(i)=(RETh(2,i)>0)*...
%          (y-D*Rh(2,i))*l*(1-l)*(gG*(-g2)-gB*(-b2))/((w(i)*l+(1-w(i))*(1-l))^2);
% cMBHw1(i)=(RETh(3,i)>0)*...
%          (y-D*Rh(3,i))*l*(1-l)*(gG*(-g3)-gB*(-b3))/((w(i)*l+(1-w(i))*(1-l))^2);
% %high by rating, second term
% aMBHw2(i)=(RETh(1,i)>0)*...
%           ((w(i)*l*gG*(g2+g3)+(1-w(i))*(1-l)*gB*(b2+b3))/(w(i)*l+(1-w(i))*(1-l)))*(-D*r*l*(1-l)*(g1+(g2+g3)*piH(i))*(b1+b3*(1-piH(i)))*(gB-gG)/((w(i)*l*(g1+(g2+g3)*piH(i))*gG+(1-w(i))*(1-l)*(b2+b3)*piH(i)*gB)^2));
% bMBHw2(i)=(RETh(2,i)>0)*...
%           ((w(i)*l*gG*(-g2)+(1-w(i))*(1-l)*gB*(-b2))/(w(i)*l+(1-w(i))*(1-l)))*(-D*r*l*(1-l)*(gB-gG)/((w(i)*l*gG+(1-w(i))*(1-l)*gB)^2));
% cMBHw2(i)=(RETh(3,i)>0)*...
%           ((w(i)*l*gG*(-g3)+(1-w(i))*(1-l)*gB*(-b3))/(w(i)*l+(1-w(i))*(1-l)))*(-D*r*l*(1-l)*g3*(1-piH(i))*(b1+b3*(1-piH(i)))*(gB-gG)/((w(i)*l*g3*(1-piH(i))*gG+(1-w(i))*(1-l)*(b1+b3*(1-piH(i)))*gB)^2));
% %total for each rating
% aMBHw(i)=aMBHw1(i)+aMBHw2(i);
% bMBHw(i)=bMBHw1(i)+bMBHw2(i);
% cMBHw(i)=cMBHw1(i)+cMBHw2(i);
% MBHw(i)=aMBHw(i)+bMBHw(i)+cMBHw(i);
% %rating terms of marginal benefit (ys cancel)
% aMBH(i)=(RETh(1,i)>0)*...
%          (y-D*Rh(1,i))*(w(i)*l*gG*(g2+g3)+(1-w(i))*(1-l)*gB*(b2+b3))/(w(i)*l+(1-w(i))*(1-l));
% bMBH(i)=(RETh(2,i)>0)*...
%          (y-D*Rh(2,i))*(w(i)*l*gG*(-g2)+(1-w(i))*(1-l)*gB*(-b2))/(w(i)*l+(1-w(i))*(1-l));
% cMBH(i)=(RETh(3,i)>0)*...
%          (y-D*Rh(3,i))*(w(i)*l*gG*(-g3)+(1-w(i))*(1-l)*gB*(-b3))/(w(i)*l+(1-w(i))*(1-l));
% MBH=aMBH+bMBH+cMBH;
     
%first deriv of MBL by omega
%low by rating, first term
aMBLw1(i)=(RETl(1,i)>0)*...
         (y-D*Rl(1,i))*l*(1-l)*(gG*(g2+g3)-gB*(b2+b3))/(((1-w(i))*l+w(i)*(1-l))^2);
bMBLw1(i)=(RETl(2,i)>0)*...
         (y-D*Rl(2,i))*l*(1-l)*(gG*(-g2)-gB*(-b2))/(((1-w(i))*l+w(i)*(1-l))^2);
cMBLw1(i)=(RETl(3,i)>0)*...
         (y-D*Rl(3,i))*l*(1-l)*(gG*(-g3)-gB*(-b3))/(((1-w(i))*l+w(i)*(1-l))^2);
%low by rating, second term
aMBLw2(i)=(RETl(1,i)>0)*...
          ((1-w(i))*l*gG*(g2+g3)+w(i)*(1-l)*gB*(b2+b3))/((1-w(i))*l+w(i)*(1-l))*(-D*r*l*(1-l)*(g1+(g2+g3)*piL(i))*(b1+b3*(1-piL(i)))*(gB-gG)/(((1-w(i))*l*(g1+(g2+g3)*piL(i))*gG+w(i)*(1-l)*(b2+b3)*piL(i)*gB)^2));
bMBLw2(i)=(RETl(2,i)>0)*...
          ((1-w(i))*l*gG*(-g2)+w(i)*(1-l)*gB*(-b2))/((1-w(i))*l+w(i)*(1-l))*(-D*r*l*(1-l)*(gB-gG)/(((1-w(i))*l*gG+w(i)*(1-l)*gB)^2));
cMBLw2(i)=(RETl(3,i)>0)*...
          ((1-w(i))*l*gG*(-g3)+w(i)*(1-l)*gB*(-b3))/((1-w(i))*l+w(i)*(1-l))*(-D*r*l*(1-l)*g3*(1-piL(i))*(b1+b3*(1-piL(i)))*(gB-gG)/(((1-w(i))*l*g3*(1-piL(i))*gG+w(i)*(1-l)*(b1+b3*(1-piL(i)))*gB)^2));
%total for each rating
aMBLw(i)=aMBLw1(i)+aMBLw2(i);
bMBLw(i)=bMBLw1(i)+bMBLw2(i);
cMBLw(i)=cMBLw1(i)+cMBLw2(i);
MBLw(i)=aMBLw(i)+bMBLw(i)+cMBLw(i);
%rating terms of marginal benefit (ys cancel)
aMBL(i)=(RETl(1,i)>0)*...
         (y-D*Rl(1,i))*((1-w(i))*l*gG*(g2+g3)+w(i)*(1-l)*gB*(b2+b3))/((1-w(i))*l+w(i)*(1-l));
bMBL(i)=(RETl(2,i)>0)*...
         (y-D*Rl(2,i))*((1-w(i))*l*gG*(-g2)+w(i)*(1-l)*gB*(-b2))/((1-w(i))*l+w(i)*(1-l));
cMBL(i)=(RETl(3,i)>0)*...
         (y-D*Rl(3,i))*((1-w(i))*l*gG*(-g3)+w(i)*(1-l)*gB*(-b3))/((1-w(i))*l+w(i)*(1-l));
MBL=aMBL+bMBL+cMBL;

end

%Eq ratings investment
figure(1)
plot(w,piH,w,piL)
legend('H','L')
title('EQ pi')

%Different parts of the der. of piH wrt omega
% piHw=(1/(alf-1))*(MBH.^((2-alf)/(alf-1))).*MBHw;
% figure(2)
% plot(w,aMBH,w,bMBH,w,cMBH,w,MBH)
% legend('a','b','c','total')
% title('MBH')
%Different parts of the der. of piL wrt omega
piLw=(1/(alf-1))*(MBL.^((2-alf)/(alf-1))).*MBLw;
figure(2)
plot(w,aMBL,w,bMBL,w,cMBL,w,MBL)
line([0.5 1],[0 0],'color','k')
title('MBL')
legend('a','b','c','total')

%Plot different terms of partial der. of MB wrt omega
% figure(3)
% plot(...
%     w,aMBHw1+bMBHw1+cMBHw1,...
%     w,aMBHw2+bMBHw2+cMBHw2,...
%     w,MBHw)
% line([0.5 1],[0 0],'color','k')
% title('d MBH / d w')
figure(3)
plot(...
    w,aMBLw1,w(1:68),aMBLw2(1:68),...
    w,bMBLw1,w,bMBLw2,...
    w,cMBLw1,w,cMBLw2,...
    w(1:68),MBLw(1:68))
line([0.5 1],[0 0],'color','k')
title('d MBL / d w')
legend('a1','a2','b1','b2','c1','c2','total','Location','best')

%%%%%%%%%%
%Ratings Distribution
%%%%%%%%%%
%EQ signal probabilities
probH=w*l+(1-w)*(1-l); 
probL=(1-w)*l+w*(1-l);
%EQ rating probabilities conditional on signal, type
probAgivGH=g1+(g2+g3)*piH;
probAgivGL=g1+(g2+g3)*piL;
probBgivGH=g2*(1-piH);
probBgivGL=g2*(1-piL);
probCgivGH=g3*(1-piH);
probCgivGL=g3*(1-piL);
probAgivBH=(b2+b3)*piH;
probAgivBL=(b2+b3)*piL;
probBgivBH=b2*(1-piH);
probBgivBL=b2*(1-piL);
probCgivBH=b1+b3*(1-piH);
probCgivBL=b1+b3*(1-piL);
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

prAgvG=probAgivGH.*w.*l+probAgivGL.*(1-w).*(1-l);
prAgvB=probAgivBH.*(1-w)+probAgivBL.*w;
prGgvA=prAgvG./(prAgvG.*l+prAgvB);
prBgvG=probBgivGH.*w.*l+probBgivGL.*(1-w).*(1-l);
prBgvB=probBgivBH.*(1-w)+probBgivBL.*w;
prGgvB=prBgvG./(prBgvG.*l+prBgvB);
prCgvG=probCgivGH.*w.*l+probCgivGL.*(1-w).*(1-l);
prCgvB=probCgivBH.*(1-w)+probCgivBL.*w;
prGgvC=prCgvG./(prCgvG.*l+prCgvB);

figure(4)
plot(w,prGgvA,w,prGgvB,w,prGgvC)

%plot ratings distribution
figure(5)
plot(w,numA,w,numB,w,numC)
legend('A','B','C')

%%%%%%%%%
%Within Rating Dispersion
%%%%%%%%%
% figure(6)
% plot(w,Rl(1,:)-Rh(1,:),'--g',w,Rl(2,:)-Rh(2,:),'-b')
% title('Interest rate spreads, within rating')
% xlabel('\omega')
% legend('A','B','Location','best')

% for i=1:lw
% test(i)=g1*(b2+b3)*(w(i)*l*(g2+g3)*gG+(1-w(i))*(1-l)*(b2+b3)*gB)/((w(i)*l*(g1+(g2+g3)*piH(i))*gG+(1-w(i))*(1-l)*(b2+b3)*piH(i)*gB)^2)...
%         -g3*b1*(w(i)*l*g3*gG+(1-w(i))*(1-l)*b3*gB)/((w(i)*l*g3*(1-piH(i))*gG+(1-w(i))*(1-l)*b3*(1-piH(i))*gB)^2);
% end
% 


hc=min(RETh(3,:))
lc=min(RETl(3,:))