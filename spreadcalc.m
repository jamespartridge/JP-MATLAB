function [piH,piL,sprdA,sprdB,sprdC,avgMBh,avgMBl,MCh,MCl]=spreadcalc(w,par)

l=0.8;        ...fraction of firms that are "good"
% w=0.7;      ...prob. signal is accurate
z=.2;
alf=1;

% Pgh=w*l/(w*l+(1-w)*(1-l));
% Pbh=(1-w)*(1-l)/(w*l+(1-w)*(1-l));
% Pgl=w*(1-l)/(w*(1-l)+(1-w)*l);
% Pbl=(1-w)*l/(w*(1-l)+(1-w)*l);

pi=0:0.001:1;
lp=length(pi);

x=zeros(3,lp);
y=zeros(3,lp);

for i=1:lp

x(1,i)=par(5)*(w*l*(par(6)+(par(7)+par(8))*pi(i))+(1-w)*(1-l)*par(8)*(1-pi(i)))/(w*l*(par(6)+(par(7)+par(8))*pi(i))*par(1)+(1-w)*(1-l)*par(8)*(1-pi(i))*par(2));
x(2,i)=par(5)*(w*l+(1-w)*(1-l))/(w*l*par(1)+(1-w)*(1-l)*par(2));
x(3,i)=par(5)*(w*l*par(8)*(1-pi(i))+(1-w)*(1-l)*(par(6)+(par(7)+par(8))*pi(i)))/(w*l*par(8)*(1-pi(i))*par(1)+(1-w)*(1-l)*(par(6)+(par(7)+par(8))*pi(i))*par(2));
% x(4,i)=par(5)*w*l*(1-w)*(1-l)*par(8)*((par(6)+(par(7)+par(8))*pi(i))+(par(7)+par(8))*(1-pi(i)))*(par(2)-par(1))/((w*l*(par(6)+(par(7)+par(8))*pi(i))*par(1)+(1-w)*(1-l)*par(8)*(1-pi(i))*par(2))^2);
% x(5,i)=par(5)*w*l*(1-w)*(1-l)*par(8)*((par(6)+(par(7)+par(8))*pi(i))+(par(7)+par(8))*(1-pi(i)))*(par(1)-par(2))/((w*l*par(8)*(1-pi(i))*par(1)+(1-w)*(1-l)*(par(6)+(par(7)+par(8))*pi(i))*par(2))^2);

y(1,i)=par(5)*((1-w)*l*(par(6)+(par(7)+par(8))*pi(i))+w*(1-l)*par(8)*(1-pi(i)))/((1-w)*l*(par(6)+(par(7)+par(8))*pi(i))*par(1)+w*(1-l)*par(8)*(1-pi(i))*par(2));
y(2,i)=par(5)*((1-w)*l+w*(1-l))/((1-w)*l*par(1)+w*(1-l)*par(2));
y(3,i)=par(5)*((1-w)*l*par(8)*(1-pi(i))+w*(1-l)*(par(6)+(par(7)+par(8))*pi(i)))/((1-w)*l*par(8)*(1-pi(i))*par(1)+w*(1-l)*(par(6)+(par(7)+par(8))*pi(i))*par(2));
% y(4,i)=par(5)*(1-w)*l*w*(1-l)*par(8)*((par(6)+(par(7)+par(8))*pi(i))+(par(7)+par(8))*(1-pi(i)))*(par(2)-par(1))/(((1-w)*l*(par(6)+(par(7)+par(8))*pi(i))*par(1)+w*(1-l)*par(8)*(1-pi(i))*par(2))^2);
% y(5,i)=par(5)*(1-w)*l*w*(1-l)*par(8)*((par(6)+(par(7)+par(8))*pi(i))+(par(7)+par(8))*(1-pi(i)))*(par(1)-par(2))/(((1-w)*l*par(8)*(1-pi(i))*par(1)+w*(1-l)*(par(6)+(par(7)+par(8))*pi(i))*par(2))^2);

end

RET=zeros(1,3,lp);
for j=1:3
for i=1:lp
    RET(1,j,i)=par(3)-par(4)*x(j,i);
    RET(2,j,i)=par(3)-par(4)*y(j,i);
end
end

V=zeros(2,lp);
M=zeros(1,lp);
B=zeros(2,lp);

for i=1:lp
    M(i)=z*(pi(i)/(1-pi(i)))^alf;

    B(1,i)=((par(5)*(w*l+(1-w)*(1-l)))^(-1))*...
         (w*l*par(1)*(par(6)+(par(7)+par(8))*pi(i))+(1-w)*(1-l)*par(2)*par(8)*(1-pi(i)))*RET(1,1,i)*(RET(1,1,i)>0)...
        +(w*l*par(1)+(1-w)*(1-l)*par(2))*par(7)*(1-pi(i))*RET(1,2,i)*(RET(1,2,i)>0)...
        +(w*l*par(1)*par(8)*(1-pi(i))+(1-w)*(1-l)*par(2)*(par(6)+(par(7)+par(8))*pi(i)))*RET(1,3,i)*(RET(1,3,i)>0);
    V(1,i)=-M(i)+B(1,i);
    
    B(2,i)=((par(5)*((1-w)*l+w*(1-l)))^(-1))*...
         ((1-w)*l*par(1)*(par(6)+(par(7)+par(8))*pi(i))+w*(1-l)*par(2)*par(8)*(1-pi(i)))*RET(2,1,i)*(RET(2,1,i)>0)...
        +((1-w)*l*par(1)+w*(1-l)*par(2))*par(7)*(1-pi(i))*RET(2,2,i)*(RET(2,2,i)>0)...
        +((1-w)*l*par(1)*par(8)*(1-pi(i))+w*(1-l)*par(2)*(par(6)+(par(7)+par(8))*pi(i)))*RET(2,3,i)*(RET(2,3,i)>0);
    V(2,i)=-M(i)+B(2,i);
end

MBh=zeros(1,lp-1);
MBl=zeros(1,lp-1);
for i=2:lp
	MBh(i-1)=(B(1,i)-B(1,i-1))/(pi(i)-pi(i-1));
	MBl(i-1)=(B(2,i)-B(2,i-1))/(pi(i)-pi(i-1));
end
avgMBh=sum(MBh(1:lp-1))/(lp-1);
avgMBl=sum(MBl(1:lp-1))/(lp-1);

[~,mh]=max(V(1,:));
piH=pi(mh);
[~,ml]=max(V(2,:));
piL=pi(ml);

MCh=z*alf*(1-pi(mh))^(-2)*(pi(mh)/(1-pi(mh)))^(alf-1);
MCl=z*alf*(1-pi(ml))^(-2)*(pi(ml)/(1-pi(ml)))^(alf-1);

SPRD=zeros(1,3);
SPRD(:)=y(:,ml)-x(:,mh);
sprdA=y(1,ml)-x(1,mh);
sprdB=y(2,ml)-x(2,mh);
sprdC=y(3,ml)-x(3,mh);

% figure(1)
% hold on
% plot(pi,M,pi,B)
% grid on
% legend('m','b')
% hold off

% figure(2)
% hold on
% plot(pi,SPRD(1,:),pi,SPRD(2,:),pi,SPRD(3,:))
% grid on
% hold off

% figure(3)
% hold on
% plot(pi,V)
% grid on
% plot(pi(mh),V(1,mh),'--gs','MarkerSize',20)
% plot(pi(ml),V(2,ml),'--ro','MarkerSize',20)
% hold off

end