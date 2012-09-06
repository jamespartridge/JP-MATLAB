function [piH,piL,Rah,Rbh,Rch,Ral,Rbl,Rcl,eq]=AAAeqcalc(w,par)

%pi grid
pi=0:0.01:1;
lp=length(pi);

%Conditional Probability of rating (1=A,2=B,3=C) given type (G or B)
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

x=zeros(3,lp);
y=zeros(3,lp);
for i=1:lp
    %High signal interest rates
    x(1,i)=par(5)*(w*par(12)*G(1,i)+(1-w)*(1-par(12))*B(1,i))/(w*par(12)*G(1,i)*par(1)+(1-w)*(1-par(12))*B(1,i)*par(2));
    x(2,i)=par(5)*(w*par(12)*G(2,i)+(1-w)*(1-par(12))*B(2,i))/(w*par(12)*G(2,i)*par(1)+(1-w)*(1-par(12))*B(2,i)*par(2));
    x(3,i)=par(5)*(w*par(12)*G(3,i)+(1-w)*(1-par(12))*B(3,i))/(w*par(12)*G(3,i)*par(1)+(1-w)*(1-par(12))*B(3,i)*par(2));
    %Low signal interest rates
    y(1,i)=par(5)*((1-w)*par(12)*G(1,i)+w*(1-par(12))*B(1,i))/((1-w)*par(12)*G(1,i)*par(1)+w*(1-par(12))*B(1,i)*par(2));
    y(2,i)=par(5)*((1-w)*par(12)*G(2,i)+w*(1-par(12))*B(2,i))/((1-w)*par(12)*G(2,i)*par(1)+w*(1-par(12))*B(2,i)*par(2));
    y(3,i)=par(5)*((1-w)*par(12)*G(3,i)+w*(1-par(12))*B(3,i))/((1-w)*par(12)*G(3,i)*par(1)+w*(1-par(12))*B(3,i)*par(2));
end

%Return if project pays off at each interest rate (y-DR)
RET=zeros(2,3,lp);
for j=1:3
for i=1:lp
    RET(1,j,i)=par(3)-par(4)*x(j,i);
    RET(2,j,i)=par(3)-par(4)*y(j,i);
end
end

V=zeros(2,lp,lp);
C=zeros(1,lp);
E=zeros(2,lp,lp);
% i is pi chosen, j is market/interest rate pi
for i=1:lp
	C(i)=par(9)*(pi(i)^par(10))/((1-pi(i))^par(11));
% 	C(i)=z*pi(i)/(1-pi(i))^alf;
% 	C(i)=z*(pi(i)^alf)/(1-pi(i));
    for j=1:lp
        E(1,i,j)=((par(5)*(w*par(12)+(1-w)*(1-par(12))))^(-1))*(...
             (w*par(12)*par(1)*G(1,i)+(1-w)*(1-par(12))*par(2)*B(1,i))*RET(1,1,j)*(RET(1,1,j)>0)...
            +(w*par(12)*par(1)*G(2,i)+(1-w)*(1-par(12))*par(2)*B(2,i))*RET(1,2,j)*(RET(1,2,j)>0)...
            +(w*par(12)*par(1)*G(3,i)+(1-w)*(1-par(12))*par(2)*B(3,i))*RET(1,3,j)*(RET(1,3,j)>0));
        V(1,i,j)=-C(i)+E(1,i,j);
    
        E(2,i,j)=((par(5)*((1-w)*par(12)+w*(1-par(12))))^(-1))*(...
             ((1-w)*par(12)*par(1)*G(1,i)+w*(1-par(12))*par(2)*B(1,i))*RET(2,1,j)*(RET(2,1,j)>0)...
            +((1-w)*par(12)*par(1)*G(2,i)+w*(1-par(12))*par(2)*B(2,i))*RET(2,2,j)*(RET(2,2,j)>0)...
            +((1-w)*par(12)*par(1)*G(3,i)+w*(1-par(12))*par(2)*B(3,i))*RET(2,3,j)*(RET(2,3,j)>0));
        V(2,i,j)=-C(i)+E(2,i,j);
    end
end

mxc=zeros(2,lp);...best choice for firms
mxe=zeros(2,lp);...best choice for investors
for i=1:lp
	[~,mxe(1,i)]=max(V(1,i,:));
    [~,mxe(2,i)]=max(V(2,i,:));
end
for j=1:lp
	[~,mxc(1,j)]=max(V(1,:,j));
    [~,mxc(2,j)]=max(V(2,:,j));
end

%find the EQ piH and piL
eqh=0; eql=0;...number of equilibria
for k=1:lp
    mh=0; ml=0;
    if mxe(1,k)==mxc(1,k)
        mh=k;
        eqh=eqh+1;
    end
    if mxe(2,k)==mxe(2,k)
        ml=k;
        eql=eql+1;
    end
end
eq=[eqh eql];

if mh~=0
    piH=pi(mh);
else
    piH=NaN;
end
if ml~=0
    piL=pi(ml);
else
    piL=NaN;
end

%EQ interest rates
if mh~=0
    Rah=x(1,mh); Rbh=x(2,mh); Rch=x(3,mh);
    Ral=y(1,ml); Rbl=y(2,ml); Rcl=y(3,ml);
else
    Rah=0; Rbh=0; Rch=0;
    Ral=0; Rbl=0; Rcl=0;
end

%marginal benefit
% MBh=zeros(1,lp);
% MBl=zeros(1,lp);
% for i=1:lp
%     if (RET(1,1,i)>=0 && RET(1,2,i)>=0 && RET(1,3,i)>=0)
%         MBh(i)=((par(5)*(w*par(12)+(1-w)*(1-par(12))))^(-1))*(...
%              (w*par(12)*par(1)*(par(7)+par(8))-(1-w)*(1-par(12))*par(2)*par(8))*RET(1,1,i)...
%             -(w*par(12)*par(1)+(1-w)*(1-par(12))*par(2))*par(7)*RET(1,2,i)...
%             +(w*par(12)*par(1)*(-par(8))+(1-w)*(1-par(12))*par(2)*(par(7)+par(8)))*RET(1,3,i));
%     elseif (RET(1,1,i)>=0 && RET(1,2,i)>=0 && RET(1,3,i)<0)
%         MBh(i)=((par(5)*(w*par(12)+(1-w)*(1-par(12))))^(-1))*(...
%              (w*par(12)*par(1)*(par(7)+par(8))-(1-w)*(1-par(12))*par(2)*par(8))*RET(1,1,i)...
%             -(w*par(12)*par(1)+(1-w)*(1-par(12))*par(2))*par(7)*RET(1,2,i));
%     elseif (RET(1,1,i)>=0 && RET(1,2,i)<0 && RET(1,3,i)<0)
%         MBh(i)=((par(5)*(w*par(12)+(1-w)*(1-par(12))))^(-1))*(...
%              (w*par(12)*par(1)*(par(7)+par(8))-(1-w)*(1-par(12))*par(2)*par(8))*RET(1,1,i));
%     elseif (RET(1,1,i)<0 && RET(1,2,i)<0 && RET(1,3,i)<0)
%         MBh(i)=0;
%     end            
%             
%     if (RET(2,1,i)>=0 && RET(2,2,i)>=0 && RET(2,3,i)>=0)
%         MBl(i)=((par(5)*((1-w)*par(12)+w*(1-par(12))))^(-1))*(...
%              ((1-w)*par(12)*par(1)*(par(7)+par(8))-w*(1-par(12))*par(2)*par(8))*RET(2,1,i)...
%             -((1-w)*par(12)*par(1)+w*(1-par(12))*par(2))*par(7)*RET(2,2,i)...
%             +((1-w)*par(12)*par(1)*(-par(8))+w*(1-par(12))*par(2)*(par(7)+par(8)))*RET(2,3,i));
%     elseif (RET(2,1,i)>=0 && RET(2,2,i)>=0 && RET(2,3,i)<0)
%         MBl(i)=((par(5)*((1-w)*par(12)+w*(1-par(12))))^(-1))*(...
%              ((1-w)*par(12)*par(1)*(par(7)+par(8))-w*(1-par(12))*par(2)*par(8))*RET(2,1,i)...
%             -((1-w)*par(12)*par(1)+w*(1-par(12))*par(2))*par(7)*RET(2,2,i));...
%     elseif (RET(2,1,i)>=0 && RET(2,2,i)<0 && RET(2,3,i)<0)
%         MBl(i)=((par(5)*((1-w)*par(12)+w*(1-par(12))))^(-1))*(...
%              ((1-w)*par(12)*par(1)*(par(7)+par(8))-w*(1-par(12))*par(2)*par(8))*RET(2,1,i));...
%     elseif (RET(2,1,i)<0 && RET(2,2,i)<0 && RET(2,3,i)<0)
%         MBl(i)=0;
%     end            
% end

%EQ marginal benefit
% eqMBh=MBh(mh); eqMBl=MBl(ml);

%EQ marginal costs
% MCh=z*((1-pi(i))^alf+alf*(1-pi(i))^(alf-1)*pi(i))/((1-pi(i))^(2*alf));
% MCl=z*((1-pi(i))^alf+alf*(1-pi(i))^(alf-1)*pi(i))/((1-pi(i))^(2*alf));
% MCh=z*alf*pi(mh)^(alf-1)*(1-pi(mh))+pi(mh)^alf/(1-pi(mh))^2;
% MCl=z*alf*pi(ml)^(alf-1)*(1-pi(ml))+pi(ml)^alf/(1-pi(ml))^2;
% MCh=par(9)*(par(10)*(pi(mh)^(par(10)-1))*((1-pi(mh))^par(11))+par(11)*(pi(mh)^par(10))*((1-pi(mh))^(par(11)-1)))/((1-pi(mh))^(2*par(11)));
% MCl=par(9)*(par(10)*(pi(ml)^(par(10)-1))*((1-pi(ml))^par(11))+par(11)*(pi(ml)^par(10))*((1-pi(ml))^(par(11)-1)))/((1-pi(ml))^(2*par(11)));

end