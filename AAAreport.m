tic
clear all

gamma_G=0.90;    ...prob. pays off if good
gamma_B=0.30;    ...prob. pays off if bad
y=1.5;           ...output
D=1;             ...investment size
r=1.01;          ...risk free rate
p1=0.5;          ...Pr(A|H)=Pr(C|L)=p1+(p2+p3)pi
p2=0.3;          ...Pr(B)=p2(1-pi)
p3=0.2;          ...Pr(C|H)=Pr(A|L)=p3(1-pi)

z=1/13;
alf=1;
delta=1.2;

l=0.7;

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

%omega grid
w=0.5:0.1:1;
% w=0.5:0.05:1;
lw=length(w);

%initialize matrix sizes
piH=zeros(1,lw);
piL=zeros(1,lw);
Rh=zeros(3,lw);
Rl=zeros(3,lw);
eq=zeros(2,lw);

%calculate EQ values at each omega
for k=1:lw
    [piH(k),piL(k),Rh(1,k),Rh(2,k),Rh(3,k),Rl(1,k),Rl(2,k),Rl(3,k),...
        eq(:,k)]=AAAeqcalc(w(k),par);
%multiple equilibria?
    if eq(2,k)==0
        display(['No Equilibrium at ' num2str(w(k))])
    elseif eq(2,k)>1
        display([num2str(eq(2,k)) ' Equilibria at ' num2str(w(k))])
    elseif eq(2,k)==1
        display(['Equilibrium at ' num2str(w(k))])
    end
end


%spread within rating at the EQ pi
% RateSprd(1,:)=Rl(1,:)-Rh(1,:);
% RateSprd(2,:)=Rl(2,:)-Rh(2,:);
% RateSprd(3,:)=Rl(3,:)-Rh(3,:);

%spread within signal at the EQ pi
% SigSprd(1,1,:)=Rh(3,:)-Rh(1,:); ...H signal, C over A
% SigSprd(1,2,:)=Rh(2,:)-Rh(1,:); ...H signal, B over A
% SigSprd(2,1,:)=Rl(3,:)-Rl(1,:); ...L signal, C over A
% SigSprd(2,2,:)=Rl(2,:)-Rl(1,:); ...L signal, B over A

%EQ rating probabilities conditional on signal
% probAgivH=par(6)+(par(7)+par(8))*piH;
% probBgivH=par(7)*(1-piH);
% probCgivH=par(8)*(1-piH);
% probAgivL=par(8)*(1-piL);
% probBgivL=par(7)*(1-piL);
% probCgivL=par(6)+(par(7)+par(8))*piL;

% numA=probAgivH*l+probAgivL*(1-l);
% numB=probBgivH*l+probBgivL*(1-l);
% numC=probCgivH*l+probCgivL*(1-l);

%Exp. return in EQ
% RETh=zeros(3,lw);
% RETl=zeros(3,lw);
% for i=1:lw
%     for j=1:3
%         RETh(j,i)=par(3)-par(4)*Rh(j,i);
%         RETl(j,i)=par(3)-par(4)*Rl(j,i);
%     end
% end

% figure(1)
% subplot(2,2,1,'replace')
% plot(w,piH,'--',w,piL)
% title('Equilibrium investment in ratings')
% xlabel('\omega')
% ylabel('\pi^*_\nu')
% legend('\pi^*_H','\pi^*_L','Location','best')
% 
% figure(3)
% subplot(2,2,2,'replace')
% plot(w,RateSprd(1,:),'--g',w,RateSprd(2,:),'-b',w,RateSprd(3,:),'r')
% title('Interest rate spreads (R(h,L)-R(h,H))')
% xlabel('\omega')
% legend('A','B','C','Location','best')

% subplot(2,2,3,'replace')
% plot(w,MBh,'--',w,MBl)
% title('Marginal Benefit')
% xlabel('\omega')
% legend('MB H','MB L','Location','best')
% % plot(w,RETh(1,:),'--',w,RETh(2,:),'-',w,RETh(3,:))
% % title('Expected Return in EQ for H')
% % xlabel('\omega')
% % legend('A','B','C','Location','SouthWest')
% 
% subplot(2,2,4,'replace')
% plot(w,MCh,'--',w,MCl)
% title('Marginal Cost')
% xlabel('\omega')
% ylim=([-0.5 1]);
% legend('MC H','MC L','Location','best')
% % plot(w,RETl(1,:),'--',w,RETl(2,:),'-',w,RETl(3,:))
% % title('Expected Return in EQ for L')
% % xlabel('\omega')
% % legend('A','B','C','Location','SouthWest')
 
% figure(2)
% hold on
% h1=plot(w,probAgivH,'b'); h1b=plot(w(1:10:end),probAgivH(1:10:end),'xb');
% h2=plot(w,probBgivH,'c'); h2b=plot(w(1:10:end),probBgivH(1:10:end),'xc');
% h3=plot(w,probCgivH,'m'); h3b=plot(w(1:10:end),probCgivH(1:10:end),'mx');
% h4=plot(w,probCgivL,'r'); h4b=plot(w(1:10:end),probCgivL(1:10:end),'xr');
% h5=plot(w,probBgivL,'g'); h5b=plot(w(1:10:end),probBgivL(1:10:end),'gx');
% h6=plot(w,probAgivL,'y'); h6b=plot(w(1:10:end),probAgivL(1:10:end),'yx');
% title('Prob of Rating, Conditional on Signal')
% xlabel('\omega')
% legend([h1 h2 h3 h4 h5 h6],'AH','BH','CH','CL','BL','AL','Location','best')
% hold off
% 
% figure(3)
% subplot(2,2,[1 2])
% plot(w,Rh(1,:),w,Rh(2,:),w,Rh(3,:))
% title('Interest Rates for H')
% xlabel('\omega')
% legend('A','B','C','Location','best')
% subplot(2,2,[3 4])
% plot(w,Rl(1,:),w,Rl(2,:),w,Rl(3,:))
% title('Interest Rates for L')
% xlabel('\omega')
% legend('A','B','C','Location','best')

% pi=0:0.001:1;
% lp=length(pi);

% for i=1:lp
%     MC(i)=par(9)*(par(10)*(pi(i)^(par(10)-1))*((1-pi(i))^par(11)))+(par(11)*(pi(i)^par(10))*((1-pi(i))^(par(11)-1)))/((1-pi(i))^(2*par(11)));
% for k=1:lw
%     if MBH(k,i)<0
%         MBH(k,i)=-1;
%     end
%     if MBL(k,i)<0
%         MBL(k,i)=-1;
%     end
% end
% end

% figure(4)
% hold on
% plot(pi(1:450),MC(1:450),pi(1:450),MBH(250,1:450),pi(1:450),MBH(350,1:450),...
%                          pi(1:450),MBL(250,1:450),pi(1:450),MBL(350,1:450))

% plot(pi(1:450),MC(1:450),pi(1:450),MBH(3,1:450),'--',pi(1:450),MBH(5,1:450),'--',pi(1:450),MBH(7,1:450),'--',pi(1:450),MBH(9,1:450),'--',...
%                          pi(1:450),MBL(3,1:450),pi(1:450),MBL(5,1:450),pi(1:450),MBL(7,1:450),pi(1:450),MBL(9,1:450))
% legend('MC','MBh w=0.6','MBh w=0.7','MBh w=0.8','MBh w=0.9',...
%             'MBl w=0.6','MBl w=0.7','MBl w=0.8','MBl w=0.9')
% hold off

% figure(5)
% plot(w,numA,w,numB,w,numC)
% legend('A','B','C')

toc