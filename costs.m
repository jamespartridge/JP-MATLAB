% Cost Functions

pi=0:0.001:0.9;
lp=length(pi);

delta=0.2:0.2:2;
leg_del=cellstr(num2str(delta'));
ld=length(delta);

alf=3;
leg_alf=cellstr(num2str(alf'));
la=length(alf);

z=1/10;
leg_z=cellstr(num2str(z'));
lz=length(z);

c_delt=zeros(lp,ld);
c_alf=zeros(lp,la,lz);
c_delt_p=zeros(lp,ld);
c_alf_p=zeros(lp,la,lz);
c=zeros(lp,ld,la,lz);
c_p=zeros(lp,ld,la,lz);

for i=1:lp
    for j=1:ld
        for k=1:la
            for l=1:lz
%                 c_delt(i,j)=z(l)*(pi(i)^delta(j))/(1-pi(i));
%                 c_delt_p(i,j)=z(l)*delta(j)*pi(i)^(delta(j)-1)*(1-pi(i))+pi(i)^delta(j)/(1-pi(i))^2;
        
                c(i,j,k,l)=z(l)*(pi(i)^alf(k))/((1-pi(i))^delta(j));
                c_p(i,j,k,l)=z(l)*(alf(k)*(pi(i)^(alf(k)-1))*((1-pi(i))^delta(j))+delta(j)*pi(i)^alf(k)*(1-pi(i))^(delta(j)-1))/((1-pi(i))^(2*delta(j)));

%                 c_alf(i,m,n)=z(l)*pi(i)/(1-pi(i))^alf(k);
%                 c_alf_p(i,m,n)=z(l)*((1-pi(i))^alf(k)+alf(k)*(1-pi(i))^(alf(k)-1)*pi(i))/((1-pi(i))^(2*alf(k)));
            end
        end
    end
end

c_pd=squeeze(c_p(:,:,1,1));
% c_pa=squeeze(c_p(:,10,:,1));
% c_pz=squeeze(c_p(:,1,1,:));

figure(1)
plot(pi,c_pd);
legend(leg_del,'Location','NorthWest')
% 
% figure(2)
% plot(pi,c_pa);
% legend(leg_alf,'Location','NorthWest')

% figure(3)
% plot(pi,c_pz);
% legend(leg_z,'Location','NorthWest')

% figure(1)
% plot(c_delt)
% legend('0.1','0.3','0.5','0.7','Location','NorthWest')
% 
% figure(2)
% plot(c_delt_p)
% legend('0.1','0.3','0.5','0.7','Location','NorthWest')
% 
% figure(3)
% plot(c_alfm)
% legend(leg_bet,'Location','NorthWest')
% 
% figure(4)
% plot(c_alfm_p)
% legend('1.3','1.5','1.7','1.9','Location','NorthWest')