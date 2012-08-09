MB=zeros(2,lp-1);
for i=2:lp
	MB(1,i-1)=(B8(1,i)-B8(1,i-1))/(pi(i)-pi(i-1));
	MB(2,i-1)=(B8(2,i)-B8(2,i-1))/(pi(i)-pi(i-1));
end
avgMBh=sum(MB(1,1:lp-1))/(lp-1)
avgMBl=sum(MB(2,1:lp-1))/(lp-1)
