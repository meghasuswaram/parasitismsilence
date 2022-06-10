close all
clear all
clc
rmin = 0.01;
rmax = 0.7;
k = 1;
ci = [0.1:0.1:1];
tmax = 50;
zsize = 100;
x = (1: 1: zsize); % Creating a vector from 0 to 100
y = makedist('Normal', 80, 10); %Making a normal distrubution with mean and std dev 
z = pdf(y,x); % Making a probability distribution 
zbar = y.mu; %Setting zbar to y.mu
h(1,:) = 0.9*z;
figure(1)
plot(x,h)
p = zeros(1,tmax);
p(1) = 0;
cioutN = zeros(1,10);
cioutzbar = zeros(1,10);
for i = 1:length(ci)
    c = ci(i)
for t = 1: tmax-1
    a = (x./(50+x));
 Rhost = rmin + (rmax - rmin)./(1 + exp(-((x - zbar))));
    h(t+1)=h(t).*exp((sum(Rhost)).*(1-h(t)./k)-(sum(a)).*p(t) );
    p(t+1)=c.*h(t).*(1-exp(-(sum(a)).*p(t)));
    
    if sum(h) ~= 0
y.mu = dot(h/sum(h),x);
else if sum(h) == 0
y.mu = 0;
    end
    end
 zbar = y.mu; %replacing the old zbar with new pop zbar
zbartime(:,t) = zbar;
    
    ht(:,t) = h;
        
%    figure(2)
%    plot(x,h)
%    title('new dist of traits')
   
   
end
S1 = sum(ht(:,end)); %N*
S2 = zbartime(:,end); %zbar*
cioutN(i) = S1;%N*
cioutzbar(i) = S2;%zbar*
end
 
figure(3)
     plot((sum(ht(:,:))),'b'),hold on, plot(p,'r'),hold off
     legend('host','parasite')
     title('host and par pop over time')
     
     figure(4)
     plot(zbartime)
     title('mean signal over time')
     
     figure(5)
     scatter(ci, cioutN)
      title('steady state of pop vs viable eggs')
      figure(6)
     scatter(ci, cioutzbar)
      title('steady state of trait vs viable eggs')
     