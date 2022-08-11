close all
clear all
clc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Initial conditions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rmax = 2;  %Max of reproductions with signaling
rmin = 0.01; %Min no. of reproductions without signaling 
sigmag = 1;
ini_pop_host = 500;
ini_pop_par = 200;
beta = 0.5; % alphaprime
k = 1000; % Carrying capacity 
tmax = 500;  % no of time steps to run 
q = 0.5;% percentage of females
fmax = 3;
zsize = 100;
nz=zeros(1,zsize);
pt=zeros(1,zsize);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setting up Initital Distribution of traits and popultaion vector
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
x = (1: 1: zsize); % Creating a vector from 0 to 100
y = makedist('Normal', 90, 10); %Making a normal distrubution with mean and std dev 
z = pdf(y,x); % Making a probability distribution 
y1 = makedist('Uniform', 0, 100);
z1 = pdf(y,x); 
plot(x,z); % Plotting to check 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% calculate zbar (this will inform R(zbar)) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
zbar = y.mu; %Setting zbar to y.mu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% N(z) i.e what is the no. of individuals with that trait (z)?
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nz(1,:) = ini_pop_host*z; %Putting 1000 individuals across the pdf
pt(1,:) = ini_pop_par*z1;
% figure(1)
% plot(x,nz)
% title('initial distribution of traits')
% xlabel('z')
% ylabel('Density distribution of population (indv/m^2)')
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
nzt = zeros(zsize, tmax); % initializing a population dataframe across time to be filled 
nzt(:,1) = nz; % the first column for t=1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tdel = 1; % if we need to slow down the function 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%time loop 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%if nz > 0   %run loop only if population sum is greater than 0
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%given R(zbar), how many individuals should be born to N(z) ie, babyz
% given S(z), how many individuals die from N(z),  deadz?
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
for   t = 1:tmax %for t = 1 to tmax
if sum(nz) > 0.00000001  % run simulation if total pop size is greater than 0 
   % Probescpar1  = probescpar( x, nz, pt, fmax);
    Rpar  = repropar( x, zbar, rmin, rmax, beta );
%     fitnesshost = (Rpar*(1-sum(nz)/k)).*Probescpar1;
%     fitnesspar = (1-Probescpar1).*q;
%     
%    hostnew = nz.*fitnesshost ;
%     parnew = nz.*fitnesspar;
%     
 nznew= nz.*exp(Rpar.*(1-sum(nz)./k)-(((x./(1+x))).*sum(pt).*fmax)/(fmax+ ((x./(1+x)).*sum(nz))));
 ptnew=nz.*(1-exp(((x./(1+x)).*sum(pt).*fmax)/(fmax+((x./(1+x)).*sum(nz)))))*q; 
% % nznew = nz.*exp(Rpar.*(1- sum(nz)./k)).*exp(-a.*sum(pt));
% ptnew = nz.*(1-exp(-a.*sum(pt)))*c;

    %nznew = ((nz.* (R* (1 - (sum(nz)/k)))) - (S.* nz))*tdel; % New population until carrying capacity
  nz =  nznew ;
  pt = ptnew;% replace old population with new population
%     
  %  nz(nz < 0) = 0; % if individuals in population falls below zero, make it zero 
  
    nzt(:,t)=nz; % storing new population vector as dataframe at every time step
    parat(:,t) = pt;
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%
    %adding mutation/offspring variance 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nzmut = zeros(1,zsize); % 
    for i = 1:zsize % for iteration from 1 to 100 of syllable rate
       indz = nz(i); % picking number of individuals at each trait
       mutdist = pdf(makedist('Normal', i, sigmag),x); %creating normal pdf at each trait
       nzmut = nzmut + (indz * mutdist);  % making a new mutation population       
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    nz = nzmut; %replacing new population with individuals who have variation 
    nzt(:,t)=nz; % storing new population vector as dataframe at every time step


%calculating zbar
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if sum(nz) ~= 0
y.mu = dot(nz/sum(nz),x);
else if sum(nz) == 0
y.mu = 0;
    end
 end
%y.mu = sum (nz.*x)/sum(nz);
%y.mu = dot(nz/sum(nz),x); %dot product to find the mean 
zbar = y.mu; %replacing the old zbar with new pop zbar
zbartime(:,t) = zbar; % Storing zbar over all timesteps 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Plotting new population after variation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%     figure(2)
%     plot(x,fitness)
%     title('per-capita/Malthusian growth')
% xlabel('z - syllable rate')
% ylabel('individual fitness')
%  % ylim([2 4])
%      xlim([0 100])
% %     pause(0.5)
figure(3) 
 plot(x,nz)
title('new distribution of traits')
xlabel('z')
ylabel('Density distribution of population (indv/m^2)')
% pause(.5)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Peak calculation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Plotting 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(4)
plot(sum(nzt(:,:)))
title('Host pop vs time')
xlabel('time')
ylabel('N')

figure(5)
plot(sum(parat(:,:)))
title('Past pop vs time')
xlabel('time')
ylabel('P')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(6)
plot(zbartime)
title('Zbar vs time')
xlabel('time')
ylabel('Zbar')

figure(7)
subplot(1,3,1)
area(nzt(:,1),'FaceColor',[0.7 0.7 0.7])

subplot(1,3,2)
area(nzt(:,10),'FaceColor',[0.7 0.7 0.7])

subplot(1,3,3)
area(nzt(:,36),'FaceColor',[0.7 0.7 0.7])



% plot1 = plot(x,nzt(:,1),'b',x,nzt(:,3),'b',x,nzt(:,10),'b',x,nzt(:,20),'b',x,nzt(:,30),'yellow',x,nzt(:,36),'yellow')
% plot1(1).LineWidth = 2;
% plot1(2).LineWidth = 2;
% plot1(3).LineWidth = 2;
% plot1(4).LineWidth = 2;
% plot1(5).LineWidth = 2;
% plot1(6).LineWidth = 2;
% plot1(5).Marker = '*';
% plot1(6).Marker = '*';