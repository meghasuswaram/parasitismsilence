close all 
clear all
clc
%Initial conditions
%Max of reproductions with signaling
rmin = 0.01; %Min no. of reproductions without signaling 
q = 0.1;
fmax = 1;
ini_pop_host = 500;
ini_pop_par = 200;
beta = 0.5; % alphaprime
k = 1000; % Carrying capacity 
tmax = 100;  % no of time steps to run 
rmaxi = [2:0.001:9];
rminoutzbar = zeros(7001,7001);
rminoutn = zeros(7001,7001);
rminout2zbar = zeros(7001,7001);
rminout2n = zeros(7001,7001);



for j = 1: length(rmaxi)
    rmax = rmaxi(j)
    
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
%plot(x,z); % Plotting to check 
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
    Probescpar1  = probescpar( x, nz, pt, fmax);
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
%

end
end

S1 = sum(nzt(:,end)); %N*
S2 = zbartime(:,end); %zbar*
rminoutn(:,j) = S1;%N*
rminoutzbar(:,j) = S2;%zbar*]
peak1a(:,j) = (sum(nzt(:,:)));
peak3a(:,j) = zbartime(:,end);
peak1b(:,j) = peak1a(end-10:end);
if length(peak3a) > 10
peak3b{j} = peak3a(end-10:end);
else
    peak3b{j} = peak3a(end);
end
peak1c(:,j) = round(peak1b(:,j));
peak3c{j} = round(peak3b{j},4);
peak1d{j} = unique(peak1c(:,j));
peak3d{j} = unique(peak3c{j});

% peak1(:,:) = findpeaks(peak1b(:,j))
% valley1(:,:) = islocalmin(peak1a(end-2:end));
% uniquex(:,:) = unique(peak1);
% uniquey(:,:) = unique(valley1);

end

rows = cellfun(@numel,peak1d);
cols = size(peak1d,2);
peak1e = zeros(max(rows),cols);
for k = 1:cols
    peak1e(1:rows(k),k) = peak1d{k};
end
peak1e(peak1e==0) = nan;

figure(4)
s1 = scatter(rmaxi, peak1e,'w','.','MarkerEdgeColor',[0.75 0.75 0.75],'DisplayName','F_{max} = 1')
hold on 

clear all
clc
%Initial conditions
%Max of reproductions with signaling
rmin = 0.01; %Min no. of reproductions without signaling 
q = 0.4;
fmax = 1;
sigmag = 1
ini_pop_host = 500;
ini_pop_par = 200;
beta = 0.5; % alphaprime
k = 1000; % Carrying capacity 
tmax = 100;  % no of time steps to run 
rmaxi = [2:0.001:7.5];
rminoutzbar = zeros(5501,5501);
rminoutn = zeros(5501,5501);
rminout2zbar = zeros(5501,5501);
rminout2n = zeros(5501,5501);



for j = 1: length(rmaxi)
    rmax = rmaxi(j)
    
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
%plot(x,z); % Plotting to check 
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
    Probescpar1  = probescpar( x, nz, pt, fmax);
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
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     nzmut = zeros(1,zsize); % 
%     for i = 1:zsize % for iteration from 1 to 100 of syllable rate
%        indz = nz(i); % picking number of individuals at each trait
%        mutdist = pdf(makedist('Normal', i, sigmag),x); %creating normal pdf at each trait
%        nzmut = nzmut + (indz * mutdist);  % making a new mutation population       
%     end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     nz = nzmut; %replacing new population with individuals who have variation 
    nzt(:,t)=nz; % storing new population vector as dataframe at every time step

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
%

end
end

S1 = sum(nzt(:,end)); %N*
S2 = zbartime(:,end); %zbar*
rminoutn(:,j) = S1;%N*
rminoutzbar(:,j) = S2;%zbar*]
peak1a(:,j) = (sum(nzt(:,:)));
peak3a(:,j) = zbartime(:,end);
peak1b(:,j) = peak1a(end-10:end);
if length(peak3a) > 10
peak3b{j} = peak3a(end-10:end);
else
    peak3b{j} = peak3a(end);
end
peak1c(:,j) = round(peak1b(:,j));
peak3c{j} = round(peak3b{j},4);
peak1d{j} = unique(peak1c(:,j));
peak3d{j} = unique(peak3c{j});

% peak1(:,:) = findpeaks(peak1b(:,j))
% valley1(:,:) = islocalmin(peak1a(end-2:end));
% uniquex(:,:) = unique(peak1);
% uniquey(:,:) = unique(valley1);

end

rows = cellfun(@numel,peak1d);
cols = size(peak1d,2);
peak1e = zeros(max(rows),cols);
for k = 1:cols
    peak1e(1:rows(k),k) = peak1d{k};
end
peak1e(peak1e==0) = nan;

s2 = scatter(rmaxi, peak1e,	'w','.','MarkerEdgeColor',[0.5 0.5 0.5],'DisplayName','F_{max} = 3')
xlabel({'Increase in reproductive incentive ==>';'r_{max}'})
ylabel('Population size of the host')
title('Bifurcation')











