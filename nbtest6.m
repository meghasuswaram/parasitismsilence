close all 
clear all
clc
%Initial conditions
rmin = 0.01; %Min no. of reproductions without signaling 
rmax=2;
ini_pop_host = 500;
ini_pop_par = 200;
beta = 0.5; % alphaprime
k = 1000; % Carrying capacity 
tmax = 100;  % no of time steps to run 
qi = [0.1:0.03:0.9];
fmaxi = [1:0.15:5];
rminoutzbarmean = zeros(27,27);
rminoutn = zeros(27,27);
rminout2zbar = zeros(27,27);
rminout2n = zeros(27,27);
temp_index = 1;
sigmag = 1;
for i = 1: length(qi)
    q = qi(i)
for j = 1: length(fmaxi)
    fmax = fmaxi(j)

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
    
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
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
peak1a(j,i) = {sum(nzt(:,:))};
S = stepinfo(nz);
S1 = sum(nzt(:,tmax)); %N*
S3 = zbartime(:,end);
rminoutn(j,i) = S1;%N*
rminoutzbarmean(j,i) = S3;%zbar*
end
end



figure(6)
contourf(qi, fmaxi, rminoutzbarmean,'EdgeColor','none')
mymap1 = [0.7 0.85 0.925
%     0.6 0.8 0.910
    0.5 0.75 0.895
    0.4 0.7 0.870
    0.3 0.6 0.850
    0.2 0.55 0.845
    0.1 0.5 0.8410
   0 0.4470 0.7410];
colormap(mymap1)
colorbar
grid off
% ylabel({'Number of viable parasitoid offspring'})
% xlabel('Proportion of parasitoid females')
% title({'Chorus mean of the population'})
% 
figure(8)
contourf(qi,fmaxi, rminoutn,'EdgeColor','none')
mymap = [0 0 0
    0.4 0.4 0.4
    0.5 0.5 0.5
    0.6 0.6 0.6
    0.7 0.7 0.7
    0.8 0.8 0.8
    0.9 0.9 0.9
    ];
colormap(mymap)
colorbar
grid off
% ylabel({'Number of viable parasitoid offspring'})
% xlabel('Proportion of parasitoid females')
% title({'Population size of the host'})


peak1b = cellfun(@(peak1a) peak1a(end-20:end),peak1a,'UniformOutput',false);
peak1c = cellfun(@(peak1b) round(peak1b,-2),peak1b,'UniformOutput',false);
peak1d = cellfun(@(peak1c) unique(peak1c),peak1c,'UniformOutput',false);
peak1e = cellfun(@(peak1d) length(peak1d),peak1d,'UniformOutput',false);
peak1f = cell2mat(peak1e);

peak1f(peak1f==1) = 0;

peak1f(peak1f>=2 & peak1f <=3) = 1;

peak1f(peak1f>=3) = 2;



figure(7)
h = heatmap(qi, fmaxi, peak1f,'CellLabelColor','none');
h.YDisplayData = flipud(h.YDisplayData);  % equivalent to 'YDir', 'Reverse'
%colorbar
ylabel({'Number of viable parasitoid offspring'})
xlabel('Proportion of parasitoid females')
title({'Population steady state';'Stability = 0, Period doubling =1, Chaos = 2'})
grid off



