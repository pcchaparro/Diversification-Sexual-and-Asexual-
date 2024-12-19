%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function definitions - Perform Individual Based Model (IBM) simulations
% of a sexual lineage for the paper "Ecological diversification in sexual
% and asexual lineages"
%
% Other m-files required: simulateSex.m, finddad.m
% Subfunctions: none
% MAT-files required: none
%
% Author: Catalina Chaparro
%
%   original version: 10.01.2024,
%   last version: 25.04.2024
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear

%Parameters----------------------------------------------------

% Habitat

n = 10;                   % Number of resources
rho=.01;                  % Resource growth rate
FTmax = 50;               % Total carrying capacity (limit around Rmax=4)
volHab = 1E3;             % Volume of the habitat (too small result in high demographic stocasticity)
DTHETA=1;                 % Trait distance between optimals to feed on different resources

% Demographic rates (assumed to be equal for both lineages)

AMAX=0.2;                 % Maximum attack rate
TAU=1/3;                  % Standard deviation of the Gaussian that determines how attack rate varies with trait
Ea=.6;                    % Efficiency to transform ingested mass into biomass
delta=.02;                % Mortality rate

% Genetic parametes

nEloci = 10;              %number of loci of ecological trait
nAloci = 5;               %number of loci of assortative mating trait
dadcont= 1;              %contribution of father to genetic material of offspring --> probability that one allele of father is inherited by offspring for each locus. When dadcont=1, the offspring receives exacty half alleles from mom and half from dad.
mutrateEl = .001;           %mutation rate ecological loci
mutrateAl = .001;           %mutation rate mating loci
sigmamutEl = .01;         % standard deviation of the mutation distribution
sigmamutAl = .01;         % standard deviation of the mutation distribution

% Simulation parameters
sim         = 1;               %number of simulations
tmax        = 1E6;             %maximum time of a simulation
tmaxafniches= 1E6;             %maximum time to run after all resources have been used

%Initial conditions

niniind = 40;       %size of initial population
traitini= 0.8;      %initial trait of individuals

%Initialize food resources (assumed to be in its carrying capacity at t=0)

Fmaxi = FTmax/n*ones(1,n);  % carrying capacity of each resource
Fi = Fmaxi;

THETAF=DTHETA*(1:n);  % vector of optimal traits to feed on food resources

%Create/initialize population matrix Pop: each row corresponds to each ind

Eloci = 1:2*nEloci;                   %columns 1:2*nEloci are ecological trait loci
Aloci = 2*nEloci+1:2*nEloci+2*nAloci; %columns 2*nEloci+1:2*nEloci+2*nAloci are assortative mating trait
sex   = 2*nEloci+2*nAloci+1;          %column  2*nEloci+2*nAloci+1 is sex
ecotr = 2*nEloci+2*nAloci+2;          %column  2*nEloci+2*nAloci+2 is ecological trait
mattr = 2*nEloci+2*nAloci+3;          %column  2*nEloci+2*nAloci+3 is assortative mating trait
Biores= 2*nEloci+2*nAloci+4;          %column  2*nEloci+2*nAloci+4 is reserve of biomass (when this reserve reaches 1 in a female, it produces 1 offspring)
aFi   = 2*nEloci+2*nAloci+5:2*nEloci+2*nAloci+5+n-1;          %other columns are the attack rate on resource ith

nfemales = floor(niniind/2);
nmales   = niniind-nfemales;

Pop          = [ones(niniind,2*nEloci)*traitini/(2*nEloci) zeros(niniind,2*nAloci)]; %genotypes
Pop(:,sex)   = [ones(nfemales,1); zeros(nmales,1)]; %1 for females, and -1 for males
Pop(:,ecotr) = sum(Pop(:,Eloci),2);
Pop(:,mattr) = mean(Pop(:,Aloci),2);
Pop(:,Biores)= 0;
for j=1:n
    Pop(:,aFi(1,j)) = AMAX*exp(-((Pop(:,ecotr)-THETAF(1,j)).^2)./(2*TAU^2));
end

% Create histogram bins in which data will be save
binsize = .1;
histedges=linspace(min(THETAF)-5,max(THETAF)+5,n*20+1);
shist=size(histedges);
nbins=shist(1,2)-1;

% ----Simulations-----------------------------------

tic

% Initialize arrays to save last resource densities and histograms of
% ecological trait
Fend       = zeros(sim,n);
timeend    = zeros(sim,1);
tfillednich= zeros(sim,1);%Time when all food resources have been used, such that its food density is depleted below half of the carrying capacity
HistEtSex  = zeros(sim,nbins);
HistEtAsex = zeros(sim,nbins);

rng shuffle % create a different seed each time

for i=1:sim 
    
    [Ffinal,distSex,time,timeniches] = simulateSex(tmax,tmaxafniches,Pop,n,rho,volHab,Fi,Fmaxi,THETAF,Ea,AMAX,TAU,delta,nEloci,nAloci,mutrateEl,mutrateAl,sigmamutEl,sigmamutAl,Eloci,Aloci,sex,ecotr,mattr,dadcont,Biores,aFi,histedges);

    Fend(i,:)       = Ffinal;
    HistEtSex(i,:)  = distSex;
    timeend(i,:)    = time;
    tfillednich(i,:)= timeniches;
end
toc

currentdirectory = pwd;
save(['sim' currentdirectory(1,end) '.mat'],'Fend','HistEtSex','HistEtAsex','timeend','tfillednich')
