%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function definitions - Make figure 2 and 3 of the paper "Ecological diversification
% in sexual and asexual lineages"
%
% Other m-files required: Eco_dynamics.m, SelectionGrad.m,
% EcoevoDynAsexual.m, EcoevoDynSexual.m, EcoEquiSexual.m, EcoEquiAsexual.m
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

% Parameters
rho=.01;                  % Resource growth rate
AMAX=0.2;                 % Maximum attack rate
TAU=1/3;                    % Standard deviation of the Gaussian that determines how attack rate varies with trait
DTHETA=1;                 % Trait distance between optimals to feed on different resources
Ea=.6;                    % Efficiency to transform ingested mass into biomass
delta=.02;                % Mortality rate
mu = .001;                % mutation rate
sigma = 1E-6;              % variance of mutational steps

figure
suptitle('Figure 2')

%Parameters
n=2;                     % Number of resources
FTmax=10;                % Total carrying capacity

%Simulate
Eco_dynamics

%Plot
plot(t1_a,x1_a(:,3),'b','LineWidth',2)
hold on
plot(t2_a,x2_a(:,3),'b--','LineWidth',2)
plot(t_s,x_s(:,3)+x_s(:,4),'r','LineWidth',1)
ylabel('Population density')
xlabel('Time')
legend('Asexual (high feeding efficiency)', 'Asexual (low feeding efficiency)','Sexual')


figure
suptitle('Figure 3')

%---Figure 3A---

% Specific parameters
n=2;                     % Number of resources
FTmax=10;                % Total carrying capacity

% Simulate

SelectionGrad

% Plot
subplot(1,3,1)
plot(traitasex,SelGradasex,'b')
hold on
plot(traitsex,SelGradsex,'r')
plot(trait,zeros(str(1,2),1),'k')
xlabel('Feeding niche trait')
ylabel('Selection gradient')
legend('Asexual','Sexual')

%---Figure 3B---

% Simulate asexual lineage

EcoevoDynAsexual

%Plot

subplot(1,3,2)
plot(1:tspeciationAsex{1},TTraiAsex{1},'b')
ylim([0 2])
xlabel('Time')
ylabel('Feeding niche trait') 

% Simulate sexual lineage

EcoevoDynSexual

%Plot
hold on
plot(1:tspeciationSex{1},TTraiSex{1},'r')

%---Figure 3C---

% Specific parameters
n=10;                     % Number of resources
FTmax=50;                 % Total carrying capacity

% Simulate asexual lineage

EcoevoDynAsexual

%Plot

subplot(1,3,3)
tfig=0;
timefig=[];
Nichesoccupied=[];
for i=1:cont
    timefig=[timefig tfig+(1:tspeciationAsex{i})];
    Nichesoccupied=[Nichesoccupied; TNPopAsex{i}];
    tfig = tfig + tspeciationAsex{i};
end
plot(timefig,Nichesoccupied,'b')
ylim([0 11])
xlabel('Time')
ylabel('Number of niches occuppied')  

% Simulate sexual lineage

EcoevoDynSexual

% Plot

hold on
tfig=0;
timefig=[];
Nichesoccupied=[];
for i=1:cont
    timefig=[timefig tfig+(1:tspeciationSex{i})];
    Nichesoccupied=[Nichesoccupied; TNPopSex{i}];
    tfig = tfig + tspeciationSex{i};
end
plot(timefig,Nichesoccupied,'r')
ylim([0 11])
xlabel('Time')
ylabel('Number of niches occuppied') 