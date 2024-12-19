% This rutine plots the demographic dynamics of an asexual and a sexual
% population for a fixed trait value


THETAF=DTHETA*(1:n); %Optimal trait values for n resources
m=1;

%Calculate attack rates on each resource
AR   = zeros(1,n);
Trai = 1.4;
for j=1:m
%Attack rates
    AR(j,:) = AMAX.*exp(-((Trai-THETAF).^2)./(2*(TAU^2))); %attack rate of morph j on preexisting resources
end

%-----Demographic dynamics using ODE solver

%find the ecological equilibrium for an asexual population
initialFood = FTmax/n*ones(n,1);
initialBiom = .1;

tmaxe = 5000;
x0    = [initialFood' initialBiom'];
fod   = @(t,x) EcoEquiAsexual(t,x,AR,n,m,rho,FTmax,Ea,delta);
[t1_a,x1_a] = ode45(fod,[0 tmaxe],x0);

%find the ecological equilibrium for an asexual population with half Ea

tmaxe = 5000;
x0    = [initialFood' initialBiom'];
fod   = @(t,x) EcoEquiAsexual(t,x,AR,n,m,rho,FTmax,Ea/2,delta);
[t2_a,x2_a] = ode45(fod,[0 tmaxe],x0);

%find the ecological equilibrium for a sexual population
initialFood = FTmax/n*ones(n,1);
initialFema = .05;
initialMale = .05;

tmaxe = 5000;
x0    = [initialFood' initialFema' initialMale'];
fod   = @(t,x) EcoEquiSexual(t,x,AR,n,m,rho,FTmax,Ea,delta);
[t_s,x_s] = ode45(fod,[0 tmaxe],x0);

