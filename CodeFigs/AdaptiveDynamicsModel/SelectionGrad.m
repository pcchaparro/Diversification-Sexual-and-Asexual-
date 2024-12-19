
THETAF=DTHETA*(1:n); %Optimal trait values for n resources
trait = min(THETAF)-DTHETA:.01:max(THETAF)+DTHETA; %range of traits to estimate selection gradient
str   = size(trait);

AR = zeros(str(1,2),n);

for j=1:str(1,2)
%Attack rates
    AR(j,:) = (AMAX.*exp(-((trait(1,j)-THETAF).^2)./(2*(TAU^2))))'; %attack rate of morph j on preexisting resources
end

%------SELECTION GRADIENT OF SEXUAL POPULATION-------

%ecological equilibrium from analytical solution of sexual population

alphasex = delta.*prod(AR,2);
betasex  = rho*delta.*sum(AR,2) - Ea*rho*FTmax/n.*prod(AR,2);
gammasex = rho^2*delta-.5*Ea*rho^2*FTmax/n.*sum(AR,2);

Nfema = max((-betasex+(betasex.^2-4*alphasex.*gammasex).^(1/2))./(4*alphasex),0);
Nmale = max((-betasex+(betasex.^2-4*alphasex.*gammasex).^(1/2))./(4*alphasex),0);
R1sex = (rho*FTmax/n)./(rho+AR(:,1).*(Nfema+Nmale));
R2sex = (rho*FTmax/n)./(rho+AR(:,2).*(Nfema+Nmale));

%selection gradient from analytical solution
SelGradsex = zeros(str(1,2),1);
for j=1:str(1,2)
    SelGradsex(j,1) = .5*(-(AMAX*Ea*R1sex(j,1)*exp(-(trait(1,j) - THETAF(1,1))^2/(2*TAU^2))*(2*trait(1,j) - 2*THETAF(1,1)))/(2*TAU^2) -(AMAX*Ea*R2sex(j,1)*exp(-(trait(1,j) - THETAF(1,2))^2/(2*TAU^2))*(2*trait(1,j) - 2*THETAF(1,2)))/(2*TAU^2));
end

%remove values where the population is extinct
maxSelS = find(SelGradsex==max(SelGradsex));
SelGradsex(1:maxSelS-1)=[];
minSelS = find(SelGradsex==min(SelGradsex));
SelGradsex(minSelS+1:end)=[];
traitsex=trait;
traitsex(1:maxSelS-1)=[];
traitsex(minSelS+1:end)=[];

%------SELECTION GRADIENT OF ASEXUAL POPULATION-------

%ecological equilibrium from analytical solution of sexual population

alphaasex = delta.*prod(AR,2);
betaasex  = rho*delta.*sum(AR,2) - 2*Ea*rho*FTmax/n.*prod(AR,2);
gammaasex = rho^2*delta - Ea*rho^2*FTmax/n.*sum(AR,2);

Nind   = max((-betaasex+(betaasex.^2-4*alphaasex.*gammaasex).^(1/2))./(2*alphaasex),0);
R1asex = (rho*FTmax/n)./(rho+AR(:,1).*Nind);
R2asex = (rho*FTmax/n)./(rho+AR(:,2).*Nind);

%selection gradient from analytical solution
SelGradasex = zeros(str(1,2),1);
for j=1:str(1,2)
    SelGradasex(j,1) = -(AMAX*Ea*R1asex(j,1)*exp(-(trait(1,j) - THETAF(1,1))^2/(2*TAU^2))*(2*trait(1,j) - 2*THETAF(1,1)))/(2*TAU^2) -(AMAX*Ea*R2asex(j,1)*exp(-(trait(1,j) - THETAF(1,2))^2/(2*TAU^2))*(2*trait(1,j) - 2*THETAF(1,2)))/(2*TAU^2);
end

%remove values where the population is extinct
maxSelA = find(SelGradasex==max(SelGradasex));
SelGradasex(1:maxSelA-1)=[];
minSelA = find(SelGradasex==min(SelGradasex));
SelGradasex(minSelA+1:end)=[];
traitasex=trait;
traitasex(1:maxSelA-1)=[];
traitasex(minSelA+1:end)=[];