function dTr = EvoDynAsexual(t,Tr,n,m,Food,Nind,Ea,AMAX,THETAF,tau,mu,sigma)

Trait=Tr;

Traitchange = zeros(m,1);
fGradF = zeros(m,1);
FitnessGrad = zeros(m,1);

for j=1:m
    %calculate fitness gradient (formula from symbolic solver
    %FitnessGrad.m)
    
    %From basal resources
    for i=1:n %number of resources
        fGradF(j,1) = fGradF(j,1) -(AMAX*Ea*Food(1,i)*exp(-(Trait(j,1) - THETAF(1,i))^2/(2*tau^2))*(2*Trait(j,1) - 2*THETAF(1,i)))/(2*tau^2);
    end
    
    FitnessGrad(j,1) = fGradF(j,1);
    Traitchange(j,1) = sigma*FitnessGrad(j,1); %*.5*mu*sigma*Nind(j,1)
end

dTr = Traitchange;
end