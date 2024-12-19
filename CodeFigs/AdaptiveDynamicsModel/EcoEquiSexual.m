function A = EcoEquiSexual(t,x,AR,n,m,rho,FTmax,Ea,delta)

    %Ecological dynamics
    % A is the matrix of attack rates of all ecomorphs on resources.
    % Rows correspond to each consumer ecomorph and columns
    % to the attacked resource.
    % x is the vector of variables to solve
    
A = zeros(n+m,1);

Food=x(1:n,1);
Fema=x(n+1:n+m,1);
Male=x(n+m+1:n+2*m,1);

ingestF = zeros(m,1);
totingest= zeros(m,1);

Femachange = zeros(m,1);
Malechange = zeros(m,1);
Foodchange = zeros(n,1);

for j=1:m %number of morphs
    %feeding on basal resources
    for i=1:n %number of resources
        ingestF(j,1) = ingestF(j,1) + AR(j,i)* Food(i,1); 
    end
    totingest(j,1) = ingestF(j,1);
    
    Femachange(j,1) = (1/2*Ea*totingest(j,1)-delta)*Fema(j,1);
    Malechange(j,1) = 1/2*Ea*totingest(j,1)*Fema(j,1)-delta*Male(j,1);
end

Fmax=FTmax/n;
grazing= zeros(n,1);
for i=1:n %number of food resources
    for j=1:m %number of morphs
        grazing(i,1) = grazing(i,1) + AR(j,i)*Food(i,1)*(Fema(j,1)+Male(j,1));
    end
    Foodchange(i,1) = rho*(Fmax-Food(i,1)) - grazing(i,1);
end

A(1:n,1)=Foodchange;
A(n+1:n+m,1)=Femachange;
A(n+m+1:n+2*m,1)=Malechange;
end