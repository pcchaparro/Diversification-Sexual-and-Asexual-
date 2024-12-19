
% This rutine computes the eco-evo dynamics of a sexually reproducing
% lineage when an ancestral population colonizes a habitat with n food resources

%-----ECO-EVO DYNAMICS SEXUAL POPULATION-----

%Initial conditions
m=1;
initialFood = FTmax/n*ones(n,1);
initialNFem = 0.05;                %Intial number of females
initialNMal = 0.05;                %Intial number of males
intialTraitSexual = 0.6;

%Initialize matrices to save info simulation

tmax1 = 1000;
FoodSex = zeros(n,tmax1);
NFem = zeros(m,tmax1);
NMal = zeros(m,tmax1);
TraiSex = zeros(m,tmax1);

FoodSex(:,1) = initialFood;
NFem(:,1) = initialNFem;
NMal(:,1) = initialNMal;
TraiSex(:,1) = intialTraitSexual;

ARsexual = zeros(m,n);
traitchange = 1;
time=1;

while (time<tmax1-1 && traitchange>1E-8)
    
    for j=1:m
    %Calculate attack rates on each resource
    ARsexual(j,1:end) = AMAX.*exp(-((TraiSex(j,time)-THETAF).^2)./(2*(TAU^2))); %attack rate of morph j on preexisting resources
    end    
    
    %find the ecological equilibrium
    tmaxe = 100;
    dif = 1;
    rep = 0;
    x0    = [FoodSex(:,time)' NFem(:,time)' NMal(:,time)'];
    while (dif>1E-6 && rep<20)
        fod = @(t,x) EcoEquiSexual(t,x,ARsexual,n,m,rho,FTmax,Ea,delta);
        [t1,x1] = ode45(fod,[0 tmaxe],x0);
        dif = abs(x1(end-1,1)-x1(end,1));
        rep = rep + 1;
        x0  = x1(end,:);
    end

    FoodSex(:,time+1)=x1(end,1:n);
    NFem(:,time+1)=x1(end,n+1:n+m);
    NMal(:,time+1)=x1(end,n+m+1:n+2*m);
    
    fod   = @(t,Tr)EvoDynSexual(t,Tr,n,m,FoodSex(:,time+1)',NFem(:,time+1)',NMal(:,time+1)',Ea,AMAX,THETAF,TAU,mu,sigma);
    [tevo,xevo] = ode23(fod,[0 1E6],TraiSex(:,time));
    TraiSex(:,time+1) = xevo(end,:);
    
    traitvar = zeros(m,1);
    for j=1:m
        traitvar(j,1) = abs(TraiSex(j,time+1)-TraiSex(j,time));
    end
    traitchange = max(traitvar);
    time = time+1;
end

TFoodSex{1} = FoodSex(:,1:time);
TNFemSex{1} = NFem(:,1:time);
TNMalSex{1} = NMal(:,1:time);
TTraiSex{1} = TraiSex(:,1:time);
TNPopSex{1} = m*ones(time,1);
tspeciationSex{1} = time;

%Calculate the 2nd derivative to determine if selection is disruptive
disruptive = zeros(1,m);
for j=1:m
    secderF = 0;
    for i=1:n
        secderF = secderF - 1/2*(AMAX*Ea*TFoodSex{1}(i,end)*exp(-(TTraiSex{1}(j,end) - THETAF(1,i))^2/(2*TAU^2))*(- TTraiSex{1}(j,end)^2 + 2*TTraiSex{1}(j,end)*THETAF(1,i) + TAU^2 - THETAF(1,i)^2))/TAU^4;
    end
    secder = secderF;
    if secder>0
        disruptive(1,j)=1;
    else
        disruptive(1,j)=0;
    end
end


%--------------------------------------------------------------------
%SPECIATION
%If selection is disruptive split the population(s) with disruptive
%selection

timemaxrad=1E8;
timerad=0;
cont = 1;

while (sum(disruptive)>0 && m<n+15) && timerad<timemaxrad

    timerad = timerad + time;
    
    initialFood = TFoodSex{cont}(:,end);
    initialNFem = [];
    initialNMal = [];
    initialTraiSex = [];
    
    %Diversification

    for j=1:m
        if disruptive(1,j)==1
            initialNFem = [initialNFem; TNFemSex{cont}(j,end)/2; TNFemSex{cont}(j,end)/2];
            initialNMal = [initialNMal; TNMalSex{cont}(j,end)/2; TNMalSex{cont}(j,end)/2];
            initialTraiSex = [initialTraiSex; TTraiSex{cont}(j,end)*1-1E-3; TTraiSex{cont}(j,end)*1+1E-3];
        else
            initialNFem = [initialNFem; TNFemSex{cont}(j,end)];
            initialNMal = [initialNMal; TNMalSex{cont}(j,end)];
            initialTraiSex = [initialTraiSex; TTraiSex{cont}(j,end)];
        end
    end
    
    m=m+sum(disruptive);
    
    % Eco-evo dynamics
    
    tmax1 = 10000;
    
    FoodSex = zeros(n,tmax1);
    NFem = zeros(m,tmax1);
    NMal = zeros(m,tmax1);
    TraiSex = zeros(m,tmax1);
    
    FoodSex(:,1) = initialFood;
    NFem(:,1) = initialNFem;
    NMal(:,1) = initialNMal;
    TraiSex(:,1) = initialTraiSex;

    ARsexual = zeros(m,n);
    traitchange = 1;
    time=1;

    while (time<tmax1-1 && traitchange>1E-8)

        for j=1:m
        %Calculate attack rates on each resource
        ARsexual(j,1:end) = AMAX.*exp(-((TraiSex(j,time)-THETAF).^2)./(2*(TAU^2))); %attack rate of morph j on preexisting resources
        end    

        %find the ecological equilibrium
        tmaxe = 100;
        dif = 1;
        rep = 0;
        x0    = [FoodSex(:,time)' NFem(:,time)' NMal(:,time)'];
        while (dif>1E-6 && rep<20)
            fod = @(t,x) EcoEquiSexual(t,x,ARsexual,n,m,rho,FTmax,Ea,delta);
            [t1,x1] = ode45(fod,[0 tmaxe],x0);
            dif = abs(x1(end-1,1)-x1(end,1));
            rep = rep + 1;
            x0  = x1(end,:);
        end

        FoodSex(:,time+1)=x1(end,1:n);
        NFem(:,time+1)=x1(end,n+1:n+m);
        NMal(:,time+1)=x1(end,n+m+1:n+2*m);

        fod   = @(t,Tr)EvoDynSexual(t,Tr,n,m,FoodSex(:,time+1)',NFem(:,time+1)',NMal(:,time+1)',Ea,AMAX,THETAF,TAU,mu,sigma);
        [tevo,xevo] = ode23(fod,[0 1E6],TraiSex(:,time));
        TraiSex(:,time+1) = xevo(end,:);

        traitvar = zeros(m,1);
        for j=1:m
            traitvar(j,1) = abs(TraiSex(j,time+1)-TraiSex(j,time));
        end
        traitchange = max(traitvar);
        time = time+1;
    end

    TFoodSex{cont+1} = FoodSex(:,1:time);
    TNFemSex{cont+1} = NFem(:,1:time);
    TNMalSex{cont+1} = NMal(:,1:time);
    TTraiSex{cont+1} = TraiSex(:,1:time);
    TNPopSex{cont+1} = m*ones(time,1);
    tspeciationSex{cont+1} = time;

    %Calculate the 2nd derivative to determine if selection is disruptive
    disruptive = zeros(1,m);
    for j=1:m
        secderF = 0;
        for i=1:n
            secderF = secderF - 1/2*(AMAX*Ea*TFoodSex{cont+1}(i,end)*exp(-(TTraiSex{cont+1}(j,end) - THETAF(1,i))^2/(2*TAU^2))*(- TTraiSex{cont+1}(j,end)^2 + 2*TTraiSex{cont+1}(j,end)*THETAF(1,i) + TAU^2 - THETAF(1,i)^2))/TAU^4;
        end
        secder = secderF;
        if secder>0
            disruptive(1,j)=1;
        else
            disruptive(1,j)=0;
        end
    end
    cont = cont+1;
end

