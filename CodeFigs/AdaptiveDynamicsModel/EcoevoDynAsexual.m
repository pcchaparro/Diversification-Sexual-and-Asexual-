% This rutine computes the eco-evo dynamics of an asexually reproducing
% lineage when an ancestral population colonizes a habitat with n food resources


%-----ECO-EVO DYNAMICS ASEXUAL POPULATION-----

%Initial conditions
THETAF=DTHETA*(1:n); %Optimal trait values for n resources
m=1;
initialFood = FTmax/n*ones(n,1);
initialNind = 0.1;
intialTraitAsexual = 0.6;

%Initialize matrices to save info simulation

tmax1 = 1000;
FoodAsex = zeros(n,tmax1);
Nind = zeros(m,tmax1);
TraiAsex = zeros(m,tmax1);

FoodAsex(:,1) = initialFood;
Nind(:,1) = initialNind;
TraiAsex(:,1) = intialTraitAsexual;

ARasexual = zeros(m,n);
traitchange = 1;
time=1;

while (time<tmax1-1 && traitchange>1E-8)
    
    for j=1:m
    %Calculate attack rates on each resource
    ARasexual(j,:) = AMAX.*exp(-((TraiAsex(j,time)-THETAF).^2)./(2*(TAU^2))); %attack rate of morph j on preexisting resources
    end    
    
    %find the ecological equilibrium
    tmaxe = 100;
    dif = 1;
    rep = 0;
    x0    = [FoodAsex(:,time)' Nind(:,time)'];
    while (dif>1E-6 && rep<20)
        fod = @(t,x) EcoEquiAsexual(t,x,ARasexual,n,m,rho,FTmax,Ea,delta);
        [t1,x1] = ode45(fod,[0 tmaxe],x0);
        dif = abs(x1(end-1,1)-x1(end,1));
        rep = rep + 1;
        x0  = x1(end,:);
    end

    FoodAsex(:,time+1)=x1(end,1:n);
    Nind(:,time+1)=x1(end,n+1:n+m);
    
    fod   = @(t,Tr)EvoDynAsexual(t,Tr,n,m,FoodAsex(:,time+1)',Nind(:,time+1)',Ea,AMAX,THETAF,TAU,mu,sigma);
    [tevo,xevo] = ode23(fod,[0 1E6],TraiAsex(:,time));
    TraiAsex(:,time+1) = xevo(end,:);
    
    traitvar = zeros(m,1);
    for j=1:m
        traitvar(j,1) = abs(TraiAsex(j,time+1)-TraiAsex(j,time));
    end
    traitchange = max(traitvar);
    time = time+1;
end

TFoodAsex{1} = FoodAsex(:,1:time);
TNindAsex{1} = Nind(:,1:time);
TTraiAsex{1} = TraiAsex(:,1:time);
TNPopAsex{1} = m*ones(time,1);
tspeciationAsex{1} = time;

%Calculate the 2nd derivative to determine if selection is disruptive
disruptive = zeros(1,m);
for j=1:m
    secderF = 0;
    for i=1:n
        secderF = secderF -(AMAX*Ea*TFoodAsex{1}(i,end)*exp(-(TTraiAsex{1}(j,end) - THETAF(1,i))^2/(2*TAU^2))*(- TTraiAsex{1}(j,end)^2 + 2*TTraiAsex{1}(j,end)*THETAF(1,i) + TAU^2 - THETAF(1,i)^2))/TAU^4;
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
    
    initialFood = TFoodAsex{cont}(:,end);
    initialNind = [];
    initialTraiAsex = [];
    
    %Diversification

    for j=1:m
        if disruptive(1,j)==1
            initialNind = [initialNind; TNindAsex{cont}(j,end)/2; TNindAsex{cont}(j,end)/2];
            initialTraiAsex = [initialTraiAsex; TTraiAsex{cont}(j,end)*1-1E-3; TTraiAsex{cont}(j,end)*1+1E-3];
        else
            initialNind = [initialNind; TNindAsex{cont}(j,end)];
            initialTraiAsex = [initialTraiAsex; TTraiAsex{cont}(j,end)];
        end
    end
    
    m=m+sum(disruptive);
    
    % Eco-evo dynamics
    
    tmax1 = 10000;
    
    FoodAsex = zeros(n,tmax1);
    Nind = zeros(m,tmax1);
    TraiAsex = zeros(m,tmax1);
    
    FoodAsex(:,1) = initialFood;
    Nind(:,1) = initialNind;
    TraiAsex(:,1) = initialTraiAsex;

    ARasexual = zeros(m,n);
    traitchange = 1;
    time=1;
    
    while (time<tmax1-1 && traitchange>1E-8)
    
        for j=1:m
        %Calculate attack rates on each resource
        ARasexual(j,1:end) = AMAX.*exp(-((TraiAsex(j,time)-THETAF).^2)./(2*(TAU^2))); %attack rate of morph j on preexisting resources
        end    

        %find the ecological equilibrium
        tmaxe = 100;
        dif = 1;
        rep = 0;
        x0    = [FoodAsex(:,time)' Nind(:,time)'];
        while (dif>1E-6 && rep<20)
            fod = @(t,x) EcoEquiAsexual(t,x,ARasexual,n,m,rho,FTmax,Ea,delta);
            [t1,x1] = ode45(fod,[0 tmaxe],x0);
            dif = abs(x1(end-1,1)-x1(end,1));
            rep = rep + 1;
            x0  = x1(end,:);
        end

        FoodAsex(:,time+1)=x1(end,1:n);
        Nind(:,time+1)=x1(end,n+1:n+m);

        fod   = @(t,Tr)EvoDynAsexual(t,Tr,n,m,FoodAsex(:,time+1)',Nind(:,time+1)',Ea,AMAX,THETAF,TAU,mu,sigma);
        [tevo,xevo] = ode23(fod,[0 1E6],TraiAsex(:,time));
        TraiAsex(:,time+1) = xevo(end,:);

        traitvar = zeros(m,1);
        for j=1:m
            traitvar(j,1) = abs(TraiAsex(j,time+1)-TraiAsex(j,time));
        end
        traitchange = max(traitvar);
        time = time+1;
    end

    TFoodAsex{cont+1} = FoodAsex(:,1:time);
    TNindAsex{cont+1} = Nind(:,1:time);
    TTraiAsex{cont+1} = TraiAsex(:,1:time);
    TNPopAsex{cont+1} = m*ones(time,1);
    tspeciationAsex{cont+1} = time;

    %Calculate the 2nd derivative to determine if selection is disruptive
    disruptive = zeros(1,m);
    for j=1:m
        secderF = 0;
        for i=1:n
            secderF = secderF -(AMAX*Ea*TFoodAsex{cont+1}(i,end)*exp(-(TTraiAsex{cont+1}(j,end) - THETAF(1,i))^2/(2*TAU^2))*(- TTraiAsex{cont+1}(j,end)^2 + 2*TTraiAsex{cont+1}(j,end)*THETAF(1,i) + TAU^2 - THETAF(1,i)^2))/TAU^4;
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