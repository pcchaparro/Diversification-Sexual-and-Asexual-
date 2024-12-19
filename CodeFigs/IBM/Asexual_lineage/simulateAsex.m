function [Ffinal,distAsex,time,timeniches] = simulateAsex(tmax,tmaxafniches,Pop,n,rho,volHab,Fi,Fmaxi,THETAF,Ea,AMAX,TAU,delta,nEloci,mutrateEl,sigmamutEl,Eloci,ecotr,Biores,aFi,histedges)

    time=1;
    tocniches=0;
    timeniches=0;
    
    while ((tocniches==0 || time-tocniches<tmaxafniches) && time<tmax)
        time = time+1;
        
        %population stats
        sPop = size(Pop);
        sp=sPop(1,1);

        %Feeding
        ingest = zeros(sp,n);
        for j=1:n
            ingest(:,j) = Pop(:,aFi(1,j)).*Fi(1,j);
        end
        Pop(:,Biores)= Pop(:,Biores) + Ea.*(sum(ingest,2)); %energy available to use in biomass

        %update density of food resources
        for j=1:n
            Fi(1,j) = max(0,(Fi(1,j) + rho*(Fmaxi(1,j)-Fi(1,j)) - sum(ingest(:,j))/volHab));
        end
        
        %This simulation runs for tmaxafniches time points after all
        %food resources have been depleted below 0.1 of their carrying
        %capacity.
        deplet=(Fi./Fmaxi>.1); 
        if (sum(deplet)==0 && tocniches==0)
            tocniches=time; %save the time at which all niches get filled for the first time
        end
        
        
        deplet5=(Fi./Fmaxi>.5); 
        if (sum(deplet5)==0 && timeniches==0)
            timeniches=time; %Save the time at which food resources have been depleted below half of their carrying capacity
        end

        %Reproduction
        Bioresfem = Pop(:,Biores);
        Probrep   = Bioresfem; %Probability to reproduce is proportional to biomass reserves
        Probrep(Probrep>1) = 1;
        repfemale = binornd(1,Probrep); %Females that reproduce
        nrepfemale=sum(repfemale);
        idrepfem=find(repfemale); %id of reproducing females
        Pop(:,Biores)= 0; %reset the reproductive buffer

        if nrepfemale>0
            %Create offspring matrix
            offspring = zeros(nrepfemale,sPop(1,2));
            baby = 0;

            for j=1:nrepfemale
                baby = baby+1;
                offspring(baby,1:2*nEloci) = Pop(idrepfem(j,1),Eloci); %Baby inherits genes from mom
            end

            %mutate the offspring randomly
            mutating=(rand(baby,2*nEloci)<mutrateEl).*sigmamutEl.*randn(baby,2*nEloci);
            offspring(:,Eloci)=offspring(:,Eloci)+mutating;

            %complete the rest of the matrix
            offspring(:,ecotr) = sum(offspring(:,Eloci),2);
            offspring(:,Biores)= 0;
            for j=1:n
                offspring(:,aFi(1,j)) = AMAX*exp(-((offspring(:,ecotr)-THETAF(1,j)).^2)./(2*TAU^2));
            end

            %add the newborns to the population
            Pop = [Pop;offspring];
        end

        %remove death individuals
        sPop = size(Pop);
        sp=sPop(1,1);
        dead=rand(sp,1)<delta;
        Pop(dead,:)=[];
        
    end
    
    %save the last time point of food densities and the population
    %histogram for each lineage
    
    Ffinal   = Fi;
    distAsex  = histcounts(Pop(:,ecotr),histedges); %count ind with trait ecotr in the bins of histedges
end
