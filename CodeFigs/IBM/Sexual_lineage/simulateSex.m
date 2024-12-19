function [Ffinal,distSex,time,timeniches] = simulateSex(tmax,tmaxafniches,Pop,n,rho,volHab,Fi,Fmaxi,THETAF,Ea,AMAX,TAU,delta,nEloci,nAloci,mutrateEl,mutrateAl,sigmamutEl,sigmamutAl,Eloci,Aloci,sex,ecotr,mattr,dadcont,Biores,aFi,histedges)

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
        Bioresfem = (Pop(:,sex)>.5).*Pop(:,Biores);
        Probrep   = Bioresfem; %Probability to reproduce is proportional to biomass reserves
        Probrep(Probrep>1) = 1;
        repfemale = binornd(1,Probrep); %Females that reproduce
        nrepfemale=sum(repfemale);
        idrepfem=find(repfemale); %id of reproducing females
        idmales = find(Pop(:,sex)<.5);
        Pop(:,Biores)= 0;      %reset the reproductive buffer

        if nrepfemale>0
            %Create offspring matrix
            offspring = zeros(nrepfemale,sPop(1,2));
            baby = 0;

            %find a dad for offspring

            ecotrmales = unique(Pop(idmales,ecotr));
            for j=1:nrepfemale
                mattrait = Pop(idrepfem(j,1),mattr);
                if mattrait>0
                    sigmaA = 1/(20*mattrait^2); %from Dieckman&Doebeli,1999
                    probdads = normpdf(ecotrmales,Pop(idrepfem(j,1),ecotr),sigmaA);
                    ecotrsel = finddad(ecotrmales,probdads); %selected eco trait
                    malesecotrsel = find((Pop(:,ecotr)>ecotrsel-1E-6 & Pop(:,ecotr)<ecotrsel+1E-6) & Pop(:,sex)<.5);
                    dad = randsample(malesecotrsel,1);
                elseif mattrait<0
                    sigmaD = 1/(mattrait^2); %from Dieckman&Doebeli,1999
                    probdads = 1-normpdf(ecotrmales,Pop(idrepfem(j,1),ecotr),sigmaD);
                    ecotrsel = finddad(ecotrmales,probdads); %selected eco trait
                    malesecotrsel = find((Pop(:,ecotr)>ecotrsel-1E-6 & Pop(:,ecotr)<ecotrsel+1E-6) & Pop(:,sex)<.5);
                    dad = randsample(malesecotrsel,1);
                else
                    dad = randsample(idmales,1);
                end

                locifromdad = binornd(1,dadcont*ones(1,nEloci+nAloci)); % loci that will receive one allele from dad
                selaldad    = locifromdad.*(unidrnd(2,1,nEloci+nAloci)+(0:2:(2*(nEloci+nAloci)-2))); %select alleles from dad
                alleledadrep= false(1,2*nEloci+2*nAloci);
                alleledadrep(selaldad(selaldad>0))=1; %These alleles from dad will replace the alleles from mom
                genfrommom  = [Pop(idrepfem(j,1),Eloci) Pop(idrepfem(j,1),Aloci)]; %The genome of the mom
                genfrommom(selaldad(selaldad>0))=0; %Remove the alleles that will be replaced by the alleles of the dad
                genfromdad = [Pop(dad,Eloci) Pop(dad,Aloci)].*alleledadrep; %The genome of the dad
                genoff = genfrommom+genfromdad;
                baby = baby+1;
                offspring(baby,1:2*nEloci+2*nAloci) = genoff;
            end

            %mutate the offspring randomly
            mutatingAl=(rand(baby,2*nAloci)<mutrateAl).*sigmamutAl.*randn(baby,2*nAloci);
            mutatedAl=offspring(:,Aloci)+mutatingAl;
            mutatedAl(mutatedAl>1)=1;
            mutatedAl(mutatedAl<-1)=-1;
            offspring(:,Aloci)=mutatedAl;
            mutatingEl=(rand(baby,2*nEloci)<mutrateEl).*sigmamutEl.*randn(baby,2*nEloci);
            offspring(:,Eloci)=offspring(:,Eloci)+mutatingEl;

            %complete the rest of the matrix
            offspring(:,sex)   = randi(2,nrepfemale,1)-1;
            offspring(:,ecotr) = sum(offspring(:,Eloci),2);
            offspring(:,mattr) = mean(offspring(:,Aloci),2);
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
    
    %save the last time point of food densities and the population histogram
    
    Ffinal   = Fi;
    distSex  = histcounts(Pop(:,ecotr),histedges); %count ind with trait ecotr in the bins of histedges
end
