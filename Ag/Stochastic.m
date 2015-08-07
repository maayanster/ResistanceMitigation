function [generations2thresh, q_array] = Stochastic(q_freq, Pref, K, ...
    WErr_ref, WErs_ref, WEss_ref, WErr_toxic, WErs_toxic, WEss_toxic,...
   gen_num)
% Stochastic simulates insecticide resistance starting 
% ------------------------------------------------------
% [generations2thresh, q_array] = Stochastic(q_freq, Pref, K, WErr_ref,...
%  WErs_ref, WEss_ref, WErr_toxic,WErs_toxic, WEss_toxic,gen_num)
% ------------------------------------------------------
% Description: Model insecticie resistance in a finite population, 
%              panmictic mating, with selection from various sources,
%              and density dependent effects
% Input:       {q_freq} initial frequency of resistant alleles
%              {Pref} proportion of area that is refuge (no Bt)
%              {K} carrying capacity
%              {WErr_ref} fitness of RR in refuge with natural enemies
%              {WErs_ref} fitness of RS in refuge with natural enemies
%              {WEss_ref} fitness of SS in refuge with natural enemies
%              {WErr_toxic} Fitness of RR in field due to natural enemies
%              {WErs_toxic} Fitness of RS in field due to natural enemies
%              {WEss_toxic} Fitness of SS in field due to natural enemies
%              {gen_num} number of generation that model runs
% Output:      {generations2thresh} Generation when q_freq is larger than 0.1
%              {q_array} array of resistant allele frequencies over number
%              of generations

% Adrian Semmelink
% Classification: Insecticide resistance project
% Last revision date: 21-July-2015

%% Initialize
population = 2900;                     % Initial pest population  
p_freq = 1 - q_freq;                  % Initial frequency of S allele
progeny = 236;                         % Number of progeny produced per female 
Wxrr_toxic = 0.207;                   % Fitness of RR in field w/o natural enemies
Wxss_toxic = 0;                       % Fitness of SS in field w/o natural enemies
Wxrs_toxic = 0;                       % Fitness of RS in field w/o natural enemies
Wxrr_ref = 0.207;                     % Fitness of RR in refuge w/o natural enemies
Wxrs_ref = 0.207;                     % Fitness of RS in refuge w/o natural enemies
Wxss_ref = 0.208;                     % Fitness of SS in refuge w/o natural enemies
MutationR = 0.00005;                  % Mutation rate of R to S or S to R (Sisterson, 2004)
i = 0;                                % Initialize generation count
q_threshold= 0.1;                     % q frequency at which program stops 
generations2thresh = 0;               % Initialize number of generations to threshold
intrinsicR = 0.15;                     % Intrinsic growth rate (Gryspeit, 2012)

%% OUTPUT
p_array=[];
q_array=[];
pop_array = [];

%% CALCULATIONS
% Calculate number of generations for loop
gen_number = gen_num - 1;

% Run stochastic model
while i <= gen_number 
    i = i+1;                      % Count generations
    %%% STEP 1: Calculate proportions of adults by genotype and sex
    % Round population each loop to convert into whole number
    population = round(population);

    % Mutation rate applied
    q_freq = q_freq + (p_freq*MutationR);
    p_freq = p_freq + (q_freq*MutationR);

    % Determine proportion of males and females by applying binomial 
    % distribution to choose number of females, assume remainder males 
    Fem = binornd(population,0.5,1);
    Male = population - Fem;
    
    % Hardy-weinberg ratios applied to find expected relative frequency  
    SS = p_freq^2;
    RS = 2*(q_freq*p_freq);
    RR = q_freq^2;

    % If error statement to check that RR/RS/SS total 1
    if ((RR+RS+SS)-1) > 0.000001 && ((RR+RS+SS)-1) < -0.000001; 
        error('RR+RS+SS do not equal 1');
    end

    %{
 Sample the number of individuals from each genotype in this generation
    % using poisson distribution - sampling "small" genotypes first
    MaleRR = poissrnd(RR*Male, 1);    
    FemRR = poissrnd(RR*Fem, 1);
    MaleRS = poissrnd(RS*Male, 1);
    FemRS = poissrnd(RS*Fem, 1);

    % Calculate remainders for SS genotypes
    MaleSS = Male - (MaleRS + MaleRR);   
    if (MaleRS + MaleRR) > Male;
        MaleSS = 0;
    end

    FemSS = Fem - (FemRS + FemRR);
    if (FemRS + FemRR) > Fem
        FemSS = 0;
    end
    %}

    MaleRR = round(RR*Male);    
    FemRR = round(RR*Fem);
    MaleRS = round(RS*Male);
    FemRS = round(RS*Fem);
    MaleSS = round(SS*Male);
    FemSS = round(SS*Fem);

    %%% Calculate pairings between genotypes
    % Tracking Female RR genotype    
    FemRRxMaleRR = poissrnd(FemRR*MaleRR/Male, 1);   
    FemRRxMaleRS = poissrnd(FemRR*MaleRS/Male, 1);
    FemRRxMaleSS = FemRR - (FemRRxMaleRR + FemRRxMaleRS);   
    if FemRRxMaleSS < 0 
        FemRRxMaleSS = 0;
        FemRRxMaleRS = FemRRxMaleRS/(FemRRxMaleRR+FemRRxMaleRS)*FemRR;
        FemRRxMaleRR = FemRRxMaleRR/(FemRRxMaleRR+FemRRxMaleRS)*FemRR;
    end

    % Tracking Female RS genotype
    FemRSxMaleRR = poissrnd(FemRS*MaleRR/Male, 1);
    FemRSxMaleRS = poissrnd(FemRS*MaleRS/Male, 1);
    FemRSxMaleSS = FemRS - (FemRSxMaleRR + FemRSxMaleRS);
    if FemRSxMaleSS < 0
        FemRSxMaleSS = 0;
        FemRSxMaleRR = FemRSxMaleRR/(FemRSxMaleRR + FemRSxMaleRS)*FemRS;
        FemRSxMaleRS = FemRSxMaleRS/(FemRSxMaleRR + FemRSxMaleRS)*FemRS;
    end

    % Tracking Female SS genotype
    FemSSxMaleRR = poissrnd(FemSS*MaleRR/Male, 1);
    FemSSxMaleRS = poissrnd(FemSS*MaleRS/Male, 1);
    FemSSxMaleSS = FemSS - (FemSSxMaleRR + FemSSxMaleRS);
    if FemSSxMaleSS < 0;
        FemSSxMaleRR = FemSSxMaleRR/(FemSSxMaleRR+FemSSxMaleRS)*FemSS;
        FemSSxMaleRS = FemSSxMaleRS/(FemSSxMaleRR+FemSSxMaleRS)*FemSS;
        FemSSxMaleSS = 0;
    end

    % Calculate 6 different types of pairings:
    RRxRR = FemRRxMaleRR; 
    RRxRS = FemRRxMaleRS + FemRSxMaleRR;
    RRxSS = FemRRxMaleSS + FemSSxMaleRR;
    SSxSS = FemSSxMaleSS;
    SSxRS = FemSSxMaleRS + FemRSxMaleSS;
    RSxRS = FemRSxMaleRS;

    %%% STEP 2: Find number of progeny produced by each genotype
    %{For each pairing type find number of progeny   
   % RRxRR_progeny = round(RRxRR*progeny);
  %  RRxRS_progeny = round(RRxRS*progeny);
    %RRxSS_progeny = round(RRxSS*progeny);
    %SSxSS_progeny = round(SSxSS*progeny);
   % SSxRS_progeny = round(SSxRS*progeny);
  %  RSxRS_progeny = round(RSxRS*progeny);
  

    %For each type of pairing assigning number of progeny for each
    % genotype (birth rate per pairing)
   
    % Find number of RR progeny produced by RRxRR pairing
    RRxRR_progeny_RR = round(RRxRR*progeny);
    
    % Find number of RS progeny produced by RRxSS pairing
    RRxSS_progeny_RS = round(RRxSS*progeny);

    % Find number of SS progeny produced by SSxSS pairing
    SSxSS_progeny_SS = round(SSxSS*progeny);

    % Find number of RR and RS progeny produced by RRxRS pairing
    RRxRS_progeny_RR = 0;
    RRxRS_progeny_RS = 0;
    for jj=1:RRxRS
        progeny_RR = round(binornd(progeny, 0.5));
        progeny_RS = progeny - progeny_RR;
        
        RRxRS_progeny_RR = RRxRS_progeny_RR + progeny_RR;
        RRxRS_progeny_RS = RRxRS_progeny_RS + progeny_RS;
    end
    
    % Find number of SS and RS progeny produced by SSxRS pairing
    SSxRS_progeny_SS = 0;
    SSxRS_progeny_RS = 0;
    for jj=1:SSxRS
        progeny_SS = round(binornd(progeny, 0.5));
        progeny_RS = progeny - progeny_SS;
        
        SSxRS_progeny_SS = SSxRS_progeny_SS + progeny_SS;
        SSxRS_progeny_RS = SSxRS_progeny_RS + progeny_RS;
    end
    
    % Find number of SS, RS, and RR progeny produced by RSxRS pairing
    RSxRS_progeny_RR = 0;
    RSxRS_progeny_SS = 0;
    RSxRS_progeny_RS = 0;
    
    for jj=1:RSxRS
        progeny_RR = round(binornd(progeny, 0.25));
        progeny_SS = round(binornd(progeny, 0.25));
        progeny_RS = progeny - progeny_SS - progeny_RR;
        
        RSxRS_progeny_RR = RSxRS_progeny_RR + progeny_RR;
        RSxRS_progeny_SS = RSxRS_progeny_SS + progeny_SS;
        RSxRS_progeny_RS = RSxRS_progeny_RS + progeny_RS;
    end

    % Find total amount of progeny for RR, SS, and RS
    tot_progeny_RR = RRxRR_progeny_RR + RRxRS_progeny_RR + RSxRS_progeny_RR;
    tot_progeny_RS = RRxRS_progeny_RS + RRxSS_progeny_RS + SSxRS_progeny_RS + RSxRS_progeny_RS;
    tot_progeny_SS = SSxSS_progeny_SS + SSxRS_progeny_SS + RSxRS_progeny_SS;

    %%% STEP 3: Divide progeny by zone
    % Dividing progeny across toxic and refuge zone 
    %RRref = binornd(tot_progeny_RR, Pref);
    %SSref = binornd(tot_progeny_SS, Pref);
    %RSref = binornd(tot_progeny_RS, Pref);
    %RRtoxic = tot_progeny_RR - RRref;
    %SStoxic = tot_progeny_SS - SSref;
    %RStoxic = tot_progeny_RS - RSref;  

    RRref = tot_progeny_RR*Pref;
    SSref = tot_progeny_SS*Pref;
    RSref = tot_progeny_RS*Pref;
    RRtoxic = tot_progeny_RR - RRref;
    SStoxic = tot_progeny_SS - SSref;
    RStoxic = tot_progeny_RS - RSref; 
    
    %%% STEP 4: Selection by zones
    % Selection not due to natural enemies 
    %RRref_x = binornd(RRref, Wxrr_ref);
    %RSref_x = binornd(RSref, Wxrs_ref);
    %SSref_x = binornd(SSref, Wxss_ref);
    %RRbt_x = binornd(RRtoxic, Wxrr_toxic);
    %SSbt_x = binornd(SStoxic, Wxss_toxic);
    %RSbt_x = binornd(RStoxic,Wxrs_toxic);

    % Selection due to natural enemies 
    %RRref = binornd(RRref_x,WErr_ref);
    %SSref = binornd(SSref_x,WEss_ref);
    %RSref = binornd(RSref_x,WErs_ref);
    %RRtoxic = binornd(RRbt_x,WErr_toxic);
    %SStoxic = binornd(SSbt_x,WEss_toxic);
    %RStoxic = binornd(RSbt_x,WErs_toxic);
    
    % Selection - fitness not due to natural enemies 
    RRref_x = RRref* Wxrr_ref;
    RSref_x = RSref*Wxrs_ref;
    SSref_x = SSref*Wxss_ref;
    RRbt_x = RRtoxic*Wxrr_toxic;
    SSbt_x = SStoxic*Wxss_toxic;
    RSbt_x = RStoxic*Wxrs_toxic;

    % Selection due to natural enemies 
    RRref = RRref_x*WErr_ref;
    SSref = SSref_x*WEss_ref;
    RSref = RSref_x*WErs_ref;
    RRtoxic = RRbt_x*WErr_toxic;
    SStoxic = SSbt_x*WEss_toxic;
    RStoxic = RSbt_x*WErs_toxic;

    % Calculate total number of RR, RS, and SS in both field types
    RR = RRtoxic + RRref;
    RS = RStoxic + RSref;
    SS = SStoxic + SSref;

    % Calculate new population 
    N_new = RR + RS + SS;
    
    % Calculate population if change in population only followed
    % logistic growth (include density dependent mortality from carrying
    % capacity and intrinsic growth rate)
    N_logistic = population + intrinsicR * population * (1 - population/K);

    % Recalculate new population including density dependent effects
    % leaving population below carrying capacity at the same proportions of
    % RR, RS, and SS genotypes
    RR = RR*(N_logistic/N_new);
    RS = RS*(N_logistic/N_new);
    SS = SS*(N_logistic/N_new); 
    
    % Calculate total population
    population = RR+RS+SS;
    
    %%% STEP 5: Calculate frequency of resistant and suceptible alleles
    % Find number of each allele after all selection stages
    q_sum = 2*RR + RS;
    p_sum = 2*SS + RS;

    % Finding p and q
    q_freq = (q_sum/(q_sum+p_sum));
    p_freq = (p_sum/(p_sum+q_sum));

    % Displaying generation at which fixation occurs
    if  q_freq >= q_threshold && generations2thresh ==0;
        generations2thresh = i;
        break
    end

    % Saving our population and fequencies info for plots (not used in 
    % all runs)
    q_array = [q_array, q_freq];
end

% If no resistance developed by the end of the run the number of
% generations to threshold is recorded as the number of generations
% simulated
if q_freq <= q_threshold;
    generations2thresh = gen_num;
end


