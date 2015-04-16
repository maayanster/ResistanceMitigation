function [generations2thresh, q_array] = Stochastic_Model(q_freq, Pref, K, ...
    WErr_ref, WErs_ref, WEss_ref, WErr_toxic, WErs_toxic, WEss_toxic,...
   gen_num)
% Stochastic_Model simulates insecticide resistance starting 
% ------------------------------------------------------
% [q_fix_percentileA, q_fix_percentileB, q_fix_percentileC, q_fix_mean]...
%  = stochastic_10March(q_freq, Pref, WErr_ref, WErs_ref, WEss_ref,...
%  WErr_toxic, WErs_toxic, WEss_toxic, gen_num, sim_num)
% ------------------------------------------------------
% Description: Model insecticie resistance in a finite population, 
%              panmictic mating, with selection from various sources,
%              and density dependent effects
% Input:       {q_freq} initial frequency of resistant alleles
%              {Pref} proportion of area that is refuge (no Bt)
%              {WErr_ref} fitness of RR in refuge with natural enemies
%              {WErs_ref} fitness of RS in refuge with natural enemies
%              {WEss_ref} fitness of SS in refuge with natural enemies
%              {WErr_toxic} Fitness of RR in field due to natural enemies
%              {WErs_toxic} Fitness of RS in field due to natural enemies
%              {WEss_toxic} Fitness of SS in field due to natural enemies
%              {gen_num} number of generation that model runs
% Output:      {generations2thresh} Generation when q_freq is larger than 0.1

% Adrian Semmelink
% Classification: Honours project
% Last revision date: 15-April-2015

%% Initialize
population = 290;                     % Initial pest population  
p_freq = 1 - q_freq;                  % Initial frequency of S allele
progeny = 69;                         % Number of progeny produced per female 
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

%% OUTPUT
p_array=[];
q_array=[];
pop_array = [];

%% CALCULATIONS
% Find intrinsic growth rate by running get_rmax function 
r_arr = get_rmax(q_freq, Pref,population, p_freq, progeny, Wxrr_ref,...
    Wxrs_ref, Wxss_ref, MutationR);
% Find the intrinsic growth rate by averaging the growth rate over a number
% of generations calculated with get_rmax
intrinsicR = mean(r_arr);

% Calculate number of generations for loop
gen_number = gen_num - 1;

% Run stochastic model
while i <= gen_number 
    i = i+1;                      % Count generations
    % Round population each loop to convert into whole number
    population = round(population);

    % Mutation rate applied
    q_freq = q_freq + poissrnd(p_freq*MutationR, 1);
    p_freq = p_freq;


    % Determine proportion of males and females by applying binomial 
    % distribution to choose number of females, assume remainder males 
    Fem = binornd(population,0.5,1);
    Male = population - Fem;

    % Hardy-weinberg ratios applied to find expected relative frequency  
    SS = p_freq^2;
    RS = 2*q_freq*p_freq;
    RR = q_freq^2;

    % If error statement to check that RR/RS/SS total 1
    if ((RR+RS+SS)-1) > 0.000001 && ((RR+RS+SS)-1) < -0.000001; 
        error('RR+RS+SS do not equal 1');
    end

    % Sample the number of individuals from each genotype in this generation
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

    % For each pairing type find number of progeny   
    RRxRR_progeny = round(RRxRR*progeny);
    RRxRS_progeny = round(RRxRS*progeny);
    RRxSS_progeny = round(RRxSS*progeny);
    SSxSS_progeny = round(SSxSS*progeny);
    SSxRS_progeny = round(SSxRS*progeny);
    RSxRS_progeny = round(RSxRS*progeny);

    %%% For each type of pairing assigning number of progeny for each
    % genotype (birth rate per pairing)
    % Find number of RR progeny produced by RRxRR pairing
    RRxRR_progeny_RR = RRxRR_progeny;

    % Find number of RR and RS progeny produced by RRxRS pairing
    RRxRS_progeny_RR = binornd(RRxRS_progeny, 0.5);
    RRxRS_progeny_RS = RRxRS_progeny - RRxRS_progeny_RR;

    % Find number of RS progeny produced by RRxSS pairing
    RRxSS_progeny_RS = RRxSS_progeny;

    % Find number of SS progeny produced by SSxSS pairing
    SSxSS_progeny_SS = SSxSS_progeny;

    % Find number of SS and RS progeny produced by SSxRS pairing
    SSxRS_progeny_RS = binornd(SSxRS_progeny, 0.5);
    SSxRS_progeny_SS = SSxRS_progeny - SSxRS_progeny_RS;

    % Find number of SS, RS, and RR progeny produced by RSxRS pairing
    RSxRS_progeny_SS = binornd(RSxRS_progeny, 0.25);
    RSxRS_progeny_RR = binornd(RSxRS_progeny, 0.25);
    RSxRS_progeny_RS = RSxRS_progeny - (RSxRS_progeny_SS+RSxRS_progeny_RR);

    % Find total amount of progeny for RR, SS, and RS
    tot_progeny_RR = RRxRR_progeny_RR + RRxRS_progeny_RR + RSxRS_progeny_RR;
    tot_progeny_RS = RRxRS_progeny_RS + RRxSS_progeny_RS + SSxRS_progeny_RS + RSxRS_progeny_RS;
    tot_progeny_SS = SSxSS_progeny_SS + SSxRS_progeny_SS + RSxRS_progeny_SS;

    % Dividing progeny across toxic and refuge zone 
    RRref = binornd(tot_progeny_RR, Pref);
    SSref = binornd(tot_progeny_SS, Pref);
    RSref = binornd(tot_progeny_RS, Pref);
    RRtoxic = tot_progeny_RR - RRref;
    SStoxic = tot_progeny_SS - SSref;
    RStoxic = tot_progeny_RS - RSref;  

    % Selection - fitness without natural enemies 
    RRref_x = RRref*Wxrr_ref;
    SSref_x = SSref*Wxss_ref;
    RSref_x = RSref*Wxrs_ref;
    RRbt_x = RRtoxic*Wxrr_toxic;
    SSbt_x = SStoxic*Wxss_toxic;
    RSbt_x = RStoxic*Wxrs_toxic;

    %Selection due to natural enemies 
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

    % Create new variable lambda 
    Lambda = intrinsicR - 1;
    
    % Calculate population if change in population only followed
    % logistic growth (include density dependent mortality from carrying
    % capacity and intrinsic growth rate)
    N_logistic = population + Lambda * population * (1 - population/K);

    % Recalculate new population including density dependent effects
    % leaving population below carrying capacity at the same proportions of
    % RR, RS, and SS genotypes
    RR = RR*(N_logistic/N_new);
    RS = RS*(N_logistic/N_new);
    SS = SS*(N_logistic/N_new);   

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
    generations2thresh = NaN;
end

% Example plot of how q_array changes over time to show where initial and
% spread stage end/begin
%figure
%x = generations2thresh - 1;
%plot(1:x,q_array)
%xlabel('Frequency of resistant allele (R)')
%ylabel('Generations')


