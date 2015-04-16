function r_arr = get_rmax(q_freq, Pref,population, p_freq, progeny,...
    Wxrr_ref,Wxrs_ref, Wxss_ref, MutationR)
% get_rmax calculates the growth rate of pest population over a number of
% generations
% ------------------------------------------------------------------------
% r_arr = get_rmax(q_freq, Pref, WErr_ref, WErs_ref, WEss_ref,...
%         population, p_freq, mean_no_progeny, std_error_progeny, ...
%         Wxrr_toxic,Wxss_toxic, Wxrs_toxic, WErr_toxic,WErs_toxic, ...
%         WEss_toxic, Wxrr_ref,Wxrs_ref, Wxss_ref, MutationR)
% ------------------------------------------------------------------------
% Description: calculates the growth rate of pest population over a ...
%              number of generations which can then be used to find ...
%              intrinsic growth rate
% Input:  {q_freq} (see other functions for details) 
%         {Pref}
%         {WErr_ref}
%         {WErs_ref}
%         {WEss_ref}
%         {population}
%         {p_freq}
%         {mean_no_progeny}
%         {std_error_progeny}
%         {Wxrr_toxic}
%         {Wxss_toxic}
%         {Wxrs_toxic}
%         {WErr_toxic}
%         {WErs_toxic}
%         {WEss_toxic}
%         {Wxrr_ref}
%         {Wxrs_ref}
%         {Wxss_ref}
%         {MutationR}
% Output: {r_arr}     an array of growth rates

% Adrian Semmelink
% Classification: Honours project
% Last revision date: 25-March-2015
%% INITIALIZE
r_arr = [];
kk=0;

%% CALCULATIONS
while kk <= 4
    kk = kk+1;                      % Count generations
    % Round population each loop to convert into whole number 
    population = round(population);

    % Mutation rate applied
    q_freq = q_freq + p_freq*MutationR;
    p_freq = p_freq + q_freq*MutationR;

    % Determine how many males/females - use binomial - assume replacement 
    % distribution to choose a number of females, assume remainder males)  
    Fem = binornd(population,0.5,1);
    Male = population - Fem;

    % Hardy-weinberg ratios are used to calculate expected relative frequency  
    SS = p_freq^2;
    RS = 2*q_freq*p_freq;
    RR = q_freq^2;

    %%% if error statement to check that RR/RS/SS total 1
    if ((RR+RS+SS)-1) > 0.000001 && ((RR+RS+SS)-1) < -0.000001; 
        error('RR+RS+SS do not equal 1');
    end

    % Sample the number of individuals from each genotype in this generation
    % using poisson distribution  - sample "small" genotypes first

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
    RRxRR_progeny = RRxRR*progeny;
    RRxRS_progeny = RRxRS*progeny;
    RRxSS_progeny = RRxSS*progeny;
    SSxSS_progeny = SSxSS*progeny;
    SSxRS_progeny = SSxRS*progeny;
    RSxRS_progeny = RSxRS*progeny;

    % For each type of pairing assigning number of progeny for each
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
    tot_progeny_RR = RRxRR_progeny_RR + RRxRS_progeny_RR + ...
        RSxRS_progeny_RR;
    tot_progeny_RS = RRxRS_progeny_RS + RRxSS_progeny_RS + ...
        SSxRS_progeny_RS + RSxRS_progeny_RS;
    tot_progeny_SS = SSxSS_progeny_SS + SSxRS_progeny_SS + ...
        RSxRS_progeny_SS;

    % Dividing progeny across treatment and refuge 
    RRref = binornd(tot_progeny_RR, Pref);
    SSref = binornd(tot_progeny_SS, Pref);
    RSref = binornd(tot_progeny_RS, Pref);
    RRtoxic = tot_progeny_RR - RRref;
    SStoxic = tot_progeny_SS - SSref;
    RStoxic = tot_progeny_RS - RSref;  

    % Selection - fitness w/o natural enemies - all zones are assumed to
    % not be toxic as we are calculating optimum growth rate 
    RRref = RRref*Wxrr_ref;
    SSref = SSref*Wxss_ref;
    RSref = RSref*Wxrs_ref;
    RRtoxic = RRtoxic*Wxrr_ref;
    SStoxic = SStoxic*Wxss_ref;
    RStoxic = RStoxic*Wxrs_ref;

    % Calculate total # of RR, RS, and SS in both field types (refuge
    % and toxic)
    RR = RRtoxic + RRref;
    RS = RStoxic + RSref;
    SS = SStoxic + SSref;

    % Calculate new population without density dependent effects
    new_population = RR + RS + SS;
    r = new_population/population;
    r_arr = [r, r_arr];
    
    population = new_population;
end


