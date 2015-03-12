q_freq_arr = [0.001, 0.001, 0.001, 0.001];
Refuge_arr = [0.3, 0.3, 0.3, 0.3];
% WE...{NE absent, NE present - no differential predation, NE present differential predation (low), NE present - differential predation (high)} 
WErr_ref_arr = [1, 0.5, 0.4, 0.2]; 
WErs_ref_arr = [1, 0.5, 0.5, 0.5]; 
WEss_ref_arr = [1, 0.5, 0.5, 0.5]; 

gen_num_arr = [100, 100, 100, 100];

q_s = zeros(1,length(q_freq_arr));
for nn=1:4  %number of trials
    [q_freq] = Deterministic_20Feb(q_freq_arr(nn),Refuge_arr(nn), WErr_ref_arr(nn), WErs_ref_arr(nn), WEss_ref_arr(nn), gen_num_arr(nn));
    q_s(nn) = q_freq;
    
end

%resistant_not_extinct = sum(q_s>0.1); % how many trials did the R allele not go to 0
%display(resistant_not_extinct)

%resistant_fixation = sum(q_s>0.1); % how many trial did the R allele get to 0.1 of pop.
%display(resistant_fixation)