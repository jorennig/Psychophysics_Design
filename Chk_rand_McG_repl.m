%% Creates stimulus matrix; randomizes; checks randomizations; retrieves files
clear all
close all
clc

%% Define variables & stimuli.
stim_dir = [pwd '/Stimuli_2/'];

s1 = 'Aud_F01_M1.mov';  % talker A, McG
s2 = 'Aud_F01_C1.mov';  % talker A, ba
s3 = 'Aud_F01_C3.mov';  % talker A, ga
s4 = 'MS1_F05_C2.mov';  % talker B, da

n_stim = 4;

%% Create matrix: STIM_NUM, TALKER,  SYLLABLE
stimuli = zeros(n_stim, 3); % empty matrix - # stim x # info columns
stimuli(:, 1) = 1:n_stim; % stimuli numbers

stimuli(:,2) = [1 1 1 2]; % id talker 
stimuli(:,3) = [1:4]; % 1 = mcg; 2 = ba; 3 = ga; 4 = da

reps_m = 20;
reps_c = 10;

%% Repetitions 
reps = [reps_m reps_c reps_c reps_c];

f = @(k) repmat(stimuli(k,:), round(reps(k)), 1);
stim_mat = cell2mat(arrayfun(f, (1:length(reps))', 'UniformOutput', false));

%% Randomize
% x dyn_reps to get to all trials
check_tot = 0;
n_stim = numel(unique(stimuli(:,1))); % vector of stimuli numbers
n_syll = numel(unique(stimuli(:,3))); % vector of syll 

check = 0;

while check == 0
    
    shuffled_stim = Shuffle(stim_mat,2);
    
    % AVOID same stimuli threepeat
     for i = 1:n_stim
        seq_s = [i i i];
        pat_s = findpattern(shuffled_stim(:,1),seq_s);
        
        if isempty(pat_s)
            check_i(i) = 0;
        else
            check_i(i) = 1;
        end
    end
    
    % AVOID same syllable threepeat
    for s = 1:n_syll 
        seq_s = [s s s];
        pat_s = findpattern(shuffled_stim(:,3),seq_s);
        
        if isempty(pat_s)
            check_s(s) = 0;
        else
            check_s(s) = 1;
        end
    end
    
    % AVOID speaker #2 threepeat
    seq_t = [2 2 2];
    pat_t = findpattern(shuffled_stim(:,2),seq_t);
    
    if isempty(pat_t)
        check_t = 0;
    else
        check_t = 1;
    end

    if sum(check_i) == 0 && sum(check_s) == 0 && check_t ==0
        check = 1;
    end

end

%% Retrieve file names
stim_file = cell(length(shuffled_stim), 1);

for j = 1:(length(shuffled_stim)) % loop through each line in stim_mat
    
    stim_num(j) = shuffled_stim(j,1);  % get stim number from shuffled matrix
    str = ['s' num2str(stim_num(j))]; % get stim string
    stim = genvarname(str); % str -> variable
    stim_file{j} = [eval(stim)]; % get file name 
end