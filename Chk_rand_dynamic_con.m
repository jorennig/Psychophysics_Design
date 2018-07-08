%% Creates stimulus matrix; randomizes; checks randomizations; retrieves files
clear all
close all
clc

%% Set variables
reps = 5; % reps 
n_talker = 5; % total talkers
n_cond = 3; % # videos per talker
n_stim = n_talker * n_cond;

%% Define stimuli
stim_dir = [pwd '/Stimuli/'];

s = {'MS1_F06_C1.mp4'; 'MS1_F06_C2.mp4'; 'MS1_F06_C3.mp4'; 'MS1_M03_C1.mp4'; 'MS1_M03_C2.mp4';'MS1_M03_C3.mp4';'MS1_F08_C1.mp4';'MS1_F08_C2.mp4';'MS1_F08_C3.mp4';'MS1_M02_C1.mp4';'MS1_M02_C2.mp4';'MS1_M02_C3.mp4';'MS1_F99_C1.mp4';'MS1_F99_C2.mp4';'MS1_F99_C3.mp4'};

%% Create matrix: STIMULI TALKER SYLLABLE MCG/CONTROL TALKER
stimuli = zeros(n_stim, 3);

stimuli(:, 1) = 1:n_stim; % stimuli numbers
talker_id = repmat(1:n_talker, n_cond, 1);
stimuli(:,2) = reshape(talker_id, n_stim, 1); % id talker (#1-5)

syllables = repmat(1:n_cond, n_talker, 1);
stimuli(:,3) = reshape(syllables', n_stim, 1); % 1 = ba, 2 = da, 3 = ga

%% Double McG talker trials (reps_m = 10)
talk = unique(stimuli(:,2)); % number of talkers - vector
stim = 1:n_stim; % vector of stimuli numbers
syll = unique(stimuli(:,3)); % vector of syll 
check_tot = 0;

while check_tot == 0
    check_v_tot = zeros(numel(talk),1);
    stim_mat = [];

    for i = 1:reps % create randomized vector x #reps

       check = 0;
       while check == 0
           check_t = zeros(numel(talk),1);
           check_s = zeros(n_stim,1);

           shuffled_set = Shuffle(stimuli,2); % shuffle a new set/repetition of stimuli

           for j = 1:numel(talk)

               seq_t = ones(3,1)*talk(j);
               pat_t = findpattern(shuffled_set(:,2),seq_t);

                if isempty(pat_t)
                    check_t(j) = 0;
                else
                    check_t(j) = 1;
                end
           end
           
           for k = 1:n_stim
               seq_s = ones(2,1)*stim(k);
               pat_s = findpattern(shuffled_set(:,1),seq_s);
               
                 if isempty(pat_s)
                    check_s(k) = 0;
                else
                    check_s(k) = 1;
                 end 
           end
                          
            if sum(check_t) == 0 && sum(check_s) == 0 
                check = 1;
            end
       end

       stim_mat = [stim_mat; shuffled_set]; % add new set to stimulus matrix

    end

   % Check overal design 
   for j = 1:numel(talk)

       seq_t = ones(3,1)*talk(j);
       pat = findpattern(stim_mat(:,2),seq_t);

        if isempty(pat)
            check_v_tot(j) = 0;
        else
            check_v_tot(j) = 1;
        end
   end

    if sum(check_v_tot) == 0
        check_tot = 1;
    end
end

%% Retrieve file names
stim_file = cell(length(stim_mat), 1);

for j = 1:(length(stim_mat)) % loop through each line in stim_mat
    
    stim_num = stim_mat(j,1);  % get stim number from shuffled matrix
    stim_file{j} = s(stim_num); % get stim string
end