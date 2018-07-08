%% Creates stimulus matrix; randomizes; checks randomizations; retrieves files
clear all
close all
clc

%% Set variables
reps = 5; % reps for control
reps_m = 10; % reps for mcgurk
mcg_talkers = 4; 
c_talkers = 1; % control talkers 
n_talker = mcg_talkers + c_talkers; % total talkers
n_cond = 3; % # videos per talker
n_stim = n_talker * n_cond;

%% Define stimuli
stim_dir = [pwd '/Stimuli/'];

s1 = 'MS1_F06_C1.mp4'; % Talker A - ba
s2 = 'MS1_F06_C3.mp4'; % Talker A - ga
s3 = 'MS1_F06_M1.mp4'; % Talker A - mcg
s4 = 'MS1_M03_C1.mp4'; % Talker B - ba
s5 = 'MS1_M03_C3.mp4'; % Talker B - ga
s6 = 'MS1_M03_M1.mp4'; % Talker B - mcg
s7 = 'MS1_F08_C1.mp4'; % Talker C - ba
s8 = 'MS1_F08_C3.mp4'; % Talker C - ga
s9 = 'MS1_F08_M1.mp4'; % Talker C - mcg
s10 = 'MS1_M02_C1.mp4'; % Talker D - ba
s11 = 'MS1_M02_C3.mp4'; % Talker D - ga
s12 = 'MS1_M02_M1.mp4'; % Talker D - mcg
s13 = 'MS1_F99_C1.mp4'; % Talker E - ba * CONTROL TALKER
s14 = 'MS1_F99_C3.mp4'; % Talker E - ga * CONTROL TALKER
s15 = 'MS1_F99_C2.mp4'; % Talker E - da * CONTROL TALKER

%% Create matrix: STIMULI TALKER SYLLABLE MCG/CONTROL TALKER
stimuli = zeros(n_stim, 4);

stimuli(:, 1) = 1:n_stim; % stimuli numbers
talker_id = repmat(1:n_talker, n_cond, 1);
stimuli(:,2) = reshape(talker_id, n_stim, 1); % id talker (#1-5)

syllables = repmat(1:n_cond, n_talker, 1);
stimuli(:,3) = reshape(syllables', n_stim, 1); % 1 = ba, 2 = ga, 3 = mcg/da
stimuli(1:mcg_talkers*n_cond,4) = 1; % 1 = mcg talker; 0 = control talker

%% Double McG talker trials (reps_m = 10)
for k = 1:(length(stimuli)) % loop through each line in stim
    if (stimuli(k, 4) == 1) == 1 % is it mcg talker?
        stimuli = [stimuli; stimuli(k,:)]; % if it is - add line again
    else
        stimuli = stimuli; % otherwise leave it the same
    end
end

% x dyn_reps to get to all trials
talk = unique(stimuli(:,2)); % number of talkers - vector
check_tot = 0;
stim = 1:n_stim; % vector of stimuli numbers
syll = unique(stimuli(:,3)); % vector of syll 

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
                          
           seq_m = [3 3 3]; % does mcgurk repeat 3 times in a row?
           pat_m = findpattern(shuffled_set(:,3),seq_m);
           
           if isempty(pat_m)
               check_m = 0;
           else
               check_m = 1;
           end

            if sum(check_t) == 0 && sum(check_s) == 0 && check_m == 0
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


% %% Randomize stimuli & check 
% % currently - checks if three of same talker in a row 
% % add so not same stimulus in a row 
% 
% new = [];
% 
% talk_order = stim_mat(:,2);
% 
% for t=1:n_talker % loop through each talker 
%     
%     flag = [t t]; % look for same talker 3x in a row
%     
%     for s = 1:(length(stim_mat(:,2))-2) % loop through each trial 
%         
%         talk1to3 = talk_order(s:s+2);
% %         test(s) = isequal(talk1to3, flag); % check if flag appears in talk order 
%                 
%         if talk1to3 == flag
% %             talk_order = Shuffle(talk_order);
%             stim_mat = Shuffle(stim_mat,2);
%             talk_order = stim_mat(:,2)
%             new(t,s) = 1;
%         else
%             new(t,s) = 0;
%             talk_order = talk_order;
%             %stim_mat = stim_mat;
%         end 
%         
%     end 
%     
% end 
% 
% 
% 
%% Retrieve file names
stim_file = cell(length(stim_mat), 1);

for j = 1:(length(stim_mat)) % loop through each line in stim_mat
    
    stim_num(j) = stim_mat(j,1);  % get stim number from shuffled matrix
    str = ['s' num2str(stim_num(j))]; % get stim string
    stim = genvarname(str); % str -> variable
    stim_file{j} = [eval(stim)]; % get file name 
end
