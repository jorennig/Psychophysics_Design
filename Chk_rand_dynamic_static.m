%% Creates stimulus matrix; randomizes; checks randomizations; retrieves files
clear all
close all
clc

%% Define variables & stimuli.
n_talker = 4; 
n_dyn = 2; % number of movies per talker
n_stat = 1; % number of imgs per talker 
n_cond = n_dyn + n_stat; % conditions 
n_stim = n_talker * n_cond; % number of stimuli without repetition

reps_d = 10; % repetition dynamic faces
reps_s = 20; % repetition static faces

stim_dir = [pwd '/Stimuli/'];

s1 = 'MS1_F06_C1.mp4';  % talker A, dynamic, ba
s2 = 'MS1_F06_C3.mp4';  % talker A, dynamic, ga
s3 = 'MS1_F06_C3.bmp';  % talker A, static, img
s4 = 'MS1_M03_C1.mp4';  % talker B, dynamic, ba
s5 = 'MS1_M03_C3.mp4';  % talker B, dynamic, ga
s6 = 'MS1_M03_M.bmp';   % talker B, static, img
s7 = 'MS1_F08_C1.mp4';  % talker C, dynamic, ba
s8 = 'MS1_F08_C3.mp4';  % talker C, dynamic, ga
s9 = 'MS1_F08_C2.bmp';  % talker C, static, img
s10 = 'MS1_M02_C1.mp4'; % talker D, dynamic, ba
s11 = 'MS1_M02_C3.mp4'; % talker D, dynamic, ga
s12 = 'MS1_M02_C3.bmp'; % talker D, static, img

%% Create in order matrix [STIMULUS, TALKER, TYPE, SYLLABLE]
stimuli = zeros(n_stim, 4); % empty matrix - # stim x # info columns
stimuli(:, 1) = 1:n_stim; % stimuli numbers

talker_id = repmat(1:n_talker, n_cond, 1);
stimuli(:,2) = reshape(talker_id, n_stim, 1); % id talker (#1-4)

type = [repmat(1, n_dyn, 1); repmat(2, n_stat, 1)];
stimuli(:,3) = repmat(type, n_talker, 1); % dynamic = 1; static = 2

syll = repmat(1:n_cond, n_talker, 1);
stimuli(:,4) = reshape(syll', n_stim, 1); % 1 = ba, 2 = ga, 3 = img

%% Repetitions 
% Add extra repetitions for static stimuli *** THERE IS PROBABLY A BETTER WAY TO DO THIS *** 
for k = 1:(length(stimuli)) % loop through each line in stim
    if (stimuli(k, 3) == 2) == 1 % is it static?
        for i = 1:((reps_s/reps_d)-1) % add a new line based on the ratio b/w reps you 
        stimuli = [stimuli; stimuli(k,:)]; % if it's an img - add the line again at the end 
        end
    else
        stimuli = stimuli; % otherwise leave it the same
    end
end

% x dyn_reps to get to all trials
talk = unique(stimuli(:,2)); % number of talkers - vector
stim = 1:n_stim; % vector of stimuli numbers
check_tot = 0;

while check_tot == 0
    check_v_tot = zeros(numel(talk),1);
    stim_mat = [];

    for i = 1:reps_d % create randomized vector x #reps

       check = 0;
       while check == 0
           check_t = zeros(numel(talk),1);
           check_s = zeros(n_stim,1);

           shuffled_set = Shuffle(stimuli, 2); % shuffle a new set/repetition of stimuli

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

   % Check overall design 
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