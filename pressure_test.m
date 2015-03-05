function data = pressure_test(trial_sequence, varargin)

% This function is for controlling the LabView program to deliver pressure
% pain and collecting ratings (continuous or one-time ratings)
%
% Usage:
% -------------------------------------------------------------------------
% data = pressure_test(trial_sequence)
%
% Inputs:
% -------------------------------------------------------------------------
% trial_sequence trial_sequence should provide information of intensity,
%                duration, repetition of simulation, rating and scale type 
%                you want to use, and cue/iti duration. 
%
% The generic form of trial sequence:
% trial_sequence{run_number}{trial_number} = ...
%          {intensity(four digits:string), duration(four digits:string),
%           repetition_number, rating_type, scale_type, cue_duration,
%           inter_trial_interval};
% For details, see example below.
%
% Optional input:
% -------------------------------------------------------------------------
% 'post_st_rating_dur'  If you are collecting continuous rating, using this 
%                       option, you can specify the duration for the 
%                       post-stimulus rating. The default is 3 seconds.
%                       (e.g., 'post_st_rating_dur', duration_in_seconds)
% 'explain_scale'       If you want to show rating scale before starting
%                       the experiment, you can use this option.
%                       (e.g., 'explain_scale', {'cont', 'overall_int'}, 'linear')
%
% Outputs:
% -------------------------------------------------------------------------
% data.
%
%
%
%
% Example:
% -------------------------------------------------------------------------
% trial_sequence{1}{1} = {'0004', '0010', '1', 'cont', 'line', '2', '5'};
%     ----------------------------
%     {1}{1}: first run, first trial
%     '0004': intensity 4kg
%     '0010': duration 10 seconds 
%     '1': one time administration (no repetition); could be any number
%     'cont': continuous rating; there are three possible options
%             'cont': continuous rating
%             'overall_int': overall intensity rating (after pressure pain stimulation)
%             'overall_pleasant': overall intensity rating (after pressure pain stimulation)
%             'all': continous + overall intensity + overall pleasantness ratings
%             You can put two things in a cell. 
%             e.g., trial_sequence{1}{1} = {'0004', '0010', '1', {'cont', 'overall_int'}, 'line', '2', '5'};
%     'line': line scale; there are three possible options
%             'line': line scale
%             'linear': right-angled triangle
%             'lms': Labeled magnitude scale (Bartoshuk)
%     '2': cue duration 2 seconds
%     '5': inter-trial interval 5 seconds
%     ----------------------------
% trial_sequence{1}{2} = {'0006', '0010', '1', 'both', 'linear', '3', '4'};
%     ----------------------------
%     {1}{2}: first run, second trial
%     Continuous and overall ratings using linear scale. 
%     ----------------------------
% trial_sequence{2}{1} = {'04.5', '07.4', '1', 'overall', 'lms', '5', '6'};
%     ----------------------------
%     {2}{1}: second run, first trial
%     Overall ratings using LMS scale. 
%     ----------------------------
% trial_sequence{2}{2} = {'0006', '0001', '8', 'cont', 'line', '3', '6'};
%     ----------------------------
%     {2}{2}: second run, second trial
%     This delivers eight quick repetitions of short 6kg pressure (1 second). 
%     ----------------------------
% data = pressure_test(trial_sequence);
%
% -------------------------------------------------------------------------
% Copyright (C) 1/10/2015, Wani Woo


%% SETUP: global
global theWindow W H; % window property
global white red orange bgcolor; % color
global t r; % pressure device udp channel
global window_rect prompt lb rb scale_W anchor_y anchor_y2 anchor promptW promptH; % rating scale

%% SETUP: Screen
bgcolor = 100;
window_rect = get(0, 'MonitorPositions'); % full screen
% window_rect = [0 0 1000 600]; % specific size
W = window_rect(3); %width of screen
H = window_rect(4); %height of screen
font = 'Helvetica';
fontsize = 24;
white = 255;
red = [158 1 66];
orange = [255 164 0];
lb = W/4; % rating scale left and right bounds
rb = (3*W)/4;
scale_W = (rb-lb).*0.1;
anchor = [0.014 0.061 0.172 0.354 0.533].*(rb-lb)+lb;
% scale_name = {'line', 'linear', 'LMS'};
post_stimulus_t = 3; % post-stimulus continuous rating seconds
doexplain_scale = false;

for i = 1:length(varargin)
    if ischar(varargin{i})
        switch varargin{i}
            % functional commands
            case {'post_st_rating_dur', 'post_st_rating'}
                post_stimulus_t = varargin{i+1};
            case {'explain_scale'}
                doexplain_scale = true;
                exp_scale.inst = varargin{i+1};
                exp_scale.scale = varargin{i+2};
        end
    end
end

%% SETUP: instructions
cont = 1; overall_int = 2; overall_pleasant = 3;
prompt{cont} = 'Please rate pain intensity continuously by moving mouse right and left.';
prompt{overall_int} = 'Please rate the overall pain intensity you experienced by clicking your mouse.';
prompt{overall_pleasant} = 'Please rate the overall unpleasantness you experienced by clicking your mouse.';

%% SETUP: DATA and Subject INFO
[fname, start_line, SID] = subjectinfo_check; % subfunction
if exist(fname, 'file'), load(fname); end

% save data using the canlab_dataset object
data.version = 'PressurePain_pilot_v1_01/09/2015_ChoongWanWoo';
data.subject = SID;
data.datafile = fname;

%% SETUP: Experiment
[run_num, trial_num, runstart, trial_starts, rating_types] = parse_trial_sequence(trial_sequence, start_line);
% cue_t = 2; % seconds for cues 

try
    % 1. SETUP: Pressure device
    
    t=udp('localhost',61557); % open udp channels
    r=udp('localhost',61158,'localport', 61556);
    
    fopen(t);
    fopen(r);
    fwrite(t, '0005,0010,o'); % open the remote channel
    
    % 2. START: Screen
    % whichScreen = max(Screen('Screens'));
	theWindow = Screen('OpenWindow', 0, bgcolor, window_rect); % start the screen
    HideCursor;
    Screen('TextFont', theWindow, font); % setting font
    Screen('TextSize', theWindow, fontsize);
    [fixW, fixH] = Screen(theWindow,'DrawText','+',0,0); % draw text in black and return new pen location; a trick to get display width and height
    promptW = Screen(theWindow, 'DrawText',prompt{1},0,0);
    promptH = 30;
    
    % anchor y location
    anchor_y = H/2+10+scale_W;
    anchor_y2 = H/2+10+scale_W+promptH;
    
    % 3. EXPLAIN SCALES
    if doexplain_scale
        explain_scale(exp_scale);
    end
    
    % 4. START: RUN
    for run_i = runstart:run_num % run starts
        
        for tr_i = trial_starts(run_i):trial_num(run_i) % trial starts
            
            if tr_i == 1
                while (1)
                    [~,~,keyCode] = KbCheck;
                    if keyCode(kbName('r'))==1
                        break
                    elseif keyCode(kbName('q'))==1
                        abort_man;
                    end
                    
                    % HERE: YOU CAN ADD MESSAGES FOR EACH RUN
                    if run_i <= run_num
                        Run_start_text{1} = 'If the participant is ready for starting the run, please press r.';
                        Run_start_text{2} = ' ';
                    end
                    
                    runtextW = 0;
                    for jj = 1:numel(Run_start_text)
                        eval(['runtextW' num2str(jj) '= Screen(''DrawText'',theWindow,Run_start_text{jj},0,0);']);
                        eval(['runtextW = max(runtextW, runtextW' num2str(jj) ');']);
                    end
                    
                    Screen(theWindow,'FillRect',bgcolor, window_rect);
                    for jj = 1:numel(Run_start_text)
                        Screen('DrawText',theWindow,Run_start_text{jj},W/2-runtextW/2,H/2+promptH*(jj-1)-150,white);
                    end
                    Screen('Flip', theWindow);
                    
                end
            end
            
            % cue or fixation cross
            Screen(theWindow,'FillRect',bgcolor, window_rect);
            Screen(theWindow,'DrawText','+',W/2-fixW/2, H/2-fixH/2,255);
            Screen('Flip', theWindow);
            cue_t = str2double(trial_sequence{run_i}{tr_i}{6});
            WaitSecs(cue_t);
            
            % SETUP: Trial stimulus
            int = trial_sequence{run_i}{tr_i}{1};
            dur = trial_sequence{run_i}{tr_i}{2};
            rep = str2double(trial_sequence{run_i}{tr_i}{3});
            
            % RECORD: Trial Info 
            data.dat{run_i}{tr_i}.intensity = int;
            data.dat{run_i}{tr_i}.duration = dur;
            data.dat{run_i}{tr_i}.repetition = rep;
            data.dat{run_i}{tr_i}.scale = rating_types.type{run_i}{tr_i};
            
            % START: Trial
            Screen(theWindow,'DrawText','+',W/2-fixW/2, H/2-fixH/2,200);
            Screen('Flip', theWindow);
            Screen(theWindow,'FillRect',bgcolor, window_rect); % Fixation cross disappear
            draw_scale(rating_types.type{run_i}{tr_i});
            
            % For continuous rating, show rating instruction before stimulus starts
            if rating_types.docont{run_i}(tr_i)
                Screen('DrawText', theWindow, prompt{cont},W/2-promptW/2,H/2-promptH/2-150,white); 
                Screen('Flip', theWindow); 
                WaitSecs(1);
            end
            
            rec_i = 0;
            for st_i = 1:rep
                
                % deliever pressure
                eval(['fwrite(t, ''' int ',' dur ',t'');']);
                
                % RECORD: Time stamp
                data.dat{run_i}{tr_i}.st_timestamp(st_i,:) = datestr(clock, 0);
                if st_i == 1
                    start_t = GetSecs; 
                    SetMouse(lb,H/2); % set mouse at the left
                end
                
                message_1 = deblank(fscanf(r));
                if strcmp(message_1,'Read Error')
                    error(message_1);
                else
                    data.dat{run_i}{tr_i}.logfile{st_i} = message_1;
                end
                
                % CONTINUOUS RATING
                if rating_types.docont{run_i}(tr_i)
                    
                    % START: Instruction and rating scale
                    deltat = 0;
                    start_t2 = GetSecs;
                    
                    while deltat <= (str2double(dur)+str2double(dur)/10) % collect data for the duration+150ms
                        deltat = GetSecs - start_t2; % message = fscanf(r);
                        % message = input('test, s:stop, others:continue? ', 's'); % for test
                        rec_i = rec_i+1; % the number of recordings
                        
                        %Track Mouse coordinate
                        x = GetMouse(theWindow);
                        
                        if x < lb
                            x = lb;
                        elseif x > rb
                            x = rb;
                        end
                        
                        cur_t = GetSecs;
                        data.dat{run_i}{tr_i}.time_from_start(rec_i,1) = cur_t-start_t;
                        data.dat{run_i}{tr_i}.cont_rating(rec_i,1) = (x-lb)./(rb-lb);
                        
                        draw_scale(rating_types.type{run_i}{tr_i}); % draw scale
                        Screen('DrawText', theWindow, prompt{cont},W/2-promptW/2,H/2-promptH/2-150,white); 
                        Screen('DrawLine', theWindow, orange, x, H/2, x, H/2+scale_W, 7);
                        Screen('Flip', theWindow); 
                    end
                else
                    WaitSecs(str2double(dur)+str2double(dur)/10);
                end
                
                message_2 = deblank(fscanf(r));
                if ~strcmp(message_2, 's') 
                    disp(message_2); 
                    error('message_2 is not s.'); 
                end % make sure if the stimulus ends
                
%                % GIVE A LITTLE INTERVAL BETWEEN REPETITION: DO WE NEED THIS?
%                 if rating_types.docont{run_i}(tr_i)
%                     
%                     % START: Instruction and rating scale
%                     deltat = 0;
%                     start_t2 = GetSecs;
%                     
%                     while deltat <= 0.1 % collect data for the duration+150ms
%                         deltat = GetSecs - start_t2; % message = fscanf(r);
%                         % message = input('test, s:stop, others:continue? ', 's'); % for test
%                         rec_i = rec_i+1; % the number of recordings
%                         
%                         %Track Mouse coordinate
%                         x = GetMouse(theWindow);
%                         
%                         if x < lb
%                             x = lb;
%                         elseif x > rb
%                             x = rb;
%                         end
%                         
%                         cur_t = GetSecs;
%                         data.dat{run_i}{tr_i}.time_from_start(rec_i,1) = cur_t-start_t;
%                         data.dat{run_i}{tr_i}.cont_rating(rec_i,1) = (x-lb)./(rb-lb);
%                         
%                         draw_scale(rating_types.type{run_i}{tr_i}); % draw scale
%                         Screen('DrawText', theWindow, prompt{cont},W/2-promptW/2,H/2-promptH/2-150,white); 
%                         Screen('DrawLine', theWindow, orange, x, H/2, x, H/2+scale_W, 7);
%                         Screen('Flip', theWindow); 
%                     end
%                 else
%                     WaitSecs(0.1);
%                 end
            end
            
            end_t = GetSecs;
            data.dat{run_i}{tr_i}.total_st_dur = end_t - start_t;
            
            % POST-STIMULUS DELAY
            if rating_types.docont{run_i}(tr_i)
                
                start_t2 = GetSecs;
                deltat = 0;
                while deltat <= post_stimulus_t %
                    deltat = GetSecs - start_t2; % message = fscanf(r);
                    rec_i = rec_i+1; % the number of recordings
                    
                    %Track Mouse coordinate
                    x = GetMouse(theWindow);
                    
                    if x < lb
                        x = lb;
                    elseif x > rb
                        x = rb;
                    end
                    
                    cur_t = GetSecs;
                    data.dat{run_i}{tr_i}.time_from_start(rec_i,1) = cur_t-start_t;
                    data.dat{run_i}{tr_i}.cont_rating(rec_i,1) = (x-lb)./(rb-lb);
                    
                    draw_scale(rating_types.type{run_i}{tr_i}); % draw scale
                    Screen('DrawText', theWindow, prompt{cont},W/2-promptW/2,H/2-promptH/2-150,white);
                    Screen('DrawLine', theWindow, orange, x, H/2, x, H/2+scale_W, 7);
                    Screen('Flip', theWindow);
                end
            else
                WaitSecs(post_stimulus_t);
            end
            
            % total_record_dur contain stimulus + post-stimulus duration
            end_t = GetSecs;
            data.dat{run_i}{tr_i}.total_record_dur = end_t - start_t;
            
            % OVERALL INTENSITY RATING
            if rating_types.dooverall_int{run_i}(tr_i)
                start_t = GetSecs;
                Screen('FillRect', theWindow, bgcolor, window_rect); % clear the screen
                Screen('Flip', theWindow);
                WaitSecs(1);
                
                SetMouse(lb,H/2); % set mouse at the left
                    
                while (1) % button
                    [x,~,button] = GetMouse(theWindow);
                    if x < lb
                        x = lb;
                    elseif x > rb
                        x = rb;
                    end
                    
                    draw_scale(rating_types.type{run_i}{tr_i}); % draw scale
                    Screen('DrawText',theWindow, prompt{overall_int},W/2-promptW/2,H/2-promptH/2-150,white);
                    Screen('DrawLine', theWindow, orange, x, H/2, x, H/2+scale_W, 7);
                    Screen('Flip', theWindow);
                    
                    if button(1), break, end
                end
                
                % freeze the screen 1 second with red line
                draw_scale(rating_types.type{run_i}{tr_i}); % draw scale
                Screen('DrawText',theWindow, prompt{overall_int},W/2-promptW/2,H/2-promptH/2-150,white); 
                Screen('DrawLine', theWindow, red, x, H/2, x, H/2+scale_W, 7);
                Screen('Flip', theWindow);
                WaitSecs(1);
                
                end_t = GetSecs;
                data.dat{run_i}{tr_i}.overall_int_rating = (x-lb)./(rb-lb);
                data.dat{run_i}{tr_i}.overall_int_rating_rt = end_t-start_t;
            end
            
            % OVERALL PLEASANTNESS RATING
            if rating_types.dooverall_int{run_i}(tr_i)
                start_t = GetSecs;
                Screen('FillRect', theWindow, bgcolor, window_rect); % clear the screen
                Screen('Flip', theWindow);
                WaitSecs(1);
                
                SetMouse(lb,H/2); % set mouse at the left
                    
                while (1) % button
                    [x,~,button] = GetMouse(theWindow);
                    if x < lb
                        x = lb;
                    elseif x > rb
                        x = rb;
                    end
                    
                    draw_scale(rating_types.type{run_i}{tr_i}); % draw scale
                    Screen('DrawText',theWindow, prompt{overall_pleasant},W/2-promptW/2,H/2-promptH/2-150,white);
                    Screen('DrawLine', theWindow, orange, x, H/2, x, H/2+scale_W, 7);
                    Screen('Flip', theWindow);
                    
                    if button(1), break, end
                end
                
                % freeze the screen 1 second with red line
                draw_scale(rating_types.type{run_i}{tr_i}); % draw scale
                Screen('DrawText',theWindow, prompt{overall_pleasant},W/2-promptW/2,H/2-promptH/2-150,white); 
                Screen('DrawLine', theWindow, red, x, H/2, x, H/2+scale_W, 7);
                Screen('Flip', theWindow);
                WaitSecs(1);
                
                end_t = GetSecs;
                data.dat{run_i}{tr_i}.overall_pleasant_rating = (x-lb)./(rb-lb);
                data.dat{run_i}{tr_i}.overall_pleasant_rating_rt = end_t-start_t;
            end
            
            % INTER-TRIAL INTERVAL
            Screen('FillRect', theWindow, bgcolor, window_rect); % basically, clear the screen
            Screen('Flip', theWindow);
            iti = str2double(trial_sequence{run_i}{tr_i}{7});
            WaitSecs(iti);
            
            save(data.datafile, 'data');
        end % trial ends
        
        % message between runs
        while (1) 
            [~,~,keyCode] = KbCheck;
            if keyCode(kbName('space'))==1
                break
            elseif keyCode(kbName('q'))==1
                abort_man;
            end
            
            % HERE: YOU CAN ADD MESSAGES FOR THE NEXT RUN
            if run_i < run_num-1
                Run_end_text{1} = ['This is the end of the run ' num2str(run_i) '.'];
                Run_end_text{2} = 'If the participant is ready for the next run, please press Space.';
            elseif run_i == run_num-1
                Run_end_text{1} = ['This is the end of the run ' num2str(run_i) '. '];
                Run_end_text{2} = 'The next run is the last run of the experiment.';
                Run_end_text{3} = 'In the next run, please remove your finger from the device as quickly as possible.';
            else
                Run_end_text{1} = 'This is the end of the experiment.';
                Run_end_text{2} = 'To finish the experiment, please press Space.';
            end
            
            runtextW = 0;
            for jj = 1:numel(Run_end_text)
                eval(['runtextW' num2str(jj) '= Screen(''DrawText'',theWindow,Run_end_text{jj},0,0);']);
                eval(['runtextW = max(runtextW, runtextW' num2str(jj) ');']);
            end
            Screen(theWindow,'FillRect',bgcolor, window_rect); 
            
            for jj = 1:numel(Run_end_text)
                Screen('DrawText',theWindow,Run_end_text{jj},W/2-runtextW/2,H/2+promptH*(jj-1)-150,white); 
            end
            Screen('Flip', theWindow);
        end
        
    end % run ends
    
    fclose(t);
    fclose(r);
    
    Screen('CloseAll');
    disp('Done');
    save(data.datafile, 'data');
    
catch err
    % ERROR 
    disp(err);
    disp(err.stack(end));
    fclose(t);
    fclose(r);
    abort_error; 
end

end

%======================================
%             SUBFUNCTIONS            
%======================================

% SUBJECT INFORMATION: file exists?
function [fname,start_line,SID] = subjectinfo_check

% Subject ID    
fprintf('\n');
SID = input('Subject ID? ','s');
    
% check if the data file exists
fname = fullfile('pressure_data', ['s' SID '.mat']);
if ~exist('pressure_data', 'dir')
    mkdir('pressure_data');
    whattodo = 1;
else
    if exist(fname, 'file')
        str = ['The Subject ' SID ' data file exists. Press a button for the following options'];
        disp(str);
        whattodo = input('1:Save new file, 2:Save the data from where we left off, Ctrl+C:Abort? ');
    else
        whattodo = 1;
    end
end
    
if whattodo == 2
    load(fname);
    start_line = 1;
    for i = 1:numel(data.dat)
        start_line = start_line + numel(data.dat{i});
    end
else
    start_line = 1;
end
   
end

% PARSING TRIAL SEQUENCE: each trial info
function [run_num, trial_num, runstart, trial_start, rating_types] = parse_trial_sequence(trial_sequence, start_line)

run_num = numel(trial_sequence);
idx = zeros(run_num+1,1);
trial_num = zeros(run_num,1);
trial_start = ones(run_num,1);

for i = 1:run_num
    trial_num(i,1) = numel(trial_sequence{i});
    idx(i+1) = idx(i) + trial_num(i);
    
    if idx(i) < start_line && idx(i+1) >= start_line
        runstart = i;
        trial_start(i,1) = start_line - idx(i);
    end
    
    for j = 1:trial_num(i)
        % (4) rating type1, continuous, overall_int, overall_pleasant, or all
        rating_types.docont{i}(j,1) = false;
        rating_types.dooverall_int{i}(j,1) = false;
        rating_types.dooverall_pleasant{i}(j,1) = false;
        
        rating_types.docont{i}(j,1) = any(strcmp(trial_sequence{i}{j}{4}, 'cont') | strcmp(trial_sequence{i}{j}{4}, 'all'));
        rating_types.dooverall_int{i}(j,1) = any(strcmp(trial_sequence{i}{j}{4}, 'overall_int') | strcmp(trial_sequence{i}{j}{4}, 'all'));
        rating_types.dooverall_pleasant{i}(j,1) = any(strcmp(trial_sequence{i}{j}{4}, 'overall_pleasant') | strcmp(trial_sequence{i}{j}{4}, 'all'));
        
        % (5) rating type2, line, linear, or lms (pseudo-log)
        rating_types.type{i}{j} = trial_sequence{i}{j}{5};
    end
end

end

% DRAWRING SCALES
function draw_scale(scale)
    
global theWindow W H; % window property
global white red orange; % color
global t r; % pressure device udp channel
global lb rb scale_W anchor_y anchor_y2 anchor; % rating scale

switch scale
    case 'line'
        xy = [lb lb lb rb rb rb; H/2 H/2+scale_W H/2+scale_W/2 H/2+scale_W/2 H/2 H/2+scale_W];
        Screen(theWindow,'DrawLines', xy, 5, 255);
        Screen(theWindow,'DrawText','Not',lb-50,anchor_y,255);
        Screen(theWindow,'DrawText','at all',lb-50,anchor_y2,255);
        Screen(theWindow,'DrawText','Worst',rb-50,anchor_y,255);
        Screen(theWindow,'DrawText','imaginable',rb-50,anchor_y2,255);
        % Screen('Flip', theWindow);
    case 'linear'
        xy = [lb H/2+scale_W; rb H/2+scale_W; rb H/2];
        Screen(theWindow, 'FillPoly', 255, xy);
        Screen(theWindow,'DrawText','Not',lb-50,anchor_y,255);
        Screen(theWindow,'DrawText','at all',lb-50,anchor_y2,255);
        Screen(theWindow,'DrawText','Worst',rb-50,anchor_y,255);
        Screen(theWindow,'DrawText','imaginable',rb-50,anchor_y2,255);
        % Screen('Flip', theWindow);
    case 'lms'
        xy = [lb H/2+scale_W; rb H/2+scale_W; rb H/2];
        Screen(theWindow, 'FillPoly', 255, xy);
        Screen(theWindow,'DrawText','Not',lb-100,anchor_y-20,255);
        Screen(theWindow,'DrawText','at all',lb-100,anchor_y2-20,255);
        Screen(theWindow,'DrawText','Worst',rb-50,anchor_y,255);
        Screen(theWindow,'DrawText','imaginable',rb-50,anchor_y2,255);
        for i = 1:5
            Screen('DrawLine', theWindow, 0, anchor(i), H/2+scale_W, anchor(i), H/2, 2);
        end
        Screen(theWindow,'DrawText','barely',anchor(1)-30,anchor_y+scale_W/2.5,255);
        Screen(theWindow,'DrawText','detectable',anchor(1)-30,anchor_y2+scale_W/2.5,255);
        Screen(theWindow,'DrawText','weak',anchor(2)-10,anchor_y,255);
        Screen(theWindow,'DrawText','moderate',anchor(3),anchor_y,255);
        Screen(theWindow,'DrawText','strong',anchor(4),anchor_y,255);
        Screen(theWindow,'DrawText','very',anchor(5),anchor_y,255);
        Screen(theWindow,'DrawText','strong',anchor(5),anchor_y2,255);
        % Screen('Flip', theWindow);
end

end            

function explain_scale(exp_scale)

global theWindow W H; % window property
global white red orange bgcolor; % color
global t r; % pressure device udp channel
global window_rect prompt lb rb scale_W anchor_y anchor_y2 anchor promptW promptH; % rating scale

prompt_ex{1} = 'Scale example: Experimenter will explain how to use the scale, and press Space.';
prompt_ex{2} = 'Scale example: Please practice rating, and when you are done, please click the mouse.';
prompt_ex{3} = 'Great job! If you are ready for the next step, please press Space.';

% some additional prompts
prompt_ex{4} = 'Welcome! The experiment will starts with an explanation of the rating scales.';
prompt_ex{5} = 'If you are ready, please press Space.';
prompt_ex{6} = 'Note that during the actual experiment, you don''t need to click the mouse for the continuous rating.';
prompt_ex{7} = 'Great! We''re done with the practice. Now, we are about to start the actual pressure pain experiment.';
prompt_ex{8} = 'If you are ready for pressure pain, please press Space.';

for i = 1:numel(prompt_ex), prompt_ex_W{i} = Screen(theWindow, 'DrawText', prompt_ex{i},0,0); end

for i = 1:numel(exp_scale.inst)
    
    switch exp_scale.inst{i}
        case 'cont'
            prompt_n = 1;
        case 'overall_int'
            prompt_n = 2;
        case 'overall_plesant'
            prompt_n = 3;
    end
    
    % first: START page
    if i == 1
        while (1) % space
            [~,~,keyCode] = KbCheck;
            
            if keyCode(kbName('space'))==1
                break
            elseif keyCode(kbName('q'))==1
                abort_man;
            end
            
            Screen('DrawText',theWindow, prompt_ex{4},W/2-prompt_ex_W{4}/2, 100, white);
            Screen('DrawText',theWindow, prompt_ex{5},W/2-prompt_ex_W{5}/2, 150, white);
            Screen('Flip', theWindow);
        end
    end
    
    % EXPLAIN
    while (1) % space
        [~,~,keyCode] = KbCheck;
        
        if keyCode(kbName('space'))==1
            break
        elseif keyCode(kbName('q'))==1
            abort_man;
        end
        
        draw_scale(exp_scale.scale); % draw scale
        Screen('DrawText',theWindow, prompt_ex{1},W/2-prompt_ex_W{1}/2,100,orange);
        Screen('DrawText',theWindow, prompt{prompt_n},W/2-promptW/2,H/2-promptH/2-150,white);
        Screen('Flip', theWindow);
    end
    
    % PRACTICE
    Screen(theWindow,'FillRect',bgcolor, window_rect);
    SetMouse(lb,H/2); % set mouse at the left
    while (1) % button
        [x,~,button] = GetMouse(theWindow);
        if x < lb
            x = lb;
        elseif x > rb
            x = rb;
        end
        
        draw_scale(exp_scale.scale); % draw scale
        Screen('DrawText',theWindow, prompt_ex{2},W/2-prompt_ex_W{2}/2,100,orange);
        Screen('DrawText',theWindow, prompt{prompt_n},W/2-promptW/2,H/2-promptH/2-150,white);
        Screen('DrawLine', theWindow, orange, x, H/2, x, H/2+scale_W, 7);
        Screen('Flip', theWindow);
        
        if button(1), break, end
    end
    
    % freeze the screen 1 second with red line
    draw_scale(exp_scale.scale); % draw scale
    Screen('DrawText',theWindow, prompt_ex{2},W/2-prompt_ex_W{2}/2,100,orange);
    Screen('DrawText',theWindow, prompt{prompt_n},W/2-promptW/2,H/2-promptH/2-150,white);
    Screen('DrawLine', theWindow, red, x, H/2, x, H/2+scale_W, 7);
    Screen('Flip', theWindow);
    WaitSecs(1);
    
    Screen(theWindow,'FillRect',bgcolor, window_rect);
    % Move to next
    if i < numel(exp_scale.inst)
        while (1) % space
            [~,~,keyCode] = KbCheck;
            
            if keyCode(kbName('space'))==1
                break
            elseif keyCode(kbName('q'))==1
                abort_man;
            end
            
            Screen('DrawText',theWindow, prompt_ex{3},W/2-prompt_ex_W{3}/2,100,orange);
            if prompt_n == 1
                Screen('DrawText',theWindow, prompt_ex{6},W/2-prompt_ex_W{6}/2,150,orange);
            end
            Screen('Flip', theWindow);
        end
    else
        while (1) % space
            [~,~,keyCode] = KbCheck;
            
            if keyCode(kbName('space'))==1
                break
            elseif keyCode(kbName('q'))==1
                abort_man;
            end
            
            Screen('DrawText',theWindow, prompt_ex{7},W/2-prompt_ex_W{7}/2,100,orange);
            Screen('DrawText',theWindow, prompt_ex{8},W/2-prompt_ex_W{8}/2,150,orange);
            
            if prompt_n == 1
                Screen('DrawText',theWindow, prompt_ex{6},W/2-prompt_ex_W{6}/2,300,orange);
            end
            
            Screen('Flip', theWindow);
        end
    end
end

end

% ABORT OPTIONS: ERROR
function abort_error
try
    fclose(t); fclose(r);
catch
end
ShowCursor; %unhide mouse
Screen('CloseAll'); %relinquish screen control
disp('Experiment aborted by error') %present this text in command window
end

% ABORT OPTIONS: MANUAL (by "q")
function abort_man
try
    fclose(t); fclose(r);
catch
end
ShowCursor; %unhide mouse
Screen('CloseAll'); %relinquish screen control
psychrethrow(psychlasterror);
disp('Experiment aborted by escape sequence'); %present this text in command window
end