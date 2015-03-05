function [ts, energy_range_trial] = generate_ts_pressure(run_number, trial_number, varargin)

% ts = generate_ts_pressure(run_number, trial_number, ...
%            rating[{'cont'}|'overall'|'both'], scale['line'|{'linear'}|'lms'], 
%            optional['post_st_rating_dur', duration_in_seconds]);
%
% trial_sequence{run_number}{trial_number} = ...
%          {intensity(four digits:string), duration(four digits:string),
%           repetition_number, rating_type, scale_type, cue_duration,
%           inter_trial_interval};
%
% Example usage
% ts = generate_ts_pressure(2, 15, 'post_st_rating_dur', 8);

if run_number < 2
    warning('Run_number should be at least two. Thus, run_number is 2 now.')
    run_number = 2;
end

run_number = run_number-1;

% rating = 'cont';
scale = 'linear';
post_stimulus_t = 8;

for i = 1:length(varargin)
    if ischar(varargin{i})
        switch varargin{i}
            % functional commands
            case {'cont', 'overall_int', 'overall_pleasant', 'all'}
                if ~exist('rating', 'var')
                    rating{1} = varargin{i};
                else
                    rating{end+1} = varargin{i};
                end
            case {'line', 'linear', 'lms'}
                scale = varargin{i};
            case {'post_st_rating_dur', 'post_st_rating'}
                post_stimulus_t = varargin{i+1};
            case {'explain_scale'}
                exp_option = true; 
        end
    end
end

if ~exist('rating', 'var'), rating = 'cont'; end

if ~exp_option
    str_suggest = ['data = pressure_test(ts, ''post_st_rating_dur'', ' num2str(post_stimulus_t) ');'];
else
    str_suggest = ['data = pressure_test(ts, ''post_st_rating_dur'', ' num2str(post_stimulus_t) ', ''explain_scale'', {''cont'', ''overall_int'', ''overall_plesant''}, ''linear'');'];
end

% default setting
energy_range= 30:10:60;
cue_mean = 2;
iti_mean = 8;

% rng('shuffle')

% possible choices
opts = {'int', 'dur', 'rep', 'rep_dur'};
int = [3, 4, 5, 6, 7];
dur = [8, 10, 12, 14, 16, 18];
rep = [0, 1];
rep_dur = [1, 2, 3];

rat = '';
for i = 1:numel(rating)
    rat = [rat ' ''' rating{i} ''''];
end
rat(1) = '{';
rat(end+1) = '}';

total_n = (run_number-1)*trial_number; % run 1 is always 5 trials
[~,idx] = sort(rand(numel(energy_range),ceil(total_n/numel(energy_range))));
energy_range_trial = energy_range(idx(:));

%% idea from 1/22 
% done: have a same sequence for everyone (3~7, 10 seconds)
% done: longer than 16-22 seconds iti
% done: then, randomize - equal numbers of trials for each energy interval 
%
% need to test longer duration and multiple peaks
% remove the finger

for i = 1
    for j = 1:numel(int)
        tr.int = int(j);
        tr.dur = 10;
        cue = cue_mean + geornd(.5);
        iti = iti_mean + geornd(.5);
        
        str = sprintf('ts{i}{j} = {''%04d'',''%04d'', ''%d'', %s, ''%s'', ''%d'', ''%d''};', ...
            tr.int, tr.dur, 1, rat, scale, cue, iti);
        eval(str);
    end
end

ii = 0;
for i = 2:run_number
    for j = 1:trial_number
        ii = ii + 1;
        tr.int = Inf; tr.dur = Inf;
        while ~((tr.int * tr.dur) >= energy_range_trial(ii) && (tr.int * tr.dur) < (energy_range_trial(ii)+10))
            for jj = 1:numel(opts)
                eval(['tr.' opts{jj} ' = ' opts{jj} '(randi(numel(' opts{jj} ')));']);
            end
        end
        
        cue = cue_mean + geornd(.5);
        iti = iti_mean + geornd(.5);
        
        if tr.rep
            str = sprintf('ts{i}{j} = {''%04d'',''%04d'', ''%d'', %s, ''%s'', ''%d'', ''%d''};', ...
                tr.int, tr.rep_dur, floor(tr.dur/tr.rep_dur), rat, scale, cue, iti);
        else
            str = sprintf('ts{i}{j} = {''%04d'',''%04d'', ''%d'', %s, ''%s'', ''%d'', ''%d''};', ...
                tr.int, tr.dur, 1, rat, scale, cue, iti);
        end
        
        eval(str);
    
    end
end

for i = run_number + 1
    for j = 1:2
        tr.int = 8;
        tr.dur = 12;
        cue = cue_mean + geornd(.5);
        iti = iti_mean + geornd(.5);
        
        str = sprintf('ts{i}{j} = {''%04d'',''%04d'', ''%d'', %s, ''%s'', ''%d'', ''%d''};', ...
            tr.int, tr.dur, 1, rat, scale, cue, iti);
        eval(str);
    end
end

disp(' ');
disp('************************************************************************');
str_info{1} = 'RUN1: 5 Trials with 3-7kg/cm^2 for the 10 seconds duration';
for i = 2:run_number
    str_info{i} = ['RUN' num2str(i) ': ' num2str(trial_number) ' Trials with a pseuro-random sequence of stimuli'];
end
str_info{i+1} = ['RUN' num2str(i+1) ': 2 Trials with 8kg/cm^2 for the 12 seconds duration.'];
str_info{i+2} = ['        RUN' num2str(i+1) ' is to test the removal of a finger from the device.'];

for i = 1:numel(str_info)
    disp(str_info{i});
end

disp(' ');
disp('************Type the following line to run the pressure_test************');
disp(str_suggest);

end