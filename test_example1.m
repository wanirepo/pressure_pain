[ts, energy_tr] = generate_ts_pressure(5, 8, 'all', 'post_st_rating_dur', 9, 'explain_scale');

% ************************************************************************
% RUN1: 5 Trials with 3-7kg/cm^2 for the 10 seconds duration
% RUN2: 8 Trials with a pseuro-random sequence of stimuli
% RUN3: 8 Trials with a pseuro-random sequence of stimuli
% RUN4: 2 Trials with 8kg/cm^2 for the 12 seconds duration.
%         RUN4 is to test the removal of a finger from the device.
%  
% ************Type the following line to run the pressure_test************
% data = pressure_test(ts, 'post_st_rating_dur', 9, 'explain_scale');
% 
clear k;
jj = 0;
for i = 2:4
    for j = 1:numel(ts{i})
        jj = jj + 1;
        k(jj) = str2num(ts{i}{j}{1}) .* str2num(ts{i}{j}{2}) .* str2num(ts{i}{j}{3});
    end
end