%% DRIVER_LOW_HIGH_ALPHA7_SWITCH_SELECTED_PLOTS.m
% Requires NE_GKS_POISSON_STRENGTH_MODEL_GKL_SCALED_CAO_CORRECTED.m
% One continuous 52-state MODEL; plot only variables from the uploaded
% original driver (Vm, Ki, Ko, Nai, Nao, Ca_o, three volumes), plus Cl.
% LOW : gKs index 1, 20 Hz, strength index 4 (=3.0), alpha7 index 3 (=1.0)
% HIGH: gKs index 3, 10 Hz, strength index 2 (=1.5), alpha7 index 2 (=0.5)
% 10 segments, five LOW -> HIGH cycles, total time 4000 s.
% This driver assumes the model's alpha7 and strength lookup tables match
% the original driver. Check the model's internal arrays before interpreting.
% Run this script from a folder containing the required model .m file.

clear; clc; close all;
rng(1);
modelName = 'NE_GKS_POISSON_STRENGTH_MODEL_GKL_SCALED_CAO_CORRECTED';
assert(exist(modelName,'file') == 2, 'Put the required model .m file on the MATLAB path.');

%% Unchanged model options
A = 2; W = 2; amp = 2;
NE_level = 1.0;
k_NE_glia = 0.5;
k_gks_mod = 0.5;
poisson_dt_ms = 1.0;
dt_out = 1.0;                 % original 1-ms output retained for Vm spikes
solver_max_step_ms = min(dt_out,poisson_dt_ms/2);
opts = odeset('RelTol',1e-6,'AbsTol',1e-9,'MaxStep',solver_max_step_ms);

%% Conditions: one alpha7 protocol, NOT seven independent alpha7 runs
cycles = 5;
segment_duration_ms = 400000;
gks_sequence = repmat([1 3],1,cycles);                  % LOW, HIGH
frequency_sequence_Hz = repmat([20 10],1,cycles);
strength_index_sequence = repmat([4 2],1,cycles);      % 3.0, 1.5
strength_value_sequence = repmat([3.0 1.5],1,cycles);
alpha7_index_sequence = repmat([3 2],1,cycles);        % 1.0, 0.5
alpha7_value_sequence = repmat([1.0 0.5],1,cycles);
segment_starts_ms = (0:2*cycles-1)*segment_duration_ms;
dur_total = 2*cycles*segment_duration_ms;

%% Memory: retain the 52 states at 250-ms spacing temporarily for extracting
% only the selected variables, not for plotting or exporting all 52.
% Each segment is solved at 1 ms, then downsampled before joining.
plot_step_ms = 250;
plot_stride = round(plot_step_ms/dt_out);
assert(abs(plot_stride*dt_out-plot_step_ms)<1e-9,'plot_step_ms must be a multiple of dt_out');
assert(mod(segment_duration_ms,plot_step_ms)==0,'Choose plot interval dividing segment duration');

root_outdir = 'RESULTS_SELECTED_PLOTS_LOW_HIGH_ALPHA7_REPEAT5';
if ~exist(root_outdir,'dir'), mkdir(root_outdir); end

%% Store selected-variable trajectory and Vm envelope
T_plot_ms = [];
Y_plot = zeros(0,52);
T_env_ms = [];
Vm_min = [];
Vm_max = [];
zoom_centers_s = [20 100 170 820 900 970];
zoom_T = cell(size(zoom_centers_s));
zoom_Vm = cell(size(zoom_centers_s));
y_current = make_initial_condition_52();
assert(numel(y_current)==52,'Initial condition must contain 52 states.');

%% Segment protocol saved for reproducibility
segment_number = (1:2*cycles).';
start_s = segment_starts_ms(:)/1000;
end_s = start_s+segment_duration_ms/1000;
phase = repmat({'LOW';'HIGH'},cycles,1);
protocol = table(segment_number,start_s,end_s,phase,gks_sequence(:), ...
    frequency_sequence_Hz(:),strength_index_sequence(:),strength_value_sequence(:), ...
    alpha7_index_sequence(:),alpha7_value_sequence(:), ...
    'VariableNames',{'Segment','Start_s','End_s','Phase','gKs_Index', ...
    'Frequency_Hz','Strength_Index','Strength','Alpha7_Index','Alpha7'});
disp(protocol);
writetable(protocol,fullfile(root_outdir,'protocol.csv'));

%% Integrate sequentially, preserving the physiological state at switches
for seg = 1:numel(gks_sequence)
    t0 = segment_starts_ms(seg);
    t1 = t0+segment_duration_ms;
    fprintf('\nSegment %d/%d: %s, %.0f-%.0f s, gKs index %d, %g Hz, strength %g (index %d), alpha7 %g (index %d)\n', ...
        seg,numel(gks_sequence),phase{seg},t0/1000,t1/1000, ...
        gks_sequence(seg),frequency_sequence_Hz(seg), ...
        strength_value_sequence(seg),strength_index_sequence(seg), ...
        alpha7_value_sequence(seg),alpha7_index_sequence(seg));

    % Reset persistent Poisson array when rate/strength changes. This does
    % NOT reset y_current: all 52 state variables remain continuous.
    eval(['clear ' modelName]);
    modelFun = str2func(modelName);
    odeHandle = @(t,y) modelFun(t,y,A,W,alpha7_index_sequence(seg),amp, ...
        gks_sequence(seg),NE_level,k_NE_glia,k_gks_mod,dur_total, ...
        frequency_sequence_Hz(seg),poisson_dt_ms,strength_index_sequence(seg));
    tspan_segment = t0:dt_out:t1;
    [Tseg,Yseg] = ode15s(odeHandle,tspan_segment,y_current,opts);
    assert(size(Yseg,2)==52,'Expected exactly 52 state variables from the model.');
    y_current = Yseg(end,:).';

    % Vm can contain spikes invisible in a 250-ms downsample: store extrema.
    Vmseg = Yseg(:,1)-Yseg(:,21);
    vm_blocks = reshape(Vmseg(1:end-1),plot_stride,[]);
    T_env_ms = [T_env_ms; (t0+plot_step_ms*(0:size(vm_blocks,2)-1).'+plot_step_ms/2)]; %#ok<AGROW>
    Vm_min = [Vm_min; min(vm_blocks,[],1).']; %#ok<AGROW>
    Vm_max = [Vm_max; max(vm_blocks,[],1).']; %#ok<AGROW>
    clear vm_blocks;

    % Retain narrow 1-ms-resolution voltage windows (first & sixth segment).
    for iz=1:numel(zoom_centers_s)
        center_ms = zoom_centers_s(iz)*1000;
        if center_ms>=t0 && center_ms<=t1
            in_window = (Tseg>=center_ms & Tseg<=center_ms+3000);
            zoom_T{iz} = Tseg(in_window);
            zoom_Vm{iz} = Vmseg(in_window);
        end
    end

    keep = 1:plot_stride:size(Yseg,1);
    if keep(end)~=size(Yseg,1), keep=[keep size(Yseg,1)]; end
    if seg>1, keep=keep(2:end); end  % omit duplicate shared boundary
    T_plot_ms = [T_plot_ms; Tseg(keep)]; %#ok<AGROW>
    Y_plot = [Y_plot; Yseg(keep,:)]; %#ok<AGROW>
    clear Tseg Yseg Vmseg odeHandle;
end

assert(size(Y_plot,2)==52 && numel(T_plot_ms)==size(Y_plot,1));
assert(T_plot_ms(1)==0 && T_plot_ms(end)==dur_total);
assert(all(diff(T_plot_ms)>0),'Trajectory times are not strictly increasing.');
% Save only selected derived signals after calculating them below.

%% Concentrations: amounts divided by CURRENT compartment volumes
% Matches the original driver's state mapping / ion_scale convention.
ion_scale = 1e-12;
Vin=Y_plot(:,14); Ving=Y_plot(:,15); Vout=Y_plot(:,16);
assert(all(Vin>0 & Ving>0 & Vout>0), ...
    'A compartment volume became nonpositive; cannot derive concentrations.');
Vm = Y_plot(:,1)-Y_plot(:,21);
Nai=Y_plot(:,5)./Vin/ion_scale;
Nao=Y_plot(:,7)./Vout/ion_scale;
Ki=Y_plot(:,8)./Vin/ion_scale;
Ko=Y_plot(:,10)./Vout/ion_scale;
Cli=Y_plot(:,11)./Vin/ion_scale;
Clig=Y_plot(:,12)./Ving/ion_scale;
Clo=Y_plot(:,13)./Vout/ion_scale;
Ca_o=Y_plot(:,52);   % Original driver calls state 52 extracellular calcium.

%% SELECTED PLOTS: same variables as the uploaded original driver, plus Cl
% Original: Vm, Ki, Ko, Nai, Nao, Ca_o, Vin, Ving, Vout.
% Added: neuronal Cli, extracellular Clo, and glial Clig.
% Exactly one continuous simulation: alpha7 alternates 1.0 -> 0.5.
t_s = T_plot_ms/1000;
plot_cycle_series(t_s,Vm,'Vm = y(1)-y(21) (mV)', ...
    'Membrane potential | all five LOW-HIGH cycles (250-ms sampled)', ...
    fullfile(root_outdir,'01_Vm_all_cycles.png'),dur_total/1000);

% At 250-ms spacing, a spike can be missed; maximum/minimum preserve it.
plot_cycle_series(T_env_ms/1000,Vm_max,'Vm maximum (mV)', ...
    'Vm maximum in each 250-ms window | spike/burst overview', ...
    fullfile(root_outdir,'01b_Vm_max_envelope.png'),dur_total/1000);
plot_cycle_series(T_env_ms/1000,Vm_min,'Vm minimum (mV)', ...
    'Vm minimum in each 250-ms window', ...
    fullfile(root_outdir,'01c_Vm_min_envelope.png'),dur_total/1000);

names = {'Neuron_volume','Glia_volume','ECS_volume', ...
    'Ki','Ko','Nai','Nao','Ca_o','Cli','Clig','Clo'};
values = {Vin,Ving,Vout,Ki,Ko,Nai,Nao,Ca_o,Cli,Clig,Clo};
labels = {'Neuron volume y(14)','Glia volume y(15)','ECS volume y(16)', ...
    'Intracellular K_i (mM)','Extracellular K_o (mM)', ...
    'Intracellular Na_i (mM)','Extracellular Na_o (mM)', ...
    'Extracellular calcium state y(52) (model units)', ...
    'Intracellular Cl_i (mM)','Glial Cl (mM)','Extracellular Cl_o (mM)'};
for k = 1:numel(names)
    plot_cycle_series(t_s,values{k},labels{k}, ...
        sprintf('%s | LOW: alpha7=1.0 -> HIGH: alpha7=0.5 | 5 cycles',names{k}), ...
        fullfile(root_outdir,sprintf('%02d_%s_all_cycles.png',k+1,names{k})), ...
        dur_total/1000);
end

%% Full-resolution Vm snippets for checking spikes / bursts
for iz=1:numel(zoom_centers_s)
    if isempty(zoom_T{iz}), continue; end
    fig=figure('Visible','off','Color','w','Position',[100 100 1200 500]);
    plot(zoom_T{iz}/1000,zoom_Vm{iz},'LineWidth',0.9);
    xlabel('Time (s)'); ylabel('Vm (mV)'); grid on;
    title(sprintf('Vm at %.0f-%.0f s | full 1-ms output', ...
        zoom_centers_s(iz),zoom_centers_s(iz)+3));
    exportgraphics(fig,fullfile(root_outdir, ...
        sprintf('ZOOM_Vm_%04d_to_%04d_sec.png',zoom_centers_s(iz), ...
        zoom_centers_s(iz)+3)),'Resolution',300);
    close(fig);
end

%% Save ONLY the variables used in the selected plots and the protocol.
save(fullfile(root_outdir,'selected_trajectories.mat'), ...
    'T_plot_ms','Vm','Vin','Ving','Vout','Ki','Ko','Nai','Nao', ...
    'Ca_o','Cli','Clig','Clo','T_env_ms','Vm_min','Vm_max', ...
    'zoom_T','zoom_Vm','protocol','-v7.3');
fprintf('\nCOMPLETE: 5 LOW-HIGH cycles with alpha7 1.0 / 0.5.\n');
fprintf('Selected plots and data saved to: %s\n',root_outdir);

%% Local functions must remain at END of this script.
function plot_cycle_series(t_s,y,ylabel_text,title_text,filename,total_s)
    fig=figure('Visible','off','Color','w','Position',[100 100 1300 500]);
    plot(t_s,y,'LineWidth',1);
    xlabel('Time (s)'); ylabel(ylabel_text,'Interpreter','tex');
    title(title_text,'Interpreter','none');
    xlim([0 total_s]); grid on;
    hold on;
    for sw=400:400:(total_s-400), xline(sw,'--'); end
    exportgraphics(fig,filename,'Resolution',300);
    close(fig);
end

function y0 = make_initial_condition_52()


    y0 = zeros(52,1);


    Vin0  = 2.09;
    Vg0   = 4.1888;
    Vout0 = 0.4;

    scale = 1e-12;


    Nai0  = 9.0;
    Naig0 = 9.0;
    Nao0  = 140.0;


    Ki0   = 113.0;
    Kig0  = 113.0;
    Ko0   = 2.8;


    Cli0  = 9.0;
    Clig0 = 9.0;
    Clo0  = 130.0;


    Vn_ref = -65.6278;
    Vg_ref = -65.3696;
    VBE0   = 0;


    y0(1) = Vn_ref;

    y0(2) = Vn_ref + Vg_ref;


    Vm = Vn_ref;


    alpha_n = ...
        0.01*(Vm+34.0)/(1.0-exp(-0.1*(Vm+34.0)));

    beta_n = ...
        0.125*exp(-(Vm+44.0)/80.0);


    alpha_h = ...
        0.07*exp(-(Vm+44.0)/20.0);

    beta_h = ...
        1.0/(1.0+exp(-0.1*(Vm+14.0)));


    y0(3) = alpha_n/(alpha_n+beta_n);

    y0(4) = alpha_h/(alpha_h+beta_h);


    y0(5)  = Nai0  * Vin0  * scale;
    y0(6)  = Naig0 * Vg0   * scale;
    y0(7)  = Nao0  * Vout0 * scale;


    y0(8)  = Ki0  * Vin0  * scale;
    y0(9)  = Kig0 * Vg0   * scale;
    y0(10) = Ko0  * Vout0 * scale;


    y0(11) = Cli0  * Vin0  * scale;
    y0(12) = Clig0 * Vg0   * scale;
    y0(13) = Clo0  * Vout0 * scale;


    y0(14) = Vin0;
    y0(15) = Vg0;
    y0(16) = Vout0;


    y0(17) = 2.0;
    y0(18) = 20.0;
    y0(19) = 2.0;
    y0(20) = 20.0;


    y0(21) = VBE0;


    y0(22) = 0;
    y0(23) = 0;
    y0(24) = 1;
    y0(25) = 0;


    y0(26) = 1;
    y0(27) = 0;
    y0(28) = 0;
    y0(29) = 0;
    y0(30) = 0;
    y0(31) = 0;


    y0(32) = 1;
    y0(33) = 0;
    y0(34) = 0;


    y0(35) = 0.01932e-4;

    y0(36) = 1.6e-4;

    y0(37) = 3.9e-4;

    y0(38) = 0.30;

    y0(39) = 0.30;


    y0(40) = 0;

    y0(41) = 0;


    y0(42) = 0;

    y0(43) = 0.05;

    y0(44) = 0.10;


    y0(45) = ...
        1/(1 + exp(-(Vm + 52)/7.4));

    y0(46) = ...
        1/(1 + exp((Vm + 80)/5));

    y0(47) = 0.05;

    y0(48) = 0;

    y0(49) = 1;

    y0(50) = 0;

    y0(51) = 0;

    % Extracellular calcium, matching the original model baseline.
    y0(52) = 2.0;

end