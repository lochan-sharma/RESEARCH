%% =========================================================================
% DRIVER FILE: DRIVER_NE_GKS_POISSON_STRENGTH_MODEL_GKL_SCALED_CAO_CORRECTED.m
% MODEL FILE REQUIRED IN THE SAME FOLDER:
%   NE_GKS_POISSON_STRENGTH_MODEL_GKL_SCALED_CAO_CORRECTED.m
%
% The driver calls the model using:
%   modelName = 'NE_GKS_POISSON_STRENGTH_MODEL_GKL_SCALED_CAO_CORRECTED';
% =========================================================================

%% =========================================================================
% DRIVER 1
%
% LOW -> HIGH cycle repeated 5 times
%
% EACH CYCLE:
%   400 sec LOW  gKs | 20 Hz | Poisson strength 3.0
%   400 sec HIGH gKs | 10 Hz | Poisson strength 1.5
%
% Repeated 5 times:
%
%   0-400 sec      LOW  | 20 Hz | 3.0
%   400-800 sec    HIGH | 10 Hz | 1.5
%
%   800-1200 sec   LOW  | 20 Hz | 3.0
%   1200-1600 sec  HIGH | 10 Hz | 1.5
%
%   1600-2000 sec  LOW  | 20 Hz | 3.0
%   2000-2400 sec  HIGH | 10 Hz | 1.5
%
%   2400-2800 sec  LOW  | 20 Hz | 3.0
%   2800-3200 sec  HIGH | 10 Hz | 1.5
%
%   3200-3600 sec  LOW  | 20 Hz | 3.0
%   3600-4000 sec  HIGH | 10 Hz | 1.5
%
% MODEL POISSON STRENGTH ARRAY:
%   [0.2, 1.5, 2.0, 3.0, 5.0]
%
% Therefore:
%   index 2 = strength 1.5
%   index 4 = strength 3.0
%
% TOTAL SIMULATION TIME = 4000 sec
% =========================================================================

clear;
clc;
close all;

rng(1);

%% =========================================================================
% MODEL NAME
% =========================================================================

modelName = 'NE_GKS_POISSON_STRENGTH_MODEL_GKL_SCALED_CAO_CORRECTED';


%% =========================================================================
% FIXED MODEL INDICES
% =========================================================================

A = 2;
W = 2;
amp = 2;      % retained for compatibility with original model


%% =========================================================================
% GKS INDICES
% =========================================================================

low_gks  = 1;
high_gks = 3;


%% =========================================================================
% PROTOCOL
%
% LOW -> HIGH repeated 5 times
% =========================================================================

number_of_cycles = 5;

segment_duration_ms = 400000;     % 400 sec


% 10 total segments:
%
% LOW HIGH LOW HIGH LOW HIGH LOW HIGH LOW HIGH

gks_sequence = repmat([low_gks, high_gks],1,number_of_cycles);


% Frequency follows gKs:
%
% LOW  = 20 Hz
% HIGH = 10 Hz

poisson_rate_sequence_Hz = ...
    repmat([20, 10],1,number_of_cycles);


% Poisson strength index:
%
% index 2 = 1.5
% index 4 = 3.0
%
% LOW  = 3.0
% HIGH = 1.5

poisson_strength_idx_sequence = ...
    repmat([4, 2],1,number_of_cycles);


% Numerical values for labels only

poisson_strength_values_for_display = ...
    repmat([3.0, 1.5],1,number_of_cycles);


% Segment durations

segment_durations_ms = ...
    segment_duration_ms * ones(1,length(gks_sequence));


%% =========================================================================
% TOTAL SIMULATION DURATION
% =========================================================================

dur_total = sum(segment_durations_ms);

% 10 segments x 400 sec = 4000 sec
% dur_total = 4,000,000 ms


%% =========================================================================
% POISSON SETTINGS
% =========================================================================

poisson_dt_ms = 1.0;


%% =========================================================================
% ALPHA7 SWEEP
% =========================================================================

alpha7_indices = 1:7;

alpha7_values = ...
    [0, 0.5, 1, 0.75, 0.9, 1.1, 1.25];


%% =========================================================================
% BASELINE NE
% =========================================================================

NE_level = 1.0;

k_NE_glia = 0.5;


%% =========================================================================
% GKS MODULATION STRENGTH
% =========================================================================
% g_KL_eff = g_KL_vol_base * ...
%            (1 + k_gks_mod*(gkk/0.14 - 1))
%
% 0.0 = no gKs modulation
% 0.5 = half-strength relative modulation (recommended starting test)
% 1.0 = full relative M2 gKs modulation

k_gks_mod = 0.5;


%% =========================================================================
% OUTPUT TIME STEP
% =========================================================================

dt_out = 1;


%% =========================================================================
% SOLVER SETTINGS
% =========================================================================

solver_max_step_ms = ...
    min(dt_out,poisson_dt_ms/2);


opts = odeset( ...
    'RelTol',1e-6, ...
    'AbsTol',1e-9, ...
    'MaxStep',solver_max_step_ms);


%% =========================================================================
% DOWNSAMPLE VOLUME VARIABLES FOR PLOTTING
% =========================================================================

plot_step_ms = 250;

plot_stride = ...
    max(1,round(plot_step_ms/dt_out));


%% =========================================================================
% OUTPUT FOLDER
% =========================================================================

root_outdir = ...
    'RESULTS_CORRECTED_CAO_LOW20Hz_STR3_HIGH10Hz_STR1p5_400sec_REPEAT5';


if ~exist(root_outdir,'dir')

    mkdir(root_outdir);

end


%% =========================================================================
% STORAGE FOR COMBINED VOLUME PLOTS
% =========================================================================

Tvol_all = cell(length(alpha7_indices),1);

Vin_all = cell(length(alpha7_indices),1);

Ving_all = cell(length(alpha7_indices),1);

Vout_all = cell(length(alpha7_indices),1);


%% =========================================================================
% ALPHA7 LOOP
% =========================================================================

for i = 1:length(alpha7_indices)


    al7 = alpha7_indices(i);


    fprintf('\n');
    fprintf('============================================================\n');

    fprintf(['Running LOW-HIGH 400-sec cycle repeated 5 times | ' ...
             'alpha7 index = %d | alpha7 = %.2f\n'], ...
             al7,alpha7_values(i));

    fprintf('Total simulation = 4000 sec\n');

    fprintf('============================================================\n');


    %% Same random starting point for every alpha7 comparison

    rng(1);


    %% Reset model before independent alpha7 run

    eval(sprintf('clear %s',modelName));

    modelFun = str2func(modelName);


    %% Initial condition

    y_current = make_initial_condition_52();


    %% Storage for complete trajectory

    T = [];

    Y = [];


    current_start_ms = 0;


    %% =====================================================================
    % SEGMENT LOOP
    % ======================================================================

    for seg = 1:length(gks_sequence)


        %% Current parameters

        gks = gks_sequence(seg);

        poisson_rate_Hz = ...
            poisson_rate_sequence_Hz(seg);

        poisson_strength_idx = ...
            poisson_strength_idx_sequence(seg);

        poisson_strength_display = ...
            poisson_strength_values_for_display(seg);


        %% Segment times

        segment_start_ms = current_start_ms;

        segment_end_ms = ...
            segment_start_ms + segment_durations_ms(seg);


        tspan_segment = ...
            segment_start_ms:dt_out:segment_end_ms;


        %% Determine name

        if gks == high_gks

            state_name = 'HIGH';

        else

            state_name = 'LOW';

        end


        %% Display progress

        fprintf('\n');

        fprintf('Segment %d of %d\n', ...
            seg,length(gks_sequence));

        fprintf('  Time = %.0f to %.0f sec\n', ...
            segment_start_ms/1000, ...
            segment_end_ms/1000);

        fprintf('  gKs = %s | index %d\n', ...
            state_name,gks);

        fprintf('  Poisson frequency = %.0f Hz\n', ...
            poisson_rate_Hz);

        fprintf('  Poisson strength = %.1f | index %d\n', ...
            poisson_strength_display, ...
            poisson_strength_idx);


        %% ================================================================
        % RESET PERSISTENT POISSON ARRAY
        %
        % Frequency changes between HIGH and LOW segments.
        % The physiological state y_current is NOT reset.
        % ================================================================

        eval(sprintf('clear %s',modelName));

        modelFun = str2func(modelName);


        %% ================================================================
        % MODEL HANDLE
        % ================================================================

        odeHandle = @(t,y) modelFun( ...
            t, ...
            y, ...
            A, ...
            W, ...
            al7, ...
            amp, ...
            gks, ...
            NE_level, ...
            k_NE_glia, ...
            k_gks_mod, ...
            dur_total, ...
            poisson_rate_Hz, ...
            poisson_dt_ms, ...
            poisson_strength_idx);


        %% ================================================================
        % RUN CURRENT 400-SEC SEGMENT
        % ================================================================

        [Tseg,Yseg] = ode15s( ...
            odeHandle, ...
            tspan_segment, ...
            y_current, ...
            opts);


        %% ================================================================
        % CONTINUE FROM FINAL STATE
        % ================================================================

        y_current = Yseg(end,:).';


        %% ================================================================
        % REMOVE DUPLICATE BOUNDARY
        % ================================================================

        if seg > 1

            Tseg = Tseg(2:end);

            Yseg = Yseg(2:end,:);

        end


        %% ================================================================
        % JOIN SEGMENTS
        % ================================================================

        T = [T; Tseg]; %#ok<AGROW>

        Y = [Y; Yseg]; %#ok<AGROW>


        current_start_ms = segment_end_ms;

    end


    %% =========================================================================
    % EXTRACT VARIABLES
    % =========================================================================

    Vn = Y(:,1);

    Vref = Y(:,21);

    Vm = Vn - Vref;


    Vin = Y(:,14);

    Ving = Y(:,15);

    Vout = Y(:,16);

    % New extracellular calcium state
    Ca_o = Y(:,52);


    %% =========================================================================
    % ION CONCENTRATIONS FROM MODEL AMOUNTS / CURRENT COMPARTMENT VOLUMES
    %
    % y(5)  = neuronal intracellular Na amount
    % y(7)  = extracellular Na amount
    % y(8)  = neuronal intracellular K amount
    % y(10) = extracellular K amount
    % y(14) = neuronal volume Vin
    % y(16) = ECS volume Vout
    %
    % The model uses coeff(46) = 1e-12, so the same conversion is used here.
    % =========================================================================

    ion_scale = 1e-12;

    Nai = Y(:,5)  ./ Y(:,14) ./ ion_scale;
    Nao = Y(:,7)  ./ Y(:,16) ./ ion_scale;

    Ki  = Y(:,8)  ./ Y(:,14) ./ ion_scale;
    Ko  = Y(:,10) ./ Y(:,16) ./ ion_scale;


    %% =========================================================================
    % STORE VOLUMES
    % =========================================================================

    idx_plot = 1:plot_stride:length(T);


    Tvol_all{i} = T(idx_plot);

    Vin_all{i} = Vin(idx_plot);

    Ving_all{i} = Ving(idx_plot);

    Vout_all{i} = Vout(idx_plot);


    alpha_tag = strrep( ...
        sprintf('alpha7_%0.2f',alpha7_values(i)), ...
        '.', ...
        'p');


    %% =========================================================================
    % FOUR ION-CONCENTRATION PLOTS
    %
    % Each alpha7 value gets its OWN PNG for each ion concentration.
    % 7 alpha7 values -> 28 ion PNG files total.
    % =========================================================================

    % 1) NEURONAL INTRACELLULAR K CONCENTRATION: Ki
    fig = figure('Visible','off','Color','w','Position',[100 100 1200 500]);
    plot(T/1000,Ki,'LineWidth',1.1);
    xlabel('Time (sec)');
    ylabel('Intracellular K concentration K_i (mM)');
    title(sprintf(['LOW 20 Hz Str 3.0 -> HIGH 10 Hz Str 1.5 | 400-sec segments repeated 5 times | alpha7 = %.2f | K_i'], ...
        alpha7_values(i)),'Interpreter','none');
    xlim([0 4000]);
    for switch_time = 400:400:3600
        xline(switch_time,'--');
    end
    grid on;
    png_name = sprintf('05_Ki_LOW_HIGH_REPEAT5_%s.png',alpha_tag);
    exportgraphics(fig,fullfile(root_outdir,png_name),'Resolution',300);
    close(fig);

    % 2) EXTRACELLULAR K CONCENTRATION: Ko
    fig = figure('Visible','off','Color','w','Position',[100 100 1200 500]);
    plot(T/1000,Ko,'LineWidth',1.1);
    xlabel('Time (sec)');
    ylabel('Extracellular K concentration K_o (mM)');
    title(sprintf(['LOW 20 Hz Str 3.0 -> HIGH 10 Hz Str 1.5 | 400-sec segments repeated 5 times | alpha7 = %.2f | K_o'], ...
        alpha7_values(i)),'Interpreter','none');
    xlim([0 4000]);
    for switch_time = 400:400:3600
        xline(switch_time,'--');
    end
    grid on;
    png_name = sprintf('06_Ko_LOW_HIGH_REPEAT5_%s.png',alpha_tag);
    exportgraphics(fig,fullfile(root_outdir,png_name),'Resolution',300);
    close(fig);

    % 3) NEURONAL INTRACELLULAR Na CONCENTRATION: Nai
    fig = figure('Visible','off','Color','w','Position',[100 100 1200 500]);
    plot(T/1000,Nai,'LineWidth',1.1);
    xlabel('Time (sec)');
    ylabel('Intracellular Na concentration Na_i (mM)');
    title(sprintf(['LOW 20 Hz Str 3.0 -> HIGH 10 Hz Str 1.5 | 400-sec segments repeated 5 times | alpha7 = %.2f | Na_i'], ...
        alpha7_values(i)),'Interpreter','none');
    xlim([0 4000]);
    for switch_time = 400:400:3600
        xline(switch_time,'--');
    end
    grid on;
    png_name = sprintf('07_Nai_LOW_HIGH_REPEAT5_%s.png',alpha_tag);
    exportgraphics(fig,fullfile(root_outdir,png_name),'Resolution',300);
    close(fig);

    % 4) EXTRACELLULAR Na CONCENTRATION: Nao
    fig = figure('Visible','off','Color','w','Position',[100 100 1200 500]);
    plot(T/1000,Nao,'LineWidth',1.1);
    xlabel('Time (sec)');
    ylabel('Extracellular Na concentration Na_o (mM)');
    title(sprintf(['LOW 20 Hz Str 3.0 -> HIGH 10 Hz Str 1.5 | 400-sec segments repeated 5 times | alpha7 = %.2f | Na_o'], ...
        alpha7_values(i)),'Interpreter','none');
    xlim([0 4000]);
    for switch_time = 400:400:3600
        xline(switch_time,'--');
    end
    grid on;
    png_name = sprintf('08_Nao_LOW_HIGH_REPEAT5_%s.png',alpha_tag);
    exportgraphics(fig,fullfile(root_outdir,png_name),'Resolution',300);
    close(fig);


    %% =========================================================================
    % EXTRACELLULAR CALCIUM PLOT: Ca_o = y(52)
    % =========================================================================

    fig = figure('Visible','off','Color','w','Position',[100 100 1200 500]);
    plot(T/1000,Ca_o,'LineWidth',1.1);
    xlabel('Time (sec)');
    ylabel('Extracellular Ca concentration Ca_o');
    title(sprintf(['LOW 20 Hz Str 3.0 -> HIGH 10 Hz Str 1.5 | 400-sec segments repeated 5 times | alpha7 = %.2f | Ca_o'], ...
        alpha7_values(i)),'Interpreter','none');
    xlim([0 4000]);
    for switch_time = 400:400:3600
        xline(switch_time,'--');
    end
    grid on;
    png_name = sprintf('09_Ca_o_LOW_HIGH_REPEAT5_%s.png',alpha_tag);
    exportgraphics(fig,fullfile(root_outdir,png_name),'Resolution',300);
    close(fig);


    %% =========================================================================
    % FULL Vm PLOT: 0-1000 sec
    % =========================================================================

    fig = figure( ...
        'Visible','off', ...
        'Color','w', ...
        'Position',[100 100 1200 500]);


    plot(T/1000,Vm,'LineWidth',1.1);


    xlabel('Time (sec)');

    ylabel('Vm = y(1)-y(21) (mV)');


    title(sprintf([ ...
        'LOW 20 Hz Str 3.0 -> HIGH 10 Hz Str 1.5 | ' ...
        '400-sec segments repeated 5 times | alpha7 = %.2f'], ...
        alpha7_values(i)), ...
        'Interpreter','none');


    xlim([0 4000]);

    grid on;


    %% Mark every switching point

    for switch_time = 400:400:3600

        xline(switch_time,'--');

    end


    png_name = sprintf( ...
        '01_Vm_LOW_HIGH_REPEAT5_%s.png', ...
        alpha_tag);


    exportgraphics( ...
        fig, ...
        fullfile(root_outdir,png_name), ...
        'Resolution',300);


    close(fig);


    clear T Y Tseg Yseg ...
          Vn Vref Vm ...
          Vin Ving Vout Ca_o ...
          Ki Ko Nai Nao ion_scale ...
          y_current;

end


%% =========================================================================
% COMBINED VOLUME PLOTS
% =========================================================================

mainTitle = [ ...
    'LOW gKs 20 Hz strength 3.0 -> ' ...
    'HIGH gKs 10 Hz strength 1.5 | ' ...
    '400-sec segments repeated 5 times'];


%% Neuron volume

plot_one_volume( ...
    Tvol_all, ...
    Vin_all, ...
    alpha7_values, ...
    root_outdir, ...
    mainTitle, ...
    'Neuron_volume', ...
    'Neuron volume y(14)', ...
    2);


%% Glia volume

plot_one_volume( ...
    Tvol_all, ...
    Ving_all, ...
    alpha7_values, ...
    root_outdir, ...
    mainTitle, ...
    'Glia_volume', ...
    'Glia volume y(15)', ...
    3);


%% ECS volume

plot_one_volume( ...
    Tvol_all, ...
    Vout_all, ...
    alpha7_values, ...
    root_outdir, ...
    mainTitle, ...
    'ECS_volume', ...
    'ECS volume y(16)', ...
    4);


fprintf('\n');
fprintf('============================================================\n');
fprintf('DRIVER 1 COMPLETE\n');
fprintf('LOW -> HIGH 400-sec sequence repeated 5 times\n');
fprintf('Total simulation time = 4000 sec\n');
fprintf('Results saved in:\n%s\n',root_outdir);
fprintf('============================================================\n');


%% =========================================================================
% LOCAL VOLUME-PLOTTING FUNCTION
% =========================================================================

function plot_one_volume( ...
    T_all, ...
    Y_all, ...
    alpha7_values, ...
    outdir, ...
    mainTitle, ...
    data_name, ...
    y_label, ...
    file_number)


    fig = figure( ...
        'Visible','off', ...
        'Color','w', ...
        'Position',[100 100 1200 500]);


    hold on;


    for i = 1:length(alpha7_values)

        plot( ...
            T_all{i}/1000, ...
            Y_all{i}, ...
            'LineWidth',1.2);

    end


    xlabel('Time (sec)');

    ylabel(y_label);


    title( ...
        [mainTitle ' | ' data_name], ...
        'Interpreter','none');


    legend( ...
        compose('\\alpha7 = %.2f',alpha7_values), ...
        'Location','best');


    xlim([0 4000]);


    %% Mark all 400-sec transitions

    for switch_time = 400:400:3600

        xline(switch_time,'--');

    end


    grid on;


    png_name = sprintf( ...
        '%02d_%s_vs_time.png', ...
        file_number, ...
        data_name);


    exportgraphics( ...
        fig, ...
        fullfile(outdir,png_name), ...
        'Resolution',300);


    close(fig);

end


%% =========================================================================
% LOCAL INITIAL-CONDITION FUNCTION
% =========================================================================

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