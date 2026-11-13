%% DRIVER: HIGH gKs -> LOW gKs, 5 cycles, one full-resolution Vm plot
% This is one continuous trajectory with parameter switches every 400 s.
% Runs with the UNMODIFIED model named below, placed in the same directory.
clear; clc; close all;
modelName = 'NE_GKS_POISSON_STRENGTH_MODEL_GKL_SCALED_CAO_CORRECTED';
assert(exist(modelName,'file')==2, '52-state model .m file not found on MATLAB path.');
A=2;
W_values=3:4;
amp=3;
NE_level=1.0; k_NE_glia=0.5; k_gks_mod=0.5;
poisson_dt_ms=1; dt_out_ms=1; plot_step_ms=250;
chunk_ms=10000;  % restart integrator every 10 s to bound output memory
cycles=5; segment_ms=400000; total_ms=2*cycles*segment_ms;
assert(mod(chunk_ms,plot_step_ms)==0 && mod(segment_ms,chunk_ms)==0);
assert(mod(plot_step_ms,dt_out_ms)==0 && mod(chunk_ms,dt_out_ms)==0);
% Model LOOKUP arrays: a7s=[0 .5 1 .75 .9 1.1 1.25],
% gkks=[.07154 .14 .22062], strength=[.2 1.5 2 3 5].
gks_indices=repmat([3 1],1,cycles);
rate_Hz=repmat([10 20],1,cycles);
strength_indices=repmat([2 4],1,cycles);
strength_values=repmat([1.5 3],1,cycles);
alpha7_indices=repmat([5 3],1,cycles);
alpha7_values=repmat([0.9 1.0],1,cycles);
assert(isequal(alpha7_values, ...
       [0.9 1.0 0.9 1.0 0.9 1.0 0.9 1.0 0.9 1.0]));
for W = W_values
    rng(1,'twister');
    outdir=sprintf('RESULTS_HIGH09_LOW10_REPEAT5_W%d_AMP%d_ONE_FULL_VM',W,amp);
    if ~exist(outdir,'dir'), mkdir(outdir); end
    opts=odeset('RelTol',1e-6,'AbsTol',1e-9, ...
                'MaxStep',min(dt_out_ms,poisson_dt_ms/2));
    % Save a reproducible segment log. No seven-condition sweep.
    segment=(1:2*cycles).';
    start_s=(segment-1)*segment_ms/1000;
    end_s=segment*segment_ms/1000;
    phase=repmat({'HIGH';'LOW'},cycles,1);
    protocol=table(segment,start_s,end_s,phase,gks_indices(:), ...
        rate_Hz(:),strength_indices(:),strength_values(:), ...
        alpha7_indices(:),alpha7_values(:), ...
        'VariableNames',{'Segment','Start_s','End_s','Phase','gKs_Index', ...
        'Frequency_Hz','Strength_Index','Strength','Alpha7_Index','Alpha7'});
    disp(protocol);
    writetable(protocol,fullfile(outdir,'protocol.csv'));
    % Keep only the states needed for volumes and ions at 250-ms intervals.
    % Column order: Na/K/Cl amounts, volumes, Ca_c, Ca_ds, Ca_o.
    cols=[5 7 8 10 11 12 13 14 15 16 36 37 52];
    Tplot=zeros(total_ms/plot_step_ms+1,1);
    Yplot=zeros(numel(Tplot),numel(cols));
    nPlot=0;
    % Store only time and Vm at full 1-ms resolution for the single Vm figure.
    % This reproduces the old plot without retaining all 52 states for 4000 s.
    Tvm_full=zeros(total_ms/dt_out_ms+1,1);
    Vm_full=zeros(numel(Tvm_full),1);
    nVm=0;
    y_current=make_initial_condition_52();
    assert(isequal(size(y_current),[52 1]) && y_current(52)==2);
    % Store every 10-s endpoint with all 52 states. These checkpoints verify that
    % no state is dropped when the solver restarts and permit later restart work.
    nChunks=total_ms/chunk_ms;
    Tcheckpoint=zeros(nChunks+1,1);
    Ycheckpoint=zeros(nChunks+1,52);
    Ycheckpoint(1,:)=y_current.';
    nCheckpoint=1;
    for seg=1:2*cycles
        t0=(seg-1)*segment_ms; t1=seg*segment_ms;
        fprintf('Segment %d/10 | %s | %g-%g s | alpha7=%g | gKs index=%d | %g Hz | strength=%g\n', ...
           seg,phase{seg},t0/1000,t1/1000,alpha7_values(seg), ...
           gks_indices(seg),rate_Hz(seg),strength_values(seg));
        % The model stores a persistent random array. Clear it exactly once per
        % segment, but DO NOT reset the 52 physiological states or reseed rng.
        eval(['clear ' modelName]);
        modelFun=str2func(modelName);
        f=@(t,y) modelFun(t,y,A,W,alpha7_indices(seg),amp, ...
            gks_indices(seg),NE_level,k_NE_glia,k_gks_mod,total_ms, ...
            rate_Hz(seg),poisson_dt_ms,strength_indices(seg));
        % Chunked 1-ms output avoids concatenating >4 million x 52 doubles.
        % Using identical 10-s chunks in both scripts makes comparisons consistent.
        for c0=t0:chunk_ms:t1-chunk_ms
            c1=c0+chunk_ms;
            [tt,yy]=ode15s(f,c0:dt_out_ms:c1,y_current,opts);
            assert(size(yy,2)==52 && all(isfinite(yy(:))), ...
                'Unexpected number of states or nonfinite solution.');
            y_current=yy(end,:).';  % retain ALL 52 states at chunk/segment changes
            vm=yy(:,1)-yy(:,21);
            nCheckpoint=nCheckpoint+1;
            Tcheckpoint(nCheckpoint)=c1;
            Ycheckpoint(nCheckpoint,:)=y_current.';
            if any(any(yy(:,[14 15 16])<=0))
                error('A volume became nonpositive near t=%g s.',c0/1000);
            end
            % Keep the full 1-ms Vm trace. Omit the duplicated first point at
            % every chunk boundary because it is the previous chunk's endpoint.
            jj=1:numel(tt);
            if c0>0, jj=jj(2:end); end
            nFull=numel(jj);
            Tvm_full(nVm+(1:nFull))=tt(jj);
            Vm_full(nVm+(1:nFull))=vm(jj);
            nVm=nVm+nFull;
            % Keep each global 250-ms sample only ONCE, including final endpoint.
            ii=1:(plot_step_ms/dt_out_ms):numel(tt);
            if c0>0, ii=ii(2:end); end
            k=numel(ii);
            Tplot(nPlot+(1:k))=tt(ii);
            Yplot(nPlot+(1:k),:)=yy(ii,cols);
            nPlot=nPlot+k;
        end
    end
    assert(nPlot==numel(Tplot));
    assert(Tplot(1)==0 && Tplot(end)==total_ms && all(diff(Tplot)>0));
    assert(nCheckpoint==numel(Tcheckpoint) && all(diff(Tcheckpoint)>0));
    assert(nVm==numel(Tvm_full) && Tvm_full(1)==0 && ...
           Tvm_full(end)==total_ms && all(diff(Tvm_full)>0), ...
           'The full 1-ms Vm trajectory was not stored correctly.');
    % Derived concentrations: ion amount / CURRENT volume / model's 1e-12 scale.
    s=1e-12; t=Tplot/1000;
    Vin=Yplot(:,8); Ving=Yplot(:,9); Vout=Yplot(:,10);
    Nai=Yplot(:,1)./Vin/s; Nao=Yplot(:,2)./Vout/s;
    Ki=Yplot(:,3)./Vin/s; Ko=Yplot(:,4)./Vout/s;
    Cli=Yplot(:,5)./Vin/s; Clig=Yplot(:,6)./Ving/s;
    Clo=Yplot(:,7)./Vout/s;
    Ca_c=Yplot(:,11); Ca_ds=Yplot(:,12); Ca_o=Yplot(:,13);
    % PLOTTING: exactly one full-resolution Vm PNG plus one PNG per other variable.
    % There is no Vm envelope and there are no separate Vm zoom figures.
    plot_label=sprintf('HIGH: alpha7 = %.1f, 10 Hz, Str 1.5 -> LOW: alpha7 = %.1f, 20 Hz, Str 3.0 | 400-s segments x 5', ...
        alpha7_values(1),alpha7_values(2));
    plot_label=sprintf('%s | W=%d | amp=%d',plot_label,W,amp);
    fig=figure('Visible','off','Color','w','Position',[100 100 1500 550]);
    plot(Tvm_full/1000,Vm_full,'b','LineWidth',0.6);
    xlabel('Time (sec)'); ylabel('Vm = y(1)-y(21) (mV)');
    title([plot_label ' | full Vm (1-ms output)'],'Interpreter','none');
    mark_switches(total_ms/1000); grid on;
    exportgraphics(fig,fullfile(outdir,'01_Vm_full_1ms.png'),'Resolution',300);
    close(fig);
    % Each neuron/glia/ECS volume has its OWN figure and PNG.
    plot_one_signal(t,Vin,'Neuron volume y(14)','Neuron volume', ...
        '02_Neuron_volume.png',outdir,plot_label,total_ms/1000);
    plot_one_signal(t,Ving,'Glia volume y(15)','Glia volume', ...
        '03_Glia_volume.png',outdir,plot_label,total_ms/1000);
    plot_one_signal(t,Vout,'ECS volume y(16)','ECS volume', ...
        '04_ECS_volume.png',outdir,plot_label,total_ms/1000);
    % One standalone figure per ion and compartment; never mix different scales.
    plot_one_signal(t,Ki,'Intracellular K_i (mM)','Intracellular potassium', ...
        '05_Ki.png',outdir,plot_label,total_ms/1000);
    plot_one_signal(t,Ko,'Extracellular K_o (mM)','Extracellular potassium', ...
        '06_Ko.png',outdir,plot_label,total_ms/1000);
    plot_one_signal(t,Nai,'Intracellular Na_i (mM)','Intracellular sodium', ...
        '07_Nai.png',outdir,plot_label,total_ms/1000);
    plot_one_signal(t,Nao,'Extracellular Na_o (mM)','Extracellular sodium', ...
        '08_Nao.png',outdir,plot_label,total_ms/1000);
    plot_one_signal(t,Cli,'Neuronal intracellular Cl_i (mM)','Neuronal chloride', ...
        '09_Cli.png',outdir,plot_label,total_ms/1000);
    plot_one_signal(t,Clig,'Glial intracellular Cl (mM)','Glial chloride', ...
        '10_Clig.png',outdir,plot_label,total_ms/1000);
    plot_one_signal(t,Clo,'Extracellular Cl_o (mM)','Extracellular chloride', ...
        '11_Clo.png',outdir,plot_label,total_ms/1000);
    plot_one_signal(t,Ca_c,'Cytosolic Ca_c (model units)','Cytosolic calcium', ...
        '12_Ca_c.png',outdir,plot_label,total_ms/1000);
    plot_one_signal(t,Ca_ds,'Submembrane Ca_ds (model units)','Submembrane calcium', ...
        '13_Ca_ds.png',outdir,plot_label,total_ms/1000);
    plot_one_signal(t,Ca_o,'Extracellular Ca_o (model units)','Extracellular calcium', ...
        '14_Ca_o.png',outdir,plot_label,total_ms/1000);
    save(fullfile(outdir,'selected_trajectories.mat'),'Tplot','Vin','Ving', ...
        'Vout','Ki','Ko','Nai','Nao','Cli','Clig','Clo','Ca_c','Ca_ds', ...
        'Ca_o','Tvm_full','Vm_full', ...
        'Tcheckpoint','Ycheckpoint','protocol','A','W','amp','-v7.3');
    fprintf(['Complete: %g s; one full-resolution Vm PNG + 13 volume/ion PNGs ' ...
             '+ protocol.csv + selected_trajectories.mat in %s\n'], ...
             total_ms/1000,outdir);
end % W loop
%% Local helper functions must be at the bottom of a MATLAB script.
function plot_one_signal(t_s,signal,ylab,signal_name,filename,outdir,protocol_text,total_s)
    fig=figure('Visible','off','Color','w','Position',[100 100 1500 550]);
    plot(t_s,signal,'b','LineWidth',1.1);
    xlabel('Time (sec)'); ylabel(ylab,'Interpreter','none');
    title([protocol_text ' | ' signal_name],'Interpreter','none');
    mark_switches(total_s); grid on;
    exportgraphics(fig,fullfile(outdir,filename),'Resolution',300);
    close(fig);
end
function mark_switches(total_s)
    xlim([0 total_s]); hold on;
    for sw=400:400:(total_s-400), xline(sw,'--','HandleVisibility','off'); end
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
