function dy = NE_GKS_POISSON_STRENGTH_MODEL_GKL_SCALED_CAO_CORRECTED(t,y,A,W,al7,amp,gks,NE_level,k_NE_glia,k_gks_mod,t_end_ms,poisson_rate_Hz,poisson_dt_ms,poisson_strength_idx)
% NE-UPDATED VERSION:
%   Adds norepinephrine (NE) as a dimensionless glial scaling factor.
%   NE selectively scales glial KIR/inward-rectifier K current, glial
%   Na/K pump, and glial NKCC/NKCC1-like cotransport.
%
%   Important update in this version:
%   NE scaling has been removed from glial Na leak and glial Cl leak.
%
% Backward-compatible calling:
%   If NE_level and k_NE_glia are not provided, the model defaults to
%   NE_level = 1.0 and k_NE_glia = 0.5.
%
%   GKS MODULATION OF THE VOLUME-MODEL K LEAK:
%   The duplicated full M2 I_K remains disabled.
%   Instead, the M2-selected gkk value modulates the ORIGINAL volume-model
%   neuronal K leak conductance by a dimensionless relative scaling:
%
%       g_KL_eff = g_KL_vol_base * ...
%                  (1 + k_gks_mod*(gkk/0.14 - 1))
%
%   k_gks_mod = 0   -> no gKs modulation; original volume-model leak retained
%   k_gks_mod = 0.5 -> half-strength relative modulation (recommended start)
%   k_gks_mod = 1   -> full relative modulation from the M2 gKs ratios
%
%   This gives one main neuronal K current:
%       gks -> g_KL_eff -> Ik_vol -> NinK/NoutK -> Ko -> osmotic volume
%   There is NO separate I_K_gks_vol term in the K amount equations.
%
%   FIXED-BIN POISSON STIMULATION UPDATE:
%   A complete binary Poisson input array is generated ONCE inside this
%   model and stored as a persistent variable. Time is divided into equal
%   bins of width poisson_dt_ms. Each bin contains either:
%       0 = no stimulation
%       1 = stimulation during that one bin
%
%   The array is generated using the professor's threshold method:
%       p_event = poisson_rate_Hz * poisson_dt_ms / 1000
%       random value > (1-p_event)  -> 1
%       otherwise                   -> 0
%
%   The adaptive ODE solver never draws new random numbers. It only reads the
%   already-generated 0/1 value for the current fixed time bin. The old
%   stml-long pulse trigger is not used by this Poisson input.
%
%   SEPARATE POISSON STRENGTH UPDATE:
%   poisson_strength_idx selects the amplitude attached to each binary
%   Poisson event. This is independent of the original amp/ampst parameter.
%   Five available Poisson pulse strengths are defined below as:
%       [1.75, 3.5, 5.0, 7.0, 9.0]
%   The applied Poisson current is:
%       I_app = poisson_strength * poisson_event
%
% Example solver call:
%   NE_level  = 1.0;   % 0.2 low/sleep-like, 1 baseline, 2 high/wake-like
%   k_NE_glia = 0.5;   % NE sensitivity
%   poisson_strength_idx = 2; % selects Poisson pulse strength 1.0
%   ode15s(@(t,y) NE_GKS_POISSON_STRENGTH_MODEL_GKL_SCALED_CAO_CORRECTED(t,y,A,W,al7,amp,gks,NE_level,k_NE_glia,k_gks_mod,t_end_ms,poisson_rate_Hz,poisson_dt_ms,poisson_strength_idx), tspan, y0);
% CLEANED VERSION: duplicated M2 Na/K/Cl/pump/KCC2 terms are not allowed
% to drive dy by default. See CLEANED MERGE CONTROL block below.
% =========================================================================
% ODE45HH427_MERGED_KEEPALL_HEAVYCOMMENTS
% =========================================================================
% STATE VARIABLE MAP
% =========================================================================
% -------------------------
% VOLUME-BACKBONE STATES
% -------------------------
% y(1)  = Vn_mem   % neuron membrane voltage-like variable
% y(2)  = Vg_mem   % glial membrane voltage-like variable
% y(3)  = n_vol    % HH-type K activation gate from volume code
% y(4)  = h_vol    % HH-type Na inactivation gate from volume code
%
% y(5)  = NinNa    % neuronal intracellular Na amount
% y(6)  = NinNag   % glial intracellular Na amount
% y(7)  = NoutNa   % extracellular Na amount
%
% y(8)  = NinK     % neuronal intracellular K amount
% y(9)  = NinKg    % glial intracellular K amount
% y(10) = NoutK    % extracellular K amount
%
% y(11) = NinCl    % neuronal intracellular Cl amount
% y(12) = NinClg   % glial intracellular Cl amount
% y(13) = NoutCl   % extracellular Cl amount
%
% y(14) = Vin      % neuronal intracellular volume
% y(15) = Ving     % glial intracellular volume
% y(16) = Vout     % extracellular volume (ECS volume)
%
% y(17) = pyrn     % neuronal pyruvate-related metabolic state
% y(18) = ATPn     % neuronal ATP-related metabolic state
% y(19) = pyrg     % glial pyruvate-related metabolic state
% y(20) = ATPg     % glial ATP-related metabolic state
%
% y(21) = VBE      % extracellular/bath electrical variable
%
% -------------------------
% ORIGINAL M2-ONLY STATES
% -------------------------
% y(22) = Po1      % RyR open state 1
% y(23) = Po2      % RyR open state 2
% y(24) = Pc1      % RyR closed state 1
% y(25) = Pc2      % RyR closed state 2
%
% y(26) = x000     % IP3 receptor Markov/state variable
% y(27) = x100     % IP3 receptor Markov/state variable
% y(28) = x110     % IP3 receptor Markov/state variable
% y(29) = x101     % IP3 receptor Markov/state variable
% y(30) = x111     % IP3 receptor Markov/state variable
% y(31) = x001     % IP3 receptor Markov/state variable
%
% y(32) = C        % alpha7 receptor closed state
% y(33) = O        % alpha7 receptor open state
% y(34) = D        % alpha7 receptor desensitized/inactive state
%
% y(35) = IP3      % intracellular IP3 concentration/state
%
% y(36) = Ca_c     % cytosolic Ca concentration
% y(37) = Ca_ds    % submembrane / dense-shell Ca concentration
% y(38) = Ca_jer   % junctional ER Ca concentration
% y(39) = Ca_ner   % network ER Ca concentration
%
% y(40) = delay    % original M2 timing/delay state
% y(41) = check    % original M2 diagnostic/helper state
%
% y(42) = M2       % M2 receptor activation state
% y(43) = ACh      % acetylcholine state
% y(44) = BF       % basal forebrain input/activity state
%
% y(45) = m_na     % fast Na activation gate from original M2 code
% y(46) = h_na     % fast Na inactivation gate from original M2 code
% y(47) = j_na     % fast K / repolarization-related gate from original M2 code
% y(48) = m_ca     % T-type Ca activation gate
% y(49) = h_ca     % T-type Ca inactivation gate
% y(50) = m_KCa    % Ca-activated K activation gate
% y(51) = m_CAN    % CAN channel activation gate
% y(52) = Ca_o     % extracellular calcium concentration
%
% =========================================================================
% QUICK INTERPRETATION OF MODEL LAYERS
% =========================================================================
% 1) Physical / compartment backbone:
%    y(1:21)
%    -> neuron-glia-ECS ion amounts, volumes, glial voltage, metabolism
%
% 2) Intracellular signaling / receptor block from M2 model:
%    y(22:39)
%    -> RyR, IP3R, alpha7, IP3, calcium stores
%
% 3) Cholinergic control block:
%    y(42:44)
%    -> M2, ACh, BF
%
% 4) Original M2 fast gating block:
%    y(45:51)
%    -> Na/Ca/KCa/CAN gating variables
%
% NOTE:
% This file is intentionally a "keep-everything" scaffold.
% Some blocks may later turn out to be redundant, but they are preserved
% here for discussion and future model refinement.
% =========================================================================
%
% LITERAL FIRST-PASS MERGE OF:
%   (1) volume / neuron-glia-ECS model
%   (2) M2 / BF / ACh / alpha7 / IP3 / RyR / Ca model
%
% -------------------------------------------------------------------------
% IMPORTANT
% -------------------------------------------------------------------------
% This file is intentionally written as a "keep everything" scaffold.
% It is NOT a reduced final model.
%
% Goal:
%   - keep the original code blocks
%   - place the M2 model inside the volume-model scaffold
%   - avoid inventing extra biology
%   - use comments to mark places that may need discussion later
%
% -------------------------------------------------------------------------
% STATE ORDER
% -------------------------------------------------------------------------
% Volume-code states:
%   y(1)  = neuron membrane voltage variable
%   y(2)  = glial membrane voltage variable
%   y(3)  = volume-code HH n gate
%   y(4)  = volume-code HH h gate
%   y(5)  = neuronal intracellular Na amount
%   y(6)  = glial intracellular Na amount
%   y(7)  = extracellular Na amount
%   y(8)  = neuronal intracellular K amount
%   y(9)  = glial intracellular K amount
%   y(10) = extracellular K amount
%   y(11) = neuronal intracellular Cl amount
%   y(12) = glial intracellular Cl amount
%   y(13) = extracellular Cl amount
%   y(14) = neuron intracellular volume
%   y(15) = glial intracellular volume
%   y(16) = extracellular volume
%   y(17) = neuronal pyruvate
%   y(18) = neuronal ATP
%   y(19) = glial pyruvate
%   y(20) = glial ATP
%   y(21) = ECS-bath voltage variable
%
% Added M2-only states:
%   y(22) = Po1
%   y(23) = Po2
%   y(24) = Pc1
%   y(25) = Pc2
%   y(26) = x000
%   y(27) = x100
%   y(28) = x110
%   y(29) = x101
%   y(30) = x111
%   y(31) = x001
%   y(32) = C
%   y(33) = O
%   y(34) = D
%   y(35) = IP3
%   y(36) = Ca_c
%   y(37) = Ca_ds
%   y(38) = Ca_jer
%   y(39) = Ca_ner
%   y(40) = delay
%   y(41) = check
%   y(42) = M2
%   y(43) = ACh
%   y(44) = BF
%   y(45) = m_na
%   y(46) = h_na
%   y(47) = j_na
%   y(48) = m_ca
%   y(49) = h_ca
%   y(50) = m_KCa
%   y(51) = m_CAN
%   y(52) = Ca_o
%
% TOTAL STATES = 52
%
% -------------------------------------------------------------------------
% NOTE ON SHARED PHYSICAL VARIABLES
% -------------------------------------------------------------------------
% The original M2 code had its own Vm and ion concentration states.
% In this merged scaffold, the SHARED physical variables come from the
% volume backbone:
%
%   Vm  <- from volume-code voltage
%   Ko, Ki, Nai, Nao, Cli, Clo <- from volume ion amounts / volumes
%
% This is the minimum necessary merge step.
% =========================================================================

% USAGE NOTE
%   This function is a keep-all merge built on the volume backbone.
%   Shared neuron voltage and shared Na/K/Cl concentrations are taken
%   from the volume model, while the alpha7/M2/BF/ACh/Ca machinery is
%   appended on top of that backbone.
%
%   Main design choice in this merge:
%   - Keep one shared neuron membrane potential: Vm = Vn_mem - VBE
%   - Keep one shared Na/K/Cl state set: derived from amount/volume vars
%   - Preserve BOTH electrophysiology blocks in the neuron equation
%   - Preserve BOTH stimulation logics for later comparison
%
dy = zeros(52,1);

%% =========================================================================
% ORIGINAL M2 PARAMETER-SELECTION BLOCK
% These arrays choose stimulation duration, delay, alpha7 scaling, leak K,
% and amplitude based on the input indices A, W, al7, amp, gks.
% =========================================================================
stms  = [20, 80, 320];
delat = [20, 40, 80, 160];
a7s   = [0, .5, 1, .75, .9, 1.1, 1.25];
gkks  = [.07154, 0.14, 0.22062];
amps  = [1.75, 3.5, 5];

% Original model stimulation/feedback amplitude. It is intentionally kept
% separate from the new Poisson pulse amplitude defined below.
stml    = stms(A);
dtarget = delat(W);
a7on    = a7s(al7);
ampst   = amps(amp); %#ok<NASGU> % retained for original-model compatibility
gkk     = gkks(gks);

%% =========================================================================
% NE MODULATION BLOCK
% =========================================================================
% NE_level is a dimensionless norepinephrine level:
%   NE_level = 0.2  -> low NE / sleep-like
%   NE_level = 1.0  -> baseline NE; original model is unchanged
%   NE_level = 2.0  -> high NE / wake-like or high-arousal condition
%
% k_NE_glia controls how strongly NE scales the selected glial mechanisms.
% The scale is written so that NE_level = 1 always gives scale = 1.
%
% Selected NE targets in this version:
%   1) glial KIR / inward-rectifier K current: Ikg
%   2) glial Na/K pump current:                itildepmpg
%   3) glial NKCC/NKCC1-like term:             Inkccg
%
% Important: NE is NOT applied to the whole glial balance equation.
% It only scales the agreed glial current/transport terms.
% NE is NOT applied to glial Na leak Inag.
% NE is NOT applied to glial Cl leak Iclg.

if nargin < 8 || isempty(NE_level)
    NE_level = 1.0;
end

if nargin < 9 || isempty(k_NE_glia)
    k_NE_glia = 0.5;
end

% -------------------------------------------------------------------------
% GKS MODULATION STRENGTH
% -------------------------------------------------------------------------
% k_gks_mod controls how strongly the M2 gKs ratios modulate the ORIGINAL
% volume-model neuronal K leak conductance.
%
%   k_gks_mod = 0.0 -> original volume leak for every gKs condition
%   k_gks_mod = 0.5 -> half of the relative M2 modulation (recommended start)
%   k_gks_mod = 1.0 -> full relative M2 modulation
%
% The middle M2 value gkk = 0.14 is the reference, so when gkk = 0.14:
%       g_KL_eff = g_KL_vol_base
%
if nargin < 10 || isempty(k_gks_mod)
    k_gks_mod = 0.5;
end

if ~isscalar(k_gks_mod) || k_gks_mod < 0 || k_gks_mod > 1
    error('k_gks_mod must be a scalar between 0 and 1.');
end

% -------------------------------------------------------------------------
% FIXED-BIN POISSON INPUT PARAMETERS
% -------------------------------------------------------------------------
% t_end_ms controls only the length of the internally generated 0/1 array.
% poisson_rate_Hz is the desired average number of events per second.
% poisson_dt_ms is the fixed width of one binary stimulation bin.
if nargin < 11 || isempty(t_end_ms)
    t_end_ms = 300000;       % fallback only; driver should pass the duration
end
if nargin < 12 || isempty(poisson_rate_Hz)
    poisson_rate_Hz = 20;    % 20 events/s on average
end
if nargin < 13 || isempty(poisson_dt_ms)
    poisson_dt_ms = 1.0;     % 1-ms Poisson bins
end

% -------------------------------------------------------------------------
% SEPARATE POISSON PULSE-STRENGTH PARAMETER
% -------------------------------------------------------------------------
% This parameter controls only the amplitude of each Poisson-driven current
% pulse. It does not alter the original amp/ampst parameter, the Poisson
% event rate, or the Poisson bin width.
%
% Index mapping:
%   1 -> 1.75
%   2 -> 3.50
%   3 -> 5.00
%   4 -> 7.00
%   5 -> 9.00
poisson_strength_values = [0.2, 1.5, 2.0, 3, 5.0];

if nargin < 14 || isempty(poisson_strength_idx)
    poisson_strength_idx = 2; % default Poisson pulse strength = 1.0
end

if ~isscalar(poisson_strength_idx) || ...
   poisson_strength_idx ~= floor(poisson_strength_idx) || ...
   poisson_strength_idx < 1 || ...
   poisson_strength_idx > numel(poisson_strength_values)
    error('poisson_strength_idx must be an integer from 1 to 5.');
end

poisson_strength = poisson_strength_values(poisson_strength_idx);

if t_end_ms <= 0
    error('t_end_ms must be positive.');
end
if poisson_rate_Hz < 0
    error('poisson_rate_Hz must be nonnegative.');
end
if poisson_dt_ms <= 0
    error('poisson_dt_ms must be positive.');
end

% Small-bin Poisson approximation used in the professor's notes.
% Example: 20 Hz with 1-ms bins gives p_event = 0.02 and threshold = 0.98.
p_event_poisson = poisson_rate_Hz * poisson_dt_ms / 1000;
if p_event_poisson > 1
    error(['poisson_rate_Hz*poisson_dt_ms/1000 must be <= 1. ' ...
           'Use a smaller Poisson bin width.']);
end

NE_glia_scale = 1 + k_NE_glia*(NE_level - 1);

% Safety clamp: prevents negative conductance/transport strength if a very
% low NE_level or large k_NE_glia is tested.
NE_glia_scale = max(0.05, NE_glia_scale);

%% =========================================================================
% ORIGINAL VOLUME-CODE PARAMETERS
% coeff(i) comments added for documentation only
% =========================================================================
coeff(1) = 6;        % neuronal delayed-rectifier K conductance g_k
coeff(2) = 20.5;     % neuronal fast Na conductance g_Na
coeff(3) = 25;       % legacy parameter from original volume code
coeff(4) = 10;       % glial inward-rectifier/KIR shape factor
coeff(5) = 5;        % legacy / not clearly used in original volume code
coeff(6) = -65.6278; % reference/resting neuronal membrane voltage
coeff(7) = -65.3696; % reference/resting glial membrane voltage
coeff(8) = 0.02;     % neuronal Na leak conductance g_naL
coeff(9) = 0.095;    % glial Na leak conductance g_naLg
coeff(10) = 0.25;    % neuronal Cl leak conductance g_clL
coeff(11) = 0.7;     % glial Cl leak conductance g_clLg
coeff(12) = 1e-6;    % current-to-mole scaling prefactor
coeff(13) = 20;      % pump maximum / ATP reference scale pmax
coeff(14) = 1;       % pump sigmoid slope pslop
coeff(15) = 0;       % ATP dependence switch used in original code
coeff(16) = 10000;   % burst repetition period
coeff(17) = 2000;    % burst duration
coeff(18) = 10;      % stimulation start time
coeff(19) = 200500;  % stimulation end time
coeff(20) = 100;     % pulse repetition period inside burst
coeff(21) = 40;      % pulse duration
coeff(22) = -0.00;   % FUJ term from original code (volume/bath related forcing term)
coeff(23) = 0;       % JRAT / coupling term used with FUJ and volume flux
coeff(24) = 0;       % legacy parameter
coeff(25) = 28.3;    % glial Na pump half-saturation term napg
coeff(26) = 1.0;     % glial K pump half-saturation term kpg
coeff(27) = 28.3;    % neuronal Na pump half-saturation term nap
coeff(28) = 1.0;     % neuronal K pump half-saturation term kp
coeff(29) = 147.8;   % immobile osmolyte concentration constant for neuron
coeff(30) = 147.8;   % immobile osmolyte concentration constant for glia
coeff(31) = 6.0;     % immobile osmolyte concentration constant for ECS/bath
coeff(32) = 0;       % optional stimulation current term in voltage equation
coeff(33) = 0*.4;    % optional Na stimulation-related conductance term
coeff(34) = .0;      % optional Cl stimulation-related conductance term
coeff(35) = 0*.4;    % optional K stimulation-related conductance term
coeff(36) = 0;       % legacy parameter
coeff(37) = 95000;   % slow relaxation timescale for volume restoration
coeff(38) = 0;       % legacy parameter
coeff(39) = 110;     % basal parameter from original code
coeff(40) = 1;       % exponent used in volume-related scaling mute2
coeff(41) = 0.4;     % reference extracellular volume
coeff(42) = 2.09;    % reference neuronal intracellular volume
coeff(43) = 4.1888;  % reference glial intracellular volume
coeff(44) = 6.02e23; % Avogadro number
coeff(45) = 1.66e-19;% electronic charge scale used in conversion
coeff(46) = 1e-12;   % pico-volume / pico-scaling factor
coeff(47) = 3.0;     % HH gating speed factor phi
coeff(48) = 1.0;     % membrane capacitance Cm
coeff(49) = 1e13;    % kappa used in bath/ECS capacitance-like term
coeff(50) = 8.85e-12;% permittivity epsilon
coeff(51) = 3.33;    % rho used in drift calculation
coeff(52) = 3;       % spatial scale / distance factor d11
coeff(53) = 5.5e-4;  % effective cell radius
coeff(54) = 2.8;     % bath potassium concentration
coeff(55) = 140.0;   % bath sodium concentration
coeff(56) = 130.0;   % bath chloride concentration
coeff(57) = 9;       % baseline neuronal intracellular chloride
coeff(58) = 2.8;     % baseline extracellular potassium
coeff(59) = 140;     % baseline extracellular sodium
coeff(60) = 130;     % baseline extracellular chloride
coeff(61) = 113;     % baseline neuronal intracellular potassium
coeff(62) = 9.0;     % baseline neuronal intracellular sodium
coeff(63) = 9;       % baseline glial intracellular chloride
coeff(64) = 113;     % baseline glial intracellular potassium
coeff(65) = 9.0;     % baseline glial intracellular sodium
coeff(66) = 1;       % neuron-glia pyruvate coupling coefficient
coeff(67) = -1;      % ATP usage coupling to pump activity
coeff(68) = 1;       % neuronal metabolic scaling
coeff(69) = 1;       % glial metabolic scaling
coeff(70) = 2;       % glycolysis threshold parameter
coeff(71) = 0.1;     % glycolysis sigmoid slope
coeff(72) = 3;       % oxidative phosphorylation threshold parameter
coeff(73) = 1;       % oxidative phosphorylation slope
coeff(74) = 3;       % pyruvate threshold parameter
coeff(75) = 1;       % pyruvate slope
coeff(76) = 1;       % neuronal oxidative metabolism scale
coeff(77) = 1;       % glial oxidative metabolism scale
coeff(78) = 1;       % neuronal glycolysis scale
coeff(79) = 1;       % glial glycolysis scale
coeff(80) = 1;       % exponent used in extracellular-volume scaling
coeff(81) = 10;      % osmotic water-coupling strength
coeff(82) = 8e-5;    % water/bath exchange coefficient dw
coeff(83) = 3e-5;    % ionic diffusion coefficient di

%% =========================================================================
% ORIGINAL M2 PARAMETERS
% These are mostly from the original M2 / Ca / receptor code.
% =========================================================================
cm    = 1.0;        % membrane capacitance used in original M2 voltage equation
g_Na  = 30.0;       % fast Na conductance in original M2 code
g_K   = 25.0;       % fast/delayed K conductance in original M2 code
g_NaL = 0.0217;     % Na leak conductance in original M2 code
g_KL  = gkk;        % original M2 K leak parameter selected by gks
% NOTE:
% The full M2 I_K is kept for reference but remains flagged out of dy.
% The value gkk is used only as the RELATIVE modulation signal for the
% volume-model K leak through g_KL_eff below.
g_ClL = 0.049;      % Cl leak conductance
g_TCa = 1.75;       % T-type Ca conductance
g_KCa = 41.322;     % Ca-activated K conductance
g_CAN = 0.826;      % CAN conductance
g_NaCAN = 0.411;    % Na component of CAN current
g_KCAN  = 0.415;    % K component of CAN current

g_GIRK = 0.30;      % GIRK conductance used in BF->ACh->M2 block

BF_tonic     = 0.10; % baseline BF activity
BF_burst_amp = 0.40; % BF increase during stimulation
tau_BF       = 300;  % BF time constant

ACh_base = 0.05;     % baseline ACh level
k_rel    = 0.80;     % BF to ACh coupling
tau_ACh  = 500;      % ACh time constant

Kd_M2    = 0.20;     % M2 half-activation constant
nHill_M2 = 2;        % Hill coefficient for M2 activation
tau_M2   = 600;      % M2 time constant

p_max     = 1.25;    % pump maximum in original M2 code
U_kcc2    = 0.3;     % KCC2 transport strength
G_glia    = 8.0;     % glial buffering strength
epsilon_k = 0.33;    % extracellular K diffusion strength
K_bath    = 4.0;     % bath K concentration in original M2 code
Na_gi     = 18.0;    % glial intracellular Na used in original M2 glial pump term
tau       = 0.001;   % time scaling used in original M2 ion/Ca equations

pi_val = 3.1415;     % pi
cell_rad = 7e-6;     % effective cell radius for M2 geometry
lam_c   = 0.70217;   % cytosolic volume fraction
lam_ds  = 4.0353e-5; % submembrane / dense shell volume fraction
lam_jer = 0.00326;   % junctional ER volume fraction
lam_ner = 0.0570;    % network ER volume fraction

V_i = (4/3)*pi_val*(cell_rad^3); % total cell volume for M2 Ca model
area_cell = 4*pi_val*(7e-4)^2;   % effective surface area used in original code
V_c   = lam_c*V_i;               % cytosol volume
V_ds  = lam_ds*V_i;              % submembrane shell volume
V_jer = lam_jer*V_i;             % junctional ER volume
V_ner = lam_ner*V_i;             % network ER volume

% Extracellular Ca compartment used only for converting membrane Ca fluxes
% into a concentration rate for the new Ca_o state.
% IMPORTANT: the original M2 model did not define a separate extracellular
% Ca compartment volume. For the first consistent implementation we use
% V_o = V_c so that the new state has a defined concentration scale without
% altering the rest of the volume backbone.
V_o = V_c;

alpha_area = 4*pi_val*(3*V_i/(4*pi_val))^(2/3); % effective area factor
Faraday = 6.02e23*1.6e-19;                     % Faraday-like charge conversion
gamma_con_c  = alpha_area/(Faraday*V_c)*0.01;  % current-to-concentration conversion in cytosol
gamma_con_ds = alpha_area/(Faraday*V_ds)*0.01; % current-to-concentration conversion in ds
gamma_con_o  = alpha_area/(Faraday*V_o)*0.01;  % current-to-concentration conversion outside

I_NCX_max = 2600; % maximum NCX current
K_mNa = 87.5;     % NCX sodium affinity
K_mCa = 1.38;     % NCX calcium affinity
k_sat = 0.1;      % NCX saturation parameter
volt  = 0.35;     % NCX voltage dependence parameter
R = 83143;        % gas constant-like scale in original code
T = 310;          % temperature
F_NCX = 964867;   % Faraday constant for NCX equation

% Ca_o is now dynamic and is unpacked below as y(52).
E_CAN = -20.0;    % reversal potential for CAN current

cmdn    = 50e-3;  % calmodulin buffer concentration
Km_cmdn = 2.38e-3;% calmodulin dissociation constant
csqn    = 15.0;   % calsequestrin buffer concentration
Km_csqn = 0.8;    % calsequestrin dissociation constant

t_xfer  = 3.1245e-3; % transfer time constant between ds and cytosol
t_tr    = 0.03448;   % transfer time constant between ER compartments
v_leak  = 0.005;     % ER leak rate
v_serca = 0.2;       % SERCA pump rate
Km_serca= 3e-4;      % SERCA half-saturation constant

v_ryr   = 0.4;       % RyR flux scale
ka_pos  = 1.5e12;    % RyR forward rate a+
ka_neg  = 0.0288;    % RyR backward rate a-
kb_pos  = 1.5e9;     % RyR forward rate b+
kb_neg  = 0.3859;    % RyR backward rate b-
kc_pos  = 0.00175;   % RyR forward rate c+
kc_neg  = 0.0001;    % RyR backward rate c-
n = 4;               % RyR calcium exponent
m = 3;               % RyR calcium exponent

v_ip3 = 24;      % IP3R flux scale
a1 = 1.672e2;    % IP3R transition parameter
a2 = 3.81;       % IP3R transition parameter
a3 = 4.0e2;      % IP3R transition parameter
a4 = 0.3000;     % IP3R transition parameter
a5 = 53.9;       % IP3R transition parameter
b1 = 0.2280;     % IP3R transition parameter
b2 = 0.000409;   % IP3R transition parameter
b3 = 0.1885;     % IP3R transition parameter
b4 = 0.000093135;% IP3R transition parameter
b5 = 0.00452;    % IP3R transition parameter

ki     = 0.0020; %#ok<NASGU> % legacy compatibility parameter from original M2 code
k_plc  = 1.1e-4;             % PLC-related IP3 production parameter
a7_scale = 100; %#ok<NASGU>  % legacy compatibility parameter in original M2 code

v_a7   = 100.0;      % alpha7 Ca influx scaling
d2     = 0.5;        % alpha7 receptor transition parameter
alpha1 = 0.05;       % alpha7 receptor transition parameter
gamma1 = 20;         % alpha7 receptor transition parameter
alpha  = 20;         % alpha7 receptor transition parameter
beta   = 0.8;        % alpha7 receptor transition parameter
d3     = 0.1;        % alpha7 receptor transition parameter
beta_IP3 = 0.3;      % IP3 decay / relaxation coefficient
IP3star  = 0.01932e-4; % baseline IP3 target

%% =========================================================================
% STATE UNPACKING: VOLUME BACKBONE
% =========================================================================
Vn_mem = y(1);  % neuron membrane-potential state from volume backbone
Vg_mem = y(2);  % glial membrane-potential state from volume backbone
n_vol  = y(3);  % HH delayed-rectifier K activation gate (volume code)
h_vol  = y(4);  % HH Na inactivation gate (volume code)

NinNa  = y(5);  % neuronal intracellular Na amount
NinNag = y(6);  % glial intracellular Na amount
NoutNa = y(7);  % extracellular Na amount

NinK   = y(8);  % neuronal intracellular K amount
NinKg  = y(9);  % glial intracellular K amount
NoutK  = y(10); % extracellular K amount

NinCl  = y(11); % neuronal intracellular Cl amount
NinClg = y(12); % glial intracellular Cl amount
NoutCl = y(13); % extracellular Cl amount

Vin    = y(14); % neuronal intracellular volume
Ving   = y(15); % glial intracellular volume
Vout   = y(16); % extracellular / ECS volume

pyrn   = y(17); % neuronal pyruvate-like metabolic pool
ATPn   = y(18); % neuronal ATP-like metabolic pool
pyrg   = y(19); % glial pyruvate-like metabolic pool
ATPg   = y(20); % glial ATP-like metabolic pool

VBE    = y(21); % bath/ECS electrical reference variable

%% =========================================================================
% STATE UNPACKING: ORIGINAL M2-ONLY STATES
% =========================================================================
Po1 = y(22);   % RyR open state 1
Po2 = y(23);   % RyR open state 2
Pc1 = y(24);   % RyR closed state 1
Pc2 = y(25);   % RyR closed state 2

x000 = y(26); % IP3R Markov state x000
x100 = y(27); % IP3R Markov state x100
x110 = y(28); % IP3R Markov state x110
x101 = y(29); % IP3R Markov state x101
x111 = y(30); % IP3R Markov state x111
x001 = y(31); % IP3R Markov state x001

C = y(32);   % alpha7 closed state
O = y(33);   % alpha7 open state
D = y(34);   % alpha7 desensitized / inactive-like state

IP3 = y(35);   % intracellular IP3 state

Ca_c   = y(36); % cytosolic calcium
Ca_ds  = y(37); % submembrane / dense shell calcium
Ca_jer = y(38); % junctional ER calcium
Ca_ner = y(39); % network ER calcium

delay = y(40);  % original M2 timing state
check = y(41); %#ok<NASGU> % original M2 diagnostic state

M2  = y(42);   % muscarinic M2 activation state
ACh = y(43);   % acetylcholine level/state
BF  = y(44);   % basal forebrain drive/activity state

m_na = y(45);
h_na = y(46);
j_na = y(47);
m_ca = y(48);
h_ca = y(49);
m_KCa = y(50);
m_CAN = y(51);
Ca_o  = y(52);   % extracellular calcium concentration

% Numerical guard used ONLY inside logarithms.
% The differential state itself remains Ca_o = y(52).
Ca_o_safe = max(Ca_o,1e-9);
Ca_c_safe = max(Ca_c,1e-12);
Ca_ds_safe = max(Ca_ds,1e-12);

%% =========================================================================
% SHARED PHYSICAL VARIABLES FROM VOLUME BACKBONE
% =========================================================================
radiusg = coeff(53);      % effective glial/cell radius scale
arean   = 4*3.14*coeff(53)^2; % neuronal membrane area
areag   = 4*3.14*radiusg^2;   % glial membrane area

E_ko  = 26.64*log((coeff(58)/coeff(61))); % baseline K reversal
E_nao = 26.64*log((coeff(59)/coeff(62))); % baseline Na reversal
E_clo = 26.64*log((coeff(57)/coeff(60))); % baseline Cl reversal

E_kog  = 26.64*log((coeff(58)/coeff(64))); % baseline glial K reversal
E_naog = 26.64*log((coeff(59)/coeff(65))); % baseline glial Na reversal
E_clog = 26.64*log((coeff(63)/coeff(60))); % baseline glial Cl reversal

Inao = coeff(8)*(coeff(6)-E_nao); % baseline neuronal Na leak current
Iclo = coeff(10)*(coeff(6)-E_clo);% baseline neuronal Cl leak current

Vm = Vn_mem - VBE;   % shared neuron membrane voltage used everywhere in merged model

% convert conserved ion amounts into concentrations
Ko   = NoutK  / Vout / coeff(46); % ECS K concentration from amount/volume
Ki   = NinK   / Vin  / coeff(46); % neuronal intracellular K concentration
Kig  = NinKg  / Ving / coeff(46); % glial intracellular K concentration
Nao  = NoutNa / Vout / coeff(46); % ECS Na concentration
Nai  = NinNa  / Vin  / coeff(46); % neuronal intracellular Na concentration
Naig = NinNag / Ving / coeff(46); % glial intracellular Na concentration
Clo  = NoutCl / Vout / coeff(46); % ECS Cl concentration
Cli  = NinCl  / Vin  / coeff(46); % neuronal intracellular Cl concentration
Clig = NinClg / Ving / coeff(46); % glial intracellular Cl concentration

Inaog = coeff(9)*(coeff(7)-E_naog); % baseline glial Na leak current
Iclog = coeff(11)*(coeff(7)-E_clog);% baseline glial Cl leak current

g_kLg = (-2*Inaog/3 - Iclog)/(coeff(7)-E_kog); % glial K leak conductance inferred from reference balance
Ikog  = (2*g_kLg/(1+exp((coeff(7)-coeff(7))/coeff(4))))*(coeff(7)-E_kog);
ipmpog = -Inaog - Iclog - Ikog; % glial pump reference current

fc2 = Iclo/(log((coeff(61)*coeff(57))/(coeff(58)*coeff(60)))); % cotransporter coefficient neuron
Icc = fc2*(log((Ki*Cli)/(Ko*Clo)));                             % neuronal cotransporter current

fc2g = (-5*Iclog-6*Ikog-4*Inaog)/(log((coeff(64)*coeff(63))/(coeff(58)*coeff(60))));
Iccg = fc2g*(log((Kig*Clig)/(Ko*Clo))); % glial cotransporter current

fc1g = (-3*Iclog-3*Ikog-2*Inaog)/(log((coeff(64)*coeff(63)*coeff(65)*coeff(63))/(coeff(58)*coeff(60)*coeff(59)*coeff(60))));
% NE-modulated glial NKCC/NKCC1-like current
Inkccg = NE_glia_scale * ( -fc1g * log((Kig*Clig*Clig*Naig)/(Ko*Clo*Clo*Nao)) ); % glial NKCC-like current

g_KL_vol_base = (-2*Inao/3-Iclo)/(coeff(6)-E_ko); % original inferred volume-model K leak

% gKs-dependent effective leak used by the ONE main neuronal K current.
% Reference gkk = 0.14 reproduces the original volume-model leak exactly.
gks_reference = 0.14;
g_KL_eff = g_KL_vol_base * ...
           (1 + k_gks_mod*(gkk/gks_reference - 1));

E_k  = 26.64*log((Ko/Ki));
E_na = 26.64*log((Nao/Nai));
E_cl = 26.64*log((Cli/Clo));

E_kg  = 26.64*log((Ko/Kig));
E_nag = 26.64*log((Nao/Naig));
E_clg = 26.64*log((Clig/Clo));

%% =========================================================================
%% VOLUME-CODE HH BLOCK
%% Possible future discussion:
%% This block and the M2 fast HH block below may overlap mechanistically.
%% Both are kept here intentionally.
%% =========================================================================
alpha_n=0.01*(Vm+34.0)/(1.0-exp(-0.1*(Vm+34.0))); % HH n-gate forward rate
beta_n=0.125*exp(-(Vm+44.0)/80.0);                % HH n-gate backward rate
alpha_m=0.1*(Vm+30.0)/(1.0-exp(-0.1*(Vm+30.0))); % HH m-gate forward rate
beta_m=4.0*exp(-(Vm+55.0)/18.0);                  % HH m-gate backward rate
alpha_h=0.07*exp(-(Vm+44.0)/20.0);               % HH h-gate forward rate
beta_h=1.0/(1.0+exp(-0.1*(Vm+14.0)));          % HH h-gate backward rate
m_inf=alpha_m/(alpha_m+beta_m);                     % instantaneous Na activation in volume block

%% =========================================================================
%% VOLUME-CODE STIMULATION BLOCK
%% =========================================================================
t2=t-coeff(18);
a=mod(t2,coeff(16));
b=mod(t2,coeff(20));

if a<coeff(17)
    f=1;
else
    f=0;
end
if b<coeff(21)
    g=1;
else
    g=0;
end
if t<coeff(18)||t>coeff(19)
    wdg=0;
else
    wdg=1;
end

rcur=0.0; % placeholder externally applied current in volume block
tgl=wdg*f*g; % binary stimulation gate from original volume code

%% =========================================================================
%% VOLUME-CODE CURRENTS
%% =========================================================================
Ina_vol=coeff(2)*(m_inf^3)*h_vol*(Vm-E_na)+(coeff(8)+tgl*coeff(33))*(Vm-E_na);
% Main merged neuronal K current:
% keep the volume-model delayed-rectifier gate/conductance,
% but use the gKs-modulated effective leak conductance.
Ik_vol=coeff(1)*(n_vol^4)*(Vm-E_k)+(g_KL_eff+tgl*coeff(35))*(Vm-E_k);
Icl_vol=(coeff(10)+tgl*coeff(34))*(Vm-E_cl);

Vg=Vg_mem-Vn_mem; % glial membrane voltage expressed relative to neuron/ECS reference

% Glial Na leak current is NOT NE-scaled in this version.
Inag = coeff(9) * (Vg - E_nag);

% NE-modulated glial KIR / inward-rectifier K current.
% coeff(4) remains the KIR rectification/shape factor.
% NE_glia_scale changes the effective KIR current strength.
Ikg = NE_glia_scale * (2*g_kLg/(1+exp((Vg-coeff(7))/coeff(4)))) * (Vg - E_kg);

% Glial Cl leak current is NOT NE-scaled in this version.
Iclg = coeff(11) * (Vg - E_clg);

%% =========================================================================
%% VOLUME-CODE PUMP / DIFFUSION / DRIFT / OSMOTIC BLOCK
%% =========================================================================
atpnn=coeff(13);
atpgn=coeff(13);

atpfacn=1/(1+exp((coeff(13)-atpnn)/coeff(14)));
ipmpo=-Inao/3;
impmp=ipmpo/(((1.0/(1.0+(coeff(27)/coeff(62))^1.5))*coeff(58)/(coeff(58)+coeff(28)))*((1-coeff(15))+coeff(15)*atpfacn));

atpfac=1/(1+exp((coeff(13)-ATPn)/coeff(14)));
Itildepump=impmp*((1.0/(1.0+(coeff(27)/Nai)^1.5))*Ko/(Ko+coeff(28)))*((1-coeff(15))+ coeff(15) *atpfac);

atpfacgn=1/(1+exp((coeff(13)-atpgn)/coeff(14)));
impmpg=ipmpog/(((1.0/(1.0+(coeff(25)/coeff(62))^1.5))*coeff(58)/(coeff(58)+coeff(26)))*((1-coeff(15))+coeff(15)*atpfacgn));

atpfacg=1/(1+exp((coeff(13)-ATPg)/coeff(14)));
% NE-modulated glial Na/K pump current.
% NE_glia_scale changes the effective glial pump strength.
itildepmpg = NE_glia_scale * impmpg * ...
    ((1.0/(1.0+(coeff(25)/Naig)^1.5))*Ko/(Ko+coeff(26))) * ...
    ((1-coeff(15))+coeff(15)*atpfacg);

mut=Vout/coeff(41);       % relative ECS volume
mute=mut^coeff(80);       % scaled ECS volume factor
mute2=mut^coeff(40);      % second ECS volume scaling factor
FUJi=mute2*coeff(22);     % FUJ-driven forcing term
areae=arean*mute;         % effective exchange area

ItildedifK=(1e-6)*coeff(83)*areae*(Ko-coeff(54))/coeff(52);  % K diffusion to bath
ItildedifN=(1e-6)*coeff(83)*areae*(Nao-coeff(55))/coeff(52); % Na diffusion to bath
ItildedifC=(1e-6)*coeff(83)*areae*(Clo-coeff(56))/coeff(52); % Cl diffusion to bath

Tcap=(1e-2)*arean*coeff(49)*coeff(50)/coeff(52); % effective capacitance-like term
drift=(1e-8)*areae*VBE/(coeff(51)*coeff(52));    % drift term driven by bath/ECS voltage

charmolcon=(1000/coeff(45))/coeff(44); % charge-to-molar conversion factor

Kave=Ko+coeff(54);
Clave=Clo+coeff(56);
Naave=Nao+coeff(55);

ItildedrK=drift*Kave/(Kave+Naave+Clave);
ItildedrN=drift*Naave/(Kave+Naave+Clave);
ItildedrC=-drift*Clave/(Kave+Naave+Clave);

ItildebK=charmolcon*ItildedrK+ItildedifK;
ItildebN=charmolcon*ItildedrN+ItildedifN;
ItildebC=charmolcon*ItildedrC+ItildedifC;

osmocAa=coeff(29)*coeff(42)/Vin;   % neuronal immobile osmolyte contribution
osmocAag=coeff(30)*coeff(43)/Ving; % glial immobile osmolyte contribution
osmocBa=coeff(31)*coeff(41)/Vout;  % ECS immobile osmolyte contribution

catc1=(5e-3)*(18); % osmotic scaling constant
osmne=coeff(81)*arean*catc1*((Ko+Nao+Clo+osmocBa)-(Cli+Nai+Ki+osmocAa));      % neuron-ECS osmotic term
osmge=coeff(81)*areag*catc1*((Ko+Nao+Clo+osmocBa)-(Clig+Naig+Kig+osmocAag));  % glia-ECS osmotic term

nconv=18;
delvout=(coeff(41)-Vout)/coeff(37); % slow relaxation toward reference ECS volume

Kad=Ko;
if delvout>(0)
    Kad=coeff(54);
end
Nad=Nao;
if delvout >0
    Nad=coeff(55);
end
Clad=Clo;
if delvout >0
    Clad=coeff(56);
end

Kaf=Ko;
if coeff(22)>0
    Kaf=coeff(54);
end
Naf=Nao;
if coeff(22)>0
    Naf=coeff(55);
end
Claf=Clo;
if coeff(22)>0
    Claf=coeff(56);
end

osmocC=coeff(31);
osmeb=coeff(81)*coeff(82)*areae*nconv*((Ko+Nao+Clo+osmocBa)-(coeff(54)+coeff(55)+coeff(56)+osmocC))/coeff(52); % ECS-bath osmotic exchange

vin=(-osmne+FUJi*coeff(23)*tgl-coeff(23)*delvout);  % neuron volume rate
ving=(-osmge+FUJi*coeff(23)*tgl-coeff(23)*delvout); % glia volume rate
vout=(osmge+osmne+osmeb+FUJi*tgl-FUJi*coeff(23)*tgl+delvout-coeff(23)*delvout); % ECS volume rate

ccon=coeff(12)*arean/coeff(44)/coeff(45); % neuron current-to-amount conversion
ccog=coeff(12)*areag/coeff(44)/coeff(45); % glia current-to-amount conversion

%% =========================================================================
%% FIXED-BIN POISSON STIMULATION BLOCK
%% =========================================================================
% The old delay-trigger and stml-long pulse logic are intentionally not used
% for this input. A value of 1 stimulates only during its own Poisson bin.
%
% The complete array is generated once inside the model. Because it is
% persistent, every alpha7 run in the same driver uses the identical random
% input, which makes the alpha7 comparison fair.

% Keep the old delay equation alive for state compatibility.
dtar = dtarget;
rec  = 1/dtar;
if t > 1000 && Vm > 0
    rst = -100;
else
    rst = 0;
end

persistent poisson_array
persistent stored_t_end_ms
persistent stored_rate_Hz
persistent stored_dt_ms

% Create the binary array only when it does not exist or its settings change.
if isempty(poisson_array) || ...
   isempty(stored_t_end_ms) || ...
   stored_t_end_ms ~= t_end_ms || ...
   stored_rate_Hz  ~= poisson_rate_Hz || ...
   stored_dt_ms    ~= poisson_dt_ms

    number_of_bins = floor(t_end_ms/poisson_dt_ms) + 1;
    random_values  = rand(1,number_of_bins);

    % Professor's rule: above threshold -> 1; otherwise -> 0.
    poisson_array = double(random_values > (1-p_event_poisson));

    stored_t_end_ms = t_end_ms;
    stored_rate_Hz  = poisson_rate_Hz;
    stored_dt_ms    = poisson_dt_ms;
end

% MATLAB indexing begins at 1. The floor operation provides a zero-order
% hold: one fixed 0/1 value is used throughout each Poisson bin.
current_poisson_bin = floor(t/poisson_dt_ms) + 1;

if current_poisson_bin >= 1 && ...
   current_poisson_bin <= numel(poisson_array)
    poisson_event = poisson_array(current_poisson_bin);
else
    poisson_event = 0;
end

% Direct binary Poisson stimulation:
%   poisson_event = 1 -> I_app = selected poisson_strength
%   poisson_event = 0 -> I_app = 0
% The original ampst parameter is not used to scale this Poisson input.
% No tstr/tstp/stml pulse is started here.
I_app = poisson_strength * poisson_event;

%% =========================================================================
%% ORIGINAL M2 BF -> ACh -> M2 BLOCK
%% =========================================================================
% The same fixed binary Poisson value drives the BF input.
sBF = poisson_event;

BF_target = BF_tonic + BF_burst_amp*sBF; % target BF activity during a Poisson event bin
dotBF = (BF_target - BF)/tau_BF;

ACh_target = ACh_base + k_rel*BF;        % ACh target driven by BF activity
dotACh = (ACh_target - ACh)/tau_ACh;

ACh_pow = ACh^nHill_M2;
Kd_pow  = Kd_M2^nHill_M2;
M2_inf  = ACh_pow / (Kd_pow + ACh_pow); % steady-state M2 activation by ACh
dotM2 = (M2_inf - M2)/tau_M2;

%% =========================================================================
%% ORIGINAL M2 ALPHA7 / HH / Ca GATING BLOCK
%% =========================================================================
if (Vm>=-78.0)
    g_alpha = 9e-9/(area_cell); % alpha7 conductance only allowed above threshold
else
    g_alpha = 0;
end

a_m = 0.32*(54+Vm)/(1-exp(-(Vm+54)/4));
b_m = 0.28*(Vm+27)/(exp((Vm+27)/5)-1);

a_h = 0.128*exp(-(50+Vm)/18);
b_h = 4/(1+exp(-(Vm+27)/5));

a_j = 0.032*(Vm+52)/(1-exp(-(Vm+52)/5));
b_j = 0.5*exp(-(Vm+57)/40);

t_m  = 0.44 + 0.15/(exp((Vm + 27)/10) + exp(-(Vm + 102)/15));
t_h  = 22.7 + 0.27/(exp((Vm + 48)/4) + exp(-(Vm + 407)/50));
m_ss = 1/(1 + exp(-(Vm + 52)/7.4));
h_ss = 1/(1 + exp((Vm + 80)/5));

t_m_KCa   = 1/(48*((Ca_ds)^2) + 0.03);
m_ss_KCa  = 48*((Ca_ds)^2)/(48*((Ca_ds)^2) + 0.03);

t_m_CAN  = 1/(20*((Ca_ds)^2) + 0.002);
m_ss_CAN = 20*((Ca_ds)^2)/(20*((Ca_ds)^2) + 0.002);

%% =========================================================================
%% ORIGINAL M2 PUMP / DIFFUSION / KCC2 BLOCK
%% Possible future discussion:
%% These terms may overlap with the volume-backbone pump / ion handling.
%% They are kept here intentionally because you asked to preserve the code.
%% =========================================================================
I_pump     = (p_max/(1.0+exp((25.0-Nai)/3.0)))*(1.0/(1.0+exp(5.5-Ko))); % neuronal Na/K pump current
I_gliapump = (p_max/3/(1.0+exp((25.0-Na_gi)/3.0)))*(1.0/(1.0+exp(5.5-Ko))); % glial Na/K pump-like term
I_glia     = G_glia/(1.0+exp((18-Ko)/2.5)); % glial K buffering / uptake term
I_diff     = epsilon_k*(Ko-K_bath);             % K diffusion from ECS toward bath

I_kcc2 = U_kcc2*log((Ki*Cli)/(Ko*Clo));      % KCC2 cotransporter flux

E_Na = 26.64*log(Nao/Nai); % Na Nernst potential
E_K  = 26.64*log(Ko/Ki);   % K Nernst potential
E_Cl = 26.64*log(Cli/Clo); % Cl Nernst potential
E_Ca = (26.64/2)*log(Ca_o_safe/Ca_c_safe); % Ca Nernst potential

I_K  = g_K*j_na^4*(Vm-E_K) + g_KL*(Vm-E_K);          % original M2 K current + K leak

% -------------------------------------------------------------------------
% DUPLICATED M2 K CURRENT HANDLING
% -------------------------------------------------------------------------
% I_K above is retained for reference only.
% It is NOT added to the ion balance, because Ik_vol is now the single
% main neuronal K current and already contains the gKs-modulated leak.
%
I_Na = g_Na*m_na^3*h_na*(Vm - E_Na) + g_NaL*(Vm-E_Na); % original M2 fast Na current + Na leak
I_Cl = g_ClL*(Vm-E_Cl);                                 % original M2 Cl leak current

I_T  = g_TCa*(m_ca^2)*h_ca*(Vm-E_Ca);                % T-type Ca current

I_KCa   = g_KCa*((m_KCa)^2)*(Vm - E_K);             % Ca-activated K current
I_CAN   = g_CAN*((m_CAN)^2)*(Vm - (E_K+E_Na)/2);          % CAN current total form
I_NaCAN = g_NaCAN*((m_CAN)^2)*(Vm - E_Na);         % Na component of CAN current
I_KCAN  = g_KCAN*((m_CAN)^2)*(Vm - E_K);            % K component of CAN current

I_NCX = I_NCX_max*((Nai)^3*Ca_o*exp((volt*Vm*F_NCX)/(R*T))-(Nao)^3*Ca_c* ...
     exp(((volt-1)*Vm*F_NCX)/(R*T)))/(((K_mNa)^3+(Nao)^3)*(K_mCa+Ca_o)* ...
     (1+k_sat*exp(((volt-1)*Vm*F_NCX)/(R*T))));

beta_c   = (1 + (cmdn*Km_cmdn)/(Km_cmdn + Ca_c^2))^-1;
beta_ds  = (1 + (cmdn*Km_cmdn)/(Km_cmdn + Ca_ds^2))^-1;
beta_jer = (1 + (csqn*Km_csqn)/(Km_csqn + Ca_jer^2))^-1;

F_TR = R*T/(2*Faraday);

v_plc = a7on*8.0e-5*D; % PLC activity driven by alpha7 state D in original code

J_xfer = (Ca_ds - Ca_c)/t_xfer;             % Ca transfer between dense shell and cytosol
J_tr   = (Ca_ner - Ca_jer)/t_tr;            % Ca transfer between ER compartments
J_leak = v_leak*(Ca_ner - Ca_c);            % passive ER-to-cytosol Ca leak
J_serca= 0.1*v_serca*(Ca_c^2/(Km_serca^2 + Ca_c^2)); % SERCA uptake into ER
J_ryr  = 20*v_ryr*(Po1 + Po2)*(Ca_jer - Ca_ds);     % RyR-mediated Ca release
J_ip3  = 100*v_ip3*(x110^4+4*x110^3*(1-x110))*(Ca_jer - Ca_ds); % IP3R-mediated Ca release
J_a7   = a7on*v_a7*O*(Ca_o - Ca_ds);       % alpha7-driven Ca influx term
J_plc  = v_plc*(Ca_ds^2/(k_plc^2 + Ca_ds^2)); % PLC/IP3 production term

I_a7 = a7on*g_alpha*(Vm - F_TR*log(Ca_o_safe/Ca_ds_safe)); % alpha7 receptor current

Vhalf_GIRK = E_K - 15;
krect_GIRK = 7;
f_rect = 1 ./ (1 + exp((Vm - Vhalf_GIRK)/krect_GIRK));
I_GIRK = g_GIRK * M2 .* f_rect .* (Vm - E_K); % M2-activated GIRK current

I_pmca_bar = 0.1;
eta_pmca   = 2.0;
Km_pmca    = 0.5;
I_pmca     = I_pmca_bar*(Ca_c^eta_pmca/(Km_pmca^eta_pmca + Ca_c^eta_pmca)); % PMCA Ca extrusion current
J_pmca     = (gamma_con_c)*I_pmca; % PMCA current converted to Ca flux

%% =========================================================================
%% CLEANED MERGE CONTROL: AVOID OVERCOUNTING IN dy EQUATIONS
%% =========================================================================
% The equations for I_Na, I_K, I_Cl, I_pump, I_kcc2, etc. are still
% calculated above.  However, by default they are NOT allowed to drive the
% main Na/K/Cl amount balances, because the volume backbone already has its
% own Na, K, Cl currents, pumps, cotransport, glial buffering, and bath
% diffusion.
%
% Default philosophy in this cleaned version:
%   1) Volume backbone controls Na/K/Cl amounts and Vin/Ving/Vout.
%   2) M2/alpha7/Ca machinery is preserved.
%   3) Only the M2-specific GIRK current is allowed to modulate voltage.
%   4) Duplicated M2 Na/K/Cl/pump/KCC2/glial/bath terms are set to zero
%      before entering dy.
%
% Later, you can turn individual terms back on by changing these flags from
% 0 to 1, but test one flag at a time.
USE_M2_NA_IN_ION_BALANCE      = 0;  % I_Na / I_NaCAN / NCX Na contribution
USE_M2_K_IN_ION_BALANCE       = 0;  % duplicated original M2 I_K only; KCa/KCAN/GIRK retained separately
USE_M2_CL_IN_ION_BALANCE      = 0;  % I_Cl contribution
USE_M2_PUMP_IN_ION_BALANCE    = 0;  % I_pump / I_gliapump contribution
USE_M2_KCC2_IN_ION_BALANCE    = 0;  % I_kcc2 contribution
USE_M2_GLIA_BATH_K_BALANCE    = 0;  %    I_glia / I_diff contribution
USE_GIRK_IN_VOLTAGE_EQUATION  = 1;  % M2-specific voltage modulation
USE_ALPHA7_IN_VOLTAGE_EQUATION= 1;  % keep off initially for stability test
USE_T_CA_IN_VOLTAGE_EQUATION  = 1;  % keep off initially for stability test
USE_KCA_CAN_NCX_IN_VOLTAGE    = 1;  % keep off initially for stability test

% Effective terms that are allowed to affect dy.
I_Na_dy     = USE_M2_NA_IN_ION_BALANCE   * I_Na;
I_NaCAN_dy  = 1  * I_NaCAN;
I_NCX_Na_dy = 1  * I_NCX;

I_K_dy      = USE_M2_K_IN_ION_BALANCE    * I_K;
I_KCa_dy    = 1    * I_KCa;
I_KCAN_dy   = 1   * I_KCAN;
I_GIRK_K_dy = 1    * I_GIRK;

I_Cl_dy     = USE_M2_CL_IN_ION_BALANCE   * I_Cl;
I_pump_dy   = USE_M2_PUMP_IN_ION_BALANCE * I_pump;
I_gliapump_dy = USE_M2_PUMP_IN_ION_BALANCE * I_gliapump;
I_kcc2_dy   = USE_M2_KCC2_IN_ION_BALANCE * I_kcc2;
I_glia_dy   = USE_M2_GLIA_BATH_K_BALANCE * I_glia;
I_diff_dy   = USE_M2_GLIA_BATH_K_BALANCE * I_diff;

I_GIRK_V_dy = USE_GIRK_IN_VOLTAGE_EQUATION   * I_GIRK;
I_a7_V_dy   = USE_ALPHA7_IN_VOLTAGE_EQUATION * I_a7;
I_T_V_dy    = USE_T_CA_IN_VOLTAGE_EQUATION   * I_T;
I_extraCaK_V_dy = USE_KCA_CAN_NCX_IN_VOLTAGE * (I_CAN + I_NCX + I_KCa + I_pmca);

%% =========================================================================
%% CLEANED MERGED VOLTAGE EQUATION
%% =========================================================================
% Volume-backbone currents remain the main voltage driver.
% Duplicated M2 Na/K/Cl/pump currents are NOT inserted here by default.
% Only selected modulatory currents, currently GIRK, are allowed through.
dVn_mem = (1.0/coeff(48))*( ...
          -Ina_vol - Ik_vol - Icl_vol - Itildepump + rcur + coeff(32)*tgl ...
          + (-1/cm)*( I_GIRK_V_dy + I_a7_V_dy + I_T_V_dy + I_extraCaK_V_dy - I_app ) ...
          );

dVg_mem = (1.0/coeff(48))*(-Inag-Ikg-Iclg-itildepmpg); % glial voltage dynamics from volume backbone currents

%% =========================================================================
%% VOLUME-CODE GATING EQUATIONS
%% =========================================================================
dn_vol = coeff(47)*(alpha_n*(1-n_vol)-beta_n*n_vol);
dh_vol = coeff(47)*(alpha_h*(1-h_vol)-beta_h*h_vol);

%% =========================================================================
%% ORIGINAL M2 GATING EQUATIONS
%% =========================================================================
dm_na = (a_m*(1-m_na)-b_m*m_na);
dh_na = (a_h*(1-h_na)-b_h*h_na);
dj_na = (a_j*(1-j_na)-b_j*j_na);

dm_ca  = (-1/t_m)*(m_ca - m_ss);
dh_ca  = (-1/t_h)*(h_ca - h_ss);
dm_KCa = (-1/t_m_KCa)*(m_KCa - m_ss_KCa);
dm_CAN = (-1/t_m_CAN)*(m_CAN - m_ss_CAN);

%% =========================================================================
%% CLEANED MERGED ION-AMOUNT EQUATIONS
%% =========================================================================
% The volume-code ion amount equations are the backbone.
% The M2 ionic flux terms remain defined above, but their effective dy terms
% are zero by default through the *_dy variables defined in the cleaned merge
% control block. This prevents double-counting Na/K/Cl currents, pumps, KCC2,
% glial K uptake, and bath diffusion.

dNinNa = (ccon)*((-Ina_vol-3.0*Itildepump)) ...
       + Vin*coeff(46)*tau*(-gamma_con_c*I_Na_dy - 3.0*I_pump_dy ...
            - gamma_con_c*3.0*I_NCX_Na_dy - gamma_con_c*I_NaCAN_dy);

dNinNag = (ccog)*((-Inag-3.0*itildepmpg+Inkccg));

dNoutNaB = -ItildebN+FUJi*tgl*Naf*coeff(46)+delvout*Nad*coeff(46);
dNoutNaC = ccon*(Ina_vol+3.0*Itildepump)+ccog*(Inag+3.0*itildepmpg-Inkccg) ...
         + Vout*coeff(46)*tau*(gamma_con_c*I_Na_dy + 3.0*I_pump_dy ...
            + gamma_con_c*3.0*I_NCX_Na_dy + gamma_con_c*I_NaCAN_dy);
dNoutNa = dNoutNaB + dNoutNaC;

% K balance:
% Ik_vol is the single main neuronal K current.
% Its leak component already contains the gKs-dependent g_KL_eff.
dNinK = ccon*((-Ik_vol - Icc + 2.0*Itildepump)) ...
      + Vin*coeff(46)*tau*( -gamma_con_c*I_K_dy + 2.0*I_pump_dy - I_kcc2_dy ...
            - gamma_con_c*I_KCa_dy - gamma_con_c*I_KCAN_dy - gamma_con_c*I_GIRK_K_dy );

dNinKg = ccog*((-Ikg-Iccg+2.0*itildepmpg+Inkccg));

% Equal/opposite neuronal K-channel contribution enters the ECS balance.
dNoutKC = ccon*(Ik_vol - 2.0*Itildepump + Icc) ...
        + ccog*(Ikg - 2.0*itildepmpg + Iccg - Inkccg) ...
        + Vout*coeff(46)*tau*( gamma_con_c*I_K_dy - 2.0*I_pump_dy + I_kcc2_dy ...
            - I_diff_dy - I_glia_dy - 2.0*I_gliapump_dy ...
            + gamma_con_c*I_KCa_dy + gamma_con_c*I_KCAN_dy + gamma_con_c*I_GIRK_K_dy );
dNoutKB = FUJi*tgl*Kaf*coeff(46)-ItildebK+delvout*Kad*coeff(46);
dNoutK = dNoutKC + dNoutKB;

dNinCl = ccon*((Icl_vol-Icc)) ...
       + Vin*coeff(46)*tau*( gamma_con_c*I_Cl_dy - I_kcc2_dy );

dNinClg = ccog*((Iclg-Iccg+2*Inkccg));

dNoutClB = -ItildebC+FUJi*tgl*Claf*coeff(46)+delvout*Clad*coeff(46);
dNoutClC = ccon*(-Icl_vol+Icc)+(-Iclg+Iccg-2*Inkccg)*ccog ...
         + Vout*coeff(46)*tau*(-gamma_con_c*I_Cl_dy + I_kcc2_dy);
dNoutCl = dNoutClB + dNoutClC;

%% =========================================================================
%% VOLUME / METABOLISM BACKBONE EQUATIONS
%% =========================================================================
glyn=coeff(78)/(1+exp((pyrn-coeff(70))/coeff(71)));
oxpn=coeff(76)/((1+exp((ATPn-coeff(72))/coeff(73)))*(1+exp((-pyrn+coeff(74))/coeff(75))));

glyg=coeff(79)/(1+exp((pyrg-coeff(70))/coeff(71)));
oxpg=coeff(77)/((1+exp((ATPg-coeff(72))/coeff(73)))*(1+exp((-pyrg+coeff(74))/coeff(75))));

dpyrn=(2*glyn-oxpn+(pyrg-pyrn)*coeff(66))*coeff(68);
dATPn=Itildepump*coeff(67)+(2*glyn+32*oxpn)*coeff(68);

dpyrg=(glyg-oxpg+(pyrn-pyrg)*coeff(66))*coeff(69);
dATPg=itildepmpg*coeff(67)+(2*glyg+32*oxpg)*coeff(68);

dVout = vout; % ECS volume rate
dVin  = vin;  % neuronal volume rate
dVing = ving; % glial volume rate

dVBE=(-dNoutCl+dNoutK+dNoutNa)/(charmolcon*Tcap); % bath/ECS electrical reference dynamics

%% =========================================================================
%% ORIGINAL M2 RYR / IP3 / ALPHA7 / Ca EQUATIONS
%% =========================================================================
dPo1 = ka_pos*Ca_c^n*Pc1 + kc_neg*Pc2 + kb_neg*Po2 ...
     - kc_pos*Po1 - kb_pos*Ca_c^m*Po1 - ka_neg*Po1;
dPo2 = kb_pos*Ca_c^m*Po1 - kb_neg*Po2;
dPc1 = ka_neg*Po1 - ka_pos*Ca_c^n*Pc1;
dPc2 = kc_pos*Po1 - kc_neg*Pc2;

dx000 = b4*x001 - a4*Ca_ds*x000 + b1*x100 - a1*IP3*x000;
dx100 = a1*IP3*x000 - b1*x100 + b5*x110 - a5*Ca_ds*x100 + b2*x101 - a2*Ca_ds*x100;
dx110 = a5*Ca_ds*x100 - b5*x110 + b2*x111 - a2*Ca_ds*x110;
dx101 = a3*IP3*x001 - b3*x101 + a2*Ca_ds*x100 - b2*x101 + b5*x111 - a5*Ca_ds*x101;
dx111 = a5*Ca_ds*x101 + a2*Ca_ds*x110 - b5*x111 - b2*x111;
dx001 = a4*Ca_ds*x000 - b4*x001 + b3*x101 - a3*IP3*x001;

dC = beta*O + d2*D - alpha*Ca_ds*C - alpha1*C;
dO = alpha*Ca_ds*C + d3*D - beta*O - gamma1*O;
dD = gamma1*O + alpha1*C - d2*D - d3*D;

dIP3 = (beta_IP3*(IP3star - IP3) + J_plc);

dCa_c  = tau*beta_c*( J_xfer + (gamma_con_c/2)*2*I_NCX - J_serca + J_leak - J_pmca );
dCa_ds = tau*beta_ds*( -(gamma_con_ds/2)*I_T - J_xfer*(V_c/V_ds) + J_ryr*(V_jer/V_ds) + J_ip3*(V_jer/V_ds) + J_a7 );
dCa_jer= tau*beta_jer*( J_tr - J_ryr - J_ip3 );
dCa_ner= tau*( J_serca*(V_c/V_ner) - J_leak*(V_c/V_ner) - J_tr*(V_jer/V_ner) );

% =========================================================================
% EXTRACELLULAR CALCIUM STATE y(52)
% =========================================================================
% Initial condition in the drivers: Ca_o(0) = 2.0 mM.
%
% IMPORTANT:
% Do NOT mix raw electrical currents directly with J_pmca.
% Each membrane current is first converted into a Ca concentration flux.
%
% Sign convention follows the existing intracellular Ca equations:
%   I_T   : Ca entering Ca_ds is removed from outside
%   J_a7  : alpha7 Ca entering Ca_ds is removed from outside
%   I_NCX : Ca added to cytosol is removed from outside
%   J_pmca: Ca removed from cytosol is added outside
%
% The volume ratios convert an intracellular concentration change into the
% corresponding extracellular concentration change.
dCa_o = tau*( ...
          +(gamma_con_o/2)*I_T ...
          - J_a7*(V_ds/V_o) ...
          - gamma_con_o*I_NCX ...
          + J_pmca*(V_c/V_o) );

ddelay = rst*delay + rec*(1-delay); % original M2 delay-variable dynamics
dcheck = (I_app - poisson_strength); % diagnostic relative to selected Poisson strength

%% =========================================================================
%% OUTPUT VECTOR
%% Each dy(i) matches the state ordering listed at the top of the file.
%% This makes it easier to trace one state from declaration -> equations -> output.
%% =========================================================================
dy(1)  = dVn_mem;
dy(2)  = dVg_mem;
dy(3)  = dn_vol;
dy(4)  = dh_vol;
dy(5)  = dNinNa;
dy(6)  = dNinNag;
dy(7)  = dNoutNa;
dy(8)  = dNinK;
dy(9)  = dNinKg;
dy(10) = dNoutK;
dy(11) = dNinCl;
dy(12) = dNinClg;
dy(13) = dNoutCl;
dy(14) = dVin;
dy(15) = dVing;
dy(16) = dVout;
dy(17) = dpyrn;
dy(18) = dATPn;
dy(19) = dpyrg;
dy(20) = dATPg;
dy(21) = dVBE;

dy(22) = dPo1;
dy(23) = dPo2;
dy(24) = dPc1;
dy(25) = dPc2;
dy(26) = dx000;
dy(27) = dx100;
dy(28) = dx110;
dy(29) = dx101;
dy(30) = dx111;
dy(31) = dx001;
dy(32) = dC;
dy(33) = dO;
dy(34) = dD;
dy(35) = dIP3;
dy(36) = dCa_c;
dy(37) = dCa_ds;
dy(38) = dCa_jer;
dy(39) = dCa_ner;
dy(40) = ddelay;
dy(41) = dcheck;
dy(42) = dotM2;
dy(43) = dotACh;
dy(44) = dotBF;
dy(45) = dm_na;
dy(46) = dh_na;
dy(47) = dj_na;
dy(48) = dm_ca;
dy(49) = dh_ca;
dy(50) = dm_KCa;
dy(51) = dm_CAN;
dy(52) = dCa_o;

end
