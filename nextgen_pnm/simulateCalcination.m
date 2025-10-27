function history = simulateCalcination(PNM, reactionSpec)
% simulateCalcination  Track CaCO3→CaO conversion with newborn micropores.
%
%   history = simulateCalcination(PNM, reactionSpec)
%   reactionSpec fields (optional):
%       .t_end_s            total simulated time [s] (default 100)
%       .dt_s               time step [s] (default 1)
%       .rate_constant      pseudo-first-order rate constant [1/s] (default 0.01)
%       .caO_porosity_frac  fraction of consumed solid volume that becomes new
%                           microporous volume (default 0.65)
%       .birth_threshold    CaCO3 volume fraction triggering new pore creation
%                           (default 0.15)
%       .max_new_pores_per_event  cap per parent pore (default 4)
%
%   Returns a struct containing time series of global conversion, porosity and
%   total BET area.  The routine modifies the supplied PNM in-place (MATLAB
%   handles structs by value, so history.finalPNM contains the updated network).

if nargin < 2
    reactionSpec = struct();
end
if ~isfield(reactionSpec, 't_end_s'), reactionSpec.t_end_s = 100; end
if ~isfield(reactionSpec, 'dt_s'), reactionSpec.dt_s = 1; end
if ~isfield(reactionSpec, 'rate_constant'), reactionSpec.rate_constant = 0.01; end
if ~isfield(reactionSpec, 'caO_porosity_frac'), reactionSpec.caO_porosity_frac = 0.65; end
if ~isfield(reactionSpec, 'birth_threshold'), reactionSpec.birth_threshold = 0.15; end
if ~isfield(reactionSpec, 'max_new_pores_per_event'), reactionSpec.max_new_pores_per_event = 4; end
if ~isfield(reactionSpec, 'new_pore_mean_nm'), reactionSpec.new_pore_mean_nm = 35; end
if ~isfield(reactionSpec, 'new_pore_sigma_nm'), reactionSpec.new_pore_sigma_nm = 12; end

[phi_init, ~] = estimatePhiSSA(PNM.meta.R, PNM.P.r_p, PNM.T.r_t, PNM.T.L_geom);
vol_total = PNM.targets.volume_m3;
solid_vol0 = vol_total * (1 - phi_init);

% Distribute solid volume among pores proportional to surface area.
pore_area = 4*pi*PNM.P.r_p.^2;
area_sum = sum(pore_area);
solid_vol_per_pore = solid_vol0 * (pore_area / area_sum);
state.caCO3_volume = solid_vol_per_pore;
state.caCO3_volume0 = solid_vol_per_pore;
state.caCO3_initial = solid_vol_per_pore;
state.class_id = PNM.P.class_id(:);
state.r_p = PNM.P.r_p(:);
state.coords = PNM.P.coords;
state.is_surface = PNM.P.is_surface;
state.alpha = PNM.meta.alpha_rt;

num_steps = ceil(reactionSpec.t_end_s / reactionSpec.dt_s) + 1;

history.time = zeros(num_steps, 1);
history.global_conversion = zeros(num_steps, 1);
history.porosity = zeros(num_steps, 1);
history.SSA = zeros(num_steps, 1);

for step = 1:num_steps
    t = (step-1) * reactionSpec.dt_s;
    history.time(step) = t;
    remaining = sum(state.caCO3_volume);
    history.global_conversion(step) = 1 - remaining / sum(state.caCO3_initial);
    [phi_est, SSA_est] = estimatePhiSSA(PNM.meta.R, PNM.P.r_p, PNM.T.r_t, PNM.T.L_geom);
    history.porosity(step) = phi_est;
    history.SSA(step) = SSA_est;

    if step == num_steps
        break;
    end

    % Reaction step: shrink CaCO3 volume following exponential decay.
    decay = reactionSpec.rate_constant * reactionSpec.dt_s;
    state.caCO3_volume = state.caCO3_volume .* exp(-decay);

    % Identify pores that crossed the birth threshold and have not yet converted.
    fraction_remaining = ones(size(state.caCO3_volume));
    mask = state.caCO3_volume0 > 0;
    fraction_remaining(mask) = state.caCO3_volume(mask) ./ state.caCO3_volume0(mask);
    newborn_candidates = find(mask & fraction_remaining < reactionSpec.birth_threshold);

    for idx = reshape(newborn_candidates, 1, [])
        volume_consumed = state.caCO3_volume0(idx) - state.caCO3_volume(idx);
        if volume_consumed <= 0
            continue;
        end
        new_volume = volume_consumed * reactionSpec.caO_porosity_frac;
        [PNM, state] = newPoreInsertion(PNM, state, idx, new_volume, reactionSpec);
        % Reset counter to avoid repeated insertions from the same pore
        state.caCO3_volume0(idx) = state.caCO3_volume(idx);
    end
end

history.finalPNM = PNM;
end

function [phi_est, SSA_est] = estimatePhiSSA(R, rp, rt, Lt)
V_sphere = 4/3*pi*R^3;
V_pores = sum(4/3*pi*rp.^3);
V_throats = sum(pi*rt.^2 .* Lt);
phi_est = (V_pores + V_throats) / V_sphere;

A_pores = sum(4*pi*rp.^2);
A_throats = sum(2*pi*rt .* Lt);
SSA_est = (A_pores + A_throats) / V_sphere;
end
