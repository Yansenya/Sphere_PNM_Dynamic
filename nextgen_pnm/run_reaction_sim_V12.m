function history = run_reaction_sim_V12()
% run_reaction_sim_V12  Demonstration driver for the next-gen CaCO3 PNM.
%
%   This script mirrors the legacy entry point name so existing automation can
%   call it, but internally it now constructs the network from experimental
%   inputs and launches the calcination simulator with micropore birth logic.

materialSpec = struct( ...
    'pellet_radius_um', 30, ...
    'SSA_total_m2_per_g', 0.57, ...
    'SSA_micro_m2_per_g', 0.39, ...
    'SSA_meso_macro_m2_per_g', 0.18, ...
    'pore_volume_cm3_per_g', 0.0017, ...
    'mean_pore_radius_nm', 8.8);

opts = struct('N_target', 6000, 'alpha_rt', 0.7, 'random_seed', 1234);
PNM = buildSphericalPNM(materialSpec, opts);

reactionSpec = struct('t_end_s', 300, 'dt_s', 5, 'rate_constant', 0.012, ...
    'caO_porosity_frac', 0.7, 'birth_threshold', 0.2, ...
    'max_new_pores_per_event', 3, 'new_pore_mean_nm', 45, ...
    'new_pore_sigma_nm', 18);

history = simulateCalcination(PNM, reactionSpec);

% Report summary metrics.
final = history.finalPNM;
QC = validatePNM(final, final.templateGeom); %#ok<NASGU>

fprintf('Final conversion: %.2f%%\n', history.global_conversion(end)*100);
fprintf('Final porosity: %.4f\n', history.porosity(end));
fprintf('Final SSA: %.4f m^2\n', history.SSA(end));
end
