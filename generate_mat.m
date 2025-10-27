clc, clear

opts = struct('R',25e-6,'N_target',5e3,'alpha_rt',0.7, ...
              'lattice','fcc','jitter_sigma',[], 'k_candidate',12, 'L_min',0.2e-6, ...
              'outdir','fig');
PNM  = buildSphericalPNM(opts, 'param_bank','out/param_bank.mat');
% PNM  = buildSphericalPNM_varN(opts, 'param_bank','out/param_bank.mat');
% mean(PNM.P.r_p)
% mean(PNM.T.r_t)