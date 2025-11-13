%% Pore network generation script based on CT-derived pore diameter distribution
% This script builds a spherical pore network without relying on the legacy
% param bank + builder pipeline. It reads EqDiameter data from CT_Pores.csv,
% estimates the pore count for a 1 mm diameter sample, matches a 12%% target
% porosity (pore volume only), connects pores, computes throat sizes using a
% fixed beta coefficient, and provides basic visualizations.

%% Configuration
sample_diameter = 1e-3;              % [m] sample diameter = 1 mm
R = sample_diameter/2;               % [m] sphere radius
phi_target = 0.12;                   % target porosity (pore volume only)
beta = 0.45;                         % throat diameter factor relative to min pore diameter
rng(42);                             % reproducibility for sampling / connections

%% Load pore diameter distribution from CT data
ct_file = fullfile('CT_Pores.csv');
if ~isfile(ct_file)
    error('CT_Pores.csv not found in the repository root.');
end
T = readtable(ct_file);
if ~ismember('EqDiameter', T.Properties.VariableNames)
    error('CT_Pores.csv must contain an ''EqDiameter'' column.');
end
eq_diam_um = T.EqDiameter;
eq_diam_um = eq_diam_um(~isnan(eq_diam_um) & eq_diam_um > 0);
if isempty(eq_diam_um)
    error('EqDiameter column is empty after filtering invalid entries.');
end
eq_diam_m = eq_diam_um * 1e-6;  % convert micrometers to meters

%% Estimate pore count to reach target porosity (initial guess)
V_sphere = (4/3) * pi * R^3;
eq_radius_samples = eq_diam_m / 2;
V_pore_samples = (4/3) * pi * (eq_radius_samples .^ 3);
mean_pore_volume = mean(V_pore_samples);
N_est = max(1, round(phi_target * V_sphere / mean_pore_volume));
fprintf('Initial pore count estimate: %d\n', N_est);

%% Sample pore diameters and rescale to match porosity exactly
idx = randi(numel(eq_diam_m), N_est, 1);
pore_diam_m = eq_diam_m(idx);
pore_radius_m = pore_diam_m / 2;
phi_initial = sum((4/3) * pi * (pore_radius_m .^ 3)) / V_sphere;
scale_factor = (phi_target / phi_initial)^(1/3);
pore_radius_m = pore_radius_m * scale_factor;
pore_diam_m = pore_radius_m * 2;
phi_final = sum((4/3) * pi * (pore_radius_m .^ 3)) / V_sphere;
fprintf('Porosity after scaling: %.4f (target %.4f)\n', phi_final, phi_target);

%% Place pores inside the sphere
opts = struct();
opts.lattice = 'fcc';
opts.k_candidate = 18;
opts.jitter_sigma = [];
[coords, is_surface, lattice_s, jitter_sigma] = placePoresInSphere(R, numel(pore_radius_m), opts);
N = size(coords, 1);
if N ~= numel(pore_radius_m)
    error('placePoresInSphere returned %d pores, expected %d.', N, numel(pore_radius_m));
end

%% Assign target coordination numbers (same PMF as legacy builder)
z_vals = 3:7;
z_pmf = [0.1, 0.25, 0.3, 0.25, 0.1];
z_pmf = z_pmf / sum(z_pmf);
u = rand(N,1);
cdf = cumsum(z_pmf);
z_target = arrayfun(@(val) z_vals(find(val <= cdf, 1, 'first')), u);
z_target(isnan(z_target)) = z_vals(end);

%% Connect pores using kNN-based algorithm
k_candidate = opts.k_candidate;
[pore1, pore2, center_distance] = connectPores(coords, z_target, k_candidate);
M = numel(pore1);

%% Derive throat properties (beta scaling and geometric length)
min_pore_diam = min(pore_diam_m(pore1), pore_diam_m(pore2));
throat_diam_m = beta * min_pore_diam;
throat_radius_m = throat_diam_m / 2;
throat_length_m = max(center_distance - (pore_radius_m(pore1) + pore_radius_m(pore2)), 1e-9);

%% Simple hydraulic conductance proxy (unit diffusivity)
D = 1.0;
G = D .* (pi .* throat_radius_m.^2) ./ max(throat_length_m, eps);

%% Assemble incidence matrix
[B, K] = assembleIncidence(N, pore1, pore2, G);

%% Package output structure and save
PNM = struct();
PNM.meta = struct('R', R, 'phi_target', phi_target, 'phi_achieved', phi_final, ...
                  'beta', beta, 'sample_diameter', sample_diameter, ...
                  'lattice_spacing', lattice_s, 'jitter_sigma', jitter_sigma);
PNM.P = struct('coords', coords, 'radius', pore_radius_m, 'diameter', pore_diam_m, ...
               'is_surface', is_surface);
PNM.T = struct('pore1', pore1, 'pore2', pore2, 'diameter', throat_diam_m, ...
               'radius', throat_radius_m, 'length', throat_length_m, 'conductance', G);
PNM.graph = struct('B', B, 'K', spdiags(K, 0, M, M));
PNM.summary = struct('num_pores', N, 'num_throats', M);

save('PNM_CT.mat', 'PNM', '-v7.3');
fprintf('Saved PNM_CT.mat with %d pores and %d throats.\n', N, M);

%% Visualization: 3D network
figure('Name','Pore Network Visualization');
hold on;
for e = 1:M
    plot3([coords(pore1(e),1), coords(pore2(e),1)], ...
          [coords(pore1(e),2), coords(pore2(e),2)], ...
          [coords(pore1(e),3), coords(pore2(e),3)], 'Color', [0.7 0.7 0.7], 'LineWidth', 0.5);
end
scatter3(coords(:,1), coords(:,2), coords(:,3), ...
         max((pore_diam_m*1e6)*2, 5), pore_diam_m*1e6, 'filled');
colormap(parula);
cb = colorbar;
cb.Label.String = 'Pore diameter [\mum]';
axis equal tight;
xlabel('x [m]'); ylabel('y [m]'); zlabel('z [m]');
title('Spherical pore network');
view(3);
grid on;
hold off;

%% Visualization: diameter distributions (log10 micrometers)
log_pore_diam = log10(pore_diam_m * 1e6);
log_throat_diam = log10(throat_diam_m * 1e6);

figure('Name','Pore Diameter Distribution');
histogram(log_pore_diam, 'BinWidth', 0.05, 'FaceColor', [0.2 0.4 0.7]);
xlabel('log_{10}(Pore diameter [\mum])');
ylabel('Count');
title('Pore diameter distribution');
grid on;

figure('Name','Throat Diameter Distribution');
histogram(log_throat_diam, 'BinWidth', 0.05, 'FaceColor', [0.7 0.3 0.2]);
xlabel('log_{10}(Throat diameter [\mum])');
ylabel('Count');
title('Throat diameter distribution');
grid on;
