function PNM = buildSphericalPNM(materialSpec, opts)
% buildSphericalPNM  Material-driven spherical PNM generator for CaCO3 pellets.
%
%   PNM = buildSphericalPNM(materialSpec, opts) consumes experimental BET/porosity
%   data (see materialSpec template below) and constructs a multiscale pore
%   network consistent with the provided targets.  Compared with the legacy
%   implementation this entry point
%     * derives geometric targets using materialToPNMParams,
%     * enforces class-aware coordination numbers, and
%     * keeps a record of pore class membership for downstream calcination logic.
%
%   Required fields in materialSpec (units in parentheses):
%       .pellet_radius_um           radius of the spherical sample (µm)
%       .SSA_total_m2_per_g         total BET surface area (m^2/g)
%       .SSA_micro_m2_per_g         micropore surface area (m^2/g)
%       .SSA_meso_macro_m2_per_g    meso+macropore surface area (m^2/g)
%       .pore_volume_cm3_per_g      total pore volume (cm^3/g)
%       .mean_pore_radius_nm        mean pore radius from nitrogen adsorption (nm)
%
%   Optional materialSpec fields:
%       .solid_density_kg_per_m3    bulk density of solid CaCO3 (default 2710)
%       .mass_scaling_g             override sample mass (g) if provided
%
%   opts fields (all optional):
%       .N_target           representative pore count (default 6000)
%       .alpha_rt           throat radius cap multiplier (default 0.7)
%       .random_seed        RNG seed for reproducibility
%       .placement          struct overriding placement heuristics
%       .outdir             visualisation output folder
%
%   The resulting structure mirrors the legacy `PNM` layout while adding
%   `PNM.P.class_id` (1:micro, 2:meso, 3:macro) and `PNM.targets` for quick
%   access to the imposed material constraints.

if nargin < 1
    error('materialSpec structure is required.');
end
if nargin < 2
    opts = struct();
end

if ~isfield(opts, 'N_target'), opts.N_target = 6000; end
if ~isfield(opts, 'alpha_rt'), opts.alpha_rt = 0.7; end
if ~isfield(opts, 'random_seed'), opts.random_seed = 42; end
if ~isfield(opts, 'outdir'), opts.outdir = fullfile(pwd, 'fig'); end
if ~exist(opts.outdir, 'dir'), mkdir(opts.outdir); end

rng(opts.random_seed, 'twister');

% Translate material measurements into geometric targets and class templates.
geom = materialToPNMParams(materialSpec, opts);
alpha = geom.alpha_rt;

% Place pores respecting class-specific radial preferences.
[coords, class_id, is_surface, placementMeta] = placePoresInSphere(geom, opts);
N = size(coords, 1);

% Assign desired coordination numbers per node.
coordProfile = geom.coordProfile;
z_target = zeros(N, 1);
for c = 1:numel(coordProfile)
    ids = find(class_id == c);
    if isempty(ids), continue; end
    range = coordProfile(c).target_range;
    z_target(ids) = randi(range, numel(ids), 1);
end

% Connect pores using class-aware degree caps.
[T_p1, T_p2, L_geom] = connectPores(coords, class_id, z_target, geom);
M = numel(T_p1);

% Draw pore radii according to class distributions.
r_p = geom.pore_radii;

% Assign throat radii/shape/length while respecting class combinations.
[T_r, T_shape, T_L] = assignThroatProps(coords, T_p1, T_p2, L_geom, r_p, class_id, geom);

% Compute initial conductances and adjust radii to match porosity/SSA targets.
[T_size] = computeSizeFactors(r_p, T_p1, T_p2, T_r, T_L, T_shape);
[~, ~, r_p, T_r] = calibratePorositySSA(geom.R_sample, r_p, T_r, T_L, ...
    geom.phi_target, geom.SSA_target_total, T_p1, T_p2, alpha);

% Recompute size factors after calibration.
[T_size] = computeSizeFactors(r_p, T_p1, T_p2, T_r, T_L, T_shape);
D = 1.0; % diffusivity placeholder
G = D .* (pi .* (T_r.^2) .* T_shape) ./ max(T_L, eps);

% Assemble incidence and diagonal conductance matrix.
[B, K] = assembleIncidence(N, T_p1, T_p2, G);

% Package output.
PNM = struct();
PNM.meta = struct('R', geom.R_sample, 'N', N, 'alpha_rt', alpha, ...
    'placement', placementMeta, 'coordProfile', coordProfile);
PNM.P = struct('coords', coords, 'r_p', r_p, 'is_surface', is_surface, ...
    'class_id', class_id);
PNM.T = struct('pore1', T_p1, 'pore2', T_p2, 'r_t', T_r, 'L_geom', T_L, ...
    'shape_factor', T_shape, 'size_factor_axial', T_size);
PNM.graph = struct('B', B, 'K', spdiags(G, 0, M, M));
PNM.targets = geom.targets;
PNM.templateGeom = geom;
PNM.qc = validatePNM(PNM, geom);

if isfield(opts, 'save_file') && ~isempty(opts.save_file)
    save(opts.save_file, 'PNM', '-v7.3');
end
end
