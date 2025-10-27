function geom = materialToPNMParams(materialSpec, opts)
% materialToPNMParams  Convert BET/porosity measurements to geometric targets.
%
%   The function embeds several modelling assumptions that differ from the
%   original implementation:
%     1. The representative pellet mass is derived from the provided radius and
%        CaCO3 density (unless overridden) so SSA/PV values convert to absolute
%        areas/volumes.
%     2. SSA fractions prescribe pore-class occupancy.  Each class receives a
%        log-normal radius distribution anchored to literature ranges.
%     3. Coordination and throat scaling templates are bundled for later use,
%        ensuring new pores created during calcination inherit compatible
%        statistics.
%
%   Output fields of `geom` include
%       .R_sample               [m] pellet radius
%       .sample_volume          [m^3]
%       .sample_mass_g          [g]
%       .phi_target             [-] porosity target
%       .SSA_target_total       [m^2]
%       .targets                struct summarising all key targets
%       .pore_classes           struct array with class meta data
%       .pore_radii             [N x 1] sampled radii (sorted by class order)
%       .class_id               [N x 1] integer class label (1:micro,2:meso,3:macro)
%       .coordProfile           coordination config per class
%       .alpha_rt               throat radius cap multiplier
%
arguments
    materialSpec struct
    opts struct
end

R = materialSpec.pellet_radius_um * 1e-6;
geom.R_sample = R;
geom.alpha_rt = opts.alpha_rt;
geom.N_total = opts.N_target;

if isfield(materialSpec, 'solid_density_kg_per_m3')
    rho = materialSpec.solid_density_kg_per_m3;
else
    rho = 2710; % kg/m^3 for CaCO3
end
vol = 4/3*pi*R^3;
geom.sample_volume = vol;

if isfield(materialSpec, 'mass_scaling_g')
    mass_g = materialSpec.mass_scaling_g;
else
    mass_g = rho * vol * 1e3; % convert kg to g
end
geom.sample_mass_g = mass_g;

SSA_total = materialSpec.SSA_total_m2_per_g * mass_g;
SSA_micro = materialSpec.SSA_micro_m2_per_g * mass_g;
SSA_meso_macro = materialSpec.SSA_meso_macro_m2_per_g * mass_g;

PV_total = materialSpec.pore_volume_cm3_per_g * mass_g * 1e-6; % cm^3→m^3
phi_target = PV_total / vol;
geom.phi_target = phi_target;
geom.SSA_target_total = SSA_total;
geom.targets = struct('SSA_total', SSA_total, 'SSA_micro', SSA_micro, ...
    'SSA_meso_macro', SSA_meso_macro, 'porosity', phi_target, ...
    'mass_g', mass_g, 'volume_m3', vol);

% Split the meso/macro SSA evenly unless the user provides a ratio.
if isfield(materialSpec, 'meso_fraction')
    meso_frac = materialSpec.meso_fraction;
else
    meso_frac = 0.7; % assume 70% of meso+macro SSA belongs to mesopores
end
SSA_meso = SSA_meso_macro * meso_frac;
SSA_macro = SSA_meso_macro * (1 - meso_frac);

% Define pore class templates (radius statistics from literature ranges).
classes = struct('name', {}, 'ssa', {}, 'radius_mean_nm', {}, ...
    'radius_sigma_nm', {}, 'coord_range', {}, 'radial_span', {});
classes(1).name = 'micropore';
classes(1).ssa = SSA_micro;
classes(1).radius_mean_nm = materialSpec.mean_pore_radius_nm; % anchored to adsorption mean
classes(1).radius_sigma_nm = 2.5; % narrow distribution
classes(1).coord_range = [8 10];
classes(1).radial_span = [0.0 0.7];

classes(2).name = 'mesopore';
classes(2).ssa = SSA_meso;
classes(2).radius_mean_nm = 40; % emphasise 10-100 nm window
classes(2).radius_sigma_nm = 12;
classes(2).coord_range = [6 8];
classes(2).radial_span = [0.5 0.9];

classes(3).name = 'macropore';
classes(3).ssa = SSA_macro;
classes(3).radius_mean_nm = 150;
classes(3).radius_sigma_nm = 45;
classes(3).coord_range = [4 6];
classes(3).radial_span = [0.8 1.0];

% Determine pore count per class via SSA weighting.
SSA_vec = [classes.ssa];
if sum(SSA_vec) <= 0
    error('SSA values must be positive.');
end
class_fraction = SSA_vec / sum(SSA_vec);
N_per_class = max(1, round(class_fraction * geom.N_total));
% adjust to ensure total matches target
N_correction = geom.N_total - sum(N_per_class);
[~, order] = sort(class_fraction, 'descend');
idx = 1;
while N_correction ~= 0
    c = order(idx);
    if N_correction > 0
        N_per_class(c) = N_per_class(c) + 1;
        N_correction = N_correction - 1;
    else
        if N_per_class(c) > 1
            N_per_class(c) = N_per_class(c) - 1;
            N_correction = N_correction + 1;
        end
    end
    idx = idx + 1;
    if idx > numel(order), idx = 1; end
end

class_id = zeros(sum(N_per_class), 1);
r_p = zeros(sum(N_per_class), 1);
ptr = 0;
for c = 1:numel(classes)
    ids = (1:N_per_class(c)) + ptr;
    class_id(ids) = c;
    mu = log((classes(c).radius_mean_nm*1e-9)^2 / sqrt((classes(c).radius_sigma_nm*1e-9)^2 + (classes(c).radius_mean_nm*1e-9)^2));
    sigma = sqrt(log(1 + (classes(c).radius_sigma_nm*1e-9)^2 / (classes(c).radius_mean_nm*1e-9)^2));
    samples = lognrnd(mu, sigma, [N_per_class(c), 1]);
    % clip to literature-supported band to avoid unrealistic pores
    lower = (classes(c).radius_mean_nm - 3*classes(c).radius_sigma_nm) * 1e-9;
    upper = (classes(c).radius_mean_nm + 3*classes(c).radius_sigma_nm) * 1e-9;
    samples = min(max(samples, max(lower, 1e-9)), min(upper, 0.5e-6));
    r_p(ids) = samples;
    ptr = ptr + N_per_class(c);
end

perm = randperm(numel(r_p));
geom.pore_classes = classes;
geom.pore_radii = r_p(perm);
geom.class_id = class_id(perm);
geom.coordProfile = struct('class', {classes.name}, 'target_range', {classes.coord_range});
geom.targets.N_per_class = N_per_class;
geom.targets.pore_classes = classes;
end
