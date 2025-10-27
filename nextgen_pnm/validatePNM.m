function QC = validatePNM(PNM, geom)
% validatePNM  Extended diagnostics for the CaCO3-specific network.
%
%   Outputs include:
%       - Degree statistics per class compared with requested ranges
%       - Porosity / SSA estimates against targets
%       - Throat radius cap violations
%
coords = PNM.P.coords;
class_id = PNM.P.class_id;
T = PNM.T;
R = geom.R_sample;

G = graph(T.pore1, T.pore2);
deg = degree(G);
QC.degree = struct('overall_mean', mean(deg), 'overall_std', std(deg));

classes = geom.pore_classes;
class_stats = struct('name', {}, 'mean_deg', {}, 'std_deg', {}, ...
    'requested', {}, 'count', {});
for c = 1:numel(classes)
    ids = find(class_id == c);
    if isempty(ids)
        class_stats(c).name = classes(c).name;
        class_stats(c).mean_deg = NaN;
        class_stats(c).std_deg = NaN;
        class_stats(c).requested = classes(c).coord_range;
        class_stats(c).count = 0;
        continue;
    end
    class_stats(c).name = classes(c).name;
    class_stats(c).mean_deg = mean(deg(ids));
    class_stats(c).std_deg = std(deg(ids));
    class_stats(c).requested = classes(c).coord_range;
    class_stats(c).count = numel(ids);
end
QC.degree.per_class = class_stats;

% Porosity / SSA check
[phi_est, SSA_est] = estimatePhiSSA(R, PNM.P.r_p, PNM.T.r_t, PNM.T.L_geom);
QC.phi_est = phi_est;
QC.SSA_est = SSA_est;
QC.target_phi = geom.phi_target;
QC.target_SSA = geom.SSA_target_total;

% Throat violations relative to alpha cap
alpha = geom.alpha_rt;
violations = PNM.T.r_t > alpha * min(PNM.P.r_p(PNM.T.pore1), PNM.P.r_p(PNM.T.pore2));
QC.throat_cap_violation_fraction = nnz(violations) / numel(violations);

% Surface pore count for sanity
radial = vecnorm(coords, 2, 2);
QC.surface_pore_fraction = nnz(radial > 0.95 * R) / numel(radial);
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
