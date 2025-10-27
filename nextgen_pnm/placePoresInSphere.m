function [coords, class_id, is_surface, meta] = placePoresInSphere(geom, opts)
% placePoresInSphere  Sample pore centres with class-specific radial spans.
%
%   Compared to the legacy jittered-FCC placement this routine:
%     * honours the class mixture delivered by materialToPNMParams,
%     * biases micropores towards the interior and macro pores towards the
%       periphery, and
%     * returns per-class placement metadata for diagnostics.
%
%   Inputs
%       geom.class_id       vector of class labels (aligned with pore radii)
%       geom.pore_classes   class templates with .radial_span entries
%       geom.R_sample       pellet radius
%
%   Outputs mirror the expectations of buildSphericalPNM.

class_id = geom.class_id(:);
classes = geom.pore_classes;
R = geom.R_sample;
N = numel(class_id);
coords = zeros(N, 3);

for i = 1:N
    c = class_id(i);
    span = classes(c).radial_span;
    r_min = span(1) * R;
    r_max = span(2) * R;
    if r_min == r_max
        r = r_min;
    else
        xi = rand();
        r = (r_min^3 + (r_max^3 - r_min^3) * xi)^(1/3);
    end
    u = randn(1,3);
    u = u / norm(u);
    coords(i, :) = r * u;
end

radial = vecnorm(coords, 2, 2);
is_surface = radial > 0.98 * R;

counts = accumarray(class_id, 1, [numel(classes), 1]);
meta = struct('counts', counts, 'radial_mean', accumarray(class_id, radial, ...
    [numel(classes), 1], @mean, NaN));
end
