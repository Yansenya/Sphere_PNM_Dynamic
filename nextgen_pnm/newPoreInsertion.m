function [PNM, state] = newPoreInsertion(PNM, state, parentIdx, new_volume, reactionSpec)
% newPoreInsertion  Spawn CaO micropores around a reacting CaCO3 pore.
%
%   The new pores inherit micropore class (1) and are positioned along random
%   directions from the parent node while keeping inside the spherical pellet.
%   Throat properties are generated using assignThroatProps_pair so the new
%   connections remain consistent with the base network statistics.

if new_volume <= 0
    return;
end

maxNew = reactionSpec.max_new_pores_per_event;
mu = log((reactionSpec.new_pore_mean_nm*1e-9)^2 / ...
    sqrt((reactionSpec.new_pore_sigma_nm*1e-9)^2 + (reactionSpec.new_pore_mean_nm*1e-9)^2));
sigma = sqrt(log(1 + (reactionSpec.new_pore_sigma_nm*1e-9)^2 / ...
    (reactionSpec.new_pore_mean_nm*1e-9)^2));

candidate_r = lognrnd(mu, sigma, [maxNew, 1]);
volumes = 4/3*pi*candidate_r.^3;
total_vol = sum(volumes);
if total_vol == 0
    return;
end
scale = (new_volume / total_vol)^(1/3);
radii = candidate_r * scale;

N0 = numel(PNM.P.r_p);
parentCoord = state.coords(parentIdx, :);
parentNeighbors = unique([PNM.T.pore2(PNM.T.pore1 == parentIdx); ...
    PNM.T.pore1(PNM.T.pore2 == parentIdx)]);
parentNeighbors(parentNeighbors == parentIdx) = [];

for k = 1:numel(radii)
    new_idx = N0 + k;
    r_new = radii(k);
    direction = randn(1,3);
    direction = direction / norm(direction);
    offset = (state.r_p(parentIdx) + r_new) * (0.6 + 0.2*rand());
    coord_candidate = parentCoord + offset * direction;
    R = PNM.meta.R;
    attempts = 0;
    while norm(coord_candidate) + r_new > R && attempts < 10
        direction = randn(1,3);
        direction = direction / norm(direction);
        coord_candidate = parentCoord + offset * direction;
        attempts = attempts + 1;
    end
    coord_candidate = min(1, R / max(norm(coord_candidate), 1e-9)) * coord_candidate;
    PNM.P.coords = [PNM.P.coords; coord_candidate];
    PNM.P.r_p = [PNM.P.r_p; r_new];
    PNM.P.is_surface = [PNM.P.is_surface; norm(coord_candidate) > 0.98 * R];
    PNM.P.class_id = [PNM.P.class_id; 1];

    state.coords = [state.coords; coord_candidate];
    state.r_p = [state.r_p; r_new];
    state.is_surface = [state.is_surface; norm(coord_candidate) > 0.98 * R];
    state.class_id = [state.class_id; 1];
    state.caCO3_volume = [state.caCO3_volume; 0];
    state.caCO3_volume0 = [state.caCO3_volume0; 0];
    state.caCO3_initial = [state.caCO3_initial; 0];

    % Connect to parent and selected neighbours.
    connect_to = parentNeighbors;
    if isempty(connect_to)
        connect_to = parentIdx;
    else
        coords_existing = PNM.P.coords(connect_to, :);
        dists = vecnorm(coords_existing - coord_candidate, 2, 2);
        [~, order] = sort(dists);
        connect_to = connect_to(order(1:min(2, numel(order))));
        connect_to = unique([parentIdx; connect_to(:)]);
    end
    if k > 1
        siblings = (N0 + (1:k-1))';
        connect_to = unique([connect_to; siblings]);
    end

    for target = reshape(connect_to, 1, [])
        L_geom = norm(PNM.P.coords(target, :) - coord_candidate);
        [r_t, shape, L_eff] = assignThroatProps_pair(PNM.P.r_p(target), r_new, ...
            PNM.P.class_id(target), 1, L_geom, state.alpha);
        PNM.T.pore1 = [PNM.T.pore1; target];
        PNM.T.pore2 = [PNM.T.pore2; new_idx];
        PNM.T.r_t = [PNM.T.r_t; r_t];
        PNM.T.L_geom = [PNM.T.L_geom; L_eff];
        PNM.T.shape_factor = [PNM.T.shape_factor; shape];
        PNM.T.size_factor_axial = [PNM.T.size_factor_axial; 0]; % placeholder
    end
end

% Recompute size factors and graph matrices after batch insertion.
M = numel(PNM.T.pore1);
PNM.T.size_factor_axial = computeSizeFactors(PNM.P.r_p, PNM.T.pore1, PNM.T.pore2, ...
    PNM.T.r_t, PNM.T.L_geom, PNM.T.shape_factor);
G = (pi .* (PNM.T.r_t.^2) .* PNM.T.shape_factor) ./ max(PNM.T.L_geom, eps);
[B, Kvec] = assembleIncidence(numel(PNM.P.r_p), PNM.T.pore1, PNM.T.pore2, G);
PNM.graph.B = B;
PNM.graph.K = spdiags(Kvec, 0, M, M);
PNM.meta.N = numel(PNM.P.r_p);
PNM.targets.N_per_class(1) = PNM.targets.N_per_class(1) + numel(radii);
if isfield(PNM, 'templateGeom') && isfield(PNM.templateGeom, 'targets') && ...
        numel(PNM.templateGeom.targets.N_per_class) >= 1
    PNM.templateGeom.targets.N_per_class(1) = ...
        PNM.templateGeom.targets.N_per_class(1) + numel(radii);
end
end
