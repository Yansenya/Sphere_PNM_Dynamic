function [r_t, shape, L_eff] = assignThroatProps_pair(rp_i, rp_j, ci, cj, L_geom, alpha)
% assignThroatProps_pair  Shared throat property sampler for existing/new edges.
%
%   This helper mirrors the heuristics used in assignThroatProps but is exposed
%   as a standalone function so calcination updates can call it directly.

r_min = min(rp_i, rp_j);
key = 10*min(ci, cj) + max(ci, cj);
switch key
    case 11 % micro-micro
        scale = betarnd(3, 6);
        shape = 0.85;
    case 12 % micro-meso
        scale = betarnd(2.5, 4.5);
        shape = 0.8;
    case 13 % micro-macro
        scale = betarnd(2.0, 5.0);
        shape = 0.75;
    case 22 % meso-meso
        scale = betarnd(2.5, 3.5);
        shape = 0.7;
    case 23 % meso-macro
        scale = betarnd(2.0, 3.0);
        shape = 0.65;
    case 33 % macro-macro
        scale = betarnd(1.8, 2.2);
        shape = 0.6;
    otherwise
        scale = betarnd(2.5, 4.0);
        shape = 0.7;
end
r_t = min(alpha * r_min, scale * r_min);
L_eff = max(L_geom - 0.5*(rp_i + rp_j), 0.2*(rp_i + rp_j));
L_eff = max(L_eff, 1e-8);
end
