function [T_r, T_shape, T_L] = assignThroatProps(~, T_p1, T_p2, L_geom, r_p, class_id, geom)
% assignThroatProps  Class-sensitive throat generation leveraging pair helper.
M = numel(T_p1);
T_r = zeros(M, 1);
T_shape = zeros(M, 1);
T_L = zeros(M, 1);
alpha = geom.alpha_rt;
for e = 1:M
    i = T_p1(e);
    j = T_p2(e);
    [T_r(e), T_shape(e), T_L(e)] = assignThroatProps_pair(r_p(i), r_p(j), ...
        class_id(i), class_id(j), L_geom(e), alpha);
end
end
