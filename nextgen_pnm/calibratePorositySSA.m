function [sc_rp, sc_rt, rp_new, rt_new] = calibratePorositySSA(R, rp, rt, Lt, phi_target, SSA_target, pore1, pore2, alpha)
% 两步缩放：先缩放 rp → 逼近 porosity；再缩放 rt → 微调 SSA（如提供）
% 当提供 pore1/pore2/alpha 时，在孔径缩放时同步缩放喉径，避免喉体积主导导致低孔隙率无法达到
if nargin < 7, pore1 = []; end
if nargin < 8, pore2 = []; end
if nargin < 9, alpha = []; end
% pore1/pore2/alpha 用于在缩放孔径时同步缩放喉道并维持 rt <= alpha*min(rp)

max_iter = 10;
rp_new = rp; rt_new = rt;
sc_rp = 1.0; sc_rt = 1.0;

[phi0, SSA0] = estimate(R, rp_new, rt_new, Lt);

for it=1:max_iter
    if ~isempty(phi_target)
        % 体积分数随 rp^3 主导，立方比例缩放
        if phi0 > 0
            sc_rp = (phi_target / phi0)^(1/3);
        else
            sc_rp = 1.0;
        end
        rp_new = rp_new * sc_rp;
        rt_new = rt_new * sc_rp;
        rt_new = enforceThroatConsistency(rt_new, rp_new, pore1, pore2, alpha);
    end

    [phi1, SSA1] = estimate(R, rp_new, rt_new, Lt);

    if ~isempty(SSA_target) && ~isnan(SSA_target)
        % 表面积 ~ rp^2 主导，平方比例缩放 throat 半径近似
        if SSA1 > 0
            sc_rt = (SSA_target / SSA1)^(1/2);
        else
            sc_rt = 1.0;
        end
        rt_new = rt_new * sc_rt;
        rt_new = enforceThroatConsistency(rt_new, rp_new, pore1, pore2, alpha);
    end

    [phi0, SSA0] = estimate(R, rp_new, rt_new, Lt);

    % 收敛检查（仅针对 phi，SSA 可选）
    if abs(phi0 - phi_target) < 0.005
        break;
    end
end

rt_new = enforceThroatConsistency(rt_new, rp_new, pore1, pore2, alpha);
end

function [phi_est, SSA_est] = estimate(R, rp, rt, Lt)
V_sphere = 4/3*pi*R^3;
V_pores  = sum(4/3*pi*rp.^3);
V_throats = sum(pi*rt.^2 .* Lt);
phi_est = (V_pores + V_throats) / V_sphere;

A_pores = sum(4*pi*rp.^2);
A_throats = sum(2*pi*rt .* Lt);
SSA_est = (A_pores + A_throats) / V_sphere;
end

function rt = enforceThroatConsistency(rt, rp, pore1, pore2, alpha)
if isempty(pore1) || isempty(pore2) || isempty(alpha)
    return;
end

min_rp = min(rp(pore1), rp(pore2));
rt = min(rt, alpha .* min_rp);
end

