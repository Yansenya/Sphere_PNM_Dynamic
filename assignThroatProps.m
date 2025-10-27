function [rt, shape, Lt] = assignThroatProps(P, p1, p2, Lgeom, rp, th, alpha, L_min)
% 从占位分布赋予喉参数并做物理一致性修正
M = numel(p1);
% 估算 E[rp] 以修正 rt 的 lognormal 均值
Erp = mean(rp);

% 先采样形状因子
shape = betarnd(th.shape_a, th.shape_b, [M,1]);
shape = max(min(shape,0.999), 0.05); % 防止极端

% 采样 rt（相对 Erp 的偏移）
rt = lognrnd( log(Erp) + th.rt_mu_offset, th.rt_sigma_log, [M,1]);
rt = min(max(rt, th.rt_minmax(1)), th.rt_minmax(2));

% 约束：rt <= alpha * min(rp_i, rp_j)
minrp = min(rp(p1), rp(p2));
rt = min(rt, alpha * minrp);

% 采样 Lt 比例并转为有效长度
beta = max(min(normrnd(th.Lt_ratio_mu, th.Lt_ratio_sigma, [M,1]), 1.1), 0.7);
Lt = max(beta .* Lgeom, L_min);
end
