function PNM = buildSphericalPNM_varN(opts, varargin)
% buildSphericalPNM_varN - 可变孔数N的球形PNM生成器（支持SSA目标开关与进度打印）
% 目标：
%   - 孔径分布近似正态（用对数正态实现），众数∈[1,10] nm（默认 5 nm）
%   - 以孔隙率 phi 为主，通过自适应 N 命中 phi（相对误差≤opts.phi_tol，默认1%）
%   - 可选：再小幅缩放（不改 N）以拉近 SSA 目标（由开关 opts.use_SSA 控制）
%
% 依赖：placePoresInSphere, connectPores, assignThroatProps, computeSizeFactors,
%       assembleIncidence, validatePNM （与你现有版本一致）

%% 参数与默认
p = inputParser;
addParameter(p, 'param_bank', '', @ischar);
parse(p, varargin{:});
param_bank_file = p.Results.param_bank;

% 几何/连通
if ~isfield(opts,'R'), opts.R = 250e-6; end
if ~isfield(opts,'alpha_rt'), opts.alpha_rt = 0.7; end
if ~isfield(opts,'lattice'), opts.lattice = 'fcc'; end
if ~isfield(opts,'jitter_sigma'), opts.jitter_sigma = []; end
if ~isfield(opts,'k_candidate'), opts.k_candidate = 12; end
if ~isfield(opts,'L_min'), opts.L_min = 5e-9; end   % 下限5 nm，避免 L=0
if ~isfield(opts,'outdir'), opts.outdir = 'fig'; end
if ~exist(opts.outdir,'dir'), mkdir(opts.outdir); end

% 可变N策略
if ~isfield(opts,'N_init'),    opts.N_init = 10000; end
if ~isfield(opts,'N_min'),     opts.N_min  = 200;  end
if ~isfield(opts,'N_max'),     opts.N_max  = 2e5;  end
if ~isfield(opts,'phi_tol'),   opts.phi_tol = 0.01; end     % 1% 相对误差
if ~isfield(opts,'maxNIter'),  opts.maxNIter = 12; end

% 孔径“近似正态”：用对数正态的众数设置
if ~isfield(opts,'rp_mode_nm'),   opts.rp_mode_nm = 500; end      % 5 nm
if ~isfield(opts,'rp_sigma_log'), opts.rp_sigma_log = 0.35; end % 对数标准差
if ~isfield(opts,'rp_min_nm'),    opts.rp_min_nm  = 0.8; end
if ~isfield(opts,'rp_max_nm'),    opts.rp_max_nm  = 50; end

% 喉参数（经验相对端孔法）
if ~isfield(opts,'rt_sigma_log'),  opts.rt_sigma_log = 0.4; end
if ~isfield(opts,'rt_scale_mean'), opts.rt_scale_mean = 0.35; end
if ~isfield(opts,'rt_scale_std'),  opts.rt_scale_std  = 0.08; end
if ~isfield(opts,'shape_a'),       opts.shape_a = 2.5; end
if ~isfield(opts,'shape_b'),       opts.shape_b = 5.0; end
if ~isfield(opts,'Lt_ratio_mu'),   opts.Lt_ratio_mu = 0.9; end
if ~isfield(opts,'Lt_ratio_sigma'),opts.Lt_ratio_sigma = 0.05; end

% 打印与随机
if ~isfield(opts,'verbose'),       opts.verbose = true; end
if ~isfield(opts,'print_every'),   opts.print_every = 1; end   % 每次 N 迭代都打印
if ~isfield(opts,'rng_seed'),      opts.rng_seed = 42; end
rng(opts.rng_seed, 'twister');

%% 读取参数库
if isempty(param_bank_file)
    error('必须提供 param_bank.mat（至少包含 global.phi_target/SSA_target）');
end
S = load(param_bank_file);
pb = S.param_bank;

phi_target  = pb.global.phi_target(1);
SSA_target  = pb.global.SSA_target(1);
if isnan(SSA_target), SSA_target = []; end
R           = opts.R;

% SSA 开关逻辑：优先遵从 opts.use_SSA；否则按是否有 SSA_target 自动决定
if ~isfield(opts,'use_SSA')
    opts.use_SSA = ~isempty(SSA_target);
end

%% 根据众数确定对数正态的 mu
% lognormal: mode = exp(mu - sigma^2) => mu = ln(mode) + sigma^2
rp_mode_m   = max(1, min(10, opts.rp_mode_nm)) * 1e-9;  % [1,10]nm
mu_log_rp   = log(rp_mode_m) + opts.rp_sigma_log^2;
rp_min      = opts.rp_min_nm * 1e-9;
rp_max      = opts.rp_max_nm * 1e-9;

%% 自适应 N：命中 φ（相对误差≤phi_tol）
alpha = opts.alpha_rt;
kCand = opts.k_candidate;
N     = opts.N_init;

best = struct('phi_err',Inf);
t_total = tic;

for itN = 1:opts.maxNIter
    t_iter = tic;

    % 1) 放置孔
    [P_coords, is_surface, lattice_s, jitter_sigma] = placePoresInSphere(R, N, opts);
    N_real = size(P_coords,1);

    % 2) 配位与连接
    z_vals = 3:7;  z_pmf = [0.1,0.25,0.3,0.25,0.1]; z_pmf = z_pmf/sum(z_pmf);
    z_target = randsample(z_vals, N_real, true, z_pmf);
    [T_p1, T_p2, L_geom] = connectPores(P_coords, z_target, kCand);
    assert(all(L_geom>0), '出现喉长<=0（检查 L_min 或连通器）');

    % 3) 孔半径（lognormal，众数锁定）
    rp = lognrnd(mu_log_rp, opts.rp_sigma_log, [N_real,1]);
    rp = min(max(rp, rp_min), rp_max);

    % 4) 喉属性（相对端孔）
    th = struct();
    th.rt_sigma_log   = opts.rt_sigma_log;
    th.shape_a        = opts.shape_a;
    th.shape_b        = opts.shape_b;
    th.Lt_ratio_mu    = opts.Lt_ratio_mu;
    th.Lt_ratio_sigma = opts.Lt_ratio_sigma;

    [T_rt, T_shape, T_L] = assignThroatProps_rel(P_coords, T_p1, T_p2, L_geom, rp, ...
        opts.rt_scale_mean, opts.rt_scale_std, th, alpha, opts.L_min);

    % 5) 估算 φ/SSA
    [phi_est, SSA_est] = estimatePhiSSA(R, rp, T_rt, T_L);
    rel_err_phi = abs(phi_est - phi_target)/max(phi_target,1e-16);

    % 打印进度
    if opts.verbose && (mod(itN, opts.print_every)==0 || itN==1)
        if isempty(SSA_target) || ~opts.use_SSA
            fprintf('[%2d/%2d] N=%7d | phi=%.6g (target=%.6g, rel.err=%.3g) | dt=%.2fs\n', ...
                itN, opts.maxNIter, N_real, phi_est, phi_target, rel_err_phi, toc(t_iter));
        else
            fprintf('[%2d/%2d] N=%7d | phi=%.6g (%.3g) | SSA=%.6g%s | dt=%.2fs\n', ...
                itN, opts.maxNIter, N_real, phi_est, rel_err_phi, SSA_est, ...
                sprintf(' (target=%.6g)', SSA_target), toc(t_iter));
        end
    end

    % 更新最佳
    if rel_err_phi < best.phi_err
        best.phi_err = rel_err_phi;
        best.data = struct('rp',rp,'T_rt',T_rt,'T_L',T_L,'T_shape',T_shape,...
            'P_coords',P_coords,'is_surface',is_surface,'T_p1',T_p1,'T_p2',T_p2,...
            'lattice_s',lattice_s,'jitter_sigma',jitter_sigma,...
            'phi_est',phi_est,'SSA_est',SSA_est,'N_real',N_real);
        if rel_err_phi <= opts.phi_tol
            break;  % 达标，停止 N 调整
        end
    end

    % 6) 调整 N（近似 φ∝N）
    if phi_est > 0
        N_new = round(N_real * phi_target / phi_est);
    else
        N_new = round(1.2*N_real);
    end
    N_new = max(opts.N_min, min(opts.N_max, N_new));
    if N_new == N
        N_new = N + sign(phi_target - phi_est)*max(5, round(0.02*N));
        N_new = max(opts.N_min, min(opts.N_max, N_new));
    end
    if opts.verbose
        fprintf('   ↳ adjust N: %d → %d (phi_est=%.3g, target=%.3g)\n', N, N_new, phi_est, phi_target);
    end
    N = N_new;
end

% 用最佳解继续
rp        = best.data.rp;
T_rt      = best.data.T_rt;
T_L       = best.data.T_L;
T_shape   = best.data.T_shape;
P_coords  = best.data.P_coords;
is_surface= best.data.is_surface;
T_p1      = best.data.T_p1;
T_p2      = best.data.T_p2;
lattice_s = best.data.lattice_s;
jitter_sigma = best.data.jitter_sigma;
phi_est   = best.data.phi_est;
SSA_est   = best.data.SSA_est;
N_real    = best.data.N_real;

% 7) （可选）SSA 微调（不改 N；优先保证 φ 误差≤phi_tol）
if opts.use_SSA && ~isempty(SSA_target) && isfinite(SSA_target) && SSA_target>0
    [sc_rp, sc_rt, rp, T_rt] = calibratePhiSSA_local(R, rp, T_rt, T_L, phi_target, SSA_target, ...
        'TolRelPhi', opts.phi_tol, 'TolRelSSA', 0.05, 'MaxIter', 40);
    [phi_est, SSA_est] = estimatePhiSSA(R, rp, T_rt, T_L);
    if opts.verbose
        fprintf('SSA tuning done: phi=%.6g (target=%.6g), SSA=%.6g (target=%.6g)\n', ...
            phi_est, phi_target, SSA_est, SSA_target);
    end
end

% 8) 导通与装配
D = 1.0;
T_size = computeSizeFactors(rp, T_p1, T_p2, T_rt, T_L, T_shape);
G = D .* (pi .* (T_rt.^2) .* T_shape) ./ T_L;

[B, K] = assembleIncidence(N_real, T_p1, T_p2, G);
QC = validatePNM(P_coords, T_p1, T_p2, rp, T_rt, R, (3:7), [], opts.alpha_rt);
[phi_check, SSA_check] = estimatePhiSSA(R, rp, T_rt, T_L);
QC.phi_est = phi_check;
QC.SSA_est = SSA_check;

% 9) 输出
PNM = struct();
PNM.meta = struct('R',R,'target_phi',phi_target,'SSA_target', iff(opts.use_SSA, SSA_target, []), ...
                  'lattice',opts.lattice,'lattice_s',lattice_s, ...
                  'jitter_sigma',jitter_sigma,'alpha_rt',opts.alpha_rt, ...
                  'k_candidate',kCand, 'N', N_real, 'use_SSA', logical(opts.use_SSA));
PNM.P = struct('coords',P_coords,'r_p',rp,'is_surface',is_surface);
PNM.T = struct('pore1',T_p1,'pore2',T_p2,'r_t',T_rt,'L_geom',T_L,'shape_factor',T_shape, ...
               'size_factor_axial',T_size);
PNM.graph = struct('B',B,'K',spdiags(K,0,numel(T_rt),numel(T_rt)),'L_builder','L_builder_placeholder');
PNM.qc = QC;

save('PNM.mat','PNM','-v7.3');
fprintf('Saved PNM.mat  |  N=%d, phi=%.6g (target=%.6g), SSA=%.6g%s  |  total %.2fs\n', ...
    N_real, phi_check, phi_target, SSA_check, ...
    iff(opts.use_SSA && ~isempty(SSA_target), sprintf(' (target=%.6g)', SSA_target), ''), ...
    toc(t_total));

end

%% ======== 辅助：喉属性（相对端孔法） ========
function [T_rt, T_shape, T_L] = assignThroatProps_rel(P_coords, T_p1, T_p2, L_geom, rp, ...
        rt_scale_mean, rt_scale_std, th, alpha, L_min)
M = numel(T_p1);
min_rp = min(rp(T_p1), rp(T_p2));
scale  = max(rt_scale_mean + rt_scale_std.*randn(M,1), 0);
T_rt   = scale .* min_rp;
T_rt   = min(T_rt, alpha*min_rp);
T_rt(T_rt<0) = 0;

T_shape = betarnd(th.shape_a, th.shape_b, [M,1]);

d_ij = vecnorm(P_coords(T_p1,:) - P_coords(T_p2,:), 2, 2);
ratio = max(th.Lt_ratio_mu + th.Lt_ratio_sigma.*randn(M,1), 0);
T_L   = max(L_min, ratio .* d_ij);
end

%% ======== 辅助：估计 φ / SSA ========
function [phi_est, SSA_est] = estimatePhiSSA(R, rp, rt, Lt)
V_sphere  = (4/3)*pi*R^3;
V_pores   = sum((4/3)*pi*rp.^3);
V_throats = sum(pi*rt.^2 .* Lt);
phi_est   = (V_pores + V_throats) / V_sphere;

A_pores   = sum(4*pi*rp.^2);
A_throats = sum(2*pi*rt .* Lt);
SSA_est   = (A_pores + A_throats) / V_sphere;
end

%% ======== 辅助：微调 φ / SSA（不改 N） ========
function [sc_rp, sc_rt, rp_new, rt_new] = calibratePhiSSA_local(R, rp, rt, Lt, phi_target, SSA_target, varargin)
p = inputParser;
addParameter(p, 'TolRelPhi', 0.01);
addParameter(p, 'TolRelSSA', 0.05);
addParameter(p, 'MaxIter', 40);
parse(p, varargin{:});
TolRelPhi = p.Results.TolRelPhi;
TolRelSSA = p.Results.TolRelSSA;
MaxIter   = p.Results.MaxIter;

rp_new = rp; rt_new = rt; sc_rp = 1.0; sc_rt = 1.0;

% φ：二分缩放 rp
f_phi = @(s) estimatePhiSSA(R, max(s*rp,0), rt_new, Lt);
s1 = bisection_target(@(s) f_phi(s), phi_target, TolRelPhi, MaxIter, 'phi');
if ~isnan(s1), rp_new = max(s1*rp,0); sc_rp = sc_rp*s1; end

% SSA：二分缩放 rt（可选）
if ~isempty(SSA_target) && isfinite(SSA_target) && SSA_target>0
    f_ssa = @(s) estimatePhiSSA(R, rp_new, max(s*rt,0), Lt);
    s2 = bisection_target(@(s) f_ssa(s), SSA_target, TolRelSSA, MaxIter, 'ssa');
    if ~isnan(s2), rt_new = max(s2*rt,0); sc_rt = sc_rt*s2; end
end
end

function s = bisection_target(fun_phi_ssa, target, tolRel, maxIter, which)
getv = @(x) (strcmp(which,'phi') * x(1) + strcmp(which,'ssa') * x(2));
v1 = getv(fun_phi_ssa(1.0));
if abs(v1 - target)/max(target,1e-16) <= tolRel, s=1.0; return; end

% bracket
if v1 > target
    s_lo = 0;    v_lo = getv(fun_phi_ssa(s_lo));
    s_hi = 1.0;  v_hi = v1;
    if v_lo > target, s = NaN; return; end
else
    s_lo = 1.0; v_lo = v1;
    s_hi = 2.0; v_hi = getv(fun_phi_ssa(s_hi));
    cnt=0;
    while v_hi < target && cnt<30
        s_lo = s_hi; v_lo = v_hi;
        s_hi = 2*s_hi; v_hi = getv(fun_phi_ssa(s_hi)); cnt=cnt+1;
        if s_hi > 1e6, s=NaN; return; end
    end
end

for it=1:maxIter
    s_mid = 0.5*(s_lo + s_hi);
    v_mid = getv(fun_phi_ssa(s_mid));
    if abs(v_mid - target)/max(target,1e-16) <= tolRel, s = s_mid; return; end
    if v_mid > target
        s_hi = s_mid;
    else
        s_lo = s_mid;
    end
end
s = 0.5*(s_lo+s_hi);
end

%% ======== 小工具 ========
function y = iff(cond, a, b)
if cond, y=a; else, y=b; end
end
