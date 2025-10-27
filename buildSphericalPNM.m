function PNM = buildSphericalPNM(opts, varargin)
% buildSphericalPNM - 主入口：读取参数→放置孔→连接→赋值→校准→装配→可视化/QC
% 用法：
%   PNM = buildSphericalPNM(opts, 'param_bank','out/param_bank.mat')

%% 参数与默认
p = inputParser;
addParameter(p, 'param_bank', '', @ischar);
parse(p, varargin{:});
param_bank_file = p.Results.param_bank;

if ~isfield(opts,'R'), opts.R = 250e-6; end
if ~isfield(opts,'N_target'), opts.N_target = 5e3; end
if ~isfield(opts,'alpha_rt'), opts.alpha_rt = 0.7; end
if ~isfield(opts,'lattice'), opts.lattice = 'fcc'; end
if ~isfield(opts,'jitter_sigma'), opts.jitter_sigma = []; end
if ~isfield(opts,'k_candidate'), opts.k_candidate = 12; end
if ~isfield(opts,'L_min'), opts.L_min = 0.2e-6; end
if ~isfield(opts,'outdir'), opts.outdir = 'fig'; end
if ~exist(opts.outdir,'dir'), mkdir(opts.outdir); end

%% 读取参数库（Python 生成的 .mat）
if ~isempty(param_bank_file)
    S = load(param_bank_file);
    pb = S.param_bank;
else
    error('必须提供 param_bank.mat（由 Python 侧生成）');
end

% 拆解
PB.global     = pb.global;
PB.pore       = pb.pore;
PB.throat     = pb.throat;
PB.flow       = pb.flow;

phi_target    = PB.global.phi_target(1);
if ~isnan(PB.global.SSA_target(1))
    SSA_target = PB.global.SSA_target(1);
else
    SSA_target = [];
end
R     = opts.R;
Ntar  = opts.N_target;
alpha = opts.alpha_rt;

%% 放置孔：抖动 FCC + 球形裁剪
[P_coords, is_surface, lattice_s, jitter_sigma] = placePoresInSphere(R, Ntar, opts);
N = size(P_coords,1);

%% 配位数目标分配（z ~ pmf）
z_vals = 3:7;
z_pmf  = PB.pore.z_pmf(:)'; z_pmf = z_pmf/sum(z_pmf);
z_target = randsample(z_vals, N, true, z_pmf);

%% 候选邻居 + 连接（度受限 + MST 保证连通）
kCand = opts.k_candidate;
[T_p1, T_p2, L_geom] = connectPores(P_coords, z_target, kCand);

%% 赋 pore 半径（lognormal 占位）
rp_mu_log    = PB.pore.rp_mu(1);
rp_sigma_log = PB.pore.rp_sigma(1);
rp_minmax    = PB.pore.rp_minmax(1,:);
rng('default');
rp = lognrnd(rp_mu_log, rp_sigma_log, [N,1]);
rp = min(max(rp, rp_minmax(1)), rp_minmax(2));

%% 赋 throat 半径/形状因子/有效长度（并做一致性修正）
th = struct();
th.rt_mu_offset = PB.throat.rt_mu(1);  % 相对于 ln(E[rp])
th.rt_sigma_log = PB.throat.rt_sigma(1);
th.rt_minmax    = PB.throat.rt_minmax(1,:);
th.shape_a      = PB.throat.shape_a(1);
th.shape_b      = PB.throat.shape_b(1);
th.Lt_ratio_mu  = PB.throat.Lt_ratio_mu(1);
th.Lt_ratio_sigma = PB.throat.Lt_ratio_sigma(1);

[T_rt, T_shape, T_L] = assignThroatProps(P_coords, T_p1, T_p2, L_geom, rp, th, alpha, opts.L_min);

M = numel(T_rt);

%% 初步尺寸因子与导通率（扩散型）
[T_size] = computeSizeFactors(rp, T_p1, T_p2, T_rt, T_L, T_shape);
% 简化：扩散导通率 G = D * A / L; D 可当作 1 占位，A ~ pi*rt^2*shape
D = 1.0; 
G = D .* (pi .* (T_rt.^2) .* T_shape) ./ max(T_L, eps);

%% 初估孔隙率/SSA 并校准
[phi0, SSA0] = estimatePhiSSA(R, rp, T_rt, T_L);
% [sc_rp, sc_rt, rp_new, rt_new] = calibratePorositySSA(R, rp, T_rt, T_L, phi_target, SSA_target);
[sc_rp, sc_rt, rp_new, rt_new] = calibratePorositySSA(R, rp, T_rt, T_L, phi_target, SSA_target, T_p1, T_p2, alpha);
rp = rp_new; T_rt = rt_new;
[T_size] = computeSizeFactors(rp, T_p1, T_p2, T_rt, T_L, T_shape);
G = D .* (pi .* (T_rt.^2) .* T_shape) ./ max(T_L, eps);

%% 装配稀疏矩阵
[B, K] = assembleIncidence(N, T_p1, T_p2, G);
% 提供 L_builder 句柄名（以便保存后读取）
L_builder_name = 'L_builder_placeholder';

%% 可视化
% vizPNM(P_coords, rp, T_p1, T_p2, T_rt, R, opts.outdir);

%% QC
QC = validatePNM(P_coords, T_p1, T_p2, rp, T_rt, R, z_vals, z_target, alpha);
[phi_est, SSA_est] = estimatePhiSSA(R, rp, T_rt, T_L);
QC.phi_est = phi_est;
QC.SSA_est = SSA_est;

%% 组装输出结构
PNM = struct();
PNM.meta = struct('R',R,'target_phi',phi_target,'SSA_target',SSA_target, ...
                  'lattice',opts.lattice,'lattice_s',lattice_s, ...
                  'jitter_sigma',jitter_sigma,'alpha_rt',alpha, ...
                  'k_candidate',kCand);
PNM.P = struct('coords',P_coords,'r_p',rp,'is_surface',is_surface);
PNM.T = struct('pore1',T_p1,'pore2',T_p2,'r_t',T_rt,'L_geom',T_L,'shape_factor',T_shape, ...
               'size_factor_axial',T_size);
PNM.graph = struct('B',B,'K',spdiags(K,0,M,M),'L_builder',L_builder_name);
PNM.qc = QC;

save('PNM.mat','PNM','-v7.3');
fprintf('Saved PNM.mat\n');
end

function [phi_est, SSA_est] = estimatePhiSSA(R, rp, rt, Lt)
% 简易估计：孔体积+喉体积 / 球体积；表面积 ~ 孔表面积+喉侧面积
V_sphere = 4/3*pi*R^3;
V_pores  = sum(4/3*pi*rp.^3);
V_throats = sum(pi*rt.^2 .* Lt);
phi_est = (V_pores + V_throats) / V_sphere;

A_pores = sum(4*pi*rp.^2);
A_throats = sum(2*pi*rt .* Lt); % 侧面积近似
SSA_est = (A_pores + A_throats) / V_sphere;
end
