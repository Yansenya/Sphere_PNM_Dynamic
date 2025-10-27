function generate_param_bank_placeholder()
% 生成与 Python 占位模型等价的 out/param_bank.mat（无需 Python）
% 字段名、单位、默认值 与你给的 DEFAULT_CONFIG 完全一致

% --- 目录 ---
if ~exist('out','dir'), mkdir('out'); end

% --- DEFAULT_CONFIG（与 Python 保持一致）---
cfg = struct();
cfg.phi_target      = 0.009;
% cfg.SSA_target      = [];          % None -> 保存为 NaN
cfg.SSA_target      = 652600;   
cfg.alpha_rt_max    = 0.7;
cfg.R               = 250e-6;
cfg.N_target        = 50000;

% pore radius lognormal
cfg.rp_mu_log       = log(5e-6);
cfg.rp_sigma_log    = 0.4;
cfg.rp_min          = 0.5e-6;
cfg.rp_max          = 30e-6;

% throat radius lognormal（注意：这是“相对 ln(E[rp]) 的偏移”，构网器内再处理）
cfg.rt_mu_log_offset= log(0.5);
cfg.rt_sigma_log    = 0.45;
cfg.rt_min          = 0.2e-6;
cfg.rt_max          = 10e-6;

% Lt 模型（长度为比例）
cfg.Lt_ratio_mu     = 0.9;
cfg.Lt_ratio_sigma  = 0.05;
cfg.Lt_min_factor   = 0.0;   % 绝对 L_min 在 MATLAB 构网器中另设

% degree PMF z in [3..7]
cfg.z_pmf           = [0.1, 0.25, 0.3, 0.25, 0.1];

% shape factor Beta(a,b)
cfg.shape_a         = 2.5;
cfg.shape_b         = 5.0;

% tortuosity proxy
cfg.tau_mu          = 1.2;
cfg.tau_sigma       = 0.15;

% 其他
cfg.seed            = 42;
cfg.note            = 'placeholder-model-config';

% --- 组装与 Python savemat 等价的数据结构 ---
param_bank = struct();

% global
param_bank.global = struct();
param_bank.global.phi_target   = double(cfg.phi_target);
param_bank.global.SSA_target   = double(ifelse(isempty(cfg.SSA_target), nan, cfg.SSA_target));
param_bank.global.alpha_rt_max = double(cfg.alpha_rt_max);
param_bank.global.R            = double(cfg.R);
param_bank.global.N_target     = double(cfg.N_target);

% pore
param_bank.pore = struct();
param_bank.pore.rp_dist   = {'logn'};                    % 与 numpy object 等价
param_bank.pore.rp_mu     = double(cfg.rp_mu_log);
param_bank.pore.rp_sigma  = double(cfg.rp_sigma_log);
param_bank.pore.rp_minmax = double([cfg.rp_min, cfg.rp_max]);
param_bank.pore.z_pmf     = double(cfg.z_pmf(:)).';      % 行向量

% throat
param_bank.throat = struct();
param_bank.throat.rt_dist        = {'logn'};
param_bank.throat.rt_mu          = double(cfg.rt_mu_log_offset); % 相对 ln(E[rp]) 偏移
param_bank.throat.rt_sigma       = double(cfg.rt_sigma_log);
param_bank.throat.rt_minmax      = double([cfg.rt_min, cfg.rt_max]);
param_bank.throat.Lt_model       = {'ratio'};
param_bank.throat.Lt_ratio_mu    = double(cfg.Lt_ratio_mu);
param_bank.throat.Lt_ratio_sigma = double(cfg.Lt_ratio_sigma);
param_bank.throat.Lt_min_factor  = double(cfg.Lt_min_factor);
param_bank.throat.shape_a        = double(cfg.shape_a);
param_bank.throat.shape_b        = double(cfg.shape_b);

% flow / tortuosity
param_bank.flow = struct();
param_bank.flow.tau_proxy = {'dist'};
param_bank.flow.tau_mu    = double(cfg.tau_mu);
param_bank.flow.tau_sigma = double(cfg.tau_sigma);

% seed / note
param_bank.seed = int32(cfg.seed);
param_bank.note = {cfg.note};

% --- 保存 ---
save(fullfile('out','param_bank.mat'), 'param_bank', '-v7');
fprintf('Saved out/param_bank.mat (placeholder, identical defaults to Python).\n');

end

% 小工具：三元
function y = ifelse(cond, a, b)
if cond, y = a; else, y = b; end
end
