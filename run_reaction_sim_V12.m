function run_reaction_sim_V12()
% ===== V11：V7 的稳健化版（奇异/病态修复）+ BET/孔体积/平均孔径记录 =====
% 新增记录：
%   BET 比表面积（m2/g）：S_internal / m_solid
%   孔体积（m3/g）：V_void_total / m_solid
%   平均孔径（m）：r_mean = 2 * V_void_total / S_internal  （几何等价直径定义）
%
% 关键稳健化（与旧V11相同）：
%   (1) Mc=diag(V_void_node)；(2) 死节点→临时Dirichlet锚到 C^n；
%   (3) 脱节分量→弱锚到 C^n；(4) 质量步对角预缩放。

clc;

%% 0) 载入 PNM
S   = load('PNM.mat');   PNM = S.PNM;

P = PNM.P;              % 节点：r_p, coords, is_surface(0/1)
Tn = PNM.T;             % 喉：pore1, pore2, r_t, L_geom, shape_factor
B = PNM.graph.B;        % 关联矩阵（M×N）
N = size(P.coords,1);
M = numel(Tn.pore1);
R_particle = PNM.meta.R;

%% 1) 常量与物性
R_gas   = 8.314462618;         % J/mol/K
sigmaSB = 5.670374419e-8;      % W/m^2/K^4

phys = struct();
phys.T_inf    = 1200;          % K
phys.ks_model = 'Arrhenius';
phys.k0       = 1.0e7;         % mol m^-2 s^-1（占位需标定）
phys.Ea       = 2.0e5;         % J/mol
phys.gammaSX  = 0.5;
phys.Dm       = 2.5e-4;        % m^2/s
phys.M_CO2    = 44.0095e-3;    % kg/mol
phys.Peq_A    = 7.079;  phys.Peq_B = 8308;
phys.deltaH   = 1.78e5;        % J/mol
phys.cp_CO2_gas = 1200;        % J/kg/K

micro = struct();
% micro.b_eff = 1.7;

therm = struct();
therm.enable_conduction = true;
therm.k_solid_eff   = 2.5;     % W/m/K
therm.A_solid_factor= 1.0;
therm.h_conv        = 120;
therm.emiss         = 0.8;

massBC = struct();
massBC.use_robin = true;
massBC.kg_mode   = 'const';
massBC.kg_const  = 0.01;       % m/s
massBC.C_inf     = 0.0;        % mol/m^3
massBC.Re        = 10;
massBC.Sc        = 0.8;
massBC.D_inf     = phys.Dm;

solid = struct();
solid.rho_CaCO3 = 2710;                 % kg/m^3
solid.M_CaCO3   = 100.0869e-3;          % kg/mol
solid.cp_CaCO3  = 900;                  % J/kg/K
solid.rho_CaO   = 3340;                 % kg/m^3
solid.M_CaO     = 56.0774e-3;           % kg/mol
solid.cp_CaO    = 800;                  % J/kg/K
solid.Omega_CaCO3 = solid.M_CaCO3 / solid.rho_CaCO3;
solid.Omega_CaO   = solid.M_CaO   / solid.rho_CaO;
solid.dOmega      = max(solid.Omega_CaCO3 - solid.Omega_CaO, 0);

%% 2) 初值与几何量
Vp  = (4/3)*pi*(P.r_p.^3);
assert(all(Tn.L_geom>0),'存在喉长<=0');
Vth = pi*(Tn.r_t.^2).*Tn.L_geom;
V_sphere = (4/3)*pi*R_particle^3;

% 初始有效反应表面积（孔面 + 半喉侧面）
S_pore    = 4*pi*(P.r_p.^2);
S_th_side = 2*pi*Tn.r_t.*Tn.L_geom;
halfS_to_end = accumarray([Tn.pore1; Tn.pore2], [S_th_side; S_th_side]/2, [N,1], @sum, 0);
S0 = S_pore + halfS_to_end;

% 初始 CaCO3 分布（按孔 + 半喉体积权重）
halfVth_to_end = accumarray([Tn.pore1; Tn.pore2], [Vth; Vth]/2, [N,1], @sum, 0);
w = Vp + halfVth_to_end;  sumw = sum(w);  assert(sumw>0,'权重和为0');
w = w / sumw;
V_void0  = sum(Vp)+sum(Vth);
V_solid0 = max(V_sphere - V_void0, 0);
V_CaCO3_0 = w * V_solid0;
V_CaCO3   = V_CaCO3_0;
X         = zeros(N,1);

% 初始温度与浓度
Tfield = phys.T_inf * ones(N,1);
p_eq0 = Peq_CO2_Pa(phys.T_inf, phys.Peq_A, phys.Peq_B);  assert(all(p_eq0>0));
C = (p_eq0 / (R_gas*phys.T_inf)) * ones(N,1);

% 表面节点
is_bnd = P.is_surface(:) ~= 0;
I_b = find(is_bnd);

% 整体转化率
Xtot = overall_conversion(V_CaCO3, V_CaCO3_0);

%% 3) 时间控制
ctrl = struct();
ctrl.dt_init   = 1e-3;
ctrl.dt_min    = 1e-6;
ctrl.dt_max    = 1000;
ctrl.grow      = 1.25;
ctrl.shrink    = 0.5;
ctrl.t_max     = 1e12;
ctrl.Xtot_stop = 0.99;
ctrl.dX_tol    = 0.01;

% ---- NEW: initial-stage limiter ----
ctrl.init.enable   = true;   % turn on/off
ctrl.init.X_end    = 0.01;   % apply while Xtot < 2%  (or use time window below)
ctrl.init.t_end    = 1e-2;   % also apply for the first 0.1 s (either condition triggers)
% ctrl.init.fp       = 0.05;   % allow per-step Δp up to 5% of p_eq  (0.02–0.1 reasonable)
ctrl.init.fp       = 1;   % allow per-step Δp up to 5% of p_eq  (0.02–0.1 reasonable)
ctrl.init.eps_q    = 1e-30;  % tiny to avoid divide-by-zero
ctrl.init.betaDiff = 0.5;    % diffusive limiter softness (0.2–1)

%% 4) 收集器（含 BET/孔体积/平均孔径）
snap_levels   = 0:100;
num_levels    = numel(snap_levels);
snap_filled   = false(1, num_levels);
snap = struct('level',[],'t',[],'Xtot',[],'phi',[], ...
              'p_pa',[],'X',[],'throat_radii',[],'T',[], ...
              'r_p',[], ...
              'BET_m2g',[],'PoreVol_m3g',[],'MeanPoreRadius_m',[]);
snap = repmat(snap,1,num_levels);

% 初始几何指标
V_void_node_ini = node_void_volume(P, Tn);
[~, ~] = node_solid_volume(V_CaCO3, V_CaCO3_0, solid);
S_internal_ini = geom_internal_area(P, Tn);              % m^2
m_solid_ini    = solid_mass_kg(V_CaCO3, V_CaCO3_0, solid); % kg
BET0   = S_internal_ini / m_solid_ini;                   % m2/kg -> m2/kg (单位一致)
PV0    = sum(V_void_node_ini) / m_solid_ini;            % m3/kg
RPM0   = 2*sum(V_void_node_ini) / max(S_internal_ini,eps); % m

% 0% 快照
snap(1).level = 0;  snap(1).t = 0.0;  snap(1).Xtot = Xtot;
snap(1).phi   = compute_porosity(P, Tn, R_particle);
snap(1).p_pa  = C_to_p(C, R_gas, Tfield);
snap(1).X     = X;
snap(1).throat_radii = Tn.r_t;
snap(1).T     = Tfield;
snap(1).r_p   = P.r_p;
snap(1).BET_m2g         = BET0 / 1e3;    % m2/g
snap(1).PoreVol_m3g     = PV0 / 1e3;     % m3/g
snap(1).MeanPoreRadius_m= RPM0;
snap_filled(1) = true;

% 时间序列
t_series    = 0.0;
Xtot_series = Xtot;
BET_series_m2g      = snap(1).BET_m2g;
PV_series_m3g       = snap(1).PoreVol_m3g;
Rmean_series_m      = snap(1).MeanPoreRadius_m;

%% 5) 主循环
dt   = ctrl.dt_init;  t = 0.0;  step = 0;

while t < ctrl.t_max
    % === 节点空隙/固体体积、局部孔隙率 ===
    V_void_node = node_void_volume(P, Tn);
    [V_solid_node, ~] = node_solid_volume(V_CaCO3, V_CaCO3_0, solid);
    eps_node = local_porosity(V_void_node, V_solid_node);


% ===== Initial-stage dt limiter (before solving the mass step) =====
if ctrl.init.enable && (Xtot < ctrl.init.X_end || t < ctrl.init.t_end)
    % Predict local source rate using current state (same as later)
    q_pred0 = source_molar_reversible(C, Tfield, X, S0, phys, R_gas);  % mol/s

    % Allowed concentration jump as a fraction of local equilibrium conc.
    p_eq_i  = Peq_CO2_Pa(Tfield, phys.Peq_A, phys.Peq_B);
    C_eq_i  = p_eq_i ./ (R_gas .* Tfield);                                   % mol/m^3
    dC_allow= ctrl.init.fp * C_eq_i;                                          % mol/m^3 per step

    % Source-driven dt cap (per node), using V_void as capacity
    dt_src_i = dC_allow .* V_void_node ./ (abs(q_pred0) + ctrl.init.eps_q);   % s
    dt_src   = max(ctrl.dt_min, min(dt_src_i(~isnan(dt_src_i) & isfinite(dt_src_i))));

    % Gentle diffusive cap: dt_diff_i ~ beta * V / Aii  (Aii later = row-sum of Ac)
    % Build a quick Ac using current eps_node (cheap, we need Aii only)
    Gm_tmp = throat_conductance_mass_structured(Tn, phys, R_gas, Tfield, eps_node, micro);
    Ac_tmp = B' * spdiags(Gm_tmp,0,M,M) * B;
    Aii    = full(diag(Ac_tmp));                             % s^-1 in "volume" units
    mask   = (Aii>0);
    if any(mask)
        dt_diff = ctrl.init.betaDiff * min( V_void_node(mask) ./ Aii(mask) );
        dt_diff = max(ctrl.dt_min, dt_diff);
    else
        dt_diff = ctrl.dt_max;
    end

    dt_cap = min([dt_src, dt_diff, ctrl.dt_max]);

    % If current dt exceeds the safe cap, shrink and restart this iteration
    if dt > dt_cap
        dt = max(ctrl.dt_min, max(dt_cap, ctrl.shrink*dt));
        % restart loop with smaller dt and rebuilt operators
        continue
    end
end




    % ===== 质量网络 =====
    Gm = throat_conductance_mass_structured(Tn, phys, R_gas, Tfield, eps_node, micro);
    Ac = B' * spdiags(Gm,0,M,M) * B;

    % 质量 Robin
    if massBC.use_robin
        [Hm_diag, bHm] = boundary_mass_robin(P, R_particle, I_b, massBC, phys);
    else
        Hm_diag = zeros(N,1);  bHm = zeros(N,1);
    end

    % 质量矩阵
    Mc = spdiags(V_void_node, 0, N, N);

    % 源项与限幅
    q_pred = source_molar_reversible(C, Tfield, X, S0, phys, R_gas);
    qn     = limit_reaction_rate(q_pred, dt, C, V_void_node, V_CaCO3, V_CaCO3_0, solid);

    % 稳健化：死节点与连通分量弱锚
    deg = effective_degree_per_node(B, Gm, N);
    dead = (V_void_node<=0) & (deg==0) & (Hm_diag==0);
    if any(dead)
        H_anchor = 1.0;
        Hm_diag(dead) = H_anchor;
        bHm(dead)     = H_anchor .* C(dead);
    end
    [Hm_diag, bHm] = ensure_component_anchors(Hm_diag, bHm, Tn, Gm, C);

    % 组装并预缩放
    LHS_c = (Mc/dt) + Ac + spdiags(Hm_diag,0,N,N);
    rhs_c = (Mc/dt)*C + qn + bHm;
    dscale = sqrt(max(abs(diag(LHS_c)), 1));
    D = spdiags(1./dscale, 0, N, N);
    LHSs = D * LHS_c * D;  rhss = D * rhs_c;

    % 解
    C_new = D * ( LHSs \ rhss );
    C_new = max(C_new, 0);

    % ===== 能量网络 =====
    Cth_vec = node_heat_capacity(C, Tfield, V_void_node, V_CaCO3, V_CaCO3_0, phys, solid);
    MT = spdiags(Cth_vec,0,N,N);
    if therm.enable_conduction
        Gth = throat_conductance_thermal(Tn, therm);
        AT  = B' * spdiags(Gth,0,M,M) * B;
    else
        AT  = sparse(N,N);
    end
    [Hdiag, bH] = boundary_heat_robin(P, R_particle, I_b, Tfield(I_b), phys.T_inf, therm, sigmaSB);
    Qreac = -phys.deltaH * qn;
    LHS_T = (MT/dt) + AT + spdiags(Hdiag,0,N,N);
    rhs_T = (MT/dt)*Tfield + Qreac + bH;
    T_new = LHS_T \ rhs_T;

    % ===== 几何演化 =====
    dn = qn * dt;                          % mol（±）
    dV_void = solid.dOmega * dn;           % m^3（±）
    [P_new, Tn_new] = apply_void_volume_and_update_geometry_signed(P, Tn, dV_void, 0.7);

    % 试算与自适应
    V_CaCO3_try = V_CaCO3 - dn*solid.Omega_CaCO3;
    V_CaCO3_try = min(max(V_CaCO3_try, 0), V_CaCO3_0);
    X_try = zeros(N,1);
    maskV0 = (V_CaCO3_0 > 0);
    X_try(maskV0) = 1 - V_CaCO3_try(maskV0)./V_CaCO3_0(maskV0);
    X_try = min(max(X_try,0),1);

    dX_max = max(abs(X_try - X));
    if dX_max > ctrl.dX_tol && dt > ctrl.dt_min
        dt = max(ctrl.dt_min, dt*ctrl.shrink);
        continue
    end

    % ===== 接受该步 =====
    C = C_new;      Tfield = T_new;
    P = P_new;      Tn     = Tn_new;
    V_CaCO3 = V_CaCO3_try;
    X = X_try;

    t    = t + dt;
    step = step + 1;

    Xtot = overall_conversion(V_CaCO3, V_CaCO3_0);

    % —— 计算并记录 BET/孔体积/平均孔径 —— 
    V_void_node = node_void_volume(P, Tn);
    S_internal  = geom_internal_area(P, Tn);           % m^2
    m_solid     = solid_mass_kg(V_CaCO3, V_CaCO3_0, solid); % kg
    BET_m2g     = (S_internal / m_solid) / 1e3;        % m2/g
    PV_m3g      = (sum(V_void_node) / m_solid) / 1e3;  % m3/g
    RP_mean     = 2*sum(V_void_node) / max(S_internal,eps); % m

    % —— 时间序列 —— 
    t_series(end+1,1)       = t;            %#ok<AGROW>
    Xtot_series(end+1,1)    = Xtot;         %#ok<AGROW>
    BET_series_m2g(end+1,1) = BET_m2g;      %#ok<AGROW>
    PV_series_m3g(end+1,1)  = PV_m3g;       %#ok<AGROW>
    Rmean_series_m(end+1,1) = RP_mean;      %#ok<AGROW>

    % —— 快照（按 level%） —— 
    lvl = floor(Xtot*100);
    for Lk = 0:lvl
        idx = Lk + 1;
        if ~snap_filled(idx)
            phi_now = compute_porosity(P, Tn, R_particle);
            p_pa    = C_to_p(C, R_gas, Tfield);
            snap(idx).level = Lk;
            snap(idx).t     = t;
            snap(idx).Xtot  = Xtot;
            snap(idx).phi   = phi_now;
            snap(idx).p_pa  = p_pa;
            snap(idx).X     = X;
            snap(idx).throat_radii = Tn.r_t;
            snap(idx).T     = Tfield;
            snap(idx).r_p   = P.r_p;
            snap(idx).BET_m2g         = BET_m2g;
            snap(idx).PoreVol_m3g     = PV_m3g;
            snap(idx).MeanPoreRadius_m= RP_mean;
            snap_filled(idx) = true;
        end
    end

    % 终止
    if Xtot > ctrl.Xtot_stop
        fprintf('Converged: Xtot=%.4f at t=%.6g s, step=%d\n', Xtot, t, step);
        break
    end

    % 放大步长
    if dX_max < 0.2*ctrl.dX_tol
        dt = min(ctrl.dt_max, dt*ctrl.grow);
    end

    if mod(step,50)==0
        fprintf('t=%.3e, Xtot=%.3f, dt=%.2e | BET=%.2e m2/g, PV=%.2e m3/g, rmean=%.2en m\n', ...
            t, Xtot, dt, BET_m2g, PV_m3g, RP_mean);
    end
end

%% 6) 保存
results = struct();
results.snap_levels   = snap_levels;
results.snapshots     = snap;
results.t_series      = t_series;
results.Xtot_series   = Xtot_series;
results.metrics.BET_series_m2g   = BET_series_m2g;
results.metrics.PV_series_m3g    = PV_series_m3g;
results.metrics.Rmean_series_m   = Rmean_series_m;

results.params.phys   = phys;
results.params.micro  = micro;
results.params.therm  = therm;
results.params.massBC = massBC;
results.params.solid  = solid;
results.params.ctrl   = ctrl;

save('reaction_results.mat','results','-v7.3');
fprintf('Saved: reaction_results.mat\n');

end % ===== end main =====


%% ==================== 辅助函数 ====================

function deg = effective_degree_per_node(B, Gm, N)
mask = (Gm > 0);
if ~any(mask), deg = zeros(N,1); return; end
Bpos = B(mask, :);
deg = full(sum(abs(Bpos),1))';
end

function [Hm_diag, bHm] = ensure_component_anchors(Hm_diag, bHm, Tn, Gm, Cn)
N = numel(Hm_diag);
i = Tn.pore1(:); j = Tn.pore2(:);
mask = (Gm > 0);
if ~any(mask)
    k = 1;
    if Hm_diag(k)==0
        Hweak = 1.0;
        Hm_diag(k) = Hweak;  bHm(k) = Hweak * Cn(k);
    end
    return;
end
Aadj = sparse([i(mask); j(mask)], [j(mask); i(mask)], 1, N, N);
Gobj = graph(Aadj);
comp = conncomp(Gobj);
nComp = max(comp);
has_anchor = accumarray(comp', Hm_diag>0, [nComp,1], @any, false);
for cc = 1:nComp
    if ~has_anchor(cc)
        k = find(comp==cc, 1, 'first');
        if ~isempty(k) && Hm_diag(k)==0
            Hweak = 1.0;
            Hm_diag(k) = Hweak;  bHm(k) = Hweak * Cn(k);
        end
    end
end
end

function G = throat_conductance_mass_structured(Tnet, phys, R_gas, Tnodes, eps_node, micro)
% 质量传输导通：G = Deff * A / L
% 这里采用 Carman–Kozeny：Deff = Dmix * ε^3 / (1-ε)^2
% 注：按你的要求，不做 ε→0 或 ε→1 的安全处理。

r  = Tnet.r_t(:);
L  = Tnet.L_geom(:);
assert(all(L>0),'存在喉长<=0');

% 几何横截面积（含 shape_factor）
A_geom = pi * (r.^2) .* max(Tnet.shape_factor(:), 0.05);

% 节点温度插值到喉
Ti = Tnodes(Tnet.pore1);
Tj = Tnodes(Tnet.pore2);
Te = 0.5*(Ti + Tj);

% 分子扩散与 Knudsen 扩散并联（调和平均）
Dm = phys.Dm;   % 如需温度相关可改为 Dm(Te)
Dk = (2/3) * r .* sqrt( 8*R_gas*Te ./ (pi*phys.M_CO2) );  % m^2/s
Dmix = 1 ./ ( 1./Dm + 1./Dk );

% 节点孔隙率插值到喉（简单平均，也可用调和/最小）
eps_i = eps_node(Tnet.pore1);
eps_j = eps_node(Tnet.pore2);
eps_e = 0.5*(eps_i + eps_j);

% ---- Carman–Kozeny 有效扩散系数 ----
Deff = Dmix .* (eps_e.^3) ./ ((1 - eps_e).^2);

% 喉道导通率
G = Deff .* (A_geom ./ L);     % m^3/s
end

function [Hm_diag, bHm] = boundary_mass_robin(P, R_particle, I_b, massBC, phys)
N = numel(P.r_p);
Hm_diag = zeros(N,1);  bHm = zeros(N,1);
if isempty(I_b) || ~massBC.use_robin, return; end
Ns = numel(I_b);
A_sphere = 4*pi*R_particle^2;
A_ext_i  = (A_sphere / Ns) * ones(Ns,1);
switch lower(massBC.kg_mode)
    case 'correlation'
        Sh = 2 + 0.6 * sqrt(max(massBC.Re,0)) * (max(massBC.Sc,0))^(1/3);
        kg = Sh * massBC.D_inf / (2*R_particle);
    otherwise
        kg = massBC.kg_const;
end
Hm_i = kg .* A_ext_i;             % m^3/s
Hm_diag(I_b) = Hm_i;
bHm(I_b)     = Hm_i .* massBC.C_inf;
end

function q = source_molar_reversible(C, Tnodes, X, S0, phys, R_gas)
p      = C_to_p(C, R_gas, Tnodes);
p_eq   = Peq_CO2_Pa(Tnodes, phys.Peq_A, phys.Peq_B);
ratio  = p ./ p_eq;
ks_T   = ks_of_T(Tnodes, phys);
n_dot  = ks_T .* (1 - ratio);                      % mol/m^2/s（±）
S_eff  = S0 .* max(0,(1 - X)) .* exp(phys.gammaSX * X);
q      = n_dot .* S_eff;
end

function q_lim = limit_reaction_rate(q_pred, dt, C, V_void_node, V_CaCO3, V_CaCO3_0, solid)
n_CaCO3_avail = V_CaCO3 ./ solid.Omega_CaCO3;
n_CaO_avail   = (V_CaCO3_0 - V_CaCO3) ./ solid.Omega_CaCO3;
n_CO2_avail   = V_void_node .* C;
dn_raw = q_pred .* dt;
s = ones(size(q_pred));
mask_pos = (dn_raw > 0);
if any(mask_pos)
    cap_fwd = n_CaCO3_avail(mask_pos);
    s(mask_pos) = min( 1, cap_fwd ./ dn_raw(mask_pos) );
end
mask_neg = (dn_raw < 0);
if any(mask_neg)
    cap_rev = min( n_CaO_avail(mask_neg), n_CO2_avail(mask_neg) );
    s(mask_neg) = min( 1, cap_rev ./ abs(dn_raw(mask_neg)) );
end
q_lim = s .* q_pred;
end

function p = C_to_p(C, R_gas, Tnodes)
p = C .* (R_gas .* Tnodes);
end

function p_eq_pa = Peq_CO2_Pa(T, A, B)
assert(all(T>0),'温度需为正');
p_eq_pa = (10.^(A - B./T)) * 101325;
end

function ks = ks_of_T(T, phys)
switch lower(phys.ks_model)
    case 'arrhenius'
        ks = phys.k0 .* exp(-phys.Ea ./ (8.314462618 .* T));
    otherwise
        ks = phys.k_s .* ones(size(T));
end
end

function [V_solid_node, V_CaO] = node_solid_volume(V_CaCO3, V_CaCO3_0, solid)
consumed_vol = (V_CaCO3_0 - V_CaCO3);  consumed_vol(consumed_vol < 0) = 0;
V_CaO = consumed_vol * (solid.Omega_CaO / solid.Omega_CaCO3);
V_solid_node = V_CaCO3 + V_CaO;
end

function eps = local_porosity(V_void_node, V_solid_node)
V_tot = V_void_node + V_solid_node;
mask = (V_tot > 0);
eps = zeros(size(V_tot));
eps(mask) = V_void_node(mask) ./ V_tot(mask);
end

function Cth_vec = node_heat_capacity(C, Tnodes, V_void_node, V_CaCO3, V_CaCO3_0, phys, solid)
consumed_vol = (V_CaCO3_0 - V_CaCO3);  consumed_vol(consumed_vol < 0) = 0;
V_CaO = consumed_vol * (solid.Omega_CaO / solid.Omega_CaCO3);
rhoCp_CaCO3 = solid.rho_CaCO3 * solid.cp_CaCO3;
rhoCp_CaO   = solid.rho_CaO   * solid.cp_CaO;
Cth_solid = rhoCp_CaCO3 * V_CaCO3 + rhoCp_CaO * V_CaO;
rho_g = C .* phys.M_CO2;
cp_g  = phys.cp_CO2_gas;
Cth_gas = rho_g .* cp_g .* V_void_node;
Cth_vec = Cth_solid + Cth_gas;
end

function Vv = node_void_volume(P, Tnet)
Vp  = (4/3)*pi*(P.r_p.^3);
Vth = pi*(Tnet.r_t.^2).*Tnet.L_geom;
halfVth_to_end = accumarray([Tnet.pore1; Tnet.pore2], [Vth; Vth]/2, [numel(P.r_p),1], @sum, 0);
Vv = Vp + halfVth_to_end;
end

function S = geom_internal_area(P, Tnet)
% 内部可接触表面积（孔面 + 喉侧面的一半分摊到两端）
S_pore    = 4*pi*(P.r_p.^2);
S_th_side = 2*pi*Tnet.r_t.*Tnet.L_geom;
halfS_to_end = accumarray([Tnet.pore1; Tnet.pore2], [S_th_side; S_th_side]/2, [numel(P.r_p),1], @sum, 0);
S = sum(S_pore + halfS_to_end);
end

function m = solid_mass_kg(V_CaCO3, V_CaCO3_0, solid)
% 当前固体质量（kg）= m_CaCO3 + m_CaO
consumed_vol = (V_CaCO3_0 - V_CaCO3);  consumed_vol(consumed_vol<0)=0;
V_CaO = consumed_vol * (solid.Omega_CaO / solid.Omega_CaCO3);
m = solid.rho_CaCO3 * sum(V_CaCO3) + solid.rho_CaO * sum(V_CaO);
end

function Gth = throat_conductance_thermal(Tnet, therm)
r  = Tnet.r_t(:);
L  = Tnet.L_geom(:);
assert(all(L>0),'存在喉长<=0');
A_solid = therm.A_solid_factor * (pi * (r.^2) .* max(Tnet.shape_factor(:), 0.05));
Gth = therm.k_solid_eff * (A_solid ./ L);
end

function [Hdiag, bH] = boundary_heat_robin(P, R_particle, I_b, Tsurf, T_inf, therm, sigmaSB)
Ns = numel(I_b);   assert(Ns>0,'未找到表面节点');
A_sphere = 4*pi*R_particle^2;
A_ext_i  = (A_sphere / Ns) * ones(Ns,1);
Tbar = 0.5*(Tsurf + T_inf);
h_rad = 4 * therm.emiss * sigmaSB .* (Tbar.^3);
h_eff = therm.h_conv + h_rad;
Hdiag = zeros(size(P.r_p));  Hdiag(I_b) = h_eff .* A_ext_i;
bH    = zeros(size(P.r_p));  bH(I_b)    = Hdiag(I_b) .* T_inf;
end

function Xtot = overall_conversion(V, V0)
S0 = sum(V0);  assert(S0>0,'初始固相总体积为0');
Xtot = 1 - sum(V)/S0;
Xtot = min(1, max(0, Xtot));
end

function phi = compute_porosity(P, Tnet, R_particle)
Vp  = (4/3)*pi*(P.r_p.^3);
Vth = pi*(Tnet.r_t.^2).*Tnet.L_geom;
V_sphere = (4/3)*pi*R_particle^3;
phi = (sum(Vp)+sum(Vth))/V_sphere;
end

function [P_out, T_out] = apply_void_volume_and_update_geometry_signed(P_in, T_in, dV_void_node, alpha_rt)
N = numel(P_in.r_p);
M = numel(T_in.pore1);
S_pore    = 4*pi*(P_in.r_p.^2);
S_th_side = 2*pi*T_in.r_t.*T_in.L_geom;
halfS_to_end  = accumarray([T_in.pore1; T_in.pore2], [S_th_side; S_th_side]/2, [N,1], @sum, 0);
S_total_node  = S_pore + halfS_to_end;

V_to_pore = zeros(N,1);
maskS = (S_total_node > 0);
V_to_pore(maskS)  = dV_void_node(maskS) .* (S_pore(maskS) ./ S_total_node(maskS));
V_to_pore(~maskS) = dV_void_node(~maskS);

coef_end = zeros(N,1);
coef_end(maskS) = dV_void_node(maskS) ./ S_total_node(maskS);
V_half_end = zeros(M,2);
V_half_end(:,1) = coef_end(T_in.pore1) .* (S_th_side/2);
V_half_end(:,2) = coef_end(T_in.pore2) .* (S_th_side/2);
V_to_throat = sum(V_half_end,2);

V_pore_old = (4/3)*pi*(P_in.r_p.^3);
V_pore_new = V_pore_old + V_to_pore;
V_pore_new(V_pore_new < 0) = 0;
rp_new = zeros(N,1);
mask_pos = (V_pore_new > 0);
rp_new(mask_pos) = ((3/(4*pi))*V_pore_new(mask_pos)).^(1/3);

A_old = pi*(T_in.r_t.^2);
dA = V_to_throat ./ T_in.L_geom;
A_new = A_old + dA;
A_new(A_new < 0) = 0;
rt_new_unclamped = zeros(M,1);
mask_Apos = (A_new > 0);
rt_new_unclamped(mask_Apos) = sqrt(A_new(mask_Apos)/pi);

rt_new = rt_new_unclamped; %#ok<NASGU> % 如需上限可加 alpha_rt

P_out = P_in;  P_out.r_p = rp_new;
T_out = T_in;  T_out.r_t = rt_new_unclamped;
end
