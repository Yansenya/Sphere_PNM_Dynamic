% function [coords, is_surface, s, jitter_sigma] = placePoresInSphere(R, N_target, opts)
% % 生成抖动的 FCC 点阵并裁剪到球体
% if ~isfield(opts,'jitter_sigma') || isempty(opts.jitter_sigma)
%     % 先预估 rp 的期望（用于 spacing 经验设定）
%     % 若不可得，则给个保守 spacing
%     E_rp = exp( (0.5) ); % 仅占位，不重要
%     s = 2.5 * 5e-6;      % 默认 12.5 μm
%     jitter_sigma = 0.15 * s;
% else
%     jitter_sigma = opts.jitter_sigma;
%     s = max(1e-6, jitter_sigma/0.15);
% end
% 
% % 构造 FCC 基元
% a = s;
% basis = [0 0 0; 0.5 0.5 0; 0.5 0 0.5; 0 0.5 0.5] * a;
% 
% % 构造包围球的立方体范围
% M = ceil((2*R)/a) + 3;
% grid = (-(M):(M))'*a;
% [X,Y,Z] = ndgrid(grid, grid, grid);
% C = [X(:), Y(:), Z(:)];
% 
% % 平移加上基元
% coords = [];
% for i=1:size(basis,1)
%     coords = [coords; C + basis(i,:)]; %#ok<AGROW>
% end
% 
% % 抖动
% coords = coords + randn(size(coords))*jitter_sigma;
% 
% % 裁剪到球
% r = sqrt(sum(coords.^2,2));
% mask = r <= R;
% coords = coords(mask,:);
% 
% % 若数量与目标差距较大，随机抽样/补点（简化处理）
% N = size(coords,1);
% if N > N_target
%     idx = randperm(N, N_target);
%     coords = coords(idx,:);
% elseif N < N_target
%     % 简易补点：在球内随机补充
%     need = N_target - N;
%     addp = randomPointsInSphere(R, need);
%     coords = [coords; addp];
% end
% 
% % 重新计算 is_surface
% N = size(coords,1);
% r = sqrt(sum(coords.^2,2));
% % 判定：与外界距离小于一格 spacing 视为表面
% is_surface = (R - r) <= (s*0.75);
% end
% 
% function pts = randomPointsInSphere(R, n)
% u = rand(n,1); v = rand(n,1); w = rand(n,1);
% theta = 2*pi*u;
% phi   = acos(2*v - 1);
% rad   = R*(w.^(1/3));
% pts = [rad.*sin(phi).*cos(theta), rad.*sin(phi).*sin(theta), rad.*cos(phi)];
% end

function [coords, is_surface, s, jitter_sigma] = placePoresInSphere(R, N_target, opts)
% 生成抖动的 FCC 点阵并裁剪到球体
% —— 关键改动：表面孔判定使用自适应壳层厚度 delta，
%    对任意 R 都稳健：delta = clip(0.8*s/sqrt(2), 0.005R, 0.08R)

if ~isfield(opts,'jitter_sigma') || isempty(opts.jitter_sigma)
    % 经验 spacing（与原版一致，保持行为不变）
    s = 2.5 * 5e-6;      % 默认 12.5 μm
    jitter_sigma = 0.15 * s;
else
    jitter_sigma = opts.jitter_sigma;
    s = max(1e-6, jitter_sigma/0.15);
end

% ---- FCC 基元 ----
a = s;
basis = [0 0 0; 0.5 0.5 0; 0.5 0 0.5; 0 0.5 0.5] * a;

% ---- 包围球的立方体范围 ----
M = ceil((2*R)/a) + 3;
grid = (-(M):(M))'*a;
[X,Y,Z] = ndgrid(grid, grid, grid);
C = [X(:), Y(:), Z(:)];

% ---- 平移加基元 ----
coords = [];
for i=1:size(basis,1)
    coords = [coords; C + basis(i,:)]; %#ok<AGROW>
end

% ---- 抖动 ----
coords = coords + randn(size(coords))*jitter_sigma;

% ---- 裁剪到球 ----
r = sqrt(sum(coords.^2,2));
mask = r <= R;
coords = coords(mask,:);

% ---- 数量与目标匹配（简化抽样/补点）----
N = size(coords,1);
if N > N_target
    idx = randperm(N, N_target);
    coords = coords(idx,:);
elseif N < N_target
    need = N_target - N;
    addp = randomPointsInSphere(R, need);
    coords = [coords; addp];
end

% ---- 自适应表面判定（对任意 R 稳健）----
N = size(coords,1);
r = sqrt(sum(coords.^2,2));

% 最近邻距离的近似：FCC d_nn ≈ a/√2 = s/√2
d_nn_est = s / sqrt(2);

% 壳层厚度：与局部点间距和 R 同时挂钩，避免 R 极端时阈值失衡
delta = 0.8 * d_nn_est;        % 以最近邻的 0.8 倍为基准
delta = max(delta, 0.005*R);   % 不低于 0.5% R
delta = min(delta, 0.08*R);    % 不高于 8%   R

% 表面孔：到外表面的距离 <= delta
is_surface = (R - r) <= delta;
end

function pts = randomPointsInSphere(R, n)
% 均匀体分布采样（半径按 r^(1/3)）
u = rand(n,1); v = rand(n,1); w = rand(n,1);
theta = 2*pi*u;
phi   = acos(2*v - 1);
rad   = R*(w.^(1/3));
pts = [rad.*sin(phi).*cos(theta), rad.*sin(phi).*sin(theta), rad.*cos(phi)];
end

