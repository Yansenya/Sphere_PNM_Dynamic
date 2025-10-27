function profile = viz_radial_pore_radius(varargin)
% 可视化：沿颗粒半径方向的孔半径 r_p(r)
% 输入（Name-Value）：
%   'PNMFile'     : 默认 'PNM.mat'
%   'ResultsFile' : 默认 'reaction_results.mat'（可选，用于选择某个level的r_p，如果结果里保存了）
%   'Level'       : 选用的整体转化率level（百分数，默认 [] 表示用初始PNM）
%   'NBins'       : 径向分桶数，默认 40
%   'OutDir'      : 输出图目录，默认 'fig_radial'
%   'MakeScatter' : 是否画散点（r/R vs r_p），默认 true
%   'MakeLines'   : 是否画统计曲线（均值+四分位带），默认 true
%   'RPUnit'      : 显示单位 'm' 或 'um' 或 'nm'（默认 'um'）
%
% 返回：
%   profile: 结构体
%       .r_mid, .r_over_R
%       .mean_rp, .q25, .q75, .count
%       .unit, .level_used, .source ('initial' 或 'snapshot')

p = inputParser;
addParameter(p,'PNMFile','PNM.mat');
addParameter(p,'ResultsFile','reaction_results.mat');
addParameter(p,'Level',[]);
addParameter(p,'NBins',40);
addParameter(p,'OutDir','fig_radial');
addParameter(p,'MakeScatter',true,@islogical);
addParameter(p,'MakeLines',true,@islogical);
addParameter(p,'RPUnit','um',@(s)ischar(s) || isstring(s));
parse(p,varargin{:});
opt = p.Results;

assert(isfile(opt.PNMFile),'未找到 PNM 文件：%s', opt.PNMFile);
S_pnm = load(opt.PNMFile);  PNM = S_pnm.PNM;

xyz = PNM.P.coords;
r   = sqrt(sum(xyz.^2,2));     % 节点半径
R   = PNM.meta.R;




% rp0 = PNM.P.r_p(:);            % 初始孔半径（m）
results = [];
if isfile(opt.ResultsFile)
    S_res = load(opt.ResultsFile);
    if isfield(S_res,'results'); results = S_res.results; end
end
S = results.snapshots(50);
rp0 = S.r_p(:);       




N   = numel(rp0);

% 尝试从结果中取指定 level 的 r_p（只有当你在快照里保存了 r_p 才会命中）
rp   = rp0;
src  = 'initial';
Luse = [];

if ~isempty(opt.Level) && isfile(opt.ResultsFile)
    S_res = load(opt.ResultsFile);
    if isfield(S_res,'results')
        snap = S_res.results.snapshots;
        idx = find(arrayfun(@(s) ~isempty(s.level) && s.level==opt.Level && isfield(s,'r_p') && ~isempty(s.r_p), snap),1);
        if ~isempty(idx)
            rp_snap = snap(idx).r_p(:);
            if numel(rp_snap)==N
                rp   = rp_snap;
                src  = 'snapshot';
                Luse = opt.Level;
            else
                warning('快照 level=%d 的 r_p 长度(%d)与节点数(%d)不符，改用初始PNM。', opt.Level, numel(rp_snap), N);
            end
        else
            warning('结果文件中未找到 level=%d 的 r_p（或未保存 r_p），改用初始PNM。', opt.Level);
        end
    end
end

% 径向分桶
NB    = opt.NBins;
edges = linspace(0, R, NB+1);
r_mid = 0.5*(edges(1:end-1)+edges(2:end));
r_over_R = r_mid / R;

[mean_rp, q25, q75, count] = radial_stats_rp(r, rp, edges);

% 单位换算
switch lower(char(opt.RPUnit))
    case 'm'
        fac = 1; ustr = 'm';
    case 'nm'
        fac = 1e9; ustr = 'nm';
    otherwise % 'um'
        fac = 1e6; ustr = 'µm';
end
mean_rp_p = mean_rp*fac;  q25_p = q25*fac;  q75_p = q75*fac;
rp_plot   = rp*fac;

% 画图
if ~exist(opt.OutDir,'dir'), mkdir(opt.OutDir); end

ttl_suffix = '';
if strcmp(src,'snapshot') && ~isempty(Luse)
    ttl_suffix = sprintf(' (level %d%%)', Luse);
end

% (1) 散点：所有孔
if opt.MakeScatter
    f1 = figure('Color','w','Position',[80 80 820 460]);
    scatter(r/R, rp_plot, 10, 'filled', 'MarkerFaceAlpha', 0.15, 'MarkerEdgeAlpha', 0.15);
    grid on; box on;
    xlabel('r / R'); ylabel(sprintf('Pore radius r_p (%s)', ustr));
    title(['Radial distribution of pore radii', ttl_suffix]);
    saveas(f1, fullfile(opt.OutDir, sprintf('rp_radial_scatter_%s.png', src)));
end

% (2) 统计曲线：均值 + 四分位
if opt.MakeLines
    f2 = figure('Color','w','Position',[100 100 820 460]);
    hold on; grid on; box on;
    % 四分位带
    fill([r_over_R fliplr(r_over_R)], [q25_p fliplr(q75_p)], [0.3 0.6 0.9], ...
        'FaceAlpha',0.18,'EdgeColor','none');
    % 均值
    plot(r_over_R, mean_rp_p, 'b-', 'LineWidth', 2.0);
    xlabel('r / R'); ylabel(sprintf('<r_p>(%s)', ustr));
    title(['Binned pore-radius profile vs radius', ttl_suffix]);
    ylim([0, max(q75_p(~isnan(q75_p)))*1.1]);
    xlim([0 1]);
    saveas(f2, fullfile(opt.OutDir, sprintf('rp_radial_profile_%s.png', src)));
end

% 输出
profile = struct('r_mid', r_mid(:), 'r_over_R', r_over_R(:), ...
                 'mean_rp', mean_rp(:), 'q25', q25(:), 'q75', q75(:), ...
                 'count', count(:), 'unit', ustr, ...
                 'level_used', Luse, 'source', src);

fprintf('完成：已输出到 %s/ （source=%s）。\n', opt.OutDir, src);

end

% ====== 逐径向桶统计孔半径 ======
function [m, q25, q75, cnt] = radial_stats_rp(r, rp, edges)
NB = numel(edges) - 1;
m   = nan(1, NB);
q25 = nan(1, NB);
q75 = nan(1, NB);
cnt = zeros(1, NB);
[~,~,bin] = histcounts(r, edges);
for k = 1:NB
    idx = (bin == k);
    cnt(k) = nnz(idx);
    if cnt(k) > 0
        rk = rp(idx);
        m(k)   = mean(rk);
        q = prctile(rk, [25 75]);
        q25(k) = q(1);
        q75(k) = q(2);
    end
end
end
