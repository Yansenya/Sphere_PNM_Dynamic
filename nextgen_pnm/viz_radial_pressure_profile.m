function profile = viz_radial_pressure_profile(varargin)
% 可视化：沿颗粒半径方向的 CO2 压强 p(r)
% 输入（Name-Value）：
%   'PNMFile'     : 默认 'PNM.mat'
%   'ResultsFile' : 默认 'reaction_results.mat'
%   'Level'       : 选用的整体转化率 level（百分数，必须提供或留空=最近可用）
%   'NBins'       : 径向分桶数，默认 40
%   'OutDir'      : 输出目录，默认 'fig_radial'
%   'MakeScatter' : 是否画散点（r/R vs p），默认 true
%   'MakeLines'   : 是否画统计曲线（均值+四分位带），默认 true
%   'PUnit'       : 'Pa'|'kPa'|'bar'|'atm'（默认 'Pa'）
%   'YScale'      : 'linear' 或 'log'（默认 'linear'）
%
% 返回 profile 结构：
%   .r_mid, .r_over_R, .mean_p, .q25, .q75, .count, .unit, .level_used, .source

p = inputParser;
addParameter(p,'PNMFile','PNM.mat');
addParameter(p,'ResultsFile','reaction_results.mat');
addParameter(p,'Level',[]);
addParameter(p,'NBins',40);
addParameter(p,'OutDir','fig_radial');
addParameter(p,'MakeScatter',true,@islogical);
addParameter(p,'MakeLines',true,@islogical);
addParameter(p,'PUnit','Pa',@(s)ischar(s) || isstring(s));
addParameter(p,'YScale','linear',@(s)any(strcmpi(s,{'linear','log'})));
parse(p,varargin{:});
opt = p.Results;

assert(isfile(opt.PNMFile),'未找到 PNM 文件：%s', opt.PNMFile);
S_pnm = load(opt.PNMFile);  PNM = S_pnm.PNM;

xyz = PNM.P.coords;
r   = sqrt(sum(xyz.^2,2));     % 节点半径
R   = PNM.meta.R;
N   = numel(PNM.P.r_p);

% 读取结果并定位 level
assert(isfile(opt.ResultsFile),'未找到结果文件：%s', opt.ResultsFile);
S_res = load(opt.ResultsFile);
assert(isfield(S_res,'results') && isfield(S_res.results,'snapshots'), '结果文件缺少 snapshots');
snap = S_res.results.snapshots;

opt.Level = 20;

% 选择 level：若未提供 Level，则选用“最后一个已填充的快照”
if isempty(opt.Level)
    idx_all = find(arrayfun(@(s) ~isempty(s.level) && ~isempty(s.p_pa), snap));
    assert(~isempty(idx_all),'results.snapshots 里没有可用的 p_pa 快照');
    idx = idx_all(end);
else
    idx = find(arrayfun(@(s) ~isempty(s.level) && s.level==opt.Level && isfield(s,'p_pa') && ~isempty(s.p_pa), snap),1);
    assert(~isempty(idx),'没有找到 level=%d 的 p_pa 快照', opt.Level);
end

p_pa = snap(idx).p_pa(:);



p_pa(p_pa>150000) = 150000;



assert(numel(p_pa)==N,'快照中的 p_pa 长度(%d)与节点数(%d)不匹配', numel(p_pa), N);
src  = 'snapshot';
Luse = snap(idx).level;

% 径向分桶
NB    = opt.NBins;
edges = linspace(0, R, NB+1);
r_mid = 0.5*(edges(1:end-1)+edges(2:end));
r_over_R = r_mid / R;

[mean_p, q25, q75, count] = radial_stats_scalar(r, p_pa, edges);

% 单位换算
switch lower(char(opt.PUnit))
    case 'pa'
        fac = 1;     ustr = 'Pa';
    case 'kpa'
        fac = 1e-3;  ustr = 'kPa';
    case 'bar'
        fac = 1e-5;  ustr = 'bar';
    case 'atm'
        fac = 1/101325; ustr = 'atm';
    otherwise
        warning('未知单位 %s，使用 Pa', opt.PUnit);
        fac = 1; ustr = 'Pa';
end
mean_p_p = mean_p*fac;  q25_p = q25*fac;  q75_p = q75*fac;  p_plot = p_pa*fac;

% 画图
if ~exist(opt.OutDir,'dir'), mkdir(opt.OutDir); end
ttl_suffix = sprintf(' (level %d%%)', Luse);

% (1) 散点：所有孔
if opt.MakeScatter
    f1 = figure('Color','w','Position',[80 80 820 460]);
    scatter(r/R, p_plot, 10, 'filled', 'MarkerFaceAlpha', 0.15, 'MarkerEdgeAlpha', 0.15);
    grid on; box on;
    set(gca,'YScale',lower(opt.YScale));
    xlabel('r / R'); ylabel(sprintf('CO_2 partial pressure (%s)', ustr));
    title(['Radial distribution of p(r)', ttl_suffix]);
    saveas(f1, fullfile(opt.OutDir, sprintf('p_radial_scatter_lvl%02d.png', Luse)));
end

% (2) 统计曲线：均值 + 四分位
if opt.MakeLines
    f2 = figure('Color','w','Position',[100 100 820 460]);
    hold on; grid on; box on;
    set(gca,'YScale',lower(opt.YScale));
    % 四分位带
    valid = ~isnan(q25_p) & ~isnan(q75_p);
    if any(valid)
        fill([r_over_R(valid) fliplr(r_over_R(valid))], ...
             [q25_p(valid)    fliplr(q75_p(valid))], [0.3 0.6 0.9], ...
             'FaceAlpha',0.18,'EdgeColor','none');
    end
    % 均值
    plot(r_over_R, mean_p_p, 'b-', 'LineWidth', 2.0);
    xlabel('r / R'); ylabel(sprintf('<p>(%s)', ustr));
    title(['Binned pressure profile vs radius', ttl_suffix]);
    xlim([0 1]);
    % 自动 y 限（线性/对数分开处理）
    if strcmpi(opt.YScale,'linear')
        yvals = q75_p(~isnan(q75_p)); 
        if ~isempty(yvals), ylim([0, max(yvals)*1.1]); end
    end
    saveas(f2, fullfile(opt.OutDir, sprintf('p_radial_profile_lvl%02d.png', Luse)));
end

% 输出
profile = struct('r_mid', r_mid(:), 'r_over_R', r_over_R(:), ...
                 'mean_p', mean_p(:), 'q25', q25(:), 'q75', q75(:), ...
                 'count', count(:), 'unit', ustr, ...
                 'level_used', Luse, 'source', src);

fprintf('完成：已输出到 %s/ （level=%d%%，unit=%s）。\n', opt.OutDir, Luse, ustr);

end

% ====== 逐径向桶统计标量 ======
function [m, q25, q75, cnt] = radial_stats_scalar(r, val, edges)
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
        vk = val(idx);
        m(k)   = mean(vk);
        q = prctile(vk, [25 75]);
        q25(k) = q(1);
        q75(k) = q(2);
    end
end
end
