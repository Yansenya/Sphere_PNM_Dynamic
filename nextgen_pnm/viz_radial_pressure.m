function profiles = viz_radial_pressure(varargin)
% 可视化：不同整体转化率下，沿颗粒半径方向的 pCO2(r) 分布
% 输入（Name-Value）：
%   'PNMFile'     : 默认 'PNM.mat'
%   'ResultsFile' : 默认 'reaction_results.mat'
%   'Levels'      : 要画的整体转化率百分数数组（默认自动挑常用档）
%   'NBins'       : 径向分桶数，默认 40
%   'OutDir'      : 输出图目录，默认 'fig_radial_p'
%   'MakeHeatmap' : 是否输出热图，默认 true
%   'MakeLines'   : 是否输出多曲线图，默认 true
%   'YScale'      : 'linear' 或 'log'（仅影响坐标轴刻度，不做取对数变换），默认 'linear'
%
% 返回：
%   profiles: 结构数组（每个 level 一项）
%       .level, .t, .Xtot
%       .r_mid, .r_over_R
%       .meanP, .q25, .q75      （单位：Pa）

% 调用：
% profiles_p = viz_radial_pressure('Levels',[0 10 20 40 60 80 90 95 99 100], 'NBins',40, 'YScale','linear');  % 或 'log'


% ---------- 解析参数 ----------
p = inputParser;
addParameter(p,'PNMFile','PNM.mat');
addParameter(p,'ResultsFile','reaction_results.mat');
% addParameter(p,'ResultsFile','reaction_results_Oct23_1.mat');
addParameter(p,'Levels',[]);
addParameter(p,'NBins',40);
addParameter(p,'OutDir','fig_radial_p');
addParameter(p,'MakeHeatmap',true,@islogical);
addParameter(p,'MakeLines',true,@islogical);
addParameter(p,'YScale','linear',@(s) any(strcmpi(s,{'linear','log'})));
parse(p,varargin{:});
opt = p.Results;

assert(isfile(opt.PNMFile),    '未找到 PNM 文件：%s', opt.PNMFile);
assert(isfile(opt.ResultsFile),'未找到结果文件：%s', opt.ResultsFile);

% ---------- 读取数据 ----------
S_pnm = load(opt.PNMFile);    PNM = S_pnm.PNM;
S_res = load(opt.ResultsFile); results = S_res.results;

xyz = PNM.P.coords;
r   = sqrt(sum(xyz.^2,2));    % 节点半径
R   = PNM.meta.R;
assert(R > 0, '颗粒半径 R 必须为正');   % 防潜在0分母

snap = results.snapshots;
has_level = arrayfun(@(s) ~isempty(s.level), snap);
avail_levels = [snap(has_level).level];

% 默认 Levels
levels = opt.Levels;
if isempty(levels)
    default_levels = [0 10 20 40 60 80 90 95 99 100];
    % default_levels = [0 1 2 3 4 5 6 7 8 9 10];
    % default_levels = [0 1 2 3 4 5 6 7 8 9 10]+10;
    levels = default_levels(ismember(default_levels, avail_levels));
end
levels = levels(ismember(levels, avail_levels));
assert(~isempty(levels), '所选 Levels 均不存在。可用 Levels: %s', mat2str(sort(avail_levels)));

N = size(xyz,1);

% ---------- 径向分桶 ----------
NB = opt.NBins;
edges = linspace(0, R, NB+1);
r_mid = 0.5*(edges(1:end-1)+edges(2:end));
r_over_R = r_mid / R;

% ---------- 热图矩阵（可选） ----------
if opt.MakeHeatmap
    heat_levels = sort(unique(avail_levels));
    H = nan(numel(heat_levels), NB);   % 行: level，列: 径向桶
end

% ---------- 逐 level 统计 ----------
profiles = struct('level',[],'t',[],'Xtot',[], ...
                  'r_mid',[],'r_over_R',[],'meanP',[],'q25',[],'q75',[]);
profiles = repmat(profiles, 1, numel(levels));

for iL = 1:numel(levels)
    L = levels(iL);
    idxSnap = find(arrayfun(@(s) ~isempty(s.level) && s.level==L, snap), 1, 'first');
    S = snap(idxSnap);

    assert(isfield(S,'p_pa') && ~isempty(S.p_pa), ...
        '快照 level=%d 中没有 p_pa（Pa）。请确认第三版模拟已保存 p_pa。', L);
    Pi = S.p_pa(:);                  % 节点压强（Pa）
    assert(numel(Pi)==N, '快照 level=%d 的 p_pa 长度(%d)与节点数(%d)不符', L, numel(Pi), N);

    [mP, q25, q75] = radial_stats(r, Pi, edges);

    profiles(iL) = struct('level',L,'t',S.t,'Xtot',S.Xtot, ...
                          'r_mid',r_mid(:),'r_over_R',r_over_R(:), ...
                          'meanP',mP(:),'q25',q25(:),'q75',q75(:));
end

% 填热图
if opt.MakeHeatmap
    for il = 1:numel(heat_levels)
        L = heat_levels(il);
        idxSnap = find(arrayfun(@(s) ~isempty(s.level) && s.level==L, snap), 1, 'first');
        S = snap(idxSnap);
        assert(isfield(S,'p_pa') && ~isempty(S.p_pa), ...
            '快照 level=%d 中没有 p_pa（Pa）。', L);
        [mP, ~, ~] = radial_stats(r, S.p_pa(:), edges);
        H(il,:) = mP;
    end
end

% ---------- 作图 ----------
if ~exist(opt.OutDir,'dir'), mkdir(opt.OutDir); end

% (1) 多曲线：p(r/R)（均值+四分位带）
if opt.MakeLines
    f1 = figure('Color','w','Position',[90 90 840 480]);
    hold on; grid on; box on;
    cmap = lines(max(7, numel(levels)));
    for iL = 1:numel(levels)
        pr = profiles(iL);
        c  = cmap(1+mod(iL-1,size(cmap,1)),:);
        % 四分位带
        fill([pr.r_over_R' fliplr(pr.r_over_R')], [pr.q25' fliplr(pr.q75')], c, ...
             'FaceAlpha',0.12,'EdgeColor','none');
        % 均值曲线
        plot(pr.r_over_R, pr.meanP, '-', 'LineWidth', 1.8, 'Color', c);
    end
    set(gca,'YScale',lower(opt.YScale));
    xlabel('r / R'); ylabel('p_{CO2} (Pa)'); xlim([0 1]);
    title('Radial pressure profiles at selected overall conversions');
    legtxt = arrayfun(@(pr) sprintf('%d%%  (t=%.3g s)', pr.level, pr.t), profiles, 'UniformOutput',false);
    legend(legtxt, 'Location','bestoutside');
    saveas(f1, fullfile(opt.OutDir, 'pressure_radial_lines.png'));
    % close(f1);
end

% (2) 热图：行=整体转化率，列=r/R，颜色=平均 p（Pa）
if opt.MakeHeatmap
    f2 = figure('Color','w','Position',[100 100 880 440]);
    imagesc(r_over_R, heat_levels, H); axis xy;
    xlabel('r / R'); ylabel('overall conversion (%)');
    cb = colorbar; ylabel(cb, 'mean p_{CO2} (Pa)');
    % 颜色范围（自动或手动都可，这里自动即可）
    title('Radial pressure heatmap');
    saveas(f2, fullfile(opt.OutDir, 'pressure_radial_heatmap.png'));
    % close(f2);
end

fprintf('完成：已输出 p(r/R) 的曲线与热图到 %s/\n', opt.OutDir);
end


% ===== 径向统计工具（均值与四分位） =====
function [mV, q25, q75] = radial_stats(r, V, edges)
NB = numel(edges) - 1;
mV  = nan(1, NB);
q25 = nan(1, NB);
q75 = nan(1, NB);
[~,~,bin] = histcounts(r, edges);
for k = 1:NB
    idx = (bin == k);
    if any(idx)
        vk = V(idx);
        mV(k)  = mean(vk);
        qq     = prctile(vk, [25 75]);
        q25(k) = qq(1);
        q75(k) = qq(2);
    end
end
end
