function profiles = viz_radial_reaction(varargin)
% 可视化：不同整体转化率下，沿颗粒半径方向的局部反应程度 X(r)
% 输入（Name-Value）：
%   'PNMFile'     : 默认 'PNM.mat'
%   'ResultsFile' : 默认 'reaction_results.mat'
%   'Levels'      : 需要绘制的整体转化率百分数数组（默认自动从结果里挑典型: [0 10 20 40 60 80 90 95 99 100]）
%   'NBins'       : 径向分桶数，默认 40
%   'OutDir'      : 输出图目录，默认 'fig_radial'
%   'MakeHeatmap' : 是否输出热图，默认 true
%   'MakeLines'   : 是否输出多曲线图，默认 true
%
% 返回：
%   profiles: 结构数组（每个 level 一项）
%       .level, .t, .Xtot
%       .r_mid          : 径向桶中心 (m)
%       .r_over_R       : 归一化半径 (r_mid/R)
%       .meanX, .q25, .q75
%       .front_r_over_R : 平均X=0.5的前沿位置 (r/R)，若不可求则 NaN
%
% 依赖：需要当前目录有 PNM.mat 与 reaction_results.mat（由你的第三版模拟脚本生成）

% 调用：
% profiles = viz_radial_reaction('Levels',[0 10 20 40 60 80 90 95 99 100], 'NBins', 40);

% ------------ 解析参数 -------------
p = inputParser;
addParameter(p,'PNMFile','PNM.mat');
addParameter(p,'ResultsFile','reaction_results.mat');
% addParameter(p,'ResultsFile','reaction_results_V8_codex_1200K.mat');
addParameter(p,'Levels',[]);
addParameter(p,'NBins',40);
addParameter(p,'OutDir','fig_radial');
addParameter(p,'MakeHeatmap',true,@islogical);
addParameter(p,'MakeLines',true,@islogical);
parse(p,varargin{:});
opt = p.Results;

assert(isfile(opt.PNMFile),'未找到 PNM 文件：%s', opt.PNMFile);
assert(isfile(opt.ResultsFile),'未找到结果文件：%s', opt.ResultsFile);

% ------------ 读取数据 -------------
S_pnm = load(opt.PNMFile);    PNM = S_pnm.PNM;
S_res = load(opt.ResultsFile); results = S_res.results;

xyz = PNM.P.coords;
r   = sqrt(sum(xyz.^2,2));     % 节点半径
R   = PNM.meta.R;
N   = size(xyz,1);

% 哪些 level 可用
snap = results.snapshots;
has_level = arrayfun(@(s) ~isempty(s.level), snap);
avail_levels = [snap(has_level).level];

% 默认 Levels
if isempty(opt.Levels)
    default_levels = [0 10 20 40 60 80 90 95 99 100];
    % default_levels = [0 1 2 3 4 5 6 7 8 9 10];
    opt.Levels = default_levels(ismember(default_levels, avail_levels));
end

% 过滤只留下确实存在的 level
levels = opt.Levels(:)';
levels = levels(ismember(levels, avail_levels));
assert(~isempty(levels),'所选 Levels 在结果中都不存在。可用的 Levels: %s', mat2str(sort(avail_levels)));

% ------------ 径向分桶 -------------
NB = opt.NBins;
edges = linspace(0, R, NB+1);
r_mid = 0.5*(edges(1:end-1)+edges(2:end));
r_over_R = r_mid / R;

% ------------ 为热图准备矩阵（可选） -------------
if opt.MakeHeatmap
    heat_levels = sort(unique(avail_levels));  % 尽量全覆盖已有快照
    H = nan(numel(heat_levels), NB);           % 行: level, 列: 径向桶
end

% ------------ 主循环：逐 level 统计并绘图数据 -------------
profiles = struct('level',[],'t',[],'Xtot',[], ...
                  'r_mid',[],'r_over_R',[],'meanX',[],'q25',[],'q75',[], ...
                  'front_r_over_R',[]);
profiles = repmat(profiles, 1, numel(levels));

for iL = 1:numel(levels)
    L = levels(iL);
    idxSnap = find(arrayfun(@(s) ~isempty(s.level) && s.level==L, snap), 1, 'first');
    S = snap(idxSnap);

    % 安全检查
    assert(numel(S.X)==N, '快照 level=%d 的 X 长度(%d)与节点数(%d)不符', L, numel(S.X), N);

    Xi = S.X(:);          % 节点转化率
    % 统计（逐桶）
    [mX, q25, q75] = radial_stats(r, Xi, edges);

    % 反应前沿（平均 X=0.5 的 r/R）
    front = front_location(r_mid, mX, 0.5);
    profiles(iL) = struct('level',L, 't',S.t, 'Xtot',S.Xtot, ...
                          'r_mid',r_mid(:), 'r_over_R',r_over_R(:), ...
                          'meanX',mX(:), 'q25',q25(:), 'q75',q75(:), ...
                          'front_r_over_R',front);
end

% 填热图矩阵（可选）
if opt.MakeHeatmap
    for il = 1:numel(heat_levels)
        L = heat_levels(il);
        idxSnap = find(arrayfun(@(s) ~isempty(s.level) && s.level==L, snap), 1, 'first');
        S = snap(idxSnap);
        [mX, ~, ~] = radial_stats(r, S.X(:), edges);
        H(il,:) = mX;
    end
end

% ------------ 作图（存图） -------------
if ~exist(opt.OutDir,'dir'), mkdir(opt.OutDir); end

% (1) 多曲线图：不同整体转化率的 X(r/R) 平均曲线 + 四分位带
if opt.MakeLines
    f1 = figure('Color','w','Position',[80 80 820 460]);
    hold on; grid on; box on;
    cmap = lines(max(7, numel(levels)));
    for iL = 1:numel(levels)
        pr = profiles(iL);
        c  = cmap(1+mod(iL-1,size(cmap,1)),:);
        % 四分位带（半透明）
        fill([pr.r_over_R' fliplr(pr.r_over_R')], [pr.q25' fliplr(pr.q75')], c, ...
             'FaceAlpha',0.12,'EdgeColor','none');
        % 均值曲线
        plot(pr.r_over_R, pr.meanX, '-', 'LineWidth', 1.8, 'Color', c);
        if ~isnan(pr.front_r_over_R)
            plot(pr.front_r_over_R, 0.5, 'o', 'Color', c, 'MarkerFaceColor', c, 'MarkerSize', 4);
        end
    end
    xlabel('r / R'); ylabel('<X>(r)'); ylim([0 1]); xlim([0 1]);
    title('Radial conversion profiles at selected overall conversions');
    legtxt = arrayfun(@(pr) sprintf('%d%%  (t=%.3g s)', pr.level, pr.t), profiles, 'UniformOutput',false);
    legend(legtxt, 'Location','SouthEast');
    saveas(f1, fullfile(opt.OutDir, 'radial_lines.png'));
    % close(f1);
end

% (2) 热图：行=整体转化率，列=r/R，颜色=桶平均 X
if opt.MakeHeatmap
    f2 = figure('Color','w','Position',[100 100 860 420]);
    imagesc(r_over_R, heat_levels, H); axis xy;
    xlabel('r / R'); ylabel('overall conversion (%)');
    cb = colorbar; ylabel(cb, 'mean X');
    caxis([0 1]);
    title('Radial conversion heatmap');
    saveas(f2, fullfile(opt.OutDir, 'radial_heatmap.png'));
    % close(f2);
end

% (3) 反应前沿位置 vs 整体转化率
f3 = figure('Color','w','Position',[100 100 560 420]);
Lv = arrayfun(@(pr) pr.level, profiles);
Fv = arrayfun(@(pr) pr.front_r_over_R, profiles);
plot(Lv, Fv, 'o-','LineWidth',1.6); grid on; ylim([0 1]); xlim([min(Lv) max(Lv)]);
xlabel('overall conversion (%)');
ylabel('front location r_{X=0.5} / R');
title('Apparent reaction front vs overall conversion');
saveas(f3, fullfile(opt.OutDir, 'front_vs_level.png'));
% close(f3);

fprintf('完成：图像已保存到 %s/ ，返回 profiles 结构用于后续分析。\n', opt.OutDir);
end


% ====== 工具函数：径向统计 ======
function [mX, q25, q75] = radial_stats(r, X, edges)
NB = numel(edges) - 1;
mX  = nan(1, NB);
q25 = nan(1, NB);
q75 = nan(1, NB);
[~,~,bin] = histcounts(r, edges);
for k = 1:NB
    idx = (bin == k);
    if any(idx)
        xk = X(idx);
        mX(k)  = mean(xk);
        q = prctile(xk, [25 75]);
        q25(k) = q(1);
        q75(k) = q(2);
    end
end
end

% ====== 工具函数：从均值曲线估计“前沿位置” ======
function r0 = front_location(r_mid, meanX, xstar)
% 在 meanX(r) 上找 X=xstar 的位置，线性插值；若不存在则 NaN
r0 = NaN;
if all(isnan(meanX)), return; end
% 找到跨越点：meanX 由 <xstar 到 >=xstar
above = meanX >= xstar;
if any(above)
    k = find(above, 1, 'first');
    if k == 1
        r0 = r_mid(1);
    else
        x1 = meanX(k-1); x2 = meanX(k);
        r1 = r_mid(k-1); r2 = r_mid(k);
        if x2 ~= x1
            r0 = r1 + (xstar - x1) * (r2 - r1) / (x2 - x1);
        else
            r0 = r1; % 平台
        end
    end
end
% 归一化在外层调用中完成
end
