function profiles = viz_radial_pore_radius_evolution(varargin)
% 可视化：不同整体转化率下，沿颗粒半径方向的孔半径 r_p(r) 演化
% 输入（Name-Value）：
%   'PNMFile'     : 'PNM.mat'（必需）
%   'ResultsFile' : 'reaction_results.mat'（建议；从中读取各 level 的 r_p）
%   'Levels'      : 需要绘制的整体转化率百分数数组（默认自动挑选常用: [0 10 20 40 60 80 90 95 99 100] 里存在的）
%   'NBins'       : 径向分桶数（默认 40）
%   'OutDir'      : 输出目录（默认 'fig_radial_rp'）
%   'MakeHeatmap' : 是否输出热图（默认 true）
%   'MakeLines'   : 是否输出多曲线图（默认 true）
%   'RPUnit'      : 半径单位：'m'|'um'|'nm'（默认 'nm'）
%
% 返回：
%   profiles: 结构数组（每个 level 一项）
%       .level, .t, .Xtot
%       .r_mid (m), .r_over_R
%       .mean_rp, .q25, .q75（单位同 RPUnit）
%       .has_rp_snapshot（该 level 是否来自快照 r_p）
%
% 用法示例：
%   profiles = viz_radial_pore_radius_evolution('Levels',[0 20 40 60 80 100],'NBins',60,'RPUnit','nm');

% ------------ 解析参数 -------------
p = inputParser;
addParameter(p,'PNMFile','PNM.mat');
addParameter(p,'ResultsFile','reaction_results.mat');
addParameter(p,'Levels',[]);
addParameter(p,'NBins',40);
addParameter(p,'OutDir','fig_radial_rp');
addParameter(p,'MakeHeatmap',true,@islogical);
addParameter(p,'MakeLines',true,@islogical);
addParameter(p,'RPUnit','nm',@(s)ischar(s) || isstring(s));
parse(p,varargin{:});
opt = p.Results;

assert(isfile(opt.PNMFile),'未找到 PNM 文件：%s', opt.PNMFile);
S_pnm = load(opt.PNMFile);    PNM = S_pnm.PNM;

xyz = PNM.P.coords;
r   = sqrt(sum(xyz.^2,2));     % 节点半径
R   = PNM.meta.R;
N   = size(xyz,1);
rp_init = PNM.P.r_p(:);        % 初始孔半径（m）

% 读取结果（如有）
results = [];
if isfile(opt.ResultsFile)
    S_res = load(opt.ResultsFile);
    if isfield(S_res,'results'); results = S_res.results; end
end

% 可用的 levels
levels_avail = [];
if ~isempty(results)
    snap = results.snapshots;
    has_level = arrayfun(@(s) ~isempty(s.level), snap);
    levels_avail = [snap(has_level).level];
end

% 默认 Levels
if isempty(opt.Levels)
    default_levels = [0 10 20 40 60 80 90 95 99 100];
    % default_levels = [0 1 2 3 4 5 6 7 8 9 10];
    if ~isempty(levels_avail)
        opt.Levels = default_levels(ismember(default_levels, levels_avail));
        if isempty(opt.Levels)  % 兜底：全取
            opt.Levels = unique(levels_avail);
        end
    else
        opt.Levels = 0;  % 没有结果文件就画初始
    end
end

levels = unique(opt.Levels(:)') ;

% 径向分桶
NB    = opt.NBins;
edges = linspace(0, R, NB+1);
r_mid = 0.5*(edges(1:end-1)+edges(2:end));
r_over_R = r_mid / R;

% 单位换算因子
switch lower(char(opt.RPUnit))
    case 'm';  fac = 1;    ustr = 'm';
    case 'um'; fac = 1e6;  ustr = 'µm';
    otherwise % 'nm'
        fac = 1e9; ustr = 'nm';
end

% 为热图准备矩阵
if opt.MakeHeatmap
    H = nan(numel(levels), NB);     % mean(r_p) heatmap
end

% 输出目录
if ~exist(opt.OutDir,'dir'), mkdir(opt.OutDir); end

% 结果容器
profiles = struct('level',[],'t',[],'Xtot',[], ...
                  'r_mid',[],'r_over_R',[], ...
                  'mean_rp',[],'q25',[],'q75',[], ...
                  'has_rp_snapshot',[]);
profiles = repmat(profiles, 1, numel(levels));

% 主循环
for iL = 1:numel(levels)
    L = levels(iL);
    rpL = rp_init; tL = NaN; XtotL = NaN; hasSnap = false;

    if ~isempty(results)
        idxSnap = find(arrayfun(@(s) ~isempty(s.level) && s.level==L, results.snapshots), 1, 'first');
        if ~isempty(idxSnap)
            S = results.snapshots(idxSnap);
            tL = S.t;  XtotL = S.Xtot;
            if isfield(S,'r_p') && ~isempty(S.r_p)
                if numel(S.r_p)==N
                    rpL = S.r_p(:);
                    hasSnap = true;
                else
                    warning('level=%d 的 r_p 长度(%d)与节点数(%d)不符→回退到初始 r_p。', L, numel(S.r_p), N);
                end
            else
                % 没有保存 r_p：回退
                warning('结果文件中 level=%d 未保存 r_p → 回退到初始 r_p。', L);
            end
        else
            warning('结果文件中不存在 level=%d 的快照 → 使用初始 r_p。', L);
        end
    end

    % 径向统计
    [m, q25, q75, ~] = radial_stats_rp(r, rpL, edges);

    % 存到输出结构（转换到目标单位）
    profiles(iL).level = L;
    profiles(iL).t     = tL;
    profiles(iL).Xtot  = XtotL;
    profiles(iL).r_mid = r_mid(:);
    profiles(iL).r_over_R = r_over_R(:);
    profiles(iL).mean_rp  = (m*fac).';
    profiles(iL).q25      = (q25*fac).';
    profiles(iL).q75      = (q75*fac).';
    profiles(iL).has_rp_snapshot = hasSnap;

    % 热图行
    if opt.MakeHeatmap
        H(iL,:) = m * fac;
    end
end

% ——作图：多曲线（均值 + 四分位带）——
if opt.MakeLines
    f1 = figure('Color','w','Position',[90 90 880 500]);
    hold on; grid on; box on;
    cmap = lines(max(7, numel(levels)));
    for iL = 1:numel(levels)
        pr = profiles(iL);
        c  = cmap(1+mod(iL-1,size(cmap,1)),:);
        % 四分位带
        fill([pr.r_over_R' fliplr(pr.r_over_R')], [pr.q25' fliplr(pr.q75')], c, ...
             'FaceAlpha',0.12,'EdgeColor','none');
        % 均值曲线
        lw = 2.0;
        ls = '-';
        if ~pr.has_rp_snapshot
            ls = '--';   % 没有保存r_p的level用虚线提示
        end
        plot(pr.r_over_R, pr.mean_rp, ls, 'LineWidth', lw, 'Color', c);
    end
    xlabel('r / R');
    ylabel(sprintf('<r_p>(%s)', ustr));
    title('Radial profiles of pore radius at different overall conversions');
    xlim([0 1]);
    % 图例文本：带时间（若有）
    le = arrayfun(@(pr) ...
        sprintf('%d%% %s', pr.level, tern(isnan(pr.t),'','')), profiles, 'UniformOutput', false);
    % 更详细：含时间（s）
    for i=1:numel(levels)
        if ~isnan(profiles(i).t)
            le{i} = sprintf('%d%% (t=%.3g s)%s', profiles(i).level, profiles(i).t, tern(profiles(i).has_rp_snapshot,'','*'));
        else
            le{i} = sprintf('%d%%%s', profiles(i).level, tern(profiles(i).has_rp_snapshot,'','*'));
        end
    end
    lg = legend(le, 'Location','bestoutside'); 
    title(lg,'levels (*=no r_p in snapshot)');
    saveas(f1, fullfile(opt.OutDir, 'rp_radial_lines.png'));
end

% ——作图：热图（行=level，列=r/R，值=<r_p>）——
if opt.MakeHeatmap
    f2 = figure('Color','w','Position',[100 100 900 420]);
    imagesc(r_over_R, levels, H); axis xy;
    xlabel('r / R'); ylabel('overall conversion (%)');
    cb = colorbar; ylabel(cb, sprintf('mean r_p (%s)', ustr));
    title('Pore radius heatmap vs radius and overall conversion');
    saveas(f2, fullfile(opt.OutDir, 'rp_radial_heatmap.png'));
end

fprintf('完成：图像已保存到 %s/ ，返回 profiles 结构用于后续分析。\n', opt.OutDir);
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

% 小工具：三目字符串
function s = tern(tf, s1, s2)
if tf, s = s1; else, s = s2; end
end
