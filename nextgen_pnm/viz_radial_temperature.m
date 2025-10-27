function profiles = viz_radial_temperature(varargin)
% 可视化：不同整体转化率下，沿颗粒半径方向的温度 T(r) 分布
% 用法示例：
%   profiles = viz_radial_temperature('Levels',[0 20 40 60 80 99], 'NBins', 40, 'Mode','T');
% 参数（Name-Value）：
%   'PNMFile'     : 默认 'PNM.mat'
%   'ResultsFile' : 默认 'reaction_results.mat'
%   'Levels'      : 需要绘制的整体转化率百分数数组（默认自动挑常用档）
%   'NBins'       : 径向分桶数（默认 40）
%   'OutDir'      : 输出图目录（默认 'fig_radial_T'）
%   'MakeHeatmap' : 是否输出热图（默认 true）
%   'MakeLines'   : 是否输出多曲线图（默认 true）
%   'Mode'        : 'T' 画绝对温度；'DeltaT' 画 T - T_inf（默认 'T'）
%
% 返回：
%   profiles: 结构数组（每个 level 一项）
%       .level, .t, .Xtot
%       .r_mid, .r_over_R
%       .meanT, .q25, .q75            （单位：K）
%       .mode                          （'T' 或 'DeltaT'）

% profiles_T = viz_radial_temperature('Levels',[0 10 20 40 60 80 90 95 99 100], 'NBins',40, 'Mode','T');     % 或 'DeltaT'

% ---------- 解析参数 ----------
p = inputParser;
addParameter(p,'PNMFile','PNM.mat');
addParameter(p,'ResultsFile','reaction_results.mat');
% addParameter(p,'ResultsFile','reaction_results_V7.mat');
addParameter(p,'Levels',[]);
addParameter(p,'NBins',40);
addParameter(p,'OutDir','fig_radial_T');
addParameter(p,'MakeHeatmap',true,@islogical);
addParameter(p,'MakeLines',true,@islogical);
addParameter(p,'Mode','T',@(s) any(strcmpi(s,{'T','DeltaT'})));
parse(p,varargin{:});
opt = p.Results;
opt.Mode = upper(opt.Mode);

assert(isfile(opt.PNMFile),    '未找到 PNM 文件：%s', opt.PNMFile);
assert(isfile(opt.ResultsFile),'未找到结果文件：%s', opt.ResultsFile);

% ---------- 读取数据 ----------
S_pnm = load(opt.PNMFile);    PNM = S_pnm.PNM;
S_res = load(opt.ResultsFile); results = S_res.results;

xyz = PNM.P.coords;
r   = sqrt(sum(xyz.^2,2));    % 节点半径
R   = PNM.meta.R;
assert(R > 0, '颗粒半径 R 必须为正');

snap = results.snapshots;
has_level = arrayfun(@(s) ~isempty(s.level), snap);
avail_levels = [snap(has_level).level];

% 外界温度（用于 DeltaT）
T_inf = 1000;
if isfield(results,'params') && isfield(results.params,'phys') && isfield(results.params.phys,'T_inf')
    T_inf = results.params.phys.T_inf;
end

% 气体常数（用于从 p 与 C 回推 T 的兜底方案）
R_gas = 8.314462618;
if isfield(results,'params') && isfield(results.params,'phys') && isfield(results.params.phys,'T_inf') %#ok<*AND2>
    % 留空即可，R_gas 就用上面的值
end

% 默认 Levels
levels = opt.Levels;
if isempty(levels)
    default_levels = [0 10 20 40 60 80 90 95 99 100];
    levels = default_levels(ismember(default_levels, avail_levels));
end
levels = levels(ismember(levels, avail_levels));
assert(~isempty(levels), '所选 Levels 均不存在。可用 Levels: %s', mat2str(sort(avail_levels)));

N  = size(xyz,1);
NB = opt.NBins;
edges = linspace(0, R, NB+1);
r_mid = 0.5*(edges(1:end-1)+edges(2:end));
r_over_R = r_mid / R;

% 热图矩阵（可选）
if opt.MakeHeatmap
    heat_levels = sort(unique(avail_levels));
    H = nan(numel(heat_levels), NB);   % 行: level，列: 径向桶
end

% ---------- 逐 level 统计 ----------
profiles = struct('level',[],'t',[],'Xtot',[], ...
                  'r_mid',[],'r_over_R',[],'meanT',[],'q25',[],'q75',[], ...
                  'mode',[]);
profiles = repmat(profiles, 1, numel(levels));

for iL = 1:numel(levels)
    L = levels(iL);
    idxSnap = find(arrayfun(@(s) ~isempty(s.level) && s.level==L, snap), 1, 'first');
    S = snap(idxSnap);

    % —— 获取温度向量（优先使用直接保存的温度）——
    Ti = [];
    if isfield(S,'T') && ~isempty(S.T)
        Ti = S.T(:);
    elseif isfield(S,'T_nodes') && ~isempty(S.T_nodes)
        Ti = S.T_nodes(:);
    elseif isfield(S,'p_pa') && isfield(S,'C') && ~isempty(S.p_pa) && ~isempty(S.C)
        % 兜底：若快照里保存了浓度 C 与压强 p，可用 p = C*R*T 推回温度
        % 注意：此路径只有你在模拟时额外保存了 C 才可用
        Ti = S.p_pa(:) ./ (R_gas .* S.C(:));
    else
        error(['快照 level=%d 中未找到温度字段（T 或 T_nodes），也无法从 p 与 C 回推。\n' ...
               '建议在第四版模拟代码的快照写入处添加：\n' ...
               '    snap(idx).T = Tfield;    %% 保存温度（K）\n' ...
               '然后重新运行得到新的 reaction_results.mat。'], L);
    end

    assert(numel(Ti)==N, '快照 level=%d 的温度长度(%d)与节点数(%d)不符', L, numel(Ti), N);

    % 根据模式选择绘制 T 或 ΔT
    if strcmpi(opt.Mode,'DELTAT')
        Vi = Ti - T_inf;
    else
        Vi = Ti;
    end

    [mT, q25, q75] = radial_stats(r, Vi, edges);

    profiles(iL) = struct('level',L,'t',S.t,'Xtot',S.Xtot, ...
                          'r_mid',r_mid(:),'r_over_R',r_over_R(:), ...
                          'meanT',mT(:),'q25',q25(:),'q75',q75(:), ...
                          'mode',opt.Mode);
end

% 填热图（若可）
if opt.MakeHeatmap
    for il = 1:numel(heat_levels)
        L = heat_levels(il);
        idxSnap = find(arrayfun(@(s) ~isempty(s.level) && s.level==L, snap), 1, 'first');
        S = snap(idxSnap);

        Ti = [];
        if isfield(S,'T') && ~isempty(S.T)
            Ti = S.T(:);
        elseif isfield(S,'T_nodes') && ~isempty(S.T_nodes)
            Ti = S.T_nodes(:);
        elseif isfield(S,'p_pa') && isfield(S,'C') && ~isempty(S.p_pa) && ~isempty(S.C)
            Ti = S.p_pa(:) ./ (R_gas .* S.C(:));
        else
            % 若某些 level 缺温度，则该行留 NaN（热图会断开）
            continue
        end
        Vi = Ti;
        if strcmpi(opt.Mode,'DELTAT'), Vi = Ti - T_inf; end
        [mT, ~, ~] = radial_stats(r, Vi, edges);
        H(il,:) = mT;
    end
end

% ---------- 作图 ----------
if ~exist(opt.OutDir,'dir'), mkdir(opt.OutDir); end

ylab = 'Temperature (K)';
ttl  = 'Radial temperature profiles';
if strcmpi(opt.Mode,'DELTAT')
    ylab = '\DeltaT = T - T_\infty (K)';
    ttl  = 'Radial \DeltaT profiles';
end

% (1) 多曲线：T(r/R)（均值 + 四分位带）
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
        plot(pr.r_over_R, pr.meanT, '-', 'LineWidth', 1.8, 'Color', c);
    end
    xlabel('r / R'); ylabel(ylab); xlim([0 1]);
    title(ttl);
    legtxt = arrayfun(@(pr) sprintf('%d%%  (t=%.3g s)', pr.level, pr.t), profiles, 'UniformOutput',false);
    legend(legtxt, 'Location','bestoutside');
    savename = 'temperature_radial_lines.png';
    if strcmpi(opt.Mode,'DELTAT'), savename = 'deltaT_radial_lines.png'; end
    saveas(f1, fullfile(opt.OutDir, savename));
    % close(f1);
end

% (2) 热图：行=整体转化率，列=r/R，颜色=平均 T 或 ΔT
if opt.MakeHeatmap
    f2 = figure('Color','w','Position',[100 100 880 440]);
    imagesc(r_over_R, heat_levels, H); axis xy;
    xlabel('r / R'); ylabel('overall conversion (%)');
    cb = colorbar; ylabel(cb, ylab);
    title(strrep([ttl ' (heatmap)'],'profiles',''));
    savename = 'temperature_radial_heatmap.png';
    if strcmpi(opt.Mode,'DELTAT'), savename = 'deltaT_radial_heatmap.png'; end
    saveas(f2, fullfile(opt.OutDir, savename));
    % close(f2);
end

fprintf('完成：已输出 %s(r/R) 的曲线与热图到 %s/\n', opt.Mode, opt.OutDir);
end


% ====== 径向统计工具：均值与四分位 ======
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
