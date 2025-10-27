% plot_radius_distribution_multi.m
% 在 0,20,40,60,80,% 与 last 转化率下绘制喉道孔径(直径)分布

clear; clc;
assert(isfile('reaction_results.mat'), '未找到 reaction_results.mat');

S = load('reaction_results.mat');
results = S.results;
snaps   = results.snapshots;

% 收集可用快照（有 throat_radii 的）
is_ok = arrayfun(@(s) ~isempty(s) && isfield(s,'throat_radii') && ~isempty(s.throat_radii), snaps);
ok_idx = find(is_ok);
assert(~isempty(ok_idx), 'snapshots 中没有可用的 throat_radii。');

levels_all = arrayfun(@(s) s.level, snaps(ok_idx));
% 目标转化率（百分比）
target_levels_req = [0 20 40 60 80]; 

% 为每个目标转化率挑选“<=该值的最大 level”，没有则退回到最小可用 level
pick_idx = zeros(1, numel(target_levels_req));
for k = 1:numel(target_levels_req)
    L = target_levels_req(k);
    le_mask = levels_all <= L;
    if any(le_mask)
        best_level = max(levels_all(le_mask));
    else
        best_level = min(levels_all);  % 如果连 0% 都没，退回最小可用（极少见）
    end
    % 找到该 level 在 snaps 中的索引（可能出现重复 level，取首个）
    idx_in_ok = find(levels_all == best_level, 1, 'first');
    pick_idx(k) = ok_idx(idx_in_ok);
end

% last：选最后一个可用快照
idx_last = ok_idx(end);
pick_idx = [pick_idx, idx_last];

% 标注用的名字
labels = arrayfun(@(L) sprintf('X = %d%%', L), target_levels_req, 'UniformOutput', false);
labels{end+1} = 'X = last';

% 取对应的直径集合
diam_sets = cell(1, numel(pick_idx));
for i = 1:numel(pick_idx)
    s = snaps(pick_idx(i));
    d = 2 * s.throat_radii(:);                % 直径 = 2 * 半径
    d = d(isfinite(d) & d > 0);
    assert(~isempty(d), '选中的某个快照 throat 直径为空或无正值。');
    diam_sets{i} = d;
end

% 统一分箱（对数刻度）
d_all = vertcat(diam_sets{:});
dmin = min(d_all); dmax = max(d_all);
edges = logspace(log10(dmin*0.95), log10(dmax*1.05), 35);  % 分箱数可调

% 计算各曲线的 pdf
pdf_sets = cell(size(diam_sets));
centers = sqrt(edges(1:end-1).*edges(2:end));
for i = 1:numel(diam_sets)
    pdf_sets{i} = histcounts(diam_sets{i}, edges, 'Normalization','pdf');
end

% 作图
figure('Color','w'); hold on;
for i = 1:numel(pdf_sets)
    plot(centers, pdf_sets{i}, 'LineWidth', 2);
end
set(gca,'XScale','log'); grid on;
xlabel('喉道孔径直径 d_t (m)');
ylabel('概率密度 (PDF)');
title('喉道孔径分布：X = 0, 20, 40, 60, 80, last');
legend(labels, 'Location','best');

% 在命令行打印每条曲线的一些统计
for i = 1:numel(diam_sets)
    d = diam_sets{i};
    s = snaps(pick_idx(i));
    fprintf('%-10s | t = %.3g s, Xt = %.3f | median = %.3g m, mean = %.3g m, N = %d\n', ...
        labels{i}, s.t, s.Xtot, median(d), mean(d), numel(d));
end

% 保存图片
outname = 'throat_diameter_distribution_X_0_20_40_60_80_last.png';
% exportgraphics(gcf, outname, 'Resolution', 300);
fprintf('已保存图像：%s\n', outname);
