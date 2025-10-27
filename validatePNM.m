function QC = validatePNM(P, p1, p2, rp, rt, R, z_vals, z_target, alpha)
% 质量控制与统计：度分布、连通分量、约束违规、径向最短通道等
N = size(P,1);
M = numel(p1);

% 度
deg = accumarray([p1; p2], 1, [N,1]);
z_hist = histcounts(deg, [z_vals-0.5, z_vals(end)+0.5]);

% 约束违规（rt <= alpha*min(rp_i,rp_j)）
viol = rt > alpha * min(rp(p1), rp(p2));
violations = sum(viol);

% 连通分量（基于图）
G = graph(p1, p2);
[bins, bin_sizes] = conncomp(G);
n_comp = max(bins);
largest = max(bin_sizes);
frac_in_gcc = largest / N;

% 径向最短路径（中心→表面任一点）
[~, cidx] = min(sum(P.^2,2));            % 近似中心点
radii = sqrt(sum(P.^2,2));
surf_nodes = find( (R - radii) <= (0.75*mean(pdist(P))) );
if isempty(surf_nodes)
    surf_nodes = find( (R - radii) <= (0.1*R) );
end
dists = distances(G, cidx, surf_nodes);
minpath = min(dists(~isinf(dists)));
if isempty(minpath) || isinf(minpath)
    minpath = NaN;
end

QC = struct();
QC.deg = deg;
QC.z_hist = z_hist;
QC.n_components = n_comp;
QC.frac_in_gcc = frac_in_gcc;
QC.violations = violations;
QC.min_center_to_surface_hops = minpath;

fprintf('[QC] N=%d, M=%d, components=%d, frac_in_GCC=%.4f, viol=%d, min_hops=%g\n',...
    N, M, n_comp, frac_in_gcc, violations, minpath);
end
