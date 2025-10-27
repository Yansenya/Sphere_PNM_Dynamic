function [p1, p2, Lgeom] = connectPores(P, z_target, kCand)
% 基于 kNN 候选 + 目标配位数（度受限）生成连边；再用 MST 保证连通
N = size(P,1);

% 候选邻：需要 knnsearch（统计工具箱）；若无，可退化 pdist2（慢）
try
    [idx, dist] = knnsearch(P, P, 'K', min(kCand, N));
catch
    D = pdist2(P,P);
    [dist, idx] = sort(D, 2, 'ascend');
    dist = dist(:,1:kCand);
    idx  = idx(:,1:kCand);
end

degree = zeros(N,1);
edges = []; Ls = [];

% 贪心/随机试探：按距离从近到远，填满每个节点的度
order = randperm(N);
for t=1:numel(order)
    i = order(t);
    nn = idx(i,2:end); % 跳过自己
    dd = dist(i,2:end);
    for k=1:numel(nn)
        j = nn(k);
        if i==j, continue; end
        if degree(i) >= z_target(i), break; end
        if degree(j) >= z_target(j), continue; end
        % 检查是否已有边
        if ~isempty(edges)
            if any( (edges(:,1)==i & edges(:,2)==j) | (edges(:,1)==j & edges(:,2)==i) )
                continue;
            end
        end
        edges = [edges; i j]; %#ok<AGROW>
        Ls    = [Ls; dd(k)]; %#ok<AGROW>
        degree(i) = degree(i)+1;
        degree(j) = degree(j)+1;
        if degree(i) >= z_target(i), break; end
    end
end

% 保证连通：构建 MST（基于所有候选的距离）
% 简化：用最近邻图的所有边（双向）做图，再 Kruskal
E_full = [];
for i=1:N
    for k=2:kCand
        j = idx(i,k);
        if i<j
            E_full = [E_full; i j]; %#ok<AGROW>
        end
    end
end
E_full = unique(sort(E_full,2),'rows');
L_full = sqrt(sum((P(E_full(:,1),:)-P(E_full(:,2),:)).^2,2));

% Kruskal MST
[~,ord] = sort(L_full,'ascend');
parent = 1:N;
findp = @(x) find_parent(x,parent);
mst = [];
for t=1:numel(ord)
    e = E_full(ord(t),:);
    a = findp(e(1)); b = findp(e(2));
    if a~=b
        parent(a) = b;
        mst = [mst; e]; %#ok<AGROW>
    end
end

% 合并现有边与 MST（去重）
allE = [edges; mst];
allE = unique(sort(allE,2),'rows');
p1 = allE(:,1);
p2 = allE(:,2);
Lgeom = sqrt(sum((P(p1,:)-P(p2,:)).^2,2));
end

function r = find_parent(x,parent)
while parent(x)~=x
    x = parent(x);
end
r = x;
end
