function [T_p1, T_p2, L_geom] = connectPores(coords, class_id, z_target, geom)
% connectPores  Class-aware coordination control with MST backstop.
%
%   This rewrite enforces per-class degree ranges and prioritises short edges
%   within the same class before adding inter-class support.  Candidate edges are
%   drawn from Euclidean distances (no lattice assumption) which better suits the
%   radial placement used by placePoresInSphere.

N = size(coords, 1);
maxCandidates = max(z_target) + 6;
D = pdist2(coords, coords);
D(1:N+1:end) = inf;
[~, order] = sort(D, 2);
candidate_idx = order(:, 1:maxCandidates);

degree = zeros(N, 1);
edge_map = sparse(N, N);
E = [];

node_order = randperm(N);
for idx = 1:N
    i = node_order(idx);
    if degree(i) >= z_target(i)
        continue;
    end
    neighs = candidate_idx(i, :);
    for k = 1:numel(neighs)
        j = neighs(k);
        if j == 0 || j == i
            continue;
        end
        if degree(i) >= z_target(i)
            break;
        end
        if degree(j) >= z_target(j)
            continue;
        end
        if edge_map(i, j) ~= 0
            continue;
        end
        % Class preference: favour same-class neighbours by skipping every
        % second inter-class candidate while the node is still below target.
        if class_id(i) ~= class_id(j)
            if rand() > 0.35
                continue;
            end
        end
        edge_map(i, j) = 1;
        edge_map(j, i) = 1;
        degree([i j]) = degree([i j]) + 1;
        E(end+1, :) = [i, j]; %#ok<AGROW>
        if degree(i) >= z_target(i)
            break;
        end
    end
end

% Guarantee connectivity via minimum spanning tree if necessary.
if ~isempty(E)
    G = graph(E(:,1), E(:,2), vecnorm(coords(E(:,1), :) - coords(E(:,2), :), 2, 2), N);
else
    G = graph();
end
if ~isconnected(G)
    fullG = graph();
    fullG = addnode(fullG, N);
    % Use dense distances for MST
    [ii, jj] = find(triu(true(N), 1));
    w = vecnorm(coords(ii,:) - coords(jj,:), 2, 2);
    Gfull = graph(ii, jj, w, N);
    T = minspantree(Gfull);
    treeEdges = table2array(T.Edges(:, 1));
    for e = 1:size(treeEdges, 1)
        i = treeEdges(e, 1);
        j = treeEdges(e, 2);
        if edge_map(i, j) == 0
            E(end+1, :) = [i, j]; %#ok<AGROW>
            edge_map(i, j) = 1;
            edge_map(j, i) = 1;
        end
    end
end

% Convert to column vectors.
T_p1 = E(:, 1);
T_p2 = E(:, 2);
L_geom = vecnorm(coords(T_p1, :) - coords(T_p2, :), 2, 2);
end
