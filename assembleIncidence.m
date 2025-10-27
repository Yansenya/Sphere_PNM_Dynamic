function [B, K] = assembleIncidence(N, p1, p2, G)
% 关联矩阵 B (M×N)，每条喉一行；导通率向量 K（对角）
M = numel(p1);
I = (1:M)'; 
B = sparse([I; I], [p1; p2], [ones(M,1); -ones(M,1)], M, N);
K = G(:);
end
