function size_factor_axial = computeSizeFactors(rp, p1, p2, rt, Lt, shape)
% 简化的尺寸因子：这里返回一个“轴向系数”占位（用于可视化/调参）
% 真实使用中可计算入口/出口/颈部三段的等效阻抗
% 目前 size_factor_axial 仅返回 shape 的轻微修正版本
size_factor_axial = shape .* (rt ./ max(Lt,eps));
end
