function vizPNM(P, rp, p1, p2, rt, R, outdir)
% 3D 几何可视化（pore: scatter3 大小映射 rp；throat: line 粗细映射 rt）

% 1) 散点 + 连线
% f1 = figure('Color','w');
f1 = figure('Color','w', 'Units','pixels', 'Position',[100 100 1000 1000]); % 正方形画布
hold on;
scatter3(P(:,1),P(:,2),P(:,3), 5 + 500*normalize01(rp), normalize01(rp), 'filled');
for i=1:numel(p1)
    plot3(P([p1(i) p2(i)],1), P([p1(i) p2(i)],2), P([p1(i) p2(i)],3), '-', ...
        'LineWidth', 0.5 + 4*normalize01(rt(i)), 'Color', [0 0 0 0.15]);
end
% 画球外框（线框）
[Xs,Ys,Zs] = sphere(24);
surf(R*Xs, R*Ys, R*Zs, 'FaceAlpha',0.03,'EdgeAlpha',0.05,'EdgeColor',[0.3 0.3 0.3],'FaceColor','none');
axis equal; grid off; view(3);
title('PNM Geometry: pores & throats');
xlabel('x'); ylabel('y'); zlabel('z');
% saveas(f1, fullfile(outdir,'PNM_geometry.png'));
exportgraphics(f1, fullfile(outdir,'PNM_geometry.png'));

% 2) 统计图
f2 = figure('Color','w');
subplot(1,2,1); histogram(rp*1e6, 40); xlabel('r_p [\mum]'); ylabel('count'); title('Pore radius PDF');
subplot(1,2,2); histogram(rt*1e6, 40); xlabel('r_t [\mum]'); ylabel('count'); title('Throat radius PDF');
saveas(f2, fullfile(outdir,'PNM_hist.png'));
end

function y = normalize01(x)
x = x(:);
mn = min(x); mx = max(x);
if mx>mn
    y = (x-mn)/(mx-mn);
else
    y = zeros(size(x));
end
end


