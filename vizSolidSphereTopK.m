function vizSolidSphereTopK(P, rp, p1, p2, rt, R, voxRes, outdir, K, alpha_solid, opts)
% 全部绘制（无裁剪、无立方体）：球表面 + 孔/喉表面
% 仅挖最大的 K 个孔及其相互连接的喉；只用 isosurface
%
% Inputs 同前；alpha_solid = 球表面透明度（0..1）
% opts 可选: .dpi(600) .smooth.box(5) .solidColor [.70 .74 .82]
%            .voidColor [.78 .80 .88] .voidAlpha(0.35)

if nargin < 9 || isempty(K), K = 100; end
if nargin < 10 || isempty(alpha_solid), alpha_solid = 0.22; end
if nargin < 11, opts = struct; end
if ~isfield(opts,'dpi'),         opts.dpi = 600; end
if ~isfield(opts,'smooth'),      opts.smooth = struct('box',5); end
if ~isfield(opts,'solidColor'),  opts.solidColor = [0.70 0.74 0.82]; end
if ~isfield(opts,'voidColor'),   opts.voidColor  = [0.78 0.80 0.88]; end
if ~isfield(opts,'voidAlpha'),   opts.voidAlpha  = 0.35; end
if ~exist(outdir,'dir'), mkdir(outdir); end

%% 1) 选择最大 K 孔 + 它们之间的喉
[~, idx] = sort(rp, 'descend');
K = min(K, numel(rp));
selPore   = false(size(rp)); selPore(idx(1:K)) = true;
selThroat = selPore(p1) & selPore(p2);
fprintf('Selected %d pores and %d throats.\n', K, nnz(selThroat));

%% 2) 体素网格（球域）
lin = linspace(-R, R, voxRes);
[x,y,z] = ndgrid(lin,lin,lin);
r2 = x.^2 + y.^2 + z.^2;
h  = lin(2)-lin(1);
solid = (r2 <= R^2);

%% 3) 挖孔与喉（仅所选）
for i = find(selPore).'
    cx = P(i,1); cy = P(i,2); cz = P(i,3); rr = rp(i);
    solid((x-cx).^2 + (y-cy).^2 + (z-cz).^2 <= rr^2) = false;
end
for e = find(selThroat).'
    A = P(p1(e),:); B = P(p2(e),:);
    v = B - A; L = norm(v); if L < eps, continue; end
    v = v / L;
    rel = [x(:)-A(1), y(:)-A(2), z(:)-A(3)];
    t   = rel * v'; t = max(0, min(L, t));
    closest = A + t .* v;
    d2 = sum(([x(:) y(:) z(:)] - closest).^2, 2);
    solid( reshape(d2 <= rt(e)^2, size(x)) ) = false;
end

void = ~solid;

%% 4) 平滑（可选）
Vs = double(solid);
Vv = double(void);
if opts.smooth.box > 0
    Vs = smooth3(Vs,'box',opts.smooth.box);
    Vv = smooth3(Vv,'box',opts.smooth.box);
end

%% 5) isosurface（无裁剪、无 isocaps）
iso = 0.5;
fv_solid = isosurface(x,y,z,Vs, iso); % 球面 + 孔壁/喉壁
fv_void  = isosurface(x,y,z,Vv, iso); % 孔/喉外表面

% 仅保留“外壳”三角面，剔除内部孔壁（孔壁交给 void 来显示）
if ~isempty(fv_solid) && ~isempty(fv_solid.vertices)
    r_vert = sqrt(sum(fv_solid.vertices.^2,2));
    r_thr  = max(R - 1.5*h, 0);
    maskF  = false(size(fv_solid.faces,1),1);
    for k = 1:size(fv_solid.faces,1)
        maskF(k) = mean(r_vert(fv_solid.faces(k,:))) > r_thr;
    end
    fv_solid.faces = fv_solid.faces(maskF,:);
end

%% 6) 绘制（仅球与孔/喉）
f = figure('Color','w','Units','pixels','Position',[100 100 1000 1000], ...
           'Renderer','opengl','GraphicsSmoothing','on');
ax = axes('Parent',f); hold(ax,'on');

if ~isempty(fv_solid) && ~isempty(fv_solid.faces)
    patch(ax, fv_solid, 'EdgeColor','none', ...
        'FaceColor', opts.solidColor, 'FaceAlpha', alpha_solid);
end
if ~isempty(fv_void) && ~isempty(fv_void.vertices)
    patch(ax, fv_void, 'EdgeColor','none', ...
        'FaceColor', opts.voidColor, 'FaceAlpha', opts.voidAlpha);
end

axis(ax,'equal'); axis off;
pbaspect(ax,[1 1 1]); view(ax,3);
grid(ax,'off'); set(ax,'XTick',[],'YTick',[],'ZTick',[]);

% camlight(ax,'headlight'); lighting(ax,'gouraud');
% 先把已有灯光删掉
delete(findall(gcf,'Type','light'));

% 柔和主光（偏上方，颜色偏暗）
light(ax,'Style','infinite','Position',[1 1 2], 'Color',[0.7 0.7 0.7]);

% 很弱的填充光，避免阴影太黑
light(ax,'Style','infinite','Position',[-2 -1 1], 'Color',[0.25 0.25 0.25]);

lighting(ax,'gouraud');   % 仍用 gouraud


exportgraphics(f, fullfile(outdir, sprintf('PNM_sphere_top%d.png',K)), ...
               'Resolution', opts.dpi);
end
