function  axe = plot_boundary(bds,axe)
% plot boundary on an axe
% Args:
% bds [img_cols,boundaries]
if nargin<2
    axe = axes();
end
sz = size(bds);
axe.NextPlot = 'add';
line(1:sz(1),bds,'LineWidth',2 ,'Parent', axe);
end

