function axe = plot_layer(layer,axe)
% plot layer masks on the image
% layer size: [img_rows, img_cols]
%             or [img_rows, img_cols, layers]
%% plot params
colormaps = [
    0 0 0;
    0 1 0;
    0 0 1;
    1 1 0;
    1 0 1;
    0 1 1;
    0.2 0.8 0;
    0.2 0.2 0.6;
    1 0 0;
    0 0 0;];
alpha = 0.2;
%% plot
if nargin<2
    axe = axes();
end
axe.NextPlot = 'add';
sz = size(layer);
if length(sz)<3
    for j = 0:max(layer(:))
        masks = layer(:,:)==j;
        alphamask(masks, colormaps(j+1,:), alpha, axe);
    end
else
    for j = 1:sz(3)
       masks = squeeze(layer(:,:,j));
       alphamask(masks, colormaps(j,:), alpha, axe);
    end
end
end
