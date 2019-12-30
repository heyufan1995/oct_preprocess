function visualizer(data, options)
% Visualize the image
% data.img [img_rows, img_cols, Bscan]
% data.bds [img_cols,Bscans,boundaries]
% data.mask [img_rows, img_cols, Bscans]
% data.lesion [img_rows, img_cols, Bscans]

sz = size(data.img);
if isfield(options,'slice')
    slice = options.slice;
else
    slice = 1:sz(3);
end
axe = axes();
for i = slice
    axe.NextPlot = 'replaceall';
    imagesc(data.img(:,:,i),'Parent',axe);axis image;colormap gray;
    title(num2str(i));
    if isfield(data,'mask')
        plot_layer(data.mask(:,:,i),axe);
    end
    if isfield(data,'bds')
        plot_boundary(squeeze(data.bds(:,i,:)),axe);
    end  
    if isfield(options,'save') && options.save
        % bug exists here
        truesize;
        ct = clock;
        savename = sprintf('%4d-%d-%2d-%2d-%1d.jpg',ct(1),ct(2),ct(3),ct(4),i);
        F = getframe(axe.Parent);
        imwrite(frame2im(F), savename)
    end
    if isfield(options,'pause') && options.pause
        pause;
    end
end

