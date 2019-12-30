function [bds_pred_sp] = dijkstra_oct(data)
% shortest path algorithms based on 
% modified from https://github.com/pangyuteng/caserel/blob/master/octGraphTheorySimplified.m
% Chiu SJ et al, Automatic segmentation of seven retinal layers in SDOCT
% images congruent with expert manual segmentation, Optics Express, 2010;18(18);19413-19428
% Author: Yufan He 
% The Jonhs Hopkins Univ.
% 1/2/2019


%% load data
bds_prob = data.bds_prob;
bds_pred_sp = zeros(size(data.bds_pred));
%% generate adjacency matrix, see equation 1.
parfor bds = 1:size(bds_prob,1)
    % process each boundary indivisually
    img = squeeze(bds_prob(bds,:,:));
    szImg = size(img);
    imgNew = zeros([szImg(1) szImg(2)+2]);
    imgNew(:,2:1+szImg(2)) = img;
    szImgNew = [szImg(1) szImg(2)+2];    
    %minimum weight
    minWeight = 1E-5;
    %arry to store weights
    adjMW = nan([numel(imgNew(:)),8]);
    %arry to store point A locations
    adjMX = nan([numel(imgNew(:)),8]);
    %arry to store point B locations
    adjMY = nan([numel(imgNew(:)),8]);

    neighborIter = [1 1  1 0  0 -1 -1 -1;...
                    1 0 -1 1 -1  1  0 -1];

    %fill in the above arrays according to Section 3.2
    szadjMW = size(adjMW);
    ind = 1; indR = 0;
    
    while ind ~= szadjMW(1)*szadjMW(2) %this step can be made more efficient to increase speed.
        [i, j] = ind2sub(szadjMW,ind);    
        [iX,iY] = ind2sub(szImgNew,i);    
        jX = iX + neighborIter(1,j);
        jY  = iY + neighborIter(2,j);
        if jX >=1 && jX <= szImgNew(1) && jY >=1 && jY <= szImgNew(2)
             %save weight
             % set to minimum if on the sides
             if jY == 1 || jY == szImgNew(2)
                adjMW(i,j) = minWeight;
             % else, calculate the actual weight based on equation 1.
             else
                adjMW(i,j) = 2 - imgNew(iX,iY) - imgNew(jX,jY) + minWeight;
             end
            %save the subscript of the corresponding nodes
            adjMX(i,j) = sub2ind(szImgNew,iX,iY);
            adjMY(i,j) = sub2ind(szImgNew,jX,jY);
        end
        ind = ind+1;
    end
    

    %assemble the adjacency matrix
    keepInd = ~isnan(adjMW(:)) & ~isnan(adjMX(:)) & ~isnan(adjMY(:)) ;
    adjMW = adjMW(keepInd);
    adjMX = adjMX(keepInd);
    adjMY = adjMY(keepInd);

    %sparse matrices, based on eq 1 with the gradient,
    adjMatrixW = sparse(adjMX(:),adjMY(:),adjMW(:),numel(imgNew(:)),numel(imgNew(:)));
                        % and the invert of gradient.

    
    %% get shortest path 

    % get layer going from dark to light
    [dist,path ] = graphshortestpath( adjMatrixW, 1, numel(imgNew(:)) );

    [pathX,pathY] = ind2sub(szImgNew,path);
    % get rid of first and last few points that is by the image borders
    pathX =pathX(gradient(pathY)~=0);
    pathY =pathY(gradient(pathY)~=0);
    
    bds_pred_sp(1,bds,:) = pathX(2:end-1);
%     %% visualize the detected boundaries, which are ilm and rpe
%     imagesc(imgNew); axis image; colormap('gray'); hold on;
%     plot(pathY,pathX,'g--','linewidth',2);
%     pause;
end
end

