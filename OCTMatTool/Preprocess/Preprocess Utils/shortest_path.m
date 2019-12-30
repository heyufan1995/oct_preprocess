%     {{Caserel}}
%     Copyright (C) {{2013}}  {{Pangyu Teng}}
% 
%     This program is free software; you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation; either version 2 of the License, or
%     (at your option) any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License along
%     with this program; if not, write to the Free Software Foundation, Inc.,
%     51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
%


%% This script demonstrates how graph theory can be used to segment 
% retinal layers in optical coherence tomography images. The method is based on 
% 
% Chiu SJ et al, Automatic segmentation of seven retinal layers in SDOCT
% images congruent with expert manual segmentation, Optics Express, 2010;18(18);19413-19428
% Section 3.2
% link(pubmed): http://goo.gl/Z8zsY
% 
% USAGE:
% run the script by pressing F5.
% 
% I am working on a more comprehensive software package for computer-aided 
% segmentation of retinal layers in optical coherence tomography images, 
% which currently includes 1. automated segmentation of 6 reitnal layers and 
% 2. GUI for examination and manual correction of the automated segmentation. 
% It is called caserel and can be downloaded at my github page. http://goo.gl/yPqhPu
%
%
% $Revision: 1.0 $ $Date: 2013/01/23 21:00$ $Author: Pangyu Teng $
% $Revision: 1.1 $ $Date: 2013/09/15 21:00$ $Author: Pangyu Teng $
%                   Comment: simplified the script to detect only ILM and RPE
%

close all;clear all;clc;
warning off;

%% load data
    
%get path of an image.
folderPath = cd;
data = load([folderPath '/test.mat']);

% get the boundary probabilities
bds_prob = data.bds_prob;
bscan = data.img;
%% ### Section 3.2 Calculatge graph weights ### 

%% generate adjacency matrix, see equation 1.
for bds = 1:size(bds_prob,1)
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

        %display progress
        if indR < round(10*ind/szadjMW(1)/szadjMW(2))
            indR = round(10*ind/szadjMW(1)/szadjMW(2));
            fprintf('progress: %1.0f%% done, this may take a while...\n\n',100*indR/10);
        end

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
    [ dist,path{1} ] = graphshortestpath( adjMatrixW, 1, numel(imgNew(:)) );

    [pathX,pathY] = ind2sub(szImgNew,path{1});

    % get rid of first and last few points that is by the image borders
    pathX =pathX(gradient(pathY)~=0);
    pathY =pathY(gradient(pathY)~=0);
    %% visualize the detected boundaries, which are ilm and rpe
    imagesc(imgNew); axis image; colormap('gray'); hold on;
    plot(pathY,pathX,'g--','linewidth',2);
    pause;
end
