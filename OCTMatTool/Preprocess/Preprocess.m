function [data, record_params] = Preprocess(filename,options)
    % Preprocess the image and its coresponding manual delineation into patches 
    % Current version processes healthy control, MS and DME data.
    % Initially written for python pipeline.
    % Written by Yufan. 1/9/2018
    % Args:
    % filename: the file path to the scan
    % options: struct, see default_options.m or string, filepath to custom_options.m
    prev_path = pwd;
    curpath =  mfilename('fullpath');
    cd(fileparts(curpath));
    cd ..;
    addpath(genpath(pwd));
    cd(prev_path);
    if nargin<2
        try
            default_options;
        catch err
            fprintf('error: could not open file!\n Error message: %s\n\n',err.message);
        end
    else
        if ischar(options)
            opt = options;
            clear options;
            run(opt);
        end
    end
    preproc_params = options.preproc_params;
    preproc_params.retinadetector_type = options.types;
    %% Reading data
    fprintf('reading data:%s\n',filename);
    volfile = filename; 
    if strcmpi(options.types,'hc') || strcmpi(options.types,'mme')
        try
            [header, img_vol, ~, ~, scanner_type] = octReader(volfile);
        catch err
            fprintf('error: could not open file!\n Error message: %s\n\n',err.message);
            return
        end        
    elseif strcmpi(options.types,'dme')
        fprintf('preprocess DME data')
        S = load(volfile);
        img_vol = S.images;
        img_vol = double(dme_preprocess(img_vol))/255;
        % The DME data don't have header information
        load('dmeheader');
        header.SizeX = size(img_vol,2);
        header.SizeZ = size(img_vol,1);
        header.NumBScans = size(img_vol,3);
        scanner_type = 'spectralis';
    elseif strcmpi(options.types,'topcon')
        fprintf('preprocess topcon')
        S = load(volfile);
        img_vol = double(S.d3)/255;
        load('topconheader');
        header.SizeX = size(img_vol,2);
        header.SizeZ = size(img_vol,1);
        header.NumBScans = size(img_vol,3);
        scanner_type = 'spectralis';     
    elseif strcmpi(options.types,'bioptigen')
        fprintf('preprocess bioptigen')
        S = load(volfile);
        img_vol = double(S.images)/255;
        load('bioheader');
        header.SizeX = size(img_vol,2);
        header.SizeZ = size(img_vol,1);
        header.NumBScans = size(img_vol,3);
        scanner_type = 'spectralis';         
    end
    sz = size(img_vol);
    sz_ori = sz;
    idx = 1:sz(3);
    %% Reading manual
    if ~strcmp(options.segfile,'')
        if strcmpi(options.types,'hc') || strcmpi(options.types,'mme')
            S = load(options.segfile);   
            if isfield(S,'clims_x') && isfield(S,'clims_y')
                img_vol = img_vol(:,S.clims_x(1):S.clims_x(2),S.clims_y(1):S.clims_y(2));
                sz = size(img_vol);
                idx = 1:sz(3);          
            end
            if ~isfield(S,'bd_pts')
                [~,bd_pts] = readControlPoints(options.segfile,sz(2));
            else
                bd_pts = S.bd_pts;
            end
            % Need to modify when segmenting different boundary schemes
            if size(bd_pts,3)==11
                bd_pts(:,:,[3,11])=[];
            end
            lesion = false(sz);
            if isfield(S,'cyst_px')           
                for j = 1:length(S.cyst_px)
                    s = lesion(:,:,j);
                    s(S.cyst_px{j}) = true;
                    lesion(:,:,j) = s;
                end        
            elseif isfield(S,'mask_vol')
                lesion = S.mask_vol;
            end
            idx = sum(sum(bd_pts,1),3)>0;        
        elseif strcmpi(options.types,'dme') 
            bd_pts = S.manualLayers1;
            bd_pts = permute(bd_pts,[2,3,1]);
            lesion = S.manualFluid1;
            lesion(lesion>0.5) = 1;
            idx = nansum(nansum(bd_pts,1),3)>0; 
        elseif strcmpi(options.types,'bioptigen') 
            bd_pts = S.layerMaps;
            bd_pts = permute(bd_pts,[2,1,3]);
            lesion = false(sz);
            idx = sum(nansum(bd_pts,3) > 0, 1) > sz(2)*2/3;
            %idx = nansum(nansum(bd_pts,1),3)>0;      
           
        end
    end
    %% Only use data with manual labels and resize
    if exist('bd_pts','var')
        bd_pts = bd_pts(:,idx,:);
    end
    if exist('lesion', 'var')
        lesion = lesion(:,:,idx);
    end
    if isfield(options,'resize') && ~isempty(options.resize) ...
        && ~strcmp(scanner_type,'spectralis')
        sz = size(img_vol);
        if isscalar(options.resize)
            options.resize = [options.resize, 1000*header.ScaleX];
        end
        new_sz = [round(sz(1)*1000*header.ScaleZ/options.resize(1)),...
                  round(sz(2)*1000*header.ScaleX/options.resize(2)),sz(3)];
        img_vol = imresize(img_vol,new_sz(1:2));
        old_scale = [header.ScaleZ, header.ScaleX];
        header.ScaleZ = sz(1)*header.ScaleZ/size(img_vol,1);   
        header.ScaleX = sz(2)*header.ScaleX/size(img_vol,2); 
        header.SizeZ = size(img_vol,1);
        header.SizeX = size(img_vol,2);
        scale_ratio = [old_scale(1)/header.ScaleZ,old_scale(2)/header.ScaleX];
        record_params.scale_ratio = scale_ratio;
        if exist('bd_pts','var')
            bd_pts = bd_pts*scale_ratio(1);
            if scale_ratio(2) ~= 1
                %interpolate for every surface and bscan
                bd_pts = interp1(bd_pts,linspace(1,size(bd_pts,1),new_sz(2)));
            end
        end 
        if exist('lesion', 'var')
            lesion = imresize(lesion,new_sz(1:2));
        end
        
    end
    %% Normalize, flatten Data
    img_vol_org = img_vol(:,:,idx);
    px_buf = round(60/(header.ScaleZ*1000));
    [img_vol,retina_mask,~,shifts] = preprocessData(img_vol,header,preproc_params,scanner_type,1);
    img_vol = img_vol(:,:,idx);
    retina_mask = retina_mask(:,:,idx);
    shifts = shifts(:,idx);
    sz = size(img_vol);
    sc = sum(sum(retina_mask,3),2); 
    bv = find(sc>0,1,'last')+px_buf;
    tv = find(sc>0,1,'first');
    if tv<bv-options.img_rows + 1 && ~strcmpi(options.types,'bioptigen') 
        fprintf('Retina too thick\n');        
        if options.crop
            fprintf('options.flatvol is disabled \n');
            options.flatvol = false;
            tv = bv-options.img_rows + 1;
        else
            fprintf('Crop adjusted\n');  
            tv = bv - 2^min(4,(ceil(log2(bv-tv+1-options.img_rows)))) ...
                    - options.img_rows + 1;
        end
    else
        tv = bv-options.img_rows + 1;     
    end
    

    %% Normalize, flatten Label    
    if exist('bd_pts','var')
        bd_pts_org = bd_pts;
        bd_pts = bd_pts - repmat(shifts,[1,1,size(bd_pts,3)]);
    end
    if exist('lesion', 'var')
        lesion_org = lesion;
        lesion = retinaFlatten(lesion,shifts,'nearest');     
    end
    %% Crop image
    if options.crop
        img_rows = options.img_rows;
        img_cols = options.img_cols;   
        if  options.segs > 2*sz(2)/options.img_cols - 1
            if ~options.train           
                options.segs = round(2*sz(2)/img_cols) - 2;
                warning(['Too many patches each B-Scan,Reconstruction algorithm need to be updated\n', ...
                         'The new patch number is %d'], options.segs);
            end
        end
        segs = options.segs;
        X_train = zeros(segs*sz(3),img_rows,img_cols);
        dilation = img_rows*ones(1,segs*sz(3));
        if exist('bd_pts','var')
            Y_train_b = zeros(segs*sz(3),img_cols,size(bd_pts,3));
        end
        if exist('lesion','var')
            Y_train_lesion = zeros(segs*sz(3),img_rows,img_cols);  
        end    
        n = 1;
        for ii = 1:sz(3)
            for jj = round(linspace(1,sz(2)-img_cols+1,options.segs))
                sc = sum(sum(retina_mask(:,jj:jj+img_cols-1,ii),3),2); 
                temp = find(sc>0,1,'first');
                if temp<tv
                    % Need to interpolate because the retina is thicker than
                    % options.img_rows
                    temp = temp - px_buf;
                    img_temp = img_vol(temp:bv,jj:jj+img_cols-1,ii);
                    % dilation records the orginal image rows.
                    dilation(n) = bv - temp + 1;                            
                    [x,y] = meshgrid(1:img_cols,1:img_rows);
                    [xx,yy] = meshgrid(1:img_cols,linspace(1,img_rows,dilation(n)));
                    img_temp = interp2(xx,yy,img_temp,x,y);
                    X_train(n,:,:) = single(img_temp); 
                    if exist('bd_pts','var')
                        bd_temp = squeeze(bd_pts(jj:jj+img_cols-1,ii,:));   
                        bd_temp = bd_temp - temp + 1;
                        bd = diff(bd_temp,1,2);
                        bd_temp(:,2:end) = bd;                                      
                        bd_temp = img_rows/dilation(n)*bd_temp; 
                        Y_train_b(n,:,:) = bd_temp;
                    end
                    if exist('lesion','var')                   
                        les_temp = lesion(temp:bv,jj:jj+img_cols-1,ii);                                                              
                        les_temp = interp2(xx,yy,les_temp,x,y,'nearest');
                        les_temp(les_temp>0.5) = 1;                                       
                        Y_train_lesion(n,:,:) = les_temp;                    
                    end
                else
                    img_temp = img_vol(tv:bv,jj:jj+img_cols-1,ii);
                    X_train(n,:,:) = single(img_temp); 
                    if exist('bd_pts','var')
                        bd_temp = bd_pts(jj:jj+img_cols-1,ii,:);
                        bd_temp = bd_temp - tv + 1;
                        bd = diff(bd_temp,1,3);
                        bd_temp(:,:,2:end) = bd; 
                        Y_train_b(n,:,:) = bd_temp;
                    end
                    if exist('lesion', 'var')                    
                        les_temp = lesion(tv:bv,jj:jj+img_cols-1,ii);                    
                        Y_train_lesion(n,:,:) = les_temp;                    
                    end
                end            
                if options.dplot
                    % plot the results
                    figure(1);
                    imagesc(squeeze(X_train(n,:,:)));colormap gray;
                    alphamask(squeeze(Y_train_lesion(n,:,:)), [1,0,0], 1);
                    line(1:img_cols,cumsum(squeeze(Y_train_b(n,:,:))'));
                    pause;
                end
                n = n+1;
            end
        end 
        % save cropped data 
        data.X_train = X_train;
        if exist('Y_train_lesion','var')        
            data.Y_train_lesion = Y_train_lesion;
        end
        if exist('Y_train_b','var')
            data.Y_train_b = Y_train_b;
        end
        record_params.dilation = dilation;
        record_params.index = round(linspace(1,sz(2)-img_cols+1,options.segs));
    end
    % record parameters for reconstruction  
    record_params.sizes_ori = sz_ori;
    record_params.sizes = sz;
    record_params.tv = tv;
    record_params.bv = bv;
    record_params.shifts = shifts;
    record_params.options = options;
    
    % return data    
    data.header = header;

    if options.flatvol
        data.flat_vol = single(img_vol(tv:bv,:,:));
        if exist('bd_pts','var')
            data.bds = bd_pts - tv + 1;
        end
        if exist('lesion', 'var')
            data.lesion = lesion(tv:bv,:,:);     
        end 
        if options.dplot
            % plot the results
            for n = 1:size(data.flat_vol,3)
                figure(1);
                imagesc(squeeze(data.flat_vol(:,:,n)));colormap gray;
                alphamask(squeeze(data.lesion(:,:,n)), [1,0,0], 1);
                line(1:size(data.flat_vol,2),squeeze(data.bds(:,n,:))');
                pause;
            end
        end
    end
    if options.org
        % return unflattened data
        % [img_rows, img_cols, bscans]
        data.img_vol = img_vol_org;
        if exist('bd_pts_org','var') 
            data.bds_org = bd_pts_org;
        end
        if exist('lesion_org','var') 
            data.lesion_org = lesion_org;
        end
    end

end

