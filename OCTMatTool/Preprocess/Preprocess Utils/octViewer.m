function octViewer(octvol,voxel_size,seg_pts)
% octViewerAdvanced.m
% 
%   GUI to display OCT data and segmentation results.
%
%   To display an OCT volume:
%       - Run octViewer without any arguments and go to File->Load Volume
%         and choose the desired format. Compatible formats include the
%         .vol format from Spectralis scanners, the .img format from
%         Cirrus (filename format should be XX_cube_z.img or
%         XX_cube_raw.img), or a saved MATLAB .mat file containing the
%         variable 'img_vol' representing the data as a 3 dimensional
%         array.
%       - Run octViewer with the first argument being a 3 dimensional array
%         representing the data (ex. a 496x1024x49 array containing 496
%         pixels per A-scan, 1024 A-scans per B-scan, and 49 B-scans per
%         volume).
%   To display a segmentation result:
%       - Run octViewer without any arguments and go to Segmentation->Load
%         File and choose the desired format. Compatible formats include a
%         .txt file from the OCT Toolbox manual segmentation tool, or a
%         .mat file stored after running the OCTLayerSegmentation function
%         in MATLAB.
%       - Run octViewer with the first argument being a 3 dimensional array
%         containing the position of each boundary at every A-scan (ex. a
%         1024x49x9 array containing the boundary position along each of
%         1024 A-scans per B-scan, 49 B-scans, and 9 boundaries).

% Update log:
% v0.4 - 9/16/2014
%    - New intensity normalization methods
%    - MME filtering - for Cirrus data
%    - Added slider for changing image
%    - Added new mode for manual segmentation (freehand)
%    - Added up and down arrow hotkeys to move the current point up or down
%    - Added hotkeys 'n' and 'p' to move to the next and previous point
%    - Added checkbox 'RB' to show automatic retina boundaries
%    - Read Bioptigen data
% v0.3 - 3/24/2014
%    - Can modify the contrast of the image based on a local area
% v0.2 - X/XX/2014
%    - Can use left and right arrow keys to change the image
%    - Can use the 'h' key to hide and show the segmentation
%    - Added pseudocyst segmentation 

% TODO:
%   Features:
%       * layer segmentation tools (including choroid):
%           * manual layer segmentation (load and save function, tag with initials of rater)
%           * automatic layer segmentation (multiple methods)
%       * cyst tools:
%           * automatic segmentation
%           * fundus density
%       * fundus image (SLO if it exists)
%       * display gradient image
%       * flatten
%       * fix aspect ratio
%
%   Bug fixes:
%       * prevent Pan tool from moving outside image
%       * setting "FontUnits 'normalized'" to all uicontrols breaks linux
%
% Andrew Lang (alang9@jhu.edu)
% Last updated: 9/16/2014
% $Id$

%% Input handling

if nargin < 3
    seg_pts = [];
end
if nargin < 2
    voxel_size = [];
end
if nargin < 1
    octvol = [];
end

% Check boundary segmentation size
if ~isempty(seg_pts)
    szv = size(octvol); % ex. [496 1024 49]
    sz = size(seg_pts); % ex. [1024 49 9]
    if (sz(1) ~= szv(2)) || (sz(2) ~= szv(3))
        error('Invalid size for segmentation points')
    end
end

%%  Initialization tasks

figure_name = 'OCT Viewer v0.3';

% Build figure
hf = figure('Visible','on',...
            'Name',figure_name,...
            'Resize','on',...
            'MenuBar','none',...
            'NumberTitle','off',...
            'Toolbar','figure',...
            'Units','normalized',... % normalized - beware - different size on different monitors
            'Position',[0.1 0.1 0.8 0.8]); % [left bottom width height]
    
control_pts = [];
bd_names = {'ILM','RNFL-GCL','GCL-IPL','IPL-INL','INL-OPL','OPL-ONL','ELM','IS-OS','OS-RPE','RPE-Choroid','Choroid-Sclera'};
active_bd = 1;
active_pt = 1;
startx = 1;
starty = 1;
fundus_proj = [];
cyst_px = [];
cyst_mask = [];
filename_vol = [];
filename_seg = [];
filename_mme = [];
savefile_mme = [];
savefile_seg = [];
bsc_header = [];
curimg = 0;
cystvals = [1 1 0]; % color of cyst overlay
minv = 0;
maxv = 1;
rpe_bdry_qf = []; % quick find boundary
ilm_bdry_qf = []; % quick find boundary
segmode = 1;

% Handles to plot objects
h_img = [];
h_mask = [];
himg_fun = [];
h_pts = zeros(11,1);
h_pts_rb = zeros(2,1);
h_ctrl = [];
h_activept = [];

% Custom pointer for segmentation


%%  Construct the components

% Axis for the image data
pax_vol = [0.05 0.2 0.9 0.75];
pax_vol_fun = [0.05 0.2 0.74 0.75];
ha = axes('Parent',hf,...
          'Position',pax_vol_fun,...
          'XTick',[],...
          'YTick',[],...
          'box','on',...
          'DrawMode','fast');
      
% Axis for the image mask overlay (MME)
hm = axes('Parent',hf,...
          'Position',pax_vol_fun,...
          'XTick',[],...
          'YTick',[],...
          'box','on',...
          'DrawMode','fast',...
          'HitTest','off');
      
linkaxes([ha, hm])
      
% Context menu for right clicking on image
% uicontextmenu

% Axis for the fundus image
pax_fun = [0.81 0.7 0.15 0.25];
ha_fun = axes('Parent',hf,...
              'Position',pax_fun,...
              'XTick',[],...
              'YTick',[],...
              'box','on',...
              'DrawMode','fast',...
              'HitTest','off');
          
% Scroll bar for images
ps_sl = [0.952 0.2 0.01 0.75];
ps_sl_fun = [0.792 0.2 0.01 0.75];
hs_sl = uicontrol('Style','slider',...
                  'Units','normalized',...
                  'Position',ps_sl_fun,...
                  'Min',1,...
                  'Max',10,...
                  'Value',5,...
                  'SliderStep',[.1 .1],...
                  'Enable','off');

% Text below image to say 'image x of y'
ptxt_img = [.87 .175 .03 .02];
ptxt_img_fun = [.71 .175 .03 .02];
htxt_img = uicontrol('Style','text',...
                     'Units','normalized',...
                     'Position',ptxt_img_fun, ...
                     'FontSize',8, ...
                     'BackgroundColor',get(hf,'Color'),...
                     'String','image');
      
% Edit box to control the current bscan index
pedit_sl = [.9 .175 .02 .02];
pedit_sl_fun = [.74 .175 .02 .02];
hedit_sl = uicontrol('Style','edit', ...
                      'Units','normalized',...
                      'Position',pedit_sl_fun, ...
                      'FontSize',8, ...
                      'BackgroundColor','white', ...
                      'String','-');

ptxt_of = [.918 .175 .01 .02];
ptxt_of_fun = [.758 .175 .01 .02];
htxt_of = uicontrol('Style','text',...
                    'Units','normalized',...
                    'Position',ptxt_of_fun, ...
                    'FontSize',8, ...
                    'BackgroundColor',get(hf,'Color'),...
                    'String','of');

% Text box displaying total number of images in the volume
ptxt_sl = [.93 .175 .02 .02];
ptxt_sl_fun = [.77 .175 .02 .02];
htxt_sl = uicontrol('Style','text', ...
                    'Units','normalized',...
                    'Position',ptxt_sl_fun, ...
                    'FontSize',8, ...
                    'BackgroundColor','white', ...
                    'String','-');
                
% BUTTON GROUP to enable different display modes (display, boundary
%   segmentation, pseudocyst segmentation)
hbg_mode = uibuttongroup('Title','Display mode',...
                         'Position',[0.05 0.03 0.13 0.14]);
                     
htb_disp = uicontrol('Style','togglebutton',...
                     'Parent',hbg_mode,...
                     'String','Display only',...
                     'Units','normalized',...
                     'Position',[0.05 0.685 0.9 0.27],...
                     'BackgroundColor','w');
                 
htb_bdry = uicontrol('Style','togglebutton',...
                     'Parent',hbg_mode,...
                     'String','Boundary segmentation',...
                     'Units','normalized',...
                     'Position',[0.05 0.365 0.9 0.27],...
                     'BackgroundColor',get(hf,'Color'),...
                     'Enable','off');
                 
htb_mme = uicontrol('Style','togglebutton',...
                     'Parent',hbg_mode,...
                     'String','Pseudocyst segmentation',...
                     'Units','normalized',...
                     'Position',[0.05 0.045 0.9 0.27],...
                     'BackgroundColor',get(hf,'Color'),...
                     'Enable','off');
                
% UIPANEL containing display options
hp_opts = uipanel('Title','Display options',...
                  'Position',[0.19 0.03 0.15 0.14]);

% Check box to show fundus projection
hcb_fun = uicontrol('Style','checkbox',...
                    'String','Show fundus/projection image',...
                    'Parent',hp_opts,...
                    'Units','normalized',...
                    'Position',[0.05 0.7 0.9 0.2],...%[0.2 0.15 0.13 0.02],...
                    'FontSize',8,...
                    'BackgroundColor','white',...
                    'Value',1,...
                    'Enable','on');
                
% Check box to hide the segmentation boundaries
hcb_hideseg = uicontrol('Style','checkbox',...
                    'String','Hide boundaries',...
                    'Parent',hp_opts,...
                    'Units','normalized',...
                    'Position',[0.05 0.4 0.9 0.2],...%[0.2 0.12 0.09 0.02],...
                    'FontSize',8,...
                    'BackgroundColor','white',...
                    'Enable','off');
                
% Check box to hide the retina boundaries
hcb_hiderb = uicontrol('Style','checkbox',...
                    'String','RB',...
                    'Parent',hp_opts,...
                    'Units','normalized',...
                    'Position',[0.7 0.4 0.2 0.2],...%[0.2 0.12 0.09 0.02],...
                    'FontSize',8,...
                    'BackgroundColor','white',...
                    'Enable','off');
                
% Check box to hide the pseudocyst segmentation
hcb_hidemme = uicontrol('Style','checkbox',...
                    'String','Hide pseudocysts',...
                    'Parent',hp_opts,...
                    'Units','normalized',...
                    'Position',[0.05 0.1 0.9 0.2],...%[0.2 0.12 0.09 0.02],...
                    'FontSize',8,...
                    'BackgroundColor','white',...
                    'Enable','off');
                 
% UIPANEL for displaying Segmentation options
hp_opts = uipanel('Title','Segmentation options',...
                  'Position',[0.35 0.03 0.3 0.14]);
              
% Popup menu for selecting the active boundary
uicontrol('Style','edit',...
          'Parent',hp_opts,...
          'Units','normalized',...
          'Position',[0.05 0.78 0.225 0.16], ...
          'FontSize',8, ...
          'BackgroundColor','w',...
          'String','Active boundary',...
          'HorizontalAlignment','left',...
          'Enable','inactive');
hp_active = uicontrol('Style','popupmenu',...
                      'Parent',hp_opts,...
                      'String',bd_names,...
                      'Units','normalized',...
                      'Position',[0.05 0.63 0.3 0.16],...
                      'FontSize',8,...
                      'BackgroundColor',get(hf,'Color'),...
                      'Enable','off');

% Segmentation mode - control points or drawing
uicontrol('Style','edit',...
          'Parent',hp_opts,...
          'Units','normalized',...
          'Position',[0.05 0.33 0.225 0.16], ...
          'FontSize',8, ...
          'BackgroundColor','w',...
          'String','Draw mode',...
          'HorizontalAlignment','left',...
          'Enable','inactive');
hp_segdraw = uicontrol('Style','popupmenu',...
                      'Parent',hp_opts,...
                      'String',{'Control points','Freehand'},...
                      'Units','normalized',...
                      'Position',[0.05 0.18 0.3 0.16],...
                      'FontSize',8,...
                      'BackgroundColor',get(hf,'Color'),...
                      'Enable','off');
                  
% Button to clear current boundary
hpb_clearseg = uicontrol('Style','pushbutton',...
                         'Parent',hp_opts,...
                         'Units','normalized',...
                         'Position',[0.45 0.55 0.5 0.35], ...
                         'FontSize',8,...
                         'String','Clear active boundary',...
                         'Enable','off');

% Button to clear current pseudocysts
hpb_clearcyst = uicontrol('Style','pushbutton',...
                          'Parent',hp_opts,...
                          'Units','normalized',...
                          'Position',[0.45 0.05 0.5 0.35], ...
                          'FontSize',8,...
                          'String','Clear active pseudocysts',...
                          'Enable','off');                  
                
% UIPANEL for displaying save options
hp_save = uipanel('Title','Save',...
                  'Position',[0.66 0.03 0.13 0.14]);

% Button to save layer segmentation
hpb_saveseg = uicontrol('Style','pushbutton',...
                         'Parent',hp_save,...
                         'Units','normalized',...
                         'Position',[0.05 0.55 0.8 0.35], ...
                         'FontSize',8,...
                         'String','Save layer segmentation',...
                         'Enable','off',...
                         'Tag','save');
                     
% Button to 'save as' layer segmentation
hpb_savesegas = uicontrol('Style','pushbutton',...
                         'Parent',hp_save,...
                         'Units','normalized',...
                         'Position',[0.86 0.55 0.12 0.35], ...
                         'FontSize',8,...
                         'String','...',...
                         'Enable','off',...
                         'Tag','saveas');
                     
% Button to save MME segmentation
hpb_savemme = uicontrol('Style','pushbutton',...
                         'Parent',hp_save,...
                         'Units','normalized',...
                         'Position',[0.05 0.05 0.8 0.35], ...
                         'FontSize',8,...
                         'String','Save MME segmentation',...
                         'Enable','off',...
                         'Tag','save');
                     
% Button to 'save as' MME segmentation
hpb_savemmeas = uicontrol('Style','pushbutton',...
                         'Parent',hp_save,...
                         'Units','normalized',...
                         'Position',[0.86 0.05 0.12 0.35], ...
                         'FontSize',8,...
                         'String','...',...
                         'Enable','off',...
                         'Tag','saveas');
                     
% Text box to display debug output
ptxt_db = [.4 .175 .2 .02];
ptxt_db_fun = [0.32 0.175 0.2 0.02];
htxt_db = uicontrol('Style','text',...
                    'Units','normalized',...
                    'Position',ptxt_db_fun, ...
                    'FontSize',8, ...
                    'BackgroundColor','w',...%get(hf,'Color'),...
                    'String','');
                
% Text box to show image quality
ptxt_snr = [0.63 0.175 0.06 0.02];
ptxt_snr_fun = [0.55 0.175 0.06 0.02];
htxt_snr = uicontrol('Style','text',...
                     'Units','normalized',...
                     'Position',ptxt_snr_fun,...
                     'FontSize',8,...
                     'BackgroundColor',get(hf,'Color'),...
                     'String','image quality');
                 
ptxt_snr_val = [0.69 0.175 0.02 0.02];
ptxt_snr_val_fun = [0.61 0.175 0.02 0.02];                 
htxt_snr_val = uicontrol('Style','text',...
                         'Units','normalized',...
                         'Position',ptxt_snr_val_fun,...
                         'FontSize',8,...
                         'BackgroundColor','w',...
                         'String','-');
                
% Toolbar - only want Zoom In, Zoom Out, Pan
htb = findall(hf,'type','UIToolbar');
hb = findobj(allchild(htb),'-not','ToolTipString','Zoom In',...
                           '-not','ToolTipString','Zoom Out',...
                           '-not','ToolTipString','Pan');
set(hb,'Visible','off')

% Menu bar
hmf = uimenu(hf,'Label','File');
hmf_load = uimenu(hmf,'Label','Load volume');
uimenu(hmf_load,'Label','Spectralis (.vol) file',...
                'Callback',@load_vol,...
                'Tag','vol')
uimenu(hmf_load,'Label','Cirrus (.img) file',...
                'Callback',@load_vol,...
                'Tag','img')
uimenu(hmf_load,'Label','Bioptigen (.oct) file',...
                'Callback',@load_vol,...
                'Tag','oct')
uimenu(hmf_load,'Label','.mat file',...
                'Callback',@load_vol,...
                'Tag','mat')

hms = uimenu(hf,'Label','Layer Segmentation','Enable','off');
hms_load = uimenu(hms,'Label','Load file');
uimenu(hms_load,'Label','OCT Toolbox file',...
                'Callback',@load_seg,...
                'Tag','otb')
uimenu(hms_load,'Label','.mat file',...
                'Callback',@load_seg,...
                'Tag','mat')
            
hmm = uimenu(hf,'Label','Pseudocyst Segmentation','Enable','off');
hmm_load = uimenu(hmm,'Label','Load file');
uimenu(hmm_load,'Label','.mat file',...
                'Callback',@load_mmeseg,...
                'Tag','mat')
            
hmip = uimenu(hf,'Label','Image Processing','Enable','off');
hmin = uimenu(hmip,'Label','Intensity normalization');
uimenu(hmin,'Label','Robust contrast stretching',...
           'Callback',@normalize_data,...
           'Tag','rcs');
uimenu(hmin,'Label','Attenuation correction',...
           'Callback',@normalize_data,...
           'Tag','ac');
uimenu(hmin,'Label','RPE/Vitreous equalization',...
           'Callback',@normalize_data,...
           'Tag','r_eq');
uimenu(hmin,'Label','RNFL/RPE equalization',...
           'Callback',@normalize_data,...
           'Tag','rr_eq');
uimenu(hmin,'Label','Full RPE equalization',...
           'Callback',@normalize_data,...
           'Tag','r_eq3');
uimenu(hmin,'Label','Smooth RPE equalization',...
           'Callback',@normalize_data,...
           'Tag','r_eq2');
uimenu(hmin,'Label','B-scan RPE equalization',...
           'Callback',@normalize_data,...
           'Tag','r_eq4');
       
uimenu(hmip,'Label','Median Filter',...
                        'Callback',@median_filter_data);
                    
hmmf = uimenu(hmip,'Label','MME Filter');
uimenu(hmmf,'Label','Single image',...
                    'Callback',@MME_filter_data,...
                    'Tag','s');
uimenu(hmmf,'Label','Whole volume',...
                    'Callback',@MME_filter_data,...
                    'Tag','v');
                    
uimenu(hmip,'Label','Experimental Filter',...
                        'Callback',@Experimental_filter);
                    
hmc = uimenu(hmip,'Label','Contrast');
uimenu(hmc,'Label','Contrast window',...
           'Callback',@contrast_window);
uimenu(hmc,'Label','Set contrast range',...
           'Callback',@contrast_range);
uimenu(hmc,'Label','Reset contrast',...
           'Callback',@reset_contrast);

%%  Initialization tasks

if ~isempty(octvol)

    initialize_vol();
    
    if ~isempty(seg_pts)
        control_pts = cell(size(seg_pts,2),11);
        for ii = 1:size(seg_pts,2)
            for jj = 1:size(seg_pts,3)
                control_pts{ii,jj} = [1:size(seg_pts,1); seg_pts(:,ii,jj)']';
            end
            control_pts{ii,10} = [];
            control_pts{ii,11} = [];
        end
        set(hp_segdraw,'Value',2);
        
        updateseg();
    end
end

%%  Callbacks for MYGUI

set(hf,'CloseRequestFcn',@figure_close);
set(hf,'WindowScrollWheelFcn',@scrollfcn);
set([hedit_sl, hs_sl],'Callback',@change_image);
set([hcb_hideseg hcb_hiderb],'Callback',@hide_seg);
set(hcb_hidemme,'Callback',@hide_mme);
set(hcb_fun,'Callback',@show_fundus);
set(hf,'KeyPressFcn',@key_press);
set(hbg_mode,'SelectionChangeFcn',@change_mode);
set(hp_active,'Callback',@set_active);
set(hp_segdraw,'Callback',@set_segmode);
set(hpb_clearseg,'Callback',@clear_active);
set(hpb_clearcyst,'Callback',@clear_mme);
set([hpb_savemme, hpb_savemmeas],'Callback',@save_mmeseg);
set([hpb_saveseg, hpb_savesegas],'Callback',@save_seg);

% t = timer('TimerFcn',@(x,y) set_debug_string(''),'StartDelay',1.5);
% uiwait(hf)


%%  Utility functions for MYGUI

    % ------------------------------------------------------------------- %
    % Callback for mouse scrollwheel to change the current image
    function scrollfcn(hfig,events)
    % Callback for mouse scrollwheel
    %   hfig: handle to figure
    %   events: event structure, containing members
    %       VerticalScrollCount  - +1 for scroll down, -1 for scroll up
    %       VerticalScrollAmount - lines to scroll per count (unused here)
    % Called when scroll wheel is adjusted in figure hfig.
        cp = get(hfig,'CurrentPoint');  % Where mouse is, in figure units
        obj_pos = get(ha,'Position');   % Axes bounds

        % Is mouse pointer within axes?
        if cp(1) > obj_pos(1) && cp(1) < obj_pos(1) + obj_pos(3) && ...
           cp(2) > obj_pos(2) && cp(2) < obj_pos(2) + obj_pos(4)

            sc = events.VerticalScrollCount;
%             sa = events.VerticalScrollAmount;
            
            if ~isempty(octvol)
                % Change current B-scan image
                curimg = curimg - sc;
                
                updateimage();
                updateseg();
                updatemme();
            end
        else                      
            % Out of axis
        end    
    end
    
    % ------------------------------------------------------------------- %
    % Update the image (display a new slice)
    function updateimage
        if ~isempty(octvol)
            if curimg < 1
                curimg = 1;
            elseif curimg > size(octvol,3)
                curimg = size(octvol,3);
            end
            
            % Show image
            if isempty(h_img)
                h_img = imagesc(octvol(:,:,curimg),'Parent',ha,'HitTest','off');
                colormap gray
                axis(ha,'off')
                axis(ha,'tight')
                hold(ha,'on')
                
                % Always have mask overlaid, alpha data makes it invisible
                maskdata = zeros([size(octvol(:,:,1)),3]);
                maskdata(:,:,1) = cystvals(1);
                maskdata(:,:,2) = cystvals(2);
                maskdata(:,:,3) = cystvals(3);                
                h_mask = imagesc(maskdata,'Parent',hm,'HitTest','off');
                axis(hm,'off')
                axis(hm,'tight')
                set(h_mask,'AlphaData',cyst_mask);
            else
                set(h_img,'CData',octvol(:,:,curimg));
            end
            
            caxis(ha,[minv maxv]);
            
            % Display current bscan index
            set(hedit_sl,'String',num2str(curimg));
            set(hs_sl,'Value',curimg);
            
            if ~isempty(bsc_header)
                set(htxt_snr_val,'String',num2str(bsc_header.Quality(curimg)));
            end
            
            updatefundus();
        end        
    end

    % ------------------------------------------------------------------- %
    % Update all boundary points on an image
    function updateseg
        hideseg = (get(hcb_hideseg,'Value') == get(hcb_hideseg,'Max'));
        hideseg_rb = (get(hcb_hiderb,'Value') ~= get(hcb_hiderb,'Max'));
        
        % Remove boundaries from plot
        delete(h_pts(h_pts>0));
        h_pts(h_pts>0) = 0;        
        delete(h_pts_rb(h_pts_rb>0));
        h_pts_rb(h_pts_rb>0) = 0;
        
        % Remove control points from screen
        if ~isempty(h_ctrl)
            delete(h_ctrl)
            h_ctrl = [];
            delete(h_activept)
            h_activept = [];
        end
        
        % Add new points
        if ~hideseg && curimg >= starty && curimg+starty-1 <= size(octvol,3)
            for i = 1:size(control_pts,2)
                ctrl_pts = control_pts{curimg-starty+1,i};
                if isempty(ctrl_pts)
                    continue
                end
                
                % Draw lines
                if size(ctrl_pts,1) > 3
                    if segmode == 1
                        [~,inds] = sort(ctrl_pts(:,1));
                        ctrl_pts = ctrl_pts(inds,:);

                        pts = interpolateCtrlPts(ctrl_pts);
                    else
                        pts = ctrl_pts;
                    end
                    
                    h_pts(i) = line(pts(:,1)+startx-1,pts(:,2),'Color','r','Parent',ha,'HitTest','off');
                end
            end
        end
        if ~hideseg_rb
            xpts = 1:size(ilm_bdry_qf,1);
            ypts = ilm_bdry_qf(:,curimg);
            h_pts_rb(1) = line(xpts,ypts,'Color','m','Parent',ha,'HitTest','off');
            ypts = rpe_bdry_qf(:,curimg);
            h_pts_rb(2) = line(xpts,ypts,'Color','m','Parent',ha,'HitTest','off');
        end
        
        % Show active boundary if necessary
        updateseg_active();
        
        % Color active line as green and move to front
        if (h_pts(active_bd) > 0) && (get(hbg_mode,'SelectedObject') == htb_bdry)
            set(h_pts(active_bd),'Color','g')
            
            % Bring to front
            c = get(ha,'Children');
            c = [h_pts(active_bd); c];
            c(find(c==h_pts(active_bd),1,'last')) = [];
            set(ha,'Children',c);
        end
    end

    % ------------------------------------------------------------------- %
    % Update only the active boundary (for better click and drag performance)
    function updateseg_active
        hideseg = (get(hcb_hideseg,'Value') == get(hcb_hideseg,'Max'));
        segedit = get(hbg_mode,'SelectedObject') == htb_bdry;

        if segedit            
            ctrl_pts = control_pts{curimg,active_bd};
            active_pt_tmp = active_pt;
            if ~isempty(ctrl_pts)
                % Draw lines
                if ~hideseg && size(ctrl_pts,1) > 3
                    if segmode == 1
                        [~,inds] = sort(ctrl_pts(:,1));
                        ctrl_pts = ctrl_pts(inds,:);
                        active_pt_tmp = find(inds==active_pt); %tmp since sorted

                        pts = interpolateCtrlPts(ctrl_pts);
                    else
                        pts = ctrl_pts;
                    end

                    if h_pts(active_bd) == 0
                        h_pts(active_bd) = line(pts(:,1),pts(:,2),...
                            'Color','g','Parent',ha,'HitTest','off');
                    else
                        set(h_pts(active_bd),'XData',pts(:,1),...
                                             'YData',pts(:,2));
                    end
                end
                if segmode == 1
                    if isempty(h_ctrl)
                        % Draw control points
                        h_ctrl = line(ctrl_pts(:,1),ctrl_pts(:,2),...
                            'Color','w','Parent',ha,'Marker','s',...
                            'MarkerSize',4,'MarkerFaceColor','k',...
                            'LineStyle','none','HitTest','off','linewidth',1);
                        h_activept = line(ctrl_pts(active_pt_tmp,1),ctrl_pts(active_pt_tmp,2),...
                            'Color','y','Parent',ha,'Marker','s',...
                            'MarkerSize',4,'MarkerFaceColor','k',...
                            'HitTest','off');
                    else
                        set(h_ctrl,'XData',ctrl_pts(:,1),...
                                   'YData',ctrl_pts(:,2),...
                                   'Color','w');
                        set(h_activept,'XData',ctrl_pts(active_pt_tmp,1),...
                                       'YData',ctrl_pts(active_pt_tmp,2),...
                                       'Color','y');
                    end
                end
            end
        end
    end

    % ------------------------------------------------------------------- %
    % Callback for clicking the "Clear active boundary" button
    function clear_active(hObject, eventdata)                 %#ok<*INUSD>
        control_pts{curimg,active_bd} = [];
        updateseg();
    end

    % ------------------------------------------------------------------- %
    % Callback for clicking the "Hide segmentation" radio button
    function hide_seg(hObject, eventdata)                
        % Update display
        updateseg();
    end

    % ------------------------------------------------------------------- %
    % Update the pseudocyst segmentation on the image
    function updatemme
        hideseg = (get(hcb_hidemme,'Value') == get(hcb_hidemme,'Max'));
        
        cyst_mask(:) = false;
        cyst_px_cur = cyst_px{curimg};
        if ~isempty(cyst_px_cur)
            if get(hbg_mode,'SelectedObject') == htb_mme
                set([hpb_clearcyst,hpb_savemme,hpb_savemmeas],'Enable','on');
            end
            
            if ~hideseg
                cyst_mask(cyst_px{curimg}) = true;
            end
        else
            set([hpb_clearcyst,hpb_savemme,hpb_savemmeas],'Enable','off')
        end
        
        set(h_mask,'AlphaData',cyst_mask)
    end

    % ------------------------------------------------------------------- %
    % Update the pseudocyst segmentation on the image at a single point
    function updatemme_pt(xpt,ypt,addpt)
        
        if addpt
            % add
            cyst_mask(ypt,xpt) = true;
        else
            % remove
            cyst_mask(ypt,xpt) = false;
        end

        set(h_mask,'AlphaData',cyst_mask)
    end

    % ------------------------------------------------------------------- %
    % Callback for clicking the "Clear active pseudocysts" button
    function clear_mme(hObject, eventdata) 
        cyst_px{curimg} = [];
        updateimage();
        updatemme();
        
        % Need to return focus to the figure, but there doesn't seem to be
        % any easy way to do it, perhaps setting the keypressfcn of a dummy
        % uicontrol that mimics the keypressfcn of the figure?
    end

    % ------------------------------------------------------------------- %
    % Callback for clicking the "Hide pseudocysts" radio button
    function hide_mme(hObject, eventdata)
        % Update display
        updatemme();
    end

    % ------------------------------------------------------------------- %
    % Update the fundus image display
    function updatefundus        
        hidefun = (get(hcb_fun,'Value') ~= get(hcb_fun,'Max'));
        
        if ~hidefun && ~isempty(octvol) && ~isempty(fundus_proj)
            fi = repmat(fundus_proj,[1 1 3]);
            fi(size(octvol,3)-curimg+1,:,2:3) = 0; % fundus is flipped

            if isempty(himg_fun)
                himg_fun = imagesc(fi,'Parent',ha_fun,'HitTest','off');
                axis(ha_fun,'square')
                set(ha_fun,'xtick',[],'ytick',[])
        %         hold(ha_fun,'on')
            else
                set(himg_fun,'CData',fi);
            end
        end
    end

    % ------------------------------------------------------------------- %
    % Callback for clicking the "Show fundus/projection image" button
    function show_fundus(hObject, eventdata)
        % Simply need to reposition things
        if (get(hObject,'Value') == get(hObject,'Max'))
            set(himg_fun,'Visible','on')
            set(ha_fun,'Visible','on')
            set([ha, hm],'Position',pax_vol_fun)
            set(htxt_sl,'Position',ptxt_sl_fun)
            set(htxt_of,'Position',ptxt_of_fun)
            set(htxt_img,'Position',ptxt_img_fun)
            set(hedit_sl,'Position',pedit_sl_fun)
            set(hs_sl,'Position',ps_sl_fun)
            set(htxt_db,'Position',ptxt_db_fun)
            set(htxt_snr,'Position',ptxt_snr_fun)
            set(htxt_snr_val,'Position',ptxt_snr_val_fun)
            updatefundus();
        else
            set(himg_fun,'Visible','off') % hide image
            set(ha_fun,'Visible','off') % hide image axis
            set([ha, hm],'Position',pax_vol)
            set(htxt_sl,'Position',ptxt_sl)
            set(htxt_of,'Position',ptxt_of)
            set(htxt_img,'Position',ptxt_img)
            set(hedit_sl,'Position',pedit_sl)
            set(hs_sl,'Position',ps_sl)
            set(htxt_db,'Position',ptxt_db)
            set(htxt_snr,'Position',ptxt_snr)
            set(htxt_snr_val,'Position',ptxt_snr_val)
        end
    end

    % ------------------------------------------------------------------- %
    % Callback for manually entering the image number to display
    function change_image(hObject,eventdata)
        
        if hObject == hedit_sl
            % Text text box
            user_entry = str2double(get(hObject,'string'));
            if isnan(user_entry)
                set(hObject,'String',curimg)
            else
                curimg = user_entry;
            end
            set(hs_sl,'Value',str2double(curimg))
        elseif hObject == hs_sl
            % Slider
            curimg = round(get(hObject,'value'));
            set(hedit_sl,'String',num2str(curimg))
        end   

        updateimage();
        updateseg();
        updatemme();
    end

    % ------------------------------------------------------------------- %
    % Callback for pressing a keyboard key
    function key_press(hObject,eventdata)
        if ~isempty(octvol)
            switch eventdata.Key
                case 'leftarrow'
                    curimg = curimg - 1;
                    updateimage();
                    updateseg();
                    updatemme();
                case 'rightarrow'
                    curimg = curimg + 1;
                    updateimage();
                    updateseg();
                    updatemme();
                case 'uparrow'
                    % Move current point up
                    incr = 0.25;
                    if ~isempty(control_pts) && ~isempty(control_pts{curimg,active_bd})
                        ctrlp = control_pts{curimg,active_bd}(active_pt,:);
                        ydat = constrain_point(ctrlp(1),ctrlp(2)-incr);
                        control_pts{curimg,active_bd}(active_pt,2) = ydat;
                        updateseg();
                    end
                case 'downarrow'
                    % Move current point down
                    incr = 0.25;
                    if ~isempty(control_pts) && ~isempty(control_pts{curimg,active_bd})
                        ctrlp = control_pts{curimg,active_bd}(active_pt,:);
                        ydat = constrain_point(ctrlp(1),ctrlp(2)+incr);
                        control_pts{curimg,active_bd}(active_pt,2) = ydat;
                        updateseg();
                    end
                case 'n'
                    % Get next control point
                    ctrlp = control_pts{curimg,active_bd};
                    [~,inds] = sort(ctrlp(:,1));
                    pt = find(inds==active_pt);
                    if pt == size(inds,1)
                        pt = 0;
                    end
                    active_pt = inds(pt+1); %tmp since sorted
                    updateseg();
                case 'p'
                    % Get previous control point
                    ctrlp = control_pts{curimg,active_bd};
                    [~,inds] = sort(ctrlp(:,1));
                    pt = find(inds==active_pt);
                    if pt == 1
                        pt = size(inds,1)+1;
                    end
                    active_pt = inds(pt-1); %tmp since sorted
                    updateseg();
                case 'h'
                    if get(htb_mme,'Value') == 1
                        if (get(hcb_hidemme,'Value') == get(hcb_hidemme,'Max'))
                            set(hcb_hidemme,'Value',get(hcb_hidemme,'Min'))
                        else
                            set(hcb_hidemme,'Value',get(hcb_hidemme,'Max'))
                        end
                        updatemme();
                    elseif get(htb_bdry,'Value') == 1                    
                        if (get(hcb_hideseg,'Value') == get(hcb_hideseg,'Max'))
                            set(hcb_hideseg,'Value',get(hcb_hideseg,'Min'))
                        else
                            set(hcb_hideseg,'Value',get(hcb_hideseg,'Max'))
                        end
                        updateseg();
                    end
            end
        end
    end

    % ------------------------------------------------------------------- %
    % Callback for clicking one of the "Display mode" toggle buttons
    function change_mode(hObject, eventdata)
        set(eventdata.OldValue,'BackgroundColor',get(hf,'Color'))
        set(eventdata.NewValue,'BackgroundColor','w')
        
        % Default cyst overlay color
        cystvals = [1 1 0];
        
        switch eventdata.NewValue
            case htb_disp
                % Update to remove active boundary
                updateseg();
                updatemme();
                
                % No callback
                set(hf,'WindowButtonDownFcn','');
                
                set([hpb_clearcyst, hpb_clearseg, hp_active,...
                     hpb_savemme, hpb_savemmeas, hpb_saveseg,...
                     hpb_savesegas, hp_segdraw],'Enable','off')
                 
            case htb_bdry
                % Enable manual layer segmentation
                
                % Update to show/hide active boundary
                updateseg();
                updatemme();
                
                % Edit callback
                set(hf,'WindowButtonDownFcn',@segclick)
                
                set([hpb_clearcyst, hpb_savemme, hpb_savemmeas],'Enable','off')
                set([hpb_clearseg, hp_active, hp_segdraw, hpb_saveseg,...
                     hpb_savesegas],'Enable','on')
                
            case htb_mme
                % Enable manual pseudocyst segmentation
                
                % Change cyst overlay color
                cystvals = [0 1 0];
                
                % Update to remove active boundary
                updateseg();
                updatemme();
                
                % Edit callback
                set(hf,'WindowButtonDownFcn',@mmeclick);
                
                set([hpb_clearseg, hp_active, hp_segdraw, hpb_saveseg,...
                     hpb_savesegas],'Enable','off')
        end
        
        maskdata = zeros([size(octvol(:,:,1)),3]);
        maskdata(:,:,1) = cystvals(1);
        maskdata(:,:,2) = cystvals(2);
        maskdata(:,:,3) = cystvals(3);

        set(h_mask,'CData',maskdata)        
    end
    
    % ------------------------------------------------------------------- %
    % Load OCT volume
    function load_vol(hObject, eventdata)
        tag = get(hObject,'Tag');
        pathname = '.';
        switch tag
            case 'vol'
                fileext = '*.vol';
                if ispref('octviewer','vol_path')
                    pathname = getpref('octviewer','vol_path');
                else
                    addpref('octviewer','vol_path','.');
                end
            case 'img'
                fileext = '*Macular Cube*_z.img';
                fileext = '*.img';
                if ispref('octviewer','img_path')
                    pathname = getpref('octviewer','img_path');
                else
                    addpref('octviewer','img_path','.');
                end
            case 'mat'
                fileext = '*.mat';
                if ispref('octviewer','mat_vol_path')
                    pathname = getpref('octviewer','mat_vol_path');
                else
                    addpref('octviewer','mat_vol_path','.');
                end
            case 'oct'
                fileext = '*.oct';
                if ispref('octviewer','oct_vol_path')
                    pathname = getpref('octviewer','oct_vol_path');
                else
                    addpref('octviewer','oct_vol_path','.');
                end
        end
        [filename_vol, pathname] = uigetfile(fileext,'Load OCT data',pathname);
        
        if isequal(filename_vol,0)
            filename_vol = [];
            return
        end
        
        set(hf,'Name',figure_name);
        drawnow;
        
        delete(h_pts(h_pts>0));
        h_pts(h_pts>0) = 0;
        delete(h_pts_rb(h_pts_rb>0));
        h_pts_rb(h_pts_rb>0) = 0;
        control_pts = [];
        bsc_header = [];
        himg_fun = [];
        
        filefull = fullfile(pathname, filename_vol);
        
        set_debug_string('Loading volume...')
        
        switch tag
            case 'vol'
                setpref('octviewer','vol_path',pathname)
                
                [header, octvol, slo, bsc_header] = octReader(filefull);
                
%                 [header, bsc_header, slo, octvol] = openVolFast(filefull,'nodisp');
%                 octvol = double(octvol).^0.25;
%                 octvol(octvol>1) = 0;
                
                voxel_size = [header.ScaleZ header.ScaleX header.Distance];
            case 'img'
                setpref('octviewer','img_path',pathname)
                
                [header, octvol, slo, bsc_header] = octReader(filefull);
                voxel_size = [header.ScaleZ header.ScaleX header.Distance];
                
%                 [octvol,vol_info] = octCirrusReader(filefull);
%                 voxel_size = vol_info.vol_res;
            case 'mat'
                setpref('octviewer','mat_vol_path',pathname)
                
                S = load(filefull);
                octvol = double(S.img_vol).^0.25;
                voxel_size = [S.header.ScaleZ S.header.ScaleX S.header.Distance];
            case 'oct'
                setpref('octviewer','oct_vol_path',pathname)
                
                [octvol,b_header] = octBioptigenReader(filefull,false,true);
                
                b_header.lineLength = round(b_header.lineLength/2);
                b_header.scanDepth = b_header.scanDepth/2;
                octvol = octvol(1:b_header.lineLength,:,:);
                
                header.ScaleX = b_header.azScanLength/b_header.lineCount / 1.4 * 6;
                header.ScaleZ = b_header.scanDepth/b_header.lineLength / 0.8 * 2 / 2;
                header.Distance = b_header.elScanLength/b_header.frameCount / 1.4 * 6;
                
                header. SizeX = b_header.lineCount;
                header.SizeZ = b_header.lineLength;
                header.NumBScans = b_header.frameCount;
                
%                 % Give values that 
%                 header.ScaleX = 6/size(octvol,1);
%                 header.ScaleZ = 0.5/size(octvol,2);
%                 header.Distance = 6/size(octvol,3);               
                
                voxel_size = [header.ScaleZ header.ScaleX header.Distance];
        end
        
        % Return focus to figure
        figure(hf);
        
        initialize_vol();
    end

    function initialize_vol
        % Set image data
        octvol = double(octvol);
        octvol = (octvol - min(octvol(:)))/(max(octvol(:))-min(octvol(:)));
        
        % Generate fundus projection image
        set_debug_string('Generating fundus projection image...')
        
        fundus_proj = [];
        startx = 1;
        starty = 1;
        if size(octvol,3) > 10 && ~isempty(voxel_size)
%             try
                [rpe_bdry_qf, ilm_bdry_qf] = quickFindRPE(octvol,voxel_size);
                fundus_proj = octFundus(octvol,rpe_bdry_qf);
                fundus_proj = (fundus_proj - min(fundus_proj(:)))/(max(fundus_proj(:))-min(fundus_proj(:)));

                % Segmentation gets centered on the data
                startx = round(size(octvol,2)/2) - floor((size(rpe_bdry_qf,1)-1)/2);
                starty = round(size(octvol,3)/2) - floor((size(rpe_bdry_qf,2)-1)/2);
%             end
        end
        
        control_pts = cell(size(octvol,3),11);
        cyst_px = cell(size(octvol,3),1);
        cyst_mask = false(size(octvol(:,:,1)));
        
        % Set current image index
        curimg = round(size(octvol,3)/2);
        set(htxt_sl,'String',num2str(size(octvol,3)));
        
        set(hs_sl,'Max',size(octvol,3),...
                  'SliderStep',[1/size(octvol,3) 5/size(octvol,3)],...
                  'Value',curimg,...
                  'Enable','on');
        
        if ~isempty(h_img)
            delete(h_img)
            h_img = [];
        end
        
        if ~isempty(h_mask)
            delete(h_mask)
            h_mask = [];
        end
        
        if~isempty(bsc_header)
            set(htxt_snr_val,'String',num2str(bsc_header.Quality(curimg)));
        end
        
        set_debug_string('Updating...')
        updateimage();
        updateseg();
        
        % Reset variables
        savefile_mme = [];
        filename_mme = [];
        savefile_seg = [];
        filename_seg = [];
        segmode = 1;
        
        % Enable buttons
        set([hcb_hideseg hcb_hiderb],'Enable','on')
        set(hcb_hidemme,'Enable','on')
        set(htb_bdry,'Enable','on')
        set(htb_mme,'Enable','on')
        set(hp_segdraw,'Value',1)
        
        set(hf,'Name',[figure_name ' - ' filename_vol]);
        set([hms,hmm,hmip],'Enable','on');
        set(hcb_fun,'Enable','on');
        set_debug_string('');
    end
 
    % ------------------------------------------------------------------- %
    % Load layer segmentation file
    function load_seg(hObject, eventdata)
        tag = get(hObject,'Tag');
        pathname = '.';
        switch tag
            case 'otb'
                fileext = '*.txt';
                if ispref('octviewer','seg_path_otb')
                    pathname = getpref('octviewer','seg_path_otb');
                else
                    addpref('octviewer','seg_path_otb','.');
                end
            case 'mat'
                fileext = '*.mat';
                if ispref('octviewer','seg_path_mat')
                    pathname = getpref('octviewer','seg_path_mat');
                else
                    addpref('octviewer','seg_path_mat','.');
                end
        end
        [filename, pathname] = uigetfile(fileext,'Load segmentation file',pathname);
        % open the directory of the folder
        
        if isequal(filename,0)
            return
        end
        filename_seg = fullfile(pathname, filename);
        
        % Return focus to figure
        figure(hf);
        
        switch tag
            case 'otb'
                setpref('octviewer','seg_path_otb',pathname)
                
                control_pts = readControlPoints(filename_seg,[],false);
                
                for i = 1:size(control_pts,1)
                    for j = 1:size(control_pts,2)
                        % OCT toolbox indexing starts at 0
                        control_pts{i,j} = control_pts{i,j} + 1;
                    end
                end
            case 'mat'
                setpref('octviewer','seg_path_mat',pathname)
                
                S = load(filename_seg);
                
                if isfield(S,'bd_pts')
                    bd_pts = S.bd_pts;
                    
                    control_pts = cell(size(bd_pts,2),11);
                    inds = [1 2 4 5 6 7 8 9 10];
                    for i = 1:size(bd_pts,2)
                        for j = 1:size(bd_pts,3)
                            control_pts{i,inds(j)} = [1:size(bd_pts,1); bd_pts(:,i,j)']';
                        end
                        control_pts{i,3} = [];
                        control_pts{i,11} = [];
                    end
                    
                    % Segmentation gets centered on the data
                    startx = round(size(octvol,2)/2) - floor((size(bd_pts,1)-1)/2);
                    starty = round(size(octvol,3)/2) - floor((size(bd_pts,2)-1)/2);
                    
                    set(hp_segdraw,'Value',2);
                    segmode = 2;
                elseif isfield(S,'control_pts')
                    control_pts = S.control_pts;
                end
                if isfield(S,'segmode')
                    segmode = S.segmode;
                else
                    segmode = 1;
                end
                set(hp_segdraw,'Value',segmode);                
        end
        
        savefile_seg = [];
        active_bd = get(hp_active,'Value');
        active_pt = 1;
        updateseg();
    end

    % ------------------------------------------------------------------- %
    % Load layer segmentation manual result to file
    function save_seg(hObject, eventdata)
        tag = get(hObject,'Tag');
        
        pathname = '.';
        if ispref('octviewer','seg_path_seg')
            pathname = getpref('octviewer','seg_path_seg');
        else
            addpref('octviewer','seg_path_seg','.');
        end
        
        % If no previously saved file or user clicked save as...
        if isempty(savefile_seg) || strcmp(tag,'saveas')
            % If previously saved file
            if ~isempty(savefile_seg)
                % Same filename as the previous saved file - let the user
                % choose to change it
                [filename, pathname] = uiputfile('*.mat','Save manual segmentation file as...',savefile_seg);
                
            elseif ~isempty(filename_seg)
                % Same file name as the one we loaded
                [filename, pathname] = uiputfile('*.mat','Save manual segmentation file as...',filename_seg);
                
            elseif ~isempty(filename_vol)                
                % Same filename as the volume that was loaded
                [~, name, ~] = fileparts(filename_vol); % remove extension
                [filename, pathname] = uiputfile('*.mat','Save manual segmentation file as...',fullfile(pathname,[name '_seg_manual.mat']));
                
            else
                % Make something up
                [filename, pathname] = uiputfile('*.mat','Save manual segmentation file as...',fullfile(pathname,'manual_segmentation.mat'));
                
            end
            
            if isequal(filename,0)
               return
            end        
            savefile_seg = fullfile(pathname, filename);
        end
        setpref('octviewer','seg_path_seg',pathname)
        
        % Return focus to figure
        figure(hf);
        
        set_debug_string('Saving...')
        save(savefile_seg,'control_pts','segmode')
        set_debug_string('File saved')
%         start(t);        
    end

    % ------------------------------------------------------------------- %
    % Load pseudocyst segmentation file
    function load_mmeseg(hObject, eventdata)        
        pathname = '.';
        if ispref('octviewer','seg_path_mme')
            pathname = getpref('octviewer','seg_path_mme');
        else
            addpref('octviewer','seg_path_mme','.');
        end
        
        [filename, pathname] = uigetfile('*.mat','Load pseudocyst segmentation file',pathname);
        
        if isequal(filename,0)
            return
        end
        filename_mme = fullfile(pathname, filename);
        
        % Return focus to figure
        figure(hf);
        
        S = load(filename_mme);
                
        if isfield(S,'cyst_px')
            cyst_px = S.cyst_px;
        elseif isfield(S,'mask_vol')
            for i = 1:size(S.mask_vol,3)
                mv = S.mask_vol(:,:,i);
                cyst_px{i} = find(mv);
            end
        else
            return            
        end
        
        setpref('octviewer','seg_path_mme',pathname)
        updatemme();
    end

    % ------------------------------------------------------------------- %
    % Save pseudocyst manual segmentation result to file
    function save_mmeseg(hObject, eventdata)
        tag = get(hObject,'Tag');
        
        pathname = '.';
        if ispref('octviewer','seg_path_mme')
            pathname = getpref('octviewer','seg_path_mme');
        else
            addpref('octviewer','seg_path_mme','.');
        end
        
        % If no previously saved file or user clicked save as...
        if isempty(savefile_mme) || strcmp(tag,'saveas')
            % If previously saved file
            if ~isempty(savefile_mme)
                % Same filename as the previous saved file - let the user
                % choose to change it
                [filename, pathname] = uiputfile('*.mat','Save pseudocyst segmentation file as...',savefile_mme);
                
            elseif ~isempty(filename_mme)
                % Same file name as the one we loaded
                [filename, pathname] = uiputfile('*.mat','Save pseudocyst segmentation file as...',filename_mme);
                
            elseif ~isempty(filename_vol)                
                % Same filename as the volume that was loaded
                [~, name, ~] = fileparts(filename_vol); % remove extension
                [filename, pathname] = uiputfile('*.mat','Save pseudocyst segmentation file as...',fullfile(pathname,[name '_mme_seg_manual.mat']));
                
            else
                % Make something up
                [filename, pathname] = uiputfile('*.mat','Save pseudocyst segmentation file as...',fullfile(pathname,'pseudocyst_manual_segmentation.mat'));
                
            end
            
            if isequal(filename,0)
               return
            end        
            savefile_mme = fullfile(pathname, filename);
        end
        setpref('octviewer','seg_path_mme',pathname)
        
        % Return focus to figure
        figure(hf);
        
        set_debug_string('Saving...')
        save(savefile_mme,'cyst_px')
        set_debug_string('File saved')
%         start(t);
    end

    % ------------------------------------------------------------------- %
    % Set active boundary for segmentation
    function set_active(hObject, eventdata)
        val = get(hObject,'Value');
        active_bd = val;
        active_pt = 1;
        
        updateseg();
    end

    function set_segmode(hObject, eventdata)
        segmode_prev = segmode;
        segmode = get(hObject,'Value');
        
        if segmode == segmode_prev
            return
        end
        
        if segmode == 1
            % Freehand -> Spline
            choice = questdlg(['Setting the segmentation mode to '...
                               'control points will add control points '...
                               'by a least squares fit changing the '...
                               'freehand result (not reversable), continue?'],...
                               'Add control points','Yes');
                           
            if strcmp(choice,'Yes')
                % # of control points to use
%                 nk = round(size(octvol,2)/75);
                nk = 20;
                
                for i = 1:size(octvol,3)
                    for j = 1:11
                        ctrlp = control_pts{i,j};
                        if ~isempty(ctrlp)
                            xpts = ctrlp(:,1);
                            ypts = ctrlp(:,2);
                            
                            % Fit a cubic spline
                            sp = spap2(nk, 4, xpts, ypts);
                            % Adjust knot points to have a better fit
                            kp = round(newknt(sp));
                            % Take the interior points (see fn2fm)
                            kp = kp(4:end-3);
                            
                            control_pts{i,j} = [kp' ypts(kp)];
                        end
                    end
                end
                active_pt = 1;
            else
                set(hObject,'Value',2)
            end 
            
        elseif segmode == 2
            % Spline -> Freehand
            choice = questdlg(['Setting the segmentation mode to '...
                               'freehand will remove all control points '...
                               '(not reversable), continue?'],...
                               'Remove control points','Yes');
                           
            if strcmp(choice,'Yes')
                % Convert spline control points to full boundaries
                pts0 = nan(size(octvol,2),2);
                pts0(:,1) = 1:size(octvol,2);
                for i = 1:size(octvol,3)
                    for j = 1:11
                        ctrlp = control_pts{i,j};
                        if ~isempty(ctrlp)
                            [~,inds] = sort(ctrlp(:,1));
                            ctrlp = ctrlp(inds,:);
                            pts = pts0;
                            
                            % Interpolate spline at every A-scan
                            pts(ctrlp(1,1):ctrlp(end,1),:) = interpolateCtrlPts(ctrlp);
                            % Extend to edge of image
                            pts(1:ctrlp(1,1),2) = ctrlp(1,2);
                            pts(ctrlp(end,1):end,2) = ctrlp(end,2);
                            
                            control_pts{i,j} = pts;
                        end                        
                    end
                end                
            else
                set(hObject,'Value',1)
            end            
        end
        updateseg();
        
    end

    % ------------------------------------------------------------------- %
    % Callback for adding and deleting points by left and right clicking 
    function segclick(src,eventdata)
        cp = get(src,'CurrentPoint');  % Where mouse is, in figure units
        obj_pos = get(ha,'Position');   % Axes bounds

        % Is mouse pointer within axes?
        if cp(1) > obj_pos(1) && cp(1) < obj_pos(1) + obj_pos(3) && ...
           cp(2) > obj_pos(2) && cp(2) < obj_pos(2) + obj_pos(4)
        
            pt = get(ha,'CurrentPoint');
            xdat = round(pt(1,1)); ydat = pt(1,2);
            
            ctrlp = control_pts{curimg,active_bd};
            
            if isempty(ctrlp) && segmode == 2
                % No freehand without points already there
                return
            end

            if strcmp(get(src,'SelectionType'),'normal')
                % Left click - add or move point
                
                % Constrain point to be between boundaries
                ydat = constrain_point(xdat,ydat);
                
                if segmode == 1
                    if isempty(ctrlp)
                        ctrlp = [xdat ydat];
                        active_pt = 1;
                    else
                        % Check for a close point
                        [ind,dist] = getClosestPoint([xdat ydat],ctrlp);

                        if dist < size(octvol,2)/128
                            % Move the close point
                            ctrlp(end+1,:) = ctrlp(ind,:);
                            ctrlp(ind,:) = [];
                        else
                            % Move the new point
                            ctrlp(end+1,:) = [xdat ydat];
                        end
                        active_pt = size(ctrlp,1);
                    end
                else
                    ctrlp(xdat,2) = ydat;
                end
                
                control_pts{curimg,active_bd} = ctrlp;
                updateseg_active();
                
                set(src,'WindowButtonMotionFcn',@wbmcb)
                set(src,'WindowButtonUpFcn',@wbucb)                
            elseif strcmp(get(src,'SelectionType'),'alt')
                % Right click
                
                if segmode == 1
                    % Check for a close point to delete                    
                    [ind,dist] = getClosestPoint([xdat ydat],ctrlp);

                    if dist < 8
                        ctrlp(ind,:) = [];
                        control_pts{curimg,active_bd} = ctrlp;
                        updateseg_active();
                    end
                else
                    ctrlp(xdat,2) = nan;
                    control_pts{curimg,active_bd} = ctrlp;
                    updateseg_active();
                end
            end
        end
        
            function wbmcb(src,evnt)
                pt = get(ha,'CurrentPoint');
                xdat2 = round(pt(1,1));
                ydat2 = pt(1,2);
                
                % Boundary check
                if xdat2 < 1, xdat2 = 1; end
                if ydat2 < 1, ydat2 = 1; end
                if xdat2 > size(octvol,2), xdat2 = size(octvol,2); end
                if ydat2 > size(octvol,1), ydat2 = size(octvol,1); end
                
                % Constrain point to be between boundaries
                ydat2 = constrain_point(xdat2,ydat2);
                
                % Draw point
                ctrlp = control_pts{curimg,active_bd};
                if segmode == 1
                    ctrlp(end,:) = [xdat2 ydat2];
                else
                    ctrlp(xdat2,2) = ydat2;
                end
                control_pts{curimg,active_bd} = ctrlp;
                updateseg_active();
            end

            function wbucb(src,evnt)
                set(src,'WindowButtonMotionFcn','')
                set(src,'WindowButtonUpFcn','')
            end
    end

    function np = constrain_point(xpt,ypt)
        yval_prev = 1;
        for i = (active_bd - 1):-1:1
            if h_pts(i) > 0
                x_prev = get(h_pts(i),'xdata');
                xind = x_prev == xpt;
                if any(xind)
                    y_prev = get(h_pts(i),'ydata');
                    yval_prev = y_prev(x_prev == xpt);
                    break
                end
            end
        end
        
        yval_next = size(octvol,1);
        for i = (active_bd + 1):length(h_pts)
            if h_pts(i) > 0
                x_next = get(h_pts(i),'xdata');
                xind = x_next == xpt;
                if any(xind)
                    y_next = get(h_pts(i),'ydata');
                    yval_next = y_next(xind);
                    break
                end
            end
        end

        if ypt < yval_prev
            np = yval_prev;
        elseif ypt > yval_next
            np = yval_next;
        else
            np = ypt;
        end
%         np = [yval_prev yval_next];
    end

    % ------------------------------------------------------------------- %
    % Callback for adding and deleting points by left and right clicking 
    function mmeclick(src,eventdata)
        cp = get(src,'CurrentPoint');  % Where mouse is, in figure units
        obj_pos = get(ha,'Position');   % Axes bounds
        
        % Is mouse pointer within axes?
        if cp(1) > obj_pos(1) && cp(1) < obj_pos(1) + obj_pos(3) && ...
           cp(2) > obj_pos(2) && cp(2) < obj_pos(2) + obj_pos(4)
        
            if strcmp(get(src,'SelectionType'),'normal')
                wbmcb(src,eventdata,true);
                if length(cyst_px{curimg}) == 1
                    set([hpb_clearcyst,hpb_savemme,hpb_savemmeas],'Enable','on')
                end
                
                set(src,'WindowButtonMotionFcn',{@wbmcb,true})
                set(src,'WindowButtonUpFcn',@wbucb)
            elseif strcmp(get(src,'SelectionType'),'alt')
                wbmcb(src,eventdata,false);
                
                set(src,'WindowButtonMotionFcn',{@wbmcb,false})
                set(src,'WindowButtonUpFcn',@wbucb)
            end
        end
        
            function wbmcb(src,eventdata,addme)
                pt = get(ha,'CurrentPoint');
                xdat = round(pt(1,1)); ydat = round(pt(1,2));
                
                % Boundary check
                if xdat < 1, xdat = 1; end
                if ydat < 1, ydat = 1; end
                if xdat > size(octvol,2), xdat = size(octvol,2); end
                if ydat > size(octvol,1), ydat = size(octvol,1); end

                pt = size(octvol,1)*(xdat-1)+ydat;
                
                px_cur = cyst_px{curimg};
                if addme
                    % Add point                    
                    if ~any(px_cur == pt)
                        px_cur(end+1) = pt;
                    end                    
                else
                    % Remove point
                    px_cur(px_cur==pt) = [];
                end
                cyst_px{curimg} = px_cur;
                
                updatemme_pt(xdat,ydat,addme);
            end

            function wbucb(src,eventdata)
                set(src,'WindowButtonMotionFcn','')
                set(src,'WindowButtonUpFcn','')
            end
    end

    % ------------------------------------------------------------------- %
    % Intensity normalization of the data
    function normalize_data(hObject, eventdata)
        set_debug_string('Intensity normalizing data')
        tag = get(hObject,'Tag');
        
        header.ScaleZ = voxel_size(1);
        header.ScaleX = voxel_size(2);
        header.Distance = voxel_size(3);
        bds = cat(3,ilm_bdry_qf,rpe_bdry_qf);
        switch tag
            case 'ac'
                % Attenuation correction - Girard et al 2011
                octvol = normalizeOCTVolume(octvol,4,header);
                octvol = (octvol - min(octvol(:)))/(max(octvol(:)) - min(octvol(:)));
            case 'rcs'
                % Median filter contrast stretching
                octvol = normalizeOCTVolume(octvol,2,header);
            case 'rr_eq'
                % Correction to make RPE and RNFL similar intensities
                octvol = normalizeOCTVolume(octvol,8,header,bds);
            case 'r_eq'
                % Correction to make RPE a value of 1 and bg a value of 0
                octvol = normalizeOCTVolume(octvol,7,header,bds);
            case 'r_eq2'
                % Only the smooth RPE intensity correction
                octvol = normalizeOCTVolume(octvol,9,header,bds);
            case 'r_eq3'
                % Full RPE intensity correction
                octvol = normalizeOCTVolume(octvol,6,header,bds);
            case 'r_eq4'
                % B-scan bias field RPE intensity correction only
                octvol = normalizeOCTVolume(octvol,10,header,bds);
        end
        
        if ~strcmp(tag,'ac')
            % Regenerate fundus image
%             rpe_bdry = quickFindRPE(octvol,voxel_size);
            fundus_proj = octFundus(octvol,rpe_bdry_qf);
            fundus_proj = (fundus_proj - min(fundus_proj(:)))/(max(fundus_proj(:))-min(fundus_proj(:)));
        end
        
        updateimage();
        set_debug_string('')
    end

    % ------------------------------------------------------------------- %
    % Median filter the data
    function median_filter_data(hObject, eventdata)
        sz_vol = size(octvol);
        dn_k = [3 3 1];
        dn_k = [7 1 1];
        octvol = permute(octvol,[2 1 3]);
        set_debug_string('Median filtering data')
        octvol = medfilt2(octvol(:,:),[dn_k(2) dn_k(1)],'symmetric');
        octvol = reshape(octvol,sz_vol(2),sz_vol(1),sz_vol(3));
        octvol = permute(octvol,[2 1 3]);
        set_debug_string('');
        
        updateimage();
    end

    % ------------------------------------------------------------------- %
    % MME Filter on current image - median filter followed by bilateral
    % filter
    function MME_filter_data(hObject, eventdata)
        
        tag = get(hObject,'Tag');
        
        if strcmp(tag,'s')
            inds = curimg;
        elseif strcmp(tag,'v')
            inds = 1:size(octvol,3);
        end
        
%         ks = 2*round([7.82/(voxel_size(1)*1000) 11.74/(voxel_size(2)*1000)]) + 1;
        ks = 2*round([8/(voxel_size(1)*1000) 8/(voxel_size(2)*1000)]) + 1;

        sigmaSpatial = round(125/(voxel_size(1)*1000));
        sigmaRange = 1/10;

        set_debug_string('Filtering data')
        
        for i = 1:length(inds)
            if length(inds) > 1
                set_debug_string(['Filtering data - ' num2str(i) '/' num2str(length(inds))])
            end
            
            img_mf = medfilt2(octvol(:,:,inds(i)),ks);
            img_mf = bilateralFilter(double(img_mf), [], 0, 1, sigmaSpatial, sigmaRange);
            img_mf = (img_mf-min(img_mf(:)))/(max(img_mf(:))-min(img_mf(:)));

            octvol(:,:,inds(i)) = img_mf;
        end
        set_debug_string('');        
        
        updateimage();
    end

    function Experimental_filter(hObject, eventdata)
        [rpe_bdry,~,shifts] = quickFindRPE(octvol,voxel_size);
        octvol2 = retinaFlatten(octvol,shifts,'linear');
        
        octvol2 = permute(octvol2,[3 1 2]);
        vs = size(octvol2);
        octvol2 = medfilt2(octvol2(:,:),[5 1]);
        octvol2 = reshape(octvol2, vs);
        octvol2 = permute(octvol2,[2 3 1]);
        
        octvol2 = retinaFlatten(octvol2,-shifts,'linear');
        
        octvol = (octvol + octvol2)/2;
        
        updateimage();
    end

    % ------------------------------------------------------------------- %
    % Allow the user to drag a rectangle over the image and change the
    % contrast such that the min and max intensity in the rectangle become
    % the min and max in the entire image
    function contrast_window(hObject, eventdata)
        pts = round(getrect(ha));
        
        if all(pts>0)
            img = octvol(pts(2):(pts(2)+pts(4)),pts(1):(pts(1)+pts(3)),curimg);
            maxv = max(img(:));
            minv = min(img(:));
        end
        
%         % Stretch intensities to be between 0 and 1
%         octvol = (octvol - minv)/(maxv - minv);
%         octvol(octvol>1) = 1;
%         octvol(octvol<0) = 0;
        
        updateimage();
    end

    % Reset image contrast
    function reset_contrast(hObject, eventdata)
        maxv = 1;
        minv = 0;
        
        updateimage();
    end

    % Set image contrast
    function contrast_range(hObject, eventdata)
        range = inputdlg({'Min','Max'},'Contrast range',2,{num2str(minv),num2str(maxv)});
        
        if ~isempty(range)
            range = str2double(range);
            if range(2) > range(1)
                maxv = range(2);
                minv = range(1);
            end
        end
        
        updateimage();
    end

    % ------------------------------------------------------------------- %
    % Close request function for the figure
    function figure_close(hObject, eventdata)
        selection = questdlg('Close This Figure?',...
            'Close Request Function',...
            'Yes','No','Yes'); 
        switch selection, 
            case 'Yes',
%                 delete(t)
                delete(gcf)
            case 'No'
                return
        end
    end

    function set_debug_string(str)
        set(htxt_db,'String',str,'ForegroundColor','r','FontWeight','bold')
        drawnow
    end
end

%%
function pts = interpolateCtrlPts(ctrl_pts)
% Interpolate a vector of control points to get values between the first and
% last point
    
    % Assume sorted
    xmin = ctrl_pts(1,1);
    xmax = ctrl_pts(end,1);
    
    xpts = (xmin:xmax)';
    
    % remove duplicate points 
    [pt1 m n] = unique(ctrl_pts(:,1));
    pt2 = ctrl_pts(m,2);
    ctrl_pts = [pt1, pt2];
    
%     pts = interp1(ctrl_pts(:,1),ctrl_pts(:,2),xpts,'cubic',nan);
    pts = cubicSplineInterp(ctrl_pts(:,1),ctrl_pts(:,2),xpts);
    pts = cat(2,xpts,pts);
end

%%
function [ind,dist] = getClosestPoint(pt,ctrlPts)
% Find the closest point in an array of control points to the current point

    d = sum(bsxfun(@minus,ctrlPts,pt).^2,2);

    [~,ind] = min(d);
    dist = sqrt(d(ind));

end

%% 
function splPts = cubicSplineInterp(xpts,ypts,xptsOut)

    xpts = xpts-xpts(1)+1;
    n = length(xpts);

    % X = calculateNaturalCubicSpline(n-1,xpts);
    Y = calculateNaturalCubicSpline(n-1,ypts);

    splPts = zeros(size(xptsOut));
    for i = 1:n-1
        splPts(xpts(i)) = ypts(i);

        t = (((xpts(i)+1):(xpts(i+1)-1))-xpts(i))/(xpts(i+1)-xpts(i));
        splPts((xpts(i)+1):(xpts(i+1)-1)) = Y(i,1) + Y(i,2)*t + Y(i,3)*t.*t + Y(i,4)*t.*t.*t;

    %     for j = (xpts(i)+1):(xpts(i+1)-1)
    %         t = (j - xpts(i))/(xpts(i+1)-xpts(i));
    %         y = getSplineValue(Y(i,:),t);
    %         splPts(j) = y;
    %     end    
    end
    splPts(xpts(n)) = ypts(n);

end

function C = calculateNaturalCubicSpline(n,x)

    gamma = zeros(n+1,1);
    delta = zeros(n+1,1);
    D = zeros(n+1,1);

    gamma(1) = 0.5;
    for i = 2:n
        gamma(i) = 1/(4-gamma(i-1));
    end
    gamma(end) = 1/(2-gamma(end-1));

    delta(1) = 3*(x(2)-x(1))*gamma(1);
    for i = 2:n
        delta(i) = (3*(x(i+1)-x(i-1))-delta(i-1))*gamma(i); 
    end
    delta(end) = (3*(x(end)-x(end-1))-delta(end-1))*gamma(end);

    D(end) = delta(end);
    for i = n:-1:1
        D(i) = delta(i) - gamma(i)*D(i+1);
    end

    C = zeros(n,4);
    for i = 1:n
        C(i,:) = [x(i),D(i),3*(x(i+1)-x(i))-2*D(i)-D(i+1),2*(x(i)-x(i+1))+D(i)+D(i+1)];
    end

end

function v = getSplineValue(c,x)

    v = c(1) + c(2)*x + c(3)*x^2 + c(4)*x^3;

end