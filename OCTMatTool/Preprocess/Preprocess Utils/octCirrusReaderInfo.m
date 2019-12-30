function vol_info = octCirrusReaderInfo(filename,cube_size,skipsize)
% hidef scan consists of 2 orthogonal bscans intersecting in the center of
% the volume. The first is in the y direction (across B-scans), the second
% is in the x direction (aligned with the center B-scan). Unfortunately,
% they are not in good alignment with the volumetric data

% Andrew Lang
% $Id: octCirrusReaderInfo.m,v 1.3 2015/03/17 19:50:35 andrew Exp $

if nargin < 2    
    cube_size = [];
    skipsize = false;
elseif nargin < 3    
    skipsize = false;
end

if isempty(cube_size)
    % assume 6mm x 6mm x 2mm cube size
    cube_size = [6 6 2];
end

% Parse file name for useful information
[~, name, ext] = fileparts(filename);
if ~strcmp(ext,'.img')
    error(['Cannot open files with extension ' ext '.'])
end

if ~isempty(strfind(name,'clinic'))
    C = textscan(name,'%s %s %s %s %s %s %s %s %s %s',...
             'Delimiter','_');
    C{1}{1} = [C{1}{1} '_' C{2}{1}];
    C(2) = [];
else
    C = textscan(name,'%s %s %s %s %s %s %s %s %s',...
             'Delimiter','_');
end

C2 = textscan(C{2}{1},'%s %s %s %s');
if strcmp(C2{1}{1},'Optic')
    C2{1}{1} = sprintf('%s %s',C2{1}{1},C2{2}{1});
    C2(2) = [];
end
                  

% patient ID
vol_info.pid = C{1}{1}; 
% Macular or Optic Disc
vol_info.scan_type = C2{1}{1}; 
% Volume size - '512x128' or '200x200' probably
vol_size = C2{3};
% Scan date - 'm-dd-yyyy'
vol_info.scan_date = datenum([C{3}{1} ' ' C{4}{1}],'mm-dd-yyyy HH-MM-SS');
% Eye side - 'OD' or 'OS' for right and left eye
vol_info.eye_side = C{5}{1};
% % Other not useful stuff     
% sn = C{6}; % 'snXXXX'
cube = C{7}; % 'cube'
% z = C{8}; % 'z' or 'raw'

if ~strcmp(cube,'cube') && ~strcmp(cube,'hidef')
    error(['Cannot read Cirrus file ' filename...
           '. Filename must end in cube_z or cube_raw.'])
end

if ~strcmp(vol_info.scan_type,'Macular') && ...
   ~strcmp(vol_info.scan_type,'Optic Disc')

    error(['Cannot read Cirrus file ' filename... 
           '. Must be Macular Cube or Optic Disc scan.'])
end

% Get volume dimensions
vol_size = regexp(vol_size,'[x]','split');
vol_size = str2double(vol_size{1});
num_ascan = vol_size(1);
num_bscan = vol_size(2);

if strcmp(cube,'hidef')
    num_bscan = 2;
    if strcmp(vol_info.scan_type,'Macular')
        num_ascan = 1024;
    elseif strcmp(vol_info.scan_type,'Optic Disc')
        num_ascan = 1000;
    end
end

if ~skipsize
    % Size of data
    fileinfo = dir(filename);
    if ~isempty(fileinfo)
        n = fileinfo.bytes;
        num_axial = n/(num_ascan*num_bscan);
    else
        num_axial = nan;
    end

    % Resolution in each direction - assume cube begins and ends at voxel
    %   center
    vol_info.vol_res = [cube_size(3)/(num_axial-1),cube_size(1)/(num_ascan-1),cube_size(2)/(num_bscan-1)];
    vol_info.vol_size = [num_axial,num_ascan,num_bscan];
end
