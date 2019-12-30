function [vol_data,vol_info,lso] = octCirrusReader(filename,cube_size)
% hidef scan consists of 2 orthogonal bscans intersecting in the center of
% the volume. The first is in the y direction (across B-scans), the second
% is in the x direction (aligned with the center B-scan). Unfortunately,
% they are not in good alignment with the volumetric data

% Added support for clinic scans

% Andrew Lang
% $Id: octCirrusReader.m,v 1.4 2015/03/17 19:53:15 andrew Exp $

if nargin < 2 || isempty(cube_size)
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Added by Yufan May/24/2017
% Solve name exception P7TMS_015_Macular Cube 512x128_6-25-2009_10-51-44_OD_sn7564_cube_z
% Unsafe, need to be changed sometime
if strncmp(C{3}{1},'Macular',6)
    C1 = C(2:end);
    C1{1}{1}=[C{1}{1} '_' C1{1}{1}];
    C = C1;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

% Read the data
fid = fopen(filename);
if fid == -1
    error(['Could not open file ' filename '! Check to make sure the ' ...
           'path is correct.']);
end
A = fread(fid,'*uchar');
fclose(fid);

% Form the volume
num_axial = numel(A)/(num_ascan*num_bscan);
if rem(num_axial,1) ~= 0 
    error(['Problem reading file' filename])
end

vol_data = reshape(A,num_ascan,num_axial,num_bscan);

% Reshape - bscans are rotated
vol_data = permute(vol_data,[2 1 3]);
vol_data = flipdim(vol_data,1); % upside down
vol_data = flipdim(vol_data,2); % flip left and right
vol_data = flipdim(vol_data,3); % flip front to back

if strcmp(cube,'hidef')
    vol_data = flipdim(vol_data,2);
end

% Resolution in each direction - assume cube begins and ends at voxel
%   center
vol_info.vol_res = [cube_size(3)/(num_axial-1),cube_size(1)/(num_ascan-1),cube_size(2)/(num_bscan-1)];

% Try to read LSO fundus image (Line Scanning Ophthalmoscope)
filename_lso = strrep(filename,'cube_z.img','lslo.bin');    
if ~exist(filename_lso,'file')
    lso = [];
else
    fid = fopen(filename_lso);
    im = fread(fid,339968,'uchar');
    fclose(fid);
    lso = fliplr(reshape(im,[664 512]));
end
