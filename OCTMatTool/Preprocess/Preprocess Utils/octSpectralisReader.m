function [header, BScanHeader, slo, BScans] = octSpectralisReader(path, options)
%OPENVOL Read Heidelberg Engineering (HE) OCT raw files (VOL ending)
% [HEADER, BSCANHEADER, SLO, BSCANS] = OPENVOL(PATH, OPTIONS)
% This function performs volume OCT data (xxxx.vol) reading. 
% HEADER: Header information as described by HE. Struct with each entry
%   named by the HE conventions. 
% BSCANHEADER: B-scan header information. Struct with each entry named by 
%   the HE conventions. The entries consist of vectors/matrices, with one 
%   field per B-Scan. 
% SLO: Slo image as unsigned integers.
% BScans: BScans. 3D matrix with floating point data of the B-Scans.
% PATH: Filename of the VOL-file to read, with ending.
%
% OPTIONS: Various output possibilites, 
%   written in the options variable as string text, i.e. 'visu writemeta'
%   Possibilities: 
%       visu: The read-in is visualized.
%       visuseg: The HE segmentation is also visualized.
%       metawrite: A metafile with the header and B-scan header information 
%           is written out (if it does not exist already)
%       metawriteforce: An existing metafile is replaced
%       header: Only the Header and BScansHeader Information is read, 
%            not the image data  
%       nodisp: nothing is diplayed during read in
%
% Originally writen by Radim Kolar, Brno University, Czech Republic
% Modified by Markus Mayer, Pattern Recognition Lab, University of
% Erlangen-Nuremberg
% Modified by Kris Sheets, Retinal Cell Biology Lab, 
% Neuroscience Center of Excellence, LSU Health Sciences Center, 
% New Orleans
% Modified by Andrew Lang, Image Analysis and Communications Lab, Johns
% Hopkins University - Modifications to increase efficiency - 5/17/2012 
%
% Currently up to Version HSF-XXX-102 is supported, without Thickness Grid
%
% First final Version: March 2010
% This version was revised and commented in January 2012, 
% and modified again in December 2011
%
% You may use this code as you want. I would be grateful if you would go to
% our homepage look for articles that you find worth citing in your next
% publication:
% http://www.vutbr.cz/lide/radim-kolar-2796/publikace
% http://www5.informatik.uni-erlangen.de/en/our-team/mayer-markus
% Thanks, Radim & Markus

% $Id: octSpectralisReader.m,v 1.1 2015/01/12 18:54:06 andrew Exp $
 
% If only one argument is defined
if nargin==1,
    options = '';
end

%--------------------------------------------------------------------------
% Open the file

visu = 0;
visuseg = 0;

if numel(strfind(options, 'visu')) ~= 0
    visu = 1;
end
if numel(strfind(options, 'visuseg')) ~= 0
    visuseg = 1;
end

fid = fopen(path);
 
%--------------------------------------------------------------------------
% File header read

header.Version = fread( fid, 12, '*int8' );
header.SizeX = fread( fid, 1, '*int32' );
header.NumBScans = fread( fid, 1, '*int32' );
header.SizeZ = fread( fid, 1, '*int32' );
header.ScaleX = fread( fid, 1, '*double' );
header.Distance = fread( fid, 1, '*double' );
header.ScaleZ = fread( fid, 1, '*double' );
header.SizeXSlo = fread( fid, 1, '*int32' );
header.SizeYSlo = fread( fid, 1, '*int32' );
header.ScaleXSlo = fread( fid, 1, '*double' );
header.ScaleYSlo = fread( fid, 1, '*double' );
header.FieldSizeSlo = fread( fid, 1, '*int32' );
header.ScanFocus = fread( fid, 1, '*double' );
header.ScanPosition = char(fread( fid, 4, '*uchar' )');
header.ExamTime = fread( fid, 1, 'int64' );
header.ScanPattern = fread( fid, 1, '*int32' );
header.BScanHdrSize = fread( fid, 1, '*int32' );
header.ID = char(fread( fid, 16, '*uchar' )');
header.ReferenceID = char(fread( fid, 16, '*uchar' )');
header.PID = fread( fid, 1, '*int32' );
header.PatientID = char(fread( fid, 21, '*uchar' )');
header.Padding = fread( fid, 3, '*int8' );
header.DOB = fread( fid, 1, '*double' );
header.VID = fread( fid, 1, '*int32' );
header.VisitID = char(fread( fid, 24, '*uchar' )');
header.VisitDate = fread( fid, 1, '*double' );
header.GridType = fread( fid, 1, '*int32');
header.GridOffset = fread( fid, 1, '*int32');
header.Spare = fread( fid, 1832, '*int8' );
 
if numel(strfind(options, 'nodisp')) == 0
    disp(['---------------------------------------------']);
    disp(['           Version: ' char(header.Version')]);
    disp(['             SizeX: ' num2str(header.SizeX)]);
    disp(['         NumBScans: ' num2str(header.NumBScans)]);
    disp(['             SizeZ: ' num2str(header.SizeZ)]);
    disp(['            ScaleX: ' num2str(header.ScaleX) ' mm']);
    disp(['          Distance: ' num2str(header.Distance) ' mm']);
    disp(['            ScaleZ: ' num2str(header.ScaleZ) ' mm']);
    disp(['          SizeXSlo: ' num2str(header.SizeXSlo)]);
    disp(['          SizeYSlo: ' num2str(header.SizeYSlo)]);
    disp(['         ScaleXSlo: ' num2str(header.ScaleXSlo) ' mm']);
    disp(['         ScaleYSlo: ' num2str(header.ScaleYSlo) ' mm']);
    disp(['FieldSizeSlo (FOV): ' num2str(header.FieldSizeSlo) 'deg']);
    disp(['         ScanFocus: ' num2str(header.ScanFocus) ' dpt']);
    disp(['      ScanPosition: ' char(header.ScanPosition)]);
    disp(['          ExamTime: ' datestr(double(header.ExamTime(1)/(1e7*60*60*24)+584755+(2/24)))]);
    disp(['       ScanPattern: ' num2str(header.ScanPattern)]);
    disp(['      BScanHdrSize: ' num2str(header.BScanHdrSize) ' bytes']);
    disp(['                ID: ' char(header.ID)]);
    disp(['       ReferenceID: ' char(header.ReferenceID)]);
    disp(['               PID: ' num2str(header.PID)]);
    disp(['         PatientID: ' char(header.PatientID)]);
    disp(['               DOB: ' datestr(header.DOB+693960)]);
    disp(['               VID: ' num2str(header.VID)]);
    disp(['           VisitID: ' char(header.VisitID)]);
    disp(['         VisitDate: ' datestr(double(header.VisitDate+693960))]);
    disp(['          GridType: ' num2str(header.GridType)]);
    disp(['        GridOffset: ' num2str(header.GridOffset)]);
    disp(['---------------------------------------------']);
end

if nargout == 1
    fclose( fid );
    return
end

status = fseek( fid, 2048, -1 );
 
%--------------------------------------------------------------------------
% SLO image read

if(~any(strcmp(options, 'header')))
    slo = fread(fid, header.SizeXSlo*header.SizeYSlo, '*uint8');
    slo = reshape(slo, header.SizeXSlo, header.SizeYSlo);
%     slo = imrotate(slo, 90);
%     slo = flipud(slo);
    slo = slo';
else
    slo = [];
end
% slo = uint8(slo);
 
%--------------------------------------------------------------------------
% Display the image

if visu==1
    scrsz = get(0,'ScreenSize');
    figure('Position',[1 0 scrsz(3) scrsz(4)-70]);
    subplot(1,2,1);
    imshow(slo);
end
 
%--------------------------------------------------------------------------
% B-scans read (including image data and scanning position)

status = fseek( fid, 2048+(header.SizeXSlo*header.SizeYSlo), -1 );
 
if(~any(strcmp(options, 'header')))
    BScans=zeros(header.SizeZ, header.SizeX ,header.NumBScans, 'single');
else
    BScans= [];
end

BScanHeader.StartX = zeros(1, header.NumBScans, 'double');
BScanHeader.StartY = zeros(1, header.NumBScans, 'double');
BScanHeader.EndX = zeros(1, header.NumBScans, 'double');
BScanHeader.EndY = zeros(1, header.NumBScans, 'double');
BScanHeader.NumSeg = zeros(1, header.NumBScans, 'int32');
BScanHeader.Quality = zeros(1, header.NumBScans, 'single');
BScanHeader.Shift = zeros(1, header.NumBScans, 'int32');
BScanHeader.ILM = zeros(header.NumBScans,header.SizeX, 'single');
BScanHeader.RPE = zeros(header.NumBScans,header.SizeX, 'single');
BScanHeader.NFL = zeros(header.NumBScans,header.SizeX, 'single');
 
for zz = 1:header.NumBScans
    ii = zz - 1;
    
    %----------------------------------------------------------------------
    %BScan position in SLO image
    
    status = fseek( fid, 16+2048+(header.SizeXSlo*header.SizeYSlo)+(ii*(header.BScanHdrSize+header.SizeX*header.SizeZ*4)), -1 );
    StartX = fread( fid, 1, '*double' ); 
    StartY = fread( fid, 1, '*double' );  
    EndX = fread( fid, 1, '*double' );
    EndY = fread( fid, 1, '*double' );  
    NumSeg = fread( fid, 1, '*int32' );
    OffSeg = fread( fid, 1, '*int32' );
    Quality = fread( fid, 1, '*float32' );
    Shift = fread( fid, 1, '*int32' );
    
    BScanHeader.StartX(zz) = StartX;
    BScanHeader.StartY(zz) = StartY;
    BScanHeader.EndX(zz) = EndX;
    BScanHeader.EndY(zz) = EndY;
    BScanHeader.NumSeg(zz) = NumSeg;
    BScanHeader.Shift(zz) = Shift;
    BScanHeader.Quality(zz) = Quality;
 
    %----------------------------------------------------------------------
    % Display the images
    
    if visu==1
        StartX_px = round(StartX/header.ScaleXSlo); %StartX in pixels
        StartY_px = round(StartY/header.ScaleYSlo); %StartY in pixels
        EndX_px = round(EndX/header.ScaleXSlo); %EndX in pixels
        EndY_px = round(EndY/header.ScaleYSlo); %EndY in pixels
        subplot(1,2,1); hold on, line([StartX_px EndX_px],[StartY_px EndY_px],'color','b');
    end
    
    %----------------------------------------------------------------------
    % BScan reading
    
    if(~any(strcmp(options, 'header')))
        status = fseek( fid, header.BScanHdrSize+2048+(header.SizeXSlo*header.SizeYSlo)+(ii*(header.BScanHdrSize+header.SizeX*header.SizeZ*4)), -1 );
        oct = fread( fid, header.SizeX*header.SizeZ, '*float32' );
        oct = reshape( oct, header.SizeX, header.SizeZ );
        BScans(:,:,zz)=oct';
    end
     
    %----------------------------------------------------------------------
    % Display the images
    
    if visu==1
        subplot(1,2,2);
        imshow( BScans(:,:,zz).^0.25,[0 1]);
        drawnow
    end
 
    %----------------------------------------------------------------------
    % Segmentation reading
    
    status = fseek( fid, 256+2048+(header.SizeXSlo*header.SizeYSlo)+(ii*(header.BScanHdrSize+header.SizeX*header.SizeZ*4)), -1 );
    seg = (fread( fid, NumSeg*header.SizeX, '*float' ))';
    BScanHeader.ILM(zz,:) = seg(1:header.SizeX);
    BScanHeader.RPE(zz,:) = seg(header.SizeX+1:header.SizeX*2);
    if NumSeg == 3
        BScanHeader.NFL(zz,:) = seg(header.SizeX*2+1:header.SizeX*3);
    end
end

%----------------------------------------------------------------------
% This is for ILM and RPE INVALID values interpolation using median of 3x3 
% neigborhood around INVALID value.

BScanHeader.ILM(BScanHeader.ILM>1e6) = nan;
BScanHeader.RPE(BScanHeader.RPE>1e6) = nan;
if NumSeg == 3
    BScanHeader.NFL(BScanHeader.NFL>1e6) = nan;
end

% max_ILM = max(max(BScanHeader.ILM));
% [N,M]=size(BScanHeader.ILM);
% for j=1:N
%     for i=1:M
%         if BScanHeader.ILM(j,i) == max_ILM
%             ILM_sec = imcrop(BScanHeader.ILM, [i-1 j-1 3 3]);
%             BScanHeader.ILM(j,i) = median(find(ILM_sec ~= max_ILM));
%         end
%         BScanHeader.ILM(j,i) = header.SizeZ-BScanHeader.ILM(j,i);
%     end
% end
% 
% max_RPE = max(max(BScanHeader.RPE));
% [N,M]=size(BScanHeader.RPE);
% for j=1:N
%     for i=1:M
%         if BScanHeader.RPE(j,i) == max_RPE
%             RPE_sec = imcrop(BScanHeader.RPE, [i-1 j-1 3 3]);
%             BScanHeader.RPE(j,i) = median(find(RPE_sec ~= max_RPE));
%         end
%         BScanHeader.RPE(j,i) = header.SizeZ-BScanHeader.RPE(j,i);
%     end
% end

%--------------------------------------------------------------------------
% Display the segmentation

if visuseg == 1,
    figure;
    mesh(double(BScanHeader.ILM));
    hold on
    mesh(double(BScanHeader.RPE));
    if NumSeg == 3
        hold on
        mesh(double(BScanHeader.NFL));
    end
end

fclose( fid );

% %--------------------------------------------------------------------------
% % Write out a .meta file (text file with the header information
% 
% if numel(strfind(options, 'metawrite')) ~= 0
%     metafilename = [path(1:end-3) 'meta'];
%     
%     if((exist(metafilename) == 0) || numel(strfind(options, 'metawriteforce')) ~= 0)
%         fidW = fopen(metafilename, 'w');
%         
%         metaWriteHelper(fidW, 'ExamTime', header.ExamTime);
%         metaWriteHelper(fidW, 'ScanFocus', header.ScanFocus);
%         metaWriteHelper(fidW, 'ScanPosition', header.ScanPosition, 'str');
%         metaWriteHelper(fidW, 'ScanPattern', header.ScanPattern());
%         metaWriteHelper(fidW, 'ID', header.ID, 'str');
%         metaWriteHelper(fidW, 'ReferenceID', header.ReferenceID, 'str');
%         metaWriteHelper(fidW, 'PID', header.PID()  );
%         metaWriteHelper(fidW, 'PatientID', header.PatientID, 'str');
%         metaWriteHelper(fidW, 'DOB', header.DOB);
%         metaWriteHelper(fidW, 'VID', header.VID);
%         metaWriteHelper(fidW, 'VisitID', header.VisitID, 'str');
%         metaWriteHelper(fidW, 'VisitDate', header.VisitDate);
%         metaWriteHelper(fidW, 'OctSize', [header.SizeX header.SizeZ header.NumBScans]);
%         metaWriteHelper(fidW, 'OctScale', [header.ScaleX header.ScaleZ header.Distance]);
%         metaWriteHelper(fidW, 'SloSize', [header.SizeXSlo header.SizeYSlo]);
%         metaWriteHelper(fidW, 'SloScale', [header.ScaleXSlo header.ScaleYSlo]);
%         metaWriteHelper(fidW, 'SloFieldSize', header.FieldSizeSlo);
%         metaWriteHelper(fidW, 'BScanStartX', BScanHeader.StartX);
%         metaWriteHelper(fidW, 'BScanStartY', BScanHeader.StartY);
%         metaWriteHelper(fidW, 'BScanEndX', BScanHeader.EndX);
%         metaWriteHelper(fidW, 'BScanEndY', BScanHeader.EndY);
%         
%         fclose( fidW );
%     end
% end
% end
% 
% function metaWriteHelper(fidW, tag, data, mode)
% if nargin < 4
%     mode = 'num';
% end
% if strcmp(mode, 'num')
%     str = sprintf('%g ', data);
% else
%     str = sprintf('%s ', data);
% end
% outputstring = [tag ' ' str];
% outputstring = deblank(outputstring);
% fprintf(fidW, '%s\n', outputstring);
% end
