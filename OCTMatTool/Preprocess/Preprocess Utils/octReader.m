function [header, oct_vol, slo, BScanHeader, scanner_type] = octReader(filename, options)
% Read OCT data. Determine scanner based on file extension
% Outputs are based on original Spectralis file reader
%
% Updates:
%   1/11/2016 - changed cirrus output scan_date to match spectralis
%               VisitDate as opposed to ExamTime since the ExamTime output
%               of the Spectralis is not right
%             - Fixed Spectralis VisitDate and DOB outputs so that datestr
%               can be used directly on them
%             - Removed trailing whitespace from header.ScanPosition
%               ('OD  ' -> 'OD')
%
%   5/26/2016 - added scan angle to the header (for posterior pole and
%               vertical scans)
%
% Andrew Lang

if nargout == 1
    % If we only want the header, then don't need to read the volumetric
    % data making it more efficient
    header_only = true;
else
    header_only = false;
end

[~, ~, ext] = fileparts(filename);

BScanHeader = [];
slo = [];
switch lower(ext)
    case '.vol'
        % Spectralis
        if header_only
            header = octSpectralisReader(filename,'-nodisp');
        else
            [header, BScanHeader, slo, oct_vol] = octSpectralisReader(filename,'-nodisp');
            oct_vol = sqrt(sqrt(oct_vol)); % 1sec faster than ^0.25
            
            % Remove infinite values that occur when the image is clipped
            % by the scanner (due to imgae registration and eye movement
            % probably)
            oct_vol(oct_vol>1) = 0;
        end
        
        % Correct VisitDate so it can be used directly by datestr
        header.VisitDate = header.VisitDate+693960;
        % Remove trailing spaces from ScanPosition
        header.ScanPosition = deblank(header.ScanPosition);
        % Correct date of birth so it can be used directly by datestr
        header.DOB = header.DOB+693960;
        
        if ~header_only
            % Get angle from the horizontal that the OCT images were acquired
            % relative to the fundus image (note: positive angles go down since
            % we use an image coordinate system with y-axis down)
            xp = [BScanHeader.StartX(1) BScanHeader.EndX(1)];
            yp = [BScanHeader.StartY(1) BScanHeader.EndY(1)];
            header.angle = atan2d(yp(2)-yp(1),xp(2)-xp(1));
        end
        
        scanner_type = 'spectralis';
        
    case '.img'
        % Cirrus
        if header_only
            vol_info = octCirrusReaderInfo(filename);
        else
            [oct_vol, vol_info, slo] = octCirrusReader(filename);
            oct_vol = im2double(oct_vol);            
            header.SizeX = size(oct_vol,2);
            header.NumBScans = size(oct_vol,3);
            header.SizeZ = size(oct_vol,1);
        end
        
        header.PID = vol_info.pid;
        
        header.ScaleX = vol_info.vol_res(2);
        header.Distance = vol_info.vol_res(3);
        header.ScaleZ = vol_info.vol_res(1);
        header.VisitDate = vol_info.scan_date;
        header.ExamTime = vol_info.scan_date;
        header.ScanPosition = vol_info.eye_side;
        header.ScanType = vol_info.scan_type;
        
        header.angle = 0;
        
        scanner_type = 'cirrus';
        
    case '.oct'
        % Bioptigen
        [oct_vol,b_header] = octBioptigenReader(filename,false,true);
        oct_vol = im2double(oct_vol);   
        oct_vol = (oct_vol - min(oct_vol(:)))/(max(oct_vol(:))-min(oct_vol(:)));
        % Don't want/need 2048 pixels!
        if b_header.lineLength > 1024
            b_header.lineLength = round(b_header.lineLength/2);
            b_header.scanDepth = b_header.scanDepth/2;
            oct_vol = oct_vol(1:b_header.lineLength,:,:);
            oct_vol = im2double(oct_vol); 
        end
        
        header.ScaleX = b_header.azScanLength/b_header.lineCount;
        header.ScaleZ = b_header.scanDepth/b_header.lineLength;
        header.Distance = b_header.elScanLength/b_header.frameCount;
        
        header.SizeX = b_header.lineCount;
        header.SizeZ = b_header.lineLength;
        header.NumBScans = b_header.frameCount;
        
        if strfind(filename,'_OD_')
            header.ScanPosition = 'OD';
        else
            header.ScanPosition = 'OS';
        end
        
        scanner_type = 'bioptigen';
        
    case '.mat'
        S = load(filename);
        
        oct_vol = S.oct_vol;
        header = S.header;
        
        header.NumBScans = size(oct_vol,3);
        header.SizeX = size(oct_vol,2);
        header.SizeZ = size(oct_vol,1);
        
        if isfield(S,'slo')
            slo = S.slo;
        end
        if isfield(S.header,'scanner_type')
            scanner_type = S.header.scanner_type;
        else
            scanner_type = 'spectralis'; % default?
        end
    case '.dcm'        
        info = dicominfo(filename);
        header.NumBScans = info.NumberOfFrames;
        header.SizeX = info.Width;
        header.SizeZ = info.Height;
        header.ScaleX = info.SharedFunctionalGroupsSequence.Item_1.PixelMeasuresSequence.Item_1.PixelSpacing(2);
        header.ScaleZ = info.SharedFunctionalGroupsSequence.Item_1.PixelMeasuresSequence.Item_1.PixelSpacing(1);
        header.Distance = info.SharedFunctionalGroupsSequence.Item_1.PixelMeasuresSequence.Item_1.SliceThickness;
        left_right = info.SharedFunctionalGroupsSequence.Item_1.FrameAnatomySequence.Item_1.FrameLaterality;
        if strcmp(left_right,'R')
            header.ScanPosition = 'OD';
        elseif strcmp(left_right,'L')
            header.ScanPosition = 'OS';
        else
            warning('Could not find laterality in DICOM header');
            header.ScanPosition = '';
        end
        scanner_type = lower(info.ManufacturerModelName);
        if contains(scanner_type, 'cirrus')
            % Dataset from Anna seems to have unstandard input
            scanner_type = 'cirrus';
            header.Distance = 0.0472;
        end
        if isfield(info.SharedFunctionalGroupsSequence.Item_1, ...
                   'PlaneOrientationSequence')
            orient = info.SharedFunctionalGroupsSequence ...
           .Item_1.PlaneOrientationSequence.Item_1.ImageOrientationPatient;
        else
            orient = [1;0;0;0;1;0];
        end
        
        if isequal(orient, [1;0;0;0;1;0])
            header.angle = 0;
        else
            warning('Angle may be non-zero, please check header information');
        end
        
        oct_vol = double(squeeze(dicomread(filename)))./255;
    otherwise
        error(['Invalid file format for reading OCT data: ' ext])
end
