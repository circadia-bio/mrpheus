function [ DATA HDR ] = ReadPhilipsScanPhysLog(filename, channels, skipprep)
%READPHILIPSSCANPHYSLOG load Philips MRI physiolog file
%
%   [ DATA HDR ] = READPHILIPSSCANPHYSLOG(FILENAME)
%   Reads the sample information (all 10 channels) from FILENAME. The 
%   samples are loaded as signed 32-bit integer values (int32). Sampling
%   rate is 496.03Hz for wireless VCG systems and 500Hz for the older wired sytems.
%
%   [ DATA HDR ] = READPHILIPSSCANPHYSLOG(FILENAME, CHANNELS)
%   Also read the log data, but only loads the signals specified by cell
%   structure CHANNELS. {v1raw, v2raw, v1, v2, ppu, resp, gx, gy, gz, mark}
%   Use {} or 'all' to load all channels (default), or an arbitrary invalid
%   name (i.e. 'none') to ignore the channel samples and only load markers.
%
%   [ DATA HDR ] = READPHILIPSSCANPHYSLOG(FILENAME, CHANNELS, SKIPPREP)
%   Also read the log data and specified channels, but skip preperation
%   fase (marked by empty comment (#) line)
%
%   DATA    is a structure with elements for the channel matrix C, 
%           marker table M and event index vectors I. 
%   DATA.C  32-bit integer matrix containing column vectors for
%           all resp. channels values. (i.e. each signal channel is stored 
%           as column vector). Note that the mark column contains
%           hexadecimal values (this will always be the last column).
%   DATA.M  is a 2-column table containing {marker,index} pairs of the
%           the markers detected in the marker channel. The first
%           column of M contains the mark values >0, and the second column 
%           contains the sample indices of the corresponding mark-value 
%           in the first column. Markers are bit-encoded values (i.e. each
%           event type corresponds to a single bit position). Currently the
%           following bits are known:
%                mark      bit pattern         description
%               0x0001 = 0000.0000.0000.0001 = Trigger ECG
%               0x0002 = 0000.0000.0000.0010 = Trigger PPU
%               0x0004 = 0000.0000.0000.0100 = Trigger Respiration
%               0x0008 = 0000.0000.0000.1000 = Measurement ('slice onset')
%               0x0010 = 0000.0000.0001.0000 = start of scan sequence (decimal 16)
%               0x0020 = 0000.0000.0010.0000 = end of scan sequence (decimal 32)
%               0x0040 = 0000.0000.0100.0000 = Trigger external
%               0x0080 = 0000.0000.1000.0000 = Calibration
%               0x0100 = 0000.0001.0000.0000 = Manual start
%               0x8000 = 1000.0000.0000.0000 = Reference ECG Trigger
%
% Copyright 2011 Academical Medical Center Amsterdam
% Created by Paul F.C. Groot

    if nargin<2
        channels = {};
    elseif ~iscell(channels)
        channels = { channels };
    end

    if nargin<3 || isempty(skipprep)
        skipprep = false;
    end
    
    bParseHeader = nargout>=2;
       
    fid=fopen(filename,'rt');
    if fid==-1
        error('Couldn''t open %s',filename);
    end
    
    if bParseHeader
        p{1} = '##\s*(.*),\sRelease\s(\w+)\s\(SWID (\d+)\)$';
        p{2} = '##\s*[\s\w]*(\d\d)-(\d\d)-(\d\d\d\d)\s(\d\d):(\d\d):(\d\d)$';
        p{3} = '##\s*([\s\-\d])+$';
        p{4} = '##\s*Dockable table\s*=\s*(\w+)$';
        iPattern = 1;
    end
    
    H = cell(4,1);
    COLUMN_NAMES = [];

    while ~feof(fid)
        str = fgetl(fid);
        if strncmp(str,'##',2)
            if bParseHeader
                T = regexp(str, p{iPattern}, 'tokens');
                if length(T)~=1 
                    warning('SCANPHYSLOG:invalidHeader','Invalid header at line 1: %s', str); 
                elseif iPattern<=length(p) 
                    H{iPattern}=T{1}; 
                    iPattern=iPattern+1;
                else
                    warning('SCANPHYSLOG:invalidHeader','Ignoring additional header line: %s', str);
                end;
            end
        elseif strncmp(str,'#',1)
            COLUMN_NAMES = regexp(str(2:end),'\w+','match'); 
            break;
        else
            error('Invalid header');
        end
    end

    if bParseHeader
        if length(H)>=1 && ~isempty(H{1})
            T = H{1};
            ID = struct('Site', T{1}, 'Release', T{2}, 'SWID', str2double(T{3}) );
        else
            ID = struct('Site', 'Unknown site', 'Release', 'unknown release', 'SWID', 0 );
        end

        if length(H)>=2 && ~isempty(H{2})
            T = H{2};
            DATETIME = struct('year', uint16(str2double(T{3})), 'month', uint16(str2double(T{2})), ...
                              'day',  uint16(str2double(T{1})), 'hour',  uint16(str2double(T{4})), ...
                              'min',  uint16(str2double(T{5})), 'sec',   uint16(str2double(T{6})));
        else
            DATETIME = struct('year',uint16(1900),'month',uint16(1),'day',uint16(1), ...
                              'hour',uint16(0),'min',uint16(0),'sec',uint16(0));
        end

        if length(H)>=3 && ~isempty(H{3})
            T = H{3};
            STATS = sscanf(T{1},'%d')';
        else
            STATS = NaN .* zeros(1,9);
        end

        if length(H)>=4 && ~isempty(H{4})
            T = H{4};
            DockableTable = strcmpi(T,'TRUE');
        else
            DockableTable = false;
        end

        HDR = struct('ID', ID, 'DATETIME', DATETIME, 'STATS', STATS, ...
            'DockableTable', DockableTable, 'COLUMN_NAMES', COLUMN_NAMES);
    end

    if skipprep
        startpos = ftell(fid);
        while ~feof(fid)
            str = fgetl(fid);
            if strcmp(strtrim(str),'#')
                break;
            end
        end
        if feof(fid)
            warning('SCANPHYSLOG:prepNotFound','Preparation fase not found: reading all samples');
            fseek(fid,startpos,'bof');
        end
    end
 
    nColumns = length(COLUMN_NAMES);
    if isempty(channels) || sum(strcmpi(channels,'all'))>0
        channels = COLUMN_NAMES;
    end
    channel_ordering = []; 
    S = logical(zeros(1,nColumns));
    for iChannel=1:length(channels)
        which_channel_vector = strcmpi( channels{iChannel}, COLUMN_NAMES );
        S = S | which_channel_vector;
        iChannelPos = find(which_channel_vector);
        if (iChannelPos>0)
            channel_ordering(end+1) = iChannelPos;
        end
    end
    column_ordering = zeros(1,length(channel_ordering));
    temp_channel_ordering = channel_ordering;
    for iChannel=1:length(temp_channel_ordering)
        iMin = find(temp_channel_ordering==min(temp_channel_ordering));
        column_ordering(iMin) = iChannel;
        temp_channel_ordering(iMin) = Inf; 
    end
 
    format = [ sprintf('%%%cd ',(S(1:end-1)==0).*'*') '%4s'];
    format = regexprep(format,'%[^*]d','%d');
    
    C=textscan(fid,format,'MultipleDelimsAsOne',1,'CommentStyle','#');
    fclose(fid);
    
    markers = uint16(hex2dec(C{end}));
    n = length(markers);
    onevec = ones(n,1,'uint16');
    I.VcgOnset      = uint32(find(bitand(markers,1*onevec)));
    I.PpuOnset      = uint32(find(bitand(markers,2*onevec)));
    I.TriggerResp   = uint32(find(bitand(markers,4*onevec)));
    I.Measurement   = uint32(find(bitand(markers,8*onevec)));
    I.ScannerStart  = uint32(find(bitand(markers,16*onevec)));
    I.ScannerStop   = uint32(find(bitand(markers,32*onevec)));
    I.TriggerExt    = uint32(find(bitand(markers,64*onevec)));
    I.Calibration   = uint32(find(bitand(markers,128*onevec)));
    I.RefTriggerVcg = uint32(find(bitand(markers,32768*onevec)));

    if S(end)
        C{end} = int32(markers);
    else
        C = C(1:end-1);
    end
    C = cell2mat(C(column_ordering));

    iMarkerIndices = uint32(find(markers>0));
    M = [uint32(markers(iMarkerIndices)), iMarkerIndices];
    
    DATA = struct('C', C, 'M', M, 'I', I);
end
