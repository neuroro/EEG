function eegMultiScaleEntropyAverage( setName )
%
% •.° Composite Multi-Scale Entropy °.•
% _________________________________________________________________________
%
% Average the composite multi-scale entropy (CMSE) of electrical brain
% activity in all *MSE.mat files named in common that are located in the
% Current Folder and sub-folders of the Current Folder and save a
% <setName>AverageCMSE.mat file for all CMSE data in the average
%
% • Usage •
% -------------------------------------------------------------------------
% Function:
% >> eegMultiScaleEntropyAverage( setName )
%
% Examples:
% >> eegMultiScaleEntropyAverage
% >> eegMultiScaleEntropyAverage( '' )
% >> eegMultiScaleEntropyAverage( 'Rest' )
% 
% Input:
%   setName: Name in common to the *MSE.mat files, for example 'Memory' or
%              'Rest', or '' for all *MSE.mat files in the Current Folder
%              and sub-folders of the Current Folder
%              (optional input, default all *MSE.mat files)
%
% Output:
%   <setName>AverageCMSE.mat file containing the CMSE of all *MSE.mat files
%   in the average
%
% • Author •
% -------------------------------------------------------------------------
% Rohan O. C. King, 2023
% @neuroro
%
% Copyright 2023 Rohan King
%
% This program is free software: you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by the
% Free Software Foundation, either version 3 of the License, or (at your
% option) any later version.
%
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
% Public License for more details (https://www.gnu.org/licenses/).
%
% Please cite this code if you use it


% Introduction
disp( ' ' )
disp( '•.° Average Composite Multi-Scale Entropy °.•' )
disp( '_________________________________________________________________________' )
disp( ' ' )
disp( 'Average the composite multi-scale entropy (CMSE) of electrical brain'      )
disp( 'activity in all *MSE.mat files named in common that are located in the'    )
disp( 'Current Folder and sub-folders of the Current Folder and save an'          )
disp( '<InputName>AverageCMSE.mat file for all CMSE data in the average'          )
disp( ' ' )


%% Input handling
% -------------------------------------------------------------------------

% Wildcard file name
aWildError = [ 'Input a name in common to all *MSE.mat files ' ...
               'located in the Current Folder and sub-folders ' ...
               'of the Current Folder to average'               ];
if ~nargin
    setName = '';
end
if ~ischar( setName ) 
    if isstring( setName )
        setName = char( setName );
    else
        error( aWildError )
    end
end
if isempty( setName )
    aWildMatFile = '*.set';
else
    aWildMatFile = [ '*' setName '*MSE.mat' ];
end


%% Individual MSE files
% -------------------------------------------------------------------------

% Search for all .mat files named in common located in the current folder 
% and sub-folders of the current folder
MatFilesStruct = dir( [ '**/' aWildMatFile ] );
nMats          = length( MatFilesStruct );
matFileList    = { MatFilesStruct(:).name   };
matFolderList  = { MatFilesStruct(:).folder };

% Separator
if any( contains( matFileList, ' ' ) )
    separator  = ' ';
elseif any( contains( matFileList, '_' ) ) && ~any( contains( matFileList, ' ' ) )
    separator  = '_';
elseif any( contains( matFileList, '-' ) ) && ~any( contains( matFileList, { ' ' '_' } ) )
    separator  = '-';
else
    separator  = '';
end


%% Average across trials and merge into a single file
% -------------------------------------------------------------------------

% MSE Struct
CMSE.Entropy    = [];
firstMatFile    = fullfile( matFolderList{1}, matFileList{1} );
CMSE.Dimensions = 'EEGs x Channels x Scales';
try     % Data definition information should be stored in the *MSE.mat files
    EntropyData                 = load( firstMatFile );
    CMSE.Scales                 = EntropyData.timeScales;
    CMSE.TimeScaleLimits        = EntropyData.timeScaleLimits;
    CMSE.SamplingRate           = EntropyData.samplingRate;
    CMSE.ChannelCoordinates     = EntropyData.chanlocs;
    clear EntropyData
catch   % Otherwise assume the first scale = 1
    disp( [ '[' 8 'Warning: Data definition information is missing.'               ']' 8 ] )
    disp( [ '[' 8 'Scales are assumed to be numbered consecutively starting at 1.' ']' 8 ] )
    load( firstMatFile, 'cmse' )
    CMSE.Scales                 = [ 1 size( cmse, 3 ) ];
    try % Load an EEGLAB dataset with the same name as the first *MSE.mat file
        mseFileSuffix           = [ separator 'MSE.mat' ];
        dataName                = extractBefore( matFileList{1}, mseFileSuffix );
        dataset                 = [ dataName '.set' ];
        EEG                     = pop_loadset( 'filename', dataset,          ...
                                               'filepath', matFolderList{1}, ...
                                               'verbose', 'off'              );
        EEG                     = eeg_checkset( EEG );
        CMSE.SamplingRate       = EEG.srate;
        CMSE.ChannelCoordinates = EEG.chanlocs;
        timeScaleLimits{1}      = [ num2str( CMSE.Scales(1) * 1000 / EEG.srate ) ' ms' ];
        timeScaleLimits{2}      = [ num2str( CMSE.Scales(2) * 1000 / EEG.srate ) ' ms' ];
        CMSE.TimeScaleLimits    = timeScaleLimits;
        clear EEG
    catch
        disp( [ '[' 8 'Warning: Sampling rate and channel co-ordinates are missing.' ']' 8 ] )
    end
    disp( ' ' )
end
clear cmse

% Join per-dataset MSE.mat files together into an average CMSE.mat file
for f = 1:nMats

    % File
    currentFile     = matFileList{f};
    currentFolder   = matFolderList{f};
    currentFullFile = fullfile( currentFolder, currentFile );

    % Load MSE for each channel x trial x time-scale
    load( currentFullFile, 'cmse' );

    % Average across trials
    cmse = mean( cmse, 2, 'omitnan' );

    % Store in MSE struct
    CMSE.Entropy(f,:,:) = squeeze( cmse ); % Files x channels x time-scales

end

% Save
if isempty( setName )
    setName = 'Pooled';
end
fileName    = [ setName separator 'Average' separator 'CMSE' ];
save( fileName, "CMSE" )

% Finish time
mseRunTime( 'Finished at' )


% _________________________________________________________________________
end



%% MSE run time
% _________________________________________________________________________

function mseRunTime( message )

% Clock
clockTime = clock;
hours     = clockTime(4);
minutes   = clockTime(5);

% 12-hour
if hours && hours < 12
    meridian = 'a.m.';
elseif hours == 12
    meridian = 'p.m.';
elseif hours > 12
    hours    = hours - 12;
    meridian = 'p.m.';
else
    hours    = 12;
    meridian = 'a.m.';
end

% Trailling zero
if minutes < 10
    tm = '0';
else
    tm = '';
end

% Time text
hoursText   = num2str( hours   );
minutesText = num2str( minutes );
timeText    = [ hoursText ':' tm minutesText ' ' meridian ];

% Print
if exist( 'message', 'var' )
    disp( [ message ' ' timeText ] );
else
    disp( timeText );
end
disp( ' ' )


% _________________________________________________________________________
end


