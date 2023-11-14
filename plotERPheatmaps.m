function plotERPheatmaps( setName, channels, smoothing )
%
% •.° Plot ERP Heatmaps °.•
% _________________________________________________________________________
%
% Plot single-trial event-related potential (ERP) heatmaps of all trials
% across multiple channels interleaved together, for all EEGLAB datasets
% named in common that are located in the Current Folder and sub-folders
% of the Current Folder, and plot a grand single-trial ERP heatmap of all 
% trials across multiple channels from all datasets interleaved together;
% with smoothing using a moving-average Gaussian window across trials
% (interleaving preserves the temporal structure of the trials)
%
% • Usage •
% -------------------------------------------------------------------------
% Set the Current Folder to a folder containing EEGLAB datasets to plot or
% to a folder of sub-folders containing EEGLAB datasets to plot
% 
% >> plotERPheatmaps( setName, channels, smoothing )
%
% Example:
% >> plotERPheatmaps( 'Memory', { 'Fz' 'FCz' }, 1.5 )
%
% •••( Function Inputs)
%
%   setName:   name in common to the datasets, for example 'Memory' or
%              'LTP', or '' for all datasets in the Current Folder and
%              sub-folders of the Current Folder
%
%   channels:  channel or channels to plot as any of the following
%              'ChannelLabel', "ChannelLabel", ChannelNumber,
%              { 'ChannelLabel1' 'ChannelLabel2' ... 'ChannelLabelN' }, or
%              [ ChannelNumber1 ChannelNumber2 ... ChannelNumberN ],
%              for example 'Fz', { 'Fz 'FCz' }, or [11 6]
%
%   smoothing: number of trials to smooth over (in standard deviations)
%              (optional input, default [] -> 1 SD trials)
%
% !!! Requires !!!
%
% EEGLAB
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
%
% King, R. (2023). plotERPheatmaps [MATLAB code]. GitHub.
%  https://github.com/neuroro/EEG/plotERPheatmaps.m


%% Startup
% -------------------------------------------------------------------------

% Initialise EEGLAB
eegStartup

% Introduction
disp( ' ' )
disp( '•.° Plot ERP Heatmaps °.•' )
disp( '_________________________________________________________________________' )
disp( ' ' )
disp( 'Plot single-trial event-related potential (ERP) heatmaps of all trials'  )
disp( 'across multiple channels interleaved together, for all EEGLAB datasets'  )
disp( 'named in common that are located in the Current Folder and sub-folders'  )
disp( 'of the Current Folder, and plot a grand single-trial ERP heatmap of all' )
disp( 'trials across multiple channels from all datasets interleaved together;' )
disp( 'with smoothing using a moving-average Gaussian window across trials'     )
disp( '(interleaving preserves the temporal structure of the trials)'           )
disp( ' ' )


%% EEGLAB dataset files
% -------------------------------------------------------------------------

% Check input
aWildError = [ 'Input a name in common to all EEGLAB datasets ' ...
               'located in the current folder and sub-folders ' ...
               'to plot ERP heatmaps' ];
if ~ischar( setName ) 
    if isstring( setName )
        setName = char( setName );
    else
        error( aWildError )
    end
end

% Wildcard dataset name
if isempty( setName )
    aWildDataset = '*.set';
else
    aWildDataset = [ '*' setName '*.set' ];
end

% Search for all datasets named in common located in the current folder and
% sub-folders of the current folder
fileStruct = dir( [ '**/' aWildDataset ] );
nFiles     = length( fileStruct );
fileList   = { fileStruct(:).name   };
folderList = { fileStruct(:).folder };


%% Channels
% -------------------------------------------------------------------------

% Check input: Channel names as a cell of characters
if iscell( channels )
    channels = cellstr( channels );
elseif ischar( channels )
    channels = { channels };
elseif isstring( channels )
    channels = { char( channels ) };

% Check input: Channel numbers as indices
elseif isnumeric( channels )
    iChannels = channels;

end

% Number of channels to plot
nChannels = length( channels );


%% Plot single-trial ERP heatmaps for each dataset
% -------------------------------------------------------------------------

% Report enumeration
reportText = 'Plotting single-trial ERP heatmaps';
if ~isnumeric( channels )
    theChannels = channels;
else
    for c = 1:nChannels
        theChannels{c} = [ 'channel ' num2str( channels(c) ) ];
    end
end
if nChannels == 1
    channelsText = 'at';
else
    channelsText = 'across';
end
for c = 1:nChannels
    channelsText = [ channelsText ' ' theChannels{c} ];
    if c < nChannels
        channelsText = [ channelsText ' &' ];
    end
end
if isempty( setName )
    filesText = [ 'in ' num2str( nFiles ) ' datasets' ];
else
    filesText = [ 'in ' num2str( nFiles ) ' ' setName ' datasets' ];
end
disp( [ reportText ' ' channelsText ' ' filesText ] )
disp( ' ' )

% Smoothing
if nargin < 2 || isempty( smoothing )
    smoothing      = 1;
end

% Pre-allocate
signals            = cell( 1, nFiles );
trialCounts        = zeros( 1, nFiles );

% Loop through: EEGLAB datasets
for f = 1:nFiles

    % File
    currentFile    = fileList{f};
    currentFolder  = folderList{f};
    currentDataset = extractBefore( currentFile, '.set' );

    % Load
    loaderargin    = { 'filename', currentFile,   ...
                       'filepath', currentFolder, ...
                       'verbose',  'off'          };
    EEG            = pop_loadset( loaderargin{:} );
    EEG            = eeg_checkset( EEG );

    % Channel indices
    if iscell( channels )
        iChannels  = eeg_chaninds( EEG, channels );
    end

    % EEG signal at the selected channels
    eegSignal      = EEG.data(iChannels,:,:);

    % Number of trials per dataset
    trialCounts(f) = size( eegSignal, 3 );

    % Interleave channels into trial epochs
    % Preserve the temporal structure of the trials
    if nChannels > 1
        tc = 0;
        for t = 1:trialCounts(f)
            for c = 1:nChannels

                % Index of interleaved trials + channels
                tc = tc + 1;

                % Current EEG signal
                signal = eegSignal(c,:,t);
                signal = squeeze( signal );

                % Store in interleaved array
                eegSignals(:,tc) = signal;

            end
        end
    else
        eegSignals = squeeze( eegSignal );
    end

    % Store in signals cell array
    signals{f}     = eegSignal;

    % Multiply the smoothing by the number of channels
    averageWidth   = smoothing * nChannels;

    % ERP image options
    plotTitle      = [ char( channels' ) currentDataset ];
    erpargin       = { 'avg_type',  'Gaussian', ...
                       'yerplabel', '\muV',     ...
                       'erp',       'on',       ...
                       'cbar',      'on'        };

    % If avg_type is set to 'Gaussian,' averageWidth is the standard
    % deviation (in units of trials, which may be non-integer) of the
    % Gaussian window used to smooth (vertically) with a moving-average.
    % Gaussian window extends three standard deviations below and three
    % standard deviations above window center (trials beyond window are
    % not incorporated into average).
    
    % Plot single-trial ERP heatmap
    figure
    erpimage( eegSignals, [], EEG.times, plotTitle, averageWidth, 0, erpargin{:} );

end % for: Files


%% Plot grand single-trial ERP heatmap from all datasets
% -------------------------------------------------------------------------

if nFiles > 1

    % Maximum number of trials across datasets
    maxTrials = max( trialCounts );

    % Interleave participants (individual datasets) and channels into trial
    % epochs to preserve the temporal structure of the trials
    tfc = 0;
    for t = 1:maxTrials
        for f = 1:nFiles
            for c = 1:nChannels

                % Handle equal and unequal number of trials across datasets
                if t <= trialCounts(f)

                    % Index of interleaved trials + participants + channels
                    tfc = tfc + 1;

                    % Current EEG signal
                    signal = signals{f};
                    signal = signal(c,:,t);
                    signal = squeeze( signal );

                    % Store in interleaved array
                    erps(:,tfc) = signal;

                end

            end
        end
    end
    
    % Multiply the smoothing by the number of participants (individual
    % datasets) and the number of channels
    averageWidth = smoothing * nFiles * nChannels;

    % ERP image options
    plotTitle    = [ 'Grand ERP heatmap at ' char( channels' ) ];
    erpargin     = { 'avg_type',  'Gaussian', ...
                     'yerplabel', '\muV',     ...
                     'erp',       'on',       ...
                     'cbar',      'on'        };
    
    % Plot grand ERP heatmap of all trials
    figure
    erpimage( erps, [], EEG.times, plotTitle, averageWidth, 0, erpargin{:} );

end


% _________________________________________________________________________
end



%%
% •.° EEGLAB Initialisation
% _________________________________________________________________________
%
function eegStartup

% Check if EEGLAB is being used
eeglabVariableUsage   = ( exist( 'EEG',    'var' ) && ~isempty( EEG )    && ~isempty( EEG.data ) )       || ...
                        ( exist( 'ALLEEG', 'var' ) && ~isempty( ALLEEG ) && ~isempty( ALLEEG(1).data ) ) || ...
                        ( exist( 'STUDY',  'var' ) && ~isempty( STUDY ) );
erplabVariableUsage   =   exist( 'ALLERP', 'var' ) && ~isempty( ALLERP );
existenceCommand      = "exist( 'globalvars', 'var' ) || exist( 'tmpEEG', 'var' )";
baseVariableExistence = evalin( "base", existenceCommand );

% Initialise EEGLAB
eeglab nogui

% Clear global variables unless they are being used
if ~eeglabVariableUsage
    clearvars -global ALLCOM ALLEEG CURRENTSET CURRENTSTUDY EEG LASTCOM PLUGINLIST STUDY
end
if ~erplabVariableUsage
    clearvars -global ALLERP
end

% Clear variables set in the Base Workspace unless they are being used
if ~baseVariableExistence
    evalin( "base", "clearvars globalvars tmpEEG" )
end

% Reset the Command Window without clearing
home


% _________________________________________________________________________
end


