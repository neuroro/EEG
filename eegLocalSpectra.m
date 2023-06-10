function eegLocalSpectra( channels )
%
% •.° Localised Event-Related Spectra °.•
% _________________________________________________________________________
%
% Extract event-related spectra at the specified electrode or average of an
% electrode cluster from all TimeFrequencyData<Participant><Condition>.mat
% files containing data centred on more than one event, which are located
% in the Current Folder and sub-folders of the Current Folder then save
% data in one TimeFrequency<Location>.mat file
%
% • Usage •
% -------------------------------------------------------------------------
% Set the Current Folder to the location of the datasets or the base folder
%
% Function:
% >> eegLocalSpectra( channels )
%
% Examples:
% >> eegLocalSpectra()
% >> eegLocalSpectra( 6 )
% >> eegLocalSpectra( 'FCz' )
% >> eegLocalSpectra( [5 6 11 12] )
% >> eegLocalSpectra( { 'Fz' 'FCz' 'F1' 'F2' } )
%
% Input:
%   Extract event-related spectra from a channel name or index; from the
%   average of a cluster of channels as a cell array of channel names or as
%   a vector of channel indices; or from all channels as (), '', [], or 0
%   (optional input, default all)
%
% Output:
%   TimeFrequency<Location>.mat file containing, for each event-related
%   windowing, pooled decompositions from all participants per condition,
%   ready for peak finding or cluster statistics, and grand averages per
%   condition ready for plotting
%
% • Requires •
% -------------------------------------------------------------------------
% Time-frequency decompositions of EEG data produced by eegTimeFrequency.m
% or eeg3TimeFreqyency.m in TimeFrequencyData<Participant><Condition>.mat
% files
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


% Introduction
disp( ' ' )
disp( '•.° Localised Event-Related Spectra °.•' )
disp( '_________________________________________________________________________' )
disp( ' ' )
disp( 'Extract event-related spectra at the specified electrode or average of an' )
disp( 'electrode cluster from all TimeFrequencyData<Participant><Condition>.mat'  )
disp( 'files containing data centred on more than one event, which are located'   )
disp( 'in the Current Folder and sub-folders of the Current Folder then save'     )
disp( 'data in one TimeFrequency<Location>.mat file'                              )
disp( ' ' )

codeRunTime( 'Started at' )
disp( ' ' )


%% Input
% -------------------------------------------------------------------------

if ~nargin || ( isscalar( channels ) && ~channels ) || isempty( channels )
    channels     = [];
end


%% Time-frequency decomposition .mat files
% -------------------------------------------------------------------------

% Wildcard file name
aWildDataset     = 'TimeFrequencyData*.mat';

% Search for all datasets named in common located in the Current Folder and
% sub-folders of the Current Folder
fileStruct       = dir( [ '**/' aWildDataset ] );
nFiles           = length( fileStruct );
fileList         = { fileStruct(:).name   };
folderList       = { fileStruct(:).folder };


%% Determine participants and conditions represented by the .mat files
% -------------------------------------------------------------------------

% Pre-allocate
participantList  = {}; % Participant number files list
participants     = []; % Unique participant numbers
iParticipant     = 0;
conditionList    = {}; % Condition files list
conditions       = {}; % Unique conditions
iCondition       = 0;

% Loop through: Files
for f = 1:nFiles

    % Current file
    currentFile        = fileList{f};

    % Current participant number
    [ iX1, iXn ]       = regexp( currentFile, '[0-9]+' );
    iParticipantNumber = iX1:iXn;
    participantNumber  = currentFile(iParticipantNumber);

    % Current condition
    currentCondition   = extractBefore( currentFile(iXn+1:end), '.mat' );


    % Build a vector of unique participant numbers
    % ---------------------------------------------------------------------

    % Check that this participant has not already been included
    if ~any( strcmp( participantList, participantNumber ) )

        % Current index in the participant number list
        iParticipant               = iParticipant + 1;

        % Add current participant number to the list if not already included
        participants(iParticipant) = str2double( participantNumber );                               %#ok

    end

    % Add current file participant to the participant files list
    participantList{f} = participantNumber;                                                         %#ok


    % Build a vector of unique conditions
    % ---------------------------------------------------------------------

    % Check that this condition has not already been included
    if ~any( strcmp( conditionList, currentCondition ) )

        % Current index in the condition list
        iCondition             = iCondition + 1;

        % Add current condition to the list given not already included
        conditions{iCondition} = currentCondition;                                                  %#ok

    end

    % Add current file condition to the condition files list
    conditionList{f} = currentCondition;                                                            %#ok


end

% Sample size N participants
N = length( participants );

% Number of conditions
nConditions = length( conditions );


%% Decomposition parameters
% -------------------------------------------------------------------------

% Determine parameters from the first file
Decomposition = load( fullfile( folderList{1}, fileList{1} ) );

% Event-centred windowings
eventCentres  = fieldnames( Decomposition );
nCentres      = length( eventCentres );

% Store conditions in struct
for w = 1:nCentres
    LocalSpectra.(eventCentres{w}).Conditions = conditions;
end

% Frequencies
for w = 1:nCentres
    LocalSpectra.(eventCentres{w}).Frequencies = Decomposition.(eventCentres{w}).Frequencies;
end
nFrequencies  = length( LocalSpectra.(eventCentres{1}).Frequencies );

% Times
%   Time points vary per window centre
for w = 1:nCentres
    LocalSpectra.(eventCentres{w}).Times = Decomposition.(eventCentres{w}).Times;
    nTimes.(eventCentres{w})             = length( LocalSpectra.(eventCentres{w}).Times );
end


% Channels
% -------------------------------------------------------------------------

% Channels in the decomposition
iChannelSet        = Decomposition.(eventCentres{1}).Channels;           % Channel indices of the EEG data (and channel co-ordinates)
ChannelCoordinates = Decomposition.(eventCentres{1}).ChannelCoordinates; % Channel co-ordinates

% Number of selected channels
nChannels          = length( channels );

% Selected channels from the decomposition
if ~isempty( channels )

    % Loop through: Channels to extract
    for ch = 1:nChannels

        % Find the index of each selected channel in the set of channels
        % (indices) that are in the decomposition
        if isnumeric( channels )
            iCurrentChannel = channels(ch);
        elseif iscell( channels )
            currentChannel  = channels{ch};
            iCurrentChannel = find( strcmp( { ChannelCoordinates(:).labels }, currentChannel ) );
        end
        iChannels(ch) = find( iChannelSet == iCurrentChannel );                                     %#ok

    end

% All channels in the decomposition
else
    iChannels = 1:length( iChannelSet );

end

% Selected channel names
channelNames  = { ChannelCoordinates(iChannelSet(iChannels)).labels };

% Store channel names and co-ordinates in struct
for w = 1:nCentres
    LocalSpectra.(eventCentres{w}).Channels           = channelNames;
    LocalSpectra.(eventCentres{w}).ChannelCoordinates = ChannelCoordinates;
end


%% Display information
% -------------------------------------------------------------------------

% Files
disp( 'Extracting localised event-related spectra...' )
disp( [ 'Windowed relative to ' num2str( nCentres ) ' experimental events' ] )
disp( [ 'From ' num2str( nFiles ) ' time-frequency decompositions' ] )
disp( [ 'Representing ' num2str( N ) ' participants and ' num2str( nConditions ) ' conditions' ] )


% Channels
% -------------------------------------------------------------------------

% Multiple selected
if nChannels > 1

    % Build list of channel names to display
    channelsText = [];
    for ch = 1:nChannels
        channelsText = [ channelsText channelNames{ch} ' ' ];                                       %#ok
    end

    disp( [ 'Localised to the average of ' num2str( nChannels ) ' channels as a cluster:' ] )
    disp( channelsText )

% One selected
elseif nChannels == 1
    disp( [ 'Localised to channel ' channelNames{1} ] )

% All selected
else
    disp( [ 'Localised to the average of all ' num2str( nChannels ) ' channels in the decompositions' ] )

end

disp( ' ' )

clear Decomposition


%% Construct time-frequency arrays
% -------------------------------------------------------------------------

% Pre-allocate
for w = 1:nCentres
    LocalSpectra.(eventCentres{w}).SpectralPower  = NaN( N, nConditions, nFrequencies, nTimes.(eventCentres{w}) );
    LocalSpectra.(eventCentres{w}).PhaseCoherence = NaN( N, nConditions, nFrequencies, nTimes.(eventCentres{w}) );
end

% Loop through: Files per participant per condition
for f = 1:nFiles

    % File
    currentFile        = fileList{f};
    currentFolder      = folderList{f};
    currentFilePath    = fullfile( currentFolder, currentFile );
    
    % Participant
    currentParticipant = participantList{f};                               % Giving e.g. '15'
    participantNumber  = str2double( currentParticipant );                 % Giving e.g. 15
    iParticipant       = find( participants == participantNumber );        % Index of the participant number (in the list of unique participant numbers and in the decomposition)
    if isempty( iParticipant )
        error( 'Participant not found = floating point error?' )
    end

    % Condition
    currentCondition   = conditionList{f};
    iCondition         = find( strcmp( conditions, currentCondition ) );   % Index of the condition (in the list of unique conditions and in the decomposition)

    % Load time-frequency data
    Decomposition      = load( currentFilePath );                          % Struct with event-centred named fields

    % Loop through: Window centres
    for w = 1:nCentres

        currentCentre  = eventCentres{w};

        % Select channels
        spectralPower  = Decomposition.(currentCentre).SpectralPower(iChannels,:,:);  % Channels x frequencies x times
        phaseCoherence = Decomposition.(currentCentre).PhaseCoherence(iChannels,:,:); % Channels x frequencies x times
        
        % Average across channels
        spectralPower  = mean( spectralPower,  1, 'omitnan' );
        phaseCoherence = mean( phaseCoherence, 1, 'omitnan' );

        % Clean up
        spectralPower  = squeeze( spectralPower  );
        phaseCoherence = squeeze( phaseCoherence );

        % Store in arrays
        LocalSpectra.(currentCentre).SpectralPower(iParticipant,iCondition,:,:)  = spectralPower;
        LocalSpectra.(currentCentre).PhaseCoherence(iParticipant,iCondition,:,:) = phaseCoherence;
        
    end % for: Window centres

end % for: Files


%% Grand average spectra per condition
% -------------------------------------------------------------------------

% Loop through: Window centres
for w = 1:nCentres

    currentCentre = eventCentres{w};

    % Mean
    LocalSpectra.(currentCentre).GrandAverage.SpectralPower  ...
         = mean( LocalSpectra.(currentCentre).SpectralPower,  1, 'omitnan' );
    LocalSpectra.(currentCentre).GrandAverage.PhaseCoherence ...
         = mean( LocalSpectra.(currentCentre).PhaseCoherence, 1, 'omitnan' );

    % Clean up
    LocalSpectra.(currentCentre).GrandAverage.SpectralPower  ...
         = squeeze( LocalSpectra.(currentCentre).GrandAverage.SpectralPower  );
    LocalSpectra.(currentCentre).GrandAverage.PhaseCoherence ...
         = squeeze( LocalSpectra.(currentCentre).GrandAverage.PhaseCoherence );

end % for: Window centres


%% Location naming
% -------------------------------------------------------------------------

% Important channel names
Fz  = 'Fz';
FCz = 'FCz';
CPz = 'CPz';
Pz  = 'Pz';
POz = 'POz';
Oz  = 'Oz';
P7  = 'P7';
P8  = 'P8';

% EGI channels
if any( contains( { ChannelCoordinates(:).labels }, 'E' ) )
    Fz  = 'E11';
    FCz = 'E6';
    CPz = 'E55';
    Pz  = 'E62';
    POz = 'E72';
    Oz  = 'E75';
end

% Important channel indices
iFz  = find( strcmp( { ChannelCoordinates(:).labels }, Fz  ) );
iFCz = find( strcmp( { ChannelCoordinates(:).labels }, FCz ) );
iCPz = find( strcmp( { ChannelCoordinates(:).labels }, CPz ) );
iPz  = find( strcmp( { ChannelCoordinates(:).labels }, Pz  ) );
iPOz = find( strcmp( { ChannelCoordinates(:).labels }, POz ) );
iOz  = find( strcmp( { ChannelCoordinates(:).labels }, Oz  ) );
iP7  = find( strcmp( { ChannelCoordinates(:).labels }, P7  ) );
iP8  = find( strcmp( { ChannelCoordinates(:).labels }, P8  ) );

% Cluster of electrodes name
if length( channels ) > 1 || ( isempty( channels ) && length( iChannelSet ) <= 7 )

    iSelectedChannels = iChannelSet(iChannels);

    % Frontal midline cluster including FCz & Fz
    if any( iSelectedChannels == iFCz ) && any( iSelectedChannels == iFz )

        locationName = 'FMCluster';

    % Frontal midline cluster including FCz
    elseif any( iSelectedChannels == iFCz ) && ~any( iSelectedChannels == iFz )
        locationName = 'FCzCluster';

    % Frontal midline cluster including Fz
    elseif ~any( iSelectedChannels == iFCz ) && any( iSelectedChannels == iFz )
        locationName = 'FzCluster';

    % Parietal cluster including Pz & CPz
    elseif any( iSelectedChannels == iPz ) && ~any( iSelectedChannels == iCPz )
        locationName = 'PcMCluster';

    % Parietal cluster including Pz
    elseif any( iSelectedChannels == iPz ) && ~any( iSelectedChannels == iCPz ) && ~any( iSelectedChannels == iPOz )
        locationName = 'PzCluster';

    % Parietal cluster including Pz & POz
    elseif any( iSelectedChannels == iPz ) && any( iSelectedChannels == iPOz )
        locationName = 'PoMCluster';

    % Occipitoparietal cluster including POz
    elseif any( iSelectedChannels == iPOz ) && ~any( iSelectedChannels == iPz ) && ~any( iSelectedChannels == iOz )
        locationName = 'POzCluster';

    % Occipital cluster including Oz & POz
    elseif any( iSelectedChannels == iOz ) && any( iSelectedChannels == iPOz )
        locationName = 'OpMCluster';

    % Occipital cluster including Oz
    elseif any( iSelectedChannels == iOz ) && ~any( iSelectedChannels == iPOz )
        locationName = 'OzCluster';

    % Bilateral parietal cluster including P7 & P8
    elseif any( iSelectedChannels == iP7 ) && any( iSelectedChannels == iP8 )
        locationName = 'BiPCluster';

    % Individual electrodes named Channel1Channel2...ChannelN
    else
        locationName = char( channelNames )';

    end

% Single electrode name
elseif length( channels ) == 1
    locationName = channelNames{1};

% All channels
else
    if length( iChannelSet ) <= 10
        locationName = 'Cluster';
    else
        locationName = 'Average';
    end

end % if Number of selected electrodes


%% Save
% -------------------------------------------------------------------------

% File name
fileName = [ 'TimeFrequency' locationName ];

% Save TimeFrequency<Location>.mat 
save( fileName, '-struct', 'LocalSpectra', '-v7.3' )

% Finish time
codeRunTime( 'Data generation completed at' )


% _________________________________________________________________________
end



%%
% •.° Code Run Time °.•
% _________________________________________________________________________
function codeRunTime( message )

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


% _________________________________________________________________________
end


