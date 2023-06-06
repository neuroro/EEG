function eegLocalSpectra( channels )
%
% •.° Localised Event-Related Spectra °.•
% _________________________________________________________________________
%
% Extract event-related spectra at the specified electrode or average of an
% electrode cluster from all TimeFrequencyDataP<#><Condition>.mat files
% containing data centred on more than one event, which are located in the
% Current Folder and sub-folders of the Current Folder then save data in
% one TimeFrequency<Location>.mat file
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
%   Extract event-related spectra from a channel name or index, from the
%   average of a cluster of channels as a cell array of channel names or as
%   a vector of channel indices, or from all channels as (), '', or []
%   (optional input, default all)
%
% Output:
%   TimeFrequency<Location>.mat file containing, for each event centring,
%   pooled decompositions from all participants and conditions ready for
%   statistics and grand averages per condition ready for plotting
%
% • Requires •
% -------------------------------------------------------------------------
% Time-frequency decompositions of EEG data produced by eegTimeFrequency.m
% or eeg3TimeFreqyency.m in TimeFrequencyDataP<#><Condition>.mat files
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
disp( '•.° Localised event-related spectra °.•' )
disp( '_________________________________________________________________________' )
disp( ' ' )
disp( 'Extract event-related spectra at the specified electrode or average of an' )
disp( 'electrode cluster from all TimeFrequencyDataP<#><condition>.mat files'     )
disp( 'containing data centred on more than one event, which are located in the'  )
disp( 'Current Folder and sub-folders of the Current Folder then save data in'    )
disp( 'one TimeFrequency<Location>.mat file'                                      )
disp( ' ' )

codeRunTime( 'Started at' )
disp( ' ' )


%% Input
% -------------------------------------------------------------------------

if ~nargin || isempty( channels ) || ~channels
    channels = [];
end


%% Participant x condition files
% -------------------------------------------------------------------------

% Wildcard file name
aWildDataset     = 'TimeFrequencyData*.mat';

% Search for all datasets named in common located in the Current Folder and
% sub-folders of the Current Folder
fileStruct       = dir( [ '**/' aWildDataset ] );
nFiles           = length( fileStruct );
fileList         = { fileStruct(:).name   };
folderList       = { fileStruct(:).folder };

% Pre-allocate
participantList  = {};
participants     = [];
iParticipant     = 0;
conditionList    = {};
conditions       = {};
iCondition       = 0;

% % Participant code
% iTitleEnd          = length( 'TimeFrequencyData' );
% iX1                = regexp( fileList{1}, '[0-9]+' );
% iParticipantPrefix = iTitleEnd+1:iX1-1;
% participantPrefix  = fileList{1}(iParticipantPrefix);


% Build vectors of participant numbers and conditions
% -------------------------------------------------------------------------

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
        participants(iParticipant) = str2double( participantNumber );

    end

    % Add current file participant to the participant file list
    participantList{f} = participantNumber;


    % Build a vector of unique conditions
    % ---------------------------------------------------------------------

    % Check that this participant has not already been included
    if ~any( strcmp( conditionList, currentCondition ) )

        % Current index in the participant number list
        iCondition             = iCondition + 1;

        % Add current participant number to the list if not already included
        conditions{iCondition} = currentCondition;

    end

    % Add current file participant to the participant file list
    conditionList{f} = currentCondition;


end

% Sample size N participants
N           = length( participants );

% Number of conditions
nConditions = length( conditions );


%% Decomposition parameters
% -------------------------------------------------------------------------

Decomposition = load( fullfile( folderList{1}, fileList{1} ) );

% Event-centred windowings
eventCentres  = fieldnames( Decomposition );
nCentres      = length( eventCentres );

% Frequencies
for w = 1:nCentres
    LocalSpectra.(eventCentres{w}).Frequencies = Decomposition.(eventCentres{w}).Frequencies;
end
nFrequencies  = length( LocalSpectra.(eventCentres{1}).Frequencies );

% Times
for w = 1:nCentres
    LocalSpectra.(eventCentres{w}).Times = Decomposition.(eventCentres{w}).Times;
    nTimes.(eventCentres{w})             = length( LocalSpectra.(eventCentres{w}).Times ); % Time points vary per window centre
end

% Channels
iChannelSet = Decomposition.(eventCentres{1}).Channels;           % Channel indices of the channel co-ordinates that are in the decomposition
chanlocs    = Decomposition.(eventCentres{1}).ChannelCoordinates; % Channel co-ordinates
if ~isempty( channels )
    for ch = 1:length( channels )                                 % Loop through: channels to extract
        % Find the index of each selected channel number
        if isnumeric( channels )
            iCurrentChannel  = channels(ch);
        elseif iscell( channels )
            currentChannel   = channels{ch};
            iCurrentChannel  = find( strcmp( { chanlocs(:).labels }, currentChannel ) );
        end
        iChannels(ch) = find( iChannelSet == iCurrentChannel );
    end
else
    iChannels = 1:length( iChannelSet );
end
channelNames = { chanlocs(iChannels).labels };
for w = 1:nCentres
    LocalSpectra.(eventCentres{w}).ChannelCoordinates = chanlocs; % Store channel co-ordinates
end

% Display channels
if length( channels ) > 1                                     % Multiple selected
    channelsText = [];
    for ch = 1:length( channels )
        channelsText = [ channelsText channelNames{ch} ' ' ]; % Build list of channe names to display
    end
    disp( [ 'Averaging channels ' channelsText ] )
elseif length( channels ) == 1
    disp( [ 'Extracting channel ' channelNames{1} ] )         % One selected
else
    disp( 'Extracting all channels' )
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


%% Save
% -------------------------------------------------------------------------

% Location naming
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
if any( contains( { chanlocs(:).labels }, 'E' ) )
    Fz  = 'E11';
    FCz = 'E6';
    CPz = 'E55';
    Pz  = 'E62';
    POz = 'E72';
    Oz  = 'E75';
end

% Important channel indices
iFz  = find( strcmp( { chanlocs(:).labels }, Fz  ) );
iFCz = find( strcmp( { chanlocs(:).labels }, FCz ) );
iCPz = find( strcmp( { chanlocs(:).labels }, CPz ) );
iPz  = find( strcmp( { chanlocs(:).labels }, Pz  ) );
iPOz = find( strcmp( { chanlocs(:).labels }, POz ) );
iOz  = find( strcmp( { chanlocs(:).labels }, Oz  ) );
iP7  = find( strcmp( { chanlocs(:).labels }, P7  ) );
iP8  = find( strcmp( { chanlocs(:).labels }, P8  ) );

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

end % if number of selected electrodes

% Save
% -------------------------------------------------------------------------

% File name
fileName = [ 'TimeFrequency' locationName ];

% Save .mat 
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


