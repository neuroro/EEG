function ...
eegTimeFrequency( setName, frequencyLimits, frequencyResolution, blending )
%
% •.° Time-Frequency Decompositions of Event-Related EEG °.•
% _________________________________________________________________________
%
% Calculate event-related spectral power and inter-trial phase coherence of
% EEG data for each trial, centred on the stimulus presentation and the
% decision response, for all EEGLAB datasets that are named in common and
% located in the Current Folder and sub-folders of the Current Folder, then
% save a .mat file for each participant x condition in the Current Folder
%
% Trials are assumed to contain a stimulus presentation event labelled per
% condition and a response-decision event
%
% EEG data may be continuous or epoched but epochs must be wider than the
% trial by half a wavelength of the lowest frequency to be decomposed
%
% • Usage •
% -------------------------------------------------------------------------
% Edit the configuration of the EEG data in the code
% Set the Current Folder to the location of the datasets or the base folder
%
% >> eegTimeFrequency( setName, frequencyLimits, frequencyResolution, blending )
%
% For example:
% >> eegTimeFrequency( 'Memory' )
% >> eegTimeFrequency( 'Memory', [ 2 30 ], 30, 'sigmoid' )
% >> eegTimeFrequency( 'Memory', [], 0, 'median' )
%
%  <-- - -  -   -     -        -             -        -     -   -  - - -->
%
% •••( Function Inputs )
%
%   setName:             Part of the file name that is common to all the 
%                         EEGLAB datasets to be decomposed that are located 
%                         in the Current Folder and sub-folders of the
%                         Current Folder, for example 'Memory' or ''
%                         (optional input, default all datasets)
%
%   frequencyLimits:     Frequencies in Hz to extract as [minimum maximum]
%                         (optional input, default [2 30])
%
%   frequencyResolution: Frequencies per octave
%                         (optional input, default 30)
%
%   blending:            Blend the edges of variable-length trials using a 
%                         'sigmoid' taper to 0 dB centred around adjacent 
%                         events or a 'median' filter then Gaussian blur of
%                         the time points with trial dropout
%                         (optional input, default sigmoid blending)
%
% [ Function Output ] = 
%
%   TimeFrequencyDataP<number><condition>.mat files containing data centred
%   on the stimulus and data centred on the response for each participant
%   and each condition, from which a localised decomposition can then be
%   extracted using eegLocalSpectra.m
%
%  <-- - -  -   -     -        -             -        -     -   -  - - -->
%
% !!! Requires !!!
% 
% waveletTransform.m by Rohan King (2023) https://github.com/neuroro/EEG/blob/main/waveletTransform.m
% spectralBlender.m  by Rohan King (2023) https://github.com/neuroro/EEG/blob/main/spectralBlender.m
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


% Introduction
disp( ' ' )
disp( '•.° EEG Time-Frequency Decomposition °.•' )
disp( '_________________________________________________________________________' )
disp( ' ' )
disp( 'Calculate event-related spectral power and inter-trial phase coherence of' )
disp( 'EEG data for each trial, centred on the stimulus presentation and the'     )
disp( 'decision response, for all EEGLAB datasets that are named in common and'   )
disp( 'located in the Current Folder and sub-folders of the Current Folder, then' )
disp( 'save a .mat file for each participant x condition in the Current Folder'   )
disp( ' ' )


%% Configuration
% -------------------------------------------------------------------------

% Participant code prefix of the participant number in the file name
% For example 'P' 'S' 'Participant' 'sub-'
participantCode = 'P';

% Standard recognition task familiarity recode prefix
% participantCode = 'Familiarity';

% !!! INPUT YOUR PARTICPANT CODE PREFIX !!!


% Stimulus event labels for each condition
% -------------------------------------------------------------------------
% Trials are assumed to contain a stimulus presentation event labelled per
% condition and a response-decision event

% Generic task
C.Task        = { '1' };

% Sternberg memory task probes
C.Probes      = { 'STI5' 'STI7' 'STI9' };

% Standard recognition task
C.Recognition = { 'Familiar' 'Recognised' 'Misidentified'     ...
                  'KnownUnverified' 'UnfamiliarUnverified'    ...
                  'ControlTrueUnfamiliar' 'ControlFalseKnown' };

% Altruism reward task
C.Altruism    = { 'AccH' 'AccM' 'AccL' 'ShaH' 'ShaM' 'ShaL' };

% !!! INPUT YOUR CONDITIONS !!!
conditions    = C.Task;
% conditions    = C.Probes;
% conditions    = C.Recognition;
% conditions    = C.Altruism;


% Events in the trial other than the stimulus
% -------------------------------------------------------------------------

% Order in the trial of the response event relative to the stimulus event,
% which is assumed to be 1st, and is 1st after removal of preceding events
responseEventOrder = 2; % 2nd after stimulus

% Non-events, such as fixation, which should not be windowed around
NE.nonEvents = { 'Fixation' 'fixation' 'FXTN' 'empty' 'boundary' };

% Non-events in the Standard recognition task
NE.nonEventsRecognition = { 'S  1' 'empty' };

% Non-events in the Altruism task
boundaryEvents   = { 'b' };
numberCharacters = cell( 1, 10 );
for n = 0:9
    numberCharacters{n+1} = num2str(n);
end
NE.nonEventsAltruism = [ boundaryEvents numberCharacters ];

% !!! INPUT YOUR NON-EVENTS !!!
nonEvents          = NE.nonEvents;
% nonEvents          = NE.nonEventsRecognition;
% nonEvents          = NE.nonEventsAltruism;


% Time parameters
% -------------------------------------------------------------------------

% !!! INPUT YOUR TIMES !!!

% Sampling rate
samplingRate     = 1000;

% Trial limits
baselineLimits   = [ -500 0    ];   % Relative to stimulus time
stimulusLimits   = [ -200 1000 ];   % Relative to stimulus time
responseLimits   = [ -800 500  ];   % Relative to response time
maxTrialWindow   = 2000;

% Minimum reaction time for a trial to be valid
minReactionTime  = 100;

% Blending parameters
blendingDuration = 80;              % Sigmoid blend time points around adjacent event
neighbourhood    = [3 30];          % Median blend neighbourhood of frequencies x time points
sigmaGauss       = 2;               % Gaussian blur standard deviations if median blending
nTrialsMinimum   = 10;              % Time points with fewer trials than this number are zeroed if median blending


% EEG cap system
% -------------------------------------------------------------------------
% For example '10-10' '10-20' 'Indices' 'Brain Products' 'Biosemi' 'EGI'

capSystem = '10-10';


% Channels of interest
% -------------------------------------------------------------------------

% International 10-10 system
ChannelSets.Frontal1010      = { 'Fz'  'F3'  'F4' 'FCz' };
ChannelSets.Parietal1010     = { 'PO7' 'PO8' 'P7' 'P8'  };
ChannelSets.Occipital1010    = {};
ChannelSets.Temporal1010     = {};

% EGI system
ChannelSets.FrontalEGI       = { 'E24' 'E19' 'E11' 'E4' 'E124' 'E12' 'E5' 'E6' };
ChannelSets.ParietalEGI      = { 'E52' 'E60' 'E61' 'E62' 'E78' 'E85' 'E92'     };
ChannelSets.OccipitalEGI     = {};
ChannelSets.TemporalEGI      = {};

% Indices in the EEG data
ChannelSets.FrontalIndices   = [ 24 19 11 4  124 5  12 6 ];
ChannelSets.ParietalIndices  = [ 52 60 61 62 78  85 92   ];
ChannelSets.OccipitalIndices = [];
ChannelSets.TemporalIndices  = [];


% !!! INPUT YOUR CAP SYSTEM AND CHANNELS !!!


%% Default inputs
% -------------------------------------------------------------------------

% Default frequency range limits
if nargin < 2 || isempty( frequencyLimits )
    frequencyLimits = [ 2 30 ]; % Minium and maximum (Hz)
end

% Defaut frequencies per octave
if nargin < 3 || isempty( frequencyResolution ) || ~frequencyResolution
    frequencyResolution = 30;
end

% Default blending method
if nargin < 4
    blending = 'sigmoid';
end


%% Derived parameters
% -------------------------------------------------------------------------

% Times
longestCycle     = samplingRate/frequencyLimits(1);
edge             = ceil( 1.5 * longestCycle );              % Edges are excised after decomposition to remove edge effects
if isempty( baselineLimits ) || baselineLimits(1) > -200
    startTime    = -200;                                    % Start no later than -200 ms so that pre-stimulus data exists
else
    startTime    = baselineLimits(1);
end
startTime        = startTime      - edge;
startSeconds     = startTime      / 1000;
endTime          = maxTrialWindow + edge;
endSeconds       = (endTime + 1)  / 1000;

% No blending for median filtering
if contains( blending, 'm', 'IgnoreCase', true ) && ~contains( blending, 's', 'IgnoreCase', true )
    blendingDuration = Inf;
end

% Blend adjustment for trials with no response 
if ~isinf( blendingDuration )
    adjustment   = blendingDuration/2;
else
    adjustment   = 0;
end

% Frequencies
octaves          = log2( frequencyLimits(2)/frequencyLimits(1) );
nFrequencies     = ceil( frequencyResolution*octaves );

% Channels
internationalSystemCaps = { 'International' '10' '20' 'Brain' 'Biosemi' 'Acti' };
if contains( capSystem, internationalSystemCaps, 'IgnoreCase', true )
    capSystem    = '1010';
elseif contains( capSystem, { 'EGI' 'Net' 'Station' 'Philips' }, 'IgnoreCase', true )
    capSystem    = 'EGI';
else
    capSystem    = 'Indices';
end
channelSetFields = fieldnames( ChannelSets );
iFieldsSystem    = contains( channelSetFields, capSystem );
nonSystemFields  = channelSetFields(~iFieldsSystem);
ChannelSets      = rmfield( ChannelSets, nonSystemFields );
channelSets      = struct2cell( ChannelSets );
channelSet       = [ channelSets{:} ];
nChannels        = length( channelSet );

% Conditions
if ischar( conditions )
    conditions   = { conditions };
end
nConditions      = length( conditions );


%% EEGLAB dataset files
% -------------------------------------------------------------------------

% Check input
aWildError = [ 'Input a name in common to all EEGLAB datasets located ' ...
               'in the Current Folder and sub-folders of the Current '  ...
               'Folder or '''' for all datasets'                        ];
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
elseif contains( setName, '.set' )
    aWildDataset = [ '*' setName ];
else
    aWildDataset = [ '*' setName '*.set' ];
end

% Search for all datasets named in common located in the current folder and
% sub-folders of the current folder
fileStruct = dir( [ '**/' aWildDataset ] );
nFiles     = length( fileStruct );
fileList   = { fileStruct(:).name   };
folderList = { fileStruct(:).folder };


%% Display parameters
% -------------------------------------------------------------------------

disp( '• Times •' )
disp( '-------------------------------------------------------------------------' )
disp( [ num2str( samplingRate ) ' Hz sampling rate' ] )
disp( [ num2str( minReactionTime ) ' ms minimum reaction time' ] )
if ~isinf( blendingDuration )
    disp( [ num2str( blendingDuration ) ' ms sigmoid blend time for variable trial durations'] )
else
    disp( [ num2str( neighbourhood(2) ) ' ms ' num2str( neighbourhood(1) ) ' Hz median blend neighbourhood for variable trial durations'] )
end
disp( 'Trial:')
disp( [ '0 ms to ' num2str( maxTrialWindow ) ' ms trial window' ] )
if ~isempty( baselineLimits )
    disp( [ num2str( baselineLimits(1) ) ' ms to ' num2str( baselineLimits(2) ) ...
            ' ms baseline window' ] )
else
    disp( 'No baseline correction' )
end
disp( [ num2str( startTime ) ' ms to ' num2str( endTime ) ...
        ' ms decomposition window' ] )
disp( 'Event-locked windows:')
disp( [ num2str( stimulusLimits(1)   ) ' ms to ' num2str( stimulusLimits(2) ) ...
        ' ms relative to the stimulus' ] )
disp( [ num2str( responseLimits(1) ) ' ms to ' num2str( responseLimits(2) ) ...
        ' ms relative to the response' ] )
disp( ' ' )
disp( '• Frequencies •' )
disp( '-------------------------------------------------------------------------' )
disp( [ num2str( frequencyLimits(1) ) ' - ' num2str( frequencyLimits(2) ) ' Hz' ] )
disp( 'logarithmically spaced in octaves' )
disp( [ num2str( frequencyResolution ) ' frequencies per octave' ] )
disp( [ num2str( octaves ) ' octaves' ] )
disp( [ num2str( nFrequencies ) ' frequencies' ] )
disp( ' ' )
disp( '• Channels •' )
disp( '-------------------------------------------------------------------------' )
ChannelList = [];
divisions   = ceil( nChannels / 10 );
for div = 1:divisions
    ChannelList.(['C' num2str(div)]) = [];
end
for ch = 1:nChannels
    if isnumeric( channelSet )
        channelName = num2str( channelSet(ch) );
    else
        channelName = char( channelSet{ch} );
    end
    currentDivision = num2str( ceil( ch/10 ) );
    ChannelList.(['C' currentDivision]) = [ ChannelList.(['C' currentDivision]) channelName ' ' ];
end
for div = 1:divisions
    disp( ChannelList.(['C' num2str(div)]) )
end
disp( ' ' )


%% Decomposition
% -------------------------------------------------------------------------

disp( '• Decomposition •')
disp( '-------------------------------------------------------------------------' )
disp( [ num2str( nFiles ) ' participant data sets' ] )

% Start time
timeFrequencyRunTime( 'Started at' )
disp( ' ' )

% Loop through: Files
for n = 1:nFiles

    % File
    currentFile    = fileList{n};
    currentFolder  = folderList{n};
    disp( [ 'Decomposing ' currentFile ] )
    disp( ' ' )

    % Load EEG dataset into an EEGLAB struct
    loaderargin    = { 'filename', currentFile,   ...
                      'filepath', currentFolder, ...
                      'verbose',  'off'          };
    CurrentEEG     = pop_loadset( loaderargin{:} );
    CurrentEEG     = eeg_checkset( CurrentEEG );

    % Sanity check: Sampling rate
    if CurrentEEG.srate ~= samplingRate
        error( [ 'Sampling rate in file not equal to ' num2str( samplingRate ) ] )
    end

    % Channel indices for named channels input
    if iscell( channelSet )
        % Compare the channel names to the channel co-ordinates labels in
        % the chanlocs.labels field
        channelSet = eeg_chaninds( CurrentEEG, channelSet );
    end

    % Loop through: Conditions
    parfor c = 1:nConditions
    % for c = 1:nConditions

        % Current condition event label
        condition     = conditions{c};

        % Parallel broadcast variables
        EEG           = CurrentEEG;             % Struct copied to each parallel worker as pop_loadeset is slower
        channels      = channelSet;             % Always indices at this point
        baselineLimit = baselineLimits;         % In ms
        stimulusLimit = stimulusLimits;         % In ms
        responseLimit = responseLimits;         % In ms
        blendDuration = blendingDuration;       % In ms
        responseEvent = responseEventOrder - 1; % Events after the stimulus in the trial

        % Pre-allocate
        Decomposition = [];
        trialCentres  = { 'Stimulus' 'Response' };
        nTrialCentres = length( trialCentres );
        for tc = 1:nTrialCentres
            Decomposition.(trialCentres{tc}).SpectralPower  = [];
            Decomposition.(trialCentres{tc}).PhaseDirection = [];
            Decomposition.(trialCentres{tc}).PhaseCoherence = [];
            Decomposition.(trialCentres{tc}).Coefficients   = [];
        end


        %% Trial data
        % -----------------------------------------------------------------

        % Deal with conditions in which there is only one trial or none at all
        eegEvents      = { EEG.event(:).type };
        iStimulusEvent = strcmpi( condition, eegEvents );
        if sum( iStimulusEvent ) == 1
            singleTrial          = true;
            disp( [ '[' 8 'Warning: Single trial for ' condition ']' 8 ] )
            iStimulus            = find( iStimulusEvent );
            stimulusLatency      = EEG.event(iStimulus).latency;
            % Duplicate single trial stimulus event so EEGLAB produces epochs
            EEG                  = pop_editeventvals( EEG, 'insert',      { iStimulus, [], [], [], [] },        ...
                                                           'changefield', { iStimulus, 'type',     condition }, ...
                                                           'changefield', { iStimulus, 'duration', 0 },         ...
                                                           'changefield', { iStimulus, 'latency',  0 });
            EEG.event(1).latency = stimulusLatency;
            EEG                  = eeg_checkset( EEG );
            EEG                  = pop_editeventvals( EEG, 'sort', { 'latency', 0 } );
            EEG                  = eeg_checkset( EEG );
            responseEvent        = responseEvent + 1; % Correction for duplicate stimulus events
        elseif ~any( iStimulusEvent )
            disp( [ '[' 8 'Warning: No trials for ' condition ']' 8 ] )
            continue
        else
            singleTrial          = false;
        end
        
        % Epoch trials
        EEG     = pop_epoch( EEG, { condition }, [startSeconds endSeconds], 'epochinfo', 'yes' );
        EEG     = eeg_checkset( EEG );
        nTrials = length( EEG.epoch );

        % Sanity check: Epoch vs. decomposition window
        if EEG.pnts ~= length( startTime:endTime )
            disp( [ num2str( length( startTime:endTime ) ) ' decomposition window time points' ] )
            disp( [ num2str( EEG.pnts ) ' epoch time points'] )
            error( 'Epoching produced the wrong number of time points' )
        end

        % Deal with conditions in which the first trials are skipped
        skipStart = true; % Set to false when a realistic trial is found
        skipCount = 0;    % Count of unrealistic trials at the start

        % Loop through: Trials
        for trial = 1:nTrials

            % Trial events
            currentTrialEvents  = EEG.epoch(trial).eventtype(:);
            nCurrentTrialEvents = length( currentTrialEvents );

            % Find the index of the first stimulus event in case any
            % preceding events not caught by nonEvents or subsequent
            % stimulus events made it into the epoch
            iEpochEventStimulus = find( matches( currentTrialEvents, condition ), 1 );

            % Index of the response event using the trial event order
            iEpochEventResponse = iEpochEventStimulus + responseEvent;

            % Non-events to exclude
            iPreEvents          = 1:nCurrentTrialEvents < iEpochEventStimulus; % Indices of events before the stimulus
            iNonEvents          = matches( currentTrialEvents, nonEvents );    % Indices of non-events such as fixation
            iExclude            = or( iPreEvents', iNonEvents );               % Exclude either (transpose was needed to keep it a vector)

            % Trial event-locked centerings
            currentTrialCentres = currentTrialEvents(~iExclude);
            nCentres            = length( currentTrialCentres );

            % Trial event times (in ms) relative to 0 ms
            stimulusTime     = 0;
            if nCentres > 1    % Stimulus and response in the epoch
                responseTime = EEG.epoch(trial).eventlatency{iEpochEventResponse};
            else               % Stimulus only, no response
                responseTime = EEG.times(end) - adjustment - 1; % Blend out by the end of the trial
            end
            centreTimes      = [ stimulusTime responseTime ];

            % Process trials with realistic reaction time (in ms)
            if responseTime >= minReactionTime

                % A realistic trial has been found
                skipStart = false;

                % Loop through: Channels
                for ch = 1:nChannels

                    % Current channel
                    channel = channels(ch);


                    %% Time-frequency decomposition
                    % -----------------------------------------------------

                    % EEG signal
                    signal = EEG.data(channel,:,trial);

                    % Continuous Wavelet Transform
                    %   See waveletTransform.m for details
                    [ spectralPower,  ...
                      phaseDirection, ...
                      coefficients,   ...
                      frequencies,    ...
                      signalIndices ] = waveletTransform( signal,                ...
                                                          samplingRate,          ...
                                                          frequencyLimits,       ...
                                                          'FrequencyResolution', ...
                                                          frequencyResolution    );

                    % Baseline mean spectrum
                    baselineLength = length( baselineLimit(1):baselineLimit(2) );
                    baselineWindow = 1:baselineLength;
                    baseline       = spectralPower(:,baselineWindow);
                    baseline       = mean( baseline, 2, 'omitnan' );

                    % No baseline
                    if ~baselineLimit(1) || isempty( baselineLimit ) || isempty( baseline )
                        baseline   = ones( nFrequencies, 1 );
                    end

                    % Convert baseline power to decibel volts^2
                    baselinePower  = 10*log10( baseline );

                    % Baseline spectrum copied across the time window
                    baseline       = repmat( baselinePower, 1, length( spectralPower ) );

                    % Sanity check: Matching size
                    if size( baseline, 1 ) ~= size( spectralPower, 1 ) || size( baseline, 2 ) ~= size( spectralPower, 2 )
                        error( 'Baseline size mismatch' )
                    end

                    % Convert spectral power to decibel volts^2
                    spectralPower  = 10*log10( spectralPower );

                    % Subtract the baseline mean spectrum in logarithmic
                    % units, giving power in decibels relative to baseline,
                    % or retaining power in decibel volts^2 if no baseline
                    spectralPower  = spectralPower - baseline;

                    % Power units
                    if baselineLimit(1)
                        powerUnits = 'Decibels relative to baseline (mean spectrum)';
                    else
                        powerUnits = 'Decibel volts^2';
                    end

                    % Times (in ms)
                    %   Input signal indices -> output time points in ms
                    timePoints     = EEG.times(signalIndices);

                    % Sanity check: Time points
                    if any( timePoints ~= startTime + signalIndices - 1 )
                        disp( [ '[' 8 'Warning: Time points calculation ' ...
                                      'is inconsistent for ' currentFile  ']' 8 ] )
                    end


                    %% Extract event-locked data windows
                    % -----------------------------------------------------

                    % Time windows (in ms)
                    %   Relative time limits adjusted by centre event times
                    %   Limit relative to centre -> limit relative to 0
                    centreLimits  = { (stimulusTime + stimulusLimit) ...
                                      (responseTime + responseLimit) };

                    % Time domains relative to centre
                    centreDomains = { stimulusLimit(1):stimulusLimit(2) ...
                                      responseLimit(1):responseLimit(2) };

                    % Blender times
                    blenderTimes  = { responseTime ...
                                      stimulusTime };

                    % Loop through: Event-locked centres
                    for centre = 1:nTrialCentres % centre = 1:nTrialCentres

                        % Current centre
                        trialCentre  = trialCentres{centre};
                        centreLimit  = centreLimits{centre};
                        centreDomain = centreDomains{centre};
                        blenderTime  = blenderTimes{centre};
                        if centre == 1
                            blenderDuration = -blendDuration; % Blend out at the response
                        else
                            blenderDuration = blendDuration;  % Blend in at the stimulus
                        end

                        % Pre-allocate open windows
                        windowLength = length( centreDomain );
                        powerWindow  = NaN( nFrequencies, windowLength );
                        phaseWindow  = NaN( nFrequencies, windowLength )*1i;
                        wavesWindow  = NaN( nFrequencies, windowLength )*1i;

                        % Process event-centerings found within the epoch
                        %   Skip too-slow reactions (not within the epoch)
                        %   leaving them as NaN
                        if centre <= nCentres

                            
                            %% Event-locked data
                            % ---------------------------------------------

                            % Time-frequency blender
                            %   See spectralBlender.m for details
                            [ blendedPower, ...
                              blendedPhase, ...
                              blendedWaves ] = spectralBlender( spectralPower,  ...
                                                                phaseDirection, ...
                                                                coefficients,   ...
                                                                frequencies,    ...
                                                                timePoints,     ...
                                                                blenderTime,    ...
                                                                blenderDuration );

                            % Centre window
                            %   Find the time point for each centre limit
                            iCentreLimits1 = find( ~( timePoints - centreLimit(1) ) );
                            iCentreLimits2 = find( ~( timePoints - centreLimit(2) ) );
                            if isempty( iCentreLimits2 )
                                iCentreLimits2 = length( timePoints );
                            end
                            centreWindow = iCentreLimits1:iCentreLimits2;
                            centreLength = length( centreWindow );

                            % Centre the decomposition
                            powerWindow(:,1:centreLength) = blendedPower(:,centreWindow);
                            phaseWindow(:,1:centreLength) = blendedPhase(:,centreWindow);
                            wavesWindow(:,1:centreLength) = blendedWaves(:,centreWindow);


                        end

                        % Sanity checks
                        if isscalar( powerWindow ) || isempty( powerWindow )
                            disp( [ '[' 8 'Warning: No decomposition for trial ' num2str( trial ) ' centred on the ' trialCentre ' in ' currentFile ']' 8 ] )
                        elseif all( ~powerWindow( ~isnan( powerWindow ) ), 'all' )
                            disp( [ '[' 8 'Warning: Decomposition of zeroes for trial ' num2str( trial ) ' centred on the ' trialCentre ' in ' currentFile ']' 8 ] )
                        end

                        % Store in struct
                        Decomposition.(trialCentre).SpectralPower(trial,ch,:,:)  = powerWindow;
                        Decomposition.(trialCentre).SpectralPowerUnits           = powerUnits;
                        Decomposition.(trialCentre).SpectralPowerDimensions      = 'Channels x Frequencies x Times'; % As it will be after averaging
                        Decomposition.(trialCentre).PhaseDirection(trial,ch,:,:) = phaseWindow;
                        Decomposition.(trialCentre).Coefficients(trial,ch,:,:)   = wavesWindow;
                        Decomposition.(trialCentre).CoefficientsDimensions       = 'Trials x Channels x Frequencies x Times';
                        Decomposition.(trialCentre).Frequencies                  = frequencies;
                        Decomposition.(trialCentre).FrequenciesUnits             = 'Hz';
                        Decomposition.(trialCentre).Times                        = centreDomain;
                        Decomposition.(trialCentre).TimesUnits                   = 'milliseconds';
                        Decomposition.(trialCentre).BaselinePowerSpectrum        = baselinePower;
                        Decomposition.(trialCentre).BaselinePowerSpectrumUnits   = 'Decibel volts^2';
                        Decomposition.(trialCentre).BaselineTimeLimits           = baselineLimit;
                        Decomposition.(trialCentre).ChannelIndices               = channels;
                        Decomposition.(trialCentre).ChannelCoordinates           = EEG.chanlocs;
                        Decomposition.(trialCentre).BlendTimes(trial,:)          = blenderTime - centreTimes(centre);


                    end % for Centres

                end % for Channels

            % Skip trials with unrealistically fast initiation time
            else

                % If the first trials are unrealistically fast set to NaN
                % after dimension sizes are realised by realistic trials
                if skipStart
                    skipCount = skipCount + 1;
                    continue
                end

                for centre = 1:nTrialCentres

                    trialCentre = trialCentres{centre};

                    Decomposition.(trialCentre).SpectralPower(trial,:,:,:)  = NaN;
                    Decomposition.(trialCentre).PhaseDirection(trial,:,:,:) = NaN*1i;
                    Decomposition.(trialCentre).Coefficients(trial,:,:,:)   = NaN*1i;

                end

            end % if Test initiation time

        end % for Trials

        % Set first trials to NaN if they were skipped
        if skipStart
            for centre = 1:nTrialCentres
            
               trialCentre = trialCentres{centre};
        
               Decomposition.(trialCentre).SpectralPower(1:skipCount,:,:,:)  = NaN;
               Decomposition.(trialCentre).PhaseDirection(1:skipCount,:,:,:) = NaN*1i;
               Decomposition.(trialCentre).Coefficients(1:skipCount,:,:,:)   = NaN*1i;
    
            end
        end

        % Store channel co-ordinates
        for centre = 1:nTrialCentres
            Decomposition.(trialCentres{centre}).ChannelCoordinates = EEG.chanlocs;
        end


        %% Power, phase coherence, and median blending
        % -----------------------------------------------------------------

        % Loop through: Event-locked centres
        for centre = 1:nTrialCentres

            trialCentre = trialCentres{centre};

            % Time-points with too-few trials for median blending
            %   Blending with a median filter won't prevent single extreme
            %   trials from skewing the average
            if isinf( blendDuration )
                missingPoints  = isnan( Decomposition.(trialCentre).SpectralPower );
                missingPoints  = squeeze( missingPoints(:,1,1,:) ); % All channels and frequencies
                trialsPerPoint = sum( ~missingPoints, 1 );          %    have an equal trial count
                iTooFew        = find( trialsPerPoint < nTrialsMinimum );
                if singleTrial
                    iTooFew    = find( trialsPerPoint < 1 );
                end
            end


            % Induced and evoked oscillations
            % -------------------------------------------------------------

            % Mean power across trials
            Decomposition.(trialCentre).SpectralPower ...
                = mean( Decomposition.(trialCentre).SpectralPower, 1, 'omitnan' );

            % Inter-trial phase coherence
            Decomposition.(trialCentre).PhaseCoherence ...
                = abs( mean( Decomposition.(trialCentre).PhaseDirection, 1, 'omitnan' ) );
            Decomposition.(trialCentre).PhaseCoherenceUnits ...
                = 'Phase alignment proportion';
            Decomposition.(trialCentre).PhaseCoherenceDimensions ...
                = 'Channels x Frequencies x Times';

            % Clean up
            Decomposition.(trialCentre).SpectralPower ...
                = squeeze( Decomposition.(trialCentre).SpectralPower  );
            Decomposition.(trialCentre).PhaseCoherence ...
                = squeeze( Decomposition.(trialCentre).PhaseCoherence );
            Decomposition.(trialCentre) ...
                = rmfield( Decomposition.(trialCentre), 'PhaseDirection' ); % Remove phase direction
            if singleTrial
            Decomposition.(trialCentre).Coefficients ...
                = Decomposition.(trialCentre).Coefficients(1,:,:,:);        % Remove duplicate trial from coefficients
            end

            % Sanity checks
            if isscalar( Decomposition.(trialCentre).SpectralPower ) || isempty( Decomposition.(trialCentre).SpectralPower )
                disp( [ '[' 8 'Warning: No decomposition centred on the ' trialCentre ' for ' currentFile ']' 8 ] )
            elseif all( not( Decomposition.(trialCentre).SpectralPower( ~isnan( Decomposition.(trialCentre).SpectralPower ) ) ), 'all' )
                disp( [ '[' 8 'Warning: Decomposition of zeroes centred on the ' trialCentre ' for ' currentFile ']' 8 ] )
            end


            % Median filter trial edges to smooth out discontinuities
            % -------------------------------------------------------------

            if isinf( blendDuration )
                times  = Decomposition.(trialCentre).Times;
                nTimes = length( times );


                % Variable trial length time points to filter
                % ---------------------------------------------------------
                % From the trial cut point closest to the centre time
                % All time points with fewer trials are filtered

                % Edges after the centre event to the end
                %   Index range of variable response times relative to stimulus
                blendTimes  = Decomposition.(trialCentre).BlendTimes;
                if size( blendTimes, 2 ) == 2
                    blendTimes = blendTimes(:,2);
                end
                earliestCut = min( abs( blendTimes ), [], 'omitnan' );
                iEdgeA1     = find( times == earliestCut );
                iEdgeA2     = nTimes;
                iEdgeA      = iEdgeA1:iEdgeA2;

                % Edges before the centre event from the start
                %   Index range of variable stimulus times relative to response
                blendTimes  = Decomposition.(trialCentre).BlendTimes;
                if size( blendTimes, 2 ) == 2
                    blendTimes = blendTimes(:,1);
                end
                iEdgeB1     = 1;
                latestCut   = min( abs( blendTimes ), [], 'omitnan' );
                iEdgeB2     = find( times == -latestCut ); % Negative relative to the centre event
                iEdgeB      = iEdgeB1:iEdgeB2;

                % Pad the edges by the neighbourhood
                padding     = neighbourhood(2) * 1000 / samplingRate;       %#ok
                iEdgeAA     = (iEdgeA1 - padding):iEdgeA2;
                iEdgeBB     = iEdgeB1:(iEdgeB2 + padding);


                % Median filter
                % ---------------------------------------------------------
                
                % Edges to filter
                switch centre
                    case 1 % Stimulus-locked
                        before = false;
                        after  = true;
                    case 2 % Response-locked
                        before = true;
                        after  = false;
                end

                % Loop through: Channels
                for ch = 1:nChannels

                    % Current channel spectra
                    channelPower = Decomposition.(trialCentre).SpectralPower(ch,:,:);
                    channelPhase = Decomposition.(trialCentre).PhaseCoherence(ch,:,:);
                    channelPower = squeeze( channelPower );
                    channelPhase = squeeze( channelPhase );

                    % Edges after the centre event to the end
                    if after

                        % Spectral power
                        edgeSpectra = channelPower(:,iEdgeAA);
                        edgeSpectra = medfilt2( edgeSpectra, neighbourhood );
                        edgeSpectra = imgaussfilt( edgeSpectra, sigmaGauss );
                        Decomposition.(trialCentre).SpectralPower(ch,:,iEdgeA) ...
                                    = edgeSpectra(:,1+padding:end);

                        % Phase coherence
                        edgeSpectra = channelPhase(:,iEdgeAA);
                        edgeSpectra = medfilt2( edgeSpectra, neighbourhood );
                        edgeSpectra = imgaussfilt( edgeSpectra, sigmaGauss );
                        Decomposition.(trialCentre).PhaseCoherence(ch,:,iEdgeA) ...
                                    = edgeSpectra(:,1+padding:end);

                    end

                    % Edges before the centre event to the end
                    if before

                        % Spectral power
                        edgeSpectra = channelPower(:,iEdgeBB);
                        edgeSpectra = medfilt2( edgeSpectra, neighbourhood );
                        edgeSpectra = imgaussfilt( edgeSpectra, sigmaGauss );
                        Decomposition.(trialCentre).SpectralPower(ch,:,iEdgeB) ...
                                    = edgeSpectra(:,1:end-padding);

                        % Phase coherence
                        edgeSpectra = channelPhase(:,iEdgeBB);
                        edgeSpectra = medfilt2( edgeSpectra, neighbourhood );
                        edgeSpectra = imgaussfilt( edgeSpectra, sigmaGauss );
                        Decomposition.(trialCentre).PhaseCoherence(ch,:,iEdgeB) ...
                                    = edgeSpectra(:,1:end-padding);

                    end

                end % for Channels
                

                % Erase time-points with too-few trials
                Decomposition.(trialCentre).SpectralPower(:,:,iTooFew)  = 0;
                Decomposition.(trialCentre).PhaseCoherence(:,:,iTooFew) = 0;
                Decomposition.(trialCentre).Coefficients(:,:,iTooFew)   = NaN*1i;


            end % if Median filter

        end % for Centres


        %% Save single participant x condition
        % -----------------------------------------------------------------

        % File name
        [ iP, iXn ] = regexp( currentFile, [ participantCode '[0-9]+' ] );
        participant = currentFile(iP:iXn);
        fileName    = [ 'TimeFrequencyData' participant condition ];

        % Save .mat file
        timeFrequencyParallelSave( fileName, Decomposition )

        % Single completion time
        timeText    = [ currentFile ' condition ' condition ' decomposition completed at' ];
        timeFrequencyRunTime( timeText )

        
    end % parfor Conditions

end % for Files


% Finish time
timeFrequencyRunTime( 'Time-frequency data generation completed at' )


% _________________________________________________________________________
end



%%
% •.° Decomposition Run Time °.•
% _________________________________________________________________________
function timeFrequencyRunTime( message )

% Time
theTimeIs = datetime( 'now' );
theTimeIs = char( theTimeIs );

% Print
if exist( 'message', 'var' )
    disp( [ message ' ' theTimeIs ] )
else
    disp( theTimeIs )
end


% _________________________________________________________________________
end



%%
% •.° Parallel Save °.•
% _________________________________________________________________________
function timeFrequencyParallelSave( fileName, structOfVariablesToSave )

save( fileName, '-struct', 'structOfVariablesToSave', '-v7.3' )


% _________________________________________________________________________
end


