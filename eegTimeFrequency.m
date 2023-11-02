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
% save a .mat file for each participant x condition in original folder
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
% >> eegTimeFrequency( 'Memory', [], 0, 'smooth' )
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
%   blending:            Blend the edges of variable-length trials
%                         'sigmoid' taper to 0 dB centred around adjacent 
%                          events
%                         'smooth' using Savitsky-Golay filtering of the 
%                          time-points then Gaussian blur with a minimum
%                          trial cut-off
%                         'median' filter in 2-D then Gaussian blur with a
%                          minimum trial cut-off
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
% waveletTransform.m by Rohan King (2023) available at
%  https://github.com/neuroro/EEG/blob/main/waveletTransform.m
% spectralBlender.m  by Rohan King (2023) available at
%  https://github.com/neuroro/EEG/blob/main/spectralBlender.m
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
% King, R. (2023). eegTimeFrequency [MATLAB code]. GitHub. 
%  https://github.com/neuroro/EEG/eegTimeFrequency.m


% Introduction
disp( ' ' )
disp( '•.° EEG Time-Frequency Decomposition °.•' )
disp( '_________________________________________________________________________' )
disp( ' ' )
disp( 'Calculate event-related spectral power and inter-trial phase coherence of' )
disp( 'EEG data for each trial, centred on the stimulus presentation and the'     )
disp( 'decision response, for all EEGLAB datasets that are named in common and'   )
disp( 'located in the Current Folder and sub-folders of the Current Folder, then' )
disp( 'save a .mat file for each participant x condition in the original folder'  )
disp( ' ' )


%% Configuration
% -------------------------------------------------------------------------

% Participant code prefix of the participant number in the file name
% For example 'P' 'S' 'Participant' 'sub-'
participantCode = 'P';

% Standard recognition task familiarity recode prefix
% participantCode = 'Familiarity';

% !!! INPUT YOUR PARTICPANT CODE PREFIX !!!


% Event-related windows in order
EventRelatedWindows.Window1 = 'Stimulus';
EventRelatedWindows.Window2 = 'Response';


% Stimulus event labels for each condition
% -------------------------------------------------------------------------
% Trials are assumed to contain a stimulus presentation event labelled per
% condition that are extracted separately and a response-decision event

% Generic task
ConditionEvents.Task        = { 'Stimulus' };

% Sternberg memory task probes
ConditionEvents.Probes      = { 'STI5' 'STI7' 'STI9' };

% Standard recognition task
ConditionEvents.Recognition = { 'Familiar' 'Recognised' 'Misidentified'     ...
                                'KnownUnverified' 'UnfamiliarUnverified'    ...
                                'ControlTrueUnfamiliar' 'ControlFalseKnown' };

% Altruism reward task
ConditionEvents.Altruism    = { 'AccH' 'AccM' 'AccL' 'ShaH' 'ShaM' 'ShaL' };

% !!! INPUT YOUR CONDITIONS !!!
% conditions                  = ConditionEvents.Task;
conditions                  = ConditionEvents.Probes;
% conditions                  = ConditionEvents.Recognition;
% conditions                  = ConditionEvents.Altruism;


% Resonse event labels
% -------------------------------------------------------------------------
% Trials are assumed to contain a stimulus presentation event and a
% response-decision event that may be labelled per response type but are
% not extracted separately

% Generic task
ResponseEvents.Task        = { 'Response' };

% Sternberg memory task responses
ResponseEvents.Probes      = { 'CIL5' 'CIL7' 'CIL9' 'COL5' 'COL7' 'COL9' ...
                               'IIL5' 'IIL7' 'IIL9' 'IOL5' 'IOL7' 'IOL9' };

% Standard recognition task
ResponseEvents.Recognition = { 'S  4' 'S  8' };

% Altruism reward task
ResponseEvents.Altruism    = { 'Resp' };

% !!! INPUT YOUR RESPONSE TYPES !!!
% responses                  = ResponseEvents.Task;
responses                  = ResponseEvents.Probes;
% responses                  = ResponseEvents.Recognition;
% responses                  = ResponseEvents.Altruism;


% Time parameters
% -------------------------------------------------------------------------

% !!! INPUT YOUR TIMES !!!

% Sampling rate
samplingRate     = 1000;

% Trial limits
baselineLimits   = [ -500 0    ];   % Relative to stimulus time
stimulusLimits   = [ -200 1000 ];   % Relative to stimulus time
responseLimits   = [ -600 500  ];   % Relative to response time

% Maximum response time for a trial to be valid
maxResponseTime  = 1500;

% Minimum reaction time for a trial to be valid
minReactionTime  = 75;

% Blending parameters
blendingDuration = 80;              % Sigmoid blend time points around adjacent event (in ms)
neighbourhood    = [3 30];          % Median blend neighbourhood of frequencies x time points
order            = 5;               % Smoothing: Savitsky-Golay filter order
                                    %   Higher  is a closer fit -> less smoothing
frames           = 85;              % Smoothing: Savitsky-Golay filter frame length (must be odd)
                                    %   Shorter is a closer fit -> less smoothing
                                    %   Longer  is a wider  fit -> smoother
sigmaGauss       = 0.5;             % Gaussian blur standard deviations
nTrialsMinimum   = 10;              % Time points with fewer trials than this number are zeroed (if median blending or smoothing)


% EEG cap system
% -------------------------------------------------------------------------
% For example '10-10' '10-20' 'Indices' 'Brain Products' 'Biosemi' 'EGI'

capSystem = 'EGI'; % 'International';


% Channels of interest
% -------------------------------------------------------------------------

% International 10-10 system
ChannelSets.Frontal1010      = { 'Fz'  'F3'  'F4' };
ChannelSets.Parietal1010     = { 'PO7' 'PO8' 'P7' 'P8' };
ChannelSets.Occipital1010    = {};
ChannelSets.Temporal1010     = {};

% EGI system
ChannelSets.FrontalEGI       = { 'E11'  'E19' 'E4'  'E24' 'E124' , 'E16'  'E18' 'E10'  'E23' 'E3' , 'E12' 'E5' , 'E6'  'E13' 'E112' };  % { 'E11' 'E24' 'E124' };
ChannelSets.ParietalEGI      = {}; % { 'E62'  'E61' 'E78'  'E60' 'E85'  'E52' 'E92' , 'E65' 'E90'  'E59' 'E91'  'E58' 'E96' };                % { 'E65' 'E90' 'E58' 'E96' };
ChannelSets.OccipitalEGI     = {};
ChannelSets.TemporalEGI      = {};

% Indices in the EEG data
ChannelSets.FrontalIndices   = [ 11 24 124 ];
ChannelSets.ParietalIndices  = [ 65 90 58 96 ];
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

% Events that trial windows are centred on
trialCentres     = struct2cell( EventRelatedWindows );

% Edge time
longestCycle     = samplingRate/frequencyLimits(1);
edge             = ceil( 1.5 * longestCycle );              % Edges are excised after decomposition to remove edge effects

% Baseline time limits
if   ~isempty( baselineLimits ) && ...
   ( isscalar( baselineLimits ) || baselineLimits(2) == 0 )
    baselineLimits = [ -abs( baselineLimits ) -1 ];         % Baseline up to the stimulus
end

% Start time
if   isempty( baselineLimits ) || ...
   ( baselineLimits(1) > -200 && baselineLimits(1) <= 0 )
    startPoint   = -200;                                    % Start no later than -200 ms so that pre-stimulus data exists
else
    startPoint   = baselineLimits(1);                       % Start at least as early as the baseline lower limit
end
earliestReaction = minReactionTime + responseLimits(1);     % Earliest possible reaction time is the response-centred window lower limit before the minimum reaction time
startTime        = min( startPoint, earliestReaction );     % Start at the baseline lower limit or at the earliest possible reaction time
startTime        = startTime - edge;
startSeconds     = startTime / 1000;

% End time
maxTrialTime     = maxResponseTime + responseLimits(2);
endTime          = maxTrialTime + edge;
endSeconds       = ( endTime + 1 ) / 1000;                  % Correction for pop_epoch

% No blending for median filtering
if  contains( blending, 'm', 'IgnoreCase', true ) && ...
   ~contains( blending, 'sig', 'IgnoreCase', true )
    blendingDuration = Inf;
end

% Blend adjustment for trials with no response (in ms)
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
conditionsNames  = cell( 1, nConditions );
for cn = 1:nConditions
    conditionsNames{cn} = char( join( conditions{cn} ) );
end


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
if contains( blending, 'sig', 'IgnoreCase', true )
    disp( [ num2str( blendingDuration ) ' ms sigmoid blend time for variable trial durations'] )
end
if contains( blending, 'smooth', 'IgnoreCase', true )
    disp( [ iptnum2ordinal( order ) ' order ' num2str( frames * 1000 / samplingRate ) ' ms Savitsky-Golay smoothing for variable trial durations' ] )
end
if contains( blending, 'median', 'IgnoreCase', true )
    disp( [ num2str( neighbourhood(1) ) ' Hz ' num2str( neighbourhood(2) ) ' ms median blend neighbourhood for variable trial durations' ] )
end
if contains( blending, { 'smooth' 'median' }, 'IgnoreCase', true )
    disp( [ num2str( sigmaGauss ) ' standard deviation Gaussian blurring for variable trial durations' ] )
end
disp( 'Trial:')
disp( [ '0 ms to ' num2str( maxResponseTime ) ' ms trial window' ] )
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
    channelsSet     = channelSet;
    if iscell( channelSet )
        % Compare the channel names to the channel co-ordinates labels in
        % the chanlocs.labels field
        iChannelSet = eeg_chaninds( CurrentEEG, channelsSet );
    elseif isnumeric( channelsSet )
        iChannelSet = channelsSet;
        channelsSet = { CurrentEEG.chanlocs(iChannelSet).labels };
    end

    % Loop through: Conditions
    parfor c = 1:nConditions
%     for c = 1:nConditions

        % Current condition event label
        condition     = conditions{c};
        if ischar( condition )
            condition = { condition };
        end

        % Parallel broadcast variables
        EEG           = CurrentEEG;             % Struct copied to each parallel worker as pop_loadeset is slower
        channels      = iChannelSet;            % Channel indices
        baselineLimit = baselineLimits;         % In ms
        stimulusLimit = stimulusLimits;         % In ms
        responseLimit = responseLimits;         % In ms
        blendDuration = blendingDuration;       % In ms

        % Pre-allocate
        Decomposition = [];
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
        eegEvents            = { EEG.event(:).type };
        iStimulusEvent       = matches( eegEvents, condition );
        if sum( iStimulusEvent ) == 1
            singleTrial      = true;
            disp( [ '[' 8 'Warning: Single trial for ' conditionsNames{c} ']' 8 ] )
            % Duplicate single trial stimulus event so EEGLAB produces epochs
            iStimulus        = find( iStimulusEvent );
            EEG.event(end+1) = EEG.event(iStimulus);
            EEG              = pop_editeventvals( EEG, 'sort', { 'latency', 0 } );
            EEG              = eeg_checkset( EEG );
        elseif ~any( iStimulusEvent )
            disp( [ '[' 8 'Warning: No trials for ' conditionsNames{c} ']' 8 ] )
            continue
        else
            singleTrial      = false;
        end
        
        % Epoch trials
        EEG     = pop_epoch( EEG, condition, [startSeconds endSeconds], 'epochinfo', 'yes' );
        EEG     = eeg_checkset( EEG );
        nTrials = length( EEG.epoch );

        % Sanity check: Epoch vs. decomposition window
        if EEG.pnts ~= length( startTime:1000/samplingRate:endTime )
            disp( [ num2str( length( startTime:1000/samplingRate:endTime ) ) ' decomposition window time points' ] )
            disp( [ num2str( EEG.pnts ) ' epoch time points'] )
            disp( [ '[' 8 'Warning: Number of epoch time points does not match ' ...
                          'start time to end time points at the sampling rate'   ']' 8 ] )
        end

        % Deal with conditions in which the first trials are skipped
        skipStart = true; % Set to false when a realistic trial is found
        skipCount = 0;    % Count of unrealistic trials at the start

        % Loop through: Trials
        fprintf( 'Trial' )
        for trial = 1:nTrials
            fprintf( [ ' ' num2str( trial ) ] )

            % Trial events
            currentTrialEvents  = EEG.epoch(trial).eventtype(:);

            % Index of the stimulus event closest to 0 ms in the epoch
            iEpochEventStimulus = find( matches( currentTrialEvents, condition ) );
            if length( iEpochEventStimulus ) > 1
                [ ~, iCurrentStimulus ] = min( abs( [ EEG.epoch(trial).eventlatency{iEpochEventStimulus} ] ) );
                iEpochEventStimulus     = iEpochEventStimulus(iCurrentStimulus);
            end

            % Index of the next response event of any type in the epoch
            iEpochEventResponse = find( matches( currentTrialEvents, responses ) );
            iEpochEventResponse = iEpochEventResponse(iEpochEventResponse > iEpochEventStimulus); % Responses after the stimulus
            if ~isempty( iEpochEventResponse )
                iEpochEventResponse = iEpochEventResponse(1);                                     % First one
            end

            % Trial event-locked centerings
            if ~isempty( iEpochEventResponse )
                nCentres = 2;
            else
                nCentres = 1;
            end

            % Trial event times (in ms) relative to the stimulus (at 0 ms)            
            switch nCentres

                % Stimulus and response
                case 2
                    stimulusTime         = EEG.epoch(trial).eventlatency{iEpochEventStimulus};
                    responseTime         = EEG.epoch(trial).eventlatency{iEpochEventResponse};

                % Stimulus but no response
                case 1
                    stimulusTime         = EEG.epoch(trial).eventlatency{iEpochEventStimulus};
                    responseTime         = EEG.times(end) - adjustment;             % Blend out by the end of the trial
                    [ ~, iResponseTime ] = min( abs( EEG.times - responseTime ) );  % Correct for sampling rate
                    responseTime         = EEG.times(iResponseTime);
            end
            centreTimes = [ stimulusTime responseTime ];

            % Sanity check
            if stimulusTime ~= 0
                disp( [ '[' 8 'Warning: Stimulus at ' num2str( stimulusTime ) ' for ' conditionsNames{c} ' trial ' num2str( trial ) ']' 8 ] )
            end

            % Process trials with realistic reaction time (in ms)
            if responseTime >= minReactionTime

                % A realistic trial has been found
                skipStart = false;

                % Loop through: Channels
%                 fprintf( '\n' )
                for ch = 1:nChannels

                    % Current channel
                    channel     = channels(ch);
%                     channelName = EEG.chanlocs(channel).labels;
%                     fprintf( [ channelName ' '  ] )


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

                    % Times (in ms)
                    %   Wavelet signal indices -> time points in ms
                    timePoints = EEG.times(signalIndices);

                    % Baseline correction
                    [ spectralPower, baseline ] = baselineCorrection( spectralPower, ...
                                                                      timePoints,    ...
                                                                      baselineLimit  );

                    % Power units
                    if ~isempty( baselineLimit ) && any( baselineLimit ~=0 )
                        powerUnits = 'Decibels relative to baseline (mean spectrum)';
                    else
                        powerUnits = 'Decibel volts^2';
                    end


                    %% Extract event-locked data windows
                    % -----------------------------------------------------

                    % Time windows (in ms)
                    %   Relative time limits adjusted by centre event times
                    %   Limit relative to centre -> limit relative to stimulus
                    centreLimits  = { ( stimulusLimit + centreTimes(1) ) ... % + 0
                                      ( responseLimit + centreTimes(2) ) };  % + responseTime

                    % Time domains relative to centre (in ms)
                    %   Centre window time points relative to the centre
                    %   event time
                    iCentreLimits = { and( timePoints >= centreLimits{1}(1), timePoints <= centreLimits{1}(2) ) ...
                                      and( timePoints >= centreLimits{2}(1), timePoints <= centreLimits{2}(2) ) };
                    centreDomains = { ( timePoints(iCentreLimits{1}) - centreTimes(1) ) ...
                                      ( timePoints(iCentreLimits{2}) - centreTimes(2) ) };

                    % Blender times
                    blenderTimes  = { responseTime ...
                                      stimulusTime };

                    % Loop through: Event-locked centres
                    for centre = 1:nTrialCentres

                        % Current centre
                        trialCentre  = trialCentres{centre};
                        centreLimit  = centreLimits{centre};  % Times in ms relative to the stimulus
                        centreDomain = centreDomains{centre}; % Times in ms relative to the centre event
                        blenderTime  = blenderTimes{centre};
                        if centre == 1
                            blenderDuration = -blendDuration; % Stimulus-locked: Blend out at the response
                        else
                            blenderDuration = blendDuration;  % Response-locked: Blend in at the stimulus
                        end

                        % Pre-allocate open windows
                        nWindowTimes = length( centreDomain );
                        powerWindow  = NaN( nFrequencies, nWindowTimes );
                        phaseWindow  = NaN( nFrequencies, nWindowTimes )*1i;
                        wavesWindow  = NaN( nFrequencies, nWindowTimes )*1i;

                        % Skip too-slow reactions found within the epoch
                        % leaving them as NaN
                        if centre == 2 && responseTime > maxResponseTime
                            continue
                        end

                        % Process event-centerings found within the epoch
                        % and skip too-slow reactions not within the epoch
                        % leaving them as NaN
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

                            % Centre window relative to the stimulus
                            %   Bounded by the centre limits relative to
                            %   the stimulus at 0 ms
                            iCentreWindow = and( timePoints >= centreLimit(1), timePoints <= centreLimit(2) );
                            nCentreWindow = sum( iCentreWindow );

                            % Centre the decomposition
                            if centre == 1                                                   % Stimulus
                                powerWindow(:,1:nCentreWindow)         = blendedPower(:,iCentreWindow);
                                phaseWindow(:,1:nCentreWindow)         = blendedPhase(:,iCentreWindow);
                                wavesWindow(:,1:nCentreWindow)         = blendedWaves(:,iCentreWindow);
                            elseif centre > 1 && responseTime <= maxResponseTime             % Response
                                powerWindow(:,end-nCentreWindow+1:end) = blendedPower(:,iCentreWindow);
                                phaseWindow(:,end-nCentreWindow+1:end) = blendedPhase(:,iCentreWindow);
                                wavesWindow(:,end-nCentreWindow+1:end) = blendedWaves(:,iCentreWindow);
                            end


                        end

                        % Sanity checks
                        if isscalar( powerWindow ) || isempty( powerWindow ) || all( isnan( powerWindow ), 'all' )
                            disp( [ '[' 8 'Warning: No decomposition for trial ' num2str( trial ) ' centred on the ' trialCentre ' in ' currentFile ']' 8 ] )
                        elseif all( ~powerWindow(~isnan( powerWindow )), 'all' )
                            disp( [ '[' 8 'Warning: Decomposition of zeroes for trial ' num2str( trial ) ' centred on the ' trialCentre ' in ' currentFile ']' 8 ] )
                        end

                        % Store in struct
                        Decomposition.(trialCentre).SpectralPower(trial,ch,:,:)  = powerWindow;
                        Decomposition.(trialCentre).PhaseDirection(trial,ch,:,:) = phaseWindow;
                        Decomposition.(trialCentre).Coefficients(trial,ch,:,:)   = wavesWindow;
                        Decomposition.(trialCentre).SpectralPowerUnits           = powerUnits;
                        Decomposition.(trialCentre).SpectralPowerDimensions      = 'Channels x Frequencies x Times'; % As it will be after averaging
                        Decomposition.(trialCentre).PhaseCoherenceUnits          = 'Phase alignment proportion';
                        Decomposition.(trialCentre).PhaseCoherenceDimensions     = 'Channels x Frequencies x Times';
                        Decomposition.(trialCentre).CoefficientsUnits            = 'Complex magnitude';
                        Decomposition.(trialCentre).CoefficientsDimensions       = 'Trials x Channels x Frequencies x Times';
                        Decomposition.(trialCentre).Frequencies                  = frequencies;
                        Decomposition.(trialCentre).FrequenciesUnits             = 'Hz';
                        Decomposition.(trialCentre).Times                        = centreDomain;
                        Decomposition.(trialCentre).TimesUnits                   = 'milliseconds';
                        Decomposition.(trialCentre).SamplingRate                 = samplingRate;
                        Decomposition.(trialCentre).SamplingRateUnits            = 'Hz';
                        Decomposition.(trialCentre).BaselinePower(trial,ch,:)    = baseline;
                        Decomposition.(trialCentre).BaselinePowerUnits           = 'Decibel volts^2';
                        Decomposition.(trialCentre).BaselinePowerDimensions      = 'Trials x Channels x Frequencies';
                        Decomposition.(trialCentre).ChannelNames                 = channelsSet;
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
        fprintf( '\n' )

        % Set first trials to NaN if they were skipped
        if skipStart
            for centre = 1:nTrialCentres
            
               trialCentre = trialCentres{centre};
        
               Decomposition.(trialCentre).SpectralPower(1:skipCount,:,:,:)  = NaN;
               Decomposition.(trialCentre).PhaseDirection(1:skipCount,:,:,:) = NaN*1i;
               Decomposition.(trialCentre).Coefficients(1:skipCount,:,:,:)   = NaN*1i;
    
            end
        end


        %% Power, phase coherence, and median blending
        % -----------------------------------------------------------------

        % Loop through: Event-locked centres
        for centre = 1:nTrialCentres

            trialCentre = trialCentres{centre};

            % Time-points with too-few trials for filter blending
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
            if singleTrial && size( Decomposition.(trialCentre).Coefficients, 1 ) == 2
            Decomposition.(trialCentre).Coefficients ...
                = Decomposition.(trialCentre).Coefficients(1,:,:,:);        % Remove duplicate trial from coefficients for single trial cases
            end

            % Sanity checks
            if isscalar( Decomposition.(trialCentre).SpectralPower ) || isempty( Decomposition.(trialCentre).SpectralPower )
                disp( [ '[' 8 'Warning: No decomposition centred on the ' trialCentre ' for ' currentFile ']' 8 ] )
            elseif all( not( Decomposition.(trialCentre).SpectralPower(~isnan( Decomposition.(trialCentre).SpectralPower )) ), 'all' )
                disp( [ '[' 8 'Warning: Decomposition of zeroes centred on the ' trialCentre ' for ' currentFile ']' 8 ] )
            end


            % Median filter (2-D) and/or Savitsky-Golay filter then
            % Gaussian blur trial edges to smooth out discontinuities
            % -------------------------------------------------------------

            if isinf( blendDuration )

                % Median, Savitsky-Golay, and Gaussian filter blending
                Decomposition = filterBlend( Decomposition, ...
                                             centre,        ...
                                             trialCentre,   ...
                                             blending,      ...
                                             neighbourhood, ...
                                             order,         ...
                                             frames,        ...
                                             sigmaGauss     );               

                % Erase time-points with too-few trials
                Decomposition.(trialCentre).SpectralPower(:,:,iTooFew)  = 0;
                Decomposition.(trialCentre).PhaseCoherence(:,:,iTooFew) = 0;
                Decomposition.(trialCentre).Coefficients(:,:,iTooFew)   = NaN*1i;

            end

        end % for Centres


        %% Save single participant x condition
        % -----------------------------------------------------------------

        % File name
        [ iP, iXn ]  = regexp( currentFile, [ participantCode '[0-9]+' ] );
        participant  = currentFile(iP:iXn);
        fileName     = [ 'TimeFrequencyData' participant conditionsNames{c} ];
        fileFullPath = fullfile( currentFolder, fileName );

        % Save .mat file
        timeFrequencyParallelSave( fileFullPath, Decomposition )

        % Single completion time
        timeText    = [ currentFile ' condition ' conditionsNames{c} ' decomposition completed at' ];
        timeFrequencyRunTime( timeText )

        
    end % parfor Conditions

end % for Files


% Finish time
timeFrequencyRunTime( 'Time-frequency data generation completed at' )


% _________________________________________________________________________
end



%%
% •.° Baseline Correction °.•
% _________________________________________________________________________
function ...
[ spectralPower, baseline ] = baselineCorrection( spectralPower, timePoints, baselineLimit )

% Baseline mean spectrum
if ~isempty( baselineLimit ) && any( baselineLimit ~=0 )
    iBaselineTimes = and( timePoints >= baselineLimit(1), timePoints <= baselineLimit(2) );
    baseline       = spectralPower(:,iBaselineTimes);
    baseline       = mean( baseline, 2, 'omitnan' );

% No baseline
elseif isempty( baselineLimit ) || all( baselineLimit == 0 )
    nFrequencies   = size( spectralPower, 1 );
    baseline       = ones( nFrequencies, 1 );

end

% Convert baseline power to decibel volts^2
baseline      = 10*log10( baseline );

% Baseline spectrum copied across the time window
baselinePower = repmat( baseline, 1, size( spectralPower, 2 ) );

% Sanity check: Matching size
if size( baselinePower, 1 ) ~= size( spectralPower, 1 ) || ...
   size( baselinePower, 2 ) ~= size( spectralPower, 2 )
    error( 'Baseline size mismatch' )
end

% Convert spectral power to decibel volts^2
spectralPower = 10*log10( spectralPower );

% Subtract the baseline mean spectrum in logarithmic units
%   This gives power in decibels relative to baseline or retains power in 
%   decibel volts^2 if no baseline
spectralPower = spectralPower - baselinePower;


% _________________________________________________________________________
end



%%
% •.° Median, Savitsky-Golay, and Gaussian filter blending °.•
% _________________________________________________________________________
function ...
Decomposition = filterBlend( Decomposition, centre, trialCentre, blending, ...
                             neighbourhood, order, frames, sigmaGauss )

% Decomposition time data
times        = Decomposition.(trialCentre).Times;
nTimes       = length( times );
samplingRate = Decomposition.(trialCentre).SamplingRate;


% Variable trial length time points to filter
% -------------------------------------------------------------------------
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
padding     = neighbourhood(2) * 1000 / samplingRate;
iEdgeAA     = (iEdgeA1 - padding):iEdgeA2;
iEdgeBB     = iEdgeB1:(iEdgeB2 + padding);


% Filtering
% -------------------------------------------------------------------------

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
for ch = 1:length( Decomposition.(trialCentre).ChannelNames )

    % Current channel spectra
    channelPower = Decomposition.(trialCentre).SpectralPower(ch,:,:);
    channelPhase = Decomposition.(trialCentre).PhaseCoherence(ch,:,:);
    channelPower = squeeze( channelPower );
    channelPhase = squeeze( channelPhase );

    % Edges after the centre event to the end
    if after

        % Spectral power
        edgeSpectra     = channelPower(:,iEdgeAA);
        if contains( blending, 'smooth', 'IgnoreCase', true )
            edgeSpectra = sgolayfilt( edgeSpectra, order, frames, [], 2 );
        end
        if contains( blending, 'median', 'IgnoreCase', true )
            edgeSpectra = medfilt2( edgeSpectra, neighbourhood );
        end
        edgeSpectra     = imgaussfilt( edgeSpectra, sigmaGauss );
        Decomposition.(trialCentre).SpectralPower(ch,:,iEdgeA) ...
                        = edgeSpectra(:,1+padding:end);

        % Phase coherence
        edgeSpectra     = channelPhase(:,iEdgeAA);
        if contains( blending, 'smooth', 'IgnoreCase', true )
            edgeSpectra = sgolayfilt( edgeSpectra, order, frames, [], 2 );
        end
        if contains( blending, 'median', 'IgnoreCase', true )
            edgeSpectra = medfilt2( edgeSpectra, neighbourhood );
        end
        edgeSpectra     = imgaussfilt( edgeSpectra, sigmaGauss );
        Decomposition.(trialCentre).PhaseCoherence(ch,:,iEdgeA) ...
                        = edgeSpectra(:,1+padding:end);

    end

    % Edges before the centre event to the end
    if before

        % Spectral power
        edgeSpectra     = channelPower(:,iEdgeBB);
        if contains( blending, 'smooth', 'IgnoreCase', true )
            edgeSpectra = sgolayfilt( edgeSpectra, order, frames, [], 2 );
        else
            edgeSpectra = medfilt2( edgeSpectra, neighbourhood );
        end
        edgeSpectra     = imgaussfilt( edgeSpectra, sigmaGauss );
        Decomposition.(trialCentre).SpectralPower(ch,:,iEdgeB) ...
                        = edgeSpectra(:,1:end-padding);

        % Phase coherence
        edgeSpectra     = channelPhase(:,iEdgeBB);
        if contains( blending, 'smooth', 'IgnoreCase', true )
            edgeSpectra = sgolayfilt( edgeSpectra, order, frames, [], 2 );
        else
            edgeSpectra = medfilt2( edgeSpectra, neighbourhood );
        end
        edgeSpectra     = imgaussfilt( edgeSpectra, sigmaGauss );
        Decomposition.(trialCentre).PhaseCoherence(ch,:,iEdgeB) ...
                        = edgeSpectra(:,1:end-padding);

    end

end % for Channels


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


