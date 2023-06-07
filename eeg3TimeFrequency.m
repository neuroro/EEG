function ...
eeg3TimeFrequency( setName, frequencyLimits, frequencyResolution, blending )
%
% •.° Time-Frequency Decompositions of Event-Related EEG °.•
% _________________________________________________________________________
%
% Calculate event-related spectral power and inter-trial phase coherence of
% EEG data for each trial, centred on the stimulus presentation, initiation
% of the response, and the response decision, for all EEGLAB datasets named
% in common that are located in the Current Folder and sub-folders of the
% Current Folder, then save a .mat file for each participant x condition to
% the Current Folder
%
% Trials are assumed to contain a stimulus presentation event labelled per
% condition, a response-initiation event, and a response-decision event
%
% EEG data may be continuous or epoched but epochs must be wider than the
% trial by half a wavelength of the lowest frequency to be decomposed
%
% • Usage •
% -------------------------------------------------------------------------
% Edit the configuration of the EEG data in the code
% Set the Current Folder to the location of the datasets or the base folder
%
% >> eeg3TimeFrequency( setName, frequencyLimits, frequencyResolution, blending )
%
% Inputs:
%   setName:             Part of the file name that is common to all the 
%                        EEGLAB datasets to be decomposed that are located 
%                        in the Current Folder and sub-folders of the
%                        Current Folder, for example 'Memory' or ''
%                          (optional input, default all datasets)
%   frequencyLimits:     Frequencies in Hz to extract as [minimum maximum]
%                          (optional input, default [2 30])
%   frequencyResolution: Frequencies per octave
%                          (optional input, default 30)
%   blending:            Blend the edges of variable-length trials using a 
%                        'sigmoid' taper to 0 dB centred around adjacent 
%                        events or a 'median' filter then Gaussian blur of
%                        the time points with trial dropout
%                          (optional input, default sigmoid blending)
%
% Output:
%   TimeFrequencyDataP<number><condition>.mat files containing data centred
%   on the stimulus, initiation, and response for each participant and each
%   condition, from which a localised decomposition can then be extracted
%   using eegLocalSpectra.m
%
% Examples:
% >> eeg3TimeFrequency( 'Memory' )
% >> eeg3TimeFrequency( 'Memory', [], [], 'median' )
% >> eeg3TimeFrequency( 'Memory', [ 2 30 ], 30, 'sigmoid' )
%
% • Requires •
% -------------------------------------------------------------------------
% waveletTransform.m and spectralBlender.m by Rohan King
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
disp( 'EEG data for each trial, centred on the stimulus presentation, initiation' )
disp( 'of the response, and the response decision, for all EEGLAB datasets named' )
disp( 'in common that are located in the Current Folder and sub-folders of the'   )
disp( 'Current Folder, then save a .mat file for each participant x condition to' )
disp( 'the Current Folder' )
disp( ' ' )


%% Configuration
% -------------------------------------------------------------------------

% Participant code in the file name that precedes the participant number
participantCode  = 'P';

% Stimulus event labels for each condition
conditions       = { 'cC-R' 'cC-S' ...
                     'iC-R' 'iC-S' ...
                     'cI-R' 'cI-S' ...
                     'iI-R' 'iI-S' };

% Order in the trial of the initiation and response events relative to the
% stimulus event, which is assumed to be 1st
initiationEvent  = 2; % 2nd after stimulus
responseEvent    = 3; % 3rd after stimulus

% Non-events, such as fixation, which should not be windowed around
nonEvents        = { 'Fixation' 'fixation' 'FXTN' 'Empty' 'empty' };

% Time parameters
samplingRate     = 1000;
baselineLimits   = [ -500 0    ]; % Relative to stimulus time
stimulusLimits   = [ -200 1000 ]; % Relative to stimulus time
initiationLimits = [ -500 500  ]; % Relative to initiation time
responseLimits   = [ -500 500  ]; % Relative to response time

% Frequency parameters
if nargin < 2 || isempty( frequencyLimits )
    frequencyLimits = [ 2 60 ]; % Minium and maximum frequency
end
if nargin < 3 || isempty( frequencyResolution ) || ~frequencyResolution
    frequencyResolution = 30;   % Frequencies per octave
end
maxTrialWindow   = 2000;

% Trial minima
minReactionTime  = 100; % Minimum initiation reaction time
nTrialsMinimum   = 10;  % Time points with fewer than this number are zeroed

% Blending
blendingDuration = 80;     % Sigmoid blend time points around adjacent event
neighbourhood    = [3 30]; % Median window frequencies x time points
sigmaGauss       = 2;      % Gaussian blur standard deviations
if nargin < 4
    blending     = 'sigmoid';
end

% Channels of interest (indices in the EEG data)
frontalChannels  = [ 5 6 11 12 ];  % [ 5  6  7  106 11 12 118 20 19 4  13 112 ];
parietalChannels = [ ];            % [ 31 80 55 54  79 62 37  87 86 53 60 85  ];


%% Derived parameters
% -------------------------------------------------------------------------

% Events
initiationEvent  = initiationEvent - 1;
responseEvent    = responseEvent   - 1;

% Times
longestCycle     = samplingRate/frequencyLimits(1);
edge             = ceil( pi/2 * longestCycle );     % Excised after decomposition to remove edge effects
if isempty( baselineLimits ) || baselineLimits(1) > -200
    startTime    = -200;                            % Start no later than -200 ms so that pre-stimulus data exists
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
channelSet       = [ frontalChannels parietalChannels ];
channelSet       = sort( channelSet );
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
disp( [ num2str( minReactionTime ) ' ms minimum initiation reaction time' ] )
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
disp( [ num2str( initiationLimits(1) ) ' ms to ' num2str( initiationLimits(2) ) ...
        ' ms relative to initiation of the response' ] )
disp( [ num2str( responseLimits(1)   ) ' ms to ' num2str( responseLimits(2) ) ...
        ' ms relative to the response decision' ] )
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
    currentFile   = fileList{n};
    currentFolder = folderList{n};
    disp( [ 'Decomposing ' currentFile ] )
    disp( ' ' )

    % Load
    loaderargin   = { 'filename', currentFile,   ...
                      'filepath', currentFolder, ...
                      'verbose',  'off'          };
    CurrentEEG    = pop_loadset( loaderargin{:} );
    CurrentEEG    = eeg_checkset( CurrentEEG );

    % Sanity check: Sampling rate
    if CurrentEEG.srate ~= samplingRate
        error( [ 'Sampling rate in file not equal to ' num2str( samplingRate ) ] )
    end

%     % Parallel loop through: Conditions
%     parpool;
%     parfor c = 1:nConditions

    % Loop through: Conditions
    for c = 1:nConditions

        % Current condition event label
        condition       = conditions{c};

        % Parallel broadcast variables
        EEG             = CurrentEEG;
        channels        = channelSet;
        baselineLimit   = baselineLimits;
        stimulusLimit   = stimulusLimits;
        initiationLimit = initiationLimits;
        responseLimit   = responseLimits;
        blendDuration   = blendingDuration;

        % Pre-allocate
        Decomposition   = [];
        trialCentres    = { 'Stimulus' 'Initiation' 'Response' };
        nTrialCentres   = length( trialCentres );
        for tc = 1:nTrialCentres
            Decomposition.(trialCentres{tc}).SpectralPower  = [];
            Decomposition.(trialCentres{tc}).PhaseDirection = [];
            Decomposition.(trialCentres{tc}).PhaseCoherence = [];
            Decomposition.(trialCentres{tc}).Coefficients   = [];
        end


        %% Trial data
        % -----------------------------------------------------------------

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
            currentTrialEvents    = EEG.epoch(trial).eventtype(:);
            nCurrentTrialEvents   = length( currentTrialEvents );

            % Stimulus event index
            %   Find the index of the event at 0 ms in case a preceding
            %   event or two made it into the epoch
            iEpochEventStimulus   = find( ~[ EEG.epoch(trial).eventlatency{:} ] );

            % Initiation and response event indices
            %   Using the stimulus event index and the trial
            % event order and the 
            iEpochEventInitiation = iEpochEventStimulus + initiationEvent;
            iEpochEventResponse   = iEpochEventStimulus + responseEvent;

            % Non-events to exclude
            iPreEvents            = 1:nCurrentTrialEvents < iEpochEventStimulus; % Indices of events before the stimulus
            iNonEvents            = contains( currentTrialEvents, nonEvents );   % Indices of non-events such as fixation
            iExclude              = or( iPreEvents', iNonEvents );

            % Trial event-locked centerings
            currentTrialCentres   = currentTrialEvents(~iExclude);
            nCentres              = length( currentTrialCentres );

            % Trial event times (in ms) relative to 0 ms
            switch nCentres

                % Stimulus only, no initiation or response
                case 1
                    stimulusTime   = 0;
                    initiationTime = EEG.times(end) - adjustment - 1; % Blend out by the end of the trial
                    responseTime   = NaN;

                % No response
                case 2
                    stimulusTime   = 0;
                    initiationTime = EEG.epoch(trial).eventlatency{iEpochEventInitiation};
                    responseTime   = EEG.times(end) - adjustment - 1; % Blend out by the end of the trial

                % All three
                case 3
                    stimulusTime   = 0;
                    initiationTime = EEG.epoch(trial).eventlatency{iEpochEventInitiation};
                    responseTime   = EEG.epoch(trial).eventlatency{iEpochEventResponse};

            end
            centreTimes = [ stimulusTime initiationTime responseTime ];

            % Process trials with realistic initiation time (in ms)
            if initiationTime >= minReactionTime

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
                        baseline   = 1;
                    end

                    % Convert power to dB relative to baseline
                    spectralPower  = 10*log10( spectralPower ./ baseline );

                    % Times (in ms)
                    %   Input signal indices -> output time points in ms
                    timePointsTest = startTime + signalIndices - 1;
                    timePoints     = EEG.times(signalIndices);
                    if samplingRate == 1000 && any( timePoints ~= timePointsTest ) % Sanity check
                        error( 'Time points calculation is inconsistent' )
                    end


                    %% Extract event-locked data windows
                    % -----------------------------------------------------

                    % Time windows (in ms)
                    %   Relative time limits adjusted by centre event times
                    %   Limit relative to centre -> limit relative to 0
                    centreLimits  = { (stimulusTime   + stimulusLimit  ) ...
                                      (initiationTime + initiationLimit) ...
                                      (responseTime   + responseLimit  ) };

                    % Time domains relative to centre
                    centreDomains = { stimulusLimit(1):stimulusLimit(2)     ...
                                      initiationLimit(1):initiationLimit(2) ...
                                      responseLimit(1):responseLimit(2) };

                    % Blender times
                    blenderTimes  = { initiationTime                ...
                                      [ stimulusTime responseTime ] ...
                                      initiationTime };

                    % Loop through: Event-locked centres
                    for centre = 1:nTrialCentres

                        % Current centre
                        trialCentre  = trialCentres{centre};
                        centreLimit  = centreLimits{centre};
                        centreDomain = centreDomains{centre};
                        blenderTime  = blenderTimes{centre};
                        if centre == 1
                            blenderDuration = -blendDuration; % Blend out at initiation
                        else
                            blenderDuration = blendDuration;
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

                        % Store in struct
                        Decomposition.(trialCentre).SpectralPower(trial,ch,:,:)  = powerWindow;
                        Decomposition.(trialCentre).PhaseDirection(trial,ch,:,:) = phaseWindow;
                        Decomposition.(trialCentre).Coefficients(trial,ch,:,:)   = wavesWindow;
                        Decomposition.(trialCentre).Frequencies         = frequencies;
                        Decomposition.(trialCentre).Times               = centreDomain;
                        Decomposition.(trialCentre).Channels            = channels;
                        Decomposition.(trialCentre).BlendTimes(trial,:) = blenderTime - centreTimes(centre);


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

            % Time-points with too-few trials
            %   Blending with a median filter won't prevent single extreme
            %   trials from skewing the average
            if isinf( blendDuration )
                missingPoints  = isnan( Decomposition.(trialCentre).SpectralPower );
                missingPoints  = squeeze( missingPoints(:,1,1,:) ); % All channels and frequencies
                trialsPerPoint = sum( ~missingPoints, 1 );          %    have an equal trial count
                iTooFew        = find( trialsPerPoint < nTrialsMinimum );
            end


            % Induced and evoked oscillations
            % -------------------------------------------------------------

            % Mean power across trials
            Decomposition.(trialCentre).SpectralPower  ...
                = mean( Decomposition.(trialCentre).SpectralPower, 1, 'omitnan' );

            % Inter-trial phase coherence
            Decomposition.(trialCentre).PhaseCoherence ...
                = abs( mean( Decomposition.(trialCentre).PhaseDirection, 1, 'omitnan' ) );

            % Clean up
            Decomposition.(trialCentre).SpectralPower  ...
                = squeeze( Decomposition.(trialCentre).SpectralPower  );
            Decomposition.(trialCentre).PhaseCoherence ...
                = squeeze( Decomposition.(trialCentre).PhaseCoherence );
            Decomposition.(trialCentre) ...
                = rmfield( Decomposition.(trialCentre), 'PhaseDirection' );


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
                padding     = neighbourhood(2) * 1000 / samplingRate;
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

        % Save mat file
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



%%
% •.° Parallel Save °.•
% _________________________________________________________________________
function timeFrequencyParallelSave( fileName, structOfVariablesToSave )

save( fileName, '-struct', 'structOfVariablesToSave', '-v7.3' )


% _________________________________________________________________________
end






%% Light median filtering on times with high trial count
% Currently bugged for trials with minimal trial count dropoff as it makes
% the time point limits outside the range, sometimes due to the padding but
% sometimes regardless of the padding = needs reworking to stay in range
% 
% % Median filter trial edges to smooth out discontinuities
%             % -------------------------------------------------------------
% 
%             if isinf( blendDuration )
%                 times  = Decomposition.(trialCentre).Times;
%                 nTimes = length( times );
% 
% 
%                 % Variable trial length time points to filter
%                 % ---------------------------------------------------------
%                 % From the trial cut point closest to the centre time
%                 % All time points with fewer trials are filtered
% 
%                 % Time points with >75% trials are lightly filtered
%                 lightPoints     = trialsPerPoint > round( nTrials*0.75 );
%                 neighboursLight = round( [ neighbourhood(1)*1/3 neighbourhood(2)*0.75 ] );
% 
%                 % Edges after the centre event to the end
%                 %   Index range of variable initiation times relative to stimulus
%                 %   Index range of variable response times relative to initiation
%                 blendTimes   = Decomposition.(trialCentre).BlendTimes;
%                 if size( blendTimes, 2 ) == 2
%                     blendTimes   = blendTimes(:,2);
%                 end
%                 earliestCut  = min( abs( blendTimes ), [], 'omitnan' );
%                 iEdgeALight1 = find( times == earliestCut );
%                 iEdgeALight2 = find( lightPoints, 1, 'last' );
%                 iEdgeALight  = iEdgeALight1:iEdgeALight2;
%                 iEdgeA1      = iEdgeALight2 + 1;
%                 iEdgeA2      = nTimes;
%                 iEdgeA       = iEdgeA1:iEdgeA2;
% 
%                 % Edges before the centre event from the start
%                 %   Index range of variable response times relative to initiation
%                 %   Index range of variable initiation times relative to response
%                 blendTimes   = Decomposition.(trialCentre).BlendTimes;
%                 if size( blendTimes, 2 ) == 2
%                     blendTimes   = blendTimes(:,1);
%                 end
%                 iEdgeBLight1 = 1;
%                 iEdgeBLight2 = find( lightPoints, 1 );
%                 iEdgeBLight  = iEdgeBLight1:iEdgeBLight2;
%                 iEdgeB1      = iEdgeBLight2 + 1;
%                 latestCut    = min( abs( blendTimes ), [], 'omitnan' );
%                 iEdgeB2      = find( times == -latestCut ); % Negative relative to the centre event
%                 iEdgeB       = iEdgeB1:iEdgeB2;
% 
%                 % Pad the edges by the neighbourhood
%                 padding      = neighbourhood(2) * 1000 / samplingRate;
%                 iEdgeAALight = (iEdgeALight1 - padding):(iEdgeALight2 + padding);
%                 iEdgeAA      = (iEdgeA1 - padding):iEdgeA2;
%                 iEdgeBBLight = iEdgeBLight1:(iEdgeBLight2 + padding);
%                 iEdgeBB      = (iEdgeB1 - padding):(iEdgeB2 + padding);
% 
% 
%                 % Median filter
%                 % ---------------------------------------------------------
%                 
%                 % Edges to filter
%                 switch centre
%                     case 1 % Stimulus-locked
%                         before = false;
%                         after  = true;
%                     case 2 % Initiation-locked
%                         before = true;
%                         after  = true;
%                     case 3 % Response-locked
%                         before = true;
%                         after  = false;
%                 end
% 
%                 % Loop through: Channels
%                 for ch = 1:nChannels
% 
%                     % Current channel spectra
%                     channelPower = Decomposition.(trialCentre).SpectralPower(ch,:,:);
%                     channelPhase = Decomposition.(trialCentre).PhaseCoherence(ch,:,:);
%                     channelPower = squeeze( channelPower );
%                     channelPhase = squeeze( channelPhase );
% 
%                     % Edges after the centre event to the end
%                     if after
% 
%                         % Spectral power
%                         edgeLightSpectra = channelPower(:,iEdgeAALight);
%                         edgeSpectra      = channelPower(:,iEdgeAA);
%                         edgeLightSpectra = medfilt2( edgeLightSpectra, neighboursLight );
%                         edgeSpectra      = medfilt2( edgeSpectra, neighbourhood );
%                         edgeLightSpectra = imgaussfilt( edgeLightSpectra, sigmaGauss );
%                         edgeSpectra      = imgaussfilt( edgeSpectra,      sigmaGauss );
%                         Decomposition.(trialCentre).SpectralPower(ch,:,iEdgeALight) ...
%                                          = edgeLightSpectra(:,1+padding:end-padding);
%                         Decomposition.(trialCentre).SpectralPower(ch,:,iEdgeA) ...
%                                          = edgeSpectra(:,1+padding:end);
% 
%                         % Phase coherence
%                         edgeLightSpectra = channelPhase(:,iEdgeAALight);
%                         edgeSpectra      = channelPhase(:,iEdgeAA);
%                         edgeLightSpectra = medfilt2( edgeLightSpectra, neighboursLight );
%                         edgeSpectra      = medfilt2( edgeSpectra, neighbourhood );
%                         edgeLightSpectra = imgaussfilt( edgeLightSpectra, sigmaGauss );
%                         edgeSpectra      = imgaussfilt( edgeSpectra,      sigmaGauss );
%                         Decomposition.(trialCentre).PhaseCoherence(ch,:,iEdgeALight) ...
%                                          = edgeLightSpectra(:,1+padding:end-padding);
%                         Decomposition.(trialCentre).PhaseCoherence(ch,:,iEdgeA) ...
%                                          = edgeSpectra(:,1+padding:end);
% 
%                     end
% 
%                     % Edges before the centre event to the end
%                     if before
% 
%                         % Spectral power
%                         edgeLightSpectra = channelPower(:,iEdgeBBLight);
%                         edgeSpectra      = channelPower(:,iEdgeBB);
%                         edgeLightSpectra = medfilt2( edgeLightSpectra, neighboursLight );
%                         edgeSpectra      = medfilt2( edgeSpectra, neighbourhood );
%                         edgeLightSpectra = imgaussfilt( edgeLightSpectra, sigmaGauss );
%                         edgeSpectra      = imgaussfilt( edgeSpectra,      sigmaGauss );
%                         Decomposition.(trialCentre).SpectralPower(ch,:,iEdgeBLight) ...
%                                          = edgeLightSpectra(:,1:end-padding);
%                         Decomposition.(trialCentre).SpectralPower(ch,:,iEdgeB) ...
%                                          = edgeSpectra(:,1+padding:end-padding);
% 
%                         % Phase coherence
%                         edgeLightSpectra = channelPhase(:,iEdgeBBLight);
%                         edgeSpectra      = channelPhase(:,iEdgeBB);
%                         edgeLightSpectra = medfilt2( edgeLightSpectra, neighboursLight );
%                         edgeSpectra      = medfilt2( edgeSpectra, neighbourhood );
%                         edgeLightSpectra = imgaussfilt( edgeLightSpectra, sigmaGauss );
%                         edgeSpectra      = imgaussfilt( edgeSpectra,      sigmaGauss );
%                         Decomposition.(trialCentre).PhaseCoherence(ch,:,iEdgeBLight) ...
%                                          = edgeLightSpectra(:,1:end-padding);
%                         Decomposition.(trialCentre).PhaseCoherence(ch,:,iEdgeB) ...
%                                          = edgeSpectra(:,1+padding:end-padding);
% 
%                     end
% 
%                 end % for Channels
%                 
% 
%                 % Erase time-points with too-few trials
%                 Decomposition.(trialCentre).SpectralPower(:,:,iTooFew)  = 0;
%                 Decomposition.(trialCentre).PhaseCoherence(:,:,iTooFew) = 0;
%                 Decomposition.(trialCentre).Coefficients(:,:,iTooFew)   = NaN*1i;
% 
% 
%             end % if Median filter