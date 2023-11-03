function eegLocalSpectra( channels, weighting, imputation )
%
% •.° Localised Event-Related Spectra °.•
% _________________________________________________________________________
%
% Extract event-related spectra at an electrode or at a local cluster of
% electrodes and form local average spectra, for each event-related window,
% from all TimeFrequencyData*.mat files located in the Current Folder and
% sub-folders of the Current Folder, then save all localised data in one
% TimeFrequency<Location>.mat file
%
% • Usage •
% -------------------------------------------------------------------------
% Set the Current Folder to the location of TimeFrequencyData*.mat files
% (generated by eegTimeFrequency.m or eeg3TimeFreqyency.m), which are
% stored in that folder or in sub-folders of that folder, ignoring other
% .mat files
%
% >> eegLocalSpectra( channels, imputation )
%
% For example:
% >> eegLocalSpectra
% >> eegLocalSpectra( 6 )
% >> eegLocalSpectra( 'FCz' )
% >> eegLocalSpectra( [11 6 5 12] )
% >> eegLocalSpectra( { 'Fz' 'FCz' 'FFC1' 'FFC2' }, 0, 10 )E
%
% .. . .  .   .     .        .             .                     .                                  .                                                       .
%
% •••( Function Inputs )
%
%   channels:   Extract event-related spectra from a channel name or index;
%                 from the average of a cluster of channels as a cell array
%                 of channel names or as a vector of channel indices; or
%                 from all channels as 'all', '', [], or 0
%                 (optional input, default all)
%
%   weighting:  Weighted or unweighted average of channels as
%                 'unweighted', false, 0
%                 'weighted',   true,  1
%                 Weighted electrode clusters are averaged with the first
%                 electrode at full-weight and other (adjacent or nearby)
%                 electrodes weighted in proportion to the ratio of the
%                 standard normal Gaussian probability density at 1 + 0.05
%                 per cluster electrode in excess of three to the standard
%                 normal Gaussian probability density at 0.
%                 (optional input, default unweighted)
%
%   imputation: Threshold for imputation as the minimum number of trials
%                 considered usable, below which imputation will be used
%                 (optional input, default 10)
%
% .. . .  .   .     .        .             .                     .                                  .                                                       .
%
% [ Function Output ] =
%
%   TimeFrequency<Location>.mat file containing, for each event-related
%   windowing, pooled decompositions from all participants per condition,
%   ready for peak finding or cluster statistics, and grand averages per
%   condition ready for plotting
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
disp( 'Extract event-related spectra at an electrode or at a local cluster of'    )
disp( 'electrodes and form local average spectra, for each event-related window,' )
disp( 'from all TimeFrequencyData*.mat files located in the Current Folder and'   )
disp( 'sub-folders of the Current Folder, then save all localised data in one'    )
disp( 'TimeFrequency<Location>.mat file'                                          )
disp( ' ' )

codeRunTime( 'Started at' )
disp( ' ' )


%% Default inputs
% -------------------------------------------------------------------------

if ~nargin || ( isscalar( channels ) && ~channels ) || isempty( channels )
    channels = [];
end

if nargin < 2 || ( isscalar( weighting ) && ~weighting ) || isempty( weighting )
    weighting = false;
end

if nargin < 3 || ~isscalar( imputation ) || ~isnumeric( imputation )
    minimumTrialCount = 10;
else
    minimumTrialCount = imputation;
end


%% Time-frequency decomposition .mat files
% -------------------------------------------------------------------------

% Wildcard file name
aWildDataset    = 'TimeFrequencyData*.mat';

% Search for all datasets named in common located in the Current Folder and
% sub-folders of the Current Folder
fileStruct      = dir( [ '**/' aWildDataset ] );
nFiles          = length( fileStruct );
fileList        = { fileStruct(:).name   };
folderList      = { fileStruct(:).folder };


%% Determine participants and conditions represented by the .mat files
% -------------------------------------------------------------------------

% Pre-allocate
participantList = {}; % Participant number files list
participants    = []; % Unique participant numbers
iParticipant    = 0;
conditionList   = {}; % Condition files list
conditions      = {}; % Unique conditions
iCondition      = 0;

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
        iParticipant = iParticipant + 1;

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
        iCondition = iCondition + 1;

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

% Determine parameters from the first file with complete data
determining = true;
d           = 1;
while determining

    % Try files in order
    try

        % Load decomposition
        Decomposition = load( fullfile( folderList{d}, fileList{d} ) );

        % Event-centred windowings
        eventCentres  = fieldnames( Decomposition );
        nCentres      = length( eventCentres );

        % Conditions
        for w = 1:nCentres
            LocalSpectra.(eventCentres{w}).Conditions = conditions;
        end

        % Frequencies
        for w = 1:nCentres
            LocalSpectra.(eventCentres{w}).Frequencies          = Decomposition.(eventCentres{w}).Frequencies;
            try
                LocalSpectra.(eventCentres{w}).FrequenciesUnits = Decomposition.(eventCentres{w}).FrequenciesUnits;
            catch
                LocalSpectra.(eventCentres{w}).FrequenciesUnits = 'Hz';
            end
        end
        nFrequencies = length( LocalSpectra.(eventCentres{1}).Frequencies );

        % Times
        %   Time points vary per window centre
        for w = 1:nCentres
            LocalSpectra.(eventCentres{w}).Times                 = Decomposition.(eventCentres{w}).Times;
            try
                LocalSpectra.(eventCentres{w}).TimesUnits        = Decomposition.(eventCentres{w}).TimesUnits;
                LocalSpectra.(eventCentres{w}).SamplingRate      = Decomposition.(eventCentres{w}).SamplingRate;
                LocalSpectra.(eventCentres{w}).SamplingRateUnits = Decomposition.(eventCentres{w}).SamplingRateUnits;
            catch
                LocalSpectra.(eventCentres{w}).TimesUnits        = 'milliseconds';
                LocalSpectra.(eventCentres{w}).SamplingRate      = [];
                LocalSpectra.(eventCentres{w}).SamplingRateUnits = 'Hz';
            end
            nTimes.(eventCentres{w}) = length( LocalSpectra.(eventCentres{w}).Times );
        end
        
        % Check decomposition fields exist
        for w = 1:nCentres
            fieldCheck = Decomposition.(eventCentres{w}).SpectralPower;                     %#ok
            fieldCheck = Decomposition.(eventCentres{w}).PhaseCoherence;                    %#ok
            fieldCheck = Decomposition.(eventCentres{w}).Coefficients;                      %#ok
        end

        % End determination
        determining = false;

    % Next file if an error occurs
    catch
        
        d = d + 1;

    end

end

% Spectra units and baseline correction
try
    spectralPowerUnits  = Decomposition.(eventCentres{1}).SpectralPowerUnits;
    phaseCoherenceUnits = Decomposition.(eventCentres{1}).PhaseCoherenceUnits;
catch
    try
        oneWave   = squeeze( Decomposition.Stimulus.Coefficients(:,1,1,:) );            % All trials
        onePower  = squeeze( Decomposition.Stimulus.SpectralPower(1,1,:) );
    catch
        oneWave   = squeeze( Decomposition.(eventCentres{1}).Coefficients(:,1,1,:) );   % All trials
        onePower  = squeeze( Decomposition.(eventCentres{1}).SpectralPower(1,1,:) );
    end
    oneWavePower  = squeeze( mean( 10*log10( oneWave .* conj( oneWave ) ), 1, 'omitnan' ) );
    maximumPower  = max( onePower, [], 'all', 'omitnan' );
    meanAbsPower  = mean( abs( onePower ), 'all', 'omitnan' );
    oneDifference = mean( abs( oneWavePower - onePower ), 'all', 'omitnan' );
%     disp( [ num2str( oneDifference ) ' decibel mean absolute difference between wavelet power and spectral power' ] )
    if oneDifference > maximumPower && oneDifference > meanAbsPower
        disp( [ '[' 8 'Baseline correction assumed from mean absolute power difference' ']' 8 ] )
        spectralPowerUnits = 'Decibels relative to baseline (mean spectrum)';
    else
        disp( [ '[' 8 'No baseline correction assumed from mean absolute power similarily' ']' 8 ] )
        spectralPowerUnits = 'Decibel volts^2';
    end
    disp( ' ' )
    phaseCoherenceUnits = 'Phase alignment proportion';
end


% Channels
% -------------------------------------------------------------------------

% Inputs
if matches( 'all', channels, 'IgnoreCase', true )
    channels = [];
end
if ischar( weighting ) || isstring( weighting )
    if matches( 'unweighted', weighting, 'IgnoreCase', true )
        weighting = false;
    elseif matches( 'weighted', weighting, 'IgnoreCase', true )
        weighting = true;
    end
end

% Channels in the decomposition
try
    iDecompositionChannels = Decomposition.(eventCentres{1}).ChannelIndices;    % Channel indices of the EEG data (and of the channel co-ordinates)
catch
    iDecompositionChannels = Decomposition.(eventCentres{1}).Channels;
end
ChannelCoordinates = Decomposition.(eventCentres{1}).ChannelCoordinates;        % Channel co-ordinates

% Selected channels from the decomposition
if ~isempty( channels )

    % Number of selected channels
    nChannels = length( channels );

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
        if isempty( iCurrentChannel )
            error( [ 'Channel ' currentChannel ' is not found in the decomposition' ] )
        end
        iChannels(ch) = find( iDecompositionChannels == iCurrentChannel, 1 );                                  %#ok

    end

% All channels in the decomposition
elseif isempty( channels ) || length( channels ) == length( iDecompositionChannels )

    % Number of channel in the decomposition
    nChannels = length( iDecompositionChannels );

    % Indices of the set of channel indices in the decomposition
    iChannels = 1:nChannels;

% Incorrect input
elseif length( channels ) > length( iDecompositionChannels )
    error( [ 'Error: You selected a greater number of channels than ' ...
             'exist in your time-frequency data files. '              ...
             'Fix: Input all or a sub-set of the channels that are '  ...
             'represented in your TimeFrequency*.mat files.'          ] )

end

% Selected channel names
channelNames  = { ChannelCoordinates(iDecompositionChannels(iChannels)).labels };

% Store channel names and co-ordinates in struct
for w = 1:nCentres
    LocalSpectra.(eventCentres{w}).Channels           = channelNames;
    LocalSpectra.(eventCentres{w}).ChannelCoordinates = ChannelCoordinates;
end

% Channel weighting
if ~weighting
    weights                 = 1;
elseif weighting
    [ weights, oneSDratio ] = electrodeClusterWeights( iChannels );
end


%% Display information
% -------------------------------------------------------------------------

% Files
disp( 'Extracting localised event-related spectra' )
disp( [ 'Windowed relative to ' num2str( nCentres ) ' experimental events' ] )
disp( [ 'From ' num2str( nFiles ) ' time-frequency decompositions' ] )
if nConditions == 1
    conditionPlurality = '';
else
    conditionPlurality = 's';
end
disp( [ 'Representing ' num2str( N ) ' participants and ' num2str( nConditions ) ' condition' conditionPlurality ] )


% Channels
% -------------------------------------------------------------------------

% One selected
if nChannels == 1
    disp( [ 'Localised to channel ' channelNames{1} ] )

% Multiple or all selected
elseif nChannels > 1

    % Multiple selected
    if ~isempty( channels )
        displayText = [ num2str( nChannels ) ' channels as a cluster:' ];

    % All selected
    else
        displayText = [ 'all ' num2str( nChannels ) ' channels in the decompositions:' ];

    end

    % Build list of channel names to display
    channelsText = [];
    for ch = 1:nChannels
        channelsText = [ channelsText channelNames{ch} ' ' ];                                       %#ok
    end

    % Display channels
    if ~weighting
        disp( [ 'Localised to the average of ' displayText ] )
        disp( channelsText )
    elseif weighting
        disp( [ 'Localised to the weighted average of ' displayText ] )
        disp( channelsText )
        disp( [ 'With ' channelNames{1} ' full-weight and other electrodes weighted at ' num2str( oneSDratio*100 ) '%' ] )
    end

end

disp( ' ' )

clear Decomposition


%% Construct time-frequency arrays
% -------------------------------------------------------------------------

% Pre-allocate
impute = false;
Imputations(nFiles).File                 = [];
Imputations(nFiles).AffectedEventWindows = [];
Imputations(nFiles).ParticipantIndex     = [];
Imputations(nFiles).ConditionIndex       = [];
for w = 1:nCentres
    LocalSpectra.(eventCentres{w}).SpectralPower            = NaN( N, nConditions, nFrequencies, nTimes.(eventCentres{w}) );
    LocalSpectra.(eventCentres{w}).SpectralPowerUnits       = spectralPowerUnits;
    LocalSpectra.(eventCentres{w}).SpectralPowerDimensions  = 'Participants x Conditions x Frequencies x Times';
    LocalSpectra.(eventCentres{w}).PhaseCoherence           = NaN( N, nConditions, nFrequencies, nTimes.(eventCentres{w}) );
    LocalSpectra.(eventCentres{w}).PhaseCoherenceUnits      = phaseCoherenceUnits;
    LocalSpectra.(eventCentres{w}).PhaseCoherenceDimensions = 'Participants x Conditions x Frequencies x Times';
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

        % Current window centre
        currentCentre  = eventCentres{w};

        % Spectral power and phase coherence at the selected channels
        % -----------------------------------------------------------------

        % Check existence
        if    isempty( Decomposition.(currentCentre).SpectralPower  )             ...
           || isempty( Decomposition.(currentCentre).PhaseCoherence )             ...
           || isempty( Decomposition.(currentCentre).Coefficients   )             ...
           || all( isnan( Decomposition.(currentCentre).SpectralPower  ), 'all' ) ...
           || all( isnan( Decomposition.(currentCentre).PhaseCoherence ), 'all' ) ...
           || all( isnan( Decomposition.(currentCentre).Coefficients   ), 'all' )
            currentExistence = false;
        else
            currentExistence = true;
        end

        % Trial count
        if currentExistence
            trialCount = size( Decomposition.(currentCentre).Coefficients, 1 );
        else
            trialCount = 0;
        end

        % Data exists in the decomposition
        if currentExistence && trialCount >= minimumTrialCount

            % Copy data from the current decomposition
            if nChannels == 1
                spectralPower(1,:,:)  = Decomposition.(currentCentre).SpectralPower;                   % Frequencies x times
                phaseCoherence(1,:,:) = Decomposition.(currentCentre).PhaseCoherence;                  % Frequencies x times
            elseif nChannels > 1
                spectralPower         = Decomposition.(currentCentre).SpectralPower(iChannels,:,:);    % Channels x frequencies x times
                phaseCoherence        = Decomposition.(currentCentre).PhaseCoherence(iChannels,:,:);   % Channels x frequencies x times
            end

        % Imputation: Set data with too-few trials or no values to NaN so
        %   that the initial grand average is true for good included data
        %   and for later replacement by the initial grand average
        else

            % Imputation is needed for this study
            impute = true;

            % Set event-related spectra to NaN
            spectralPower  = NaN( nChannels, nFrequencies, nTimes.(currentCentre) );
            phaseCoherence = NaN( nChannels, nFrequencies, nTimes.(currentCentre) );
            disp( [ '[' 8 'Warning: ' currentCentre '-related data in ' currentFile ' contains too-few or no trials' ']' 8 ] )
        
            % Store the information needed to impute the correct data
            Imputations(f).File                    = currentFile; 
            Imputations(f).AffectedEventWindows(w) = true;
            Imputations(f).ParticipantIndex        = iParticipant;
            Imputations(f).ConditionIndex          = iCondition;
        
        end
        
        % Average across channels
        spectralPower  = mean( spectralPower  .* weights, 1, 'omitnan' ); % An entirety of NaNs will still average to NaN
        phaseCoherence = mean( phaseCoherence .* weights, 1, 'omitnan' );

        % Clean up
        spectralPower  = squeeze( spectralPower  );
        phaseCoherence = squeeze( phaseCoherence );

        % Store in arrays
        LocalSpectra.(currentCentre).SpectralPower(iParticipant,iCondition,:,:)  = spectralPower;
        LocalSpectra.(currentCentre).PhaseCoherence(iParticipant,iCondition,:,:) = phaseCoherence;

        clear spectralPower phaseCoherence
        
    end % for Window centres

end % for Files


%% Grand average spectra per condition (initial Grand Average if imputing)
% -------------------------------------------------------------------------

% Loop through: Window centres
for w = 1:nCentres

    currentCentre = eventCentres{w};

    % Spectral Power
    LocalSpectra.(currentCentre).GrandAverage.SpectralPower                 ...
         = mean( LocalSpectra.(currentCentre).SpectralPower,  1, 'omitnan' );
    LocalSpectra.(currentCentre).GrandAverage.SpectralPowerUnits            ...
         = spectralPowerUnits;
    if nConditions == 1
        LocalSpectra.(currentCentre).GrandAverage.SpectralPowerDimensions   ...
         = 'Frequencies x Times';
    elseif nConditions > 1
        LocalSpectra.(currentCentre).GrandAverage.SpectralPowerDimensions   ...
         = 'Conditions x Frequencies x Times';
    end

    % Phase Coherence
    LocalSpectra.(currentCentre).GrandAverage.PhaseCoherence                ...
         = mean( LocalSpectra.(currentCentre).PhaseCoherence, 1, 'omitnan' );
    LocalSpectra.(currentCentre).GrandAverage.PhaseCoherenceUnits           ...
         = phaseCoherenceUnits;
    if nConditions == 1
        LocalSpectra.(currentCentre).GrandAverage.PhaseCoherenceDimensions  ...
         = 'Frequencies x Times';
    elseif nConditions > 1
        LocalSpectra.(currentCentre).GrandAverage.PhaseCoherenceDimensions  ...
         = 'Conditions x Frequencies x Times';
    end

    % Clean up
    LocalSpectra.(currentCentre).GrandAverage.SpectralPower  ...
         = squeeze( LocalSpectra.(currentCentre).GrandAverage.SpectralPower  );
    LocalSpectra.(currentCentre).GrandAverage.PhaseCoherence ...
         = squeeze( LocalSpectra.(currentCentre).GrandAverage.PhaseCoherence );

end % for Window centres


%% Imputation
% -------------------------------------------------------------------------

% Check imputation is needed for this study
if impute

    % Loop through: Participant x condition files
    for f = 1:nFiles
    
        % Check imputation is needed for the current file
        if     ~isempty( Imputations(f).File )                 ...
            && ~isempty( Imputations(f).AffectedEventWindows ) ...
            &&      any( Imputations(f).AffectedEventWindows )
    
            % Affected participant and condition indices in LocalSpectra
            % for the current file
            iParticipant = Imputations(f).ParticipantIndex;
            iCondition   = Imputations(f).ConditionIndex;
    
            % Loop through: Event-related windows
            for w = 1:nCentres
    
                % Check imputation is needed for the current window
                if Imputations(f).AffectedEventWindows(w)
    
                    % Affected event-related window
                    affectedCentre = eventCentres{w};


                    % Impute missing or insufficient data
                    % -----------------------------------------------------

                    disp( [ '[' 8 'Imputing ' affectedCentre '-related data for the ' ...
                                  iptnum2ordinal(iParticipant) ' participant in the ' ...
                                  iptnum2ordinal(iCondition)   ' condition'           ']' 8 ] )
                    
                    % Replace individual data with the initial grand average
                    if nConditions > 1

                        LocalSpectra.(affectedCentre).SpectralPower(iParticipant,iCondition,:,:)  ...
             = squeeze( LocalSpectra.(affectedCentre).GrandAverage.SpectralPower(iCondition,:,:)  );

                        LocalSpectra.(affectedCentre).PhaseCoherence(iParticipant,iCondition,:,:) ...
             = squeeze( LocalSpectra.(affectedCentre).GrandAverage.PhaseCoherence(iCondition,:,:) );

                    else

                        LocalSpectra.(affectedCentre).SpectralPower(iParticipant,iCondition,:,:)  ...
                      = LocalSpectra.(affectedCentre).GrandAverage.SpectralPower;

                        LocalSpectra.(affectedCentre).PhaseCoherence(iParticipant,iCondition,:,:) ...
                      = LocalSpectra.(affectedCentre).GrandAverage.PhaseCoherence;

                    end


                end % if Impute this window
    
            end % for Event-related windows
    
        end % if Impute this file
    
    end % for Files


    % Grand average spectra per condition (complete imputed)
    % ---------------------------------------------------------------------
    
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
    
    end % for Window centres


end % Imputation


%% Save
% -------------------------------------------------------------------------

% Location naming
locationName = locationNaming( iChannels,              ...
                               channelNames,           ...
                               iDecompositionChannels, ...
                               ChannelCoordinates      );

% File name
fileName = [ 'TimeFrequency' locationName ];

% Save TimeFrequency<Location>.mat 
save( fileName, '-struct', 'LocalSpectra', '-v7.3' )

% Finish time
codeRunTime( 'Data generation completed at' )


% _________________________________________________________________________
end



%%
% •.° Electrode Cluster Weighting °.•
% _________________________________________________________________________
function [ w, oneSDratio ]  = electrodeClusterWeights( electrodeCluster )

% Number of selected electrodes
nElectrodes = length( electrodeCluster );

% Default multipes of the standard deviation (SD)
if ~exist( 'SD', 'var' ) || nargin < 2
    SD = 1 + 0.05 * max( 0, nElectrodes - 3 );
end

% Weighting as a standard Gaussian 
peakGaussian     = exp( -0^2 / 2 ) / sqrt( 2*pi );

% Primary electrode weighted 1
w(1)             = 1;

% Adjacent electrodes weighted at 1+ SD in a Gaussian
oneSDGaussian    = exp( -SD^2 / 2 ) / sqrt( 2*pi );
oneSDratio       = oneSDGaussian / peakGaussian;
w(2:nElectrodes) = oneSDratio;

% Normalise weights
w = w ./ sum( w ) * nElectrodes;

% Channels in the first dimension of w
if size( w, 1 ) == 1
    w = w';
end


% _________________________________________________________________________
end



%%
% •.° Location Naming °.•
% _________________________________________________________________________
function ...
locationName = locationNaming( iChannels, channelNames, ...
                               iDecompositionChannels, ChannelCoordinates )

% Important channel names
Fz  = 'Fz';
FCz = 'FCz';
CPz = 'CPz';
Pz  = 'Pz';
POz = 'POz';
Oz  = 'Oz';
P7  = 'P7';
P8  = 'P8';
PO7 = 'PO7';
PO8 = 'PO8';

% EGI channels
if any( contains( { ChannelCoordinates(:).labels }, 'E' ) )
    Fz  = 'E11';
    FCz = 'E6';
    CPz = 'E55';
    Pz  = 'E62';
    POz = 'E72';
    Oz  = 'E75';
    P7  = 'E58';
    P8  = 'E96';
    PO7 = 'E65';
    PO8 = 'E90';
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
iPO7 = find( strcmp( { ChannelCoordinates(:).labels }, PO7 ) );
iPO8 = find( strcmp( { ChannelCoordinates(:).labels }, PO8 ) );

% Cluster of electrodes name
if length( iChannels ) > 1 || ( isempty( iChannels ) && length( iDecompositionChannels ) <= 7 )

    iSelectedChannels = iDecompositionChannels(iChannels);

    % Frontal midline cluster including FCz & Fz
    if any( iSelectedChannels == iFCz ) && any( iSelectedChannels == iFz )
        locationName = 'FMidCluster';

    % Frontal midline cluster including FCz
    elseif any( iSelectedChannels == iFCz ) && ~any( iSelectedChannels == iFz )
        locationName = 'FCzCluster';

    % Frontal midline cluster including Fz
    elseif any( iSelectedChannels == iFz ) && ~any( iSelectedChannels == iFCz )
        locationName = 'FzCluster';

    % Parietal midline cluster including Pz & CPz
    elseif any( iSelectedChannels == iPz ) && ~any( iSelectedChannels == iCPz )
        locationName = 'PcMidCluster';

    % Parietal midline cluster including Pz
    elseif any( iSelectedChannels == iPz ) && ~any( iSelectedChannels == iCPz ) && ~any( iSelectedChannels == iPOz )
        locationName = 'PzCluster';

    % Parietal midline cluster including Pz & POz
    elseif any( iSelectedChannels == iPz ) && any( iSelectedChannels == iPOz )
        locationName = 'PoMidCluster';

    % Occipitoparietal midline cluster including POz
    elseif any( iSelectedChannels == iPOz ) && ~any( iSelectedChannels == iPz ) && ~any( iSelectedChannels == iOz )
        locationName = 'POzCluster';

    % Occipital midline cluster including Oz & POz
    elseif any( iSelectedChannels == iOz ) && any( iSelectedChannels == iPOz )
        locationName = 'OMidCluster';

    % Occipital midline cluster including Oz
    elseif any( iSelectedChannels == iOz ) && ~any( iSelectedChannels == iPOz )
        locationName = 'OzCluster';

    % Bilateral parietal cluster including P7 & P8
    elseif any( iSelectedChannels == iP7 ) && any( iSelectedChannels == iP8 )
        locationName = 'PBiCluster';

    % Occipitoparietal bilateral cluster including PO7 & PO8
    elseif any( iSelectedChannels == iPO7 ) && any( iSelectedChannels == iPO8 )
        locationName = 'OPBiCluster';

    % Individual electrodes named Channel1Channel2...ChannelN
    else
        locationName = char( channelNames' );

    end

% Single electrode name
elseif length( iChannels ) == 1
    locationName = channelNames{1};

% All channels
else
    if length( iDecompositionChannels ) <= 10
        locationName = 'Cluster';
    else
        locationName = 'Average';
    end

end % if Number of selected electrodes


% _________________________________________________________________________
end



%%
% •.° Code Run Time °.•
% _________________________________________________________________________
function codeRunTime( message )

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


