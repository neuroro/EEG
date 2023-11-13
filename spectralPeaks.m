function ...
spectralPeaks( fileNamePart, method, frequencyBand, timeLimit, peakSize )
%
% •.° Find Oscillatory Peaks and Sinks °.•
% _________________________________________________________________________
%
% Find oscillatory bursts or suppression in neural event-related spectra
% 
% Extract event-related spectral peaks or sinks within a frequency band per
% event-related window, per condition, per participant, and in the grand
% average(s), enumerated in magnitude as the mean value in the peak or sink
% neighbourhood and by the frequency and time of the peak or sink
%
% Save all oscillatory extrema in a .mat file and save separate .csv files 
% per metric per enumeration per event-related window
%
% • Usage •
% ------------------------------------------------------------------------- 
% Set the Current Folder to the location of the localised event-related 
% spectra file TimeFrequency<Cluster>.mat
% 
% >> spectralPeaks( fileNamePart, method, frequencyBand, timeLimit, peakSize )
% 
% For example:
% >> spectralPeaks
% >> spectralPeaks( 'FMCluster' )
% >> spectralPeaks( '', 'localmax', 'theta', 600 )
% >> spectralPeaks( '', 'min', [2.5 8.75], [-200 300; 75 600] )
% >> spectralPeaks( '', '', '', 0, [1.5 20] )
%
% .. . .  .   .     .        .             .                     .                                  .
%
% •••( Function Inputs )
%
%   fileNamePart:  Name of the localised time-frequency decomposition .mat
%                   file generated by eegLocalSpectra.m, or a unique part
%                   of the file name, for example 'FMCluster'
%                   (optional input: default TimeFrequency*Cluster.mat)
%
%   method:        Peak finding method or sink finding method
%                    'localmax' local maximum with max prominence in 2-D
%                    'localmin' local minimum with max prominence in 2-D
%                    'max'      max of the max across times per frequency
%                    'min'      min of the min across times per frequency
%                   (optional input: default local maxima)
%
%   frequencyBand: Frequency band (named) or frequency range (limits)
%                   within which to find peaks or sinks
%                    'Delta'        1 to <4 Hz
%                    'Theta'        4 to <8 Hz
%                    'Theta+'     2.5 to 8.5 Hz
%                    'Theta++'      2 to 10 Hz
%                    'Alpha'        8 to <13 Hz
%                    'Beta'        13 to 30 Hz
%                    'Gamma'      >30 to 80 Hz
%                    'Low Gamma'  >30 to 50 Hz
%                    'High Gamma' >50 to 80 Hz
%                    [f1 f2]      f1 to f2 Hz frequency range
%                   (optional input: default theta extended to 2.5-8.5 Hz)
%
%   timeLimit:     Time limit or limits (in ms) to find peaks within as a
%                   scalar of the maximum time for stimulus-related spectra
%                    with other time limits dervied in proportion to it; a
%                   [t1 t2] vector of stimulus-related time limits with
%                    other time limits dervied in proportion to t2; or a
%                   w x t matrix of time limits (t) for each event-related
%                    window (w), semicolon-separated in alphbetical order, 
%                    for example [t1Resp t2Resp; t1Stim t2Stim]
%                   (optional input: default 75-600 ms after the stimulus, 
%                    and -third to half the maximum around the response)
%
%   peakSize:      Size of the peak or sink mean magnitude neighbourhood in
%                   frequencies (in Hz) and times (in ms) on either side of
%                   the peak as [f t]
%                   (optional input: default 1.5 Hz 20 ms)
%
% .. . .  .   .     .        .             .                     .                                  .
%
% [ Function Outputs ] =
%
%   Separate .csv files per metric per enumeration per event-related window
%   and a .mat file containing all peaks or sinks data
%
% • Author •
% -------------------------------------------------------------------------
% Rohan O. C. King, 2023
%
% GitHub @neuroro
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
% King, R. (2023). spectralPeaks [MATLAB code]. GitHub. 
% https://github.com/neuroro/EEG/spectralPeaks.m


% Introduction
disp( ' ' )
disp( '•.° Find Oscillatory Peaks and Sinks °.•' )
disp( '_________________________________________________________________________' )
disp( ' ' )
disp( 'Find oscillatory bursts or suppression in neural event-related spectra'    )
disp( ' ' )
disp( 'Extract event-related spectral peaks or sinks within a frequency band per' )
disp( 'event-related window, per condition, per participant, and in the grand'    )
disp( 'average(s), enumerated in magnitude as the mean value in the peak or sink' )
disp( 'neighbourhood and by the frequency and time of the peak or sink'           )
disp( ' ' )


%% Inputs and defaults
% -------------------------------------------------------------------------

% Local decomposition
if ~exist( 'fileNamePart', 'var' ) || isempty( fileNamePart )
    fileNamePart = 'Cluster';
end
if ~ischar( fileNamePart )
    fileNamePart = char( fileNamePart ); % Convert to character if needed
end

% Peak finding method
if ~exist( 'method', 'var' ) || ~contains( method, { 'Max' 'Min' 'Local' }, 'IgnoreCase', true )
    method   = 'LocalMax';
end
method       = lower( method );
if contains( method, 'max' )
    extremum = 'Peak';
    minima   = false;
elseif contains( method, 'min' )
    extremum = 'Sink';
    minima   = true;
end

% Frequency band
if ~exist( 'frequencyBand', 'var' ) || isempty( frequencyBand )
    frequencyBand = 'ThetaExtended';
end
if isstring( frequencyBand )
    frequencyBand = char( frequencyBand );
end

% Time limit
if ~exist( 'timeLimit', 'var' ) || ( isscalar( timeLimit ) && ~timeLimit )
    timeLimit = 600;
end

% Peak size
if ~exist( 'peakSize', 'var' ) || isempty( peakSize )
    peakSize  = [1.5 20];
end


%% Localised time-frequency decomposition
% -------------------------------------------------------------------------

% Wildcard file name
wildSpectra  = [ '*TimeFrequency*' fileNamePart '*.mat' ];

% Search for fileName.mat files named in common located in the Current
% Folder and sub-folders of the Current Folder
fileStruct   = dir( [ '**/' wildSpectra ] );
nFiles       = length( fileStruct );

% Sanity check
if nFiles > 1
    disp( [ '[' 8 'WARNING: Multiple files found, using the first one.' ']' 8 ] )
    disp( ' ' )
end

% File to use
fileName     = fileStruct(1).name;
folderpath   = fileStruct(1).folder;
fullFilePath = fullfile( folderpath, fileName );

% Cluster
if contains( fileName, 'TimeFrequency' ) && contains( fileName, 'Cluster' )
    cluster  = extractAfter( extractBefore( fileName, 'Cluster.mat' ), 'TimeFrequency' );
else
    cluster  = 'Neural';
end

% Load localised event-related spectra
LocalSpectra = load( fullFilePath );


%% Decomposition properties
% -------------------------------------------------------------------------

% Event-centred windowings
eventCentres   = fieldnames( LocalSpectra );
nCentres       = length( eventCentres );

% Time-frequency metrics
grandFields    = fieldnames( LocalSpectra.(eventCentres{1}).GrandAverage );
iMiscFields    = contains( grandFields, { 'Dimension' 'Unit' }, 'IgnoreCase', true );
metrics        = grandFields(~iMiscFields);
nMetrics       = length( metrics );
iPower         = and( contains( metrics, 'Power', 'IgnoreCase', true ), ~contains( metrics, 'Baseline', 'IgnoreCase', true ) );
iCoherence     = contains( metrics, 'Coherence', 'IgnoreCase', true );

% Conditions
conditions     = LocalSpectra.(eventCentres{1}).Conditions;
nConditions    = length( conditions );

% Sample size N participants
N = size( LocalSpectra.(eventCentres{1}).SpectralPower, 1 ); % Participants x conditions x frequencies x times

% Frequencies in the decomposition
frequenciesAll = LocalSpectra.(eventCentres{1}).Frequencies;

% Metric units
metricUnits    = cell( 1, nMetrics );
for m = 1:nMetrics
    try
        metricUnits{m} = LocalSpectra.(eventCentres{1}).GrandAverage.([ metrics{m} 'Units' ]);
    catch
        if m == find( iPower )
            metricUnits{m} = 'Decibels';
        elseif m == find( iCoherence )
            metricUnits{m} = 'Proportion';
        end
    end
end

% Short metric names
measures             = metrics;
measures(iPower)     = { 'Power'     };
measures(iCoherence) = { 'Coherence' };
for m = 1:nMetrics
    if ~minima
        measures{m}  = [ extremum measures{m} ];
    else
        measures{m}  = [ measures{m} extremum ];
    end
end


%% Peak or sink finding limits
% -------------------------------------------------------------------------

% Frequency band limits
[ frequencyBand, bandName ] = aprioriFrequencyBands( frequencyBand );

% Frequency band indices
% Index of the frequency closest to each frequency band limit
[ ~, iBand1 ]   = min( abs( frequenciesAll - frequencyBand(1) ) );
[ ~, iBand2 ]   = min( abs( frequenciesAll - frequencyBand(2) ) );

% Frequencies in the frequency band
frequenciesBand = frequenciesAll(iBand1:iBand2);
nFrequencies    = length( frequenciesBand );

% Time limit scalar or vector of time limits
if isscalar( timeLimit ) || isvector( timeLimit )
    
    timeLimits = zeros( nCentres, 2 );

    try
        earliestTime = timeLimit(1);
        maximumTime  = timeLimit(2);
    catch
        earliestTime = max( 0, 90 - timeLimit/40 );
        maximumTime  = timeLimit;
    end

    for w = 1:nCentres

        % Stimulus-related window
        if contains( eventCentres{w}, 'Stim', 'IgnoreCase', true )
            timeLimits(w,:) = [ earliestTime maximumTime ];

        % Initiation-related window
        elseif contains( eventCentres{w}, 'Init', 'IgnoreCase', true )
            timeLimits(w,:) = [ -maximumTime/2 maximumTime/3 ];

        % Response-related window
        elseif contains( eventCentres{w}, 'Resp', 'IgnoreCase', true )
            timeLimits(w,:) = [ -maximumTime/3 maximumTime/2 ];

        % General event-related window
        else
            timeLimits(w,:) = [ -maximumTime/2 maximumTime ];

        end

    end

% Matrix of time limits
elseif size( timeLimit, 1 ) == nCentres && size( timeLimit, 1 ) == 2
    timeLimits = timeLimit;

else
    errorText  = [ 'Input time limits (in ms) to find peaks within as a '     ...
                   'scalar of the maximum time for stimulus-related spectra ' ...
                   'with other time limits dervied in proportion to it; a '   ...
                   '[t1 t2] vector of stimulus-related time limits with '     ...
                   'other time limits dervied in proportion to t2; or a '     ...
                   'w x t matrix of time limits (t) for each event-related '  ...
                   'window (w), semicolon-separated in alphbetical order, '   ...
                   'for example [t1Resp t2Resp; t1Stim t2Stim]'               ];
    error( errorText )

end


%% Display parameters
% -------------------------------------------------------------------------

disp( '• Conditions •' )
disp( '-------------------------------------------------------------------------' )
disp( 'Peaks or sinks are extracted per condition in order' )
disp( char( conditions ) )
disp( ' ' )
disp( '• Event-centred windows •' )
disp( '-------------------------------------------------------------------------' )
disp( char( eventCentres ) )
disp( ' ' )
disp( '• Time-frequency metrics •' )
disp( '-------------------------------------------------------------------------' )
disp( char( metrics ) )
disp( ' ' )
disp( '• Oscillatory peaks or sinks •' )
disp( '-------------------------------------------------------------------------' )
disp( [ extremum 's are extracted within frequencies of ' num2str( frequencyBand(1) ) ' - ' num2str( frequencyBand(2) ) ' Hz'] )
disp( 'Times in ms are limited per window centre in order' )
disp( num2str( timeLimits ) )
if contains( method, 'local' )
    if ~minima
        disp( 'Peaks are 2-D joint local maxima with the largest joint prominence' )
    else
        disp( 'Sinks are 2-D joint local minima with the largest joint prominence' )
    end
else
    if ~minima
        disp( 'Peaks are the maximum of the maximum across times per frequency' )
    else
        disp( 'Sinks are the minimum of the minimum across times per frequency' )
    end
end
disp( ' ' )


%% Set up
% -------------------------------------------------------------------------

% SpectralPeaks struct peak field names and field pre-allocation
dimensions      = { '' 'Frequency' 'Time'         };
dimensionsUnits = { '' 'Hz'        'milliseconds' };
for w = 1:nCentres
    for m = 1:nMetrics
        dimensions{1} = metrics{m};
        if ~isempty( metricUnits{m} )
            dimensionsUnits{1} = metricUnits{m};
        end
        for d = 1:length( dimensions )
            if ~minima
                peakFields{d} = [ extremum dimensions{d} ];                 %#ok
            else
                peakFields{d} = [ dimensions{d} extremum ];                 %#ok
            end
            SpectralPeaks.(eventCentres{w}).(metrics{m}).(peakFields{d}) = NaN( N, nConditions );
            SpectralPeaks.(eventCentres{w}).(metrics{m}).([peakFields{d} 'Units']) = dimensionsUnits{d};
        end
    end
end

% Loop through: Event-related windows
for w = 1:nCentres

    % Current event centre
    currentCentre = eventCentres{w};

    % Current event-centred times
    timesAll      = LocalSpectra.(eventCentres{w}).Times;
    
    % Current event-centred time limit indices
    [ ~, iTime1 ] = min( abs( timesAll - timeLimits(w,1) ) );
    [ ~, iTime2 ] = min( abs( timesAll - timeLimits(w,2) ) );

    % Current event-centred window times
    timesLimited  = timesAll(iTime1:iTime2);

    % Loop through: Time-Frequency Metrics
    for m = 1:nMetrics

        % Current time-frequency metric
        currentMetric = metrics{m};

        % Metric fields
        dimensions    = { '' 'Frequency' 'Time' };
        dimensions{1} = metrics{m};
        for d = 1:length( dimensions )
            if ~minima
                peakFields{d} = [ extremum dimensions{d} ];
            else
                peakFields{d} = [ dimensions{d} extremum ];
            end
        end

        % Loop through: Conditions
        for c = 1:nConditions

            % Loop through: Participants
            % Participant 0 is the grand average
            for p = 0:N


                %% Time-frequency window in which to find a peak or sink
                % ---------------------------------------------------------

                % Current participant x condition time-frequency window
                if p == 0
                    timeFrequency   = LocalSpectra.(currentCentre).GrandAverage.(currentMetric);
                    if nConditions > 1
                    timeFrequency   = timeFrequency(c,:,:);
                    end
                elseif p > 0
                    timeFrequency   = LocalSpectra.(currentCentre).(currentMetric);
                    timeFrequency   = timeFrequency(p,c,:,:);
                end
                timeFrequency       = squeeze( timeFrequency );

                % Apply the time limits
                timeFrequencyWindow = timeFrequency(:,iTime1:iTime2);

                % Select the frequency band
                timeFrequencyWindow = timeFrequencyWindow(iBand1:iBand2,:);


                %% Find the peak or sink
                % ---------------------------------------------------------

                % Sink minimum sign correction: min = -max( -data )
                % Find minima using max( -data ) followed by -max below
                if minima
                    timeFrequencyWindow = -timeFrequencyWindow;
                end

                % Most prominent local maximum across frequencies and times
                % ---------------------------------------------------------
                if contains( method, 'local' )

                    % Metric local maxima across frequencies and times
                    [ fLocalMaxima, fProminences ] = islocalmax( timeFrequencyWindow, 1 ); % 2-D max indices and prominence values
                    [ tLocalMaxima, tProminences ] = islocalmax( timeFrequencyWindow, 2 ); % 2-D max indices and prominence values

                    % Frequency jointness tolerance
                    %   0-2 frequencies on either side of the max are
                    %   considered joint depending on frequency resolution
                    fTolerance1 = 0;
                    fTolerance2 = 2;
                    fTolerance  = max( fTolerance1, round( nFrequencies*0.05 ) ); % 5% of the number of frequencies in the frequency band
                    fTolerance  = min( fTolerance2, fTolerance );

                    % Time jointness tolerance
                    %   1-3+ samples on either side of the max are
                    %   considered joint depending on the sampling rate
                    sampleTime  = mean( diff( timesLimited ) );
                    tTolerance  = max( 1, round( 2 / sampleTime ) );


                    % Joint local maxima across frequencies and times with
                    % jointness tolerance
                    % -----------------------------------------------------
                    iLocalMaxima = false( size( timeFrequencyWindow ) );
                    nTimes       = length( timesLimited );

                    % Loop through: Frequencies and times
                    for f = 1:nFrequencies
                        for t = 1:nTimes

                            % Test if the current time-frequency point is
                            % locally maximal across frequencies or times
                            if fLocalMaxima(f,t) || tLocalMaxima(f,t)

                                % Tolerated frequency indices
                                fRange1 = max( f - fTolerance, 1 );
                                fRange2 = min( f + fTolerance, nFrequencies );
                                fRange  = fRange1:fRange2;

                                % Tolerated time indices
                                tRange1 = max( t - tTolerance, 1 );
                                tRange2 = min( t + tTolerance, nTimes );
                                tRange  = tRange1:tRange2;

                                % Tolerated window around the current local
                                % maximum that is considered joint
                                fWindow = fLocalMaxima(fRange,tRange);
                                tWindow = tLocalMaxima(fRange,tRange);

                                % Test if the current local maximum is a
                                % joint time-frequency local maximum
                                if any( fWindow, 'all' ) && any( tWindow, 'all' )
                                    iLocalMaxima(f,t) = true;
                                end

                            end

                        end
                    end

                    % No joint local maximum anywhere
                    if ~any( iLocalMaxima, 'all' )

                        % Warning
                        if ~minima
                            minOrMax = 'maximum';
                        else
                            minOrMax = 'minimum';
                        end
                        disp( [ '[' 8 'Warning: No local ' minOrMax ' in ' currentMetric ' found jointly in 2-D'             ']' 8 ] )
                        if p
                            disp( [ '[' 8 'for the ' iptnum2ordinal(p) ' participant in the ' iptnum2ordinal(c) ' condition' ']' 8 ] )
                        else
                            disp( [ '[' 8 'for the grand average of the ' iptnum2ordinal(c) ' condition'                     ']' 8 ] )
                        end
                        disp( [ '[' 8 'Increasing the neighbourhood aroud each ' minOrMax ' considered singularly joint'     ']' 8 ] )
                        disp( ' ' )

                        % Increase the tolerances
                        fTolerance = 2 * fTolerance + 1;
                        tTolerance = 2 * tTolerance + 1;

                        % Loop through: Frequencies and times (again)
                        for f = 1:nFrequencies
                            for t = 1:nTimes
    
                                % Test if the current time-frequency point
                                % is locally maximal in time and frequency 
                                if fLocalMaxima(f,t) || tLocalMaxima(f,t)
    
                                    % Tolerated frequency indices
                                    fRange1 = max( f - fTolerance, 1 );
                                    fRange2 = min( f + fTolerance, nFrequencies );
                                    fRange  = fRange1:fRange2;
    
                                    % Tolerated time indices
                                    tRange1 = max( t - tTolerance, 1 );
                                    tRange2 = min( t + tTolerance, nTimes );
                                    tRange  = tRange1:tRange2;
    
                                    % Tolerated window around the current
                                    % local maximum to be considered joint
                                    fWindow = fLocalMaxima(fRange,tRange);
                                    tWindow = tLocalMaxima(fRange,tRange);
    
                                    % Test if the current local maximum is
                                    % a joint time-frequency local maximum
                                    if any( fWindow, 'all' ) && any( tWindow, 'all' )
                                        iLocalMaxima(f,t) = true;
                                    end
    
                                end
    
                            end
                        end

                    end % if No joint local maxima

                    % Still no joint local maximum anywhere
                    if ~any( iLocalMaxima, 'all' )
                        disp( [ '[' 8 'Warning: Still no local ' minOrMax ' in ' currentMetric ' found jointly in 2-D'   ']' 8 ] )
                        if p
                            disp( [ '[' 8 'for the ' iptnum2ordinal(p) ' participant in the ' iptnum2ordinal(c) ' condition' ']' 8 ] )
                        else
                            disp( [ '[' 8 'for the grand average of the ' iptnum2ordinal(c) ' condition' ']' 8 ] )
                        end
                        disp( [ '[' 8 'Using 1-D local ' minOrMax(1:5) 'a across times instead, for this case'           ']' 8 ] )
                        disp( ' ' )
                        iLocalMaxima = tLocalMaxima;                                    % 2-D logical indices
                    end

                    % Joint sum prominences
                    prominences     = fProminences + tProminences;
                    prominentMaxima = prominences .* iLocalMaxima;                      % 2-D prominence of joint maxima

                    % Most prominent local maximum
                    maxProminent    = max( prominentMaxima, [], 'all', 'omitnan' );     % max will give the first of multiple equal max values
                    iLocalMaximum   = prominentMaxima == maxProminent;                  % 2-D logical indices allowing multiple equal maxima

                    % Most powerful local maximum of multiple equally
                    % jointly prominent maxima
                    if sum( iLocalMaximum, 'all' ) > 1
                        localMaxPower = timeFrequencyWindow .* iLocalMaximum;           % 2-D power of joint maxima with max prominence
                        maxPower      = max( localMaxPower, [], 'all', 'omitnan' );     % max will give the first of multiple equal max values
                        iLocalMaximum = localMaxPower == maxPower;                      % 2-D logical indices allowing multiple equal maxima
                    end


                    % Peak frequency and time
                    % -----------------------------------------------------

                    % Peak indices
                    [ iFrequencyPeak, iTimePeak ] = find( iLocalMaximum );
                    if numel( iFrequencyPeak ) > 1
                        iFrequencyPeak            = iFrequencyPeak(1);                  % First in case multiple remain
                    end
                    if numel( iTimePeak ) > 1
                        iTimePeak                 = iTimePeak(1);                       % First in case multiple remain
                    end

                    % Peak frequency
                    peakFrequency = frequenciesBand(iFrequencyPeak);

                    % Peak time
                    peakTime = timesLimited(iTimePeak);


                % Maximum of the maximum across times per frequency
                % ---------------------------------------------------------
                elseif contains( method, { 'max' 'min' } ) && ~contains( method, 'local' )

                    % Pre-allocate
                    outputPeaks = NaN( 1, nFrequencies ); % Peak value for each frequency
                    iTimePeaks  = NaN( 1, nFrequencies ); % Time index of the peak at each frequency

                    % Loop through: Frequencies in the frequency band
                    for f = 1:nFrequencies

                        % Maximum across times at the current frequency 
                        [ outputPeaks(f), iTimePeaks(f) ] = max( timeFrequencyWindow(f,:), [], 'omitnan' );

                    end

                    % Peak frequency index at the grand maximum peak
                    [ ~, iFrequencyPeak ] = max( outputPeaks, [], 'omitnan' );

                    % Peak frequency
                    peakFrequency = frequenciesBand(iFrequencyPeak);

                    % Peak time index for the peak frequency of the grand maximum peak
                    iTimePeak = iTimePeaks(iFrequencyPeak);

                    % Peak time
                    peakTime = timesLimited(iTimePeak);


                end % if method

                % No peak
                if isempty( peakFrequency ) || isempty( peakTime )
                    peakFrequency = NaN;
                    peakTime      = NaN;
                end


                %% Determine the peak window
                % ---------------------------------------------------------

                % Peak window frequency limits
                peakFrequencyLimits         = [ ( peakFrequency - peakSize(1) ) ...
                                                ( peakFrequency + peakSize(1) ) ];

                % Peak window frequency limit indices
                [ ~, iFrequencyPeakLimit1 ] = min( abs( frequenciesAll - peakFrequencyLimits(1) ) );
                [ ~, iFrequencyPeakLimit2 ] = min( abs( frequenciesAll - peakFrequencyLimits(2) ) );
                iFrequenciesPeak            = iFrequencyPeakLimit1:iFrequencyPeakLimit2;

                % Peak window time limits
                peakTimeLimits              = [ ( peakTime - peakSize(2) ) ...
                                                ( peakTime + peakSize(2) ) ];

                % Peak window time limit indices
                [ ~, iTimePeakLimit1 ]      = min( abs( timesAll - peakTimeLimits(1) ) );
                [ ~, iTimePeakLimit2 ]      = min( abs( timesAll - peakTimeLimits(2) ) );
                iTimesPeak                  = iTimePeakLimit1:iTimePeakLimit2;


                %% Extract the mean oscillatory peak or sink
                % ---------------------------------------------------------

                % Select the peak or sink window
                peakWindow = timeFrequency(iFrequenciesPeak,iTimesPeak);

                % Unweighted average across the peak or sink window
                peak       = mean( peakWindow(:), 'omitnan' );

                % Sink minimum sign correction: min = -max( -data )
                if minima
                    peak   = -peak;
                end


                %% Error checks
                % ---------------------------------------------------------

                % Participant error reporting text
                if p == 0
                    participantText = 'grand average';
                else
                    participantText = [ iptnum2ordinal(p) ' participant' ];
                end

                % Extremum frequency outside the frequency range
                if peakFrequency < frequenciesBand(1) || peakFrequency > frequenciesBand(end)
                    disp( [ '[' 8 '•••••••••••••••••••••••••' ']' 8 ] )
                    disp( [ '[' 8 '•      •         •      •' ']' 8 ] )
                    disp( [ '[' 8 '•  :(  •  ERROR  •  :(  •' ']' 8 ] )
                    disp( [ '[' 8 '•      •         •      •' ']' 8 ] )
                    disp( [ '[' 8 '•••••••••••••••••••••••••' ']' 8 ] )
                    disp( [ '[' 8 extremum ' found outside ' num2str( frequenciesBand(1) ) ' Hz - ' num2str( frequenciesBand(end) ) ' Hz' ']' 8 ] )
                    disp( [ '[' 8 'In ' currentMetric ' relative to the ' currentCentre ']' 8 ] )
                    disp( [ '[' 8 'For the ' participantText ']' 8 ] )
                    disp( [ '[' 8 'In the '  iptnum2ordinal(c) ' condition' ']' 8 ] )
                end

                % Extremum time outside the time limits
                if peakTime < timesLimited(1) || peakTime > timesLimited(end)
                    disp( [ '[' 8 '•••••••••••••••••••••••••' ']' 8 ] )
                    disp( [ '[' 8 '•      •         •      •' ']' 8 ] )
                    disp( [ '[' 8 '•  :(  •  ERROR  •  :(  •' ']' 8 ] )
                    disp( [ '[' 8 '•      •         •      •' ']' 8 ] )
                    disp( [ '[' 8 '•••••••••••••••••••••••••' ']' 8 ] )
                    disp( [ '[' 8 extremum ' found outside ' num2str( timesLimited(1) ) ' ms - ' num2str( timesLimited(end) ) ' ms' ']' 8 ] )
                    disp( [ '[' 8 'In ' currentMetric ' relative to the ' currentCentre ']' 8 ] )
                    disp( [ '[' 8 'For the ' participantText ']' 8 ] )
                    disp( [ '[' 8 'In the '  iptnum2ordinal(c) ' condition' ']' 8 ] )
                end

                % Non-existent data
                if all( isnan( timeFrequency ), 'all' )
%                     disp( [ '[' 8 ':S  Missing data  :S' ']' 8 ] )
%                     disp( [ '[' 8 '••••••••••••••••••••' ']' 8 ] )
                    disp( [ '[' 8 'Missing data' ']' 8 ] )
                    disp( [ '[' 8 '••••••••••••' ']' 8 ] )
                    disp( [ '[' 8 'No ' currentMetric ' relative to the ' currentCentre ']' 8 ] )
                    disp( [ '[' 8 'For the ' participantText ']' 8 ] )
                    disp( [ '[' 8 'In the '  iptnum2ordinal(c) ' condition' ']' 8 ] )
                    disp( ' ' )
                    peak          = NaN;
                    peakFrequency = NaN;
                    peakTime      = NaN;

                % Non-numeric extremum
                elseif ~all( isnan( timeFrequency ), 'all' ) && isnan( peak )
                    disp( [ '[' 8 '•••••••••••••••••••••••••••' ']' 8 ] )
                    disp( [ '[' 8 '•      •           •      •' ']' 8 ] )
                    disp( [ '[' 8 '•  :/  •  Warning  •  :/  •' ']' 8 ] )
                    disp( [ '[' 8 '•      •           •      •' ']' 8 ] )
                    disp( [ '[' 8 '•••••••••••••••••••••••••••' ']' 8 ] )
                    disp( [ '[' 8 measures{m} ' relative to the ' currentCentre ' is not-a-number' ']' 8 ] )
                    disp( [ '[' 8 'For the ' participantText ']' 8 ] )
                    disp( [ '[' 8 'In the '  iptnum2ordinal(c) ' condition' ']' 8 ] )
                    disp( [ '[' 8 'This is not due to missing data and may be an error' ']' 8 ] )
                    disp( ' ' )
                    peakFrequency = NaN;
                    peakTime      = NaN;

                end


                %% Store in struct
                % ---------------------------------------------------------

                % Peaks or sinks at the grand average per condition
                if p == 0
                    SpectralPeaks.(currentCentre).(currentMetric).GrandAverage.(peakFields{1})(c)        = peak;
                    if ~isempty( metricUnits{m} )
                    SpectralPeaks.(currentCentre).(currentMetric).GrandAverage.([peakFields{1} 'Units']) = metricUnits{m};
                    end
                    SpectralPeaks.(currentCentre).(currentMetric).GrandAverage.(peakFields{2})(c)        = peakFrequency;
                    SpectralPeaks.(currentCentre).(currentMetric).GrandAverage.([peakFields{2} 'Units']) = 'Hz';
                    SpectralPeaks.(currentCentre).(currentMetric).GrandAverage.(peakFields{3})(c)        = peakTime;
                    SpectralPeaks.(currentCentre).(currentMetric).GrandAverage.([peakFields{3} 'Units']) = 'milliseconds';

                % Peaks or sinks per participant per condition
                elseif p > 0
                    SpectralPeaks.(currentCentre).(currentMetric).(peakFields{1})(p,c) = peak;
                    SpectralPeaks.(currentCentre).(currentMetric).(peakFields{2})(p,c) = peakFrequency;
                    SpectralPeaks.(currentCentre).(currentMetric).(peakFields{3})(p,c) = peakTime;
                    
                end


            end % for: Participants 

        end % for: Conditions

        % Conditions
        SpectralPeaks.(currentCentre).(currentMetric).GrandAverage.Conditions = conditions;
        SpectralPeaks.(currentCentre).(currentMetric).Conditions              = conditions;


%% Save peaks or sinks
% -------------------------------------------------------------------------

        % Export a comma separated matrix for each enumerated dimension of
        % the peak or sink of each oscillatory metric in each window centre
        for d = 1:length( peakFields )
            currentDimension   = peakFields{d};
            if d == 1
                currentMeasure = measures{m};
            else
                currentMeasure = currentDimension;
            end
            outputPeaks        = SpectralPeaks.(currentCentre).(currentMetric).(currentDimension);
            outputFileName     = [ cluster bandName currentMeasure 'In' currentCentre 'Related' currentMetric '.csv' ];
            writematrix( outputPeaks, outputFileName )
        end


    end % for: Time-frequency metrics

end % for: Event-related windows


% Save matlab matrix
save( [ cluster bandName extremum 's.mat' ], '-struct', 'SpectralPeaks', '-v7.3' )


% _________________________________________________________________________
end



%% 
% •.° Frequency Band Definitions °.•
% _________________________________________________________________________
%
function [ frequencyLimits, bandName ] = aprioriFrequencyBands( frequencyBand )

% Frequency band names
delta      = { 'Delta' 'D' };
theta      = { 'Theta' 'T' };
alpha      = { 'Alpha' 'A' };
beta       = { 'Beta'  'B' };
gamma      = { 'Gamma' 'G' };
lowGamma   = { 'Low Gamma'  'LowGamma'  'Gamma1' 'G1' };
highGamma  = { 'High Gamma' 'HighGamma' 'Gamma2' 'G2' };

% Organisation for Human Brain Mapping definitions (Pernet et al., 2018)
deltaBand  = [ 1       3.999  ];
thetaBand  = [ 4       7.999  ];
alphaBand  = [ 8       12.999 ];
betaBand   = [ 13      30     ];
gammaBand  = [ 30.001  80     ];

% Gamma-band divisions
gammaBand1 = [ 30.001  50     ];
gammaBand2 = [ 50.001  80     ];

% Extended theta+ peak finding window
thetaX     = { 'Theta+'  'T+'  'Tp'    'ThetaPlus'  'ThetaExtended'  'Tx' };
thetaXF    = [ 2.5 8.5 ];

% 2-10 Hz theta++ peak finding window (Gyurkovics & Levita, 2021) 
thetaX2    = { 'Theta++' 'T++' 'Tpp'   'ThetaPlus2' 'ThetaPlusPlus' ...
               'ThetaGyurkovicsLevita' 'ThetaGL'    'TGL'           };
thetaX2F   = [ 2 10 ];

% Frequency limits
if ischar( frequencyBand )
    switch lower( frequencyBand )
        case lower( delta )
            frequencyLimits = deltaBand;
            bandName        = delta{1};
        case lower( theta )
            frequencyLimits = thetaBand;
            bandName        = theta{1};
        case lower( alpha )
            frequencyLimits = alphaBand;
            bandName        = alpha{1};
        case lower( beta )
            frequencyLimits = betaBand;
            bandName        = beta{1};
        case lower( gamma )
            frequencyLimits = gammaBand;
            bandName        = gamma{1};
        case lower( lowGamma )
            frequencyLimits = gammaBand1;
            bandName        = lowGamma{2};
        case lower( highGamma )
            frequencyLimits = gammaBand2;
            bandName        = highGamma{2};
        case lower( thetaX )
            frequencyLimits = thetaXF;
            bandName        = thetaX{1};
        case lower( thetaX2 )
            frequencyLimits = thetaX2F;
            bandName        = thetaX2{1};
    end
elseif isnumeric( frequencyBand )
            frequencyLimits = frequencyBand;
            bandName        = [ num2str( frequencyBand(1) ) 'Hz-' num2str( frequencyBand(2) ) 'Hz' ];
else
    error( 'Specify a named frequency band or frequency limits in Hz as [f1 f2]' )
end

% • References •
% -------------------------------------------------------------------------
%
% Gyurkovics, M. & Levita, L. (2021). Dynamic Adjustments of Midfrontal
%   Control Signals in Adults and Adolescents. Cerebral Cortex, 31(2), 
%   795–808. https://doi.org/10.1093/cercor/bhaa258
%
% Pernet, C. R., Garrido, M., Gramfort, A., Maurits, N., Michel, C., Pang,
%   E., … Puce, A. (2018, August 9). Best Practices in Data Analysis and
%   Sharing in Neuroimaging using MEEG.
%   https://doi.org/10.31219/osf.io/a8dhx


% _________________________________________________________________________
end



%% Alternative rank sum maxima selection (incomplete)
% % Minimum rank sum of prominence and power (first ordinal)
% powerfulMaxima     = timeFrequencyWindow .* iLocalMaxima;
% maxProminences     = prominentMaxima(prominentMaxima ~= 0);
% maxPowers          = powerfulMaxima(powerfulMaxima ~= 0);
% rankProminence     = !!!!!
% rankPower          = !!!!!
% rankSumProPow      = rankProminence + rankPower;
% [ ~, iMinRankSum ] = min( rankSumProPow );
% maxRankProminent   = maxProminences(iMinRankSum);
% iLocalMaximum      = prominentMaxima == maxRankProminent;    % 2-D logical indices allowing multiple equal maxima