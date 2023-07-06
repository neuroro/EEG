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
% >> spectralPeaks( fileNamePart, method, frequencyBand, timeLimit, peakSize )
%
% Inputs:
%   fileNamePart:  Name of the localised time-frequency decomposition .mat
%                    file generated by eegLocalSpectra.m, or a unique part
%                    of the file name, for example 'FMCluster'
%                    (optional input: default TimeFrequency*Cluster.mat)
%   method:        Peak finding method
%                    'localmax' local maximum with max prominence in 2-D
%                    'localmin' local minimum with max prominence in 2-D
%                    'max'      max of the max across times per frequency
%                    'min'      min of the min across times per frequency
%                    (optional input: default local maxima)
%   frequencyBand: Frequency band to extract peaks from as frequency limits
%                    [f1 f2] in Hz or as pre-set named frequency bands
%                    'delta', 'theta', 'theta+', 'alpha', 'beta', 'gamma'
%                    (optional input: default 2-10 Hz)
%   timeLimit:     Time limit or limits (in ms) within which to extract
%                    peaks or sinks as a single reference time limit from 
%                    which the rest are calculated or as a matrix of time 
%                    limits per (semicolon-separated) event-related window 
%                    (optional input: default 100 ms to 600 ms relative to 
%                     the stimulus and -300 to 200 ms relative to response)
%   peakSize:      Frequencies (in Hz) and times (in ms) on either
%                    side of the peak to average across as [f t]
%                    (optional input: default 1.5 Hz 25 ms)
% Outputs:
%   Separate .csv files per metric per enumeration per event-related window
%   and a .mat file containing all peaks or sinks data
%
% Examples:
% >> spectralPeaks
% >> spectralPeaks( 'FMCluster' )
% >> spectralPeaks( '', 'localmax', [2 10], 600 )
% >> spectralPeaks( '', 'min', 'alpha', [100 600; -300 200] )
% >> spectralPeaks( '', '', '', 0, [1.5 25] )
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
    fileNamePart = 'TimeFrequency*Cluster';
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
    frequencyBand = 'Theta+';
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
    peakSize  = [2 25];
end


%% Localised time-frequency decomposition
% -------------------------------------------------------------------------

% Wildcard file name
wildSpectra  = [ '*' fileNamePart '*.mat' ];

% Search for fileName.mat files named in common located in the Current
% Folder and sub-folders of the Current Folder
fileStruct   = dir( [ '**/' wildSpectra ] );
nFiles       = length( fileStruct );

% Sanity check
if nFiles > 1
    disp( [ '[' 8 'WARNING: Multiple files found, using the first one.' ']' 8 ] )
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
eventCentres = fieldnames( LocalSpectra );
nCentres     = length( eventCentres );

% Time-frequency metrics
grandFields  = fieldnames( LocalSpectra.(eventCentres{1}).GrandAverage );
iMiscFields  = contains( grandFields, { 'Dimension' 'Unit' }, 'IgnoreCase', true );
metrics      = grandFields(~iMiscFields);
nMetrics     = length( metrics );

% Conditions
conditions   = LocalSpectra.(eventCentres{1}).Conditions;
nConditions  = length( conditions );

% Sample size N participants
N = size( LocalSpectra.(eventCentres{1}).SpectralPower, 1 ); % Participants x conditions x frequencies x times

% Frequencies in the decomposition
frequenciesAll = LocalSpectra.(eventCentres{1}).Frequencies;


%% Peak finding limits
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

% Time limits
if isscalar( timeLimit )
    for w = 1:nCentres
        if w == 1
            timeLimits(w,:) = round( [ 0            timeLimit   ] );            %#ok
        else
            timeLimits(w,:) = round( [ -timeLimit/2 timeLimit/3 ] );            %#ok
        end
    end
else
    timeLimits = timeLimit;
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

% Short metric names
iPower               = contains( metrics, 'Power',     'IgnoreCase', true );
iCoherence           = contains( metrics, 'Coherence', 'IgnoreCase', true );
measures             = metrics;
measures(iPower)     = { 'Power'     };
measures(iCoherence) = { 'Coherence' };
for m = 1:nMetrics
    if ~minima
        measures{m} = [ extremum measures{m} ];
    else
        measures{m} = [ measures{m} extremum ];
    end
end

% SpectralPeaks struct peak field names and field pre-allocation
dimensions = { '' 'Frequency' 'Time' };
for w = 1:nCentres
    for m = 1:nMetrics
        dimensions{1} = metrics{m};
        for d = 1:length( dimensions )
            if ~minima
                peakFields{d} = [ extremum dimensions{d} ];                 %#ok
            else
                peakFields{d} = [ dimensions{d} extremum ];                 %#ok
            end
            SpectralPeaks.(eventCentres{w}).(metrics{m}).(peakFields{d}) = NaN( N, nConditions );
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
                    [ fLocalMaxima, fProminences ] = islocalmax( timeFrequencyWindow, 1 );
                    [ tLocalMaxima, tProminences ] = islocalmax( timeFrequencyWindow, 2 );

                    % Frequency jointness tolerance
                    %   0-2 frequencies on either side of the max are
                    %   considered joint depending on frequency resolution 
                    fTolerance1 = 0;
                    fTolerance2 = 2;
                    fTolerance  = max( fTolerance1, round( nFrequencies*0.05 ) );
                    fTolerance  = min( fTolerance2, fTolerance );

                    % Time jointness tolerance
                    %   2-4 samples on either side of the max are
                    %   considered joint depending on the sampling rate
                    if timesLimited(2) - timesLimited(1) < 1      % > 1000 Hz
                        tTolerance = 3;
                    elseif timesLimited(2) - timesLimited(1) == 1 % = 1000 Hz
                        tTolerance = 2;
                    elseif timesLimited(2) - timesLimited(1) > 1  % < 1000 Hz
                        tTolerance = 1;
                    end

                    % Joint local maxima across frequencies and times with
                    % jointness tolerance
                    % -----------------------------------------------------
                    localMaxima = false( size( timeFrequencyWindow ) );
                    nTimes      = length( timesLimited );

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
                                    localMaxima(f,t) = true;
                                end

                            end

                        end
                    end

                    % No joint local maximum anywhere
                    if ~any( localMaxima, 'all' )

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

                        % Double the tolerances
                        fTolerance = fTolerance * 2;
                        tTolerance = tTolerance * 2;

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
                                        localMaxima(f,t) = true;
                                    end
    
                                end
    
                            end
                        end

                    end % if No joint local maxima

                    % Still no joint local maximum anywhere
                    if ~any( localMaxima, 'all' )
                        disp( [ '[' 8 'Warning: Still no local ' minOrMax ' in ' currentMetric ' found jointly in 2-D'   ']' 8 ] )
                        if p
                            disp( [ '[' 8 'for the ' iptnum2ordinal(p) ' participant in the ' iptnum2ordinal(c) ' condition' ']' 8 ] )
                        else
                            disp( [ '[' 8 'for the grand average of the ' iptnum2ordinal(c) ' condition' ']' 8 ] )
                        end
                        disp( [ '[' 8 'Using 1-D local ' minOrMax(1:5) 'a across times instead, for this case'           ']' 8 ] )
                        disp( ' ' )
                        localMaxima = tLocalMaxima;
                    end

                    % Joint sum prominences
                    prominences     = fProminences + tProminences;
                    prominentMaxima = prominences .* localMaxima;

                    % Most prominent local maximum
                    maxProminent    = max( prominentMaxima, [], 'all', 'omitnan' ); % max will give the first of multiple equal max values
                    iLocalMaximum   = prominentMaxima == maxProminent;
                    
                    % Multiple equally jointly prominent maxima
                    % Most powerful local maximum
                    if sum( iLocalMaximum, 'all' ) > 1
                        localMaxPower = timeFrequencyWindow(iLocalMaximum);
                        [ ~, iMax ]   = max( localMaxPower, [], 'omitnan' ); % max will give the first of multiple equal max values
                        iLocalMaximum = iLocalMaximum(iMax);
                    end


                    % Peak frequency and time
                    % -----------------------------------------------------

                    % Peak indices
                    [ iFrequencyPeak, iTimePeak ] = find( iLocalMaximum );

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

                % Sanity check: Non-numeric peak or sink
                if isnan( peak )
                    disp( [ '[' 8 '•••••••••••••••••••••••••' ']' 8 ] )
                    disp( [ '[' 8 '•      •         •      •' ']' 8 ] )
                    disp( [ '[' 8 '•  :(  •  ERROR  •  :(  •' ']' 8 ] )
                    disp( [ '[' 8 '•      •         •      •' ']' 8 ] )
                    disp( [ '[' 8 '•••••••••••••••••••••••••' ']' 8 ] )
                    disp( [ '[' 8 measures{m} ' not-a-number' ']' 8 ] )
                    disp( [ '[' 8 'For the ' iptnum2ordinal(p) ' participant' ']' 8 ] )
                    disp( [ '[' 8 'In the ' iptnum2ordinal(c) ' condition'    ']' 8 ] )
                    disp( ' ' )
                end


                %% Store in struct
                % ---------------------------------------------------------

                % Peaks or sinks at the grand average per condition
                if p == 0
                    SpectralPeaks.(currentCentre).(currentMetric).GrandAverage.(peakFields{1})(c) = peak;
                    SpectralPeaks.(currentCentre).(currentMetric).GrandAverage.(peakFields{2})(c) = peakFrequency;
                    SpectralPeaks.(currentCentre).(currentMetric).GrandAverage.(peakFields{3})(c) = peakTime;

                % Peaks or sinks per participant per condition
                elseif p > 0
                    SpectralPeaks.(currentCentre).(currentMetric).(peakFields{1})(p,c) = peak;
                    SpectralPeaks.(currentCentre).(currentMetric).(peakFields{2})(p,c) = peakFrequency;
                    SpectralPeaks.(currentCentre).(currentMetric).(peakFields{3})(p,c) = peakTime;
                    
                end


            end % for: Participants 

        end % for: Conditions


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
save( [ cluster bandName extremum 's' ], '-struct', 'SpectralPeaks', '-v7.3' )


% _________________________________________________________________________
end



%% 
% •.° Frequency Band Definitions °.•
% _________________________________________________________________________
%
function [ frequencyLimits, bandName ] = aprioriFrequencyBands( frequencyBand )

% Frequency band names
delta     = { 'Delta' 'delta' 'd' };
theta     = { 'Theta' 'theta' 't' };
alpha     = { 'Alpha' 'alpha' 'a' };
beta      = { 'Beta'  'beta'  'b' };
gamma     = { 'Gamma' 'gamma' 'g' };

% Organisation for Human Brain Mapping definitions (Pernet et al., 2018)
deltaBand = [ 1       3.999  ];
thetaBand = [ 4       7.999  ];
alphaBand = [ 8       12.999 ];
betaBand  = [ 13      30     ];
gammaBand = [ 30.001  80     ];

% Gyurkovics and Levita (2021) theta peak finding window
thetaPlus = { 'ThetaPlus' 'Theta+' 'theta+' 'tp' 'tw' };
thetaWide = [ 2 10 ];

% Frequency limits
if ischar( frequencyBand )
    switch frequencyBand
        case delta
            frequencyLimits = deltaBand;
            bandName        = delta{1};
        case theta
            frequencyLimits = thetaBand;
            bandName        = theta{1};
        case alpha
            frequencyLimits = alphaBand;
            bandName        = alpha{1};
        case beta
            frequencyLimits = betaBand;
            bandName        = beta{1};
        case gamma
            frequencyLimits = gammaBand;
            bandName        = gamma{1};
        case thetaPlus
            frequencyLimits = thetaWide;
            bandName        = thetaPlus{1};
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

