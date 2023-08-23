function eegMultiScaleEntropyData( timeScales, setName )
%
% •.° Composite Multi-Scale Entropy °.•
% _________________________________________________________________________
%
% Measure the composite multi-scale entropy (CMSE) of electrical brain
% activity in all EEGLAB datasets named in common that are located in the
% Current Folder and sub-folders of the Current Folder and save
% <Dataset>MSE.mat files for each dataset at the single-trial level
% 
% MSE needs a large number of continuous time points to accurately
% determine entropy, particularly at higher time-scales:
%     10000 points is the default recommendation
%   > 30000 points is recommended for accuracy at higher time-scales
%      4000 points is advised as an absolute minimum
%
% • Usage •
% -------------------------------------------------------------------------
% Function:
% >> eegMultiScaleEntropyData( timeScales, setName )
%
% Examples:
% >> eegMultiScaleEntropyData
% >> eegMultiScaleEntropyData( 20, '' )
% >> eegMultiScaleEntropyData( [1 20], 'Rest' )
% 
% Inputs:
%   timeScales: Estimate entropy for scales 1 to the specified scale or
%               specify a vector of scale limits or of individual scales
%               (optional input, default scales 1 - 20)
%   setName:    Name in common to the datasets, for example 'Memory' or
%               'Rest', or '' for all datasets in the Current Folder and
%               sub-folders of the Current Folder
%               (optional input, default all datasets)
%
% Outputs:
%   <Dataset>MSE.mat files for each dataset containing individual CMSE at
%   the single-trial level
%
% • Author •
% -------------------------------------------------------------------------
% Rohan O. C. King, 2023
% @neuroro
%
% Copyright 2023 Rohan King
%
% Composite multi-scale entropy (CMSE) was elucidated and defined by Wu et 
% al. (2013), who demonstrated it to have superior convergence.
%
% The compositeMultiScaleEntropy, coarseGrain, and sampleEntropy
% sub-functions that calculate CMSE in this code were written with
% reference to Wu et al. (2013), Lee (2012), and Lu and Wang (2021).
%
% Estimated reliability of CMSE as a function of signal length is in
% reference to Wu et al. (2013). Estimated accuracy of sample entropy
% values as a function of signal length is inferred from Richman and
% Moorman (2000).
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
% • References •
% -------------------------------------------------------------------------
% Lee, K. (2012). Sample Entropy. MATLAB Central File Exchange.
%   https://www.mathworks.com/matlabcentral/fileexchange/35784-sample-entropy
% Lu, J., & Wang, Z. (2021). The Systematic Bias of entropy calculation in
%   the Multi-scale entropy algorithm. Entropy, 23(6), 659.
% Richman, J. S., & Moorman, J. R. (2000). Physiological time-series
%   analysis using approximate entropy and sample entropy. American Journal
%   of Physiology-Heart and Circulatory Physiology, 278(6), H2039-H2049.
%   https://doi.org/10.1152/ajpheart.2000.278.6.H2039
% Wu, S.-D., Wu, C.-W., Lin, S.-G., Wang, C.-C., & Lee, K.-Y. (2013). Time
%   series analysis using composite multiscale entropy. Entropy, 15(3),
%   1069-1084.


% Introduction
disp( ' ' )
disp( '•.° Composite Multi-Scale Entropy °.•' )
disp( '_________________________________________________________________________' )
disp( ' ' )
disp( 'Measure the composite multi-scale entropy (CMSE) of electrical brain'      )
disp( 'activity in all EEGLAB datasets named in common that are located in the'   )
disp( 'Current Folder and sub-folders of the Current Folder and save'             )
disp( '<Dataset>MSE.mat files for each dataset at the single-trial level'         )
disp( ' ' )


%% Input handling
% -------------------------------------------------------------------------

% Default scales
if ~nargin
    timeScales = 20;
end

% Scalar input
if isscalar( timeScales )
    timeScales = [1 timeScales];
end

% Wildcard dataset name
aWildError = [ 'Input a name in common to all EEGLAB datasets ' ...
               'located in the Current Folder and sub-folders ' ...
               'to measure the composite multi-scale entropy '  ...
               'that is present within them'                    ];
if nargin < 2
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
    aWildDataset = '*.set';
else
    aWildDataset = [ '*' setName '*.set' ];
end


%% EEGLAB dataset files
% -------------------------------------------------------------------------

% Search for all datasets named in common located in the current folder and
% sub-folders of the current folder
FilesStruct = dir( [ '**/' aWildDataset ] );
nFiles      = length( FilesStruct );
fileList    = { FilesStruct(:).name   };
folderList  = { FilesStruct(:).folder };


%% Parameters
% -------------------------------------------------------------------------

% Load first EEG dataset
disp( 'Determining parameters from the first dataset...' )
disp( ' ' )
EEG                = pop_loadset( 'filename', fileList{1},   ...
                                  'filepath', folderList{1}, ...
                                  'verbose', 'off'           );
EEG                = eeg_checkset( EEG );

% Determine parameters
samplingRate       = EEG.srate;
nChannels          = EEG.nbchan;
nTrials            = EEG.trials;
nPoints            = EEG.pnts;
clear EEG

% Time-scales
timeScaleLimits{1} = [ num2str( timeScales(1)   / samplingRate * 1000 ) ' ms' ];
timeScaleLimits{2} = [ num2str( timeScales(end) / samplingRate * 1000 ) ' ms' ];


%% Stability - do you have enough samples?
% -------------------------------------------------------------------------

% 30000+ points
if nPoints >= 30000
    disp( 'You have enough time samples to accurately estimate MSE at all time-scales' )
    disp( [ 'from 1 to ' num2str( round( nPoints / 300 ) ) ' with good approximation thereafter'] )

% 10000 - 30000 points
elseif ( nPoints > 10000 && nPoints < 30000 )
    disp( 'You have enough time samples to estimate MSE with good accuracy up to' )
    disp( [ 'high time-scales. You have ' num2str( ( nPoints / 10000 - 1 ) * 100 ) '% more time samples than 10000,' ] )
    disp( 'which has been shown to give reliable estimates (Wu et al., 2013)' )

% 10000 points
elseif nPoints == 10000
    disp( 'You have enough time samples to estimate MSE with good accuracy up to' )
    disp( 'high time-scales. You have 10000 time samples, which has been shown to' )
    disp( 'give reliable estimates (Wu et al., 2013)' )

% 4000 - 10000 points
elseif ( nPoints >= 4000 && nPoints < 10000 )
    if nPoints >= 7000
        alert   = 'Caution';
        qualify = 'decent';
        caveat  = 'is decent.';
    elseif ( nPoints >= 5000 && nPoints < 7000 )
        alert   = 'Warning';
        qualify = 'acceptable';
        caveat  = 'is acceptable.';
    else
        alert   = 'WARNING';
        qualify = 'limited';
        caveat  = 'is just sufficient.';
    end
    warningText = [ '[' 8 alert ': You have enough time samples to estimate CMSE with ' qualify '\n'   ...
                          'accuracy, but estimates are likely to be unstable at higher time-scales.\n' ...
                          'Averaging across trials is necessary (done for you in the last step).\n'    ...
                          'Averaging across a scientifically-valid electrode cluster may help.\n'      ']' 8 ];
    warningTxt2 = [ '[' 8 'You have ' num2str( ( nPoints / 4000 - 1 ) * 100 ) '%' ...
                          ' more points than the minimum of 4000, which ' caveat  ']' 8 ];
    fprintf( warningText )
    disp( warningTxt2 )

% 4000 points
elseif nPoints == 4000
    warningText = [ '[' 8 'WARNING: You have enough time samples to estimate CMSE with limited\n'   ...
                          'accuracy, but estimates may be unstable at higher time-scales.\n'        ...
                          'Averaging across trials is necessary (done for you in the last step).\n' ...
                          'Averaging across a scientifically-valid electrode cluster may help.\n'   ...
                          'You have 4000 time samples, which is the minimum shown to give\n'        ...
                          'adequate convergence (Wu et al., 2013).\n'                               ']' 8 ];
    fprintf( warningText )

% < 4000 points
elseif nPoints < 4000
    if nPoints >= 3000
        qualify = { 'TOO' 'may' };
        caveat  = [ 'means CMSE estimates are likely to be inaccurate, '            ...
                    'but time-scales up to ' num2str( round( nPoints / 200 ) ) '\n' ...
                    'might be usable with extensive averaging.'                     ];
    else
        qualify = { 'FAR TOO' 'will' };
        caveat  = 'means CMSE estimates, if they exist, are likely to be highly inaccurate.';
    end
    warningText = [ '[' 8 'WARNING: YOUR SIGNAL IS ' qualify{1} ' SHORT! CMSE estimates ' qualify{2} ' be inaccurate!\n' ...
                          'Averaging across trials is necessary (done for you in the last step).\n'                    ...
                          'Averaging across a scientifically-valid electrode cluster is probably\n'                    ']' 8 ];
    warningTxt2 = [ '[' 8 'necessary, but may not be sufficient. You have ' num2str( nPoints/4000*100 ) '% of 4000 points, which' ']' 8 ];
    fprintf( warningText )
    disp( warningTxt2 )
    fprintf( [ '[' 8 caveat '\n' ']' 8 ] )

end

disp( ' ' )


%% Predicted accuracy
% -------------------------------------------------------------------------

% Richman and Moorman (2000) show confidence intervals for sample entropy
% with m = 2 and r = 0.2 for different simulated signal lengths.
% The following is inferred from their work, for a range of numbers of 
% signal or scaled-signal points:
% < 135 - 150 entropy estimates have extreme confidence intervals and so
%             are likely to be highly inexact, but may be indicative of the
%             trend in the context of multiple scales
%   150 - 200 entropy estimates have large confidence intervals and so are 
%             likely to be inexact but may indicate the trend across scales
%   200 - 300 entropy estimates have medium confidence intervals and so
%             are likely to be usable approximations
%   300 - 500 entropy estimates have small confidence intervals and so are
%             likely to be fairly accurate
%   > 500     entropy estimates have tiny confidence intervals and so are
%             likely to be quite exact

% Signal points needed for sample entropy to be accurate to heuristically 
% different degrees or thresholds
pointThresholds = [ 150 200 300 500 ];
pointThresholds = flip( pointThresholds );

% Equivalent scale thresholds for the EEG data
scaleThresholds = round( nPoints ./ pointThresholds );
scaleThresholds = [ 1 scaleThresholds Inf ];

% Scale zones
scaleZones = cell(1,5);
for z = 1:5
    scaleZones{z} = [ num2str( scaleThresholds(z) ) ' - ' num2str( scaleThresholds(z+1) ) ];
end

% Confidence zones
confidenceZones = { [ '  Sample entropy may be exact, with a tiny\n'      ...
                      '  confidence interval of the estimate\n'           ] ...
                    [ '  Sample entropy may be fairly accurate, with\n'   ...
                      '  a small confidence interval of the estimate\n'   ] ...
                    [ '  Sample entropy may be approximate, with a \n'    ...
                      '  medium confidence interval of the estimate\n'    ] ...
                    [ '  Sample entropy may be fairly inaccurate, with\n' ...
                      '  a large confidence interval of the estimate;\n' ...
                      '  but indicative of the entropy trend across\n'    ...
                      '  scales\n'                                        ] ...
                    [ '  Sample entropy may be highly inexact, with an\n' ...
                      '  extreme confidence interval of the estimate;\n'  ...
                      '  but potentially suggestive of the entropy\n'     ...
                      '  trend across scales\n'                           ] };


%% Display CMSE information
% -------------------------------------------------------------------------

disp( '• Parameters •' )
disp( '-------------------------------------------------------------------------' )
if isempty( setName )
    disp( [ num2str( nFiles ) ' datasets' ] )
else
    disp( [ num2str( nFiles ) ' ' setName ' datasets' ] )
end
disp( [ num2str( samplingRate ) ' Hz sampling rate'] )
disp( [ num2str( nChannels ) ' channels'] )
disp( [ num2str( nTrials ) ' trials' ] )
disp( [ num2str( nPoints ) ' time samples at the sampling rate' ] )
disp( ' ' )
disp( '• Entropy in the EEG •' )
disp( '-------------------------------------------------------------------------' )
disp( [ 'Measured at scales ' num2str( timeScales(1) ) ' - ' num2str( timeScales(end) ) ] )
disp( [ 'Equivalent to time-scales of ' timeScaleLimits{1} ' - ' timeScaleLimits{2} ] )
disp( [ 'Probably realisable at time-scales of ' num2str( 1000/samplingRate ) ' ms - ' num2str( scaleThresholds(end-2)*1000/samplingRate ) ' ms' ] )
disp( ' ' )
disp( '• Predicted accuracy •')
disp( '-------------------------------------------------------------------------' )
disp( 'The following predictions are approximate' )
disp( 'They are numerically derived from simulated data but not validated' )
for z = 1:length( confidenceZones )
    disp( [ 'Scale ' scaleZones{z} ':' ] )
    fprintf( confidenceZones{z} )
end
disp( ' ' )


%% Calculate composite multi-scale entropy (CMSE) for each dataset
% -------------------------------------------------------------------------

% Start time
mseRunTime( 'Started at' )

% Loop through: EEGLAB datasets
for f = 1:nFiles

    % File
    currentFile   = fileList{f};
    currentFolder = folderList{f};
    currentName   = extractBefore( currentFile, '.set' );
    disp( [ 'Measuring the entropy in ' currentName ] )
    disp( ' ' )

    % Pre-allocate
    cmse          = NaN( nChannels, nTrials, timeScales(end) );

    % Load
    EEG           = pop_loadset( 'filename', currentFile, 'filepath', currentFolder, 'verbose', 'off' );
    EEG           = eeg_checkset( EEG );
    data          = EEG.data;
    chanlocs      = EEG.chanlocs;

    % Sanity checks
    if EEG.nbchan > nChannels
        warning( [ 'The number of channels in the current dataset is greater '     ...
                   'than in the first dataset. Excess channels are not processed!' ] )
        disp( ' ' )
    elseif EEG.nbchan < nChannels
        warning( [ 'The number of channels in the current dataset is '      ...
                   'fewer than in the first dataset. This usually matters!' ] )
        disp( ' ' )
    end

    % Loop through: Channels
    parfor ch = 1:nChannels

        % Loop through: Trials
        for t = 1:nTrials

            % EEG signal
            eegSignal = squeeze( data(ch,:,t) );

            % CMSE for each channel x trial x time-scale
            cmse(ch,t,:) = compositeMultiScaleEntropy( eegSignal, timeScales ); %#ok

        end

    end

    % Save
    if contains( currentName, ' ' )
        fileName = [ currentName ' MSE' ];
    elseif contains( currentName, '_' )
        fileName = [ currentName '_MSE' ];
    elseif contains( currentName, '-' ) && ~contains( currentName, '_' )
        fileName = [ currentName '-MSE' ];
    else
        fileName = [ currentName 'MSE' ];
    end
    fullFilePath = fullfile( currentFolder, fileName );
    save( fullFilePath, 'cmse', 'timeScales', 'timeScaleLimits', 'samplingRate', 'chanlocs' )

    % File finish time
    mseRunTime( [ currentName ' finished at' ] )

end


% _________________________________________________________________________
end



%% Composite Multi-Scale Entropy
% _________________________________________________________________________
%
% Wu, et al. (2013), Lu and Wang (2021)
%
function mse = compositeMultiScaleEntropy( signal, timeScales, distanceThreshold )

% Crete vector of scales if input is scale limits
if length( timeScales ) == 2
    timeScales = timeScales(1):timeScales(end);
end

% Default distance threshold as a percentage of the standard deviation
if nargin < 3
    distanceThreshold = 0.2; % Default set to 20% to align with Richman and Moorman (2000)
end

% Pre-allocate
mse = zeros( 1, length( timeScales ) );

% Loop through: Time-scales
for s = timeScales

    % Loop through: Coarse grained time-series
    for cg = 1:s

        % Coarse grained signal
        scaledSignal = coarseGrain( signal(cg:end), s );

        % Tolerance
        % Percentage of the standard deviation of the coarse-grained signal
        r = distanceThreshold * std( scaledSignal , 0, 'omitnan' );

        % Composite multi-scale entropy
        mse(s) = mse(s) + sampleEntropy( scaledSignal, r ) / s;

    end

end


% _________________________________________________________________________
end



%% Coarse-grained signal
% _________________________________________________________________________
%
% Wu, et al. (2013)
%
function scaledSignal = coarseGrain( signal, scale )

% Signal time points
nTimes = length( signal );

% Coarse-grained limit
cgMax  = floor( nTimes / scale );

for cg = 1:1:cgMax

    coarseGrainIndices = (cg - 1) * scale + 1:cg * scale;

    scaledSignal(cg)   = mean( signal(coarseGrainIndices) );                %#ok

end


% _________________________________________________________________________
end



%% Sample entropy
% _________________________________________________________________________
%
% For embedding dimension, m = 2, and a tolerance, r
%
function se = sampleEntropy( signal, r )

% Data points
nPoints = length( signal );

% Initialise
Nn = 0;
Nd = 0;

% Loop through: All pairs of points with m successive pairs of points
for i = 1:nPoints - 2
    for j = i + 1:nPoints - 2

        % Point pair differences
        pointDifference   = abs( signal(i)   - signal(j)   );
        m1pointDifference = abs( signal(i+1) - signal(j+1) );
        m2pointDifference = abs( signal(i+2) - signal(j+2) );

        % Test point pair and next point pair difference
        if pointDifference < r && m1pointDifference < r

            % Count of close point pairs & close next point pairs
            % (m successive pairs of close points)
            Nn = Nn + 1;

            % Test m-next point pair difference
            if m2pointDifference < r

                % Count of close m-next point pairs
                % (m+1 successive pairs of close points)
                Nd = Nd + 1;

            end

        end

    end
end

% Sample entropy
% SE = - ln( m+1 successive pairs of close points 
%            / m successive pairs of close points )
se = -log( Nd / Nn );


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


