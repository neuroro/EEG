function Ct = globalFieldPowerComponents( setName, smoothing )
%
% •.° ERP Component Times In Global Field Power °.•
% _________________________________________________________________________
%
% Calculate global field power (GFP), the spatial standard deviation, of
% event-related EEG data, with and without smoothing using Savitsky-Golay
% filtering, and find local maxima and minima in GFP and smoothed GFP, to
% obtain candidate peak latencies for components of the event-related
% potential (ERP)
%
% Savitsky-Golay filtering may help resolve candidate ERP components in
% data with a lower signal-to-noise ratio
%
% • Usage •
% -------------------------------------------------------------------------
% Set the Current Folder to the location of the EEGLAB datasets (including
% sub-folders) from which to extract global field power and candidate ERP
% component latencies
%
% >> Ct = globalFieldPowerComponents( setName, smoothing );
%
% For example:
% >> globalFieldPowerComponents;
% >> Ct = globalFieldPowerComponents( 'LTP' );
% >> Ct = globalFieldPowerComponents( '', [5 85] );
%
% ~  ~ ~\/~ ~  ~   ~  ~ ~/\~ ~  ~   ~  ~ ~\~.~/~ ~  ~   ~  ~ ~/~'~\~ ~  ~
%
% •••( Function Inputs )
%
%   setName:   Part of the file name that is common to all the EEGLAB
%               datasets to be analysed that are located in the Current
%               Folder and sub-folders of the Current Folder, for example
%               'LTP' or '' for all datasets
%               (optional input, default all datasets)
%
%   smoothing: Vector of Savitsky-Golay filter polynomial order (positive
%               integer) and window length (odd integer) as [order length],
%               which control the degree of smoothing: smaller and longer
%               are smoother, larger and shorter are a closer fit
%               (optional input, default [5 85])
%
% [ Function Outputs ] = 
%
%   GlobalFieldPowerComponents.mat file containing GFP and smoothed GFP,
%    local maxima and minima in each, and estimated ERP component latencies
%    for all datasets
%
%   Plots of GFP and smoothed GFP with candidate ERP components labelled
%    for all datasets as <dataset>GlobalFieldPower.fig and <>SmoothG*r.fig
%
% ~  ~ ~\/~ ~  ~   ~  ~ ~/\~ ~  ~   ~  ~ ~\~.~/~ ~  ~   ~  ~ ~/~'~\~ ~  ~
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


%% Startup
% -------------------------------------------------------------------------

% Initialise EEGLAB
eegStartup

% Introduction
disp( ' ' )
disp( '•.° ERP Component Times In Global Field Power °.•' )
disp( '_________________________________________________________________________' )
disp( ' ' )
disp( 'Calculate global field power (GFP), the spatial standard deviation, of'  )
disp( 'event-related EEG data, with and without smoothing using Savitsky-Golay' )
disp( 'filtering, and find local maxima and minima in GFP and smoothed GFP, to' )
disp( 'obtain candidate peak latencies for components of the event-related'     )
disp( 'potential (ERP)'                                                         )
disp( ' ' )
disp( 'Savitsky-Golay filtering may help resolve candidate ERP components in'   )
disp( 'data with a lower signal-to-noise ratio'                                 )
disp( ' ' )


%% Default inputs
% -------------------------------------------------------------------------

% Dataset naming
if ~exist( 'setName', 'var' )
    setName   = '';
end

% Savitsky-Golay polynomial order and window length
%   5 gives peaks that are sufficiently sharp
%   85 is a low-frequency state estimate based on 100 ms delta brain states
%   and 20 ms gamma brain states (Sadaghiani, 2023) sampled at 1000 Hz
if ~exist( 'smoothing', 'var' )
    smoothing = [5 85];
end

% • Reference •
% -------------------------------------------------------------------------
% Sadaghiani, S. (2023, July 24). Timescale-overarching principles shape
%   brain state dynamics along multiple concurrent infraslow-to-fast
%   trajectories. In Temporal Organizing Principles of the Connectome
%   [Symposium]. Organization for Human Brain Mapping 2023 Annual Meeting,
%   Montréal, Quebec, Canada. https://www.humanbrainmapping.org


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


%% Event-related or evoked potential component determination
% -------------------------------------------------------------------------

% Loop through: Pairs of EEGLAB datasets
for f = 1:nFiles

    % File
    aFile       = fileList{f};
    aFolder     = folderList{f};
    aDataset    = extractBefore( aFile, '.set' );

    % Load
    loaderargin = { 'filename', aFile,   ...
                    'filepath', aFolder, ...
                    'verbose',  'off'          };
    EEG         = pop_loadset( loaderargin{:} );
    EEG         = eeg_checkset( EEG );

    % Sanity check: Trial-level (epoched) data
    if EEG.trials == 1 || ismatrix( EEG.data )
        disp( [ '[' 8 'Warning: Skipping ' aDataset ' as it has not been split into trials (not epoched)' ']' 8 ] )
        continue
    end

    % Double precision (if needed)
    if ~strcmpi( class( EEG.data ), 'double' )
        EEG.data = double( EEG.data );
    end

    % Description
    Ct.Contents  = [ 'Global field power (GFP) '             ...
                     'and estimated Components of the ERP '  ...
                     'using GFP local maxima '               ...
                     '(minima are also found)'               ];
    Ct.Contentss = [ 'Smoothed global field power (GFPs) '   ...
                     'using Savitsky-Golay filtering, '      ...
                     'and estimated Componentss of the ERP ' ...
                     'using GFPs local maxima '              ...
                     '(minima are also found)'               ];

    % Sampled times
    Ct.Times{f} = EEG.times;
    Ct.TimeUnit = 'milliseconds';


    %% Global field power and Savitsky-Golay smoothing
    %  --------------------------------------------------------------------

    % Generate ERP at every channel by averaging across trials
    erps        = mean( EEG.data, 3, 'omitnan' );
    erps        = squeeze( erps );

    % Global field power
    gfp         = std( erps, 0, 1, 'omitnan' );
    Ct.GFP{f}   = gfp;

    % Smoothed global field power
    order       = smoothing(1); % Higher is closer fit so less smoothing
    frames      = smoothing(2); % Must be odd; shorter is closer fit so less smoothing, longer is wider fit so smoother
    sggfp       = sgolayfilt( gfp, order, frames );
    Ct.GFPs{f}  = sggfp;

    Ct.VoltageUnit = 'microVolts';


    %% Find candidate ERP component times as local maxima and minima
    %  --------------------------------------------------------------------

    % Local maxima
    Ct.Maxima.iMaximaGFP{f}  = islocalmax( gfp, 2 );
    Ct.Maxima.iMaximaGFPs{f} = islocalmax( sggfp, 2 );
    Ct.Maxima.iUnit          = 'indices of times';
    Ct.Maxima.MaximaGFP{f}   = EEG.times(Ct.Maxima.iMaximaGFP{f});
    Ct.Maxima.MaximaGFPs{f}  = EEG.times(Ct.Maxima.iMaximaGFPs{f});
    Ct.Maxima.Unit           = 'milliseconds';

    % Local minima
    Ct.Minima.iMinimaGFP{f}  = islocalmax( -gfp, 2 );
    Ct.Minima.iMinimaGFPs{f} = islocalmax( -sggfp, 2 );
    Ct.Minima.iUnit          = 'indices of times';
    Ct.Minima.MinimaGFP{f}   = EEG.times(Ct.Minima.iMinimaGFP{f});
    Ct.Minima.MinimaGFPs{f}  = EEG.times(Ct.Minima.iMinimaGFPs{f});
    Ct.Minima.Unit           = 'milliseconds';

    % Candidate ERP components
    iCandidateMaxima  = and( Ct.Maxima.MaximaGFP{f} > 0, Ct.Maxima.MaximaGFP{f} < 700 );
    Ct.Components{f}  = Ct.Maxima.MaximaGFP{f}(iCandidateMaxima);
    iCandidateMaxima  = and( Ct.Maxima.MaximaGFPs{f} > 0, Ct.Maxima.MaximaGFPs{f} < 700 );
    Ct.Componentss{f} = Ct.Maxima.MaximaGFPs{f}(iCandidateMaxima);
    Ct.ComponentUnit  = 'milliseconds';


    %% Plotting
    %  --------------------------------------------------------------------

    figure

    % Plot GFP
    p = plot( EEG.times, gfp );
    title( [ aDataset ' Global Field Power' ] )
    xlabel( 'Time (ms)' )
    ylabel( 'Voltage Spatial Standard Deviation (\muV)' ) % Standard deviation across electrodes of the potential difference versus the reference channel

    % Data tips
    for dt = 1:length( Ct.Components{f} )
        x = Ct.Components{f}(dt);
        t = EEG.times == x;
        y = gfp(t);
        datatip( p, x, y );
    end

    % Save figure
    savefig( [ aDataset 'GlobalFieldPower' ] )

    clf

    % Plot smoothed GFP
    p = plot( EEG.times, sggfp );
    title( [ aDataset ' Smoothed Global Field Power' ] )
    xlabel( 'Time (ms)' )
    ylabel( 'Voltage Spatial Standard Deviation (\muV)' ) % Standard deviation across electrodes of the potential difference versus the reference channel

    % Data tips
    for dt = 1:length( Ct.Componentss{f} )
        x = Ct.Componentss{f}(dt);
        t = EEG.times == x;
        y = sggfp(t);
        datatip( p, x, y );
    end

    % Save figure
    savefig( [ aDataset 'SmoothGlobalFieldPower' ] )
    

end % for Dataset files


% Save data
save GlobalFieldPowerComponents Ct


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


