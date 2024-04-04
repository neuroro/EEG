function fig = plotTimeFrequency( eventRelatedSpectra, frequencies, times, varargin )
%
% •.° Plot Time-Frequency Decompositions or Event-Related Spectra °.•
% _________________________________________________________________________
%
% Plot event-related spectra as a heatmap with colours scaled to an input
% range or with a fixed colour for zero positioned according to the sign of
% the data, which can facilitate visual comparison of plots
%
% • Usage •
% -------------------------------------------------------------------------
% >> fig = plotTimeFrequency( eventRelatedSpectra, frequencies, times, ...
%            'OptionName1', optionValue1, 'OptionName2', optionValue2, ... )
%
% For example:
% >> fig = plotTimeFrequency( eventRelatedSpectra, frequencies, times );
% >> plotTimeFrequency( eventRelatedSpectra, frequencies, times, 'colourscale', [z1 z2], 'cmap', 'jet' );
% >> plotTimeFrequency( Stimulus.GrandAverage.SpectralPower, Stimulus.Frequencies, Stimulus.Times );
%
%  <-- - -  -   -     -        -             -        -     -   -  - - -->
%
% •••( Function Inputs )
%
%   eventRelatedSpectra: Time-frequency data to plot using colour values
%
%   frequencies:         Vector of frequencies in Hz
%
%   times:               Vector of time point values in milliseconds
%
%  <-- - -  -   -     -        -             -        -     -   -  - - -->
%
% •••( 'Optional', Inputs )
%
%   'colourvalue' Colour values to plot
%    or 'cv'        'exact'   exact time-frequency data values
%                   'ratio'   power ratio relative to baseline
%                             (for time-frequency data values of power in
%                              decibels relative to baseline)
%                   'percent' percent change in power relative to baseline
%                             (for time-frequency data values of power in
%                              decibels relative to baseline)
%                   (default exact)
%
%   'colourscale' Scale of the colour space
%    or 'cs'        [] or ''  scale to the data
%                   [z1 z2]   scale from z1 to z2
%                   '0' or    scale to the absolute maximum so that zero is
%                   'abs'     a fixed colour, which is useful for comparing
%                             plots side-by-side
%                   '0+'      zero minimum for positive data, as above
%                   '0-'      zero maximum for negative data, as above
%                   '0='      zero centred for positive and negative data, 
%                             as above
%                   (default scale to the data if jetzeroed, colourszeroed,
%                    or bluewhitered exist, otherwise scale to the absolute
%                    maximum)
%
%   'colourmap'   Colour map used to draw the data
%    or 'cmap'      Example pre-set colour spaces:
%                   'jet'
%                   'parula'
%                   'hsv'
%                   'turbo'
%                   Nx3 matrix of N gradients x RGB values
%                   or
%                   'jetzeroed' by Rohan King (2023) https://github.com/neuroro/EEG/blob/main/jetzeroed.m
%                   'colourszeroed' by Rohan King (2023) https://github.com/neuroro/EEG/blob/main/colourszeroed.m
%                   'bluewhitered' by Nathan Childress (2008) https://au.mathworks.com/matlabcentral/fileexchange/4058-bluewhitered
%                   (default jetzeroed, colourszeroed, or bluewhitered if 
%                    they exist, otherwise jet)
%
%   'colourbar'   Draw a colour bar at the specified location
%    or 'cbar'      'T'   top outside the plot box
%                   'B'   bottom outside the plot box
%                   'R'   right outside the plot box
%                   'L'   left outside the plot box
%                   'off' none
%                   (default on the right)
%
% • Author •
% -------------------------------------------------------------------------
% Rohan O. C. King, 2024
% @neuroro
%
% Copyright 2024 Rohan King
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


% Optional inputs
% -------------------------------------------------------------------------

if isempty( varargin )
    vars = [];
else
    vars = plotTFvars( varargin{:} );
end


% Relative power options
% -------------------------------------------------------------------------

ratioColour   = false;
percentColour = false;
if isfield( vars, 'colourvalue' )
    switch vars.colourvalue

        % Power ratio relative to baseline
        case 'ratio'
            eventRelatedSpectra = 10 .^ ( eventRelatedSpectra / 10 ) - 1;   % Scale colours relative to 1 by subtracting 1 from the ratio then adding 1 to the colour bar labels
            ratioColour         = true;

        % Percent change in power relative to baseline
        case 'percent'
            eventRelatedSpectra = ( 10 .^ ( eventRelatedSpectra / 10 ) - 1 ) * 100;
            percentColour       = true;

    end
end


% Colour scale
% -------------------------------------------------------------------------

% Defaults: jetzeroed, colourszeroed, bluewhitered, or zero +/- abs max
if isfield( vars, 'colourscale' )
    colourScale = vars.colourscale;
else
    if exist( 'jetzeroed.m',     'file' ) || ...
       exist( 'colourszeroed.m', 'file' ) || ...
       exist( 'bluewhitered.m',  'file' )
        colourScale = [];
    else
        colourScale = '0';
    end
end

% Fixed zero colour
if ~isempty( colourScale ) && ischar( colourScale )

    dataMax = max( eventRelatedSpectra, [], 'all', 'omitnan' );
    dataMin = min( eventRelatedSpectra, [], 'all', 'omitnan' );

    switch colourScale

        % Zero positioned depending on the sign of the data
        case { 'abs' '0' }

            % All positive -> min always 0
            if dataMin > 0
                colourScale = [ 0 dataMax ];

            % All negative -> max always 0
            elseif dataMax < 0
                colourScale = [ dataMin 0 ];

            % Positive & negative -> centre always zero
            else
                absoluteMax = max( abs( eventRelatedSpectra ), [], 'all', 'omitnan' );
                colourScale = [ -absoluteMax absoluteMax ];

            end

        % Zero minimum for positive data
        case '0+'
            colourScale = [ 0 dataMax ];

        % Zero maximum for negative data
        case '0-'
            colourScale = [ dataMin 0 ];

        % Zero centred for positive & negative data
        case '0='
            absoluteMax = max( abs( eventRelatedSpectra ), [], 'all', 'omitnan' );
            colourScale = [ -absoluteMax absoluteMax ];

    end

end


% Plot
% -------------------------------------------------------------------------

% Draw pixel image
fig = figure;
if isempty( colourScale )
    imagesc( times, frequencies, eventRelatedSpectra )
else
    imagesc( times, frequencies, eventRelatedSpectra, colourScale )
end

% Colour bar
if isfield( vars, 'colourbar' )
    switch vars.colourbar
        case 'T'
            location = 'NorthOutside';
        case 'B'
            location = 'SouthOutside';
        case 'L'
            location = 'WestOutside';
        case 'R'
            location = 'EastOutside';
    end
else
    location = 'EastOutside';
end
if exist('location','var')
    colorbar( gca, location )
end

% Colour bar labels for power ratio relative to baseline
if ratioColour
    cb = colorbar;
    for c = 1:length( cb.TickLabels )
        cb.TickLabels{c} = string( str2double( cb.TickLabels{c} ) + 1 );    % Add 1 to correct for earlier subtraction
    end
end

% Colour bar labels for percent change in power relative to baseline
if percentColour
    cb = colorbar;
    for c = 1:length( cb.TickLabels )
        cb.TickLabels{c} = [ cb.TickLabels{c} '%' ];
    end
end

% Colour map
if isfield( vars, 'colourmap' )
    colormap( vars.colourmap );
else
    try
        colormap( jetzeroed() )
    catch
        try
            colormap( colourszeroed() )
        catch
            try
                colormap( bluewhitered() )
            catch
                colormap jet
            end
        end
    end
end

% White background
set( gcf, 'color', 'w' );


% Axes
% -------------------------------------------------------------------------

xticklabels( 'auto' )
yticklabels( 'auto' )

ax = gca;

% Frequency spacing
if std( diff( frequencies ) ) < 1e-12
    frequencySpacing = 'linear';
else
    frequencySpacing = 'log';
end
ax.YScale = frequencySpacing;

% Frequencies increasing from the origin
ax.YDir = 'normal';

% Logarithmically spaced frequency labels
if strcmp( frequencySpacing, 'log' )
    fMin            = frequencies(1);
    fMax            = frequencies(end);
    nOctaves        = log2( fMax / fMin );
    nTicks          = floor( 3 * nOctaves );
    fTicks          = 2 .^ ( linspace( log2( fMin ), log2( fMax ), nTicks ) );
    fTicks(2:end-1) = round( fTicks(2:end-1), 0 );  % Round the middle ticks
    fTicks          = unique( fTicks );             % Keep unique ticks after rounding
    ax.YTick        = fTicks;                       % Set frequency ticks
    ylim( [ fMin fMax ] )                           % Remove white space
end

% Labels
xlabel( 'Time (ms)'      )
ylabel( 'Frequency (Hz)' )


% _________________________________________________________________________
end



%% 
% •.° Variable Input Arguments -> Variables Struct °.•
% _________________________________________________________________________
%
function vars = plotTFvars( varargin )

[ varoptions, varoops ] = plotTFcases;

for i = 1:2:length( varargin )-1

    varname  = varargin{i};
    varvalue = varargin{i+1};

    switch varname

        case varoptions.cv
            if ischar( varvalue )
                vars.colourvalue = varvalue;
            else
                error( varoops.cv )
            end

        case varoptions.cs
            if isvector( varvalue ) || ischar( varvalue ) || isempty( varvalue )
                vars.colourscale = varvalue;
            else
                error( varoops.cs )
            end

        case varoptions.cb
            if ischar( varvalue )
                vars.colourbar = varvalue;
            else
                error( varoops.cb )
            end

        case varoptions.cm
            if ischar( varvalue ) || isstring( varvalue ) ...
               || ( ismatrix( varvalue ) && size( varvalue, 2 ) == 3 )
                vars.colourmap = varvalue;
            else
                error( varoops.cm )
            end

    end

end


% _________________________________________________________________________
end



%% 
% •.° Variable Input Cases °.•
% _________________________________________________________________________
%
function [ varoptions, varoops ] = plotTFcases

% Named variable input cases
varoptions.cv = { 'colourvalue' 'colorvalue' 'colour value' 'color value' ...
                  'ColourValue' 'ColorValue' 'Colour Value' 'Color Value' ...
                  'cv' 'cvalue' 'CV' 'CValue'                             };
varoptions.cs = { 'colourscale' 'colorscale' 'colour scale' 'color scale' ...
                  'ColourScale' 'ColorScale' 'Colour Scale' 'Color Scale' ...
                  'cs' 'cscale' 'CS' 'CScale'                             };
varoptions.cb = { 'cbar' 'colourbar' 'colorbar' 'colour bar' 'color bar' ...
                  'CBar' 'ColourBar' 'ColorBar' 'Colour Bar' 'Color Bar' ...
                  'cb' 'CB'                                              };
varoptions.cm = { 'cm' 'cmap' 'colour' 'color' 'colours' 'colors' ...
                  'colourmap' 'colormap' 'colour map' 'color map' ...
                  'CM' 'CMap' 'Colour' 'Color' 'Colours' 'Colors' ...
                  'ColourMap' 'ColorMap' 'Colour Map' 'Color Map' };

% Error messages for each case
varoops.cv = "Specify colour values as 'exact', 'ratio', or 'percent'";
varoops.cs = "Specify colour scaling according to the position of " + ...
             "zero as '0+' '0-' or '0=', or as a [z1 z2] range";
varoops.cb = "Specify the location to draw the colour bar as 'L', " + ...
             "'R', 'T', or 'B'";
varoops.cm = "Specify the colour map as a named pre-set or as " + ...
             "a matrix of gradients x RGB";


% _________________________________________________________________________
end


