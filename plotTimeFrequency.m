function plotTimeFrequency( eventRelatedSpectra, frequencies, times, varargin )
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
% >> plotTimeFrequency( eventRelatedSpectra, frequencies, times, ...
%                       'OptionName1', optionValue1, 'OptionName2', optionValue2, ... )
%
% For example:
% >> plotTimeFrequency( eventRelatedSpectra, frequencies, times )
% >> plotTimeFrequency( eventRelatedSpectra, frequencies, times, 'colourscale', [], 'cmap', 'jet' )
% >> plotTimeFrequency( squeeze(mean(Stimulus.SpectralPower,1)), Stimulus.Frequencies, Stimulus.Times )
% >> plotTimeFrequency( squeeze(Stimulus.GrandAverage.SpectralPower(1,:,:)), Stimulus.Frequencies, Stimulus.Times )
%
%  <-- - -  -   -     -        -             -        -     -   -  - - -->
%
% •••( Function Inputs )
%
%   eventRelatedSpectra: Time-frequency data to plot using colour values
%
%   frequencies:         Vector of logarithmically spaced frequencies in Hz
%
%   times:               Vector of time point values in milliseconds
%
%  <-- - -  -   -     -        -             -        -     -   -  - - -->
%
% :( 'Optional', Inputs )
%
%   'fticks'      Number of frequencies displayed on the y-axis
%                   (default 10)
%
%   'colourscale' Scale of the colour space
%                   [] or '' scale to the data
%                   [z1 z2]  scale from z1 to z2
%    or 'cs'        '0' or   scale to the absolute maximum so that zero is
%                   'abs'    a fixed colour, which is useful for comparing
%                            plots side-by-side
%                   '0+'     zero minimum for positive data, as above
%                   '0-'     zero maximum for negative data, as above
%                   '0='     zero centred for positive and negative data, 
%                            as above
%                   (default scale to the absolute maximum or scale to the 
%                    data if jetzeroed or bluewhitered exist)
%
%   'colourmap'   Colour map used to draw the data by MATLAB's colormap
%    or 'cmap'      Example pre-set colour spaces:
%                   'jet'
%                   'parula'
%                   'hsv'
%                   'turbo'
%                   Nx3 matrix of N gradients x RGB values
%                   or
%                   'jetzeroed' by Rohan King (2023) https://github.com/neuroro/EEG/blob/main/jetzeroed.m
%                   'bluewhitered' by Nathan Childress (2008) https://au.mathworks.com/matlabcentral/fileexchange/4058-bluewhitered
%                   (default jet or bluewhitered if it exists)
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


% Inputs
% -------------------------------------------------------------------------

if isempty( varargin )
    vars = [];
else
    vars = plotTFvars( varargin{:} );
end
vars.XL  = 'Time (ms)';
vars.X   = times;


% Frequencies spaced logarithmically
% -------------------------------------------------------------------------

vars.YL  = 'Frequency (Hz)';
vars.YS  = 'log';

% Default f-ticks
if ~isfield( vars, 'fticks' )
    vars.fticks = 10;
end

% Flip if needed
if frequencies(1) > frequencies(end)
    frequencies         = flip( frequencies );
	eventRelatedSpectra = flip( eventRelatedSpectra, 1 );
end

% Base 2 logarithmic limits
fmin = log2( frequencies(1)   );
fmax = log2( frequencies(end) );

% Octave space
F = round( 2.^( linspace( fmin, fmax, vars.fticks ) ), 1 );

% Any matching log base and exponent will give the same spacing except
% for floating point errors
% In 2022a logarithmic spaces for bases 2 vs. e vs. 10 are identical
% to 1e-13 and nearly identical to 1e-14

% Rounding resolution
df0 = find( ~diff( F ), 1 );
if ~isempty( df0 )
    F = round( 2.^( linspace( fmin, fmax, vars.fticks ) ), 2 );
end

% Frequencies
vars.Y = F;


% Colour scale
% -------------------------------------------------------------------------

% Defaults: jetzeroed or bluewhitered or zero +/- abs max
if isfield( vars, 'colourscale' )
    colourscale = vars.colourscale;
else
    if exist( 'jetzeroed.m', 'file' ) || exist( 'bluewhitered.m', 'file' )
        colourscale = [];
    else
        colourscale = '0';
    end
end

% Fixed zero colour
if ~isempty( colourscale ) && ischar( colourscale )

    datamax = max( eventRelatedSpectra, [], 'all', 'omitnan' );
    datamin = min( eventRelatedSpectra, [], 'all', 'omitnan' );

    switch colourscale

        % Zero positioned depending on the sign of the data
        case { 'abs' '0' }

            % All positive -> min always 0
            if datamin > 0
                colourscale = [ 0 datamax ];

            % All negative -> max always 0
            elseif datamax < 0
                colourscale = [ datamin 0 ];

            % Positive & negative -> centre always zero
            else
                absolutemax = max( abs( eventRelatedSpectra ), [], 'all', 'omitnan' );
                colourscale = [ -absolutemax absolutemax ];

            end

        % Zero minimum for positive data
        case '0+'
            colourscale = [ 0 datamax ];

        % Zero maximum for negative data
        case '0-'
            colourscale = [ datamin 0 ];

        % Zero centred for positive & negative data
        case '0='
            absolutemax = max( abs( eventRelatedSpectra ), [], 'all', 'omitnan' );
            colourscale = [ -absolutemax absolutemax ];

    end

end


% Plot
% -------------------------------------------------------------------------

% Set up pixel data
C = { 'XData' vars.X 'YData' vars.Y 'CData' eventRelatedSpectra };

% Draw pixel image
figure
if isempty( colourscale )
    imagesc( C{:} )
else
    imagesc( C{:}, colourscale )
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

% Colour map
if isfield( vars, 'colourmap' )
    colormap( vars.colourmap );
else
    try
        colormap( jetzeroed() )
    catch
        try
            colormap( bluewhitered() )
        catch
            colormap jet
        end
    end
end

% White background
set( gcf, 'color', 'w' );


% Axes
% -------------------------------------------------------------------------

% Logarithmic scale
set( gca, 'YScale', vars.YS )

% Values
xticklabels( 'auto' )
yticklabels( 'auto' )
ax       = gca;
ax.XTick = vars.X;
ax.YTick = vars.Y;
xticks( 'auto' )

% Limits
xlim( [ vars.X(1) vars.X(end) ] )
ylim( [ vars.Y(1) vars.Y(end) ] )

% Labels
xlabel( vars.XL )
ylabel( vars.YL )


% _________________________________________________________________________
end


%% Variable input arguments -> variables struct
% -------------------------------------------------------------------------
function vars = plotTFvars( varargin )

[ varoptions, varoops ] = plotTFcases;

for i = 1:2:length( varargin )-1

    varname  = varargin{i};
    varvalue = varargin{i+1};

    switch varname

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

        case varoptions.ft
            if isscalar( varvalue ) && isnumeric( varvalue )
                vars.fticks = varvalue;
            else
                error( varoops.ft )
            end

    end

end

end


%% Variable input cases
% -------------------------------------------------------------------------
function [ varoptions, varoops ] = plotTFcases

% Named variable input cases
varoptions.cs = {'colourscale' 'colorscale' 'colour scale' 'color scale' ...
                 'ColourScale' 'ColorScale' 'Colour Scale' 'Color Scale' ...
                 'cs' 'cscale' 'CS' 'CScale'};
varoptions.cb = {'cbar' 'colourbar' 'colorbar' 'colour bar' 'color bar' ...
                 'CBar' 'ColourBar' 'ColorBar' 'Colour Bar' 'Color Bar' ...
                 'cb' 'CB'};
varoptions.cm = {'cm' 'cmap' 'colour' 'color' 'colours' 'colors' ...
                 'colourmap' 'colormap' 'colour map' 'color map' ...
                 'CM' 'CMap' 'Colour' 'Color' 'Colours' 'Colors' ...
                 'ColourMap' 'ColorMap' 'Colour Map' 'Color Map'};
varoptions.ft = {'ft' 'ftick' 'fticks' 'tick' 'ticks' 'f tick' 'f ticks' ...
                 'FT' 'FTick' 'FTicks' 'Tick' 'Ticks' 'F tick' 'F ticks' ...
                 'F Tick' 'F Ticks'};

% Error messages for each case
varoops.cs = "Specify colour scaling according to the position of " + ...
             "zero as '0+' '0-' or '0=', or as a [z1 z2] range";
varoops.cb = "Specify the location to draw the colour bar as 'L', " + ...
             "'R', 'T', or 'B'";
varoops.cm = ['Specify the colour map as a named pre-set or as ', ...
              'a matrix of gradients x RGB'];
varoops.ft = 'Specify the number of frequencies to display on the y-axis';

end

