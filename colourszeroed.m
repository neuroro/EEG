function cmap = colourszeroed( limits, resolution )
%
% •.° Zeroed Colour Map °.•
% _________________________________________________________________________
%
% Colour map with a fixed colour for zero and positive colours scaled
% separately to negative colours using specified or automatic colour scale
% limits and colour map resolution
%
% • Usage •
% -------------------------------------------------------------------------
% >> cmap = colourszeroed( limits, resolution );
% >> colormap( colourszeroed( limits, resolution ) );
%
% Examples:
% >> colormap colourszeroed
% >> colormap( colourszeroed( [] ) )
% >> colormap( colourszeroed( [ minimum maximum ] ) )
%
% <-- - -  -   -     -        -             -        -     -   -  - - -->
%
% •••( Function Inputs )
%
%   limits:     Colour scale limits as a vector or scalar
%                [] or 0             image data minimum and maximum
%                [ minimum maximum ] vector of fixed limits
%                +number             fixed maximum and image data minimum
%                -number             fixed minimum and image data maximum
%            
%                (optional input, default limits scaled to the data)
%
%   resolution: Colour map resolution in number of RGB triplets as the 
%                specified +scalar or 0 to match the current colour map
%                (optional input, default current colour map or 256)
%
% <-- - -  -   -     -        -             -        -     -   -  - - -->
% 
% [ Function Output ] = 
%   cmap:       Zeroed colormap in a resolution x 3 matrix of RGB values
%
% • Attribution •
% -------------------------------------------------------------------------
% COLOURSZEROED written by Rohan O. C. King, 2023
% @neuroro
%
% Copyright 2023 Rohan King
%
% COLOURSZEROED was inspired by bluewhitered (Childress, 2008).
% 
% Parts of the COLOURSZEROED code were adapted & generalised from
% bluewhitered (Childress, 2008).
% 
% Reference:
%   Nathan Childress (2008). bluewhitered [MATLAB code]. MATLAB Central
%    File Exchange. Retrieved July 11, 2023.
%    https://www.mathworks.com/matlabcentral/fileexchange/4058-bluewhitered
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


%% Inputs
%  ------------------------------------------------------------------------

% Colour scale limits
if ~exist( 'limits', 'var' )
    limits = [];
end

% Colour map resolution
if ~exist( 'resolution', 'var' ) || ~resolution
    resolution = size( get( gcf, 'colormap' ), 1 );
end
if resolution <= 0
    error( 'Colour map resolution in number of RGB triplets must be a positive scalar.' )
end


%% Colours
%  ------------------------------------------------------------------------

% Current colour map
currentMap       = colormap;
nColours         = size( currentMap, 1 );

% Positive value colours (RGB)
positiveMaximum  = currentMap(end,:);
iPositiveMain    = round( nColours * 7/8 );
positiveMain     = currentMap(iPositiveMain,:);
iPositiveMiddle  = round( nColours * 6/8 );
positiveMiddle   = currentMap(iPositiveMiddle,:);
iPositiveMinimum = round( nColours * 5/8 );
positiveMinimum  = currentMap(iPositiveMinimum,:);

% Zero colour (RGB)
iMiddleColour    = round( nColours * 4/8 );
fixedZero        = currentMap(iMiddleColour,:);

% Negative value colours (RGB)
iNegativeMinimum = round( nColours * 3/8 );
negativeMinimum  = currentMap(iNegativeMinimum,:);
iNegativeMiddle  = round( nColours * 2/8 );
negativeMiddle   = currentMap(iNegativeMiddle,:);
iNegativeMain    = round( nColours * 1/8 );
negativeMain     = currentMap(iNegativeMain,:);
negativeMaximum  = currentMap(1,:);

% Fixed colours
fixedPositives   = [ fixedZero; positiveMinimum; positiveMiddle; positiveMain; positiveMaximum ];
fixedNegatives   = [ negativeMaximum; negativeMain; negativeMiddle; negativeMinimum; fixedZero ];
nFixedPositives  = length( fixedPositives );
nFixedNegatives  = length( fixedNegatives );


%% Colour limits
%  ------------------------------------------------------------------------

% Current colour limits
currentLimits = clim;
clim auto

% Fixed - or + limit
if isscalar( limits )
    if limits < 0
        cLimits = [ limits currentLimits(2) ];
    elseif limits > 0
        cLimits = [ currentLimits(1) limits ];
    else
        cLimits = currentLimits; % 0 = automatic limits
    end

% Fixed limits
elseif isvector( limits )
        cLimits = limits;

% Automatic limits
elseif isempty( limits )
    cLimits = currentLimits;

end

% Adjust colour limits if needed
if cLimits ~= currentLimits
    clim( cLimits )
end

% Positive and negative existence
if cLimits(2) > 0
    positives = true;
else
    positives = false;
end
if cLimits(1) < 0
    negatives = true;
else
    negatives = false;
end


%% Positive and negative colour enumeration
%  ------------------------------------------------------------------------

if positives && negatives

    positiveProportion = cLimits(2) / diff( cLimits );              % (Childress, 2008)
    nPositive          = round( positiveProportion * resolution );  % (Childress, 2008)
    nNegative          = resolution - nPositive;                    % (Childress, 2008)
    nColours           = [ nNegative nPositive ];
    fixedColours       = { fixedNegatives fixedPositives };
    nFixedColours      = [ nFixedNegatives nFixedPositives ];

elseif positives && ~negatives

    nColours           = [ resolution resolution ];
    fixedColours       = { fixedPositives fixedPositives };
    nFixedColours      = [ nFixedPositives nFixedPositives ];

elseif negatives && ~positives

    nColours           = [ resolution resolution ];
    fixedColours       = { fixedNegatives fixedNegatives };
    nFixedColours      = [ nFixedNegatives nFixedNegatives ];

end

%% Create colour map by interpolation from fixed colours (Childress, 2008)
%  ------------------------------------------------------------------------

% Positive colour map 
if positives
    
    fixedColourDomain = linspace( 0, 1, nFixedColours(2) );
    colourDomain      = linspace( 0, 1, nColours(2) );
    positiveColourMap = zeros( nColours(2), 3 );
    for rgb = 1:3
        positiveColourMap(:,rgb) = ...
            min( max( interp1( fixedColourDomain, fixedColours{2}(:,rgb), colourDomain )', 0 ), 1 );
    end
    cmap = positiveColourMap;

end

% Negative colour map
if negatives

    fixedColourDomain = linspace( 0, 1, nFixedColours(1) );
    colourDomain      = linspace( 0, 1, nColours(1) );
    negativeColourMap = zeros( nColours(1), 3 );
    for rgb = 1:3
        negativeColourMap(:,rgb) = ...
            min( max( interp1( fixedColourDomain, fixedColours{1}(:,rgb), colourDomain )', 0 ), 1 );
    end
    cmap = negativeColourMap;

end

% Positive & negative colour map
if positives && negatives

    cmap = [ negativeColourMap; positiveColourMap ];

end


% _________________________________________________________________________
end



%% Use of bluewhitered code with adaptation and generalisation
% _________________________________________________________________________
% 
% Parts of the COLOURSZEROED code were adapted & generalised from
% bluewhitered (cited as Childress, 2008), which is covered in the
% bluewhitered licence (copied below) as use in source form with
% modifications.
%
% Licence
% -------------------------------------------------------------------------
%
% Copyright (c) 2009, Nathan Childress
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
% 
% * Redistributions of source code must retain the above copyright notice, this
%   list of conditions and the following disclaimer.
% 
% * Redistributions in binary form must reproduce the above copyright notice,
%   this list of conditions and the following disclaimer in the documentation
%   and/or other materials provided with the distribution
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
% DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE
% FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
% DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
% SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
% CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
% OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
% OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

