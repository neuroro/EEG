function convertFigures( fileFormat, resolution, figurePath, outputPath )
%
% •.° Convert Figures °.•
% _________________________________________________________________________
%
% Convert all matlab figure files in a folder to an image file format
%
% Vector images with transparent backgrounds can be saved as pdf files
% Flat images at the specified resolution can be saved as jpg, tiff, or png
%
% • Usage •
% -------------------------------------------------------------------------
% >> convertFigures( fileFormat, resolution, figurePath, outputPath )
%
% For example:
% >> convertFigures
% >> convertFigures( 'pdf' )
% >> convertFigures( 'jpg', 600, 'C:\Users\<user>\Documents\MATLAB' )
%
% .........................................................................
%
% •••( Function Inputs )
%
%   fileFormat: 'pdf', 'jpg', 'tiff', or 'png' 
%                (optional input, default pdf)
%
%   resolution: image resolution in dots per inch
%                (optional input, default 600, redundant for pdf)
%
%   figurePath: folder containing matlab figures to convert 
%                (optional input, default Current Folder)
%
%   outputPath: folder in which to save figure images
%                (optional input, default same folder)
%
% .........................................................................
%
% [ Function Outputs ] =
%
%   Figure images saved in the specified file format
%
% • Authors •
% -------------------------------------------------------------------------
% Rohan O. C. King & Chris C. King, 2023
% @neuroro
%
% Copyright 2023 Rohan King & Chris King
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

% Format
if ~exist( 'fileFormat', 'var' )
    fileFormat = 'pdf';
end
if ~matches( fileFormat(1), '.' )
    fileFormat = [ '.' fileFormat ];
end

% Image resolution
if ~exist( 'resolution', 'var' )
    resolution = 600;
end

% Figure folder path
if ~exist( 'figurePath', 'var' )
    figurePath = cd;
end

% Output folder path
if ~exist( 'outputPath', 'var' )
    outputPath = figurePath;
end


% Convert figures
% -------------------------------------------------------------------------

% Listing of all figures in the chosen folder
wildFigurePath = fullfile( figurePath, '*.fig' );
figureFiles    = dir( wildFigurePath );
nFigures       = length( figureFiles );

% Load and export
for f = 1:nFigures

    % Figure name and path
    fileName     = figureFiles(f).name;
    figureName   = extractBefore( fileName, '.fig' );
    fileFullPath = fullfile( figurePath, fileName );

    % Load figure object
    fig = openfig( fileFullPath );

    % Output image name and path
    imageFileName = [ figureName fileFormat ];
    imageFilePath = fullfile( outputPath, imageFileName );

    % Image options
    if matches( fileFormat, '.pdf', 'IgnoreCase', true )
        imargin = { 'ContentType', 'vector', 'BackgroundColor', 'none' };
    else
        imargin = { 'ContentType', 'image', 'Resolution', resolution };
    end

    % Export image
    exportgraphics( fig, imageFilePath, imargin{:} )

    % Clean up
    clf
    close( fig )

end


% _________________________________________________________________________
end

