function [ map ] = HueOverLum(hueMap, lumMap, colorMap, hueLim, lumLim, nanColor,...
    normalize)
  %HueOverLum(hueMap, lumMap, colorMap, hueLim, lumLim, nanColor, normalize) Create an RGB matrix representing the combined luminosity map and the hue map, with the hue given by the supplied colormap
  % Inputs:
  %        hueMap: matrix
  %        lumMap: matrix
  %        colorMap: 3 column matrix
  %        hueLim: vector of lower/upper limits for hueMap, empty for automatic
  %        lumLim: vector of lower/upper limits for lumMap, empty for automatic
  %        nanColor: nonempty 3-element vector to make nans this specific RGB color
  %        normalize: if true, keep singles/doubles from 0 to 1, otherwise output is uint8
  % Outputs:
  %        map: what we want.
  
  % This Functions follow the coding style that can be
  % sumarized in:
  % * Variables have lower camel case
  % * Functions upper camel case
  % * Constants all upper case
  % * Spaces around operators
  %
  % Authors:  Néstor Uribe-Patarroyo
  %
  % NUP: 
  % 1. Wellman Center for Photomedicine, Harvard Medical School, Massachusetts
  % General Hospital, 40 Blossom Street, Boston, MA, USA;
  % <uribepatarroyo.nestor@mgh.harvard.edu>

  % MGH Flow Measurement project
  %
  % Changelog:
  %
  % V2.0 (2017-06-01): Fixed bug (ind2rgb does not accept more than 256 colors) and added memmory management for large matrices
  %
  % Copyright Néstor Uribe-Patarroyo (2017)

  if nargin < 3 || isempty(colorMap)
    colorMap = ametrine(256);
  end
  
  if nargin < 4 || isempty(hueLim)
    hueLim(1) = nanmin(hueMap(:));
    hueLim(2) = nanmax(hueMap(:));
  end
  
  if nargin < 5 || isempty(lumLim)
    lumLim(1) = nanmin(lumMap(:));
    lumLim(2) = nanmax(lumMap(:));
  end

  if nargin < 6 || isempty(nanColor)
    nanColor = colorMap(1, :);
  end

  if nargin < 7 || isempty(normalize)
    % if normalize = true, then we don't transform into uint8 and keep
    % everything from 0 to 1
    normalize = false;
  end
  
  % Save map of nans before changing it to uint8!
  nanMap = isnan(hueMap);
  if ~isa(hueMap, 'uint8')
    % If hueMap is an 8-bit integer, don't do any normalization!
    hueMap = (hueMap - hueLim(1)) / range(hueLim);
    hueMap = Coerce(hueMap);
    if ~normalize
      hueMap = uint8(round(255 * hueMap));
    end
  end
  
  if ~isa(lumMap, 'uint8')
    % If lumMap is an 8-bit integer, don't do any normalization!
    lumMap = (lumMap - lumLim(1)) / range(lumLim);
    lumMap = Coerce(lumMap);
  else
    % If lumMap is an 8-bit integer, I need to do normalization converting to
    % single!
    lumMap = single(lumMap) / 255;
  end

  % ind2rgb supports max 256 colors!
  if size(colorMap, 1) > 256
    colorMap = interp1(linspace(0, 1, size(colorMap, 1)).', colorMap, linspace(0, 1, 256).');
  end
  
  if 0
    % I need HSL, not HSV
    hueMapRGB = ind2rgb(uint8(round(255 * hueMap)), colorMap);
    hueMapHSV = rgb2hsv(hueMapRGB);
    
    map = hsv2rgb(cat(3, hueMapHSV(:, :, 1), ones(size(lumMap), 'like', lumMap), lumMap));
  else
    if ~normalize
      hueMapRGB = single(ind2rgb(hueMap, colorMap));
    else
      hueMapRGB = ind2rgb(hueMap, colorMap);
    end
    % Especial handling of nans, make them white or nanColor
    % But only if nanColor is different from colorMap(1, :)!
    if nanColor ~= colorMap(1, :)
      for k = 1:3
        thisChannel = hueMapRGB(:, :, k);
        thisChannel(nanMap) = nanColor(k);
        hueMapRGB(:, :, k) = thisChannel;
      end
    end
    % Save memory
    clear hueMap
    if ~normalize
      % We cannot mix integer vectors with anything else
      map = uint8(round(255 * bsxfun(@times, hueMapRGB, lumMap)));
       %map = uint8(round(255* bsxfun(@times, lab2rgb(bsxfun(@times, rgb2lab(hueMapRGB), [2, 1, 1])), lumMap)));
    else
      map = bsxfun(@times, hueMapRGB, lumMap);
       %map = bsxfun(@times, lab2rgb(bsxfun(@times, rgb2lab(hueMapRGB), [2, 1, 1])), lumMap);
    end
  end
  
end

