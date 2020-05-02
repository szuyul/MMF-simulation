%complex_imagesc plots the complex image on the current figure
%
% input: 
% img is an 2D image either flattened or not
% amp_range is the range of amplitude for colormap
%
% output: in the current figure, plot the complex field without ticks
%
%
% 2017-2019 Szu-Yu Lee
% Bouma Lab - The Wellman Center for Photomedicine

function [  ] = complex_imagesc(img, amp_range, varargin)

if ismember(1, size(img))
    img_dim = sqrt(length( img ));
    oimg = reshape( img , [img_dim img_dim]);
else 
    oimg = img;
end

imgmax = max(max(abs(oimg)));
field = oimg/imgmax;
if nargin == 2    % if the amplitude range is given
    image(HueOverLum(angle(field), abs(field), colormap(gca, cmap('C6')), [-pi, pi], amp_range)); 
else              % otherwise normalize the amplitude
    image(HueOverLum(angle(field), abs(field), colormap(gca, cmap('C6')), [-pi, pi], [0 1])); 
end

set(gca,'XTick',[], 'YTick',[]);    axis image

end

