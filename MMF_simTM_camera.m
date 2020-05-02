%MMF_simTM_camera  simulates a transmission matrix in a focal grid representation (e.g., 2D camera pixels) 
% considering both horizontal(H) and vertical(V) polarizations of a MMF given the fiber specifications and geometry.
%
%[ T_HH, T_HV, T_VH, T_VV, mode_to_x, mode_to_y, x_to_mode, y_to_mode, ill_E, propconst ]...
%            = MMF_simTM_camera( lambda, D, NA, Length, Rho, Theta, N, input_num )
% 
% outputs:
% T_HH ... T_VV (N^2 by n) are the four quadrants of full transmission matrix considering both H and V polarizations. 
%   An column vector of each quadrant is a flattend N by N 2D image of MMF spatial response to a correspondingly coupled input channel.
%   e.g., T_VH is the transmission maxtrix with H polarization inputs and V polarization outputs 
% mode_to_H, mode_to_V are the basis transformation matrices from PIM mode basis to H or V polarization output
% H_to_mode, V_to_mode are the basis transformation matrices from H or V polarization input coupled to PIM mode basis 
% ill_E is a N^2 by input_num image stack of input coupling to the MMF, each column vector is an flattened N by N 2D image of a focal spot 
%   at a specific spatial channel.
%
% inputs:
% lambda is the operating wavelength (unit: m)
% D is the fiber core diameter of the MMF (unit: m)
% NA is the numerical aperture of the MMF 
% Length is an array of the segmental length of the MMF (unit: m)
% Rho is an array of the bending radius of curvature (unit: m) of each defined segment. If the segment is straight, Rho = inf
% Theta is an array of the orientation of the bending of the MMF (unit: rad.)
% N is the image dimension
% input_num is the number of coupling focal spots to the MMF, should be a square of an integer, e.g., 400
%
%
% 2017-2019 Szu-Yu Lee
% Bouma Lab - The Wellman Center for Photomedicine

function [ T_HH, T_HV, T_VH, T_VV, mode_to_H, mode_to_V, H_to_mode, V_to_mode, ill_E, propconst ]...
            = MMF_simTM_camera( lambda, D, NA, Length, Rho, Theta, N, input_num )
%% calculate the TM in PIM representation
[ T, ~, ~, ~, ~, propconst, Er, Ep, ~, ~, ~, ~, img_size ] = MMF_simTM_PIM( lambda, D, NA, Length, Rho, Theta, N );

%% excitation of MMF modes with delta spikes LP input modes
% initialization
focus_waist = 0.6*lambda/NA;                                                % effective focal spot coupled into the fiber, depends on NA of MMF
                                                                            % using higher NA of Obj. causes loss due to limited acceptance angle of MMF
NMode = size(T,1);
input_dim = sqrt(input_num);

% define the input focal spots 
% assume the input is polarized along x(y) direction => Ex(Ey)
ill_E = zeros(N,N, input_num);                                              % illumination Ex(Ey) field
[x, y] = meshgrid(linspace(-img_size/2, img_size/2, N));
[theta, ~] = cart2pol(x,y);

% 2D XY raster scanning 
for xx = 1:input_num
    [fyi, fxi] = ind2sub([input_dim input_dim], xx);
    fxi = fxi - (input_dim + 1)/2;
    fyi = fyi - (input_dim + 1)/2;
    ill_E(:,:,xx) = exp( -( (x - fxi*(img_size/input_dim)).^2 + (y - fyi*(img_size/input_dim)).^2 ) / focus_waist^2 );
    ill_E_amp = sqrt(sum(sum(ill_E(:,:,xx).*conj(ill_E(:,:,xx)))));
    ill_E(:,:,xx) = ill_E(:,:,xx) / ill_E_amp;
end
%{
% 2D rotational scanning
R_rot = 15;
for xx = 1:input_num
    fxi = R_rot*cos(2*pi*(xx/input_num));
    fyi = R_rot*sin(2*pi*(xx/input_num));
    inEx(:,:,xx) = exp( -( (x - fxi*(a/input_dim)).^2 + (y - fyi*(a/input_dim)).^2 ) / focus_waist^2 );
    amplitude = sqrt(dx*dy*abs(sum(sum(inEx(:,:,xx)*(1/imp0).*conj(inEx(:,:,xx))))));
    inEx(:,:,xx) = inEx(:,:,xx) / amplitude;
end
%}  

%% excitaion of MMF modes due to energy coupling
ill_E = reshape(ill_E, [N^2, input_num]);
Er  = reshape( Er, [N^2, NMode]);
Ep  = reshape( Ep, [N^2, NMode]);

sin_M = sparse( diag(sin(theta(:))) );
cos_M = sparse( diag(cos(theta(:))) );

% transform Ex(Ey) into Er and Ep, then use inner product to calculate coupling strength

%%%% H pol. excitation %%%%
% calculate the amplitude, or the coupling coefficient of a mode
% dot product of Ex and modal fields (if to be more stringent, should consider Ez also)
H_to_mode = bsxfun( @rdivide, (Er'*cos_M - Ep'*sin_M), sqrt(diag(Er'*Er) + diag(Ep'*Ep)) );
in_vector_H = H_to_mode * ill_E;
%%%%%%%%%%%%%%%%%%%%%%%%

%%%% V pol. excitation %%%%
% calculate the amplitude, or the coupling coefficient of a mode
% dot product of Ey and modal fields (if to be more stringent, should consider Ez also)
V_to_mode = bsxfun( @rdivide, (Er'*sin_M + Ep'*cos_M), sqrt(diag(Er'*Er) + diag(Ep'*Ep)) );
in_vector_V = V_to_mode * ill_E;   
%%%%%%%%%%%%%%%%%%%%%%%%


%% transform the TM from PIM to focal spot representations  
% coupled modes after MMF propagation
out_vector_H = T*in_vector_H;                                               % this is the output image in PIM representation subject to H polarized input
out_vector_V = T*in_vector_V;                                               % this is the output image in PIM representation subject to V polarized input

mode_to_H = cos_M*Er - sin_M*Ep;
mode_to_V = sin_M*Er + cos_M*Ep;

%%%% H pol. excitation %%%%
img_out_HH = mode_to_H*out_vector_H;
% detect H polarization 
img_out_VH = mode_to_V*out_vector_H;
% detect V polarization


%%%% V pol. excitation %%%%
img_out_HV = mode_to_H*out_vector_V;
% detect H polarization 
img_out_VV = mode_to_V*out_vector_V;
% detect V polarization

%% reshape the output images into matrices
T_HH = img_out_HH;                                                          % the same as mode_to_x * T * x_to_mode * ill_E
T_VH = img_out_VH;
T_HV = img_out_HV;
T_VV = img_out_VV;

end
