%% define the fiber specification, geometry, and operating wavelength 
lambda = 1.55e-6;                                                           % wavelength in air
D = 50e-6;                                                                  % fiber core size
NA = 0.15;                                                                  % NA of fiber
Length = 1;                                                                 % total length of MMF
Rho = inf;                                                                  % radius of curvature of the bending (m)
Theta = 0;                                                                  % orientation of the bending projected on x-y plane
N = 30;                                                                 
input_dim = 30;
input_num = input_dim^2;

[ T_HH, T_HV, T_VH, T_VV, mode_to_H, mode_to_V, H_to_mode, V_to_mode, ill_E ]...
            = MMF_simTM_camera( lambda, D, NA, Length, Rho, Theta, N, input_num );

%% simulate experimental TM measurement of the MMF (consider H polarized inputs)
% To measure the TM of a MMF, a straightforward method is probing individual MMF input spatial channels
%   with a scanning focus and recording the output response subject to the illumination realization.
%   The red circular supports outline the core parts on the MMF facet
inE = abs(reshape( ill_E, [N, N, input_num] )).^2;
output_img_H = reshape(abs(T_HH).^2, [N N input_num]);
output_img_V = reshape(abs(T_VH).^2, [N N input_num]);

close all 
figure('Position', [200, 200, 1000, 400]);
for ii = 1:input_num
    subplot(131)
    imagesc( output_img_H(:,:,ii), [0 1e12] ); colormap gray; axis image
    title('H polarized ligh at output');
    viscircles( [N/2 N/2], 0.94*round(N/2) );
    
    subplot(132)
    imagesc( output_img_V(:,:,ii), [0 1e12] ); colormap gray; axis image
    title('V polarized light at output');
    viscircles( [N/2 N/2], 0.94*round(N/2) );
    
    subplot(133)
    imagesc( inE(:,:,ii) ); colormap gray; axis image
    title('H polarized scanning focus at input');
    viscircles( [N/2 N/2], 0.94*round(N/2) );
    drawnow
end

%% generate a H polarized focus through the MMF using measured T_HH
% To generate a focus, we first compute the regularized inversion of T_HH. The input wave-fronts necessary to
%   create foci at defined output positions correspond to column vectors of the approximately inverted matrix.
invT_HH = Tikinv( T_HH ); invT_HH = invT_HH/max(invT_HH(:));
output_img_H = reshape( abs(T_HH*invT_HH).^2, [N N N^2] );
output_img_V = reshape( abs(T_VH*invT_HH).^2, [N N N^2] );

close all
figure('Position', [200, 200, 1000, 400]);
for ii = 1:N^2
    subplot(131)
    imagesc( output_img_H(:,:,ii), [0 4e14] ); colormap gray; axis image
    title('H polarized ligh at output');
    viscircles( [N/2 N/2], 0.94*round(N/2) );
    
    subplot(132)
    imagesc( output_img_V(:,:,ii), [0 4e14] ); colormap gray; axis image
    title('V polarized light at output');
    viscircles( [N/2 N/2], 0.94*round(N/2) );
    
    subplot(133)
    complex_imagesc( reshape(invT_HH(:,ii), [N N]) ); colormap gray; axis image
    title('H polarized wavefront at input');
    viscircles( [N/2 N/2], 0.94*round(N/2) );
    drawnow
end

%% compute the new TM of the MMF when the fiber shape is deformed
Length = 1e-1*ones(1, 10);                                                  % 10 segments, total length adds up to 1m, the same as the original setting
Rho = 0.25 + 0.1 * (rand(1, numel(Length)) - 0.5);                           % radius of curvature of the bending between 0.2 - 0.3 m
Theta = (2*pi) * rand(1, numel(Length));                                    % orientation of the bending projected on x-y plane

[ T_HH_new, T_HV_new, T_VH_new, T_VV_new ] = MMF_simTM_ccd( lambda, D, NA, Length, Rho, Theta, N, input_num );

%% simulate distorted focus when the fiber shape is deformed
output_img_H_new = reshape( abs(T_HH_new*invT_HH).^2, [N N N^2] );
output_img_V_new = reshape( abs(T_VH_new*invT_HH).^2, [N N N^2] );

close all
figure('Position', [200, 200, 1000, 400]);
for ii = 1:N^2
    subplot(131)
    imagesc( output_img_H_new(:,:,ii), [0 3e14] ); colormap gray; axis image
    title('H polarized ligh at output');
    viscircles( [N/2 N/2], 0.94*round(N/2) );
    
    subplot(132)
    imagesc( output_img_V_new(:,:,ii), [0 3e14] ); colormap gray; axis image
    title('V polarized light at output');
    viscircles( [N/2 N/2], 0.94*round(N/2) );
    
    subplot(133)
    complex_imagesc( reshape(invT_HH(:,ii), [N N]) ); colormap gray; axis image
    title('H polarized wavefront at input');
    viscircles( [N/2 N/2], 0.94*round(N/2) );
    drawnow
end

