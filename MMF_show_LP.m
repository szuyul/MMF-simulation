%% define the fiber specification, geometry, and operating wavelength 
lambda = 1.55e-6;                                                           % wavelength in air
D = 25e-6;                                                                  % fiber core size
NA = 0.2;                                                                  % NA of fiber
Length = 1;                                                                 % total length of MMF
Rho = inf;                                                                  % radius of curvature of the bending (m)
Theta = 0;                                                                  % orientation of the bending projected on x-y plane
N = 30;                                                                 

[ T, NMode, lmap, mmap, LPxymap, propconst, Ex, Ey, Ez, Hx, Hy, Hz, img_size ]...
            = MMF_simTM_LP( lambda, D, NA, Length, Rho, Theta, N );

%% show the propagation constants
LPlmax = max(abs( lmap ));

close all
scatter_size = 20;
figure('Position', [100, 100, 800, 300]);                                   % plot the propagation constant for EH
for l = 0:LPlmax                                                            % starts from (l, m) = (0, 1)
    temp = propconst( lmap == l );
    LPm = numel(temp);                                                      % # of modes with that l
    scatter(l*ones(1,LPm),temp,20,'black', 'filled');
    hold on
end
title('LP modes');  grid on
xlabel('orbital angular momentum ( l index)');  ylabel('propagation const. (\beta in m^-1)')

%% plot the 2D modal E fields of LP modes (x polarization for example)

close all
figure('Position', [100, 100, 1100, 400]);

for ii = 1:(NMode/2)
    temp = Ex(:,:,ii);                                                      
    pind = (mmap(ii)-1)*(2*LPlmax+1) + lmap(ii) + LPlmax + 1;
    
    subplot(max(mmap(LPxymap==1)), 2*LPlmax+1, pind);
    complex_imagesc( temp );
    title(['LPx l = ',num2str(lmap(ii)),', m = ',num2str(mmap(ii))]);
    axis off
end

%% plot the 2D Laguerre-Gaussian modes by superposing LP modes (x polarization for example)
close all
figure('Position', [100, 100, 1200, 600]);

for l = 0:max(lmap)
    idx_pm = find( lmap == l & LPxymap == 1 );
    idx_nm = find( lmap ==-l & LPxymap == 1 );
    
    for ii = 1:numel(idx_pm)
        temp = Ex(:,:,idx_pm(ii)) + Ex(:,:,idx_nm(ii));                                                      
        imgmax = max(max(abs( temp )));
        pind = (ii-1)*(max(lmap)+1) + (l+1);
    
        subplot(max(mmap), max(lmap)+1, pind);
        complex_imagesc( temp/imgmax );
        title(['LG (',num2str(ii-1),',',num2str(l), ')']);
        axis off
    end
end

