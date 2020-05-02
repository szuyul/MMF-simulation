%% define the fiber specification, geometry, and operating wavelength 
lambda = 1.55e-6;                                                           % wavelength in air
D = 25e-6;                                                                  % fiber core size
NA = 0.1;                                                                  % NA of fiber
Length = 1;                                                                 % total length of MMF
Rho = inf;                                                                  % radius of curvature of the bending (m)
Theta = 0;                                                                  % orientation of the bending projected on x-y plane
N = 30;                                                                 

[ T, NMode, lmap, mmap, EHHEmap, propconst, Er, Ep, Ez, Hr, Hp, Hz, img_size ]...
            = MMF_simTM_PIM( lambda, D, NA, Length, Rho, Theta, N );

%% show the propagation constants
EHlmax = max(abs( lmap(EHHEmap == 1) ));
HElmax = max(abs( lmap(EHHEmap == 0) ));

close all
figure('Position', [100, 100, 800, 300]);                                   % plot the propagation constant for EH
subplot(1,2,1)
for l = 0:EHlmax                                                            % starts from (l, m) = (0, 1)
    temp = propconst( lmap == l & EHHEmap == 1 );
    EHm = numel(temp);                                                      % # of modes with that l
    scatter(l*ones(1,EHm),temp, 20,'black', 'filled');
    hold on
end
title('EH modes');  grid on
xlabel('orbital angular momentum ( l index)');  ylabel('propagation const. (\beta in m^-1)')

subplot(1,2,2)                                                              % plot the propagation constant for HE
for l = 0:HElmax
    temp = propconst( lmap == l & EHHEmap == 0 );
    HEm = numel(temp);
    scatter(l*ones(1,HEm),temp, 20,'black', 'filled');
    hold on
end
title('HE modes');  grid on
xlabel('orbital angular momentum ( l index)');  ylabel('propagation const. (\beta in m^-1)')

%% 2D plot the modal E fields of PIMs
close all
figure('Position', [100, 500, 1100, 400]);
figure('Position', [100, 100, 1100, 400]);

for ii = 1:NMode
    temp = Er(:,:,ii);                                                      % choose a scalar field to plot
    if EHHEmap(ii) == 1
        figure(1)
        pind = (mmap(ii)-1)*(2*EHlmax+1) + lmap(ii) + EHlmax + 1;
        subplot(max(mmap(EHHEmap==1)), 2*EHlmax+1, pind);
        complex_imagesc( temp );
        title(['EH (Er) l = ',num2str(lmap(ii)),', m = ',num2str(mmap(ii))]);
    else
        figure(2)
        pind = (mmap(ii)-1)*(2*HElmax+1) + lmap(ii) + HElmax + 1;
        subplot(max(mmap(EHHEmap==0)), 2*HElmax+1, pind);
        complex_imagesc( temp );
        title(['HE (Er) l = ',num2str(lmap(ii)),', m = ',num2str(mmap(ii))]); 
    end
    axis off
end

