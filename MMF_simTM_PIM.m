%MMF_simTM_PIM  simulates the vectorial propagation-invariant modes (PIMs), a set of Bessel functions, 
% of a MMF given the specifications and geometry. The code is the implementation of theoretical MMF model,
% which can be found in the 3rd chapter of "Optical Electronics in Modern Communications, 1997" by Prof. Amnon Yariv.
% The fiber transmission of a deformed MMF is simulated based on a perturbation theory developed by Ploschner et al., 
% "Seeing through chaos in multimode fibres," Nat. Photonics 9, 2015.
%
% [ T, NMode, lmap, mmap, EHHEmap, propconst, Er, Ep, Ez, Hr, Hp, Hz, img_size ]...
%            = MMF_simTM_PIM( lambda, D, NA, Length, Rho, Theta, N ) 
% 
% outputs:
% T is the transmission matrix of a straight or deformed MMF in PIM representation at a certain wavelength
% NMode is the number of total PIM modes
% lmap is the orbital angular momentum, or l index, of each PIM mode
% mmap is the m index of each PIM mode
% EHHEmap specifies whether a mode is either a EH or HE mode, 1 for EH and 0 for HE
% propconst is the propagation constant (unit: m) of each PIM mode
% Er, Ep, Ez (N by N by NMode) is the transverse electric field distribution of each PIM mode in radial, azimuthal, and axial directions
% Hr, Hp, Hz (N by N by NMode) is the transverse magnetic field distribution of each PIM mode
% img_size is the physical size (unit: m) of the transverse field images
%
% inputs:
% lambda is the operating wavelength (unit: m)
% D is the fiber core diameter of the MMF (unit: m)
% NA is the numerical aperture of the MMF 
% Length is an array of the segmental length of the MMF (unit: m)
% Rho is an array of the bending radius of curvature (unit: m) of each defined segment. If the segment is straight, Rho = inf
% Theta is an array of the orientation of the bending of the MMF (unit: rad.)
% N is the image dimension
%
%
% 2017-2019 Szu-Yu Lee
% Bouma Lab - The Wellman Center for Photomedicine

function [ T, NMode, lmap, mmap, EHHEmap, propconst, Er, Ep, Ez, Hr, Hp, Hz, img_size ]...
            = MMF_simTM_PIM( lambda, D, NA, Length, Rho, Theta, N )
%% MMF parameter initialization (unit in m)
c = 3e8;
epsilon = 8.85e-12;
k0 = 2*pi/lambda;                                                           % wavenumber
w = 2*pi*c/lambda;
mu = 4*pi*1e-7;
%imp0 = w*mu/k0;

a = D/2;                                                                    % radius of the fiber core

% refractive index of the MMF core based on fused silica, please see:
% https://refractiveindex.info/?shelf=glass&book=fused_silica&page=Malitson                                      
n_core = 1.4469;                                                            % refractive index of the core
n_clad = sqrt(n_core^2 - NA^2);                                             % refractive index of the cladding 
%delta_n = n_core- n_clad;

epsilon_core = epsilon*n_core^2;
epsilon_clad = epsilon*n_clad^2;
k_core = k0*n_core;
%k_clad = k0*n_clad;
imp_core = w*mu/k_core;
%imp_clad = w*mu/k_clad;

%V = k0*a*NA;                                                               % V number of the fibre (also known as waveguide or fibre parameter)
nndif = (n_core^2 - n_clad^2)/2/n_core^2;                                   % normalized refractive index difference
nnsum = (n_core^2 + n_clad^2)/2/n_core^2;                                   % normalized refractive index summation

%% Find roots of characteristic equation to define propagation constants of the EH and HE modes
dbn = 8e4; 
% how fine the root will be found (the # of modes should saturate when dbn large enough and all roots found)                                                                 
   
Beta = linspace(n_core*k0, n_clad*k0, dbn);                                 % from large to small prop. const (low to high order)
ha = @(x) a*sqrt( (n_core*k0)^2 - x.^2 );
qa = @(x) a*sqrt( x.^2 - (n_clad*k0)^2 );

L = 50;                                                                     
% L is the maximal possible absolute of orbital angular momentum. Increase this if the MMF has a very large # of modes 

EHbeta = zeros(L + 1, L);
EHm = zeros(L,1);                                                           % # of modes of each l index 
HEbeta = zeros(L + 1, L);
HEm = zeros(L,1);

for l=0:L                                                                   % for each l index, find roots for beta       
    %Jd = @(ri) besselj(l-1,ri) - (l./ri).*besselj(l,ri);
    Kd = @(ro) (l./ro).*besselk(l,ro) - besselk(l+1,ro);
    R  = @(x) sqrt( (nndif^2)*(Kd(qa(x))./(qa(x).*besselk(l,qa(x)))).^2 ...
                  + (l*x/(n_core*k0)).^2 .* ( qa(x).^(-2) + ha(x).^(-2) ) .^2 );
    %%% EH modes %%%
    EH_transcend_func = @(x) besselj(l+1,ha(x))./(ha(x).*besselj(l,ha(x)))...
                            - nnsum*(Kd(qa(x))./(qa(x).*besselk(l,qa(x)))) - (l./(ha(x).^2) - R(x));
    EH = EH_transcend_func( Beta );
    EH_transcend_func_abs = @(x) abs( EH_transcend_func(x) ).^2;
    
    EHind = intersect( find(abs(diff(sign(EH))) > 0), find(abs(EH) < 3), 'sorted');
    % approximate the position of roots of the characteristic equation
    % coarse search for the solution with criterion: changed sign with small magnitude to avoid abrupt jump
    
    beta = propconst_fine_search( EH_transcend_func_abs, Beta(EHind), Beta(EHind+2), Beta(EHind-2));
    % fine search within the interval bounded by consecutive roots found in coarse search
    
    EHbeta(l+1,1:length(EHind)) = beta;   
    % collect beta of m modes for a l 
    EHm(l+1) = sum(beta ~= 0);                                              % keep track of # of modes for each l index
    %%%%%%%%%%%%%%%%
    
    %%% HE modes %%%
    HE_transcend_func = @(x) besselj(l-1,ha(x))./(ha(x).*besselj(l,ha(x)))...
                            + nnsum*(Kd(qa(x))./(qa(x).*besselk(l,qa(x)))) - (l./(ha(x).^2) - R(x));
    HE = HE_transcend_func( Beta );
    HE_transcend_func_abs = @(x) abs( HE_transcend_func(x) ).^2;
    
    HEind = intersect( find(abs(diff(sign(HE))) > 0), find(abs(HE) < 3), 'sorted');
    % coarse search for the solution with criterion: changed sign with small magnitude to avoid abrupt jump
    beta = propconst_fine_search( HE_transcend_func_abs, Beta(HEind), Beta(HEind+2), Beta(HEind-2));
    % fine search within the interval bounded by consecutive roots found in coarse search
    HEbeta(l+1,1:length(HEind)) = beta;
    HEm(l+1) = sum(beta ~= 0);
    %%%%%%%%%%%%%%%%
end

%% sort the modes and index them
EHlmax = find(EHm ~= 0, 1, 'last' ) - 1;                                    % EH mode largest l number, -1 since l starts from 0
HElmax = find(HEm ~= 0, 1, 'last' ) - 1;

NEH = sum(EHbeta(1,:) ~= 0) + 2*sum(sum(EHbeta(2:end,:)~=0));               % the number of EH modes, including negative l 
NHE = sum(HEbeta(1,:) ~= 0) + 2*sum(sum(HEbeta(2:end,:)~=0));               % the number of HE modes, including negative l 
NMode = NEH + NHE;                    
% keep track of mode indices (EH or HE, l = -m...,-3,-2,-1,0,1,2,3...,m)
EHHEmap = vertcat( ones(NEH,1), zeros(NHE,1) );                             % 1 for EH modes, 0 for HE modes
lmap = [];
mmap = [];
for l = -EHlmax:EHlmax
    ind = abs(l) + 1;
    lmap = vertcat(lmap, l*ones(EHm(ind),1));
    mmap = vertcat(mmap, (1:EHm(ind)).');
end

for l = -HElmax:HElmax
    ind = abs(l) + 1;
    lmap = vertcat(lmap, l*ones(HEm(ind),1));
    mmap = vertcat(mmap, (1:HEm(ind)).');
end

%% derive "vector" E field (EH & HE modes) for each orbital angular momentum in polar coordinates
img_size = 1.06*D;                                                           
% physical size of the image covering fiber facet, leaving some margin for evanescent field in the cladding part

[x, y] = meshgrid(linspace(-img_size/2, img_size/2, N),linspace(-img_size/2, img_size/2, N));
dx = x(1,2) - x(1,1);
dy = y(2,1) - y(1,1);
[theta,rho] = cart2pol(x,y);
% create a N x N 2D polar coordinate

% preallocate the EM fields for each mode
Er = zeros(N, N, NMode);                                                    % radial E in a 2D grid, mode index
Ep = zeros(N, N, NMode);                                                    % azimuthal E
Ez = zeros(N, N, NMode);                                                    % axial E
Hr = zeros(N, N, NMode);                   
Hp = zeros(N, N, NMode);                   
Hz = zeros(N, N, NMode); 

% calculate the EM field for each mode given (l,m) 
% when l ~= 0, the negative l is also a degeneracy
for ii = 1:NMode
    l = lmap(ii);
    m = mmap(ii);
    if EHHEmap(ii) == 1
        mode_beta = EHbeta( abs(l)+1, m );
    else
        mode_beta = HEbeta( abs(l)+1, m );
    end
    [ Er(:,:,ii), Ep(:,:,ii), Ez(:,:,ii), Hr(:,:,ii), Hp(:,:,ii), Hz(:,:,ii) ] ...
         = MMFEHHEfields( a, n_core, epsilon_core, n_clad, epsilon_clad, k0, w, EHHEmap(ii), mode_beta, l, rho, theta, 0 );
    
    amplitude = sqrt( dx*dy*abs(sum(sum( Er(:,:,ii).*conj(Hp(:,:,ii)) - Ep(:,:,ii).*conj(Hr(:,:,ii)) ))) );
    % normalize the propagating power along z axis of each mode based on Poynting theorem
        
    Er(:,:,ii) = Er(:,:,ii) / amplitude;
    Ep(:,:,ii) = Ep(:,:,ii) / amplitude;
    Ez(:,:,ii) = Ez(:,:,ii) / amplitude;
    Hr(:,:,ii) = Hr(:,:,ii) / amplitude;
    Hp(:,:,ii) = Hp(:,:,ii) / amplitude;
    Hz(:,:,ii) = Hz(:,:,ii) / amplitude;
end

Er_vec = reshape(Er, [N^2, NMode]);
Ep_vec = reshape(Ep, [N^2, NMode]);

%% construction of transmission matrix in PIM representation based on EH, HE solutions
% each EH(l,m) and HE(l,m) has corresponding propagation constant                                                                                           
propconst = zeros(NMode, 1);
% propagation constants ordered in (l,m) = (0,1), (0,2), (0,3)...(1,1), (-1,1), (1,2), (-1,2)...

for ii = 1:NMode
    l = lmap(ii);
    m = mmap(ii);
    if EHHEmap(ii) == 1
        propconst(ii) = EHbeta( abs(l)+1 , m );
    else
        propconst(ii) = HEbeta( abs(l)+1 , m );
    end
end
B0 = diag(propconst);

% allocate the T matrix
n = numel(Length);                                                          % how many pieces of MMF segment
T = eye(NMode);

%%%%%%% mode dependent loss %%%%%%%%
% define the mode dependent loss, otherwise make it all 0
MDL = zeros(NMode, 1);
loss = 0;
%loss = linspace(0.1/n, 0.1/n, NMode/2).';                                   % low to high loss mode dependent loss (MDL)
[~, ind] = sort( propconst(1:NMode/2),        'descend' );
MDL(ind) = (1-loss);                                        
[~, ind] = sort( propconst(NMode/2 + 1: end), 'descend' );
MDL(NMode/2 + ind) = (1-loss);  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sigma = 0.17;
xi = 1 - (1 - 2*sigma)*((n_core - 1)/n_core);

for pp = 1:n                                                                % iterate over each piece of MMF segment
    if Rho(pp) == inf                                                       % if there is no bending
        T = diag( MDL .* exp(1i*propconst*Length(pp)) ) * T;
    else
%% build the perturbed T matrix, Tp, by calculated modes
% for details, please see Ploschner et al., "Seeing through chaos in multimode fibres," Nat. Photonics 9, 2015.
        Rho_eff = Rho(pp)/xi;                                               % effective radius of curvature considering refractive index variation
        X = reshape( cos(Theta(pp))*x + sin(Theta(pp))*y, [N^2, 1]);        % perturbation term for that segment of MMF
        A = dx*dy*(1/imp_core)* ( Er_vec'*bsxfun(@times, X, Er_vec)...      % implement the spatial operator applied to individual modes
                                + Ep_vec'*bsxfun(@times, X, Ep_vec) );      % pay attention to the normalization
                            
        B = B0 - (n_core*k0/Rho_eff)*A;
        T = (diag(MDL) * expm(1i*B*Length(pp))) * T;                        % note the difference between exp and expm, the former is element-wise exponential, 
                                                                            % the latter is matrix exponential 
    end
end



end
