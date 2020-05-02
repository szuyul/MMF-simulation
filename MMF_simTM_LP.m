%MMF_simTM_LP  simulates the vectorial linearly polarized modes (LPs), a set of Bessel functions, 
% of a MMF given the specifications and geometry. The LPs are the approximation of PIMs under the 
% assumption of weakly guided modes. The code is the implementation of theoretical MMF model,
% which can be found in the 3rd chapter of "Optical Electronics in Modern Communications, 1997" by Prof. Amnon Yariv.
% The fiber transmission of a deformed MMF is simulated based on a perturbation theory developed by Ploschner et al., 
% "Seeing through chaos in multimode fibres," Nat. Photonics 9, 2015.
%
% [ T, NMode, lmap, mmap, LPxymap, propconst, Ex, Ey, Ez, Hx, Hy, Hz, img_size ]...
%            = MMF_simTM_LP( lambda, D, NA, Length, Rho, Theta, N )
% 
% outputs:
% T is the transmission matrix of a straight or deformed MMF in LP representation at a certain wavelength
% NMode is the number of total LP modes
% lmap is the orbital angular momentum, or l index, of each LP mode
% mmap is the m index of each LP mode
% LPxymap specifies whether a mode is either a horizontally(x) or vertically(y) polarized mode, 1 for x and 0 for y
% propconst is the propagation constant (unit: m) of each LP mode
% Ex, Ey, Ez (N by N by NMode) is the transverse electric field distribution of each LP mode in x, y, and z directions
% Hx, Hy, Hz (N by N by NMode) is the transverse magnetic field distribution of each LP mode
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

function [ T, NMode, lmap, mmap, LPxymap, propconst, Ex, Ey, Ez, Hx, Hy, Hz, img_size ]...
            = MMF_simTM_LP( lambda, D, NA, Length, Rho, Theta, N )
%% MMF parameter initialization (unit in m)
c = 3e8;
%epsilon = 8.85e-12;
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

%epsilon_core = epsilon*n_core^2;
%epsilon_clad = epsilon*n_clad^2;
k_core = k0*n_core;
%k_clad = k0*n_clad;
imp_core = w*mu/k_core;
%imp_clad = w*mu/k_clad;

%V = k0*a*NA;                                                                % V number of the fibre (also known as waveguide or fibre prameter)
%nndif = (n_core^2 - n_clad^2)/2/n_core^2;                                   % normalized refractive index difference
%nnsum = (n_core^2 + n_clad^2)/2/n_core^2;                                   % normalized refractive index summation

%% Find roots of characteristic equation to define propagation constants of the weakly guided LP modes
dbn = 8e4;                                                                  
% how fine the root will be found (the # of modes should saturate when dbn large enough and all roots found)

Beta = linspace(n_core*k0, n_clad*k0, dbn);                                 % from large to small prop. const (low to high order)
ha = @(x) a*sqrt( (n_core*k0)^2 - x.^2 );
qa = @(x) a*sqrt( x.^2 - (n_clad*k0)^2 );

L = 50;                                                                     
% L is the maximal possible absolute of orbital angular momentum. Increase this if the MMF has a very large # of modes 

LPbeta = zeros(L + 1, L);
LPm = zeros(L,1);                                                           % # of modes of each l index 

for l=0:L                                                                   % for each orbital angular momentum, find roots for beta     
    LP_transcend_func = @(x) ha(x).*( besselj(l+1,ha(x))./besselj(l,ha(x)) )...
                           - qa(x).*( besselk(l+1,qa(x))./besselk(l,qa(x)) );
    LP = LP_transcend_func( Beta );
    LP_transcend_func_abs = @(x) abs( LP_transcend_func(x) ).^2;
    % approximate position of roots of the characteristic equation
    % coarse search: sign changed and magnitude smaller than an abrupt jump
    
    LPind = intersect( find(abs(diff(sign(LP))) > 0), find(abs(LP) < 3), 'sorted');                                                  
    % fine search, for each approx. root, use prior and later element as lower and upper search bounds
    
    beta = propconst_fine_search( LP_transcend_func_abs, Beta(LPind), Beta(LPind+2), Beta(LPind-2)); 
    LPbeta(l+1,1:length(LPind)) = beta;                                     % collect beta of m modes for a l 
    LPm(l+1) = sum(beta ~= 0);                                              % keep track of # of m modes for l
end

%% sort the modes and index them
LPlmax = find(LPm ~= 0, 1, 'last' ) - 1;                                    % LP mode largest l number, -1 since l starts from 0

NLP = sum(LPbeta(1,:) ~= 0) + 2*sum(sum(LPbeta(2:end,:)~=0));               % the number of LP modes per polarization, including negative l 
NMode = 2*NLP;                    
% keep track of mode indices (EH or HE, l = -m...,-3,-2,-1,0,1,2,3...,m)
LPxymap = vertcat( ones(NLP,1), zeros(NLP,1) );                             % 1 for LPx modes, 0 for LPy modes
lmap = [];
mmap = [];
for l = -LPlmax:LPlmax                                                      % for LPx modes
    ind = abs(l) + 1;
    lmap = vertcat(lmap, l*ones(LPm(ind),1));
    mmap = vertcat(mmap, (1:LPm(ind)).');
end

for l = -LPlmax:LPlmax                                                      % for LPy modes
    ind = abs(l) + 1;
    lmap = vertcat(lmap, l*ones(LPm(ind),1));
    mmap = vertcat(mmap, (1:LPm(ind)).');
end

%%  derive "vector" E field (LPx & LPy) for each orbital angular momentum in polar coordinates
img_size = 1.06*D;
% physical size of the image covering fiber facet, leaving some margin for evanescent field in the cladding part

[x, y] = meshgrid(linspace(-img_size/2, img_size/2, N),linspace(-img_size/2, img_size/2, N));
dx = x(1,2) - x(1,1);
dy = y(2,1) - y(1,1);
[theta,rho] = cart2pol(x,y);
% create a N x N 2D polar coordinate

% preallocate the LP fields for each mode
Ex = zeros(N, N, NMode);                                                   
Ey = zeros(N, N, NMode);                                                    
Ez = zeros(N, N, NMode);                                                    
Hx = zeros(N, N, NMode);                   
Hy = zeros(N, N, NMode);                   
Hz = zeros(N, N, NMode); 

% calculate the LP field for each mode given (l,m)
% when l ~= 0, the negative l is also a degeneracy
for ii = 1:NMode
    l = lmap(ii);
    m = mmap(ii);
    mode_beta = LPbeta( abs(l)+1, m );
    
    [ Ex(:,:,ii), Ey(:,:,ii), Ez(:,:,ii), Hx(:,:,ii), Hy(:,:,ii), Hz(:,:,ii) ] ...
        = MMFLPfields( a, n_core, n_clad, k0, w, mode_beta, l, LPxymap(ii), rho, theta, 0 );

    amplitude = sqrt( dx*dy*abs(sum(sum( Ex(:,:,ii).*conj(Hy(:,:,ii)) - Ey(:,:,ii).*conj(Hx(:,:,ii)) ))) );  
    % normalize the propagating power along z axis of each mode based on Poynting theorem
    
    Ex(:,:,ii) = Ex(:,:,ii) / amplitude;
    Ey(:,:,ii) = Ey(:,:,ii) / amplitude;
    Ez(:,:,ii) = Ez(:,:,ii) / amplitude;
    Hx(:,:,ii) = Hx(:,:,ii) / amplitude;
    Hy(:,:,ii) = Hy(:,:,ii) / amplitude;
    Hz(:,:,ii) = Hz(:,:,ii) / amplitude;
end

Ex_vec = reshape(Ex, [N^2, NMode]);
Ey_vec = reshape(Ey, [N^2, NMode]);


%% construction of TransM based on EH, HE solutions
% each LP(l,m) has corresponding propagation constant                                                                                        
propconst = zeros(NMode, 1);
% propagation constants ordered in (l,m) = (0,1), (0,2), (0,3)...(1,1), (-1,1), (1,2), (-1,2)...

for ii = 1:NMode
    l = lmap(ii);
    m = mmap(ii);
    if LPxymap(ii) == 1
        propconst(ii) = LPbeta( abs(l)+1 , m );
    else
        propconst(ii) = LPbeta( abs(l)+1 , m );
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
        X = reshape( cos(Theta(pp))*x + sin(Theta(pp))*y ,[N^2, 1]);        % perturbation term for that segment of MMF
        A = dx*dy*(1/imp_core)* ( Ex_vec'*bsxfun(@times, X, Ex_vec)...
                                + Ey_vec'*bsxfun(@times, X, Ey_vec) );
        B = B0 - (n_core*k0/Rho_eff)*A;
        T = (diag(MDL) * expm(1i*B*Length(pp))) * T;                        % note the difference between exp and expm
    end
end



end