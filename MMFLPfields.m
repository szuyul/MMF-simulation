%MMFLPfields calculates the vectorial fields of MMF linearly polarized (LP) modes
% The code is the implementation of theoretical MMF model, which can be found in the 3rd chapter of 
% "Optical Electronics in Modern Communications, 1997" by Prof. Amnon Yariv.
%
% outputs:
% Ex, Ey, Ez (N by N by NMode) is the transverse electric field distribution of each LP mode in horizontal, vertical, and axial directions
% Hx, Hy, Hz (N by N by NMode) is the transverse magnetic field distribution of each LP mode
%
% inputs:
% a is the radius of the fiber core
% n_core, n_clad are the refractive indices of core and cladding parts
% k0, w is the wavenumber and frequency in vacuum
% beta is the propagation constant (unit: m)
% l is the orbital angular momentum
% LPxymap is an indicator of whether the mode is x(1) or y(0) polarized
% rho, theta are the 2D polar coordinates
% z is the axial coordinate
%
% 
% 2017-2019 Szu-Yu Lee
% Bouma Lab - The Wellman Center for Photomedicine

function [ Ex, Ey, Ez, Hx, Hy, Hz ] = MMFLPfields( a, n_core, n_clad, k0, w, beta, l, LPxymap, rho, theta, z )
mu = 4*pi*1e-7;
h = sqrt( (n_core*k0)^2 - beta^2 );
q = sqrt( beta^2 - (n_clad*k0)^2 );
ha = h*a;
qa = q*a;

N = size(rho,1);
Ex = zeros(N, N);                   
Ey = zeros(N, N);                   
Ez = zeros(N, N);                   
Hx = zeros(N, N);                   
Hy = zeros(N, N);                   
Hz = zeros(N, N);

B = besselj(l, ha)./besselk(l, qa);
%%%%%%%%%%%%%%%%%
if LPxymap == 1                     % for LPx modes
% calculate the fields within the core
Excore = @(r,pha) besselj(l,h*r).*exp(1i*(l*pha - beta*z));
Eycore = @(r,pha) 0;
Ezcore = @(r,pha) (1i*h/(2*beta)).*(...
                                besselj(l+1,h*r).*exp(1i*pha) - besselj(l-1,h*r).*exp(-1i*pha) ).*...
                                        exp(1i*(l*pha - beta*z));
Hxcore = @(r,pha) 0;
Hycore = @(r,pha) (beta/(w*mu)).*besselj(l,h*r).*exp(1i*(l*pha - beta*z));
Hzcore = @(r,pha) (h/(2*w*mu)).*(...
                                besselj(l+1,h*r).*exp(1i*pha) + besselj(l-1,h*r).*exp(-1i*pha) ).*...
                                        exp(1i*(l*pha - beta*z));
% calculate the fields within the cladding
Exclad = @(r,pha) B.*besselk(l,q*r).*exp(1i*(l*pha - beta*z));
Eyclad = @(r,pha) 0;
Ezclad = @(r,pha) (1i*q*B/(2*beta)).*(...
                                besselk(l+1,q*r).*exp(1i*pha) + besselk(l-1,q*r).*exp(-1i*pha) ).*...
                                        exp(1i*(l*pha - beta*z));
Hxclad = @(r,pha) 0;
Hyclad = @(r,pha) (beta*B/(w*mu)).*besselk(l,q*r).*exp(1i*(l*pha - beta*z));
Hzclad = @(r,pha) (q*B/(2*w*mu)).*(...
                                besselk(l+1,q*r).*exp(1i*pha) - besselk(l-1,q*r).*exp(-1i*pha) ).*...
                                        exp(1i*(l*pha - beta*z));

else                                % for LPy modes
% calculate the fields within the core
Excore = @(r,pha) 0;
Eycore = @(r,pha) besselj(l,h*r).*exp(1i*(l*pha - beta*z));
Ezcore = @(r,pha) (h/(2*beta)).*(...
                                besselj(l+1,h*r).*exp(1i*pha) + besselj(l-1,h*r).*exp(-1i*pha) ).*...
                                        exp(1i*(l*pha - beta*z));
Hxcore = @(r,pha) -(beta/(w*mu)).*besselj(l,h*r).*exp(1i*(l*pha - beta*z));
Hycore = @(r,pha) 0;
Hzcore = @(r,pha) -(1i*h/(2*w*mu)).*(...
                                besselj(l+1,h*r).*exp(1i*pha) - besselj(l-1,h*r).*exp(-1i*pha) ).*...
                                        exp(1i*(l*pha - beta*z));
% calculate the fields within the cladding
Exclad = @(r,pha) 0;
Eyclad = @(r,pha) B.*besselk(l,q*r).*exp(1i*(l*pha - beta*z));
Ezclad = @(r,pha) (q*B/(2*beta)).*(...
                                besselk(l+1,q*r).*exp(1i*pha) - besselk(l-1,q*r).*exp(-1i*pha) ).*...
                                        exp(1i*(l*pha - beta*z));
Hxclad = @(r,pha) -(beta*B/(w*mu)).*besselk(l,q*r).*exp(1i*(l*pha - beta*z));
Hyclad = @(r,pha) 0;
Hzclad = @(r,pha) -(1i*q*B/(2*w*mu)).*(...
                                besselk(l+1,q*r).*exp(1i*pha) + besselk(l-1,q*r).*exp(-1i*pha) ).*...
                                        exp(1i*(l*pha - beta*z));
end
%%%%%%%%%%%%%%%%%
% no division of zero, so do not have to worry of origin singularity
% assigning the fields within the core
region = find( rho < a );
Ex( region ) = Excore( rho( region ), theta( region ) );
Ey( region ) = Eycore( rho( region ), theta( region ) );
Ez( region ) = Ezcore( rho( region ), theta( region ) );
Hx( region ) = Hxcore( rho( region ), theta( region ) );
Hy( region ) = Hycore( rho( region ), theta( region ) );
Hz( region ) = Hzcore( rho( region ), theta( region ) );

% assigning the fields within the cladding
region = find( rho > a );
Ex( region ) = Exclad( rho( region ), theta( region ) );
Ey( region ) = Eyclad( rho( region ), theta( region ) );
Ez( region ) = Ezclad( rho( region ), theta( region ) );
Hx( region ) = Hxclad( rho( region ), theta( region ) );
Hy( region ) = Hyclad( rho( region ), theta( region ) );
Hz( region ) = Hzclad( rho( region ), theta( region ) );

end

