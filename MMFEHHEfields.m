%MMFEHHEfields calculates the vectorial fields of MMF propagation invariant modes (PIMs)
% The code is the implementation of theoretical MMF model, which can be found in the 3rd chapter of 
% "Optical Electronics in Modern Communications, 1997" by Prof. Amnon Yariv.
%
% outputs:
% Er, Ep, Ez (N by N by NMode) is the transverse electric field distribution of each PIM mode in radial, azimuthal, and axial directions
% Hr, Hp, Hz (N by N by NMode) is the transverse magnetic field distribution of each PIM mode
%
% inputs:
% a is the radius of the fiber core
% n_core, n_clad are the refractive indices of core and cladding parts
% epsilon, epsilon2 are the permittivity of the fiber core and cladding parts
% k0, w is the wavenumber and frequency in vacuum
% EHorHE is an indicator of whether the mode is EH(1) or HE(0)
% beta is the propagation constant (unit: m)
% l is the orbital angular momentum
% rho, theta are the 2D polar coordinates
% z is the axial coordinate
%
% 
% 2017-2019 Szu-Yu Lee
% Bouma Lab - The Wellman Center for Photomedicine

function [ Er, Ep, Ez, Hr, Hp, Hz ] = MMFEHHEfields( a, n_core, epsilon, n_clad, epsilon2, k0, w, EHorHE, beta, l, rho, theta, z )
mu = 4*pi*1e-7;
h = sqrt( (n_core*k0)^2 - beta^2 );
q = sqrt( beta^2 - (n_clad*k0)^2 );
ha = h*a;
qa = q*a;

N = size(rho,1);
Er = zeros(N, N);                   % radial E 2D grid sizes, mode index
Ep = zeros(N, N);                   % azimuthal E
Ez = zeros(N, N);                   % axial E
Hr = zeros(N, N);                   
Hp = zeros(N, N);                   
Hz = zeros(N, N);

Jd = @(r) besselj(l-1,r) - (l./r).*besselj(l,r);
Kd = @(r) (l./r).*besselk(l,r) - besselk(l+1,r);

if l == 0                           % corresponds to TE01, TE02... or TM01, TM02... 
    if EHorHE == 1                  % EH l = 0 special case --> TM modes
        A = 1;
        B = 0;
        C = A * (besselj(0,ha)/besselk(0,qa));
        D = 0;
    else                            % HE l = 0 special case --> TE modes
        A = 0;
        B = 1;
        C = 0;
        D = B * (besselj(0,ha)/besselk(0,qa));
    end
else
    A = 1;
    B = (1i*l*beta/(w*mu))*( qa^(-2) +  ha^(-2) )*( Jd(ha)/(ha*besselj(l,ha)) + Kd(qa)/(qa*besselk(l,qa)) )^(-1)*A;
    C = besselj(l,ha)/besselk(l,qa)*A;
    D = B*C;
end
%%%%%%%%%%%%%%%%%
% calculate the fields within the core
Ercore = @(r,pha) -((1i*beta)/h.^2)*(...            
                                A*h*Jd(h*r) + (1i*w*mu*l./(beta.*r)).*B.*besselj(l,h*r) ).*...
                                        exp(1i*(l*pha - beta*z));
Epcore = @(r,pha) -((1i*beta)/h.^2)*(...            
                                (1i*l./r).*A.*besselj(l,h*r) - (w*mu./beta).*B.*h.*Jd(h*r) ).*...
                                        exp(1i*(l*pha - beta*z));
Ezcore = @(r,pha) A*besselj(l,h*r).*exp(1i*(l*pha - beta*z));

Hrcore = @(r,pha) -((1i*beta)/h.^2)*(...            
                                B*h*Jd(h*r) - (1i*w*epsilon *l./(beta.*r)).*A.*besselj(l,h*r) ).*...
                                        exp(1i*(l*pha - beta*z));
Hpcore = @(r,pha) -((1i*beta)/h.^2)*(...            
                                (1i*l./r).*B.*besselj(l,h*r) + (w*epsilon ./beta).*A.*h.*Jd(h*r) ).*...
                                        exp(1i*(l*pha - beta*z));
Hzcore = @(r,pha) B*besselj(l,h*r).*exp(1i*(l*pha - beta*z));
% calculate the fields within the cladding
Erclad = @(r,pha) ((1i*beta)/q.^2)*(...            
                                C*q*Kd(q*r) + (1i*w*mu*l./(beta.*r)).*D.*besselk(l,q*r) ).*...
                                        exp(1i*(l*pha - beta*z));
Epclad = @(r,pha) ((1i*beta)/q.^2)*(...            
                                (1i*l./r).*C.*besselk(l,q*r) - (w*mu./beta).*D.*q.*Kd(q*r) ).*...
                                        exp(1i*(l*pha - beta*z));
Ezclad = @(r,pha) C*besselk(l,q*r).*exp(1i*(l*pha - beta*z));

Hrclad = @(r,pha) ((1i*beta)/q.^2)*(...            
                                D*q*Kd(q*r) - (1i*w*epsilon2*l./(beta.*r)).*C.*besselk(l,q*r) ).*...
                                        exp(1i*(l*pha - beta*z));
Hpclad = @(r,pha) ((1i*beta)/q.^2)*(...            
                                (1i*l./r).*D.*besselk(l,q*r) + (w*epsilon2./beta).*C.*q.*Kd(q*r) ).*...
                                        exp(1i*(l*pha - beta*z));
Hzclad = @(r,pha) D*besselk(l,q*r).*exp(1i*(l*pha - beta*z));
%%%%%%%%%%%%%%%%%
% assigning the fields within the core
region = find( rho < a & rho ~= 0 );
Er( region ) = Ercore( rho( region ), theta( region ) );
Ep( region ) = Epcore( rho( region ), theta( region ) );
Ez( region ) = Ezcore( rho( region ), theta( region ) );
Hr( region ) = Hrcore( rho( region ), theta( region ) );
Hp( region ) = Hpcore( rho( region ), theta( region ) );
Hz( region ) = Hzcore( rho( region ), theta( region ) );

% assigning the fields within the cladding
region = find( rho > a );
Er( region ) = Erclad( rho( region ), theta( region ) );
Ep( region ) = Epclad( rho( region ), theta( region ) );
Ez( region ) = Ezclad( rho( region ), theta( region ) );
Hr( region ) = Hrclad( rho( region ), theta( region ) );
Hp( region ) = Hpclad( rho( region ), theta( region ) );
Hz( region ) = Hzclad( rho( region ), theta( region ) );

% assigning the fields at the rho == 0 singular point
% interpolation by the amp. and phase of the four neighboring pixels
region = find( rho == min(rho(rho~=0)));                                    
Er( rho==0 ) = mean( abs(Er( region )) ) * exp(1i* mean( angle(Er( region )) )) ;
Ep( rho==0 ) = mean( abs(Ep( region )) ) * exp(1i* mean( angle(Ep( region )) )) ;
Ez( rho==0 ) = mean( abs(Ez( region )) ) * exp(1i* mean( angle(Ez( region )) )) ;
Hr( rho==0 ) = mean( abs(Hr( region )) ) * exp(1i* mean( angle(Hr( region )) )) ;
Hp( rho==0 ) = mean( abs(Hp( region )) ) * exp(1i* mean( angle(Hp( region )) )) ;
Hz( rho==0 ) = mean( abs(Hz( region )) ) * exp(1i* mean( angle(Hz( region )) )) ;

end

