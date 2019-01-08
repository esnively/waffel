%%%% Initialize frequency components of THz field %%%%

param.w0 = 2*pi*c/param.lambda0;                                                 % angular frequency in free space

field_refine_factor = 20;
if(param.nfreq == 1)
    omega = param.w0;
    delta_omega = omega;
    lengthfactor =  param.lambda0/c;
    cq0 = 1;
    deltat_nyquist = c*pi*2/omega/field_refine_factor; % made smaller for genesis case, not sampling enough points
    z0 = -4*param.lambda0:deltat_nyquist:4*param.lambda0;
else 
    omega = linspace(param.w0*(1-param.deltaBW),param.w0*(1+param.deltaBW),param.nfreq); %param.w0*(1-param.deltaBW+2*param.deltaBW*(0:param.nfreq-1)/(param.nfreq-1));
    delta_omega = diff(omega(1:2));
    lengthfactor = (param.nfreq-1)*param.w0/range(omega)*param.lambda0/c; %(param.nfreq-1)/2/param.deltaBW*param.lambda0/c;
    cq0 = exp(-(omega-param.w0).^2/2/param.sw^2);
    
    deltat_nyquist = c*pi*2/(max(omega)-min(omega)); 
    deltat_tot = 2*pi*c/delta_omega;
  
    z0 = -deltat_tot/2:deltat_nyquist/field_refine_factor:deltat_tot/2;
end

k_01 = (pi+atan(param.b/sqrt(2.*param.Rppwg*param.b-param.b^2)))/param.b;   % PPWG transverse wavenumber
k_zq = sqrt(omega.^2/c^2 - k_01^2);                                            % dispersion relation

k_zq(omega.^2 <= c^2*k_01^2) =[];
cq0(omega.^2 <= c^2*k_01^2) =[]; 
omega(omega.^2 <= c^2*k_01^2) =[];
param.nfreq = length(omega);
k_zq(imag(k_zq)>0) = -k_zq(imag(k_zq)>0);

kKl_temp = delta_omega/2/pi*cq0*exp(-1i*omega(:)*z0/c);
[peakvalue,peakposition] = max(abs(real(kKl_temp)));
cq0 = cq0.*param.E0/peakvalue.*exp(-1i*omega(:)*z0(peakposition)/c).';
kKl0 = delta_omega/2/pi*cq0*exp(-1i*omega(:)*z0/c);

param.zoffset = param.timingTHz*c;
cq0 = cq0.*exp(1i*k_zq(:).'*param.zoffset);
cqprime = zeros(param.nfreq,1);

param.k_z = k_zq;
param.omega = omega;
param.delta_omega = delta_omega;
param.z0 = z0;
param.k_z0 = sqrt(param.w0^2/c^2 - k_01^2); 

%% effective mode area TE_01 (based on fit to mathematica mode model in "CPPWG mode coupling.nb")
A_mode = pi*param.waist^2/2; %original approximation
%Ar = param.b/2;
%A_mode = 6.37927e-7-0.00106325*Ar+0.32453*Ar^2+338.599*Ar^3-70241.9*Ar^4;

ETHz = A_mode*sqrt(eps0/mu0)*(omega(2)-omega(1))/2/pi*sum(abs(cq0).^2);
