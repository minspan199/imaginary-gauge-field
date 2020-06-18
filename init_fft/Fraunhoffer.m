function [U_det] = Fraunhoffer(U_ap, k, z, xs, ys, lambda)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fraunhoffer
% Translates Near Field (U_ap) to Far Field (U) 
% (assuming various Fraunhoffer assumptions are valid)
% INPUTS:
% U_ap          (Field at aperture)
% 
% OUTPUTS:
% U             (Field at detector)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Angular_Spectrum = fftshift(fft2(fftshift(U_ap))); %Near Field
U_det = exp(i*k*z)*exp(i*k*(xs.^2+ys.^2)/(2*z)).*Angular_Spectrum...
    /(i*lambda*z);
end

