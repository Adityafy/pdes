function N1h = N1hat(p,psihmat,zetahmat)
% function for the nonlinear part of equation 1
    Kx = p.smesh.Kx;
    Ky = p.smesh.Ky;
    N1h = -fft2(real((ifft2(psihmat)).^3)) ...
            - fft2(real( ifft2(1i*Ky.*zetahmat) .* ifft2(1i*Kx.*psihmat) )) ...
                + fft2(real( ifft2(1i*Kx.*zetahmat) .* ifft2(1i*Ky.*psihmat) ));
end