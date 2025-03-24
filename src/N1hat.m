function N1h = N1hat(psih,zetah,qx,qy)
% N1h = N1hat(psih,zetah,qx,qy)
% yields the fft of the nonlinear part of gsh equation 1
% for the pseudospectral solve
    N1h = -fft(real((ifft(psih)).^3)) ...
            - fft(real( ifft(1i*qy.*zetah) .* ifft(1i*qx.*psih) )) ...
                + fft(real( ifft(1i*qx.*zetah) .* ifft(1i*qy.*psih) ));
end