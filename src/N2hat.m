function N2h = N2hat(gm,psih,qx,qy)
% N2h = N2hat(gm,psih,qx,qy)
% yields the fft of the nonlinear part of gsh equation 2
% for the pseudospectral solve

    % derivatives in the fourier domain
    psi_y = ifft(1i*qy.*psih);
    psi_xxx = ifft(-1i*qx.^3.*psih);
    psi_xyy = ifft(-1i*qx.*qy.^2.*psih);
    psi_x = ifft(1i*qx.*psih);
    psi_yxx = ifft(-1i*qy.*qx.^2.*psih);
    psi_yyy = ifft(-1i*qy.^3.*psih);
    
    % fft of the nonlinear part
    N2h = fft(real(-gm * psi_y.*(psi_xxx+psi_xyy))) ...
         + fft(real(gm * psi_x.*(psi_yxx+psi_yyy)));
end