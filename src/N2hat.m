function N2h = N2hat(p,psih)
% function for the nonlinear part of equation 2

frgamma = p.sim.frgamma;
gm = p.con.gm;
Kx = p.smesh.Kx;
Ky = p.smesh.Ky;

    psi_y = ifft2(1i*Ky.*psih);
    psi_xxx = ifft2(-1i*Kx.^3.*psih);
    psi_xyy = ifft2(-1i*Kx.*Ky.^2.*psih);
    psi_x = ifft2(1i*Kx.*psih);
    psi_yxx = ifft2(-1i*Ky.*Kx.^2.*psih);
    psi_yyy = ifft2(-1i*Ky.^3.*psih);

    if frgamma == 0
        % fft of the nonlinear part
        N2h =  fft2(real(-gm *  psi_y.*(psi_xxx+psi_xyy))) ...
            + fft2(real(gm * psi_x.*(psi_yxx+psi_yyy))) ;

    elseif frgamma > 0
        %%% taking gaussian filter of rhs (seems to work)
        rhs_wo_gm = -psi_y.*(psi_xxx+psi_xyy) + psi_x.*(psi_yxx+psi_yyy);
        rhs_wo_gm_matrix = real(rhs_wo_gm);
        filtered_rhs_wo_gm = imgaussfilt(rhs_wo_gm_matrix,frgamma, ...
            'FilterDomain','spatial','Padding','circular');
        filtered_rhs_matrix = gm * filtered_rhs_wo_gm;
        N2h = fft2(real(filtered_rhs_matrix));
        % imagesc(gm*[filtered_rhs_wo_gm filtered_rhs_wo_gm]); colorbar; colormap jet;
        % drawnow;
    end
    

    %%% fft of the nonlinear part if F_gamma is applied to the N2h (this
    %%% seems to be working in a weird way)
    % N2h =  F_gamma.*( fft(real(-gm *  psi_y.*(psi_xxx+psi_xyy))) ...
    %      + fft(real(gm * psi_x.*(psi_yxx+psi_yyy))) ) ;

    %%% if F_gamma is added multiplied within the fft brackets (doesn't
    %%% work)
    % N2h =  fft( F_gamma.*real(-gm *  psi_y.*(psi_xxx+psi_xyy)) ...
    %      + F_gamma.*real(gm * psi_x.*(psi_yxx+psi_yyy))) ;

    %%% if F_gamma is added to the derivatives (doesn't work)
    % psi_y = ifft(1i*F_gamma.*qy.*psih);
    % psi_xxx = ifft((-1i)*F_gamma.*qx.^3.*psih);
    % psi_xyy = ifft((-1i)*F_gamma.*qx.*qy.^2.*psih);
    % psi_x = ifft(1i*F_gamma.*qx.*psih);
    % psi_yxx = ifft((-1i)*F_gamma.*qy.*qx.^2.*psih);
    % psi_yyy = ifft((-1i)*F_gamma.*qy.^3.*psih);
end