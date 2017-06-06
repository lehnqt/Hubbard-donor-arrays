%version 2 takes into account the case when N_up or N_down is zero
function rho_jk = spdm2d(N_size, N_up, N_down, w, j1,k1,j2,k2,sigma)
rho_jk=w'*cjck2d(N_size, N_up, N_down, w, j1,k1,j2,k2,sigma);
end

