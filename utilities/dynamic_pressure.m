function [pbar_s, pbar] = dynamic_pressure(alti,veloc,pbar_lower_bnd,pbar_upper_bnd)

pbar = 0.5*rho_fcn(alti)*veloc^2;   % dynamic pressure

pbar_s = (pbar-pbar_lower_bnd)/(pbar_upper_bnd-pbar_lower_bnd)*2-1;  % scaled dynamic pressure, [-1,1]