function rho = rho_fcn(alt)
    r0 = 2.37764e-3;  % const, sea-level density, slugs/ft^3
    tfac = 1.0 - 0.703e-5 * alt;
    rho = r0.* ( tfac.^ 4.14 );  % density
end