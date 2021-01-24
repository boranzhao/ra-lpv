%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                adc.m                       %%
%%                                            %%
%%  Author : Ying Huo                         %%
%%                                            %%
%% standard atmosphere model used in F-16     %%
%% model  ( air data computer )               %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [tfac, t, rho, amach, qbar, ps ] = adc ( vt, alt )

r0 = 2.37764e-3;  % const, sea-level density, slugs/ft^3

tfac = 1.0 - 0.703e-5 * alt;
t = 519.0 * tfac;  % temperature, Rankine scale
if ( alt >= 35000.0 )
    t = 390.0;
end
rho = r0 * ( tfac ^ 4.14 );  % density
amach = vt / sqrt ( 1.4 * 1716.3 * t ); % Mach number
qbar = 0.5 * rho * vt * vt;  % dynamic pressure, psf
ps = 1715.0 * rho * t;  % static pressure


