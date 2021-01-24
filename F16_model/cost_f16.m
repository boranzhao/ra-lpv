%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%             cost_f16.m                     %%
%%                                            %%
%%  Author : Ying Huo                         %%
%%                                            %%
%% This is to define & compute the cost of    %%
%% the current simplex where the states are   %%
%% constrained by different flight conditions %%
%% as defined in constraint.m                 %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function cost = cost_f16( state, control, s, xcg )

% state   - states before constrained
% control - current control
%  s      - current simplex

control = [ s(1); s(2); s(4); s(5) ];
state(2) = s( 3 );
state(3) = s( 6 );
state(13) = tgear( control(1) );
%---- constraints -------------------------------
state_constr = constraint_f16( state );
time = 0;
% global xcg;
[ xd, an, alat, qbar, amach, q, alpha ] = f16_dynam ( time, state_constr, control, xcg );
% original 
% cost = xd(1)^2 + 100 * ( xd(2)^2 + xd(3)^2 ) + ...
%                  10 * ( xd(7)^2 + xd(8)^2 + xd(9)^2 );
             
             
cost = xd(1)^2 + 100 * ( xd(2)^2 + xd(3)^2 + xd(4)^2 + xd(5)^2+xd(6)^2)+... 
                 10 * ( xd(7)^2 + xd(8)^2 + xd(9)^2 ) +xd(11)^2;

