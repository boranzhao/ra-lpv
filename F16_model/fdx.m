%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                   fdx.m                          %%
%%                                                  %%
%%  Author : Ying Huo                               %%
%%                                                  %%
%% This is to approximate the partial derivative    %%
%% of f w.r.t. state by computing the difference    %% 
%% of the state derivative, i.e. f(state, control)  %%
%% w.r.t. state at (state+dx) and (state-dx).       %% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function  dfdx = fdx( x, control, xcg, i, j, dx )

time = 0.0;
t = x(j);
x(j) = t - dx;  
[ x_dot, an, alat, qbar, amach, q, alpha ] = f16_dynam ( time, x, control, xcg );
xd1 = x_dot(i); 
x(j) = t + dx;  
[ x_dot, an, alat, qbar, amach, q, alpha ] = f16_dynam ( time, x, control, xcg );
xd2 = x_dot(i); 
dfdx = ( xd2 - xd1 ) / ( dx + dx ); 
x(j) = t;       % return to original x

