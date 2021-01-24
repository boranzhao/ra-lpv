%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                   fdu.m                          %%
%%                                                  %%
%%  Author : Ying Huo                               %%
%%                                                  %%
%% This is to approximate the partial derivative    %%
%% of f w.r.t. control by computing the difference  %% 
%% of the state derivative, i.e. f(state, control)  %%
%% w.r.t.control at (control+du) and (control-du).  %% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function  dfdu = fdu( x, control, xcg, i, j, du )

time = 0.0;
t = control(j);
control(j) = t - du;  
[ x_dot, an, alat, qbar, amach, q, alpha ] = f16_dynam ( time, x, control, xcg );
xd1 = x_dot(i);       
control(j) = t + du;  
[ x_dot, an, alat, qbar, amach, q, alpha ] = f16_dynam ( time, x, control, xcg );
xd2 = x_dot(i);       
dfdu = ( xd2 - xd1 ) / ( du + du ); 
control(j) = t;       % return to original control(j)

