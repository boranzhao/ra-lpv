%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                  CM.m                       %%
%%                                             %%
%%  Author : Ying Huo                          %%
%%                                             %%
%%  pitching moment coefficients in F-16 model %% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function cm_value = cm( alpha, el )

A = [ .205   .168   .186   .196   .213   .251   .245   .238   .252   .231   .198   .192;
      .081   .077   .107   .110   .110   .141   .127   .119   .133   .108   .081   .093;
     -.046  -.020  -.009  -.005  -.006   .010   .006  -.001   .014   .000  -.013   .032;
     -.174  -.145  -.121  -.127  -.129  -.102  -.097  -.113  -.087  -.084  -.069  -.006;
     -.259  -.202  -.184  -.193  -.199  -.150  -.160  -.167  -.104  -.076  -.041  -.005 ];
A = A';  
row = 3;
col = 3;

s = 0.2 * alpha;
k = fix ( s );
if ( k <= -2 )
    k = -1;
end
if ( k >= 9 )
    k = 8;
end
da = s - k;
l = k + fix ( sign ( da ) * 1.1 );
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% added as boundary condition
if l < -2
    l = -2;
elseif l > 9
        l = 9;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
s = el / 12.0;
m = fix ( s );
if ( m <= -2 )
   m = -1;
end
if ( m >= 2)
   m = 1;
end
de = s - m;
n = m + fix ( sign ( de ) * 1.1 );
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% added as boundary condition
if n < -2
    n = -2;
elseif n > 2
        n = 2;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t = A( k+row, m+col );
u = A( k+row, n+col );
v = t + abs( da ) * ( A( l+row, m+col ) - t );
w = u + abs( da ) * ( A( l+row, n+col ) - u );
cm_value = v + ( w - v ) * abs( de );
