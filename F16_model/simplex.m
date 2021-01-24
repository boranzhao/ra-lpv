%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                  simplex.m                            %%
%%                                                       %%
%%  Author : Ying Huo                                    %%
%%                                                       %%
%% This is to set up a function minimization algorithm   %%
%% ( simplex search )to determine a steady-state trim    %%
%% condition for the 6-DOF F-16 aircraft model.          %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ s_trim, y_low ] = simplex( s, ds, x0, u0, no_step, epsilon, xcg )
 
% ---- construct the simplex 'vertices' ----
vertices = []; 
for i = 1 : 6
    for j = 1 : 7
        vertices( i, j ) = s( i );
        vertices( i, i+1 ) = s( i ) + ds( i );
    end
end
% -------------------------------------------------
%    compute the cost at the vertices of the 
%    initial simplex and sort them to get the 
%    'best' and the 'worst' vertices.
% -------------------------------------------------
y = []; % Here y represents the cost value
y0 = cost_f16( x0, u0, s, xcg ); % initial cost
y(1) = y0;
for j = 2 : 7
    y(j) = cost_f16( x0, u0, vertices( :,j ), xcg );
end
k = 7; 
y_high = y(1);
y_low = y(1);
no_high = 1;
no_low = 1;
for j = 2 : 7
    if y(j) > y_high
        y_high = y( j );
        no_high = j;
    elseif y(j) < y_low
        y_low = y( j );
        no_low = j;
    end
end
% -------------------------------------------------
%    compute the second 'worst' vertices at which
%    f = max(f(s(i))) for i ~= no_high. f here is
%    computed for later comparison with the 
%    reference point s_ref.
% -------------------------------------------------
if no_high == 1
    y_temp = y( 2 : 7 );
else
    y_temp_1 = y( 1 : no_high-1 );
    y_temp_2 = y( no_high+1 : 7 );
    y_temp = [ y_temp_1'; y_temp_2' ];
end
y_second = y_temp(1);
for j = 2 : 6
    if y_temp(j) > y_second
        y_second = y_temp( j );
    end
end

% -------------------------------------------------
%    compute the standard deviation to  
%    set up the iteration terminate criteria.
% -------------------------------------------------
y_sum = y(1);
for j = 2 : 7
    y_sum = y_sum + y(j);
end
y_mean = y_sum / 7;
deviation = 0.0;
for i = 1 : 7
    deviation = deviation + ( y(j) - y_mean )^2;
end
devia_stand = sqrt( deviation / 7 );

%-----------------------------------------------------
%  simplex search algorithm
%-----------------------------------------------------
while k <= no_step & devia_stand > epsilon 
    s_centroid = [];
    s_ref = [];
    %---- compute the centroid point ----
    for i = 1 : 6
        s_centroid(i) = 0.0;
        for j = 1 : 7
            if j ~= no_high
                s_centroid( i ) = s_centroid( i ) + vertices( i, j ); 
            end
        end
    s_centroid(i) = s_centroid(i) / 6;    
    end
    %---- compute the reflection point ----
    for i = 1 : 6
        s_ref(i) = 2 * s_centroid( i )- vertices( i, no_high ); 
    end
    y_ref = cost_f16( x0, u0, s_ref, xcg );   
    
    %---- compute the new point for the simplex ------------------
    if y_ref < y_low         %---- case I   ----
        s_exp = 2 * s_ref - s_centroid;
        y_exp = cost_f16( x0, u0, s_exp, xcg );
        if y_exp < y_ref
            s_new = s_exp;
        else 
            s_new = s_ref;
        end
    elseif y_ref < y_second    %---- case II  ----
        s_new = s_ref;
    else                       %---- case III ----
        if y_ref >= y_high   
            s_new = 0.5 * ( vertices( :, no_high )' + s_centroid );
        else
            s_new = 0.5 * ( s_ref + s_centroid );
        end
    end
    %---- form the new simplex with the max point no_high
    %     replaced by the new point derived above ------------------
    vertices( :, no_high ) = s_new';
    %---- upgrade step k -----------
    k = k + 6;
    % -------------------------------------------------
    %    compute the cost at the vertices of the 
    %    new simplex and sort them to get the 
    %    'best' and the 'worst' vertices.
    % -------------------------------------------------
    y = []; % Here y represents the cost value
    y0 = cost_f16( x0, u0, vertices( :, 1 ), xcg ); % new cost
    y(1) = y0;
    for j = 2 : 7
        y(j) = cost_f16( x0, u0, vertices( :,j ), xcg );
    end
    y_high = y(1);
    y_low = y(1);
    no_high = 1;
    no_low = 1;
    for j = 2 : 7
        if y(j) > y_high
            y_high = y( j );
            no_high = j;
        elseif y(j) < y_low
            y_low = y( j );
            no_low = j;
        end
    end
    % -------------------------------------------------
    %    compute the second 'worst' vertices 
    % -------------------------------------------------
    if no_high == 1
        y_temp = y( 2 : 7 );
    else
        y_temp_1 = y( 1 : no_high-1 );
        y_temp_2 = y( no_high+1 : 7 );
        y_temp = [ y_temp_1'; y_temp_2' ];
    end
    y_second = y_temp(1);
    for j = 2 : 6
        if y_temp(j) > y_second
            y_second = y_temp( j );
        end
    end
    %---- calculate the standard deviation again ----
    y_sum = y(1);
    for j = 2 : 7
        y_sum = y_sum + y(j);
    end
    y_mean = y_sum / 7;
    deviation = 0.0;
    for i = 1 : 7
        deviation = deviation + ( y(j) - y_mean )^2;
    end
    devia_stand = sqrt( deviation / 7 );
end

epsilon = devia_stand;
no_step = k;
% terminate_text = [ ' The iteration is terminated successfully at step = ',...
%             num2str(k-1), ' and the standard deviation epsilon = ' num2str(epsilon) ];
% display( terminate_text );
y_low = y( no_low );
s_trim = vertices( :, no_low );
