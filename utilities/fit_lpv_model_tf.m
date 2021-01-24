% version history
% 01/03/2020: Change the interpolation method into 
%  interpolating the poles and zeros of linearized models instead of the
%  coefficients 

clear;
clc;
% close all;



scheduleParaNum = 2;

syms pbar_s veloc_s % scaled dynamic pressure, [-1, 1]

% Fit the Bm matrix
load('lin_results.mat','pbars','As','Bs','Cs','Ds','altitudes','velocities','state_trims','control_trims','xcg')

[m,n] = size(altitudes);

% sort the dynamic pressure
pbar_vec = reshape(pbars,m*n,1);
[pbar_vec,Index] = sort(pbar_vec);

veloc_vec = reshape(velocities,m*n,1);
% veloc_vec = veloc_vec(Index)/1000; % sort and nomalize
veloc_min = min(veloc_vec);
veloc_max = max(veloc_vec);

veloc_svec = (veloc_vec-veloc_min)/(veloc_max-veloc_min)*2-1;


% state_trims = reshape(state_trims,size(state_trims,1),m*n);
% state_trims = state_trims(:,I);
% 
% control_trims = reshape(control_trims,size(control_trims,1),m*n);
% control_trims = control_trims(:,I);


pbar_lower_bnd = pbar_vec(1);
pbar_upper_bnd = pbar_vec(end);
% test of scaling
pbar_svec = (pbar_vec-pbar_lower_bnd)/(pbar_upper_bnd-pbar_lower_bnd)*2-1;
% pbar_vec2 = (pbar_s_vec+1)*(pbar_upper_bnd-pbar_lower_bnd)/2+pbar_lower_bnd
% figure(2);
% plot(pbar_vec,pbar_vec2);

s = tf('s');
tf1_nums = zeros(2,m,n);
tf2_nums = zeros(2,m,n);
tf1_dens = zeros(3,m,n);

for i=1:m
    for j=1:n       
        tf1 = tf(ss(As(:,:,i,j),Bs(:,:,i,j),[1 0],0));
        tf2 = tf(ss(As(:,:,i,j),Bs(:,:,i,j),[0 1],0));
        
        [tf1_num, tf1_den] = tfdata(tf1);
        [tf2_num, tf2_den] = tfdata(tf2);
        if norm(tf1_den{1}-tf2_den{1})>1e-10
            error('The denominators of the two transfer functions should be the same');
        end
        if norm([tf1_num{1}(1) tf2_num{1}(1)])>1e-10 
            error('The order of at least one numerator is large than 1');
        end
        tf1_nums(:,i,j) = (tf1_num{1}([2 3]))';
        tf2_nums(:,i,j) = (tf2_num{1}([2 3]))';
        tf1_dens(:,i,j) = (tf1_den{1})';
    end
end


% Fit the tf1_num and tf2_num vectors
tf1_num1_vec = reshape(tf1_nums(1,:,:),m*n,1); % get the first element of tf1_num
tf1_num2_vec = reshape(tf1_nums(2,:,:),m*n,1); % get the second element of tf1_num
tf2_num1_vec = reshape(tf2_nums(1,:,:),m*n,1); % get the first element of tf1_num
tf2_num2_vec = reshape(tf2_nums(2,:,:),m*n,1); % get the second element of tf1_num

% Fit the tf_den vectors (not that the first element is always 1
tf1_den1_vec = reshape(tf1_dens(2,:,:),m*n,1);
tf1_den2_vec = reshape(tf1_dens(3,:,:),m*n,1);


tf1_num1_vec = tf1_num1_vec(Index);
tf1_num2_vec = tf1_num2_vec(Index);
tf2_num1_vec = tf2_num1_vec(Index);
tf2_num2_vec = tf2_num2_vec(Index);

tf1_den1_vec = tf1_den1_vec(Index);
tf1_den2_vec = tf1_den2_vec(Index);


order = 1;

if scheduleParaNum == 1
% use only the pbar_s as the independent variable
    p_coef_tf1_num1 =  polyfit(pbar_svec,tf1_num1_vec,order);
    p_coef_tf1_num2 =  polyfit(pbar_svec,tf1_num2_vec,order);
    p_coef_tf2_num1 =  polyfit(pbar_svec,tf2_num1_vec,order);
    p_coef_tf2_num2 =  polyfit(pbar_svec,tf2_num2_vec,order);

    tf1_num1_fit = polyval(p_coef_tf1_num1,pbar_svec);
    tf1_num2_fit = polyval(p_coef_tf1_num2,pbar_svec);
    tf2_num1_fit = polyval(p_coef_tf2_num1,pbar_svec);
    tf2_num2_fit = polyval(p_coef_tf2_num2,pbar_svec);

    % fit a first order polynomial curve
    p_coef_tf1_den1 =  polyfit(pbar_svec,tf1_den1_vec,order);
    p_coef_tf1_den2 =  polyfit(pbar_svec,tf1_den2_vec,order);

    tf1_den1_fit = polyval(p_coef_tf1_den1,pbar_svec);
    tf1_den2_fit = polyval(p_coef_tf1_den2,pbar_svec);

    figure(1);
    title('Fit for numerator');
    subplot(2,2,1);
    plot(pbar_svec,tf1_num1_vec,'o',pbar_svec,tf1_num1_fit,'-'); hold on;
    legend('data','linear fit');
    xlabel('dynamic pressure');
    ylabel('tf1_num1');
    subplot(2,2,2);
    plot(pbar_svec,tf1_num2_vec,'o',pbar_svec,tf1_num2_fit,'-'); hold on;
    legend('data','linear fit');
    xlabel('dynamic pressure');
    ylabel('tf1_num2');
    subplot(2,2,3);
    plot(pbar_svec,tf2_num1_vec,'o',pbar_svec,tf2_num1_fit,'-'); hold on;
    legend('data','linear fit');
    xlabel('dynamic pressure');
    ylabel('tf2_num1');
    subplot(2,2,4);
    plot(pbar_svec,tf2_num2_vec,'o',pbar_svec,tf2_num2_fit,'-'); hold on;
    legend('data','linear fit');
    xlabel('dynamic pressure');
    ylabel('tf2_num2');

    figure(2);
    title("Fit for denominator");
    subplot(2,1,1);
    plot(pbar_svec,tf1_den1_vec,'o',pbar_svec,tf1_den1_fit,'-'); hold on;
    legend('data','linear fit');
    xlabel('dynamic pressure');
    ylabel('tf1_den1');
    subplot(2,1,2);
    plot(pbar_svec,tf1_den2_vec,'o',pbar_svec,tf1_den2_fit,'-'); hold on;
    legend('data','linear fit');
    xlabel('dynamic pressure');
    ylabel('tf1_den2');
    
    tf1_num1_sym = p_coef_tf1_num1(1)*pbar_s+p_coef_tf1_num1(2);
    tf1_num2_sym = p_coef_tf1_num2(1)*pbar_s+p_coef_tf1_num2(2);
    tf2_num1_sym = p_coef_tf2_num1(1)*pbar_s+p_coef_tf2_num1(2);
    tf2_num2_sym = p_coef_tf2_num2(1)*pbar_s+p_coef_tf2_num2(2);
    
    tf1_den1_sym = p_coef_tf1_den1(1)*pbar_s+p_coef_tf1_den1(2);
    tf1_den2_sym = p_coef_tf1_den2(1)*pbar_s+p_coef_tf1_den2(2);
    
    
    T = [tf1_num2_sym tf1_num1_sym; tf2_num1_sym tf2_num1_sym];
    Tinv = inv(T);
    T = simplify(T);
    Tinv = simplify(Tinv);
    
     
    Ap = simplify(T*[0 1;-tf1_den2_sym -tf1_den1_sym]*Tinv);
    Bp = T*[0;1];
    
    % extract the constant and linear parts
    Ap0 = double(subs(Ap,0));
    Ap1 = double(subs(Ap,1))-Ap0;
    Bp0 = double(subs(Bp,0));
    Bp1 = double(subs(Bp,1))-Bp0;
    
    Ap2 = zeros(size(Ap0));
    Bp2 = zeros(size(Bp0));    

elseif scheduleParaNum == 2

    % use both pbar_s and velocities
    p_coef_tf1_num1 =  polyfitn([pbar_svec,veloc_svec],tf1_num1_vec,order);
    p_coef_tf1_num2 =  polyfitn([pbar_svec,veloc_svec],tf1_num2_vec,order);
    p_coef_tf2_num1 =  polyfitn([pbar_svec,veloc_svec],tf2_num1_vec,order);
    p_coef_tf2_num2 =  polyfitn([pbar_svec,veloc_svec],tf2_num2_vec,order);

    tf1_num1_fit = polyvaln(p_coef_tf1_num1,[pbar_svec,veloc_svec]);
    tf1_num2_fit = polyvaln(p_coef_tf1_num2,[pbar_svec,veloc_svec]);
    tf2_num1_fit = polyvaln(p_coef_tf2_num1,[pbar_svec,veloc_svec]);
    tf2_num2_fit = polyvaln(p_coef_tf2_num2,[pbar_svec,veloc_svec]);

    % fit a first order polynomial curve
    p_coef_tf1_den1 =  polyfitn([pbar_svec,veloc_svec],tf1_den1_vec,order);
    p_coef_tf1_den2 =  polyfitn([pbar_svec,veloc_svec],tf1_den2_vec,order);

    tf1_den1_fit = polyvaln(p_coef_tf1_den1,[pbar_svec,veloc_svec]);
    tf1_den2_fit = polyvaln(p_coef_tf1_den2,[pbar_svec,veloc_svec]);
    

    figure(1);
    title('Fit for numererators');
    subplot(2,2,1);
    plot3(pbar_svec,veloc_svec, tf1_num1_vec,'o',pbar_svec,veloc_svec,tf1_num1_fit,'-'); hold on;
    legend('data','linear fit');
    xlabel('dynamic pressure');
    ylabel('velocity');    
    zlabel('tf1_num1');
    subplot(2,2,2);
    plot3(pbar_svec,veloc_svec,tf1_num2_vec,'o',pbar_svec,veloc_svec,tf1_num2_fit,'-'); hold on;
    legend('data','linear fit');
    xlabel('dynamic pressure');
    ylabel('velocity');
    zlabel('tf1_num2');
    subplot(2,2,3);
    plot3(pbar_svec,veloc_svec,tf2_num1_vec,'o',pbar_svec,veloc_svec,tf2_num1_fit,'-'); hold on;
    legend('data','linear fit');
    xlabel('dynamic pressure');
    ylabel('velocity');
    zlabel('tf2_num1');
    subplot(2,2,4);
    plot3(pbar_svec,veloc_svec,tf2_num2_vec,'o',pbar_svec,veloc_svec,tf2_num2_fit,'-'); hold on;
    legend('data','linear fit');
    xlabel('dynamic pressure');
    ylabel('velocity');
    zlabel('tf2_num2');
    
    figure(2);
    title("Fit for denominator");
    subplot(2,1,1);
    plot3(pbar_svec,veloc_svec,tf1_den1_vec,'o',pbar_svec,veloc_svec,tf1_den1_fit,'-'); hold on;
    legend('data','linear fit');
    xlabel('dynamic pressure');
    ylabel('velocity');
    zlabel('tf1_den1');
    subplot(2,1,2);
    plot3(pbar_svec,veloc_svec,tf1_den2_vec,'o',pbar_svec,veloc_svec,tf1_den2_fit,'-'); hold on;
    legend('data','linear fit');
    xlabel('dynamic pressure');
    ylabel('velocity');
    zlabel('tf1_den2');
    
    % generate the parameter-dependent matrices
    tf1_num1_sym = p_coef_tf1_num1.Coefficients(1)*pbar_s+p_coef_tf1_num1.Coefficients(2)*veloc_s+p_coef_tf1_num1.Coefficients(3);
    tf1_num2_sym = p_coef_tf1_num2.Coefficients(1)*pbar_s+p_coef_tf1_num2.Coefficients(2)*veloc_s+p_coef_tf1_num2.Coefficients(3);
    tf2_num1_sym = p_coef_tf2_num1.Coefficients(1)*pbar_s+p_coef_tf2_num1.Coefficients(2)*veloc_s+p_coef_tf2_num1.Coefficients(3);
    tf2_num2_sym = p_coef_tf2_num2.Coefficients(1)*pbar_s+p_coef_tf2_num2.Coefficients(2)*veloc_s+p_coef_tf2_num2.Coefficients(3);
    
    tf1_den1_sym = p_coef_tf1_den1.Coefficients(1)*pbar_s+p_coef_tf1_den1.Coefficients(2)*veloc_s+p_coef_tf1_den1.Coefficients(3);
    tf1_den2_sym = p_coef_tf1_den2.Coefficients(1)*pbar_s+p_coef_tf1_den2.Coefficients(2)*veloc_s+p_coef_tf1_den2.Coefficients(3);
    
    
    T = [tf1_num2_sym tf1_num1_sym; tf2_num1_sym tf2_num1_sym];
    Tinv = inv(T);
    T = simplify(T);
    Tinv = simplify(Tinv);
    
     
    Ap = simplify(T*[0 1;-tf1_den2_sym -tf1_den1_sym]*Tinv);
    Bp = T*[0;1];
    
    % extract the constant and linear parts
    Ap0 = double(subs(Ap,{pbar_s,veloc_s},{0,0}));
    Ap1 = double(subs(Ap,{pbar_s,veloc_s},{1,0}))-Ap0;
    Ap2 = double(subs(Ap,{pbar_s,veloc_s},{0,1}))-Ap0;
    
    Bp0 = double(subs(Bp,{pbar_s,veloc_s},{0,0}));
    Bp1 = double(subs(Bp,{pbar_s,veloc_s},{1,0}))-Bp0;
    Bp2 = double(subs(Bp,{pbar_s,veloc_s},{0,1}))-Bp0;
end

return;
%% verify the LPV model at a numer of points

alti_0 = 10000;
veloc_0 = 900; 

idx_row = find(velocities(:,1)==veloc_0,1);
idx_col = find(altitudes(1,:)==alti_0,1);

Alti = As(:,:,idx_row,idx_col);
Blti = Bs(:,:,idx_row,idx_col);

[pbar_s,~] = dynamic_pressure(alti_0,veloc_0,pbar_lower_bnd,pbar_upper_bnd);

veloc_s = veloc_0/900;

Alpv = double(subs(Ap));
Blpv = double(subs(Bp));



figure; bode(ss(Alti,Blti,eye(2),0),ss(Alpv,Blpv,eye(2),0));
legend('lti','lpv');


% fprintf(1,'A11 error norm: %.3f,omega: %.3f pole bounds: [%.3f  %.3f],actual poles: %.3f, %.3f\n',theta1, abs(poles(1)), -0.7*(theta1+4+offset_upp), -0.7*(theta1+4+offset_low),poles); 


%  desired dynamics
% note that pbar_s =
% (pbar-pbar_lower_bnd)/(pbar_upper_bnd-pbar_lower_bnd)*2-1

 
% pbar = (pbar_s+1)*(pbar_upper_bnd-pbar_lower_bnd)/2+pbar_lower_bnd;
% Am = [0 1; -omega_p^2 2*zeta*omega_p];
