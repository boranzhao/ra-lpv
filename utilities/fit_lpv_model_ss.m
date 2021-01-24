%% Fit a state-space LPV model

% clear;
clc;
close all;



scheduleParaNum = 2;

syms pbar_s veloc_s % scaled dynamic pressure, [-1, 1]

% Fit the Bm matrix
load('lin_results.mat','pbars','As','Bs','Cs','Ds','altitudes','velocities','state_trims','control_trims','xcg')

[m,n] = size(altitudes);

% sort the dynamic pressure
pbar_vec = reshape(pbars,m*n,1);
[pbar_vec,Index] = sort(pbar_vec);

veloc_vec = reshape(velocities,m*n,1);

veloc_min = min(veloc_vec);
veloc_max = max(veloc_vec);

veloc_svec = (veloc_vec-veloc_min)/(veloc_max-veloc_min)*2-1;
% veloc_vec = veloc_vec(Index)/1000; % sort and nomalize

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

% Fit the B matrix
Bs1_vec = reshape(Bs(1,1,:,:),m*n,1); % get the first element 
Bs2_vec = reshape(Bs(2,1,:,:),m*n,1); % get the second element 
Bs1_vec = Bs1_vec(Index);
Bs2_vec = Bs2_vec(Index);

% Fit the As matrix
A11_vec = reshape(As(1,1,:,:),m*n,1); % get the second element 
A12_vec = reshape(As(1,2,:,:),m*n,1); % get the second element 
A21_vec = reshape(As(2,1,:,:),m*n,1); % get the second element 
A22_vec = reshape(As(2,2,:,:),m*n,1); % get the second element 

A11_vec = A11_vec(Index);
A12_vec = A12_vec(Index);
A21_vec = A21_vec(Index);
A22_vec = A22_vec(Index);

order = 1;

if scheduleParaNum == 1
% use only the pbar_s as the independent variable
    p_coef_A11 =  polyfit(pbar_svec,A11_vec,order);
    p_coef_A12 =  polyfit(pbar_svec,A12_vec,order);
    p_coef_A21 =  polyfit(pbar_svec,A21_vec,order);
    p_coef_A22 =  polyfit(pbar_svec,A22_vec,order);

    A11_fit = polyval(p_coef_A11,pbar_svec);
    A12_fit = polyval(p_coef_A12,pbar_svec);
    A21_fit = polyval(p_coef_A21,pbar_svec);
    A22_fit = polyval(p_coef_A22,pbar_svec);

    % fit a first order polynomial curve
    p_coef_B1 =  polyfit(pbar_svec,Bs1_vec,order);
    p_coef_B2 =  polyfit(pbar_svec,Bs2_vec,order);

    Bs1_fit = polyval(p_coef_B1,pbar_svec);
    Bs2_fit = polyval(p_coef_B2,pbar_svec);

    figure(1);
    title('Fit for A matrix');
    subplot(2,2,1);
    plot(pbar_svec,A11_vec,'o',pbar_svec,A11_fit,'-'); hold on;
    legend('data','linear fit');
    xlabel('dynamic pressure');
    ylabel('A11');
    subplot(2,2,2);
    plot(pbar_svec,A12_vec,'o',pbar_svec,A12_fit,'-'); hold on;
    legend('data','linear fit');
    xlabel('dynamic pressure');
    ylabel('A12');
    subplot(2,2,3);
    plot(pbar_svec,A21_vec,'o',pbar_svec,A21_fit,'-'); hold on;
    legend('data','linear fit');
    xlabel('dynamic pressure');
    ylabel('A21');
    subplot(2,2,4);
    plot(pbar_svec,A22_vec,'o',pbar_svec,A22_fit,'-'); hold on;
    legend('data','linear fit');
    xlabel('dynamic pressure');
    ylabel('A22');

    figure(2);
    title("Fit for B matrix");
    subplot(2,1,1);
    plot(pbar_svec,Bs1_vec,'o',pbar_svec,Bs1_fit,'-'); hold on;
    legend('data','linear fit');
    xlabel('dynamic pressure');
    ylabel('B1');
    subplot(2,1,2);
    plot(pbar_svec,Bs2_vec,'o',pbar_svec,Bs2_fit,'-'); hold on;
    legend('data','linear fit');
    xlabel('dynamic pressure');
    ylabel('B2');
    
    
    % generate the parameter-dependent matrices
    Ap = [p_coef_A11(1)*pbar_s+p_coef_A11(2) 0.94; p_coef_A21(1)*pbar_s+p_coef_A21(2) p_coef_A22(1)*pbar_s+p_coef_A22(2)];
    Bp = [p_coef_B1(1)*pbar_s+p_coef_B1(2) p_coef_B2(1)*pbar_s+p_coef_B2(2)]';
    Cp = [1 0];

    % extract the constant and linear parts
    Ap0 = double(subs(Ap,0));
    Ap1 = double(subs(Ap,1))-Ap0;
    Bp0 = double(subs(Bp,0));
    Bp1 = double(subs(Bp,1))-Bp0;
    
    Ap2 = zeros(size(Ap0));
    Bp2 = zeros(size(Bp0));    

elseif scheduleParaNum == 2

    % use both pbar_s and velocities
    p_coef_A11 =  polyfitn([pbar_svec,veloc_svec],A11_vec,order);
    p_coef_A12 =  polyfitn([pbar_svec,veloc_svec],A12_vec,order);
    p_coef_A21 =  polyfitn([pbar_svec,veloc_svec],A21_vec,order);
    p_coef_A22 =  polyfitn([pbar_svec,veloc_svec],A22_vec,order);

    A11_fit = polyvaln(p_coef_A11,[pbar_svec,veloc_svec]);
    A12_fit = polyvaln(p_coef_A12,[pbar_svec,veloc_svec]);
    A21_fit = polyvaln(p_coef_A21,[pbar_svec,veloc_svec]);
    A22_fit = polyvaln(p_coef_A22,[pbar_svec,veloc_svec]);

    % fit a first order polynomial curve
    p_coef_B1 =  polyfitn([pbar_svec,veloc_svec],Bs1_vec,order);
    p_coef_B2 =  polyfitn([pbar_svec,veloc_svec],Bs2_vec,order);

    Bs1_fit = polyvaln(p_coef_B1,[pbar_svec,veloc_svec]);
    Bs2_fit = polyvaln(p_coef_B2,[pbar_svec,veloc_svec]);

    figure(1);
    title('Fit for A matrix');
    subplot(2,2,1);
    plot3(pbar_svec,veloc_svec, A11_vec,'o',pbar_svec,veloc_svec,A11_fit,'-'); hold on;
    legend('data','linear fit');
    xlabel('dynamic pressure');
    ylabel('velocity');    
    zlabel('A11');
    subplot(2,2,2);
    plot3(pbar_svec,veloc_svec,A12_vec,'o',pbar_svec,veloc_svec,A12_fit,'-'); hold on;
    legend('data','linear fit');
    xlabel('dynamic pressure');
    ylabel('velocity');
    zlabel('A12');
    subplot(2,2,3);
    plot3(pbar_svec,veloc_svec,A21_vec,'o',pbar_svec,veloc_svec,A21_fit,'-'); hold on;
    legend('data','linear fit');
    xlabel('dynamic pressure');
    ylabel('velocity');
    zlabel('A21');
    subplot(2,2,4);
    plot3(pbar_svec,veloc_svec,A22_vec,'o',pbar_svec,veloc_svec,A22_fit,'-'); hold on;
    legend('data','linear fit');
    xlabel('dynamic pressure');
    ylabel('velocity');
    zlabel('A22');
    
    figure(2);
    title("Fit for B matrix");
    subplot(2,1,1);
    plot3(pbar_svec,veloc_svec,Bs1_vec,'o',pbar_svec,veloc_svec,Bs1_fit,'-'); hold on;
    legend('data','linear fit');
    xlabel('dynamic pressure');
    ylabel('velocity');
    zlabel('B1');
    subplot(2,1,2);
    plot3(pbar_svec,veloc_svec,Bs2_vec,'o',pbar_svec,veloc_svec,Bs2_fit,'-'); hold on;
    legend('data','linear fit');
    xlabel('dynamic pressure');
    ylabel('velocity');
    zlabel('B2');
    
    % generate the parameter-dependent matrices
    Ap = [p_coef_A11.Coefficients(1)*pbar_s+p_coef_A11.Coefficients(2)*veloc_s+p_coef_A11.Coefficients(3) ...
        p_coef_A12.Coefficients(1)*pbar_s+p_coef_A12.Coefficients(2)*veloc_s+p_coef_A12.Coefficients(3);...
        p_coef_A21.Coefficients(1)*pbar_s+p_coef_A21.Coefficients(2)*veloc_s+p_coef_A21.Coefficients(3)...
        p_coef_A22.Coefficients(1)*pbar_s+p_coef_A22.Coefficients(2)*veloc_s+p_coef_A22.Coefficients(3)];
    Bp = [p_coef_B1.Coefficients(1)*pbar_s+p_coef_B1.Coefficients(2)*veloc_s+p_coef_B1.Coefficients(3);
          p_coef_B2.Coefficients(1)*pbar_s+p_coef_B2.Coefficients(2)*veloc_s+p_coef_B2.Coefficients(3)];
    Cp = [1 0];

    % extract the constant and linear parts
    Ap0 = double(subs(Ap,{pbar_s,veloc_s},{0,0}));
    Ap1 = double(subs(Ap,{pbar_s,veloc_s},{1,0}))-Ap0;
    Ap2 = double(subs(Ap,{pbar_s,veloc_s},{0,1}))-Ap0;
    
    Bp0 = double(subs(Bp,{pbar_s,veloc_s},{0,0}));
    Bp1 = double(subs(Bp,{pbar_s,veloc_s},{1,0}))-Bp0;
    Bp2 = double(subs(Bp,{pbar_s,veloc_s},{0,1}))-Bp0;
end

% fprintf(1,'A11 error norm: %.3f,omega: %.3f pole bounds: [%.3f  %.3f],actual poles: %.3f, %.3f\n',theta1, abs(poles(1)), -0.7*(theta1+4+offset_upp), -0.7*(theta1+4+offset_low),poles); 

%  desired dynamics
% note that pbar_s =
% (pbar-pbar_lower_bnd)/(pbar_upper_bnd-pbar_lower_bnd)*2-1

 
% pbar = (pbar_s+1)*(pbar_upper_bnd-pbar_lower_bnd)/2+pbar_lower_bnd;
% Am = [0 1; -omega_p^2 2*zeta*omega_p];
