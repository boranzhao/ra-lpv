
clear;
% close all;
load('parameters_twoparas_ss_large_omega_range.mat','X','L');

lmi_formulation = 2; %1 for standard version based on [Apk 1998], 2 for refined version that allows for parameter-dependent X & Y
include_filter = 1;  % 1 for including a filter, 0 for not including the filter 

K =30; % filter bandwidth;
xi =1; % tuning parameter for the trade off between unmatched uncertainty compensation and stability satisfaction 

syms theta  theta2 %scheduling parameter
%%%%%%%%%%%%%%%%% for just testing %%%%%%%%%%%%%
% theta = -1; theta2 =0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
A = [-0.966 0.941; -3.441 -1.302]+theta*[-0.695 -0.022; -2.988 -0.884]+theta2*[-0.00377  -0.000767; 0.0864  -0.00373]; % Note that currently this is the open-loop A matrix
Bm = [-0.00202-0.00144*theta-6.13e-7*theta2; -0.265-0.241*theta-1.35e-4*theta2];
% formulate the closed-loop A matrix
Xtheta = X(:,:,1) +theta*X(:,:,2) +theta2*X(:,:,3);
Ltheta = L(:,:,1) +theta*L(:,:,2) +theta2*L(:,:,3);
A = A+Bm*Ltheta/Xtheta;

% Bm = [1;3];
% Bum = [0.264+0.241*theta; -0.002+0.001*theta];
Bum  = [Bm(2);-Bm(1)]; %
% Cm = [1 0];
Cm = [1 0;0 1];

% We = makeweight(10,1000,0.1)/10; 

% We = ss(tf([1 1],[1 1]));
nz = size(Cm,1);
% We = ss(tf(1))*eye(size(Cm,1)); 
We = ss(tf(1))*[1 0; 0 xi]; 
nx = size(A,1);
nu = 1;
ny = 1;

nw = 1; % dimension of the d/isturbance
nWe = order(We);


if ~include_filter
    Gsym.A = [A        zeros(nx) zeros(nx,nWe); 
              zeros(nx) A         zeros(nx,nWe);
              We.b*Cm   We.b*Cm   We.a];
    Gsym.B1 = [Bum;zeros(nx,1); zeros(nWe,1)];
    Gsym.B2 = [zeros(nx,nu); Bm; zeros(nWe,nu)];

    Gsym.C1 = [We.d*Cm We.d*Cm We.c];
    Gsym.D11 = zeros(nz,nw);
    Gsym.D12 = zeros(nz,nu);
    Gsym.C2 = zeros(nw,nx+nx+nWe);
    Gsym.D21 = 1;
    Gsym.D22 = zeros(ny,nu);
else
    sysF = ss(-K,1,K,0); nf=1;
    Gsym.A = [A            zeros(nx)    zeros(nx,nf) zeros(nx,nWe); 
              zeros(nx)    A            Bm*sysF.c   zeros(nx,nWe);
              zeros(nf,nx) zeros(nf,nx) sysF.a      zeros(nf,nWe);
              We.b*Cm       We.b*Cm     zeros(nWe,nf) We.a];
    Gsym.B1 = [Bum;zeros(nx+nf+nWe,1)];
    Gsym.B2 = [zeros(nx*2,nu); sysF.b; zeros(nWe,nu)];

    Gsym.C1 = [We.d*Cm We.d*Cm zeros(nz,nf) We.c];
    Gsym.D11 = zeros(nz,nw);
    Gsym.D12 = zeros(nz,nu);
    Gsym.C2 = zeros(nw,nx+nx+nf+nWe);
    Gsym.D21 = 1;
    Gsym.D22 = zeros(ny,nu);
end


XY = 3; %1: X is parameter-dependent and Y is constant; 2: X is constant and Y is parameter-dependent

%%%%%%%%%%%%%%%%%%%%%% for testing %%%%%%%%%%%%%%%%%%%
% ThetaGrid = 1;   % range of scaled dynamic pressure)
% Theta2Grid = 1;  % range of scaled aircraft speed. 
% d_Theta = [0];
% d_Theta2 = [0];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ThetaGrid = -1:0.2:1;   % range of scaled dynamic pressure)
Theta2Grid = -1:0.4:1;  % range of scaled aircraft speed. 
d_Theta = [-0.02 0.02];
d_Theta2 = [-0.05 0.05];
%%% original value
% d_Theta = [-0.01 0.01];
% d_Theta2 = [-0.05 0.05];
%%% 
F_theta = @(x) [1 x(1) x(2)]; %function for PD matrices
d_F_theta = @(x) [0 1 1];
FthetaNum = [1 1 1];

% F_theta = @(x) [1]; %function for PD matrices
% d_F_theta = @(x) [0];
% FthetaNum = [1 0 0];
% 
% F_theta = @(x) [1 x(1)]; %function for PD matrices
% d_F_theta = @(x) [0 1];
% FthetaNum = [1 1 0];

% F_theta = @(x) [1 x(1) x(1)^2 x(2)]; %function for PD matrices
% d_F_theta = @(x) [0 1 2*x(1) 1];
% FthetaNum = [1 2 1];

% Am = A;
% Hm = ss(Am,Bm,Cm,0);
% Hum = ss(Am,Bum,Cm,0);
% Kinv = -Hm\Hum;
% zero(Hm)

%% induced L2 norm minimization: no problem
% ga = ss(Gsym.A, [Gsym.B1 Gsym.B2], [Gsym.C1;Gsym.C2],[Gsym.D11 Gsym.D12; Gsym.D21 Gsym.D22]);
% NMEAS = size(Gsym.C2,1);
% NCONS = size(Gsym.B2,2);
% opts = hinfsynOptions('Method','LMI','Display','on');
% [Khinf,CL,gam] = hinfsyn(ga,NMEAS,NCONS,opts);
% 
% figure('name','controller'); bode(Khinf,Kinv);legend('Infinity norm min','Dynamic inversion');
% figure('name','Closed loop system'); bode(Hum+Hm*Khinf,Hum+Hm*Kinv);legend('Infinity norm min','Dynamic inversion');
% return;
%% peak-to-peak gain minimization. 
LAMBDA = 3.6; 3.2:0.2:4;  0.1:.1:2; %1.3;  0.1:0.1:2; %  0.25; 
ZETAOPT = inf*ones(1,length(LAMBDA));
% lamopt = LAMBDA(index2)
% X = Xopt(:,:,1);
% Y = Yopt(:,:,1);
% M = eye(size(X));
% N = eye(size(X)) - Y*X;
% Akhat = Khatopt.A(:,:,1);
% Bkhat = Khatopt.B(:,:,1);
% Ckhat = Khatopt.C(:,:,1);
% Dkhat = Khatopt.D(:,:,1);

plant_paras.ThetaGrid = ThetaGrid;
plant_paras.Theta2Grid = Theta2Grid;
plant_paras.d_Theta ={d_Theta,d_Theta2};
design_paras.XY_PD = XY;
% parameter-independent
% design_paras.F_theta = @(theta) [1 0];
% design_paras.d_F_theta =  @(theta) [0 0];
% design_paras.FthetaNum = [1 1 0];

%%% affine
design_paras.F_theta = F_theta;
design_paras.d_F_theta =  d_F_theta;
design_paras.FthetaNum = FthetaNum;
%%% quadratic
% % design_paras.F_theta = @(theta) [1 theta(1) theta(1)^2];
% % design_paras.d_F_theta =  @(theta) [0 1 2*theta(1)];
% % design_paras.FthetaNum = [1 2 0];
opt_yalmip = sdpsettings('verbose',1,'solver','mosek'); %do not display the optimization process if verbose is 0 
for i = 1:length(LAMBDA)
    lambda = LAMBDA(i);
    design_paras.lambda = lambda;
    if lmi_formulation == 1
        opt = lpv_ppg(Gsym,plant_paras,design_paras,opt_yalmip);
    else
        opt = lpv_ppg2(Gsym,plant_paras,design_paras,opt_yalmip);
    end
    ZETAOPT(i) = opt.zeta;
end
[zetamin,index] = min(ZETAOPT)
lambda = LAMBDA(index)

We = We.d;
filename = sprintf('unmatched_comp_filter_%0.2f_K30.mat',xi) ;      % 1.618
save(filename,'opt','lambda','zetamin','index','We')
return;
%% check the stability of the controller
for rou  = -1:0.01:1
    Frou = design_paras.F_theta([rou 0]);
    X = opt.X(:,:,1);
    Y = opt.Y(:,:,1);
    Akhat = opt.Ahat(:,:,1);
    Bkhat = opt.Bhat(:,:,1);
    Ckhat = opt.Chat(:,:,1);
    Dkhat = opt.Dhat(:,:,1);
    for i=2:length(Frou)
        X = X + Frou(i)*opt.X(:,:,i);
        Y = Y + Frou(i)*opt.Y(:,:,i);
        Akhat = Akhat + Frou(i)*opt.Ahat(:,:,i);
        Bkhat = Bkhat + Frou(i)*opt.Bhat(:,:,i);
        Ckhat = Ckhat + Frou(i)*opt.Chat(:,:,i);
        Dkhat = Dkhat + Frou(i)*opt.Dhat(:,:,i);
    end
    S = X-Y;
    Dk = Dkhat;
    Bk = Bkhat;
    Ck = Ckhat/S;
    Ak = Akhat/S;
    
    if max(real(eig(Ak)))>0
       fprintf(1,'controller unstable at theta =%.3f',rou);
    end
end


%% Find the LTI controller
X = opt.X(:,:,1);
Y = opt.Y(:,:,1);

Akhat = opt.Ahat(:,:,1);
Bkhat = opt.Bhat(:,:,1);
Ckhat = opt.Chat(:,:,1);
Dkhat = opt.Dhat(:,:,1);

Am = double(vpa(A,5));
if lmi_formulation == 1
    M = eye(size(X));
    N = eye(size(X)) - Y*X;
    Dk = Dkhat;
    Ck = (Ckhat - Dk*Gsym.C2*X)/M';
    Bk = N\(Bkhat - Y*Gsym.B2*Dk);
    Ak = N\(Akhat - Bkhat*Gsym.C2*X - Y*Gsym.B2*Ckhat - Y*(double(Gsym.A)-Gsym.B2*Dk*Gsym.C2)*X)/M';
else
    S = X-Y;
    Dk = Dkhat;
    Bk = Bkhat;
    Ck = Ckhat/S;
    Ak = Akhat/S;
end

Kff = ss(Ak,Bk,Ck,Dk);
figure('name','Controller'); bode(Kinv,Kff);legend('dynamic inversion','lmi');
figure('name','Closed loop system'); bode(Hum+Hm*Kinv, Hum+Hm*Kff);legend('dynamic inversion','lmi');
