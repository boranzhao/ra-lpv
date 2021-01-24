 clear;
load('parameters_twoparas_ss_large_omega_range.mat','X','L');
syms theta  theta2 %scheduling parameter
%%%%%%%%%%% should comment this %%%%%%%%%
% theta = 1;theta2 = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
A = [-0.966 0.941; -3.441 -1.302]+theta*[-0.695 -0.022; -2.988 -0.884]+theta2*[-0.00377  -0.000767; 0.0864  -0.00373]; % Note that currently this is the open-loop A matrix
Bm = [-0.00202-0.00144*theta-6.13e-7*theta2; -0.265-0.241*theta-1.35e-4*theta2];
% formulate the closed-loop A matrix
Xtheta = X(:,:,1) +theta*X(:,:,2) +theta2*X(:,:,3);
Ltheta = L(:,:,1) +theta*L(:,:,2) +theta2*L(:,:,3);
A = A+Bm*Ltheta/Xtheta;

% Bm = [1;3];
% Bum = [0.264+0.241*theta; -0.002+0.001*theta];
Bum  = [Bm(2);-Bm(1)]; %
Cm = [1 0];
n=size(A,1);
nu = size(Bm,2);
ny = size(Cm,1);
nf=1;
n1 = 2*n+nf;
m = nu;

ThetaGrid = -1:0.2:1;   % range of scaled dynamic pressure)
Theta2Grid = -1:0.4:1;  % range of scaled aircraft speed. 
d_Theta = [-0.02 0.02];
d_Theta2 = [-0.05 0.05];

%%%%%%%%%%%%%%%%%%%%%% for testing %%%%%%%%%%%%%%%%%%%
% ThetaGrid = 1;   % range of scaled dynamic pressure)
% Theta2Grid = 1;  % range of scaled aircraft speed. 
% d_Theta = [0];
% d_Theta2 = [0];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% A = A;
% Hm = ss(A,Bm,Cm,0);
% Hum = ss(A,Bum,Cm,0);
% Kinv = -Hm\Hum;
% zero(Hm)
% return;

% F_theta = @(x) [1 x(1) x(2)]; %function for PD matrices
% d_F_theta = @(x) [0 1 1];
% FthetaNum = [1 1 1];

F_theta = @(x) [1 x(1) x(1)^2 x(2)]; %function for PD matrices
d_F_theta = @(x) [0 1 2*x(1) 1];
FthetaNum = [1 2 1];

% F_theta = @(x) [1 x(1) x(1)^2 x(2) x(2)^3]; % does not improve
% d_F_theta = @(x) [0 1 2*x(1) 1 3*x(2)^2];
% FthetaNum = [1 2 2];

plant_paras.ThetaGrid = ThetaGrid;
plant_paras.Theta2Grid = Theta2Grid;
plant_paras.d_Theta ={d_Theta,d_Theta2};

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
% opt_yalmip = sdpsettings('verbose',1,'solver','sedumi','sedumi.maxiter',200); %do not display the optimization process if verbose is 0 
% opt_yalmip = sdpsettings('verbose',1,'solver','sdpt3','sdpt3.maxit',200); %do not
opt_yalmip = sdpsettings('verbose',1,'solver','mosek'); %do not
% opt_yalmip = sdpsettings('verbose',1,'solver','lmilab'); 

% Control settings
w = 1; % unknown input gain, assume to be 1
K= 30; % filter bandwidth 
C = ss(-w*K,eye(m),w*K,0);
%% compute L1 norm bound of Gum and Gxum
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n=size(A,1);n1 = 2*n+nf;
opt_yalmip = sdpsettings('verbose',0,'solver','mosek'); %do not
% opt_yalmip = sdpsettings('verbose',1,'solver','sedumi','sedumi.maxiter',200); %do not display the op
A_bH = ones(n1,n1)*nan;
B_bH = ones(n1,n-m)*nan;
C_bH = ones(m,n1)*nan;
D_bH = ones(m,n-m)*nan;

%%%%%%%%%%%%%%%%%%%%%%%%% for test only %%%%%%%%%
% A_bH = rand(n1,n1); 
% A_bH = A_bH + A_bH';
% A_bH = A_bH-(max(abs(eig(A_bH))))*eye(n1); 
% B_bH = rand(n1,n-m);
% C_bH = rand(m,n1);
% D_bH = rand(m,n-m);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Gxum
for id_test = 1:2
    % load unmatched_comp_filter_09.mat
    load unmatched_comp_filter_1.00_K30.mat
    
    Psym.A = [A zeros(n,n+m+n1);
              zeros(n) A Bm*w*K zeros(n,n1);
              zeros(m,2*n) -w*K C_bH;
              zeros(n1,2*n+m) A_bH];
    Psym.B = [Bum; zeros(n,n-m); D_bH; B_bH];
    if id_test == 1
        Psym.C = [eye(n) eye(n) zeros(n,m+n1)];
        Psym.D = zeros(n,n-m);
    else
        %%%%%%% for evaluting Gum= Hum-HmCHbar %%%%%%%
        Psym.C = [Cm Cm zeros(1,m+n1)]; 
        Psym.D = zeros(1,n-m);
    end

    Psym.X = opt.X;
    Psym.Y = opt.Y;
    Psym.Ahat = opt.Ahat;
    Psym.Bhat = opt.Bhat;
    Psym.Chat = opt.Chat;
    Psym.Dhat = opt.Dhat;
    LAMBDA = 3.6; 3:0.2:4; 2:0.2:3.0; %  1:.2:2; 2:0.2:3; 2.6; %:0.2:3; .001:1:6;%:0.5:6; %1.3;  0.1:0.1:2; 1.8; 
    ZETAOPT = inf*ones(1,length(LAMBDA));

    %%%%%%%%%%%%%%%%%%%%%%% verify whether the state space matrices are correct: no problem
    % sys_P= ss(Psym.A,Psym.B,Psym.C,Psym.D);
    % bH = ss(A_bH,B_bH,C_bH,D_bH);
    % Hxm = ss(A,Bm,eye(n),0);
    % C = ss(-w*K,eye(m),w*K,0);
    % Hxum = ss(A,Bum,eye(n),0);
    % sys_P2= Hxum+Hxm*C*bH;
    % figure;
    % % bode(sys_P,sys_p2);
    % bode(sys_P);hold on;
    % bode(sys_P2)
    % return;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for i = 1:length(LAMBDA)
        lambda = LAMBDA(i);
        design_paras.lambda = lambda;   
        opt1 = lpv_ppg_analysis(Psym,plant_paras,design_paras,opt_yalmip);
        ZETAOPT(i) = opt1.zeta;
    end
    [zetamin_Gum,id] = min(ZETAOPT);    
    lambda_Gum = LAMBDA(id);
    if id_test == 1
        fprintf(1,'||Gxum||_ppg: %.3f\n', zetamin_Gum);
    else
        fprintf(1,'||Gum||_ppg: %.3f\n', zetamin_Gum);
    end
end
return;






%% compute L1 norm bound of Gxm
Psym.A = [A -Bm*w*K; zeros(m,n) -w*K];
Psym.B = [Bm;eye(m)];
Psym.C = [eye(n) zeros(n,m)];
Psym.D = zeros(n,m);

%%%%%%%%%%%% compute the L1 norm for comparison
% model = ss(Psym.A,Psym.B,Psym.C,Psym.D);
% l1norm_Gm = L1norm(model)
% return;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

LAMBDA = 3.4:.1:3.8; 3.5; %3.2:.1:3.7; %3.5:.1:3.7; %1.3;  0.1:0.1:2; 1.8; 
ZETAOPT = inf*ones(1,length(LAMBDA));
for i = 1:length(LAMBDA)
    lambda = LAMBDA(i);
    design_paras.lambda = lambda;   
    opt = lpv_ppg_analysis(Psym,plant_paras,design_paras,opt_yalmip);
    ZETAOPT(i) = opt.zeta;
end
[zetamin_Gm,index] = min(ZETAOPT)
lambda_Gm = LAMBDA(index)
% save('Gum_bound.mat','opt','lambda','zetamin','index')
return;

%% compute L1 norm bound of HxmCkg
kg = -inv(Cm/A*Bm);
Psym.A = [A Bm*w*K;
          zeros(m,n) -w*K];
Psym.B = [zeros(n,m); kg];
Psym.C = [eye(n) zeros(n,m)];
Psym.D = zeros(n,m);
LAMBDA = 3.1:.1:3.3; %2.2:.1:3; %1.3;  0.1:0.1:2; 1.8; 
ZETAOPT = inf*ones(1,length(LAMBDA));
for i = 1:length(LAMBDA)
    lambda = LAMBDA(i);
    design_paras.lambda = lambda;   
    opt = lpv_ppg_analysis(Psym,plant_paras,design_paras,opt_yalmip);
    ZETAOPT(i) = opt.zeta;
end
[zetamin_HCKg,index] = min(ZETAOPT)
lambda_HCKg = LAMBDA(index)
return;






%% Compute the Lipschitz constant
b0 = 3;
L= 20*pi;
norm_Bum = ones(1,10000)*nan;
norm_Bm = ones(1,10000)*nan;
norm_Buminv = ones(1,10000)*nan;
norm_Bminv = ones(1,10000)*nan;
i = 0;
for theta = -1:0.1:1
    for theta2 = -1:0.1:1
        i = i+1;
        Bm0 = double(subs(Bm)); 
        Bum0 = double(subs(Bum));  
        Bminv0 = pinv(Bm0);
        Buminv0 = pinv(Bum0);
        norm_Bm(i) =norm(Bm0);
        norm_Bum(i) =norm(Bum0);
        norm_Bminv(i) =norm(Bminv0);
        norm_Buminv(i) =norm(Buminv0);
    end
end
%%
norm_Bm(isnan(norm_Bm))=[];
norm_Bum(isnan(norm_Bum))=[];
norm_Bminv(isnan(norm_Bminv))=[];
norm_Buminv(isnan(norm_Buminv))=[];
B10 = max(norm_Bminv)*b0
B20 = max(norm_Buminv)*b0
L1= max(norm_Bminv)*L
L2= max(norm_Buminv)*L

%% compute rho_in
load('parameters_twoparas_ss_large_omega_range.mat','X'); % inverse of Lyapunov matrix
max_eig_P = ones(1,1000)*nan;
min_eig_P = ones(1,1000)*nan;
i=0;
for theta = -1:0.1:1
    for theta2 = -1:0.1:1
        i = i+1;
        Ftheta = F_theta([theta;theta2]);
        
        Xeval = X(:,:,1)*Ftheta(1) +  X(:,:,2)*Ftheta(2)+X(:,:,3)*Ftheta(3);
%         P = inv(Xeval);
        eig_vec = eig(Xeval);
        max_eig_P(i) = 1/min(eig_vec);
        min_eig_P(i) = 1/max(eig_vec);        
    end
end
max_eig_P(isnan(max_eig_P)) = [];
min_eig_P(isnan(min_eig_P)) = []; 
lambda2 = max(max_eig_P);
lambda1 = min(min_eig_P);

rho0 = 0.1;
rhoin = rho0*sqrt(lambda2/lambda1)

%% solve the stability condition for rho_r
Gm_l1 = 0.0219;
HCkg_l1 = 5.3115;
Gum_l1= 5.617;
l0 = L2/L1;
b0 = max(B10,B20/l0);

r_max = 20/180*pi;

syms rhor
% L1 = 6.2832

assume(rhor>0);
% solve(Gm_l1+Gum_l1*l0 < (rhor-HCkg_l1*r_max-rhoin)/(L1*rhor+b0),rhor)
solve(Gm_l1 < (rhor-HCkg_l1*r_max-rhoin)/(L1*rhor+b0),rhor)

