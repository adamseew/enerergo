
% Energy-aware ergodic search, i.e., an ergodic controller that uses the
% controller from econtrol.m to derive a control action that respects the
% battery constraints while maximazing the coverage


answer=questdlg('Would you like to clean the environment?',...
    'Clean',...
	'Yes','No','No');
switch answer
case 'Yes'
    clear;
end
clear("answer");


%% definitions

D=2; % dimension, e.g., 2D, 3D, etc.
L=2; % period
epsilon=0.05; % tollerance interval
approx_f=1; % approximation factor, i.e., 1 no approximation. Use powers
            % of 10

if D~=2
    error("myComponent:notImplemented",strcat("Error. \nlinear transf",...
          "ormation matrices Am, indeces K are not yet implemented for",...
          " dimensions other than 2 (dimension is %d)"), D);
end

disp("starting gauss app")
gauss_app=gauss("on"); % starting Guassian mixture model designer app
try
    waitfor(gauss_app,"closed","yes"); % wait for Guassian mixture model 
                                       % designer app to terminate
catch
    delete(gauss_app.ErgodiccontrollerdesignerUIFigure)
    error("Data from gauss.mlapp are required")
end
Mu=gauss_app.Mu;
Sigma=gauss_app.Sigma;
alpha=gauss_app.alpha;
k=gauss_app.k; % number of freq. in Fourier transform
[K(1,:,:) K(2,:,:)]=ndgrid(0:1:k,0:1:k); % set of indices
N=gauss_app.N; % horizon
x0=gauss_app.x0; % initial guesses
xf=gauss_app.xf; % desired final point, e.g., recharge station

delete(gauss_app.ErgodiccontrollerdesignerUIFigure)
if length(alpha)>1
    disp(strcat("there are ",string(length(alpha)),...
        " gaussians in the model"))
else
    disp("there is one gaussian in the model")
end
clear("gauss_app")

dt=1e-2;

xlim_=[-L/2 L/2]; % limits (square)
ulim=[0 4];

args.L=L; % wrapping arguments for AUX functions
args.D=D;
args.Mu=Mu;
args.Sigma=Sigma;
args.K=K;
args.N=N;

x=x0;
debug.x=nan(2,N);
debug.u=nan(2,N-1);
debug.x(:,1)=x;
debug.alpha=nan(length(Mu),1);

wt=zeros(length(K)^D,1);
Lambda_k=zeros(length(K)^D);


%% battery

% battery model from: 

min_b=.05; % minimum battery constraint

args.C1=4.78*1e2;
args.R1=2.85*1e-2; % first RC element in the ECM
args.C2=1.83*1e4;
args.R2=4.44*1e-2;
args.Q=.250; % charge
args.V=3; % voltage
args.Rs=5.55*1e4; % R in the ECM

I=3; % load (fixed)
v1=3;
v2=3;
z=100; % initial SoC
x0_b=[v1;v2;z];
debug.z=z;
Voc=args.V+x0_b(1)+x0_b(2)+I*args.Rs;
debug.Voc=Voc;


%% CasADi

addpath('~/casadi');

import casadi.*;
opti=casadi.Opti();

PHI_K_VAL=opti.variable(length(K)^D,1); % aux variables
f_k_x_val=[];
df_k_x_val=[];

X=opti.variable(2,N/approx_f); % state
U=opti.variable(2,N/approx_f-1); % input
ALPHA=opti.variable(length(alpha),1); % control variable
X0_B3=opti.variable(1,1); % battery at final time step (just to add it as a 
                          % constraint)

opti.set_initial(ALPHA,alpha); % setting initial guess

opti.subject_to(sum(ALPHA)<=1); % constraints on control variable
opti.subject_to(ALPHA>0);
opti.subject_to(X(:,1)==x); % state initial guess

args.alpha=ALPHA;

disp('entering main loop, i.e., transcription')
wb=waitbar(0,''); % initializing progressbar

for t=1:N-1
    time_start=tic;
    utilde=0;

    for k=1:length(K)
       
        f_k_x_val=[f_k_x_val;...
            f_k_x(K(:,:,k),x,args)];
        df_k_x_val=[df_k_x_val ...
            df_k_x(K(:,:,k),x,args)*L^D];

        if t==1 % phi_k is not dependent on x, i.e., it is the same at
                % each t
            PHI_K_VAL((k-1)*length(K)+1:k*length(K))=...
                phi_k(K(:,:,k),args);
            for j=1:length(args.K)
                Lambda_k((k-1)*length(K)+j,(k-1)*length(K)+j)=...
                    (sum(args.K(:,j,k).^2)+1).^(-(args.D+1)/2);        
            end
        end
    end
    wt=wt+f_k_x_val;
    utilde=utilde-1*df_k_x_val*Lambda_k*(wt/t-PHI_K_VAL);
    
    u=utilde*max(ulim)/(norm(utilde)+1E-1);
    x=x+u*dt;
    
    if mod(t,approx_f)==0
        opti.subject_to(U(:,t/approx_f)==u);
        opti.subject_to(X(:,t/approx_f+1)==x);
    end

    x0_b=x0_b+dt*thevenin(x0_b,1,args); % battery model
    Voc=args.V+x0_b(1)+x0_b(2)+I*args.Rs;
    debug.z=[debug.z;x0_b(3)];
    debug.Voc=[debug.Voc;Voc];

    f_k_x_val=[];
    df_k_x_val=[];

    remaining=double(toc(time_start))*(N-1-t)/60; % data for progressbar
    waitbar(t/(N-1),wb,...
        strcat(string(round(remaining))," minutes remaining"));
    
end

close(wb);
disp('exiting main loop')

opti.subject_to(X0_B3==x0_b(3));
opti.subject_to(norm(X(:,end)-xf)<=epsilon); % final constraint (over-
                                             % writes previous one from the 
                                             % loop)
opti.subject_to(X0_B3>=min_b); % battery constraint
opti.minimize(-1*sum(ALPHA)); % cost (best coverage)

opti.solver('ipopt');

disp('starting the solver')

try
    sol=opti.solve();

    disp('optimal solution found')
    debug.x=sol.value(X);
    debug.u=sol.value(U);
    debug.alpha=sol.value(ALPHA);

    debug.phi_k_val=sol.value(PHI_K_VAL);
catch
    disp('optimal solution NOT found')

    debug.x=opti.debug.value(X);
    debug.u=opti.debug.value(U);
    debug.alpha=opti.debug.value(ALPHA);

    debug.phi_k_val=opti.debug.value(PHI_K_VAL);
end

fprintf("alphas: ")
fprintf('%d ', debug.alpha); fprintf('\n');


%% visualization

% to visualize the probability distribution in th plot

sampl=linspace(xlim_(1),xlim_(2),length(K)^D); % build samples per each 
                                               % point in space
[probx,proby]=ndgrid(sampl,sampl);
probx=probx(:);
proby=proby(:);
f_k_val_prob=nan((length(K)^D)^2,1);
f_k_val_prob_val=nan(length(K)^D,1);

for j=1:(length(K)^D)^2
    for k=1:length(K)
        f_k_val_prob_val((k-1)*length(K)+1:k*length(K))=...
            f_k_x(K(:,:,k),[probx(j);proby(j)],args);
    end
    f_k_val_prob(j)=f_k_val_prob_val'*debug.phi_k_val;
end

figure;
h=surf(...
       reshape(probx,length(K)^D,length(K)^D),...
       reshape(proby,length(K)^D,length(K)^D),...
       reshape(f_k_val_prob,length(K)^D,length(K)^D),...
       'FaceColor','interp','EdgeAlpha',0,'FaceAlpha',.7...
      );
colormap(flipud(bone(30)));
z=get(h,'ZData');
set(h,'ZData',z-10);
view(2);
grid on;
hold on;
for j=1:length(alpha)
    plot(Mu(1,j),Mu(2,j),'g^');
end
xlim([0 1]);
ylim([0 1]);
axis square;
set(gcf,'color','w');
ax1=gca;
ax1.YGrid='on';
ax1.Layer='top';
ax1.GridLineStyle=':';
ax1.GridAlpha=.25;
pl=plot(debug.x(1,:),debug.x(2,:),'blue');
plot(debug.x(1,1),debug.x(2,1),'rx');
plot(debug.x(1,end),debug.x(2,end),'rs');
% remove the comments bellow for an animation
% for j=2:N
    % pause(.05);
    % delete(pl);
    % pl=plot(debug.x(1,1:j),debug.x(2,1:j),'blue');
% end
set(ax1,'XTick',get(ax1,'YTick'));

