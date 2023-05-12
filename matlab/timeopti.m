
% Time-optimal ergodic search with IPOPT (CasADi framework for
% transcription, RK4 for numerical simulation)

%% definitions

k=9; % number of freq. in Fourier transform
D=2; % dimension, e.g., 2D, 3D, etc.
[K(1,:,:) K(2,:,:)]=ndgrid(0:1:k,0:1:k); % set of indices
L=2; % period
Ai={[1 0;0 1],[1 0;0 -1],[-1 0;0 1],[-1 0;0 -1]};

if D~=2
    error("myComponent:notImplemented",strcat("Error. \nlinear transf",...
          "ormation matrices Am, indeces K are not yet implemented for",...
          " dimensions other than 2 (dimension is %d)"), D);
end
Am=@(i) cell2mat(Ai(i)); % linear transformation matrices

Mu(:,1)=[.5;.7];
Mu(:,2)=[.6;.3];
Sigma(:,:,1)=[.3;.1]*[.3;.1]'*5e-1+eye(D)*5e-3;
Sigma(:,:,2)=[.1;.2]*[.1;.2]'*3e-1+eye(D)*1e-2;
alpha=[.5;.5];

N=200;
dt=1e-3;
tf=N*dt;

xlim=[-L/2 L/2;-L/2 L/2]; % limits
ulim=[-6 6];

x0=[.1;.3];%[.5;.1]; % initial guesses
xf=[.9;-.1];%[2;3.2];

args.L=L; % wrapping arguments for AUX functions
args.D=D;
args.alpha=alpha;
args.Am=Am;
args.Mu=Mu;
args.Sigma=Sigma;
args.K=K;
args.N=N;

fdot=@(x,u) eye(2)*u; % craft's dynamics


%% CasADi
addpath('~/casadi');

import casadi.*;
opti=casadi.Opti();

X=opti.variable(2,N); % states and controls
U=opti.variable(2,N-1);
T=opti.variable(); % time variable (i.e., final time)

opti.subject_to(T>=0); % setting time constraint
opti.subject_to(X(:,end)==xf); % setting final guess
opti.subject_to(X(:,1)==x0); % setting initial guess


opti.subject_to(ulim(1)<=U(:,:)<=ulim(2)); % setting state/control bounds
opti.subject_to(xlim(1,1)<=X(1,:)<=xlim(1,2));
opti.subject_to(xlim(2,1)<=X(2,:)<=xlim(2,2));

opti.set_initial(T,dt); % setting initial guess T into CasADi

for k=1:N-1
    
    k1=fdot(X(:,k),        U(:,k));
    k2=fdot(X(:,k)+dt/2*k1,U(:,k));
    k3=fdot(X(:,k)+dt/2*k2,U(:,k));
    k4=fdot(X(:,k)+dt*k3,  U(:,k));
    xdot=X(:,k)+dt/6*(k1+2*k2+2*k3+k4); % Runge-Kutta 4th order integration
    opti.subject_to(X(:,k+1)==xdot); % dynamics constraint
end

epsilon_xt_val=epsilon_xt(X,args);
opti.subject_to(epsilon_xt_val<=.0035); % ergodicity metrics constraint

opti.minimize(T); % objective, i.e., minimum time

opti.solver('ipopt');
sol=opti.solve();

%% visualization

figure;
xres=opti.debug.value(X);
plot(xres(1,:)',xres(2,:)');
clear xlim;
xlim([0 1]);
ylim([0 1]);
hold on;
plot(Mu(1,1),Mu(2,1),'o');
plot(Mu(1,2),Mu(2,2),'o');

