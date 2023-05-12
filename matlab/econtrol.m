
% Ergodic controller for a spatial distribution (two gaussians in the
% mixture model)

% Same parameters as in (see calinon.m):
% Calinon, S. (2020). Mixture Models for the Analysis, Edition, and 
% Synthesis of Continuous Time Series. In: Bouguila, N., Fan, W. (eds) 
% Mixture Models and Applications. Unsupervised and Semi-Supervised 
% Learning. Springer, Cham. 
% https://doi.org/10.1007/978-3-030-23876-6_3

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

N=2000;
dt=1e-2;

xlim=[-L/2 L/2;-L/2 L/2]; % limits
ulim=[0 10];

x0=[.1;.3];%[.5;.1]; % initial guesses

args.L=L; % wrapping arguments for AUX functions
args.D=D;
args.alpha=alpha;
args.Am=Am;
args.Mu=Mu;
args.Sigma=Sigma;
args.K=K;
args.N=1;

phi_k_val=nan(length(K),1);
x=x0;
debug.x=nan(2,N);
debug.u=nan(2,N-1);
debug.x(:,1)=x;

wt=zeros(length(K)^D,1);
phi_k_val=nan(length(K)^D,1);
f_k_x_val=nan(length(K)^D,1);
df_k_x_val=nan(2,length(K)^D);
Lambda_k=zeros(length(K)^D);

for t=1:N-1
    utilde=0;


    for k=1:length(K)
        f_k_x_val((k-1)*length(K)+1:k*length(K))=...
            f_k_x(K(:,:,k),x,args);
        phi_k_val((k-1)*length(K)+1:k*length(K))=...
            phi_k(K(:,:,k),args);
        df_k_x_val(:,(k-1)*length(K)+1:k*length(K))=...
            df_k_x(K(:,:,k),x,args)*L^D;

        for j=1:length(args.K)
            Lambda_k((k-1)*length(K)+j,(k-1)*length(K)+j)=...
                (sum(args.K(:,j,k).^2)+1).^(-(args.D+1)/2);        
        end
    end
    wt=wt+f_k_x_val;
    utilde=utilde-1*df_k_x_val*Lambda_k*(wt/t-phi_k_val);

    u=utilde*max(ulim)/(norm(utilde)+1E-1);
    x=x+u*dt;

    debug.x(:,t+1)=x;
    debug.u(:,t)=u;
end


%% visualization

figure;
plot(debug.x(1,:),debug.x(2,:))
clear xlim;
xlim([0 1]);
ylim([0 1]);
hold on;
plot(Mu(1,1),Mu(2,1),'*');
plot(Mu(1,2),Mu(2,2),'*');
