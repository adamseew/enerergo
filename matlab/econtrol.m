
% Ergodic controller for a spatial distribution (four gaussians in the
% mixture model)

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
Mu(:,1)=[.33333;.14285];
Mu(:,2)=[.83332;.1    ];
Mu(:,3)=[.4    ;.57142];
Mu(:,4)=[.83332;.85713];
Sigma(:,:,1)=[.1;.1]*[.1;.1]'*1e-3+eye(D)*1e-3;
Sigma(:,:,2)=Sigma(:,:,1);
Sigma(:,:,3)=Sigma(:,:,1);
Sigma(:,:,4)=Sigma(:,:,1);
alpha=[.25;.25;.25;.25];

N=2000;
dt=1e-2;

xlim=[-L/2 L/2]; % limits (square)
ulim=[0 2];

x0=[.05;.05];%[.5;.1]; % initial guesses

args.L=L; % wrapping arguments for AUX functions
args.D=D;
args.alpha=alpha;
args.Am=Am;
args.Mu=Mu;
args.Sigma=Sigma;
args.K=K;
args.N=1;

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
        df_k_x_val(:,(k-1)*length(K)+1:k*length(K))=...
            df_k_x(K(:,:,k),x,args)*L^D;

        if t==1 % phi_k is not dependent on x, i.e., it is the same at
                % each t
            phi_k_val((k-1)*length(K)+1:k*length(K))=...
                phi_k(K(:,:,k),args);
            for j=1:length(args.K)
                Lambda_k((k-1)*length(K)+j,(k-1)*length(K)+j)=...
                    (sum(args.K(:,j,k).^2)+1).^(-(args.D+1)/2);        
            end
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


% to visualize the probability distribution in th plot
sampl=linspace(xlim(1),xlim(2),length(K)^D); % build samples per each point
                                             % in space
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
    f_k_val_prob(j)=f_k_val_prob_val'*phi_k_val;
end

figure(1);
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
clear xlim;
xlim([0 1]);
ylim([0 1]);
axis square;
set(gcf,'color','w');
ax1=gca;
ax1.YGrid='on';
ax1.Layer='top';
ax1.GridLineStyle=':';
ax1.GridAlpha=.25;
pl=plot(debug.x(1,1),debug.x(2,1),'blue');
for j=2:N
    pause(.1);
    delete(pl);
    pl=plot(debug.x(1,1:j),debug.x(2,1:j),'blue','LineWidth',2);
end
set(ax1,'XTick',get(ax1,'YTick'));

csvwrite('traj.csv',(debug.x(:,1:10:N).*3)');

return

%% battery
args.C1=4.78*1e2;
args.C2=1.83*1e4;
args.R1=2.85*1e-2;
args.R2=4.44*1e-2;
args.Q=.250;
args.V=3;
args.Rs=5.55*1e4;

I=3;
v1=3;
v2=3;
z=100;
x0=[v1;v2;z];
debug.z=z;
Voc=args.V+x0(1)+x0(2)+I*args.Rs;
debug.Voc=Voc;
for j=1:N-1
    x0=x0+dt*thevenin(x0,1,args);
    Voc=args.V+x0(1)+x0(2)+I*args.Rs;
    debug.z=[debug.z;x0(3)];
    debug.Voc=[debug.Voc;Voc];
end

