
% Calinon, S. (2020). Mixture Models for the Analysis, Edition, and 
% Synthesis of Continuous Time Series. In: Bouguila, N., Fan, W. (eds) 
% Mixture Models and Applications. Unsupervised and Semi-Supervised 
% Learning. Springer, Cham. 
% https://doi.org/10.1007/978-3-030-23876-6_3

nbData = 2000; % number of datapoints
nbFct = 10; % number of basis functions along x and y
nbVar = 2; % dimension of datapoint
nbStates = 2; % number of Gaussians to represent the spatial distribution
sp = (nbVar + 1) / 2; % sobolev norm parameter
dt = 1E-2; % time step
xlim = [0; 1]; % domain limit for each dimension (considered to be 1 for each dimension in this implementation)
L = (xlim(2) - xlim(1)) * 2; % size of [-xlim(2),xlim(2)]
u_max = 1E1; % maximum speed allowed 

% desired spatial distribution represented as a mixture of Gaussians
Mu(:,1) = [.5; .7]; 
Sigma(:,:,1) = [.3;.1]*[.3;.1]' *5E-1 + eye(nbVar)*5E-3; %eye(nbVar).*1E-2; 
Mu(:,2) =  [.6; .3]; 
Sigma(:,:,2) = [.1;.2]*[.1;.2]' *3E-1 + eye(nbVar)*1E-2;

Priors = ones(1,nbStates) / nbStates; % mixing coefficients


% compute Fourier series coefficients phi_k of desired spatial distribution

[xx, yy] = ndgrid(1:nbFct, 1:nbFct);
rg = 0:nbFct-1;
[KX(1,:,:), KX(2,:,:)] = ndgrid(rg, rg);
Lambda = (KX(1,:).^2 + KX(2,:).^2 + 1)'.^-sp; % weighting vector (Eq.(15))

% explicit description of phi_k by exploiting the Fourier transform properties of Gaussians (optimized version by exploiting symmetries),
% by enumerating symmetry operations for 2D signal ([-1,-1],[-1,1],[1,-1] and [1,1]), and removing redundant ones -> keeping ([-1,-1],[-1,1])
op = hadamard(2^(nbVar-1));
op = op(1:nbVar,:);
% compute phi_k
kk = KX(:,:) * 2 * pi / L;
w_hat = zeros(nbFct^nbVar, 1);
for j=1:nbStates
	for n=1:size(op,2)
		SigmaTmp = diag(op(:,n)) * Sigma(:,:,j) * diag(op(:,n))'; % Eq.(20)
		w_hat = w_hat + Priors(j) * cos(kk' * diag(op(:,n)) * Mu(:,j)) .* exp(diag(-.5 * kk' * SigmaTmp * kk)); %Eq.(21)
	end
end
w_hat = w_hat / L^nbVar / size(op,2);

% fourier basis functions (for a discretized map)
nbRes = 100;
xm1d = linspace(xlim(1), xlim(2), nbRes); % spatial range for 1D
[xm(1,:,:), xm(2,:,:)] = ndgrid(xm1d, xm1d); % spatial range
phim = cos(KX(1,:)' * xm(1,:) .* 2 * pi / L) .* cos(KX(2,:)' * xm(2,:) * 2 * pi / L) * 2^nbVar; %Fourier basis functions
hk = [1; 2*ones(nbFct-1,1)];
HK = hk(xx(:)) .* hk(yy(:)); 
phim = phim .* repmat(HK,[1,nbRes^nbVar]);

% desired spatial distribution 
g = w_hat' * phim;

% ergodic control 
x = [.1; .3]; % initial position
wt = zeros(nbFct^nbVar, 1);
for t=1:nbData
	r.x(:,t) = x; % log data
	
	% fourier basis functions and derivatives for each dimension (only cosine part on [0,L/2] is computed since the signal is even and real by construction) 
	phi1 = cos(x * rg * 2 * pi / L); % Eq.(18)
	dphi1 = -sin(x * rg * 2 * pi / L) .* repmat(rg,nbVar,1) * 2 * pi / L;
	
	dphi = [dphi1(1,xx) .* phi1(2,yy); phi1(1,xx) .* dphi1(2,yy)]; % gradient of basis functions
	wt = wt + (phi1(1,xx) .* phi1(2,yy))' / L^nbVar;	% wt./t are the Fourier series coefficients along trajectory (Eq.(17))

	
	% controller with constrained velocity norm
	u = -dphi * (Lambda .* (wt/t - w_hat)); % Eq.(24)
	u = u * u_max / (norm(u)+1E-1); % velocity command
	
	x = x + u * dt; % update of position
	r.g(:,t) = (wt/t)' * phim; % reconstructed spatial distribution (for visualization)
	r.w(:,t) = wt/t; % fourier coefficients along trajectory (for visualization)
end

% plot

% x
% subplot(1,3,1); 

figure;
G = reshape(g,[nbRes,nbRes]);
G([1,end],:) = max(g); % add vertical image borders
G(:,[1,end]) = max(g); % add horizontal image borders
surface(squeeze(xm(1,:,:)), squeeze(xm(2,:,:)), zeros([nbRes,nbRes]), G, 'FaceColor','interp','EdgeColor','interp');
hold on;
plot(r.x(1,:), r.x(2,:), '-','linewidth',1,'color',[0 0 0]);
plot(r.x(1,1), r.x(2,1), '.','markersize',15,'color',[0 0 0]);
%axis([xlim(1),xlim(2),xlim(1),xlim(2)]); axis equal;
clear xlim;
xlim([0 1]);
ylim([0 1]);

% w
% subplot(1,3,2); hold on; axis off; title('$w$','interpreter','latex','fontsize',20);
% imagesc(reshape(wt./t,[nbFct,nbFct]));
% axis tight; axis equal; axis ij;

% w_hat
% subplot(1,3,3); hold on; axis off; title('$\hat{w}$','interpreter','latex','fontsize',20);
% imagesc(reshape(w_hat,nbFct,nbFct));
% axis tight; axis equal; axis ij;

