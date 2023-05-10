function phi_k_val = phi_k(k,args)
%PHI_K fourier series coefficients of the spatial distribution
%   Function which returns fourier series coefficients of the spatial 
%   distribution passed as gaussian mixture model in args
% Inputs:
%   k   : current index from the set of indices K
%   args: additional arguments. You must pass at least alpha (mixing 
%         coefficients), D (dimension), Am (linear transformation matrices
%         function), Mu (centers of the gaussians), L (period), Sigma 
%         (covariance matrices of the gaussians)
% Outputs: 
%   phi_k_val: length(k)x1 coefficients of the spatial distribution (in
%              Rlength(k))

    if or(size(args.Sigma,3)~=length(args.alpha),...
          length(args.alpha)~=size(args.Mu,2)...
         )
        error("myComponent:invalidArgument",strcat("Error. \nThe numbe",... 
              "r of matrices Sigma must match the number of mixing coe",...
              "fficients alpha and the number of centers Mu (for the g",...
              "aussians)"));
    end

    phi_k_val=0;

    % iterating gaussians one-by-one in the gaussian mixutre model
    for j=1:length(args.alpha) 

        % phi_k_val are real/even, i.e., iterating up to 2^(args.D-1) is
        % sufficient to characterize the signal (/fully)
        for m=1:2^(args.D-1)

            phi_k_val=phi_k_val+(args.alpha(j)/(2^(args.D-1)))*...
                cos(2*pi*k'*args.Am(m)*args.Mu(:,j)/args.L).*...
                exp(diag(-2*pi^2*k'*args.Am(m)*args.Sigma(:,:,j)*...
                    args.Am(m)'*k/args.L^2)...
                   );
        end
    end

    phi_k_val=phi_k_val/(args.L^args.D);
end

