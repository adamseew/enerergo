function f_k_x_val = f_k_x(k,x,args) 
%F_K_X fourier series coefficients of the basis functions with complex 
%      exponentials for a spatial signal x
%   Fuction which returns the fourier series coefficients of the basis 
%   functions with complex exponentials for a spatial signal x (a points) 
%   on the interval [-L/2,L/2], where L is the period
% Inputs:
%   k   : current index from the set of indices K
%   x   : point as a row vector (works with CasADi variables)
%   args: additional arguments. You must pass at least D (dimension), and 
%         L (period)
% Outputs:
%   f_k_x_val: length(k)x1 fourier series coefficients for a spatial 
%              signal x

    f_k_x_val=ones(length(k),1)/args.L^args.D;

    for d=1:args.D
            % ?                .*
            f_k_x_val=f_k_x_val.*cos((2*pi*x(d)*k(d,:))/args.L)';%-...  
            %                 1i*sin((2*pi*x(d)*k(d,:))/args.L)';
            % ?                 .*
            % ?               '.*
            %f_k_x_val=f_k_x_val.*...
            %    (-2*pi*k(d,:)'.*sin((2*pi*x(d)*k(d,:))/args.L)'/args.L);%-...  
            %   2*1i*pi*k(d,:)'.*cos((2*pi*x(d)*k(d,:))/args.L)'/args.L);
    end
end

