function [dx, Voc] = thevenin(x, I, args)

    if size(x,1)~=3
        error("myComponent:invalidArgument",strcat("Error. \nThe state",...
              " must be a row vector of three values. v1, v2, and z (s",...
              "tate of charge)"));
    end
    
    Voc=args.V+x(1)+x(2)+I*args.Rs;
    dx=[-1/(args.R1*args.C1) 0 0; 0 -1/(args.R2*args.C2) 0; 0 0 0]*x+...
        [1/args.C1;1/args.C2;-1/args.Q]*I;
end
