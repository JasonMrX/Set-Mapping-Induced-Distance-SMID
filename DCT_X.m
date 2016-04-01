function T = DCT_X(Dim)
    N=Dim;
    M=Dim;
    T = zeros(Dim);
    for u=0:Dim-1
        for i=0:Dim-1
            T(u+1,i+1) = cos(pi*u*(2*i+1)/2/N) * Gam(u) * sqrt(2/N) ;
        end
    end
end