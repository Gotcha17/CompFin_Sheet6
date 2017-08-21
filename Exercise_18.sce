clc; clear; //clears the console and all previously stored variables
function V0 = BS_EuCall_FFT (S0, r, sigma, T, K, R, N, M)
    //defining initial parameters
    delta=M/N;
    m=1:N;
    kappa_1=0;
    
    //f_tilde_0 for European Call for K=e^kappa
    function y = f_tilde(z)
        y = 1/ (z*(z-1));       //since kappa=0
    endfunction
    
    //Characteristic function for the BS model
    function x = chi(u)
        x=exp(%i*u*(log(S0)+r*T)-(%i*u+u^2)*sigma^2/2*T);
    endfunction
    
    //defining the g function
    function y = g(u)
        y=f_tilde(R+%i*u) * chi(u-%i*R);
    endfunction
    
    //computing x_n vector for n=1,...,N
    for n=m
        x_n(n)=g((n-1/2)*delta)*delta*exp(-%i*(n-1)*delta*kappa_1)
        kappa_m(n)=kappa_1+(n-1)*2*%pi/M
    end
    
    //applying Fast Fourier Transform to compute x_hat_m from x_n for m=1,...,N
    x_hat_all=fft(x_n)
    
    //defining function for computation of option Values for different kappas (4.16)
    function y = V0_k(x)
        y = exp(-r*T+(1-R).*kappa_m)/%pi.*real(x.*exp(-%i/2*delta.*kappa_m))
    endfunction
    
    //initializing matrix with dimension of Nx2, to store all results
    //computed for option value and strike
    V0_km=zeros(N,2)
    //Storing option price values for a range of strike prices (kappas)
    V0_km(:,2)=V0_k(x_hat_all)
    //Storing strike prices 
    V0_km(:,1)=exp(kappa_m)
    //matrix transpose for later use in the interpln function (needed in this form)
    V0_km=V0_km'
    //linear interpolation of the results, compute option prices for needed strike prices
    V0=interpln(V0_km, K)
endfunction

S0=100; r=0.05; sigma=0.2; T=1; K=80:130; R=1.1; N=2^11; M=50;

V0 = BS_EuCall_FFT (S0, r, sigma, T, K, R, N, M)
disp("V0 = " + string(V0) + " K = " + string(K))
