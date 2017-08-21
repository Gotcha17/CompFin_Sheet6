clc; clear; //clears the console and all previously stored variables

function Vt = Heston_EuCall_Laplace (St, r, gamma0, kappa, lambda, sigma_tilde, T, t, K, R)
    //Laplace transform for the call is defined
    function y = f_tilde(z)
        y = K^(1-z) / (z*(z-1));
    endfunction
    
    //Conditional characteristic function for the Heston model
    function x = chi_t(u)
        //function d(u) is defined
        d=sqrt(lambda^2+sigma_tilde^2*(%i*u+u^2));
        //gamma is a function of t in the characteristic function of the Heston model,
        //since in our case gamma(0) is given and t=0, simple implementation of gamma(t) is used
        gamma_t=gamma0
        //first part of the chi_t(u) function
        chi_1_t=exp(%i*u*(log(St)+r*(T-t)))*(exp(lambda*(T-t)/2)/(cosh(d*(T-t)/2)+lambda*sinh(d*(T-t)/2)/d))^(2*kappa/sigma_tilde^2);
        //second part of the chi_t(u) function
        chi_2_t=exp(-gamma_t*((%i*u+u^2)*sinh(d*(T-t)/2)/d)/(cosh(d*(T-t)/2)+lambda*sinh(d*(T-t)/2)/d));
        //finally chi_t is "glued together"
        x=chi_1_t*chi_2_t;
    endfunction
    
    //Laplace transform for the call and the conditional characteristic function
    //for the Heston model chi_t
    function integ=integrand(u)
        integ=real(f_tilde(R+%i*u) * chi_t(u-%i*R));
    endfunction
    
    //Computing value of the call option through integration in the range of 
    //0 and 28 (28 is chosen instead of infinityb, ecause the displayed value
    //does not get more precise after this bound)
    Vt=(exp(-r*(T-t))/%pi) * intg(0, 28, integrand);
endfunction

St=100; r=0.05; gamma0=0.2^2; kappa=0.5; lambda=2.5; sigma_tilde=1; T=1; t=0; K=100; R=3;

Vt = Heston_EuCall_Laplace (St, r, gamma0, kappa, lambda, sigma_tilde, T, t, K, R)
disp("Price of the European Call in the Heston model via the Laplace transform approach is "+string(Vt))
