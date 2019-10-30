function [xopt,fopt,min_val_initial] = epi_sph_grad_descent(x,p1,p2,ptsnum,diff_amplifier,tol,maxiter,alpha,beta)
if length(x) ~= 4
    error('Incorrect number of input arguments.')
end
fun=@(x)calc_phi3(x,p1,p2,ptsnum,diff_amplifier);
% termination tolerance
% tol = 1e-6;
% maximum number of allowed iterations
% maxiter = 100000;
% step size ( 0.33 causes instability, 0.2 quite accurate)
% alpha = 0.001;
% step direction;
% angle_div=10;
% beta = linspace(0,2*pi*(1-1/angle_div),angle_div);
% initialize gradient norm, optimization vector, iteration counter, perturbation
length_beta=length(beta);
niter = 0;
min_val_initial=fun(x);
min_val_old=min_val_initial;

flag=1;
x_candis=cell(length_beta,length_beta);
f_candis=zeros(length_beta,length_beta);
while flag
    e1s=repmat(x(1:2)',1,length_beta)+alpha.*[cos(beta);sin(beta)];
    e2s=repmat(x(3:4)',1,length_beta)+alpha.*[cos(beta);sin(beta)];
    for m=1:length_beta
        for n=1:length_beta
            x_candis{m,n}=[e1s(:,m);e2s(:,n)]';
            f_candis(m,n)=fun(x_candis{m,n});
        end
    end
    min_val=min(min(f_candis));
    if min_val < min_val_old-tol
        [min_indm,min_indn]=find(f_candis==min_val);
        x=x_candis{min_indm,min_indn};
        min_val_old=min_val;
    else
        if min_val<min_val_old && min_val >= min_val_old-tol
            x=x_candis{min_indm,min_indn};
            flag = 0;
            disp('min_val change smaller than tol.');
        else 
            flag = 0;
            disp('min_val larger than min_val_old.');
        end
    end
    niter=niter+1;
    if niter >= maxiter
        flag = 0;
        disp('niter reaches maxiter.');
    end
end
xopt=x;
fopt=min_val;
            
        