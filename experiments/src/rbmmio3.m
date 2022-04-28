% code for using RBM as a warmstart for MIO3
% q is the number of points for LQS
% m_normal is the number of non-outliers (at the top of the dataset)
% dep_var=true means there is a dependent variable 
%

function rbmmio3(iteration, datafname, q, m_normal, dep_var, resloc, timelimit, lqs_beta)
formulation="rbm-mio3"

lqs_beta

X = readtable(datafname);
X = X{:,:};
[m,n] = size(X); 
if dep_var == true
    X = [X(:,1) ones(m,1) X(:,2:n)]; % add column of 1s
else
    X = [ones(m,1) X]; % add column of 1s at beginning
end
[m,n] = size(X); 

RBMparam = struct();
RBMparam.feaTol=1.0e-6;
RBMparam.minFrac=0.5;
RBMparam.maxRemoveFrac=0.75;
RBMparam.maxDist=-16;
RBMparam.stopCondn='noBetter';
RBMparam.reinstate=1;


[errs,solTime,nHPs] = RBMpcaV1(X, m_normal, q, RBMparam,dep_var)


% get sum of squared error on non outliers
if dep_var == true % get error along response direction
    beta_star = [-1.0; errs.weights(:,3)];
    dist = abs(X*beta_star);
    tot_err = sum(dist(1:m_normal,1).*dist(1:m_normal,1));
    errs.totSqDistTru(1,3) = tot_err;
else
    beta_star = [-errs.RHS(1,3); errs.weights(:,3)];
    beta_star = (beta_star/beta_star(1,1))*(n-1);
    tot_err = errs.totSqDistTru(1,3);
end


[beta_star, f_beta_star] = mio3(iteration, datafname, q, lqs_beta, m_normal, dep_var, formulation, resloc, timelimit, beta_star, errs)

return

end




