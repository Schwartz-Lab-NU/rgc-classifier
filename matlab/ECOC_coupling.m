function [p_hat,r_hat,n_iter] = ECOC_coupling(M,r,n)
%M is the K-by-L coding matrix
%r is the N-by-L calibrated ECOC response matrix
%n is the L-by-1 vector of the number of training examples
%p_hat is the N-by-K probability estimates
%r_hat is the coupled approximation to r

tol=1e-4;

%From Zadrozny, 2001
[K,L] = size(M);
N = size(r,1);
r = reshape(r,[N,1,L]);
n = shiftdim(n,-2); %size 1-by-1-by-L
mP = shiftdim(M>0,-1); %size 1-by-K-by-L
mN = shiftdim(M<0,-1);

%guess p randomly
p_hat = rand(N,K);

%normalize p <--- p = p./sum(p)
%p_hat = p_hat-min(p_hat,[],'all'); %all positive
p_hat = p_hat./sum(p_hat,2); %sum to 1
p_hat(p_hat<=0) = eps;
old_p = p_hat;

get_r_hat();

delta=inf;
n_iter=0;
while delta>tol
    %update p
    p_hat = p_hat .* ( sum(n.*r.*mP,3) + sum(n.*(1-r).*mN,3) ) ./ ( sum(n.*r_hat.*mP,3) + sum(n.*(1-r_hat).*mN,3) );
    %lots of unnecessary computation here, but probably fast for matlab
    
    %renormalize p
    %p_hat = p_hat-min(p_hat,[],'all');
    p_hat = p_hat./sum(p_hat,2);
    p_hat(p_hat<=0) = eps;
    
    %recompute r_hat
    get_r_hat();
    
    %calculate delta (max squared error)
    delta = max((p_hat-old_p).^2,[],'all');
    old_p = p_hat;
    n_iter = n_iter+1;
end
r_hat = squeeze(r_hat);

    function get_r_hat()
        %compute r_hat = sum(p+) / [sum(p+) + sum(p-)]
       
        r_pos = sum(p_hat .* mP,2); % N-by-1-by-L
        r_neg = sum(p_hat .* mN,2);
        r_hat = r_pos ./ (r_pos+r_neg);
    end

end