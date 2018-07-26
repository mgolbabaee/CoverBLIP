function res = algo_IPA(data, param)
% Iterative projection (matched-filtering) algorithm. Depending on D.match
% this algorithm performs exact/inexact matched filtering. i.e. BLIP/CoverBLIP
% Inputs: data
%                 y: k-space data
%                 D: dictionary/matching structure
%                 F: forward/adjoint operator
%                 N,M,L: mrf data size; N*M spatial resolution, L temporal resolution
% Inputs: param
%                 tol: progress tolerence for stoppage 
%                 max iter: maximum iteration
%                 epsilon: approximation level (in case of inexact iteration)
%
% Outputs: res
%                 dm: indices of chosen fingerprints (map)
%                 pd: mproton density map
%                 err: meaasurement fidelity error
%                 comp_cost: total projection cost
%                 time_proj_cost: total projection time
%                 iter: number of performed iterations before stoppage
%
% (c) Mohammad Golbabaee, 2017
%%
N = data.N; M = data.M; L = data.L;
y = data.y; F = data.F; D = data.D;
% Vk = data.Vk;
% Vkt = data.Vkt;

proposed_step = param.step; tol = param.tol;
max_iter = param.iter; 
if ~isfield(param,'epsilon'), epsilon=0.2; else epsilon = param.epsilon;end
if ~isfield(param,'current_est'), current_est=0; else current_est = param.current_est;end

% --initialization
if param.svd==1, x = zeros(N,M,param.k); else x = zeros(N,M,L);end
Fx = zeros(size(y));

% Initial residual
err = -y;
err_norm_old = norm(err(:), 2);
%initial gradient
grad = F.adjoint(err);

% Indices of D coefficient
ind = zeros(1,N*M);
% pseudo-Proton density
pseudo_density = zeros(1, N*M);

comp_cost = zeros(1, max_iter);
comp_svd_cost = zeros(1, max_iter);
time_proj_cost = zeros(1, max_iter);
err_iter = zeros(1, max_iter);
% % sol_err = nan(1,max_iter);

stepsize_tol = 1;%0.99;%9999;
run_iter = 0;
dumm=0;
for iter = 1:max_iter
    
    oldx = x;
    oldFx = Fx;
    step = proposed_step;
    done = 0;
    
    while ~done
        % Gradient step
        x = x - step*grad;

        x = transpose(reshape(x,N*M,[]));
       
        ts = tic;
        [pseudo_density, ind, x, comp] = D.match(x, epsilon, ind);
        time_proj_cost(iter) = time_proj_cost(iter)+ toc(ts);

        comp_cost(iter) = comp_cost(iter) + comp;
        
        % Taking svd dim reduction cost to account
        if param.svd==1, comp_svd_cost(iter) = comp_svd_cost(iter) + 2*param.k*L*M*N ;end
        
        % scale up x
        x = x.*repmat(pseudo_density,[size(x,1),1]);
        x = reshape(transpose(x), N, M, []);
        
        Fx = F.forward(x);
        
        %initial rescaling
        if dumm == 0
            sc = norm(y(:))/norm(Fx(:));
            x = x*sc;
            Fx = Fx*sc;
            dumm=1;
        end
        
        omega = (norm(x(:)-oldx(:),2)/norm(Fx(:)-oldFx(:),2))^2;
        
        
        if (step > omega*stepsize_tol) && (max_iter>1)
            disp('reducing step size');
            step = step/2;
            x = oldx;
        else
            err = Fx - y;
            err_norm = norm(err(:), 2)^2;
            grad = F.adjoint(err);
            err_iter(iter+1) = err_norm;
            done=1;
        end
        run_iter = run_iter + 1;
    end
    
    if abs(err_norm_old - err_norm) < tol * err_norm
        fprintf('Status: BLIP convergenced\n');
        break;
    end
    err_norm_old = err_norm;
    
    % Verbose
    fprintf('Iter %i, obj= %e, |y-F.forward(x)|/|y|= %e\n',...
        iter, err_norm, sqrt(err_norm)/norm(y(:)) );
    
end

ind = reshape(ind, N, M);
density = reshape(pseudo_density, N, M)./D.scale(ind);

x=reshape(x,N*M,[]);
x=transpose( data.Vkt(transpose(x)) );

res.dm = ind;
res.pd = density;
res.err = err_iter;
res.comp_cost = comp_cost;
res.time_proj_cost = time_proj_cost;
res.iter = run_iter;
res.comp_svd_cost = comp_svd_cost;
end