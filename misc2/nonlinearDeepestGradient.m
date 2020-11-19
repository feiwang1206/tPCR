function [z, z_iter] = nonlinearDeepestGradient(f,df,z,tol,maxit,cg_method,verbose_flag,reference_norm)

% usage:
%   [z, z_iter] = nonlinearConjugateGradient(f,df,z0,tol,maxit,cg_method,verbose_flag,reference_norm)
%
% Finds the minimum of the function specified by the function handle "f"
% using a nonlinear conjugate gradient algorithm. This algorithm is based on
% the freely available code from Michael Overton "NLCG 1.0".
%
% f = function handle of the cost function
% df = function handle of the derivative of the cost function
% z0 = set a starting point
% tol = stopping tolerance
% maxit = maximum iterations
% verbose = use 0 to switch to silent mode
% cg_method = set the CG-Method
%
% 30.09.2011
% Thimo Hugger
% 21.10.2019
% Fei Wang


persistent h;

if nargin<=3 || isempty(tol)
    tol = 1e-6;
end
if nargin<=5 || isempty(cg_method)
    cg_method = 'fr+pr';
end
if nargin<=6 || isempty(verbose_flag)
    verbose_flag = 2;
end
if nargin<=7 || isempty(reference_norm)
    reference_norm = 1;
end


c1 = 0.05;
backtrack_step = 0.8;

% strongwolfe = 1;
% wolfe1 = 0;
% wolfe2 = 0.5; 
frec = nan; % so defined in case of immediate return
alpharec = nan;

[fz,re] = f(z);
dfz = df(z,re);

% if size(dfz,2) > size(dfz,1) % error return is appropriate here
%     error('gradient must be returned as a column vector, not a row vector');
% end
gnorm = l2norm(dfz);
if fz == inf % better not to generate an error return
    if verbose_flag > 0
        fprintf('nlcg: f is infinite at initial iterate\n');
    end
    return
elseif isnan(fz)
    if verbose_flag > 0
        fprintf('nlcg: f is nan at initial iterate\n');
    end
    return
elseif gnorm < tol
    if verbose_flag > 0
       fprintf('nlcg: tolerance on gradient satisfied at initial iterate\n');
    end
    return
end

p = -dfz;  % start with steepest descent

if nargout==2
    z_iter = cell(1,maxit(1));
end
for iter = 1:maxit(1)
    
    gtp = real(dotprod(dfz,p));
    if  gtp >= 0 || isnan(gtp)
        if verbose_flag > 0
            fprintf('Not descent direction, taking gradient descent at iteration %d, f(z) = %1.5e, gnorm = %1.5e\n', iter, fz, gnorm);
        end
        p = -dfz;
    end
    
    % backtracking line search
    alpha = 1;
    if (fz + c1*alpha*gtp) > 0
        fy = f(z,alpha,p,1);
        counter = 1;
        while ( (fy > (fz + c1*alpha*gtp)) ) % sufficient decrease condition aka first Wolfe condition or Armijo rule, second Wolfe condition is too expensive
            alpha = backtrack_step * alpha;
%             fy_old = fy;
            fy = f(z,alpha,p,0);
            counter = counter + 1;
            if counter>=200
                break;
                alpha = 0;
            end
        end
        
    else % rare case
        fy_old = Inf;
        fy = f(z,alpha,p,1);
        counter = 1;
        while fy < fy_old
            alpha = backtrack_step * alpha;
            fy_old = fy;
            fy = f(z,alpha,p,0);
            counter = counter + 1;
            if counter>=200
                break;
                alpha = 0;
            end
        end
        alpha = alpha / backtrack_step;
        
    end
    
    z = z + alpha*p;
    
    [fz,re] = f(z);
    dfz = df(z,re);
    
    
% original line search method by Michael Overton (backtracking seams to work as well and is faster)
%     if strongwolfe % strong Wolfe line search is usually essential (default)
%         [alpha, z, fz, g, fail] = ...
%             linesch_sw(z, fz, g, p, pars, wolfe1, wolfe2, fvalquit, verbose_flag);
%     else  % still, leave weak Wolfe as an option for comparison
%         [alpha, z, fz, g, fail] = ...
%             linesch_ww(z, fz, g, p, pars, wolfe1, wolfe2, fvalquit, verbose_flag);
%     end


    gnorm = l2norm(dfz);
    z_iter{3}(iter)=alpha*gnorm/reference_norm;
    z_iter{7}(iter) = fz;
    alpharec(iter) = alpha;
    
    if verbose_flag==1
        if isempty(h) || ~ishandle(h) || ~strcmp(get(h,'Tag'),'cg_figure')
            h = figure('Tag','cg_figure');
        end
        if gcf~=h
            set(0,'CurrentFigure',h);
        end
        if isvector(z)
            plot(abs(z));
        elseif is2Darray(z)
            imagesc(abs(z));
            colormap gray;colorbar
        elseif length(size(z))==3
            imagesc(abs(array2mosaic(z)));
            colormap gray;colorbar
        end
        axis off;
        title(sprintf('iter %d: step = %1.5e, f(z) = %1.5e, gnorm = %1.5e\n', iter, alpha*gnorm/reference_norm, fz, gnorm));
        drawnow;
%         fprintf('iter %d: step = %1.5e, f(z) = %1.5e, gnorm = %1.5e\n', iter, alpha*gnorm/reference_norm, fz, gnorm);
    elseif verbose_flag==2
        fprintf('iter %d: step = %1.5e, f(z) = %1.5e, gnorm = %1.5e\n', iter, alpha*gnorm/reference_norm, fz, gnorm);
    end
    
    for ii=1:length(maxit)
        if iter==maxit(ii)
            z_keep{1}{ii}=z;
            z_keep{2}(ii)=maxit(ii);
        end
    end
    if iter>2 && (z_iter{3}(iter-1)+z_iter{3}(iter))/2 < tol
        if verbose_flag > 0
            fprintf('step length below tolerance, quit at iteration %d, f(z) = %1.5e\n', iter, fz);
        end
        break
    end
    
    
    p =  - dfz;

    
end % for loop
if iter<maxit
    z_keep{1}{1}=z;
    z_keep{2}(1)=iter;
end
z_iter{1} = z_keep;

if verbose_flag > 0
    fprintf('%d iterations reached, f(z) = %1.5e, gnorm = %1.5e\n', maxit, fz, gnorm);
end