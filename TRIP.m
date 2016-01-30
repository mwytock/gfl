%out=TRIP(fgx,reg,opt)        Trust Region Proximal solver
%
% General proximal solver for an optimization problem in the form
%
%   min_x  f(x)  +  r(x)
%
% where f is smooth and convex and r is non-smooth convex.
% TRIP finds the solution by means of alternating gradient steps
% for f and proximity steps for r. 
% 
% The length of the steps is determined by the spectral
% formulae of Barzilai-Borwein, which is non-monotonous.
% This generally ensures good performance in practice but can lead 
% to non-convergent behavior in some setings. If such
% situation is detected, TRIP performs a monotonic step
% following the minimun norm subgradient. Such subgradient is
% computed efficiently through the proximity operator.
%
% Inputs:
%
%   - fgx:  computes FX and GFX 
%   - reg:  computes R(X) and if needed, prox of R(X)
%   - opt: optimization options.
%
% Outputs:
%
%   - out: solution structure with fields: 
%
% (c) 2010, Dongmin Kim and Suvrit Sra
% (c) 2013, adapted for TVprox by Alvaro Barbero.
% (c) 2014, simplified for new TV by Suvrit Sra
% 
%

function out=TRIP(FX, REG, x0, opt)
   out.name = 'TRIP';
   out.version = 'Version 3: March 2014';
    
    % Check input parameters
    if ~exist('opt','var'), opt = options('trip'); end

    out.startTime = tic;
    out.obj = nan*ones(opt.maxit,1);
    out.iterTime = nan*ones(opt.maxit,1);
    out.iter = 0;
    out.term_reason=0;
    out.delta = 1;

    % TRIP initialization
    [x out.oldx out.g out.oldg out.o out.oldo] = initTrip(FX, REG, x0);

    %% Misc. init.
    % safe bet at the first iter.
    alpha = opt.alpmax;
    lazy = 1;         % just start the BB dance,
    refx = x;
    refo = out.o;
    refg = out.g;
    out.term_reason = 0;
    out.tmpo = refo;
    out.delta = opt.delmax;

    % -------------------------------------------------------------------------
    %  The main iterative loop
    % -------------------------------------------------------------------------
    while 1
        out.iter = out.iter + 1;
        
        % ---------------------------------------------------------------------
        % Internal Solver: TRIP proximity step
        % ---------------------------------------------------------------------
        [newx s model_diff newo newg term_reason] = solver(FX, REG, x, out.g, out.delta, alpha, opt, lazy);
        
        
        % ---------------------------------------------------------------------
        % check termination
        % ---------------------------------------------------------------------
        term_reason = commonCheckTermination(out, opt);
        if term_reason > 0
           out.X=newx;
            break;
        end

        % Gather statistics
        out.oldo = out.o;
        out.o = newo;
        out.tmpo = out.o;
        out.obj(out.iter)  = out.o;
        out.iterTime(out.iter) = toc(out.startTime);
        commonShowStatus(out, opt);

        pass = 0;
        
        % Check if we can avoid updating the Trust-Region.
        % We can if we are still in the *null* iterations or if objective function has improved
        if lazy < opt.maxnull || out.o < refo
            % Keep lazy until finishing null iterations
            if lazy < opt.maxnull, lazy = lazy + 1;
            % In not-null iterations update reference objective
            % TODO: should we update reference as well in null iterations if objective improved?
            else
                lazy = 1;
                refx = newx;
                refo = newo;
                refg = newg;
            end
            
            % Recompute BB step and go for a new iteration
            out.oldx = x;
            x = newx;
            out.oldg = out.g;
            out.g = newg;
            alpha = update_alpha(opt, s, out.oldg, out.g, out.iter);
            continue;
        end
        
        %
        % Bad case: should go through all the nasty updates.
        %
        
        % Reset the current iter. to the reference
        x = refx;
        out.o = refo;
        out.g = refg;
        %alpha = 1; %TODO: necessary?
        
        % Enforced monotonicity in TRIP proximity step
        lazy = 0;
        [newx s model_diff newo newg out.term_reason] = ...
            solver(L, reg, x, out.g, out.delta, alpha, opt, lazy);

        % measure the model accuracy
        rho = (out.o - newo) / model_diff;

        % ---------------------------------------------------------------------
        % Update trust-region
        % ---------------------------------------------------------------------
        update = 0;
        if rho >= opt.eta2
            out.delta = min(opt.delmax, opt.gam3 * out.delta);
            update = 1;
        elseif rho >= opt.eta1
            out.delta = opt.gam2 * out.delta;
            update = 1;
        else
            out.delta = opt.gam1 * out.delta;
            %fprintf('        delta=%f iter=%d\n', out.delta, out.iter);
        end

        % Update BB stepsize if needed
        if update > 0
            out.oldx = x;
            x = newx;
            out.oldg = out.g;
            out.g = newg;
            alpha = update_alpha(opt, s, out.oldg, out.g, out.iter);

            refx = newx;
            refo = newo;
            refg = newg;
        end
            
    end % of while

    % -------------------------------------------------------------------------
    %  Final statistics and wrap up
    % -------------------------------------------------------------------------
    out.time = toc(out.startTime);
    out.status = 'Success';
    out.term_reason = tripTermReason(term_reason);
    out.tmpo = out.oldo;
    out.X=newx;
end

%
%  Internal solver: performs TRIP proximity step
%
function [newx s model_diff newo newg term_reason] = solver(FX, REG, x, g, delta, alpha, opt, lazy)
    term_reason = 0;        % no reason for term.

    % Compute displacement within trust region s
    %fprintf(1,'alpha=%g\n',alpha); %FIXME
    s = -g / alpha;
    % TRIP proximity step
    newx = get_newx(REG, x, s, delta, alpha, opt);
    % Displacement by proximity
    s = newx - x;
    %fprintf(1,'||s||=%g\n',norm(s)); %FIXME
    %fprintf(1,'||x||=%g\n',norm(x)); %FIXME
    

    % New objective values
    [newo newg] = FX.fgx(newx);
    newr = REG.fval(newx);
    newo = newo + newr;

    % If in null iterations, we are done
    if lazy > 0
        model_diff = 0;
        return;
    end
    
    %
    % Else we need to enforce monotonicity in the proximal step
    %
    % TODO: IMPL THIS PHASE!!
    % Compute threshold of the descent condition: norm of the minimum norm subgradient
    threshold = thresholding(FX, REG, x);
    % Account for the reduction and trust-region size parameters
    threshold = threshold * opt.beta * min(opt.delmin, delta);
    r = REG.fval(x);
    %fprintf(1,'r=%g\n',r); %FIXME
    %fprintf(1,'x=[ '); fprintf(1,'%f ',x); fprintf(1,']\n'); %FIXME
    %fprintf(1,'s=[ '); fprintf(1,'%f ',s); fprintf(1,']\n'); %FIXME
    %fprintf(1,'g=[ '); fprintf(1,'%f ',g); fprintf(1,']\n'); %FIXME
    model_diff = r - (s(:)' * g(:)) - .5 * (s(:)' * s(:)) * alpha - newr;

    % -------------------------------------------------------------------------
    %  Line search until improvement above threshold
    % -------------------------------------------------------------------------
    iiter = 0;

    fprintf(1,'model_diff=%g, threshold=%g\n',model_diff,threshold); %FIXME
    while model_diff <= threshold
        fprintf(1,'THRESHOLDING LOOP'); %FIXME
        if iiter > opt.maxls
            term_reason = 1;
            break;
        elseif model_diff <= 0
            term_reason = 2;
            break;
        end
        
        % Compute proximity scaling following linear search
        alpha = opt.tau * alpha;
        s = -g / alpha;
        
        fprintf(1,'alpha=%g\n',alpha); %FIXME

        % Proximity step with scaling
        newx = get_newx(REG, x, lambda, s, delta, alpha, opt);
        s = newx - x;

        % Current improvement
        model_diff = r - (s(:)' * g(:)) - .5 * (s(:)' * s(:)) * alpha - ...
            reg.val(newx); 
        iiter = iiter + 1;
    end

    % New objective values
    [newo newg] = FX.fgx(newx);
    newo = newo + REG.fval(newx);
end

% 
%  Our friend, BB
%
function alpha = update_alpha(opt, s, oldg, g, iter)
    v = g(:) - oldg(:);

    if mod(iter, 2) == 0
        alpha = (s(:)' * v(:)) / (s(:)' * s(:));
    else
        alpha = (v(:)' * v(:)) / (s(:)' * v(:));
    end

    % Bound alpha
    alpha = min(opt.alpmax, max(opt.alpmin, alpha));
end

%
% TRIP proximity step
%
function newx = get_newx(REG, x, s, delta, alpha, opt)
    scale = 1/alpha;
    newx = REG.prox(x+s,scale);

    % -------------------------------------------------------------------------
    %  delta inf-norm projection
    % -------------------------------------------------------------------------

    upper = x + delta;
    lower = x - delta;

    uidx = find(newx > upper);
    lidx = find(newx < lower);

    newx(uidx) = upper(uidx);
    newx(lidx) = lower(lidx);
end

%
% TRIP reasons for termination
%
function r=tripTermReason(tr)
    switch tr
      case 1
        r = 'Search for step-size fails';
      case 2
        r = 'Fail for sufficient reduction in the model';
      case 3
        r = 'Maximum number of iterations reached';
      case 4
        r = 'Too small trust-region';
      case 5
        r = 'Absolute change in obj below threshold';
      case 6 
        r = 'Target LCP achieved';
      otherwise
        r = 'Undefined';
    end
end

%
% TRIP initialization: find two initial points
%
function [x oldx g oldg o oldo] = initTrip(FX, REG, x0)
    % Find another random initial point.
    %
    x1             = ones(size(x0))/1e+2;
    [c0 g0]        = FX.fgx(x0);
    [c1 g1]        = FX.fgx(x1);
    r0             = REG.fval(x0);
    r1             = REG.fval(x1);

    o0 = c0 + r0;
    o1 = c1 + r1;

    if o0 < o1
        x         = x0;
        oldx      = x1;
        g         = g0;
        oldg      = g1;
        c         = c0;
        oldc      = c1;
        r         = r0;
        oldr      = r1;
        o         = o0;
        oldo      = o1;
    else
        x         = x1;
        oldx      = x0;
        g         = g1;
        oldg      = g0;
        c          = c1;
        oldc       = c0;
        r         = r1;
        oldr      = r0;
        o         = o1;
        oldo      = o0;
    end
end

%
% commonCheckTermination -- test various term criteria
%
function flag=commonCheckTermination(out, opt)
    flag=0;

    if out.iter >= opt.maxit
        flag = 1;
        return;
    elseif ( abs(out.tmpo-out.oldo) < opt.tolo || abs(out.tmpo-out.oldo)/out.oldo < opt.tolro) && out.iter > 3
        flag = 3;
        return;
    elseif out.iter > 1 && out.iterTime(out.iter-1) > opt.maxt
        flag = 7;
        return;
    end
end

%
% commonShowStatus -- display running info.
%
function st=commonShowStatus(out, opt)
    st='Running';

    if (opt.verbose)
        fprintf('Iter: %d, Obj: %g\n', out.iter, out.obj(out.iter));
    end
end

%
% thresholding operator: TODO
%
function th = thresholding(FGX, reg, x) 
   disp('THRESHOLDING INVOKED');
   th=x;
end
 