function out = plb(fx, stats, x0, lb, ub, opt)
% PLB   --  Optimizes f(x) s.t., lb <= x <= ub
%
% The implementation follows a 'reverse-communication' interface wherein
% the function f(x) and its gradient f'(x) are computed via function
% handles.
%
% Usage:
%       OUT = PLB(FX, x0, LB, UB, OPT)
%
% FX -- function handle to compute f(x) and f'(x)
%       it will be invoked as FX(x)
%
%
% x0 -- Starting vector (useful for warm-starts) -- *can* be zero.
%
% OPT -- This structure contains important options that control how the
% optimization procedure runs. To obtain a default structure the user can
% use 'options = plboptions'. Use 'help plboptions' to get a description of
% what the individual options mean.
%
% OUT contains the solution and other information about the optimization or
% info indicating whether the method succeeded or failed.
%
% See also: options
%
% Version 1.2 (c) 2008, 2010  Dongmin Kim  and Suvrit Sra
% minor modifications: 2014, Suvrit Sra
%
% Reference:
%
% "Tackling box-constrained convex optimization via a new projected quasi-Newton approach"
%    by D. Kim, S. Sra, I. S. Dhillon
%    in SIAM Journal on Scientific Computing (SISC) Oct. 2010;
%
% NOTES: Works well for large-scale problems; Experiment with opt.maxmem



   FX = @(x) fx(x);

   % ------------------------------------------------------
   %  INITIALIZATION
   %  ------------------------------------------------------
   out.iter = 0;
   out.time = 0;
   out.algo = 'PLB';
   out.startTime = cputime;
   out.status = 'Failure';

   %  INITIALIZATION
   %  ------------------------------------------------------
   rho = ones(opt.maxmem, 1);
   alp = rho;
   delx = zeros(length(x0), opt.maxmem);
   delg = delx;
   last = 1;

   out.x = x0;
   out.oldx = x0;
   [out.obj, out.grad] = FX(x0);
   out.oldobj = out.obj;
   out.oldgrad = out.grad;
   out.srch = -out.grad;

   startTime=tic;
   [out.x, flag,step] = pqnLineSearch(out, FX, lb, ub, opt);
   %out.x = pqnFastLS(out, FX, lb, ub, opt);
   lst=toc(startTime);

   [obj out.grad] = FX(out.x);
   out.obj(1) = obj;
   out.iterTime(1) = cputime;
   out.stats(1) = stats(out.x);
   out.statsTime(1) = cputime;

   % -----------------------------------------------------
   %  The main iterative loop
   %  -----------------------------------------------------
   while 1
      out.iter = out.iter + 1;
      commonShowStatus(out, opt);

      gp = find(out.x==0 & out.grad > opt.eps);

      if isfield(opt, 'k_max') & length(gp) > opt.k_max
         [val, ind] = sort(abs(out.grad(gp)), 'descend');
         % TODO(mwytock): invert_tridiag assumes this is in sorted order!!
         gp = sort(gp(ind(1:opt.k_max)));
      end

      % Compute L-BFGS.
      delx(:, last) = out.x - out.oldx;
      delg(:, last) = out.grad - out.oldgrad;
      delx(gp, :) = 0;
      delg(gp, :) = 0;

      out.oldx = out.x;
      out.oldgrad = out.grad;
      out.oldobj = obj;

      out.grad(gp) = 0;
      out.srch = out.grad;
      rho(last) = 1 / (delx(:, last)' * delg(:, last));
      pt = last;

      for i = 1 : min(out.iter, opt.maxmem)
         alp(pt) = rho(pt) * delx(:, pt)' * out.srch;
         out.srch = out.srch - alp(pt) * delg(:, pt);
         pt = opt.maxmem - mod(-pt + 1, opt.maxmem);
      end

      out.srch = 1 / rho(last) / (delg(:, last)' * delg(:, last)) * out.srch;

      for i = 1 : min(out.iter, opt.maxmem)
         pt = mod(pt, opt.maxmem) + 1;
         b = rho(pt) * delg(:, pt)' * out.srch;
         out.srch = out.srch + (alp(pt) - b) * delx(:, pt);
      end

      last = mod(last, opt.maxmem) + 1;

      out.srch = -out.srch;
      out.srch(gp) = 0;

      % ----------------------------------------------
      tmp=tic;
      [out.x, obj, flag, step] = pqnLineSearch(out, FX, lb, ub, opt);
      lst=lst+toc(tmp);

      if (~isempty(opt.grad))
         out.grad = feval(opt.grad, out.x, aux);
      else
         [obj out.grad] = FX(out.x);
      end
      out.obj(out.iter+1) = obj;
      out.iterTime(out.iter+1) = cputime;
      out.stats(out.iter+1) = stats(out.x);
      out.statsTime(out.iter+1) = cputime;

      % termination
      if flag < 0
         term_reason = 6;
      else
         term_reason = commonCheckTermination(out, opt);
      end
      if (term_reason > 0)
         break;
      end

      if cputime-out.startTime > opt.max_time
        break
      end
   end % of while

   % ------------------------------------------------
   %  Final statistics and wrap up
   %  ------------------------------------------------
   out.time = toc(startTime);
   out.lineSearchTime=lst;
   out.status = 'Success';
   %out.term_reason = pqnSetTermReason(term_reason);
end


function st=commonShowStatus(out, opt)
% commonShowStatus -- display running info.

   st='Running';

   if (opt.verbose)
      fprintf('Iter: %d, Obj: %g\n', out.iter, out.obj(out.iter));
   end
end

function flag=commonCheckTermination(out, opt)
% commonCheckTermination -- test various term criteria

   flag=0;

   if out.iter >= opt.maxit
      flag = 1;
      return;
   end

   % elseif out.delta < options.toldel
   %     v = 2;
   %     return;
   % elseif abs(out.tmpo-out.oldo)/out.tmpo < options.tolo
   %     v = 3;
   %     return;
   % else
   %     LCP = [min(options.REG.lambda + out.oldg(:), max(out.x(:), 0)); ...
   %            min(options.REG.lambda - out.oldg(:), max(-out.x(:), 0))];
   %     LCP = norm(LCP(:), inf);
   %     LCP = LCP / max(1e-6, norm(out.x(:), inf));
   %     if LCP < options.tollcp
   %         v = 0;
   %         return;
   %     end
   % end
end

function [x, obj, flag, step] = pqnLineSearch(out, FX, lb, ub, opt)
%  Armijo along projection arc.

   step = 1;
   x = out.x + step * out.srch;
   flag = -1;

   for i = 1 : opt.max_func_evals
      al = find(x < lb);
      x(al) = lb(al);
      au = find(x > ub);
      x(au) = ub(au);
      delx = out.x - x;
      fc = FX(x);

      if out.obj - fc >= opt.sigma * out.grad' * delx
         obj = fc;
         flag = 1;
         return;
      end

      step = step * opt.beta;
      x = out.x + step * out.srch;
   end
   obj = fc;
   x = out.x;
end