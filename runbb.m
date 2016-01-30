function out = runbb(GFX, stats, x0, opt)
% runBB   -- Solve BB based nonneg problem
%
% NOTE --- guaranteed convergence phase: *REMOVED*
%
%
% x0 -- Starting vector (useful for warm-starts).
%
% OPT -- This structure contains important opt that control how the
% optimization procedure runs. To obtain a default structure the user can
% use 'opt = solopt'. Use 'help solopt' to get a description of
% what the individual opt mean.
%
% Most important options to tune as: opt.tolg, opt.maxit
%
%
% OUT contains the solution and other information about the optimization or
% info indicating whether the method succeeded or failed.
%
% See also: solopt, pbbnnls
%
% Version 1.0 (c) 2010 Suvrit Sra
% Version 1.1 (c) 2014 Suvrit Sra
%

   out.iter = 0;
   out.startTime = cputime;
   out.status = 'Failure';

   out.x      = x0;
   [x out.oldx out.grad out.oldg out.o out.oldo] = initBB(GFX, x0);


   %% Begin the main algorithm
   fprintf('Running: **** SBB- ****\n\n');
   fprintf('Iter   \t     Obj\t\t  ||pg||_inf\t\t ||x-x*||\n');
   fprintf('-------------------------------------------------------\n');

   % There is some ugly code-duplication below for efficiency reasons...
   while 1
      out.iter = out.iter + 1;
      [termReason, out.pgTimes(out.iter)] = checkTermination(opt, out);
      if (termReason > 0), break; end
      if cputime-out.startTime > opt.max_time
        break
      end
      xd=out.x-out.oldx;
      gd=out.grad-out.oldg;
      step = computeBBStep(out,xd,gd);
      out.oldx=out.x;
      out.x = out.x - step * out.grad;
      out.oldg = out.grad;
      out.x(out.x < 1e-8) = 0;
      [out.dd_obj(out.iter) out.grad] = GFX(out.x);
      out.iterTime(out.iter) = cputime;
      out.obj(out.iter) = stats(out.x);
      out.statsTime(out.iter) = cputime;
      if (opt.truex), out.trueError(out.iter) = norm(opt.xt-out.x); end
      if (opt.verbose)
         fprintf('%04d\t %E\t%E\t%E\n', out.iter, out.dd_obj(out.iter), out.pgTimes(out.iter), nan);
      end
   end % of while

   %%  Final statistics and wrap up
   out.status = 'Success';
   out.termReason = setTermReason(termReason);
end

% Compute BB step; for SBB also modifies out.oldg
function step = computeBBStep(out,xd,gd)
   gp = find(out.x==0 & out.grad > 1e-8);

   xd(gp)=0;
   gd(gp)=0;

   % update x & gradient
  if (mod(out.iter, 2) == 0)
     nr=xd'*xd;
     if (nr==0), step=0; return; end
     dr=xd'*gd;
     step = nr/dr;
  else
     nr=xd'*gd;
     if (nr==0), step=0; return; end
     dr=gd'*gd;
     step = nr/dr;
  end


  %if (step <= 0 || step > 1e8), step=1; end

end

% check various termination criteria; return norm of pg
function [v pg] = checkTermination(options, out)
   % pgnorm limit -- need to check this first of all
   gp = find( (out.x ~= 0 | out.grad < 0));
   pg = norm(out.grad(gp), 'inf');
   if (pg < options.tolg), v=8; return; end

   if (out.iter >= options.maxit)
      v = 4;
      return;
   end

   % All is ok...
   v = 0;
end

%% Prints status
function showStatus(out, options)
   if (options.verbose)
      fprintf('.');
      if (mod(out.iter, 30) == 0)
         fprintf('\n');
      end
   end
end

% String representation of termination
function r = setTermReason(t)
   switch t
     case 1
       r = 'Exceeded time limit';
     case 2
       r = 'Relative change in x small enough';
     case 3
       r = 'Relative change in objvalue small enough';
     case 4
       r = 'Maximum number of iterations reached';
     case 5
       r = '|x_t+1 - x_t|=0 or |grad_t+1 - grad_t| < 1e-9';
     case 6
       r = 'Line search faild';
     case 7
       r = '|x^T * grad| < opt.pbb_gradient_norm';
     case 8
       r = '|| grad ||_inf < opt.tolg';
     case 100
       r = 'The active set converged';
     otherwise
       r = 'Undefined';
   end
end


%
% TRIP initialization: find two initial points
%
function [x oldx g oldg o oldo] = initBB(FX, x0)
    % Find another random initial point.
    %
    x1             = ones(size(x0))/1e+2;
    [c0 g0]        = FX(x0);
    [c1 g1]        = FX(x1);

    o0 = c0;
    o1 = c1;

    if o0 < o1
        x         = x0;
        oldx      = x1;
        g         = g0;
        oldg      = g1;
        c         = c0;
        oldc      = c1;
        o         = o0;
        oldo      = o1;
    else
        x         = x1;
        oldx      = x0;
        g         = g1;
        oldg      = g0;
        c          = c1;
        oldc       = c0;
        o         = o1;
        oldo      = o0;
    end
end
