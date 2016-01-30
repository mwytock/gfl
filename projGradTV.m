function out=projGradTV(Y,lam,opts)
% PROJGRADTV -- run projected gradient for GFL-TV op
  t_start = cputime;
  UET = @(X) rightMultByEtrans(X);
  U = initU(Y,lam,opts);               % to allow unified initialization
  % TODO(mwytock): Fix to take lambda as vector
  lam = lam(1);
  [d,n]=size(Y);

  for i=1:opts.maxit
    T = UET(U)-Y;                % X = -T, i.e., X = Y-UE'
    XE=T(:,2:n)-T(:,1:n-1);      % X*E
    U = projInfTwoBall(lam,U+0.25*XE);
    out.iterTime(i) = cputime;
    [out.obj(i) out.gap(i)] = primal_objval2(Y,U,lam);
    out.statsTime(i) = cputime;
    if cputime-t_start > opts.max_time
      break
    end
  end
  out.U=U;
  out.time = out.statsTime - cumsum(out.statsTime-out.iterTime) - t_start;
  out.gap(1) = inf;
end
