function out=runMethod(Y,lam,opts)
% RUNMETHOD --- invokes GFL TV operator
%
%
   switch opts.alg
     case 'pg'
       out=projGradTV(Y,lam,opts);
     case 'og'
       out=optGradTV(Y,lam,opts);
     case 'spg'
       out=spgTV(Y,lam,opts);
     case 'dr'
       out=drTV(Y,lam,opts);
     case 'lb'
       out=lbfgsbTV(Y,lam,opts);
     case 'aspn'
       out=aspnTV(Y,lam,opts);
     case 'pn'
       out=pnTV(Y,lam,opts);
     case 'gflasso'
       out=gflassoTV(Y,lam,opts);
     case 'plb'
       out=plbTV(Y,lam,opts);
     case 'sbb'
       out=sbbTV(Y,lam,opts);
     otherwise
       error('ERROR: Unknown TV algorithm requested!');
   end

end
