function [varargout] = odeparareal (ode_c, solver_c, opt_c, ode_f, solver_f, opt_f, trange, y0, varargin)
%ODEPARAREAL Solves an evolution problem using Parareal.
%   SOL = ODEPARAREAL(ode_c, solver_c, opt_c, ode_f, solver_f, opt_f, trange, y0, varargin) 
%   returns a structure array SOL with fields:
%   SOL.t is fine discretization in time, SOL.y is converged fine solution,
%   SOL.err, SOL.f, SOL.c are error, fine and coarse solutions at each iteration
%   [t,y] = ODEPARAREAL(ode_c, solver_c, opt_c, ode_f, solver_f, opt_f, trange, y0, varargin) 
%   returns the fine discretization in time t, and the converged fine solution y.
%
%   Parallelization requires Octave and the parallel package 
%   https://octave.sourceforge.io/parallel/

% The source code is distributed under the terms of the GNU General
% Public License (GPL) (version 2 or later).

  if (nargin < 8)
    print_usage ();
  end

  if (nargin >= 9)
    if (~isstruct (varargin{1}))
      %# varargin{1:len} are parameters for fun
      odeopts = odeparaset ();
      funarguments = varargin;
    elseif (numel (varargin) > 1)
      %# varargin{1} is an ODE odeopts structure opt
      odeopts = varargin{1};
      funarguments = {varargin{2:numel (varargin)}};
    else  %# if (isstruct (varargin{1}))
      odeopts = varargin{1};
      funarguments = {};
    end
  else  %# nargin == 8
    odeopts = odeparaset ();
    funarguments = {};
  end

  %# Check time span vector
  if (~isnumeric (trange) || ~isvector (trange))
    error ('odeparareal: TRANGE must be a numeric vector');
  end

  if (numel (trange) < 2)
    error ('odeparareal: TRANGE must contain at least 2 elements');
  elseif (trange(1) == trange(end))
    error ('odeparareal: invalid time span, TRANGE(1) == TRANGE(2)');
  else
    odeopts.direction = sign (trange(end) - trange(1));
  end
  trange = trange(:);

  %# FIXME: make it work when integrating backwards in time
  if (odeopts.direction==-1)
    error ('odeparareal: invalid time span, TRANGE(1) > TRANGE(2)');
  end

  %# Check initial condition vector
  if (~isnumeric (y0) || ~isvector (y0))
    error ('odeparareal: Y0 must be a numeric vector');
  end
  y0 = y0(:);

  %# Check coarse solver input
  if ~( isa (ode_c, 'function_handle') || exist (ode_c, 'file'))
    error ('odeparareal: ODE_C must be a valid function handle');
  end

  if (ischar (solver_c))
    try
      solver_c = str2func (solver_c);
    catch
      %warning (lasterr);
    end
  end
  if (~isa (solver_c, 'function_handle'))
    error ('odeparareal: SOLVER_C must be a function handle or file');
  end

  if (isempty (opt_c))
    opt_c = odeset ();
  end

  %# Check fine solver input
  if ~( isa (ode_f, 'function_handle') || exist (ode_f, 'file'))
    error ('odeparareal: ODE_F must be a function handle or file');
  end

  if (ischar (solver_f))
    try
      solver_f = str2func (solver_f);
    catch
      %warning (lasterr);
    end
  end
  if (~isa (solver_f, 'function_handle'))
    error ('odeparareal: SOLVER_F must be a valid function handle');
  end

  if (isempty (opt_f))
    opt_f = odeset ();
  end

  odeopts.funarguments = funarguments;
  
  %# If no ouput arguments use odeplot
  if (isempty (odeopts.OutputFcn) && nargout == 0)
    odeopts.OutputFcn = @odeplot;
    odeopts.haveoutputfunction = true;
  else
    odeopts.haveoutputfunction = ~isempty (odeopts.OutputFcn);
  end
  
  % TODO: add Prolongator
  % isempty (odeopts.Prolongator)

  if (strcmp (odeopts.Verbose, 'on'))
    verbosity = 1;
    tim_p = tic;
  else
    verbosity = 0;
  end

  %# Set time range subdivision according to processors
  %# use nonequidistant partitioning if given in trange
  idx = [(1:odeopts.NProc);(2:odeopts.NProc+1)]';
  if (length(trange) == 2)
    trange  = linspace (trange(1), trange(end), odeopts.NProc+1);
  elseif (length(trange) ~= odeopts.NProc+1)
    error ('odeparareal: TRANGE and odeopts.NProc do not match');
  end
  trange_proc  = trange(idx);

  if (isempty (odeopts.NIter))
    odeopts.NIter = odeopts.NProc;
  end

  %# Initialize the OutputFcn
  if (odeopts.haveoutputfunction)
    if (~isempty (odeopts.OutputSel))
      init_vals = y0(odeopts.OutputSel,end);
    else
      init_vals = y0;
    end
    feval (odeopts.OutputFcn, trange, init_vals, 'init', odeopts.funarguments{:});
    %end
  end

  %# check for parallel package
  try
    pkg load parallel
    has_parallel = true;
  catch
    has_parallel = false;
    warning('parareal: sequential time stepping since parallel package not available.');
  end

  %# initialize solutions
  err          = zeros (odeopts.NIter, 1);              %# error per iteration
  y_cold       = 0;                                     %# old solution from coarse propagator
  y_cnew       = zeros (odeopts.NProc+1, length(y0));   %# matrix containing the new coarse solution
  y_f0         = zeros (odeopts.NProc+1, length(y0));   %# matrix containing the fine solution for proc i+1 at t_0
  y_fend       = zeros (odeopts.NProc+1, length(y0));   %# matrix containing the fine solution for proc i at t_end
  y_tot        = zeros (odeopts.NProc+1, length(y0));   %# matrix containing the initial values
  y_tot(1,:)   = y0;

  %# initialize solution structure
  %# save solutions only if verbosity = 1
  sol.c = cell (1+(odeopts.NIter-1)*verbosity, odeopts.NProc);
  sol.f = cell (1+(odeopts.NIter-1)*verbosity, odeopts.NProc);

  %# start parareal
  k = 0;
  while (k < 1 || err(k) > 1)
   k  = k+1;
    kk = 1+(k-1)*verbosity;

    %# solve coarse grid
    if (verbosity)
      fprintf ('parareal: start %d/%d iteration\n', k, odeopts.NIter);
      fprintf ('parareal: start coarse grid\n');
      tim_c = tic;
    end
    for i = 1:odeopts.NProc
      sol.c{kk,i}   = feval (solver_c, ode_c, trange_proc(i,:), y_tot(i,:), opt_c, funarguments{:});
      %# work around for transpose bug of odepkg < 0.9
      %# see https://savannah.gnu.org/bugs/?49402
      if size(sol.c{kk,i}.y,1)~=length(y0)
        sol.c{kk,i}.y = sol.c{kk,i}.y.';
      end
      y_cold = y_cnew(i+1,:);
      %# if the coarse solution uses a lower fidelity
      %# upscaling by a prolongator is necessary
      if ~isfield(odeopts,'Prolongator') || isempty(odeopts.Prolongator)
        y_cnew(i+1,:) = sol.c{kk,i}.y(:,end);
      else
        y_cnew(i+1,:) = feval(odeopts.Prolongator,sol.c{kk,i}.y(:,end));
      end
      y_tot(i+1,:)  = y_cnew(i+1,:) + y_fend(i+1,:) - y_cold;
      if any(isnan(y_tot)) 
        fprintf ('parareal: isnan\n');
        keyboard
      end
    end
    
    if (verbosity)
      fprintf('parareal: coarse grid took %fs\n', toc (tim_c));
    end

    %# define fine propagator
    fun_f = @(i) feval (solver_f, ode_f, trange_proc(i,:), y_tot(i,:), opt_f, funarguments{:});
    
    %# solve (in parallel) using fine propagator
    if (verbosity)
      fprintf ('parareal: start fine grid\n');
      tim_f = tic;
    end

    if has_parallel % parallel fine-scale time stepping
      sol_f = pararrayfun (odeopts.NProc, fun_f, 1:odeopts.NProc,...
                           'UniformOutput', false, 'VerboseLevel', verbosity);
    else %# sequential fine-scale time stepping
      sol_f = {};
      for i = 1:odeopts.NProc
       sol_f{i} = fun_f(i);
      end
    end

    %# work around for transpose bug of odepkg < 0.9
    %# see https://savannah.gnu.org/bugs/?49402
    for i=1:odeopts.NProc
      if size(sol_f{i}.y,1)~=length(y0)
        sol_f{i}.y = sol_f{i}.y.';    % FIXME: is there a better way than a for?
      end
      if any(isnan(sol_f{i}.y)) 
        fprintf ('parareal: isnan\n');
        keyboard
      end
    end

    %# reconstruct global t and y vectors
    t = [trange(1)];
    y = [y0(:)'];
    
    temp_t = cellfun(@(c) c.x(2:end), sol_f, 'UniformOutput', false);
    t = [t; [temp_t{:}]'];
    temp_y = cellfun(@(c) c.y(:,2:end), sol_f, 'UniformOutput', false);
    y = [y; [temp_y{:}]'];
    [sol.f{kk,:}] = deal (sol_f);

    if (odeopts.haveoutputfunction)
      if (~isempty (odeopts.OutputSel))
        vals = y(:,odeopts.OutputSel);
      else
        vals = y;
      end
      stop_solve = feval (odeopts.OutputFcn, t, vals', [],...
                          odeopts.funarguments{:});
    end
    
    if (verbosity)
      fprintf ('parareal: fine grid took %fs\n', toc (tim_f));
    end

    %# compare solutions at start and end of each time interval
    %# use 2norm in space and inf-norm in time
    for i = 2:odeopts.NProc
      y_fend(i,:) = sol_f{i-1}.y(:,end);
      y_f0(i,:)   = sol_f{i}.y(:,1);
    end
    y_fend(end,:) = sol_f{end}.y(:,end);        
    err(k) = para_AbsRel_norm(y_fend(1:end-1,:),y_f0(1:end-1,:), odeopts.AbsTol, odeopts.RelTol);

    %# error message if the impossible happens and we never converge
    if (k > odeopts.NIter)
      warning (['parareal: not converge after odeopts.NIter=', ...
              num2str(odeopts.NIter), ' iterations']);
      break;
    end
  end % while parareal   

  %# Postprocessing
  if (odeopts.haveoutputfunction)  %# Cleanup plotter
    feval (odeopts.OutputFcn, [], [], 'done', odeopts.funarguments{:});
  end
  
  if nargout == 1
    sol.x = t.';            %# return time vector (as row)
    sol.y = y.';            %# return solution vector (as rows)
    sol.err = err(1:k);     %# return errors
    sol.f = sol.f(1:kk,:);  %# truncate sol.f 
    sol.c = sol.c(1:kk,:);  %# truncate sol.c
    varargout{1} = sol;
  elseif nargout == 2
    varargout{1} = t;   %# return time vector (as column)
    varargout{2} = y;   %# return solution vector (as columns)
  end

  if (verbosity)
    fprintf ('parareal: total execution time %fs\n', toc (tim_p));
  end
end


function infNorm = para_AbsRel_norm (x, y, AbsTol, RelTol)
%# inspired by Hairer, Norsett, Wanner 
%# Solving Ordinary Equations 1, page 167  
  sc = AbsTol + abs (x) .* RelTol;
  M = abs (x - y) ./ sc;
  twoNorm = sqrt(sum(M.^2,2))/sqrt(size(x,2));      %# 2norm of each row
  infNorm  = max(twoNorm);          %# infinity norm of the 2norm
end %function

%# COMMON TO ALL TESTS
%!function ydot = fpol (t, y)  # The Van der Pol ODE
%!  ydot = [y(2); (1 - y(1)^2) * y(2) - y(1)];
%!endfunction
%!function ref = fref ()       # The computed reference sol
%!  ref = [0.32331666704577, -1.83297456798624];
%!endfunction
%!
%!demo
%! vopt_c = odeset ('RelTol',1e-1,'AbsTol',1e-1,'NormControl','on');
%! vopt_f = odeset ('RelTol',1e-10,'AbsTol',1e-10,'NormControl','on');
%! vopt_p = odeparaset('Nproc',10,'NIter',10,'Verbose','on','OutputSel',1:2,'OutputFcn',@odeplot);
%! fvdb = @(vt,vy) [vy(2); (1 - vy(1)^2) * vy(2) - vy(1)];
%! sol = odeparareal (fvdb,'ode23',vopt_c,fvdb,'ode23',vopt_f,[0 40],[2 0],vopt_p);
%!
%!test  # for dimension of output and correctness
%! vopt_c = odeset ('RelTol',1e-1,'AbsTol',1e-1,'NormControl','on');
%! vopt_f = odeset ('RelTol',1e-10,'AbsTol',1e-10,'NormControl','on');
%! vopt_p = odeparaset('Nproc',10,'NIter',10,'Verbose','off');
%! [t y] = odeparareal (@fpol,'ode23',vopt_c,@fpol,'ode23',vopt_f,[0 2],[2 0],vopt_p);
%! assert (size(t,2)==1);
%! assert (size(t,1)==size(y,1) & size(y,2)==2);
%! assert ([t(end), y(end,:)], [2, fref], 1e-3);
