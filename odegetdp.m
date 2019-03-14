function varargout = odegetdp (fun, trange, init, varargin)
%ODEGETDP Calls getdp with a syntax similar to odeXY solvers.
%   Download getdp from http://getdp.info

% The source code is distributed under the terms of the GNU General
% Public License (GPL) (version 2 or later).

  if (nargin >= 4)
    if (~ isstruct (varargin{1}))
      getdpopts = struct ();
      funarguments = varargin;
    elseif (numel (varargin) > 1)
      %# varargin{1} is an options structure opt
      getdpopts = varargin{1};
      funarguments = {varargin{2:numel (varargin)}};
    else  % if (isstruct (varargin{1}))
      getdpopts = varargin{1};
      funarguments = {};
    end %if
  else  % nargin == 3
    getdpopts = struct ();
    funarguments = {};
  end %if

  if (~isnumeric (trange) || ~isvector (trange))
    error ('Octave:invalid-input-arg', ...
           'odegetp: TRANGE must be a numeric vector');
  end %if

  if (numel (trange) < 2)
    error ('Octave:invalid-input-arg', ...
           'odegetp: TRANGE must contain at least 2 elements');
  elseif (trange(1) >= trange(2))
    error ('Octave:invalid-input-arg', ...
           'odegetp: invalid time span, TRANGE(1) >= TRANGE(2)');
  end %if
  trange = trange(:);

  if (~ isnumeric (init) || ~ isvector (init) || any(isnan(init)) )
    error ('Octave:invalid-input-arg', ...
           'odegetp: INIT must be a numeric vector');
  end %if
  init = init(:);

  % Check if model exists
  if (~ischar(fun))
      fun = func2str(fun);
  end
  [fdir,fname,fext] = fileparts(fun);
  if ~exist(fun,'file') || ~strcmp(fext,'.pro')
    error ('Octave:invalid-input-arg', ...
           'odegetdp: FUN must be a valid getdp FUN.pro file');
  else % get absolute path
    fdir = cd(cd(fdir));
  end %if

  % Check for getdp executable
  if ~isfield(getdpopts,'Executable')
    getdpopts.Executable = 'getdp';
  end %if
  [status, output] = system([getdpopts.Executable ' --version' ' 2>&1']);
  if status
    error ('Octave:invalid-input-arg', ...
           'odegetdp: getdp not found, you may use "Executable" in the option structure.');
  else
    getdpver = regexp(output,'\d\.\d\d\.\d','match','once');
  end %if

  if ~isfield(getdpopts,'Verbose')
    getdpopts.Verbose=0;
  end %if

  % Choose TimeStep
  if ~isfield(getdpopts,'TimeStep')
    getdpopts.TimeStep = (trange(end)-trange(1))/100;
    fprintf(['odegetdp: choosing default time step dt=' num2str(getdpopts.TimeStep) 's\n'])
  elseif getdpopts.TimeStep<=0
    error ('Octave:invalid-input-arg', ...
           'odegetdp: Time step must be positive.');
  end %if

  % Process optional funarguments
  funargstr = '';
  for i=1:2:length(funarguments)
    if isnumeric(funarguments{i+1})
      funargstr = [funargstr ' -setnumber ' funarguments{i} ' ' num2str(funarguments{i+1})]; 
    else
      funargstr = [funargstr ' -setstring ' funarguments{i} ' ' funarguments{i+1}]; 
    end %if
  end %if

  % Choose PreProcessing
  if ~isfield(getdpopts,'PreProcessing')
    getdpopts.PreProcessing = '#1';
    %fprintf('odegetdp: choosing default PreProcessing "#1"\n')
  end %if

  % Choose first Resolution
  if ~isfield(getdpopts,'Resolution')
    getdpopts.Resolution = '#1';
    %fprintf('odegetdp: choosing default Resolution "#1"\n')
  end %if

  % Check if mesh exists
  mshfile =  fullfile(fdir,[fname '.msh']);
  if ~exist(mshfile,'file')
    error(['odegetdp: run gmsh on beforehand to create a msh-file: ' mshfile]);
  end %if

  % Create a file with the given initial value
  tmpdir  = tempname();
  mkdir(tmpdir);
  tmpname = fullfile(tmpdir,fname);
  resdir  = fullfile(tmpdir,'res');
  prefile = fullfile(tmpdir,[fname '.pre']);
  resfile = fullfile(tmpdir,[fname '.res'])

  % Do preprocessing step
  exe_string = [getdpopts.Executable ' "' fun '"' ...
                ' -pre "' getdpopts.PreProcessing '"' ...
                ' -msh "' mshfile '"' ...
                ' -name "' tmpname '"' ...
                ' -res "' resfile '"' ...
                ' -setnumber timemax ' num2str(trange(2)) ...
                ' -setnumber dtime ' num2str(getdpopts.TimeStep) ...
                ' -setstring ResDir "' resdir '"' ...
                funargstr ' 2>&1'];

  % Error and debug output
  if getdpopts.Verbose==1
    [status] = system(exe_string);
  else
    [status,output] = system(exe_string);
  end %if
  if status
    error('odegetdp: preprocessing failed');
  end %if

  % check if init is correct
  numdof  = length(init);
  numpres = getdp_preresolution(prefile);
  if (numdof~=sum(numpres))
    warning (['odegetdp: INIT has wrong size: ' int2str(numdof) ' instead of ' int2str(numpres) ': ' prefile]);
  end

  % Create initial data
  set_resolution(resfile,trange(1),init,numdof);

  % Choose PostOperation
  if isfield(getdpopts,'PostOperation')
    funargstr = [funargstr ' -pos ' getdpopts.PostOperation];
  end %if

  % Launch GetDP to get solution and postproc
  exe_string = [getdpopts.Executable ' "' fun '"'...
                ' -restart' ... %"' getdpopts.Resolution '" ' ...
                ' -msh "' mshfile '"' ...
                ' -name "' tmpname '"' ...
                ' -res "' resfile '"' ...
                ' -setnumber timemax ' num2str(trange(2)) ...
                ' -setnumber dtime ' num2str(getdpopts.TimeStep) ...
                ' -setstring ResDir "' resdir '"' ...
                funargstr ' 2>&1'];
  
  if getdpopts.Verbose==1
    [status] = system(exe_string);
  else
    [status,output] = system(exe_string);
  end %if
  if status
    error('odegetdp: solving failed');
  end %if

  % Read results
  [t,y] = getdp_read_resolution(resfile,numdof);

  % Return solution
  if nargout==2
    varargout{1} = t';
    varargout{2} = y';
  else
    varargout{1} = struct('x',t,'y',y,'solver','odegetdp',...
                          'version',getdpver,'TempName',tmpname);
  end %if
  %[status] = rmdir(tmpdir,'s');

end %function

function set_resolution(file,t,x,numdofs)
  fid = fopen(file,'w');
  % get positions of dofdata in x vector
  dofpos = cumsum([0;numdofs]);
  % start writing res file
  fprintf(fid, '$ResFormat /* GetDP 2.10.0, ascii */\n');
  fprintf(fid, '1.1 0\n');
  fprintf(fid, '$EndResFormat\n');
  for j=1:length(t)
    for k=1:length(numdofs)
      fprintf(fid, '$Solution  /* DofData #%d */\n',num2str(k-1));
      fprintf(fid, '%d %.16f 0 %d\n',k-1,t(j),j-1);
      y(1,:)=x([1+dofpos(k):dofpos(k+1)]);
      fprintf(fid, '%.16f %.16f\n',[real(y);imag(y)]);
      fprintf(fid, '$EndSolution\n');
    end %for k
  end %for j
  fclose(fid);
  if any(isnan(x)) || any(isnan(t))
    warning('odegetdp: something went wrong')
  end
end %function

function [t,x] = getdp_read_resolution(file,numdofs)
  fid = fopen(file);
  % init solution vector, may contain several dofdata sets
  x = zeros(sum(numdofs),0);
  % init vector of time steps
  t = [];
  % init vector of time step numbers
  j = 0;
  oldstep = 0;
  % get positions of dofdata in x vector
  dofpos = cumsum([0;numdofs]);
  while feof(fid) == 0
    line = fgetl(fid);
    val  = str2double(line);
    if findstr(line,'$Solution')
      line = fgetl(fid);
      % dofdata-number time-value time-imag-value time-step-number
      tmp  = sscanf(line,'%d %f %f %d');
      % check if the time-step-number has advanced
      if oldstep<1+tmp(4)
        j = j + 1;
        oldstep = 1 + tmp(4);
      elseif oldstep>1+tmp(4)
        error('odegetdp: error reading file %s.\ntime step %d is stored after %d',...
              file, tmp(4),oldstep-1);
      end
      k = 1 + tmp(1);
      t(1,j) = tmp(2);
      % read complex dofdata set into solution vector
      xtmp  = fscanf(fid,'%f',[2 numdofs(k)])';
      x([1+dofpos(k):dofpos(k+1)],j) = xtmp(:,1)+imag(xtmp(:,2));
    elseif findstr(line,'$ResFormat')
      line = fgetl(fid);
      if ~strcmp(line(1:3),'1.1')
        error('odegetdp: unknown file format version');
      end %if wrong version
    end %if solution
  end %while
  fclose(fid);
  if any(any(isnan(x))) || any(isnan(t))
    error('getdp_read_resolution: file contains NaN')
  end
end %function