function numdof = getdp_numdof (fun, trange, varargin)
%GETDP_NUMDOF extract the information of degrees of freedom from .
%   Download GetDP from http://getdp.info

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
  resfile = fullfile(tmpdir,[fname '.res']);

  % Do preprocessing step
  exe_string = [getdpopts.Executable ' "' fun  '"'...
                ' -pre "' getdpopts.PreProcessing '"' ...
                ' -msh "' mshfile '"'...
                ' -name "' tmpname '"'...
                ' -res "' resfile '"'...
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
  numdof = getdp_preresolution(prefile);

  % delete the temporaray directory
  [status] = rmdir(tmpdir,'s');

end

