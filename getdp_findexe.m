function p = getdp_findexe
%GETDP_FINDEXE search some common places for getdp
%   Download GetDP from http://getdp.info

% The source code is distributed under the terms of the GNU General
% Public License (GPL) (version 2 or later).

  % windows binary name is different 
  if ispc
    filename = 'getdp.exe';
  else
    filename = 'getdp';
  end

  p = [];

  % get default path
  syspath  = getenv('PATH');
  dirs = regexp(syspath,pathsep,'split');
  
  % add some common places for binaries
  dirs = [dirs(:)',...
          {'/opt/bin'},...
          {'/usr/local/bin'},...
          {[getenv('HOME') '/bin']},...
          {pwd},...
          {[pwd '/getdp']},...
         ];
  
  for mydir = 1:numel(dirs)
      tp = fullfile(dirs{mydir},filename);
      if exist(tp,'file')
          p = tp;
          break
      end
  end

  if isempty(p)
      error('getdp not found, please install from http://getdp.info')
  end

end