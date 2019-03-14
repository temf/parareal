function numdofs = getdp_preresolution(file)
%GETDP_PRESOLUTION read GetDP's pre file.
%   Download GetDP from http://getdp.info

% The source code is distributed under the terms of the GNU General
% Public License (GPL) (version 2 or later).

  fid = fopen(file);
  numdofs = [];

  while feof(fid) == 0
    line = fgetl(fid);
    if findstr(line,'$DofData')
      % TODO: read the other values
      line1 = fgetl(fid);
      line2 = fgetl(fid);
      line3 = fgetl(fid);
      line4 = fgetl(fid);
      % parse only numdof for now
      line5 = fgetl(fid);
      tmp  = sscanf(line5,'%d %d');
      numdofs(end+1,1) = tmp(2);
    end %if
  end %while

  fclose(fid);

end %function