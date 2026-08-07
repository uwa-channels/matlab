function install()
%INSTALL Adds this toolbox to the MATLAB or Octave search path permanently.
%
% Run INSTALL once from the root of a cloned or downloaded copy of this
% repository.  It adds the src and examples folders to the search path and
% saves the path so the functions stay available in later sessions.
%
% MATLAB users can instead install the packaged toolbox from the Add-On
% Explorer or from a .mltbx file; see the README.  This script is the
% supported route for Octave, which cannot read .mltbx files.

here = fileparts(mfilename("fullpath"));

folders = {fullfile(here, "src"), fullfile(here, "examples")};
for k = 1:numel(folders)
    if ~isfolder(folders{k})
        error("install:missingFolder", "Expected folder is missing: %s", folders{k});
    end
    addpath(folders{k});
end

if savepath() ~= 0
    warning("install:savepathFailed", ...
        ["Could not save the path.  The toolbox works in this session, but " ...
         "you will need to run install again next time, or add the folders " ...
         "to your startup file."]);
    return
end

fprintf("Installed.  Try: help replay\n");
end
