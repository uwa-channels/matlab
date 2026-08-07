function outputFile = buildToolbox(version)
%BUILDTOOLBOX Package this repository into an installable .mltbx file.
%
% BUILDTOOLBOX() reads the version string from the VERSION file at the
% repository root, stages the toolbox files under <repo>/build, and writes
% uwa-channels-<version>.mltbx next to that staging folder.
%
% BUILDTOOLBOX(VERSION) overrides the version string, for example "1.2.0".
% The release workflow passes the Git tag this way.
%
% OUTPUTFILE = BUILDTOOLBOX(...) returns the path of the generated file.
%
% Running this function requires R2023a or later.  The toolbox it produces
% installs on R2021a and later.
%
% Example:
%    buildToolbox
%    matlab.addons.install('build/uwa-channels-1.0.0.mltbx')

% Keep this identifier fixed forever.  MATLAB uses it to recognize that a
% new .mltbx is an upgrade of an installed toolbox rather than a separate
% add-on, so changing it would leave users with two copies on the path.
% It need not be a GUID; a readable name shows up in the Identifier column
% of matlab.addons.installedAddons and makes the uninstall command legible.
identifier = "uwa-channels";

repoRoot = fileparts(fileparts(mfilename("fullpath")));

if nargin < 1 || isempty(version)
    version = strtrim(fileread(fullfile(repoRoot, "VERSION")));
end
version = string(version);
if isempty(regexp(version, "^\d+(\.\d+){1,3}$", "once"))
    error("buildToolbox:badVersion", ...
        "Version must be 2 to 4 dot-separated numbers, got '%s'.", version);
end

buildDir = fullfile(repoRoot, "build");
stageDir = fullfile(buildDir, "uwa-channels");
if isfolder(buildDir)
    rmdir(buildDir, "s");
end
mkdir(stageDir);
mkdir(fullfile(stageDir, "examples"));

copyfile(fullfile(repoRoot, "src", "*.m"), stageDir);
copyfile(fullfile(repoRoot, "examples", "*.m"), fullfile(stageDir, "examples"));
copyfile(fullfile(repoRoot, "LICENSE"), stageDir);
copyfile(fullfile(repoRoot, "README.md"), stageDir);

opts = matlab.addons.toolbox.ToolboxOptions(stageDir, identifier);
opts.ToolboxName = "Underwater Acoustic Channel Toolbox";
opts.ToolboxVersion = version;
opts.Summary = "Replay signals through measured underwater acoustic " + ...
    "channels and generate realistic ocean noise.";
opts.Description = join([ ...
    "Replays passband signals through measured underwater acoustic channels,"
    "generates realistic ocean noise (pink Gaussian, spatially correlated,"
    "or impulsive alpha-stable), and reconstructs full time-varying impulse"
    "responses from the compressed channel representation."
    ""
    "Channel files: https://doi.org/10.5281/zenodo.21287414"
    "Documentation: https://uwa-channels.github.io/"], " ");
opts.AuthorName = "Underwater Acoustic Channels Group";
opts.AuthorCompany = "Underwater Acoustic Channels Group";
opts.MinimumMatlabRelease = "R2021a";
opts.ToolboxMatlabPath = [stageDir; fullfile(stageDir, "examples")];
opts.OutputFile = fullfile(buildDir, "uwa-channels-" + version + ".mltbx");

matlab.addons.toolbox.packageToolbox(opts);

% The release also ships an unversioned copy.  GitHub serves the assets of
% the newest release under a fixed /releases/latest/download/<name> URL, so
% a stable file name is what lets the documented one-line install command
% never mention a version.  The version still travels inside the package.
stableFile = fullfile(buildDir, "uwa-channels.mltbx");
copyfile(opts.OutputFile, stableFile);

outputFile = opts.OutputFile;
fprintf("Packaged %s\n", outputFile);
fprintf("Packaged %s\n", stableFile);
end
