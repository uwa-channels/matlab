# Packaging and releasing

MATLAB has no command-line package index like PyPI, so there is no literal `pip install uwa-channels`.  The closest equivalent is a **toolbox package** (`.mltbx`) published as a GitHub release asset under a fixed name, which users install with one command:

```matlab
matlab.addons.install(websave([tempname '.mltbx'], 'https://github.com/uwa-channels/matlab/releases/latest/download/uwa-channels.mltbx'))
```

Three things make that command work, and each one is easy to break:

* **`websave` returns the path it wrote**, so it nests directly inside `matlab.addons.install` with no temporary variable.
* **`matlab.addons.install` requires an absolute path.**  Given a relative one it reports "The file is not a valid toolbox file", which is misleading.  `tempname` returns an absolute path, so the command is safe.
* **The asset name carries no version.**  GitHub resolves `/releases/latest/download/<name>` to the newest non-prerelease release, so the command in the README never has to change.  This is why the build emits `uwa-channels.mltbx` alongside `uwa-channels-<version>.mltbx`.

## Cut a release

1. Bump the version in `VERSION` and commit it.
2. Tag and push:

   ```
   git tag v1.0.1
   git push origin v1.0.1
   ```

That is the whole release.  The `Release` workflow then:

1. checks that the tag matches `VERSION`, failing rather than shipping a mislabeled package;
2. builds `build/uwa-channels-<version>.mltbx` and the unversioned copy;
3. installs the built file and asserts `replay`, `noisegen`, and `unpack` land on the path;
4. creates the GitHub release if it does not exist yet, or uploads to it if it does, attaching both `.mltbx` files;
5. runs the README's install command against the published release and asserts the installed version equals the tag.

Step 5 means a broken install URL fails the release instead of reaching users.

## Build locally

Requires R2023a or later to run; the package it produces installs on R2021a and later.

```matlab
addpath('packaging')
buildToolbox            % version from the VERSION file
buildToolbox('1.2.0')   % or an explicit version
```

`buildToolbox` wipes and recreates `build/`, so it does not accumulate stale packages.  Test the result with

```matlab
t = matlab.addons.install(fullfile(pwd, 'build', 'uwa-channels.mltbx'));
which replay
matlab.addons.uninstall(t.Identifier)
```

## Do not change the identifier

`buildToolbox.m` hard-codes the identifier `uwa-channels`.  MATLAB uses it to tell that a new `.mltbx` upgrades an installed toolbox instead of being a separate add-on.  Changing it would leave existing users with two copies of the functions on their path.

The packaging dialog defaults to a GUID, but the identifier is a free-form string; MathWorks' own documentation example uses `com-mathworks-guilayout`.  A readable one is worth having because it is what appears in the `Identifier` column of `matlab.addons.installedAddons` and what users type to uninstall:

```matlab
matlab.addons.uninstall('uwa-channels')
```

## Dependencies

`.mltbx` cannot express "requires Signal Processing Toolbox".  `ToolboxOptions` has no product-requirement property at all, and its `RequiredAddons` property is for other add-ons rather than licensed MathWorks products: each entry needs `Name`, `Identifier`, `EarliestVersion`, `LatestVersion`, and `DownloadURL` fields, and the installer tries to *download and install* them, which is the wrong behavior for a product the user must license.

So the toolbox declares only what the format supports, `MinimumMatlabRelease = "R2021a"`, and the product requirement stays in the README.  The authoritative list comes from MATLAB's own dependency analysis:

```matlab
[~, prods] = matlab.codetools.requiredFilesAndProducts( ...
    {'src/replay.m','src/noisegen.m','src/unpack.m'});
{prods.Name}    % MATLAB, Signal Processing Toolbox
```

Re-run that after adding any function call, and update the README if the list grows.  Note that the alpha-stable noise generator in `noisegen.m` is implemented by hand precisely so it does not pull in Statistics and Machine Learning Toolbox.

## Optional: File Exchange listing

Publishing to [MATLAB File Exchange](https://www.mathworks.com/matlabcentral/fileexchange) additionally puts the toolbox in the in-product Add-On Explorer (**Home > Add-Ons > Get Add-Ons**), which some users prefer over pasting a command.  It is a one-time manual step on the MathWorks side, not something this workflow can do: sign in with the group account, choose **Publish > From GitHub**, give `https://github.com/uwa-channels/matlab`, and title it `Underwater Acoustic Channel Toolbox` to match `ToolboxName` in `buildToolbox.m`.  After that, File Exchange re-syncs on each new GitHub release.
