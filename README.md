# SETUP


## Submodule Revisions

The src submodule is linked to a specific commit of the https://lex-gitlab.aer.com/AQAC/common_ASP repository, and modifications may have been made at numerous points. This is fine and indeed what Git intends (see flow charts in http://paynedigital.com/articles/2011/10/introduction-to-git-submodules), but one may want to keep track of what was modified. We can do this by examining the working copy, the latest commit into the local repository (HEAD), or the copy in the central repository (i.e., the Gitlab host).

In the next few examples, let's assume I made a small change to the EulerianBoxDriver.f90 file.

### Working Copy vs. Local Repository

```
% git diff EulerianBoxDriver.f90
```

The output will be the standard output from a call to the Linux "diff" command.

I prefer the "tkdiff" approach, which is just a graphical "diff" that nicely visualizes (e.g., color codes) the differences between two files. Assuming this (or something similar) is installed on the user's system, they can set the difftool setting:

```
% git config --global diff.tool tkdiff
% git config --global merge.tool tkdiff
% git config --global --add difftool.prompt false
% git config --list | grep diff
diff.tool=tkdiff
difftool.prompt=false
```

Then we can use the command:

```
% git difftool EulerianBoxDriver.f90
```

If called like this verbatim, the HEAD version is in the left tkdiff panel, and the working copy is in the right panel.

### Local Repository vs. Central Repository

Should the simple change to the working copy of EulerianBoxDriver.f90 be committed into the user's local repository, the differences can be compared with what is in the central repository in Gitlab with:

```
% git difftool HEAD...origin/master
```

In tkdiff, the most recent commit in the central repo is on the right, and my local repo copy (HEAD) is on the left. these comparisons assume we just wanna compare the HEAD copies with the master branch -- origin/master could presumably just be changed to origin/whateverBranch if branches for the repo exist (which they don't, so i have not tested).

## Updating the Submodules

To update the src submodule to the current state of the master branch of the common_ASP repository:

```
% pwd
/rd47/scratch/ASP/GitASP
% git submodule foreach git pull origin master
```

One can also likewise push any changes in src/ to https://lex-gitlab.aer.com/AQAC/common_ASP.

## Building and Running the Model

See `How_to_run_ASP.pdf` at the top of this repo. Also, from [Rick's Docker Cloud ASP Page](<https://cloud.docker.com/repository/docker/pernak18/asp>):

To run a container, try:

```
% docker run -v ${YOUR_OUTPUT_DIR}/userOutputFiles:/asp/OutputFiles -v ${YOUR_INPUT_DIR}/userInputFiles:/asp/userInputFiles pernak18/asp
```

It is assumed that `${YOUR_INPUT_DIR}` is a directory with the *.in files that ASP expects. Typically, one can just do a copy of the defaults:

```
% cp -r InorganicInputDecks/ ${YOUR_INPUT_DIR}
% rm -rf ${YOUR_INPUT_DIR}/defaults ${YOUR_INPUT_DIR}/VBS_draft
```

This method is recommended because then the user inputs are not symbolic links (so there is no dependence on the host paths), but it's important to execute the `rm` command because the entrypoint script does an `ls` on `${YOUR_INPUT_DIR}`, removes the files in the container's `InorganicInputDecks/` directory, then replaces the container's original input files with the user-provided input files. If the `rm` is not performed, an attempt at removing and copying directories is done, and the entrypoint script will crash.

This command should work as well but has not been tested:

```
% docker run -v {ANY_OUTDIR}:/asp/OutputFiles -v {ANY_INDIR}:/asp/userInputFiles asp -i {ANY_INDIR}
```
