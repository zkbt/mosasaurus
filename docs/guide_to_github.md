# `github` cheatsheet

For those us new to collaborative coding with github, here are a bunch of commands.

Create a new branch:

    git branch issue#10-bjd

Switch over to that branch, for testing and developing:

    git checkout issue#10-bjd

Add files to the staging area for this branch:

    git add .

Commit changes to this branch (for staged files):

    git commit -m "replace Irwin's lfa with astropy"

Push this branch's changes to the repository:

    git push --set-upstream origin issue#10-bjd

Fetch all branches from remote, and remove dead ones (?):

    git fetch --all --prune
