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
