#!/bin/bash
# ----------------------------------------------------------------------------
# USAGE: bash update_git.sh 'some_bullshit_comment_about_whats_been_changed'
# ----------------------------------------------------------------------------

# Update/ pull the .git repository from the Master version,
# i.e. check there has been no local conflicting changes before pulling
git fetch && git log ..origin/master
git fetch && git log origin/master..
git pull

# Add local changes to the .git repository
git add *
git commit -m $1

# Commit/upload local changes to the Master
git push -u origin master
