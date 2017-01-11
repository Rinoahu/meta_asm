#!/bin/bash

rm canopies.npy

git config --global user.email xiaohu@iastate.edu
git config --global user.name Rinoahu


git remote rm origin

git add -A .
git commit -m 'first_commit'
git remote add origin https://github.com/Rinoahu/meta_asm
git pull origin master
git push origin master

git checkout master
