#!/bin/bash

ls -ltd */ | sed "1,$1d" | rev | cut -d ' ' -f1 | rev &> dir_list
if [ -s dir_list ]; then
  echo "Removing:"
  cat dir_list
  srun -N 1 -n 20 -p pdebug -t 30 drm `< dir_list`
else
  echo "No directories to remove at this time."
fi
rm dir_list
