#!/bin/sh

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/

module load matlab

matlab -nodesktop $* -r "addpath('compute_persistence', 'matpcl');"
