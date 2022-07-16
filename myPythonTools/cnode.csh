#!/bin/csh -f
#PBS -q long
# Number of nodes (CHANGE THIS if needed)
#PBS -l walltime=1:00:00,nodes=1:ppn=16
# output file base name
#PBS -N interactive
# Put standard error and standard out in same file
#PBS -j oe
# Export all Environment variables
#PBS -V

