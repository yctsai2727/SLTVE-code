#!/bin/bash
mkdir -p ParallelSol
for id in {0..7}
do
    octave -qf driveCN.m $id
done