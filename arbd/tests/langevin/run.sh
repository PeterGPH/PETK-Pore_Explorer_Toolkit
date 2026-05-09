#! /bin/bash

mkdir -p output
nvprof ../../src/arbd_lin_lin $@ Langevin.bd output/BrownDyn

