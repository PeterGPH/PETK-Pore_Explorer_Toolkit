#! /bin/bash

mkdir output
nvprof ../../src/arbd_lin_lin $@ BrownDyn.bd output/BrownDyn

