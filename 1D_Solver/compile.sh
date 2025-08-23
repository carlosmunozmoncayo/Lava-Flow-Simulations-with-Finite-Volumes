#!/bin/bash

f2py -c newton_solver.f90 -m newton_solver -llapack -lblas
f2py -c rp1_geoclaw_swe_rs.f90 -m geoclaw_swe_rs