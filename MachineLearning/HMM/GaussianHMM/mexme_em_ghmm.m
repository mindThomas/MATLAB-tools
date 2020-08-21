echo on

mex  em_ghmm.c

mex -Dbetanormalize  forward_backward.c

mex  likelihood_mvgm.c

mex  ndellipse.c

mex  -DranSHR3 sample_ghmm.c


echo off