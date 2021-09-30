MCCINC=-I ~/matlab/poly -I ~/matlab/fit -I ~/matlab/filter -I ~/matlab/functions -I ~/matlab/misc -I ~/matlab/io -I ~/matlab/util -I ~/matlab/Dedalus -I ~/matlab/graphic_util -I ~/matlab/statistics -I ~/matlab/fortran -I ~/matlab/struct -I ~/matlab/cell  -I ~/matlab/video

default:ded_poseidon

ded_poseidon:ded_poseidon.m 
	mcc $(MCCINC) -m -R -nodesktop -o $@ ded_poseidon


