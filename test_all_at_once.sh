t1=4096
t2=65536
t3=17367040
max=10
param="params.txt"
program="BlackScholes"
tag=$1
current="~/GPU/Project/BlackSholesCUDA"
locate=$2
thr1=4
thr2=8

function test_cuda_all {
	echo generating performance test result files: test_mpi_4...$tag
	echo program: $program
	echo nthread: 128 
	echo mode:    0
	echo location:$locate
	make cuda_4096_run > $locate/test_cuda4096_$t1\_$tag
	make cuda_65536_run > $locate/test_cuda65536_$t2\_$tag
	make cuda_17367040_run > $locate/test_cuda17367040_$t3\_$tag
}

test_cuda_all
