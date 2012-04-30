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
	echo generating performance test result files: test_cuda_...$tag
	echo program: $program
	echo nthread: 128 
	echo mode:    0
	echo location:$locate
	make cuda_4096_run > $locate/test_cuda_$t1\_$tag
	make cuda_65536_run > $locate/test_cuda_$t2\_$tag
	make cuda_131072_run > $locate/test_cuda_131072\_$tag
	make cuda_8388608_run > $locate/test_cuda_8388608\_$tag
	make cuda_17367040_run > $locate/test_cuda_17367040\_$tag
	make cuda_exceed_run > $locate/test_cuda_67108864\_$tag
}

test_cuda_all
