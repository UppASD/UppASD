set -e
SP="$1"
PRE=./build_rcg08_base/bin/gpu_adaptive_runtime_benchmark
POST=./build_rcg08_cuda_fp64/bin/gpu_adaptive_runtime_benchmark
ARGS="--blocks 2048 --atoms-per-block 4 --repetitions 5 --iterations 10"
: > $SP/ab_pre.txt; : > $SP/ab_post.txt
for round in 1 2 3; do
  echo "== round $round ==" | tee -a $SP/ab_pre.txt >> $SP/ab_post.txt
  nvidia-smi --query-gpu=index,utilization.gpu,clocks.sm,temperature.gpu --format=csv,noheader >> $SP/ab_env.txt
  if [ $((round % 2)) -eq 1 ]; then
    $PRE $ARGS >> $SP/ab_pre.txt;  $POST $ARGS >> $SP/ab_post.txt
  else
    $POST $ARGS >> $SP/ab_post.txt; $PRE $ARGS >> $SP/ab_pre.txt
  fi
done
echo AB_DONE
