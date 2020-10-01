echo $1
echo $2
echo $3
p="$(pwd)"
docker run --rm -it \
--gpus all \
--mount type=bind,source=$p/weights,target=/weights,readonly \
--mount type=bind,source=$p/dataset,target=/human_trainer_dataset,readonly \
--mount type=bind,source=$p/output,target=/output \
joms/trainer:latest bash
#-e INI="$1" \
#-e npersons="$2" \
#-e nscans="$3" \
#hypermod/trainer:latest bash
