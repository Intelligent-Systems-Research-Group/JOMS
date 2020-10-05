DEBUG=""
#DEBUG="gdb --args"

DATASET_FOLDER=/human_trainer_dataset
backup=/output/$(date +"%y-%m-%d_%T_train_hello_world")

train=1

mkdir ${backup}

rm -rf tmp
mkdir tmp

export OMP_NUM_THREADS=4
q=""
q+=" -b $DATASET_FOLDER -p 10 -a 999 -e 90"

device=0
ndeformjoints=0
label_term=.5
lambda=3
robust_radius=1.5
m2d=0
shape_reg=.5
model_reg=2
dropoff=5
cont_iter=2  #300
npersons=3  #15
nscans=3
q+=" -j ${cont_iter} -r ${cont_iter} -d "
q+="experiments/data_config_hello_3p3s.ini"
q+=" -D ${ndeformjoints} "
q+=" -T template_hello_world "
muscle=1.5
pose=0.015
rigid_deform=1
npoints=1000
ENABLE_DATA_TERM=0
q+=" -i ${ENABLE_DATA_TERM} "
#ncomp is the number of shape blend-shapes
ncomp=1
continuous_step_size=1
#no joints are used in the model, this just a limitation of the joms trainer
njoints=4
q+=" -w ${continuous_step_size} -J ${njoints}"
MAX_VERTEX_RING=64
q+=" -V ${MAX_VERTEX_RING} "
q+=" -C 0"
mkdir ${backup}/newout$i
echo "Init: " > ${backup}/newout$i/progress.txt
echo "njoints=${njoints}.0" > ${backup}/newout$i/config.t
echo "MAX_VERTEX_RING = ${MAX_VERTEX_RING}" >> ${backup}/newout$i/config.t
echo "nverts=20.0" >> ${backup}/newout$i/config.t
echo "nfaces=20.0" >> ${backup}/newout$i/config.t
echo "nlocal=20.0" >> ${backup}/newout$i/config.t
echo "nscans=${nscans}" >> ${backup}/newout$i/config.t 
echo "npersons=${npersons}" >> ${backup}/newout$i/config.t
echo "ncomp=$ncomp" >> ${backup}/newout$i/config.t
echo "npoints=$npoints" >> ${backup}/newout$i/config.t
echo "TRAIN=true" >> ${backup}/newout$i/config.t
echo "RIGID_DEFORM_FACTOR=${rigid_deform}" >> ${backup}/newout$i/config.t
echo "ENABLE_DATA_TERM=${ENABLE_DATA_TERM}" >>${backup}/newout$i/config.t 
echo "MODEL_TO_DATA_TERM=${m2d}" >>${backup}/newout$i/config.t 
echo "LABEL_TERM=${label_term}" >>${backup}/newout$i/config.t 
echo "ROBUST_RADIUS=1" >> ${backup}/newout$i/config.t 
echo "MEAN_DROPOFF=${robust_radius}" >> ${backup}/newout$i/config.t
echo "INIT_MEAN_REG=${shape_reg}" >> ${backup}/newout$i/config.t
echo "ndeformjoints=${ndeformjoints}" >> ${backup}/newout$i/config.t
echo "MUSCLE_DEFORM_FACTOR=${muscle}" >> ${backup}/newout$i/config.t
echo "REST_MODEL_REG=${model_reg}" >> ${backup}/newout$i/config.t
echo "POSE_REG=${pose}" >> ${backup}/newout$i/config.t
echo "lambda_sparse=${lambda}" >> ${backup}/newout$i/config.t
echo "COEFF_REG=.0075" >> ${backup}/newout$i/config.t
echo "ENABLE_GROUND=65" >> ${backup}/newout$i/config.t #161
echo "continuous_step_size=${continuous_step_size}" >> ${backup}/newout$i/config.t
echo "REDUCE_POSE_REG=85" >> ${backup}/newout$i/config.t
echo "USE_TEMPORAL=0" >> ${backup}/newout$i/config.t
echo "USE_SYMMETRY=0" >> ${backup}/newout$i/config.t
cat experiments/config4.t >> ${backup}/newout$i/config.t
CUDA_VISIBLE_DEVICES=$device ${DEBUG} ./trainer $p $q -c $ncomp -s ${backup}/newout$i/config.t -o ${backup}/newout$i -n ${npoints} #> /dev/null 2> ${backup}/newout$i/progress.txt &
