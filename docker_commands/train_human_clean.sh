DEBUG=""
#DEBUG="gdb --args"

DATASET_FOLDER=/human_trainer_dataset
backup=/output/$(date +"%y-%m-%d_%T_train_human")

nprocesses=1
train=1
i=1
mkdir ${backup}

rm -rf tmp
mkdir tmp

export OMP_NUM_THREADS=4

q=" -b $DATASET_FOLDER -p 30 -a 5 "
q+=" -T template_human "

poses=("0.015" "0.015" "0.015" "0.015") 
device=0

ndeformjoints=2
label_term=0.5
lambda=3
m2d=0
robust_radius=1.5 #IGNORED!
shape_reg=1
model_reg=2

dropoff=5
cont_iter=125  #300
npersons=510  #15
nscans=1000
q+=" -j ${cont_iter} -r ${cont_iter} -d "
q+="experiments/data_config_male2d_5p250s.ini"
q+=":experiments/data_config_female2d_5p250s.ini" 
q+=":experiments/data_config_caesar_male2d_mixed_250p250s.ini"
q+=":experiments/data_config_caesar_female2d_mixed_250p250s.ini"
q+=" -D ${ndeformjoints} "
q+=" -E m_lhand,m_rhand "
muscle=1.5
pose=0.015
rigid_deform=1

npoints=17000

ENABLE_DATA_TERM=50
#q+=" -i ${ENABLE_DATA_TERM} "
q+=" -i 40 "
ncomp=10
continuous_step_size=0.01
njoints=17
q+=" -w ${continuous_step_size} -J ${njoints} -S "
MAX_VERTEX_RING=16
q+=" -V ${MAX_VERTEX_RING} "
q+=" -C 259 "
mkdir ${backup}/newout$i
echo "Init: " > ${backup}/newout$i/progress.txt
echo "njoints=${njoints}.0" > ${backup}/newout$i/config.t
echo "MAX_VERTEX_RING = ${MAX_VERTEX_RING}" >> ${backup}/newout$i/config.t
echo "nverts=1362.0" >> ${backup}/newout$i/config.t
echo "nfaces=1324.0" >> ${backup}/newout$i/config.t
echo "nlocal=10000.0" >> ${backup}/newout$i/config.t
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
echo "continuous_step_size=0.01" >> ${backup}/newout$i/config.t
echo "REDUCE_POSE_REG=85" >> ${backup}/newout$i/config.t
echo "USE_TEMPORAL=0" >> ${backup}/newout$i/config.t
echo "USE_SYMMETRY=1" >> ${backup}/newout$i/config.t
cat experiments/config4.t >> ${backup}/newout$i/config.t
CUDA_VISIBLE_DEVICES=$device ${DEBUG} ./trainer $p $q -c $ncomp -s ${backup}/newout$i/config.t -o ${backup}/newout$i -n ${npoints} #> /dev/null 2> ${backup}/newout$i/progress.txt &

