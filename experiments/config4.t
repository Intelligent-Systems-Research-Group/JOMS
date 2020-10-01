--ndeformjoints = 2
nbasis = 1
lambda = 100
nlabels = 25.0 --63.0
njointlabels = 0 --13
njoints = 33
newWeight = 10.0
basis_max = 16
local_size = 28
pose_size = 3*njoints
THRESH = 30
THRESH_MUL = 1
ZERO_JOINT = 2
scaling = 1--1000
--lambda_sparse = 5 --5
LABEL_ITER_MUL = 85.0/30.0
SHAPE_PRIOR_PCS = 2
--INIT_MEAN_REG = 30 --30
--MEAN_DROPOFF=5
WEAKEN_LATENT_CONSTRAINT=3000
USE_TEMPORAL=0

local file = assert(loadfile("experiments/shapefit.t"))
env = getfenv(1)
setfenv(file,env)
file()
