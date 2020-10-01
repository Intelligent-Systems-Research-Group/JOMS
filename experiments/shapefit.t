--MAX_VERTEX_RING = 16
ndeformshapes = 9

cfg = 0
function inc()
   cfg = cfg+1
   return cfg-1
end

function pretty(x)
   for i=3,table.getn(x),3 do
      print(x[i])
   end
end
--@TODO To config

function softmax(w)
    local out = {}
    for i=1,table.getn(w) do
        --out[i] = Select(mask[i],1.0 / (1 + ad.exp(-lambda*w[i])) ,0)
        --out[i] = Select(mask[i], .5 + .5*w[i] / sqrt(1 + w[i]*w[i]) ,0)
        out[i] = .5 + .5*w[i] / sqrt(1 + w[i]*w[i])
    end
    return out
end

function SEnergy(e)
	
    Energy(scaling*e)
    
end

N, MW, M3, M2, M1 = 
	Dim("W",0), Dim("MW",1), Dim("MP",2), Dim("M3",3), Dim("M2",4)
--JOINT_DTYPE = opt_float51

K3 = TRAIN and N or M3
K2 = TRAIN and N or M2
K1 = TRAIN and N or M1

X3 = Unknown("X3",opt_float3,{N},inc())
X2 = Unknown("X2",opt_float2,{N},inc())
X1 = Unknown("X1",opt_float,{N},inc())
DW = Array("DW",opt_float28,{MW},inc())
D3 = Array("D3",opt_float3,{M3},inc())
D2 = Array("D2",opt_float2,{M2},inc())
D1 = Array("D1",opt_float,{M1},inc())

Y3 = TRAIN and X3 or D3
Y2 = TRAIN and X2 or D2
Y1 = TRAIN and X1 or D1

local w_fitSqrt = Param("w_fitSqrt", opt_float, inc())+1
local w_surface = Param("w_surface", opt_float, inc())

--lambda_sparse = (nscans/50.0)*5

lambda_sparse = Select(lesseq(w_fitSqrt,130),lambda_sparse,4*lambda_sparse)
lambda_sparse = Select(lesseq(w_fitSqrt,140),lambda_sparse,4*lambda_sparse)
lambda_sparse = Select(lesseq(w_fitSqrt,150),lambda_sparse,4*lambda_sparse)
lambda_sparse = Select(lesseq(w_fitSqrt,160),lambda_sparse,2*lambda_sparse)



function decay(input,base,exponent)
    for j=1,exponent do
        input = Select(lesseq(w_fitSqrt,j*THRESH),input,base*input)
    end
    return input
end

function disable()
    --return Select(less(sin(2*3.141*w_fitSqrt/10),0),0,1)
    return w_surface
end

function basis_bspline_deriv(u)
    local t = u
    local t2 = t*t
    local base =  {}
    base[0] = -0.5*t2 +     t    -.5
    base[1] =  1.5*t2 - 2.0*t
    base[2] = -1.5*t2 +     t    +.5
    base[3] =  0.5*t2
    return base
end

function patch_bspline_deriv(u,v,sweights,tweights)
    local derivU = {}
    local derivV = {}
    local dsWeights = basis_bspline_deriv(u)
    local dtWeights = basis_bspline_deriv(v)
    for i=0,3 do
        for j=0,3 do
            derivU[4*i+j] = dsWeights[j]*tweights[i]
            derivV[4*i+j] = sweights[j]*dtWeights[i]
        end
    end
    return derivU,derivV
end

function data()
    local GDS = Graph("GDS",inc(),
        "d",{M3},inc(),
        "T",{N},inc())

    local cw = Select(lesseq(w_fitSqrt,4*THRESH),1,sqrt(1-w_fitSqrt/4*THRESH))
    cw = Select(lesseq(w_fitSqrt,4*THRESH-1),cw,0)

    SEnergy(cw*(D3(GDS.d)-X3(GDS.T)))
end

function bilinear()
    local GD = Graph("GD",inc(),
        "d",{M3},inc(),
        "u",{N},inc(),
        "m11",{N},inc(),
        "m12",{N},inc(),face_idx[i*ip.npoints+p]
        "m22",{N},inc(),
        "m21",{N},inc()
    )

    local uin = {}
    local umask = {}
    uin[1] = D2(GD.u)(0)
    uin[2] = D2(GD.u)(1)
    umask[1] = 1
    umask[2] = 1
    --local u = softmax(uin,umask)
    local u = {}
    u[1] = uin[1]

    local m11 = X3(GD.m11)
    local m12 = X3(GD.m12)
    local m22 = X3(GD.m22)
    local m21 = X3(GD.m21)
    local bx1 = (1-u[1])*m11 + u[1]*m21
    local bx2 = (1-u[1])*m12 + u[1]*m22;
    local qp = bx1*(1-u[2])+bx2*u[2]; 
    
    local cw = Select(eq(w_fitSqrt,0),1,0)

    SEnergy(cw*(D3(GD.d)-qp), GD) --0.00001
end

function label()
    local params = {"GL", inc(), "d", {M3}, inc(), "m", {M1}, inc(),
		"s", {M1}, inc(), "u",{M2},inc()}
        
    for i=0,basis_max-1 do
        table.insert(params,"T"..i)
        table.insert(params, {N})
        table.insert(params, inc())
    end

    local GD = Graph(unpack(params))

    local u = D2(GD.u)
    
    local us = {}
    local vs = {}
    
    local t = u(0)
    local t2 = t*t
    local t3 = t2*t
    us[0] = 1.0/6.0 * (1.0 - 3.0*(t - t2) - t3)
    us[1] = 1.0/6.0 * (4.0 - 6.0*t2 + 3.0*t3)
    us[2] = 1.0/6.0 * (1.0 + 3.0*(t + t2 - t3))
    us[3] = 1.0/6.0 * (t3)

    t = u(1)
    t2 = t*t
    t3 = t2*t
    vs[0] = 1.0/6.0 * (1.0 - 3.0*(t - t2) - t3)
    vs[1] = 1.0/6.0 * (4.0 - 6.0*t2 + 3.0*t3)
    vs[2] = 1.0/6.0 * (1.0 + 3.0*(t + t2 - t3))
    vs[3] = 1.0/6.0 * (t3)

    local qpbase = 0
    for i=0,3 do
        for j=0,3 do
            qpbase = qpbase + vs[i]*us[j]*X3(GD["T"..(4*i+j)])
        end
    end
    
   --local cw = Select(lesseq(w_fitSqrt,0*THRESH),1,sqrt(1-w_fitSqrt/(4*THRESH)))
   local cw = Select(lesseq(w_fitSqrt,ENABLE_DATA_TERM),1,1.0/8.0)
   cw = Select(lesseq(w_fitSqrt,LABEL_ITER_MUL*THRESH),cw,0)
   cw = sqrt(npoints/nlabels) * cw
   cw = LABEL_TERM*cw
   
   local frac = nlabels/1326.0
   local cr = Select(lesseq(w_fitSqrt,ENABLE_DATA_TERM),0,(7.0/8.0)*frac) --ENABLE_DATA_TERM

   local mask = D1(GD.m)

   local isReverseDataTerm = D1(GD.s)
   cr = MODEL_TO_DATA_TERM*sqrt(cr)

   cw = Select(lesseq(isReverseDataTerm,.5),cw,cr)

   local error = Vector(
       mask*cw*(D3(GD.d)(0)  - qpbase(0)),
       mask*cw*(D3(GD.d)(1)  - qpbase(1)),
       mask*cw*(D3(GD.d)(2)  - qpbase(2))
   )

   Energy(error)
end

function label2D()
	local params = {"G2D", inc(), "d", {M2}, inc(),
		"w",{M1}, inc(),
        "j", {N}, inc(),
        "t", {N}, inc(),
		"s",{M1}, inc(),
		"ct",{M3},inc(),
		"cr",{M3},inc()
        }

	local G2D = Graph(unpack(params)) 

	local s = D1(G2D.s)
	local j = X3(G2D.j) + (1-s)*X3(G2D.t)
	local w = D1(G2D.w)
	local p = D2(G2D.d)
	p = Vector(p(0),p(1),1)

	--local fx = 1498.2242623
	--local fy = 1498.2242623
	--local cx = 790.263706
	--local cy = 578.90334

	--local R = Vector(0.99924469, -0.00488005,  0.03855168,
    --  				-0.01071995, -0.98820433,  0.15276545,
    --  				0.03735144, -0.15306333, -0.98751025)

	local fx = 2*800
	local fy = 2*800
	local cx = 2*395.131
	local cy  =2*310.548

	local C = Vector(1/fx,0,-cx/fx,
                   0,1/fy,-cy/fy,
                   0,0,1)

	--R = Transpose3(R)
	--local ct = Vector(-0.03609917,0.43416458, 2.37101226)
	--local ct = Vector(-0.04783459,0.79178219,3.77646524)
	local ct = D3(G2D.ct)
	local R = Rotate3D( D3(G2D.cr) ) 

	local u = Matrix3x3Mul(C,p)
	--u = Matrix3x3Mul(R,u-ct)
	--local x1 = Matrix3x3Mul(R,-ct)
	local u = Matrix3x3Mul(R,u) + ct
	local x1 = ct 
	local x2 = u
	local x0 = j 
	local nom = Cross(x2-x1,x0-x2)
	local denom = x2-x1
	--local dist = sqrt(Dot3(nom,nom))/sqrt(Dot3(denom,denom))
	local dist = nom/sqrt(Dot3(denom,denom))
	--local cw = Select(lesseq(w_fitSqrt,ENABLE_DATA_TERM),1.0/8.0,0)
	--30*
	
	local cw = Select(lesseq(w_fitSqrt,ENABLE_DATA_TERM),1,1.0/8.0)
    cw = Select(lesseq(w_fitSqrt,LABEL_ITER_MUL*THRESH),cw,0)
    cw = sqrt(npoints/nlabels) * cw
    cw = LABEL_TERM*cw

	Energy(cw*w*dist)
end

function catmulSimple()
    local params = {"GD", inc(), "d", {M3}, inc(),
        "dn", {M3}, inc(),
		"du", {N}, inc(),
        "u",{M2},inc(),
		"su",{M2},inc(),
        "w",{N}, inc()
        }
        
    for i=0,basis_max-1 do
        table.insert(params,"T"..i)
        table.insert(params, {N})
        table.insert(params, inc())
    end

    local GD = Graph(unpack(params))

    local u = D2(GD.u)
    local w = X1(GD.w)    

	local active = disable()
	u = u  + continuous_step_size*active*X2(GD.du)*D2(GD.su) --0.00001

    local us = {}
    local vs = {}
    
    local t = u(0)
    local t2 = t*t
    local t3 = t2*t
    us[0] = 1.0/6.0 * (1.0 - 3.0*(t - t2) - t3)
    us[1] = 1.0/6.0 * (4.0 - 6.0*t2 + 3.0*t3)
    us[2] = 1.0/6.0 * (1.0 + 3.0*(t + t2 - t3))
    us[3] = 1.0/6.0 * (t3)
    
    local dus =  {}
    dus[0] = -0.5*t2 +     t    -.5
    dus[1] =  1.5*t2 - 2.0*t
    dus[2] = -1.5*t2 +     t    +.5
    dus[3] =  0.5*t2

    t = u(1)
    t2 = t*t
    t3 = t2*t
    vs[0] = 1.0/6.0 * (1.0 - 3.0*(t - t2) - t3)
    vs[1] = 1.0/6.0 * (4.0 - 6.0*t2 + 3.0*t3)
    vs[2] = 1.0/6.0 * (1.0 + 3.0*(t + t2 - t3))
    vs[3] = 1.0/6.0 * (t3)
    
    local dvs =  {}
    dvs[0] = -0.5*t2 +     t    -.5
    dvs[1] =  1.5*t2 - 2.0*t
    dvs[2] = -1.5*t2 +     t    +.5
    dvs[3] =  0.5*t2
    

    local derivU
    local derivV
    --derivU,derivV = patch_bspline_deriv(u(0),u(1),us,vs)

    local qpbase = 0
    local du = 0
    local dv = 0
    for i=0,3 do
        for j=0,3 do
            local T = X3(GD["T"..(4*i+j)])
            qpbase = qpbase + us[j]*vs[i]*T
            du = du + dus[j]*vs[i]*T
            dv = dv + us[j]*dvs[i]*T
        end
    end

	local dn = D3(GD.dn)    
    local dp =  (D3(GD.d) - qpbase)

	local n = Vector(
        du(1)*dv(2) - du(2)*dv(1),
        du(2)*dv(0) - du(0)*dv(2),
        du(0)*dv(1) - du(1)*dv(0)
    )
	n = n / sqrt(n(0)*n(0) + n(1)*n(1) + n(2)*n(2))
	du = du / sqrt(du(0)*du(0) + du(1)*du(1) + du(2)*du(2))
	dv = dv / sqrt(dv(0)*dv(0) + dv(1)*dv(1) + dv(2)*dv(2))
	o1 = dn(0)*du(0)+dn(1)*du(1)+dn(2)*du(2)
	o2 = dn(0)*dv(0)+dn(1)*dv(1)+dn(2)*dv(2)

    local frac = 1 --nlabels/npoints
    --local cw = Select(lesseq(w_fitSqrt,ENABLE_DATA_TERM),0,sqrt(1.0/8.0)*frac) --.5
	local cw = Select(lesseq(w_fitSqrt,ENABLE_DATA_TERM),0,sqrt(7.0/8.0)*frac)
    --local cw = Select(lesseq(w_fitSqrt,0*THRESH),0,0)
	cw = .1*Select(lesseq(w_fitSqrt,LABEL_ITER_MUL*THRESH),cw,frac) --7* multiplier???
    --local cw = Select(lesseq(w_fitSqrt,1*THRESH),0,.5*(nlabels/npoints)) --.5
    --cw = Select(lesseq(w_fitSqrt,2*THRESH),cw,.75*nlabels/npoints)
    --cw = Select(lesseq(w_fitSqrt,3*THRESH),cw,(1-.125)*nlabels/npoints)
    --cw = Select(lesseq(w_fitSqrt,4*THRESH),cw,nlabels/npoints)  
    --local cw = nlabels/npoints
    
    local tau = Select(lesseq(w_fitSqrt,70), 100, 50) --100: 100 80
    tau = Select(lesseq(w_fitSqrt,80), tau, 25) --110: 60
	tau = Select(lesseq(w_fitSqrt,90), tau, 10) --120: 40
	tau = Select(lesseq(w_fitSqrt,100), tau, 2) --130: 20
    local wreg = cw*(tau/sqrt(2)) * (w*w - 1)
    --wreg = Select(lesseq(w_fitSqrt,1*THRESH),(w*w - 1),wreg)
    
    local uw = .01 --Select(lesseq(w_fitSqrt,1*THRESH),.1,0)
    
    --cw = Select(lesseq(w_fitSqrt,1*THRESH),cw,sqrt(nlabels/npoints)) --5
    local nw = 1 - active
    local r = dp --Select(lesseq(w_fitSqrt,60),dp,eprime)
	
	
    --sigma = Select(lesseq(w_fitSqrt,2*THRESH),cw,sqrt(.75*nlabels/npoints))
	
    --SEnergy(cw*(r*r/(r*r+sigma*sigma)))
    --,(tau/sqrt(2)) * (w*w - 1)
    --SEnergy(cw*r, GD)
	--local re = sqrt(20*cw)*r*w
    local re = sqrt(20*cw)*r*w
	local ren = sqrt(20*cw)*w*0.01*(dn - n)
	local reu = sqrt(20*cw)*w*0.01*o1
	local rev = sqrt(20*cw)*w*0.01*o2
    --local result = Vector(re(0),re(1),re(2),
	--	0.1*(X2(GD.du)(0)), 0.1*(X2(GD.du)(1)), wreg)
    SEnergy(re)
	--SEnergy(reu)
	--SEnergy(rev)
	--SEnergy(ren)
    
	--SEnergy(re, GD)
    SEnergy(nw*X2(GD.du)) --nw*
    
	SEnergy(wreg)
    --SEnergy(cw*nw*(dn*dn*l2 - n*n))
    --Energy(cw*(tau/sqrt(2)) * (w*w - 1))
    --SEnergy(cw*nw*o1,cw*nw*o2)
    
    --SEnergy(cw*w*(dn - n))
    --SEnergy(cw*w*(Dot3(dn,n)/L_2_norm(n) - 1))
    
    --SEnergy(cw*w*(Dot3(dn,n) - 1))
    
    
    --SEnergy(cw*nw*(dn(0)*n(0)+dn(1)*n(1)+dn(2)*n(2) - 1))
    
    --local unit = 1 - (dn(0)*n(0)+dn(1)*n(1)+dn(2)*n(2))
    --SEnergy(cw*nw*(unit))
end


function meanShapePrior()
    local params = {"GS", inc(), 
        "u1", {N}, inc(),
        "u2", {N}, inc(),
        "u3", {N}, inc(),
        "u4", {N}, inc(),
        "v1",{M3},inc(),
        "v2",{M3},inc(),
        "v3",{M3},inc(),
        "v4",{M3},inc(), 
        "m1",{M1},inc(),
        "m2",{M1},inc(),
        "m3",{M1},inc(),
        "m4",{M1},inc(),
        "mask1",{M3},inc(),
        "mask2",{M3},inc(),
        "mask3",{M3},inc(),
        "mask4",{M3},inc(),
        "m",{M3}, inc()}
    local GS = Graph(unpack(params))

    local u1 = X3(GS.u1)
    local u2 = X3(GS.u2)
    local u3 = X3(GS.u3)
    local u4 = X3(GS.u4)
        
    local mask1 = D3(GS.mask1)
    local mask2 = D3(GS.mask2)
    local mask3 = D3(GS.mask3)
    local mask4 = D3(GS.mask4)

    local v1 = mask1*D3(GS.v1)
    local v2 = mask2*D3(GS.v2)
    local v3 = mask3*D3(GS.v3)
    local v4 = mask4*D3(GS.v4)

    local m1 = D1(GS.m1)
    local m2 = D1(GS.m2)
    local m3 = D1(GS.m3)
    local m4 = D1(GS.m4)
        
    local m = D3(GS.m)(0)

    local cw = 1--Select(eq(v1,v2),SHAPE_PRIOR_PCS,1)

    local cost = m1*((u2-u1)-(v2-v1))
               + m2*((u3-u2)-(v3-v2))
               + m3*((u4-u3)-(v4-v3))
               + m4*((u1-u4)-(v1-v4))

	local cw = Select(lesseq(w_fitSqrt,-1),10*INIT_MEAN_REG,INIT_MEAN_REG)    
    cw = Select(lesseq(w_fitSqrt,50),cw,cw/MEAN_DROPOFF) --REMOVED a 10 factor! 10 -> .2
	cw = Select(lesseq(w_fitSqrt,75),cw,cw/2)

	local cq = Select(lesseq(w_fitSqrt,85),1000*INIT_MEAN_REG,2.5*cw) --10m currently 2.5, now 1.0, .3 
	cw = Select(lesseq(m,1.5),cw,cq)
    --cn = sqrt(m2*cn*nscans) --10c--remove nscans!
    
	cn = (npoints*nscans) / (nfaces*(ncomp+1))
	cn = sqrt(cn)
--cw = Select(lesseq(w_fitSqrt,1*THRESH),cw,sqrt(nlabels/npoints)) --5
    
    Energy(cn*cw*cost)
end

function coeffRegTerm()
    local cn = Select(lesseq(w_fitSqrt,15),6*3,3) --.1
	cn = Select(lesseq(w_fitSqrt,50),cn,.5) --.1
	cn = Select(lesseq(w_fitSqrt,75),cn,.25)
	local GCR = Graph("GCR",inc(),"c",{N}, inc(), "w", {M1},inc())
    wperson = D1(GCR.w)

    local cn = cn*sqrt((npoints*nscans*wperson) / (npersons*ncomp))

	--cn = decay(cn,sqrt(.5),10)
    --m = Select(lesseq(wperson,0.5),1,.1) --20
	--m = sqrt(wperson/30.0)

	SEnergy(COEFF_REG*cn*X1(GCR.c))
end


function ground()
    local params = {"GG", inc(), "u",{M2},inc(),
		"a",{M1},inc(),
		"d",{M1},inc()
	}

    for i=0,basis_max-1 do
        table.insert(params,"T"..i)
        table.insert(params, {N})
        table.insert(params, inc())
    end

    local GD = Graph(unpack(params))

	local a = D1(GD.a)
	local d = D1(GD.d)
    local u = D2(GD.u)

    local us = {}
    local vs = {}

    local t = u(0)
    local t2 = t*t
    local t3 = t2*t
    us[0] = 1.0/6.0 * (1.0 - 3.0*(t - t2) - t3)
    us[1] = 1.0/6.0 * (4.0 - 6.0*t2 + 3.0*t3)
    us[2] = 1.0/6.0 * (1.0 + 3.0*(t + t2 - t3))
    us[3] = 1.0/6.0 * (t3)

    t = u(1)
    t2 = t*t
    t3 = t2*t
    vs[0] = 1.0/6.0 * (1.0 - 3.0*(t - t2) - t3)
    vs[1] = 1.0/6.0 * (4.0 - 6.0*t2 + 3.0*t3)
    vs[2] = 1.0/6.0 * (1.0 + 3.0*(t + t2 - t3))
    vs[3] = 1.0/6.0 * (t3)

    local qpbase = 0
    for i=0,3 do
        for j=0,3 do
            qpbase = qpbase + vs[i]*us[j]*X3(GD["T"..(4*i+j)])
        end
    end


   local frac = nlabels/1326.0
   local ground = d   --  -.984

   local groundAxis = qpbase(0)
   groundAxis = Select(lesseq(a,0.5),groundAxis, qpbase(1))
   groundAxis = Select(lesseq(a,1.5),groundAxis, qpbase(2))

   ground = groundAxis-ground
   ground = Select(lesseq(ground,0),ground,0)
   ground = ground*Select(lesseq(w_fitSqrt,ENABLE_GROUND),0*frac,100*frac) --100
   local error = Vector(
       ground
   )
   Energy(error)
end



function pjointTerm()
    local params = {"GPJ", inc(), "PJ",{N},inc(), "meanJoint", {K3}, inc()}
    for j=1,ncomp do
        table.insert(params,"PC"..j)
        table.insert(params, {K3})
        table.insert(params, inc())

        table.insert(params,"c"..j)
        table.insert(params, {N})
        table.insert(params, inc())
    end

    local GPJ = Graph(unpack(params))

    local joint = Y3(GPJ["meanJoint"])
    for i=1,ncomp do
        joint = joint + X1(GPJ['c'..i]) * Y3(GPJ["PC"..i]) 
    end
	local cw = sqrt((nscans*npoints)/(njoints*npersons))
	Energy(cw*(X3(GPJ["PJ"]) - joint)) 
    --SEnergy(5*lambda_sparse*sqrt(nscans)*(X3(GPJ["PJ"]) - joint)) --#nscans 50
end

function weightPriorTerm()
    local GWP = Graph("GWP", inc(),"w",{N},inc(),"wp",{M1},inc())
    local w = X1(GWP.w) --/ lambda
    local wp = D1(GWP.wp) / lambda

    w = .5 + .5*w / sqrt(1 + w*w)
    --local wreg = {}
    --local wmask = {}
    --wmask[1] = 1
    --wreg[1] = w
    --local myw = softmax(wreg,wmask)[1]
    --local cn = Select(lesseq(w_fitSqrt,4*THRESH),1000,.0001)
	local cn = (nscans*npoints)/(2*nverts)
    cn = cn*Select(lesseq(w_fitSqrt,85),1,.00000001) --0.0001
    SEnergy(10*cn*(w-wp))
end

function weightNormTerm()
    local params = {"GW", inc()}
    for j=1,4 do
        table.insert(params,"w"..j)
        table.insert(params, {N})
        table.insert(params, inc())
    end
    local GW = Graph(unpack(params))

    local EW = -1
    local wreg = {}
    local wmask = {}
    for j=1,4 do
        wreg[j] = X1(GW["w"..j]) / lambda
        wmask[j] = 1
    end
    local myw = softmax(wreg,wmask)
    for j=1,4 do
        EW = EW + myw[j]
    end
	local cn = (nscans*npoints)/(nverts)
    --local cn = Select(lesseq(w_fitSqrt,3*THRESH),1000,1000)
    SEnergy(10*cn*EW) --1000
end

function jointRingPriorTerm()
	local params = {"GJT", inc(), "d", {N}, inc()}
	for j=1,MAX_VERTEX_RING do
		table.insert(params,"j"..j)
        table.insert(params, {N})
        table.insert(params, inc())

		table.insert(params,"w"..j)
        table.insert(params, {M1})
        table.insert(params, inc())
	end

	local GJT = Graph(unpack(params))

	local e = Vector(0,0,0)
	local w = 0
	for i=1,MAX_VERTEX_RING do
		e = e + D1(GJT["w"..i]) * X3(GJT["j"..i])
		w = w + D1(GJT["w"..i])
	end

	e = e / w
	e = e - X3(GJT.d)

	--cn = Select(lesseq(w_fitSqrt,2*THRESH),1,.1)
	cn = 5*sqrt((nscans*npoints)/((ncomp+1)*njoints))
	--SEnergy(sqrt(nscans)*cn*40*e) --nscans
	cn = Select(lesseq(w_fitSqrt,REDUCE_POSE_REG),cn,0.1*cn) 
	Energy(cn*e)
end

function jointPriorTerm()
     local GJT = Graph("GJT", inc(),
        "jointPrior",{M3},inc(),
        "joint",{N},inc(),
        "vertexPrior",{M3},inc(),
        "vertex",{N},inc()
    )

    local jointPrior = D3(GJT.jointPrior)
    local joint = X3(GJT.joint)
    local vertexPrior = D3(GJT.vertexPrior)
    
    local cw = Select(eq(jointPrior,vertexPrior),SHAPE_PRIOR_PCS,1)

    local vertex = X3(GJT.vertex)
    local cost = (vertexPrior-jointPrior) - (vertex - joint)

    local cn = 1 * cw * Select(lesseq(w_fitSqrt,2*THRESH),1,.1) --.1 --10
    
    SEnergy(sqrt(nscans)*cn*cost) --nscans
end

function advancedPosePrior()

	local params = {"GAP", inc()}
	
	table.insert(params,"meanbase")
    table.insert(params, {M3})
    table.insert(params, inc())

	if TRAIN then 
    	table.insert(params,"mean")
    	table.insert(params, {N})
    	table.insert(params, inc())
	end


    table.insert(params,"thetabase")
    table.insert(params, {M3})
    table.insert(params, inc())
    table.insert(params,"theta")
	table.insert(params, {N})
	table.insert(params, inc())
    
    
	local G = Graph(unpack(params))
	
	local thetabase = D3(G.thetabase)
	local meanbase = D3(G.meanbase)
	
	local mean = Vector(0,0,0)
	if TRAIN then
		mean = X3(G["mean"])
	end
    local theta = X3(G["theta"])
        
	local R_a = Matrix3x3Dot(Rotate3D(thetabase), Rod(theta))
	local R_b = Matrix3x3Dot(Rotate3D(meanbase), Rod(mean))
	local EAnim = R_a - R_b --Vector(1,0,0,0,1,0,0,0,1)
	local Eprior = R_b - Vector(1,0,0,0,1,0,0,0,1)
	--local m = Select(lesseq(w_fitSqrt,7),1000,POSE_REG)	
	local m = POSE_REG
	local cn = sqrt((nscans*npoints)/(nscans*njoints))
	cn = Select(lesseq(w_fitSqrt,ENABLE_DATA_TERM),cn,cn/3.0)
	cn = Select(lesseq(w_fitSqrt,REDUCE_POSE_REG),cn,cn/1000.0) --cn,0
	SEnergy(m*cn*EAnim) --0.01 .05
	--SEnergy((1.0/sqrt(nscans))*m*Eprior)
	--SEnergy(m*Eprior)
end


function centerTerm()
    local GC = Graph("GC", inc(),"w",{N},inc())
    local w = X3(GC.w)
    SEnergy(1000*(w))
end

function shapeRegTerm()
    local GSRT = Graph("GSRT", inc(),"w",{N},inc())
    local w = X3(GSRT.w)
    SEnergy(0.0001*(w))
end


function localTerm()
    local params = {"GBA", inc(),"d",{N},inc(),"w",{MW},inc()}
    for j=0,(local_size-1) do
        table.insert(params,"T"..j)
        table.insert(params, {N})
        table.insert(params, inc())
    end
    local GBA = Graph(unpack(params))

    local resi = -X3(GBA.d)
    local w = DW(GBA.w)
    --table.insert(resi,XP(GBA.b))
    for j=0,(local_size-1) do
        resi = resi+w(j)*X3(GBA["T"..j])
    end
    --local cn = Select(lesseq(w_fitSqrt, THRESH),1,10)
    --SEnergy(lambda_sparse * resi)
	local cn = 1*sqrt((nscans*npoints)/(nscans*nlocal))
	Energy(lambda_sparse * cn*resi)
end

function temporalPriorTerm()
	local params = {"GTP", inc()}
	
	--table.insert(params,"trans_a")
    --table.insert(params, {N})
    --table.insert(params, inc())

    --for i=1,njoints do
    table.insert(params,"theta_a")
    table.insert(params, {N})
    table.insert(params, inc())
    --end
    
    table.insert(params,"thetabase_a")
    table.insert(params, {M3})
    table.insert(params, inc())
    
    --table.insert(params,"trans_b")
    --table.insert(params, {N})
    --table.insert(params, inc())

    --for i=1,njoints do
    table.insert(params,"theta_b")
    table.insert(params, {N})
    table.insert(params, inc())
    --end
    
    table.insert(params,"thetabase_b")
    table.insert(params, {M3})
    table.insert(params, inc())
    
    local GTP = Graph(unpack(params))
    
	local cn = sqrt((nscans*npoints)/(nscans*njoints))
    local tw = 0.000001
	tw = Select(lesseq(w_fitSqrt,12),tw,5)
	tw = Select(lesseq(w_fitSqrt,REDUCE_POSE_REG),tw,tw/12) --/12 48
	tw = tw*cn*POSE_REG
    --tw = decay(tw,sqrt(.5),20)
    --if i == 1 then
    --    Energy(tw*(X3(GTP.trans_a)-X3(GTP.trans_b)))
	--end
    local thetabase_a = D3(GTP.thetabase_a)
    local thetabase_b = D3(GTP.thetabase_b)
    
    
    --for i=1,njoints do
	local theta_a = X3(GTP.theta_a)
	local theta_b = X3(GTP.theta_b)
	--local thetabase = D3(G["thetabase"..i])
	--local thetabase = Slice(thetabases,3*i,3*(i+1))
	local relA = Matrix3x3Dot(Rotate3D(thetabase_a), Rod(theta_a))
	local relB = Matrix3x3Dot(Rotate3D(thetabase_b), Rod(theta_b))
	--local relB = Rotate3D(-thetabase_b)
	--local c = Matrix3x3Dot(relA, relB)
	
	--c = Matrix3x3Dot(c,Rod(theta_b))
	--local EAnim = c - Vector(1,0,0,0,1,0,0,0,1)
	Energy(tw*(relA-relB))
	--relativesR[i] = Rod(theta)1
    --end
end

function symmetryTerm()
	local GSym = Graph("GSym",inc(),
    "u",{N},inc(),
    "v",{N},inc())
    
    local u = X3(GSym.u)
    local v = X3(GSym.v)

	local cn = 1--Select(lesseq(w_fitSqrt,85),1,.00001)
    --local cn = Select(lesseq(w_fitSqrt,2.5*THRESH),100,10)
    local wsym = cn*10*sqrt((nscans*npoints)/((ncomp+1)*nverts))
    Energy(Vector(
		wsym*(u(0)+v(0)),
		wsym*(u(1)-v(1)),
		wsym*(u(2)-v(2))
	))
    
end

function vertex_vars(params)
    table.insert(params,"v")
    table.insert(params, {N})
    table.insert(params, inc())

    for i=1,4 do
        table.insert(params, "weight"..i)
        table.insert(params, {K1})
        table.insert(params, inc())
    end

    
end

function vertex_vars_rest(params)
    table.insert(params,"v")
    table.insert(params, {K3})
    table.insert(params, inc())

    for i=1,ncomp do 
        table.insert(params, "PC"..i)
        table.insert(params, {K3})
        table.insert(params, inc())
    end

    for i=1,ndeformjoints do
    for j=1,ndeformshapes do
        table.insert(params,"def"..i.."d"..j)
        table.insert(params, {K3})
        table.insert(params, inc())
    end
    end

    for i=1,ndeformjoints do
	table.insert(params,"influence_base"..i)
        table.insert(params, {M3})
        table.insert(params, inc())
    
        table.insert(params,"influence"..i)
        table.insert(params, {N})
        table.insert(params, inc())
        
        table.insert(params,"ellbow"..i)
        table.insert(params, {M3})
        table.insert(params, inc())
    end
end

function transform_basis(GR,resultR, resultT)
    --M = X3(GR["def".."1".."d".."1"])
    local M = X3(GR["v"]) 
    --M = Vector(0,0,0)

    
    local win = {}
    local mask = {}
    for i=1,4 do
        win[i] = Y1(GR["weight"..i]) / lambda
        --win[i] = 1
    end   
    
    local wout = softmax(win)    
    --local skinned = Vector(0,0,0)
     
    local transformR = Vector(0,0,0,0,0,0,0,0,0)
    local transformT = Vector(0,0,0)
    
    for i=1,4 do
        local w = wout[i]
        transformR = transformR + w*resultR[i]
        transformT = transformT + w*resultT[i]
    end
    local skinned = Matrix3x3Mul(transformR,M) + transformT
    skinned = skinned + X3(GR["trans"])
    return skinned
end

function transform_rest_pose(GR)
    local M = Y3(GR["v"]) 

    for i=1,ncomp do  
        M = M + X1(GR['c'..i])*Y3(GR["PC"..i])
    end

    for i=1,ndeformjoints do
		local thetabase = D3(GR["influence_base"..i])
        local theta = X3(GR["influence"..i])
        local ellbow = D3(GR["ellbow"..i])(0)
        SEnergy((1-ellbow)*theta)
        theta = ellbow*theta
        --theta[ZERO_JOINT] = theta[ZERO_JOINT]*ellbow
        
        --theta[0] = theta[0]*ellbow
        --theta[1] = theta[1]*ellbow
        --theta[2] = theta[2]*ellbow
        
        
        local R = Matrix3x3Dot(Rotate3D(thetabase), Rod(theta))
		R = R - Vector(1,0,0,0,1,0,0,0,1)
    for j=1,ndeformshapes do
        M = M + R(j-1)*Y3(GR["def"..i.."d"..j])
    end
    end
    
    return M
end



function anc(name)
 for i=1,njoints do
    if ancestors[i][1] == name then
        return ancestors[i][2]
    end
 end
end

function jidx(name)
  for i=1,njoints do
    if ancestors[i][1] == name then
        return i
    end
 end
end

function modelTerm(name)
    --local relativesT = {}
    --local fixedJointTransformT = {}
    local resultR = {}
    local resultT = {}
    local joints = {}
    local posedJoints = {}
    local params = {name, inc(), "d", {N}, inc()}

    vertex_vars(params)

    for i=1,4 do
    table.insert(params,"meanJoint"..i)
    table.insert(params, {N})
    table.insert(params, inc())
    end

    for i=1,4 do
    table.insert(params,"posedJoint"..i)
    table.insert(params, {N})
    table.insert(params, inc())
    end
    
    table.insert(params,"trans")
    table.insert(params, {N})
    table.insert(params, inc())

    for i=1,4 do
    table.insert(params,"theta"..i)
    table.insert(params, {N})
    table.insert(params, inc())
    end
    for i=1,4 do
    table.insert(params,"thetabase"..i)
    table.insert(params, {M3})
    table.insert(params, inc())
    end
    local G = Graph(unpack(params))
    pretty(params)
    --print(G.thetabase)
    
    UsePreconditioner(true)

    for i=1,4 do
       table.insert(joints,X3(G["meanJoint"..i]))
    end

    for i=1,4 do
       table.insert(posedJoints,X3(G["posedJoint"..i]))
    end

    local absolutesR = {}

    thetabases = {}
    for i=1,4 do
       table.insert(thetabases,D3(G["thetabase"..i]) )
    end
    
    ----local garbage = 0
    local j = 1
    for i=1,4 do
        
        local k = i - 1
        ----print("joint id: ",j)
        local thetabase = Vector(thetabases[i](0),
                                 thetabases[i](1),
                                 thetabases[i](2))
        
        --if isRigid(ancestors[i][1]) then
	    -- TODO: THIS IS INCORRECT WITH THE NEW SETUP!!
       
        print("var joint id: ", j, " real joint ", i)
        local theta = X3(G["theta"..j])
        --relativesR[i] = Matrix3x3Dot(Rotate3D(thetabase), Rod(theta))
        absolutesR[i] = Matrix3x3Dot(Rotate3D(thetabase), Rod(theta))
        j = j + 1     
    end

    local absolutesT = {}
    local absolutesFixedT = {}

    --absolutesR[1] = relativesR[1]
    --absolutesT[1] = relativesT[1]
    absolutesT[1] = posedJoints[1]
    --absolutesFixedR[1] = Rotate3D(Vector(0,0,0))
    absolutesFixedT[1] = joints[1]

    for i=1,4 do
        --absolutesR[i] = Matrix3x3Dot(absolutesR[j],relativesR[i])
        --absolutesT[i] = Matrix3x3Mul(absolutesR[j],relativesT[i]) + absolutesT[j]
	absolutesT[i] = posedJoints[i]

        absolutesFixedT[i] = joints[i]
        --absolutesFixedT[i] = fixedJointTransformT[i] + absolutesFixedT[j]
    end

    for i=1,4 do
        resultR[i] = absolutesR[i]
        resultT[i] = Matrix3x3Mul(-resultR[i], absolutesFixedT[i]) + absolutesT[i]
    end

    local output = Vector(0,0,0)
    output = output + transform_basis(G,resultR, resultT)
    local modelError = X3(G.d)-output
    --SEnergy(D3(G.d)-output) -- +X3(G["trans"]) OLD!!
    --SEnergy(sqrt(1) * modelError) -- +X3(G["trans"])
    
    --local cn = (2*w_fitSqrt/THRESH + 1)
    --cw = Select(lesseq(w_fitSqrt,2*THRESH),cw,sqrt(.01)) --5!
    --cw = Select(lesseq(w_fitSqrt,3*THRESH),cw,sqrt(.01)) --5!
    --cw = Select(lesseq(w_fitSqrt,4*THRESH),cw,sqrt(0.1)) --5!
    --cw = Select(lesseq(w_fitSqrt,4*THRESH),cw,sqrt(1))
    --local cn = Select(lesseq(w_fitSqrt,4*THRESH),1,5)
    --cn = Select(lesseq(w_fitSqrt,6*THRESH),cn,10)  

	local cn = sqrt((nscans*npoints)/(nscans*nverts))
    SEnergy(cn*lambda_sparse*modelError)    
end



function restModelTerm(name)
	local relativesT = {}
	local fixedJointTransformT = {}
	local resultR = {}
	local resultT = {}
	local joints = {}

    local params = {name, inc(), "d", {N}, inc()}


    vertex_vars_rest(params)

    
    for i=1,ncomp do 
    table.insert(params,"c"..i)
    table.insert(params, {N})
    table.insert(params, inc())
    end
    
    table.insert(params,"r")
    table.insert(params, {M3})
    table.insert(params, inc())

    local G = Graph(unpack(params))
    --for index, value in ipairs(params) do
    --    print(index, value)
    --end
    UsePreconditioner(true)

    local output = transform_rest_pose(G)
    local modelError = X3(G.d)-output
    
    --local cn = 1 --Select(lesseq(D3(G.r)(0),.5),1,100) --.1
    --SEnergy(lambda_sparse*modelError)
	local cn = sqrt((nscans*npoints)/(nscans*nverts))
	Energy(lambda_sparse*REST_MODEL_REG*modelError)    
end

function relativeToAbsoluteTerm()
     local params = {"GRTA", inc()}

    table.insert(params,"theta_a")
    table.insert(params, {N})
    table.insert(params, inc())

    table.insert(params,"thetabase_a")
    table.insert(params, {M3})
    table.insert(params, inc())

    table.insert(params,"theta_r")
    table.insert(params, {N})
    table.insert(params, inc())

    table.insert(params,"thetabase_r")
    table.insert(params, {M3})
    table.insert(params, inc())

    table.insert(params,"theta_b")
    table.insert(params, {N})
    table.insert(params, inc())

    table.insert(params,"thetabase_b")
    table.insert(params, {M3})
    table.insert(params, inc())

    table.insert(params,"j_a")
    table.insert(params, {N})
    table.insert(params, inc())

    table.insert(params,"j_b")
    table.insert(params, {N})
    table.insert(params, inc())

    table.insert(params,"J_b")
    table.insert(params, {N})
    table.insert(params, inc())

    table.insert(params,"J_0")
    table.insert(params, {N})
    table.insert(params, inc())

    table.insert(params,"root")
    table.insert(params, {M3})
    table.insert(params, inc())

    table.insert(params,"rigid")
    table.insert(params, {M3})
    table.insert(params, inc())

    local G = Graph(unpack(params))
    local thetabase_a = D3(G.thetabase_a)
    local thetabase_b = D3(G.thetabase_b)
    local thetabase_r = D3(G.thetabase_r)
    local root = D3(G.root)(0)
    local rigid = D3(G.rigid)(0)

    local theta_a = X3(G.theta_a)
    local theta_b = X3(G.theta_b)
    local theta_r = (1-rigid)*X3(G.theta_r)
    local j_a = X3(G.j_a)
    local j_b = X3(G.j_b)
    local J_b = X3(G.J_b)
    local J_0 = X3(G.J_0)

    local R_a = Matrix3x3Dot(Rotate3D(thetabase_a), Rod(theta_a))
    R_a = root*Vector(1,0,0,0,1,0,0,0,1) + (1-root)*R_a
    local R_b = Matrix3x3Dot(Rotate3D(thetabase_b), Rod(theta_b))
    local R_r = Matrix3x3Dot(Rotate3D(thetabase_r), Rod(theta_r))
    local R_c = Matrix3x3Dot(R_a,R_r)
    local diff = R_b - R_c

	local cn = sqrt((nscans*npoints)/(nscans*njoints))
    Energy(cn * lambda_sparse*diff)

    J_0 = (1-root)*J_0
    j_a = (1-root)*j_a
    local jnew = Matrix3x3Mul(R_a,J_b-J_0) + j_a
    local jdiff = jnew - j_b

	cn = sqrt((nscans*npoints)/(nscans*njoints))
    Energy(cn * lambda_sparse * jdiff)
end

function linkVariableTerm()
    local params = {"GLV", inc()}

    table.insert(params,"a")
    table.insert(params, {N})
    table.insert(params, inc())

    table.insert(params,"b")
    table.insert(params, {N})
    table.insert(params, inc())

    local G = Graph(unpack(params))
    local a = X3(G.a)
    local b = X3(G.b)

	local cn = sqrt((nscans*npoints)/(nscans))
    Energy(cn * lambda_sparse*(a-b))    
end


function filler1()
	local params = {"GJLT", inc(), "P", {MP}, inc()}
end

catmulSimple()
label()
label2D()
restModelTerm("G")
modelTerm("G2")
pjointTerm()
if TRAIN then
	weightPriorTerm()
	weightNormTerm()
	meanShapePrior()
end
coeffRegTerm()
ground()
if TRAIN then
	jointRingPriorTerm()
	--jointPriorTerm()
	centerTerm()
end
localTerm()
if USE_TEMPORAL > 0 then
	temporalPriorTerm()
end
if TRAIN and USE_SYMMETRY>0 then
	symmetryTerm()
end
advancedPosePrior()
relativeToAbsoluteTerm()
linkVariableTerm()




--[[
local GAnim = Graph("GAnim",inc(),
    "theta_a",{N},inc(),
    "base_a",{N},inc(),
    "theta_b",{N},inc(),
    "base_b",{N},inc())

local theta_a = X3(GAnim["theta_a"])
local base_a = D3(GAnim["base_a"])
local theta_b = X3(GAnim["theta_b"])
local base_b = D3(GAnim["base_b"])

local R_a = Matrix3x3Dot(Rotate3D(base_a), Rod(theta_a))
local R_b = Matrix3x3Dot(Rotate3D(base_b), Rod(theta_b))

local R_diff = Matrix3x3Dot(Transpose(R_a),R_b)
local EAnim = R_diff - Vector(1,0,0,0,1,0,0,0,1)
SEnergy(EAnim)
]]--


print("nparams in lua: ",cfg)
print("bytes used before gc: ", collectgarbage("count")/1024)
collectgarbage("collect")
print("bytes used after gc: ", collectgarbage("count")/1024)

