import paths
import sys

for p in paths.external:
    sys.path.append(p)

import caffe
from caffe import layers as L
from caffe import params as P

def conv_param(k,n):
    return  dict(   num_output=n,
                    pad=k/2,
                    kernel_size=k,
                    stride=1,
                    weight_filler=dict(type="xavier"),
                    bias_filler=dict(type="constant"))
    

def mynet(batch,steps,dep=False):
    n=caffe.NetSpec()

    sp=dict(bias_term=True,filler=dict(value=1.0))

    if dep:
        n.source = L.Input(input_param=dict(shape=[dict(dim=[1,1,128,128])]))
    else:
        n.data = L.Data(data_param=dict(source="db",batch_size=batch,backend=P.Data.LMDB))
        n.expected, n.source = L.Slice(n.data,slice_param=dict(axis=1,slice_point=1),ntop=2)

    n.conv1 = L.Convolution(n.source,convolution_param=conv_param(5,16),param=[dict(lr_mult=1),dict(lr_mult=2)])
    n.bn1   = L.BatchNorm(n.conv1,batch_norm_param=dict(eps=0.001))
    n.scale1= L.Scale(n.bn1,scale_param=sp)
    n.scale1= L.ReLU(n.scale1)
    inp = "scale1"
    for m in range(steps):
        k=m+1
        cid1="step%d/conv1" % k
        cid2="step%d/conv2" % k
        bid1="step%d/bn1"   % k
        bid2="step%d/bn2"   % k
        eid ="step%d/elt"   % k
        
        n[cid1] = L.Convolution(n[inp],convolution_param=conv_param(3,16),param=[dict(lr_mult=1),dict(lr_mult=2)])
        n[bid1] = L.BatchNorm(n[cid1],batch_norm_param=dict(eps=0.001))
        n[bid1] = L.Scale(n[bid1],scale_param=sp)
        n[bid1] = L.ReLU(n[bid1])

        n[cid2] = L.Convolution(n[bid1],convolution_param=conv_param(3,16),param=[dict(lr_mult=1),dict(lr_mult=2)])
        n[bid2] = L.BatchNorm(n[cid2],batch_norm_param=dict(eps=0.001))
        n[bid2] = L.Scale(n[bid2],scale_param=sp)
        n[bid2] = L.ReLU(n[bid2])

        n[eid]  = L.Eltwise(n[bid2],n[inp])
        inp = eid

    if dep:
        outname="generated"
    else:
        outname="topconv"
    n[outname]=L.Convolution(n[inp],convolution_param=conv_param(5,1),param=[dict(lr_mult=1),dict(lr_mult=2)])
    if not dep:
        #n.cross_entropy_loss=L.SigmoidCrossEntropyLoss(n.topconv,n.expected,name="loss",loss_weight=1)
        #n.generated=L.Sigmoid(n.topconv)
        n.l2_loss=L.EuclideanLoss(n.expected,n.topconv,name="loss",loss_weight=1)
    return n


ns=mynet(30,3,True)
print ns.to_proto()

