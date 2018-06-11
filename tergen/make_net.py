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

def ip_param(n):    
    return  dict(   num_output=n,
                    weight_filler=dict(type="xavier"),
                    bias_filler=dict(type="constant"))

def mynet(batch,steps,dep=False,descr=False):

    conv_lr = [dict(lr_mult=1,decay_mult=1),dict(lr_mult=2,decay_mult=1)]
    scale_lr = [dict(lr_mult=1,decay_mult=1),dict(lr_mult=1,decay_mult=1)]
    bn_param = dict(eps=0.001,use_global_stats=False)

    fr_lr = [dict(lr_mult=0,decay_mult=0),dict(lr_mult=0,decay_mult=0)]
    fr_bn = dict(eps=0.001,use_global_stats=True)

    gen_conv_lr = conv_lr
    gen_scale_lr = scale_lr
    gen_bn_param = bn_param

    dsc_conv_lr = fr_lr

    n=caffe.NetSpec()

    sp=dict(bias_term=True,filler=dict(value=1.0))

    if dep:
        n.source = L.Input(input_param=dict(shape=[dict(dim=[1,1,128,128])]))
    else:
        if descr:
            n.data,n.label = L.Data(data_param=dict(source="db",batch_size=batch,backend=P.Data.LMDB),ntop=2)
        else:
            n.data = L.Data(data_param=dict(source="db",batch_size=batch,backend=P.Data.LMDB))
        n.expected, n.source = L.Slice(n.data,slice_param=dict(axis=1,slice_point=1),ntop=2)

    n.conv1 = L.Convolution(n.source,convolution_param=conv_param(5,16),param=gen_conv_lr)
    n.bn1   = L.BatchNorm(n.conv1,batch_norm_param=gen_bn_param)
    n.scale1= L.Scale(n.bn1,scale_param=sp,param=gen_scale_lr)
    n.scale1= L.ReLU(n.scale1)
    inp = "scale1"
    for m in range(steps):
        k=m+1
        cid1="step%d/conv1" % k
        cid2="step%d/conv2" % k
        bid1="step%d/bn1"   % k
        bid2="step%d/bn2"   % k
        eid ="step%d/elt"   % k
        
        n[cid1] = L.Convolution(n[inp],convolution_param=conv_param(3,16),param=gen_conv_lr)
        n[bid1] = L.BatchNorm(n[cid1],batch_norm_param=gen_bn_param)
        n[bid1] = L.Scale(n[bid1],scale_param=sp,param=gen_scale_lr)
        n[bid1] = L.ReLU(n[bid1])

        n[cid2] = L.Convolution(n[bid1],convolution_param=conv_param(3,16),param=gen_conv_lr)
        n[bid2] = L.BatchNorm(n[cid2],batch_norm_param=gen_bn_param)
        n[bid2] = L.Scale(n[bid2],scale_param=sp,param=gen_scale_lr)
        n[bid2] = L.ReLU(n[bid2])

        n[eid]  = L.Eltwise(n[bid2],n[inp])
        inp = eid

    if dep:
        outname="generated"
    else:
        outname="topconv"
    n[outname]=L.Convolution(n[inp],convolution_param=conv_param(5,1),param=gen_conv_lr)
    if not dep:
        #n.cross_entropy_loss=L.SigmoidCrossEntropyLoss(n.topconv,n.expected,name="loss",loss_weight=1)
        n.generated=L.Sigmoid(n.topconv)
        n.l2_loss=L.EuclideanLoss(n.expected,n.generated,name="loss",loss_weight=1)
    if descr:
        n.d_conv1 = L.Convolution(n[outname],convolution_param=conv_param(5,32),param=dsc_conv_lr) 
        n.d_pool1 = L.Pooling(n.d_conv1,pooling_param=dict(kernel_size=3,stride=2,pool=P.Pooling.MAX)) 
        n.d_pool1 = L.ReLU(n.d_pool1)
        
        n.d_conv2 = L.Convolution(n.d_pool1,convolution_param=conv_param(5,32),param=dsc_conv_lr) 
        n.d_pool2 = L.Pooling(n.d_conv2,pooling_param=dict(kernel_size=3,stride=2,pool=P.Pooling.MAX)) 
        n.d_pool2 = L.ReLU(n.d_pool2)
        
        n.d_conv3 = L.Convolution(n.d_pool2,convolution_param=conv_param(5,64),param=dsc_conv_lr) 
        n.d_pool3 = L.Pooling(n.d_conv3,pooling_param=dict(kernel_size=3,stride=2,pool=P.Pooling.MAX)) 
        n.d_pool3 = L.ReLU(n.d_pool3)
        
        n.d_conv4 = L.Convolution(n.d_pool3,convolution_param=conv_param(5,64),param=dsc_conv_lr) 
        n.d_pool4 = L.Pooling(n.d_conv4,pooling_param=dict(kernel_size=3,stride=2,pool=P.Pooling.MAX)) 
        n.d_pool4 = L.ReLU(n.d_pool4)

        n.d_ip1 = L.InnerProduct(n.d_pool4,param=dsc_conv_lr,inner_product_param=ip_param(512))
        n.d_ip1 = L.ReLU(n.d_ip1)
        n.d_ip2 = L.InnerProduct(n.d_ip1,param=dsc_conv_lr,inner_product_param=ip_param(2))
        n.softmax_with_loss = L.SoftmaxWithLoss(n.d_ip2,n.label,name="loss",loss_weight=0)

    return n


ns=mynet(16,3,False,False)
print ns.to_proto()

