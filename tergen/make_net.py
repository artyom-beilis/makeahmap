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

def mynet(batch,steps,loss_type,dep=False,descr=False,part='gen'):

    conv_lr = [dict(lr_mult=1,decay_mult=1),dict(lr_mult=2,decay_mult=1)]
    scale_lr = [dict(lr_mult=1,decay_mult=1),dict(lr_mult=1,decay_mult=1)]
    bn_param = dict(eps=0.001,use_global_stats=False)

    fr_lr = [dict(lr_mult=0,decay_mult=0),dict(lr_mult=0,decay_mult=0)]
    fr_bn = dict(eps=0.001,use_global_stats=True)

    if part=='gen':
        gen_conv_lr = conv_lr
        gen_scale_lr = scale_lr
        gen_bn_param = bn_param
        dsc_conv_lr = fr_lr
    else:
        gen_conv_lr = fr_lr
        gen_scale_lr = fr_lr
        gen_bn_param = fr_bn
        dsc_conv_lr = conv_lr

    n=caffe.NetSpec()

    sp=dict(bias_term=True,filler=dict(value=1.0))

    if dep:
        n.source = L.Input(input_param=dict(shape=[dict(dim=[1,1,64,64])]))
    else:
        if descr:
            if part == 'gen':
                n.data,n.label = L.Data(data_param=dict(source="db",batch_size=batch,backend=P.Data.LMDB),ntop=2)
            else:
                n.data,n.label = L.Data(data_param=dict(source="db",batch_size=batch/2,backend=P.Data.LMDB),ntop=2)
        else:
            n.data = L.Data(data_param=dict(source="db",batch_size=batch,backend=P.Data.LMDB))
        
        n.expected, n.source = L.Slice(n.data,slice_param=dict(axis=1,slice_point=1),ntop=2)
        if part!='gen':
            #n.data_ref,n.label_ref = L.Data(data_param=dict(source="db_ref",batch_size=batch/2,backend=P.Data.LMDB),ntop=2)
            n.data_ref = L.Split(n.expected)
            n.label_ref = L.DummyData(shape=[dict(dim=[batch/2])],data_filler=dict(value=1.0))

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

    outname="topconv"
    n[outname]=L.Convolution(n[inp],convolution_param=conv_param(5,1),param=gen_conv_lr)
    n.generated=L.Sigmoid(n.topconv)
    if not dep:
        lw = 1 if part == 'gen' else 0
        if loss_type == 'euc':
            n.l2_loss=L.EuclideanLoss(n.generated,n.expected,name="loss",loss_weight=lw)
        else:
            n.l2_loss=L.EuclideanLoss(n.generated,n.expected,name="loss",loss_weight=0)
            n.cross_entropy_loss=L.SigmoidCrossEntropyLoss(n.topconv,n.expected,name="loss",loss_weight=lw)
    if descr:
        if part!='gen':
            n.desc_inp = L.Concat(n[outname],n.data_ref,concat_param=dict(axis=0))
            #n.desc_lbl = L.Concat(n.label,n.label_ref,concat_param=dict(axis=0))
            n.desc_lbl = L.Concat(n.label_ref,n.label,concat_param=dict(axis=0))
            outname = "desc_inp"
            label_name = "desc_lbl"
        else:
            label_name = 'label'
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
        n.d_ip2 = L.InnerProduct(n.d_ip1,param=dsc_conv_lr,inner_product_param=ip_param(1))

            
        n.sigmoid_loss = L.SigmoidCrossEntropyLoss(n.d_ip2,n[label_name],name="loss",loss_weight=100)
        n.score = L.Sigmoid(n.d_ip2)
        n.lbl_flat=L.Reshape(n[label_name],reshape_param=dict(shape=dict(dim=[-1,1])))
        n.diff = L.Eltwise(n.score,n.lbl_flat,eltwise_param=dict(coeff=[1.0/batch,-1.0/batch]))
        n.error = L.Reduction(n.diff,reduction_param=dict(operation=P.Reduction.ASUM))
        #n.score_out = L.Split(n.score)

    return n

def print_to_file(name,n):
    s=str(n.to_proto())
    f=open(name,'w')
    f.write(s)
    f.close()


batch=128
res_steps=2
loss='euc'
#loss='cross_entropy'

print_to_file('gan_dep.prototxt',       mynet(batch,res_steps,loss,True))
print_to_file('gan_pretrain.prototxt',  mynet(batch,res_steps,loss,False))
print_to_file('gan_train_gen.prototxt', mynet(batch,res_steps,loss,False,True,'gen'))
print_to_file('gan_train_dsc.prototxt', mynet(batch,res_steps,loss,False,True,'dsc'))

