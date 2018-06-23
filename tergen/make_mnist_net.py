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

def conv_param_nb(k,n):
    return  dict(   num_output=n,
                    bias_term=False,
                    pad=k/2,
                    kernel_size=k,
                    stride=1,
                    weight_filler=dict(type="xavier"))

def ip_param(n):    
    return  dict(   num_output=n,
                    weight_filler=dict(type="xavier"),
                    bias_filler=dict(type="constant"))

def mynet(batch,steps,loss_type,dep=False,descr=False,part='gen'):

    conv_lr = [dict(lr_mult=1,decay_mult=1),dict(lr_mult=2,decay_mult=1)]
    bcnv_lr = [dict(lr_mult=1,decay_mult=1)]
    scale_lr = [dict(lr_mult=1,decay_mult=1),dict(lr_mult=1,decay_mult=1)]
    bn_param = dict(eps=0.001,use_global_stats=False)
    #bn_param = dict(eps=0.001,use_global_stats=True)

    fr_lr = [dict(lr_mult=0,decay_mult=0),dict(lr_mult=0,decay_mult=0)]
    fr_clr = [dict(lr_mult=0,decay_mult=0)]
    #fr_bn = dict(eps=0.001,use_global_stats=True)
    fr_bn = dict(eps=0.001,use_global_stats=False)

    if part=='gen':
        gen_conv_lr = conv_lr
        gen_bcnv_lr = bcnv_lr
        gen_scale_lr = scale_lr
        gen_bn_param = bn_param
        dsc_conv_lr = fr_lr
    else:
        gen_conv_lr = fr_lr
        gen_bcnv_lr = fr_clr
        gen_scale_lr = fr_lr
        gen_bn_param = fr_bn
        dsc_conv_lr = conv_lr

    n=caffe.NetSpec()

    sp=dict(bias_term=True,filler=dict(value=1.0))

    if dep:
        n.source = L.Input(input_param=dict(shape=[dict(dim=[1,1,32,32])]))
    else:
        if descr:
            if part == 'gen':
                bs = batch
            else:
                bs = batch / 2
        else:
            bs = batch
        n.data = L.Data(data_param=dict(source="db",batch_size=bs,backend=P.Data.LMDB))
        
        n.expected, n.source = L.Slice(n.data,slice_param=dict(axis=1,slice_point=1),ntop=2)
        if descr:
            if part!='gen':
                #n.data_ref = L.Split(n.expected)
                n.data_ref = L.Data(data_param=dict(source="db_ref",batch_size=batch/2,backend=P.Data.LMDB))
                n.label_0 = L.DummyData(shape=[dict(dim=[batch/2])],data_filler=dict(value=0.0))
                n.label_1 = L.DummyData(shape=[dict(dim=[batch/2])],data_filler=dict(value=1.0))
                n.label = L.Concat(n.label_0,n.label_1,concat_param=dict(axis=0))
            else:
                n.label = L.DummyData(shape=[dict(dim=[batch])],data_filler=dict(value=1.0))

    n.conv1 = L.Convolution(n.source,convolution_param=conv_param_nb(3,16),param=gen_bcnv_lr)
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
        
        n[cid1] = L.Convolution(n[inp],convolution_param=conv_param_nb(3,16),param=gen_bcnv_lr)
        n[bid1] = L.BatchNorm(n[cid1],batch_norm_param=gen_bn_param)
        n[bid1] = L.Scale(n[bid1],scale_param=sp,param=gen_scale_lr)
        n[bid1] = L.ReLU(n[bid1])

        n[cid2] = L.Convolution(n[bid1],convolution_param=conv_param_nb(3,16),param=gen_bcnv_lr)
        n[bid2] = L.BatchNorm(n[cid2],batch_norm_param=gen_bn_param)
        n[bid2] = L.Scale(n[bid2],scale_param=sp,param=gen_scale_lr)
        n[bid2] = L.ReLU(n[bid2])

        n[eid]  = L.Eltwise(n[bid2],n[inp])
        inp = eid

    outname="topconv"
    n[outname]=L.Convolution(n[inp],convolution_param=conv_param(3,1),param=gen_conv_lr)
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
            n.desc_inp = L.Concat(n.generated,n.data_ref,concat_param=dict(axis=0))
            cinp = "desc_inp"
        else:
            cinp = "generated"
        n.d_conv1 = L.Convolution(n[cinp],convolution_param=conv_param(5,32),param=dsc_conv_lr) 
        n.d_pool1 = L.Pooling(n.d_conv1,pooling_param=dict(kernel_size=3,stride=2,pool=P.Pooling.MAX)) 
        n.d_pool1 = L.ReLU(n.d_pool1)
        
        n.d_conv2 = L.Convolution(n.d_pool1,convolution_param=conv_param(5,32),param=dsc_conv_lr) 
        n.d_pool2 = L.Pooling(n.d_conv2,pooling_param=dict(kernel_size=3,stride=2,pool=P.Pooling.MAX)) 
        n.d_pool2 = L.ReLU(n.d_pool2)
        
        #n.d_conv3 = L.Convolution(n.d_pool2,convolution_param=conv_param(5,64),param=dsc_conv_lr) 
        #n.d_pool3 = L.Pooling(n.d_conv3,pooling_param=dict(kernel_size=3,stride=2,pool=P.Pooling.MAX)) 
        #n.d_pool3 = L.ReLU(n.d_pool3)
        
        n.d_conv4 = L.Convolution(n.d_pool2,convolution_param=conv_param(3,64),param=dsc_conv_lr) 
        n.d_pool4 = L.Pooling(n.d_conv4,pooling_param=dict(kernel_size=3,stride=2,pool=P.Pooling.MAX)) 
        n.d_pool4 = L.ReLU(n.d_pool4)

        #n.d_ip1 = L.InnerProduct(n.d_pool4,param=dsc_conv_lr,inner_product_param=ip_param(512))
        #n.d_ip1 = L.ReLU(n.d_ip1)
        n.d_ip2 = L.InnerProduct(n.d_pool4,param=dsc_conv_lr,inner_product_param=ip_param(1))

            
        n.sigmoid_loss = L.SigmoidCrossEntropyLoss(n.d_ip2,n.label,name="loss",loss_weight=100)
        n.score = L.Sigmoid(n.d_ip2)
        n.lbl_flat=L.Reshape(n.label,reshape_param=dict(shape=dict(dim=[-1,1])))
        n.diff = L.Eltwise(n.score,n.lbl_flat,eltwise_param=dict(coeff=[1.0/batch,-1.0/batch]))
        n.error = L.Reduction(n.diff,reduction_param=dict(operation=P.Reduction.ASUM))
        #n.output = L.Split(n[cinp])
        #n.output_labels = L.Split(n.score)
        #n.inputs = n.source 

    return n

def print_to_file(name,n):
    s=str(n.to_proto())
    f=open(name,'w')
    f.write(s)
    f.close()


batch=60
res_steps=2
loss='euc'
#loss='cross_entropy'

print_to_file('gan_dep.prototxt',       mynet(batch,res_steps,loss,True))
print_to_file('gan_pretrain.prototxt',  mynet(batch,res_steps,loss,False))
print_to_file('gan_train_gen.prototxt', mynet(batch,res_steps,loss,False,True,'gen'))
print_to_file('gan_train_dsc.prototxt', mynet(batch,res_steps,loss,False,True,'dsc'))

