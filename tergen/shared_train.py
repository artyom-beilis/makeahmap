import paths
import sys

for p in paths.external:
    sys.path.append(p)

import caffe
import numpy as np
import signal

stop=False

def handler(sig,frame):
    global stop
    stop=True

signal.signal(signal.SIGINT,handler)

caffe.set_mode_gpu()
caffe.set_device(0)

slv_dsc = caffe.get_solver('solver_dsc.prototxt')
slv_gen = caffe.get_solver('solver_gen.prototxt')

slv_dsc.net.copy_from('snapshots/pretrain_iter_1000.caffemodel')
slv_dsc.share_weights(slv_gen.net)

def run_steps(slv,steps):
    hist = np.zeros([2,steps])
    for k in range(steps):
        slv.step(1)
        hist[0,k] = 1 - slv.net.blobs['error'].data
        hist[1,k] = slv.net.blobs['l2_loss'].data
    #print np.round(hist[1,:]*100) / 100
    accuracy = np.mean(hist[0,:])
    l2_loss = np.mean(hist[1,:])
    return (accuracy,hist[0,:],l2_loss)


gen_steps = 0
dsc_steps = 0
th = 0.7

while gen_steps < 200000 and not stop:
    while not stop:
        dsc_step=1
        accuracy,hist,l2 = run_steps(slv_dsc,dsc_step)
        print 'Dsc accuracy[%d] = %f l2=%f %s' % (dsc_steps,accuracy,l2,str(np.round(hist*100)/100))
        dsc_steps += dsc_step
        if accuracy > th:
            break
    while not stop:
        gen_step=1
        accuracy,hist,l2 = run_steps(slv_gen,gen_step)
        gen_steps+=gen_step
        print 'Gen accuracy[%d] = %f l2=%f %s' % (gen_steps,accuracy,l2,str(np.round(hist*100)/100))
        if accuracy > th: # and l2 < 30:
            #slv_gen.snapshot()
            break
slv_gen.snapshot() 
