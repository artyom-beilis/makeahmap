import paths
import sys

for p in paths.external:
    sys.path.append(p)

import caffe

caffe.set_mode_gpu()
caffe.set_device(0)

s1 = caffe.get_solver('solver.prototxt')
s2 = caffe.get_solver('solver.prototxt')
s1.share_weights(s2.net)

for step in range(3):
    for k in range(50):
        s1.step(10)
        print "Step 1",k
    for k in range(50):
        s2.step(10)
        print "Step 2",k    


