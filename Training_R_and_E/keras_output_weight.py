import scipy.io as sio
import numpy as np
from keras.models import load_model

model = load_model('your_trained_NN_model.h5')

with open('your_MATLAB_NN_model_assuming_elu_activation.mat','wb') as f:
    j = 0
    for i in range(0,len(model.layers)):
        if len(model.layers[i].get_weights()) > 0:              # workaround for dropout layers
            Wname = 'W' + str(j)
            print(Wname)
            Bname = 'B' + str(j)
            Wmat = model.layers[i].get_weights()[0]
            Wmat = np.transpose(Wmat)
            print(Wmat.shape)
            Bmat = model.layers[i].get_weights()[1]
            Bmat = Bmat.reshape((-1,1))
            print(Bmat.shape)
            sio.savemat(f, {Wname: Wmat})
            sio.savemat(f, {Bname: Bmat})
            j += 1

f.close()