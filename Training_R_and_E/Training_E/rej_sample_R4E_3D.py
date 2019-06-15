import scipy.io as sio
import numpy as np
from keras.models import load_model

model = load_model('path_to_your_trained_python_NN_model_for_R_3D.h5')

samSize = 100000000
q1_init = np.random.uniform(-0.9, 1.5, (samSize,1))
q23_init = np.random.uniform(-0.6, 0.6, (samSize,2))
q4_init = np.random.uniform(-2.1, 0.1, (samSize,1))
q5_init = np.random.uniform(-0.9, 0.9, (samSize,1))
q_init = np.concatenate((q1_init, q23_init, q4_init, q5_init), axis=1)

dq123_init = np.random.uniform(-5.0, 5.0, (samSize,3))
dq4_init = np.random.uniform(-10.0, 10.0, (samSize,1))
dq5_init = np.random.uniform(-15.0, 15.0, (samSize,1))
dq_init = np.concatenate((dq123_init, dq4_init, dq5_init), axis=1)

tau1_init = np.random.uniform(-250, 250, (samSize,1))
tau23_init = np.random.uniform(-80, 80, (samSize,2))
tau45_init = np.random.uniform(-250, 250, (samSize,2))
tau_init = np.concatenate((tau1_init, tau23_init, tau45_init), axis=1)

x_init = np.concatenate((q_init, dq_init, tau_init), axis=1)

x_nn = np.copy(x_init)
x_nn[:, 5:9] /= 2.5
x_nn[:, 9:10] /= 7.5
x_nn[:, 10:15] /= 100.0

residues = model.predict(x_nn, batch_size=1000)
norm_res = np.linalg.norm(residues, axis=1)

x_sampling = x_init[norm_res < 5,:]
print(x_sampling.shape)

with open('R4E_3D.mat','wb') as f:
	sio.savemat(f, {'valid_qdqtau': x_sampling})