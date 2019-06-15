import scipy.io as sio
import numpy as np
from keras.models import load_model

model = load_model('path_to_your_trained_python_NN_model_for_R_2D.h5')

samSize = 12000000
q1_init = np.random.uniform(-1.2, 2.4, (samSize,1))
q4_init = np.random.uniform(-2.5, 0.1, (samSize,1))
q5_init = np.random.uniform(-1.2, 1.1, (samSize,1))
q_init = np.concatenate((q1_init, q4_init, q5_init), axis=1)

dq1_init = np.random.uniform(-5.0, 5.0, (samSize,1))
dq4_init = np.random.uniform(-15.0, 15.0, (samSize,1))
dq5_init = np.random.uniform(-15.0, 15.0, (samSize,1))
dq_init = np.concatenate((dq1_init, dq4_init, dq5_init), axis=1)

tau1_init = np.random.uniform(-200, 200, (samSize,1))
tau45_init = np.random.uniform(-250, 250, (samSize,2))
tau_init = np.concatenate((tau1_init, tau45_init), axis=1)

x_init = np.concatenate((q_init, dq_init, tau_init), axis=1)

x_nn = np.copy(x_init)
x_nn[:, 3:4] /= 2.5
x_nn[:, 4:6] /= 7.5
x_nn[:, 6:9] /= 100.0

residues = model.predict(x_nn, batch_size=1000)
norm_res = np.linalg.norm(residues, axis=1)

x_sampling = x_init[norm_res < 5,:]
print(x_sampling.shape)

with open('R4E_2D.mat','wb') as f:
	sio.savemat(f, {'valid_qdqtau': x_sampling})