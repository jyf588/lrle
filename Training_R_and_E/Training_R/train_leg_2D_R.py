from keras.models import Sequential
from keras.layers import Dense
from keras.layers import Dropout
from keras import utils
from keras import optimizers
import scipy.io as sio
import numpy as np


mat_contents = sio.loadmat('InOutR_300k_q-1224&-2501&-1211_dq5&15&15_tau200250250_2392.mat')
xall = mat_contents['InR']
# whitening input a bit
xall[:, 3:4] /= 2.5
xall[:, 4:6] /= 7.5
xall[:, 6:9] /= 100.0

yall = mat_contents['OutR']

# # Reverse columns (state and u orders), for the jumping example (Traj optimization)
# xall[:,[0, 1, 2, 3, 4, 5, 6, 7, 8]] = xall[:, [2, 1, 0, 5, 4, 3, 8, 7, 6]]
# yall[:,[0, 1, 2]] = yall[:,[2, 1, 0]] 

length = yall.shape[0]

print(xall)
print(yall)

x_train = xall[0:int(length*0.9),:]
x_test = xall[int(length*0.9):,:]

y_train = yall[0:int(length*0.9),:]
y_test = yall[int(length*0.9):,:]

print(x_train.shape)
print(y_train.shape)
print(x_test.shape)
print(y_test.shape)

model = Sequential()
model.add(Dense(128, input_dim = 9, activation='elu'))
# model.add(Dropout(0.05))
model.add(Dense(128, activation='elu'))
# model.add(Dropout(0.05))
model.add(Dense(128, activation='elu'))
# model.add(Dropout(0.05))
model.add(Dense(3))

n_epoch = 600
model.compile(loss='mean_squared_error', optimizer='adam')
model.fit(x_train, y_train, epochs=n_epoch, validation_split = 0.05, batch_size=256)
score = model.evaluate(x_test, y_test, batch_size=256)
print(score)


model.save('InOutR_300k_2392_q-1224&-2501&-1211_dq5&15&15_tau200250250_elu.h5')

print(y_test)
y_pred = model.predict(x_test, batch_size=256)
print(y_pred)
print(np.linalg.norm(y_pred[:,0]-y_test[:,0])/np.linalg.norm(y_test[:,0]))
print(np.linalg.norm(y_pred[:,1]-y_test[:,1])/np.linalg.norm(y_test[:,1]))
print(np.linalg.norm(y_pred[:,2]-y_test[:,2])/np.linalg.norm(y_test[:,2]))
print(np.linalg.norm(y_pred-y_test)/np.linalg.norm(y_test))
