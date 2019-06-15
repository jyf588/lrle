from keras.models import Sequential
from keras.layers import Dense
from keras.layers import Dropout
from keras import utils
from keras import optimizers
import scipy.io as sio
import numpy as np


mat_contents = sio.loadmat('InOutR_2_4M_3D_q-0915&-0606&-0606&-2101&-0909_dq555&10&15_tau250&80&80&250&250_2392.mat')
xall = mat_contents['InR']

# whitening input a bit
xall[:, 5:9] /= 2.5
xall[:, 9:10] /= 7.5
xall[:, 10:15] /= 100.0

yall = mat_contents['OutR']

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
model.add(Dense(256, input_dim = 15, activation='elu'))
# model.add(Dropout(0.01))
model.add(Dense(180, activation='elu'))
# model.add(Dropout(0.01))
model.add(Dense(64, activation='elu'))
# model.add(Dropout(0.01))
model.add(Dense(5))

n_epoch = 700
model.compile(loss='mean_squared_error', optimizer='adam')
model.fit(x_train, y_train, epochs=n_epoch, validation_split = 0.05, batch_size=512)
score = model.evaluate(x_test, y_test, batch_size=512)
print(score)


model.save('InOutR_2_4M_3D_q-0915&-0606&-0606&-2101&-0909_dq555&10&15_tau250&80&80&250&250_2392_elu.h5')

y_pred = model.predict(x_test, batch_size=512)
print(y_pred)
print(np.linalg.norm(y_pred[:,0]-y_test[:,0])/np.linalg.norm(y_test[:,0]))
print(np.linalg.norm(y_pred[:,1]-y_test[:,1])/np.linalg.norm(y_test[:,1]))
print(np.linalg.norm(y_pred[:,2]-y_test[:,2])/np.linalg.norm(y_test[:,2]))
print(np.linalg.norm(y_pred[:,3]-y_test[:,3])/np.linalg.norm(y_test[:,3]))
print(np.linalg.norm(y_pred[:,4]-y_test[:,4])/np.linalg.norm(y_test[:,4]))
print(np.linalg.norm(y_pred-y_test)/np.linalg.norm(y_test))
