from keras.models import Sequential
from keras.layers import Dense
from keras.layers import Dropout
from keras import utils
from keras import optimizers
import scipy.io as sio
import numpy as np

mat_contents = sio.loadmat('InOutRE_ValidR4E_2D_q-1224&-2501&-1211_dq5&15&15_tau200250250_2392.mat')
xall = mat_contents['InR']
xall[:, 3:4] /= 2.5
xall[:, 4:6] /= 7.5
xall[:, 6:9] /= 100.0
eall = mat_contents['OutE']
rall = mat_contents['OutR']

# exclude tail too-large values
xall = xall[(eall[:,0]<1700),:]
rall = rall[(eall[:,0]<1700),:]
eall = eall[(eall[:,0]<1700),:]

yall = eall / 10.0

xyall = np.concatenate((xall,yall), axis=1)
np.random.shuffle(xyall)
xall = xyall[:,:9]
yall = xyall[:,9:10]

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
model.add(Dense(256, input_dim = 9, activation='elu'))
# model.add(Dropout(0.02))
model.add(Dense(180, activation='elu'))
# model.add(Dropout(0.02))
model.add(Dense(64, activation='elu'))
# model.add(Dropout(0.02))
model.add(Dense(1))

n_epoch = 700
model.compile(loss='mean_squared_error', optimizer='adam')
model.fit(x_train, y_train, epochs=n_epoch, validation_split = 0.05, batch_size=512)
score = model.evaluate(x_test, y_test, batch_size=512)
print(score)

model.save('InOutE_tilde_2D_q-1224&-2501&-1211_dq5&15&15_tau200250250_2392_512360180elu_mse_ediv10_r4e.h5')

print(y_test)
y_pred = model.predict(x_test, batch_size=512)
print(y_pred)
print(np.linalg.norm(y_pred-y_test)/np.linalg.norm(y_test))
