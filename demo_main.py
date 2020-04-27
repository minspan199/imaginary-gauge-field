import cv2
from phase_retrieval_GS import *
import matplotlib.pyplot as plt
import numpy as np

imgNear = cv2.imread("D:/Google Drive/Non-Hermitian Gauge Field/non-hermitian-gauge-field/Phase_retrieval_single_constraint-master/A_Near.png", cv2.IMREAD_GRAYSCALE)
imgNear = imgNear.astype(float)
imgNear = np.asarray(imgNear, float)


imgFar = cv2.imread('A_Far.png', cv2.IMREAD_GRAYSCALE)
imgFar = imgFar.astype(float)
imgFar = np.asarray(imgFar, float)



max_iters = 1000
phase_mask = Ger_Sax_algo(imgNear, imgFar, max_iters)
plt.figure(1)
plt.subplot(131)
plt.imshow(imgFar)
plt.title('Desired image')
plt.subplot(132)
plt.imshow(phase_mask)
plt.title('Phase mask')
plt.subplot(133)
recovery = np.fft.ifft2(np.exp(phase_mask * 1j))
plt.imshow(np.absolute(recovery)**2)
plt.title('Recovered image')
plt.show()


