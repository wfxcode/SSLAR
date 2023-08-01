
import numpy as np

from ssnmf import SSNMF
import os
import torch
import copy
os.environ["KMP_DUPLICATE_LIB_OK"]="TRUE"
# os.environ["CUDA_DEVICE_ORDER"] = "PCI_BUS_ID"
os.environ["CUDA_VISIBLE_DEVICES"] = "-1"

# # Numpy不遗漏数据的有监督
def ssNMF(X, k, Y, A=None, S=None, N=100, lam=100, mult=True):
    if A is None or S is None:
      # lam * X矩阵的2范式
        nmf_mod = SSNMF(X=X, k=k, Y=Y, lam=lam * np.linalg.norm(X), modelNum=3)
    else:
        nmf_mod = SSNMF(X=np.array(X), k=int(k), Y=np.array(Y), lam=lam * np.linalg.norm(X), A=np.array(A), S=np.array(S), modelNum=3)
    if not mult:
        return nmf_mod
    else:
        nmf_mod.mult(numiters=N, saveerrs=True)
    return nmf_mod
  
def getTensor(m, device):
  	mt = torch.from_numpy(copy.deepcopy(m))
  	mt = mt.type(torch.FloatTensor)
  	mt = mt.to(device)
  	return mt

def torch_ssNMF(X, k, Y, A, S, N=100, lam=100, mult=True):
  	devise = torch.device("cuda") if torch.cuda.is_available() else torch.device("cpu")
  	# print(1)
  	Xt = getTensor(X, devise)
  	Yt = getTensor(Y, devise)
  	At = getTensor(A, devise)
  	Bt = torch.eye(k)
  	St = getTensor(S, devise)
  	# print(torch.cuda.is_available())
  	# print(lam * torch.norm(Xt))
  	# nmf_mod = SSNMF(Xt, int(k), Y = Yt, lam=lam * torch.norm(Xt), A=At, B=Bt , S=St, modelNum=5)
  	nmf_mod = SSNMF(Xt, int(k), Y = Yt, lam=lam * torch.norm(Xt), A=At, S=St, modelNum=6)
  	arr = nmf_mod.B.numpy()
  	# print(arr.shape)
  	arr = nmf_mod.W.numpy()
  	# print((arr == 1).all())
  	arr = nmf_mod.L.numpy()
  	# print((arr == 1).all())
  	nmf_mod.mult(numiters=N, saveerrs=True)
  	nmf_mod.A = nmf_mod.A.numpy()
  	nmf_mod.S = nmf_mod.S.numpy()
  	return nmf_mod
