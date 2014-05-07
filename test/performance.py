from time import time
import numpy as np
import cvxpy as cp
def func():
  N= 200
  start = time()
  x=cp.Variable(1)
  y=cp.Variable(N, N)
  X=cp.Variable(N, N)
  c = np.ones((N, 1))
  obj = cp.Minimize(c.transpose() * (y+X) * c)
  constr = [x>=3, y <= X, 2*y>=0]
  prob = cp.Problem(obj, constr)
  prob.solve()
  print time() - start
