from time import time
import numpy as np
import cvxpy as cp
def func():
  N= 200
  start = time()
  x=cp.Variable(1)
  y=cp.Variable(N, N)
  c = np.ones((N, 1))
  obj = cp.Minimize(c.transpose() * (y+x) * c)
  constr = [x>=3, y <= x, 2*y>=0];
  prob = cp.Problem(obj, constr);
  print time() - start
  prob.solve();
  print time() - start
