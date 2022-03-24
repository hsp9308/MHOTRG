import time
import numpy as np
from scipy.optimize import minimize_scalar
import subprocess
import sys


def func(x):
	params = x
	s = subprocess.check_output(['./TRG_it'] + ['%d' % int(sys.argv[1])] + [str(params)])
	print(s)
	s1 = complex(s)
	sr = np.real(s1)
	si = np.imag(s1)
	print(sr)
	print(params)
	return (sr)



start_time=time.time()
res = minimize_scalar(func, bounds=(0.0002250,0.0002341),method='bounded',options={'disp':True,'xatol':1e-10})

print(res.x)
print("--------  %s seconds elapsed  --------" % (time.time() - start_time))
