from scipy.integrate import odeint
import numpy as np
from numpy.linalg import norm
import matplotlib.pyplot as plt
from scipy.integrate import RK45

####0-0---0-0---0-0---0-0---0-0---0-0---0-0---0-0---0-0---0-0---0-0---0-0
def u_diff(u, v, A, B):
	res = A + u * u * v - (B + 1) * v
	return res
####0-0---0-0---0-0---0-0---0-0---0-0---0-0---0-0---0-0---0-0---0-0---0-0
def v_diff(u, v, B):
	res = B * u - u * u * v
	return res
####0-0---0-0---0-0---0-0---0-0---0-0---0-0---0-0---0-0---0-0---0-0---0-0	
def func(y, t, A, B):
		u_d = u_diff(y[0], y[1], A, B)
		v_d = v_diff(y[0], y[1], B)
		
		return (u_d, v_d)
####0-0---0-0---0-0---0-0---0-0---0-0---0-0---0-0---0-0---0-0---0-0---0-0
def runge_kut_1(A, B, y0, t0, e, h):	
	y_cur = y0[len(y0) - 1]
	
	k1 = func(y_cur, t0, A, B)
	
	y_new = (y_cur[0] + k1[0] * h,  y_cur[1] + k1[1] * h)
	y_res = np.append(y0, [y_new], axis = 0)
	
	#print(norm(y0 - y_new, 2))
	e = e - 1
	
	if e > 0:
		return runge_kut_1(A, B, y_res, t0 + h, e, h)
	
	else:
		return y_res
####0-0---0-0---0-0---0-0---0-0---0-0---0-0---0-0---0-0---0-0---0-0---0-0		
def plot_runge_kut_t(y_res, t, name):	
	
	y_res = y_res.transpose()
	
	plt.figure(name) 	
	plt.plot(t, y_res[0], '-', t, y_res[1], '--')
	plt.ylabel("- u, -- v")
	plt.xlabel("t")
	plt.grid()
####0-0---0-0---0-0---0-0---0-0---0-0---0-0---0-0---0-0---0-0---0-0---0-0
def runge_kut_4(A, B, y0, t0, h, n):

	time = t0 + h * n
	
	t = np.array(t0)
	y_res = y0
	
	def fun(t, y):
		return func(y, t, A, B)
	
	sol = RK45(fun, t0, y0[0], time, h)

	while(sol.status != "finished" and sol.status != "failed"):
		sol.step()
		y_res = np.append(y_res, [sol.y], axis = 0)
		t = np.append(t, sol.t)	
	
	return y_res, t	
####0-0---0-0---0-0---0-0---0-0---0-0---0-0---0-0---0-0---0-0---0-0---0-0
def razn_resh_plot(A, B, y0, t0):
	y_res_1 = runge_kut_1(A, B, y0, t0, 900, 0.01)
	y_res_4, t_4 = runge_kut_4(A, B, y0, t0, 0.01, 900 - 1)
	
	y_res_1 = y_res_1.transpose()
	y_res_4 = y_res_4.transpose()
	
	u = y_res_1[0] - y_res_4[0]
	v = y_res_1[1] - y_res_4[1]
	
	y_res = np.array([u, v])
	
	plot_runge_kut_t(y_res.transpose(), t_4, "ЯРМК(1) - ЯРМНК(4)")
####0-0---0-0---0-0---0-0---0-0---0-0---0-0---0-0---0-0---0-0---0-0---0-0
def faz_traect_plot(A, B, y0, t0, n, h):
	y_res_1 = runge_kut_1(A, B, y0, t0, n, h)
	y_res_4, t = runge_kut_4(A, B, y0, t0, h, n - 1)

	
	y_res_4 = y_res_4.transpose()
	y_res_1 = y_res_1.transpose()
	
	plt.figure("RK1 faz") 	
	plt.plot(y_res_1[0], y_res_1[1], '-')
	plt.ylabel("u")
	plt.xlabel("v")
	plt.grid()
	
	plt.figure("RK4 faz") 	
	plt.plot(y_res_4[0], y_res_4[1], '-')
	plt.ylabel("u")
	plt.xlabel("v")
	plt.grid()

####0-0---0-0---0-0---0-0---0-0---0-0---0-0---0-0---0-0---0-0---0-0---0-0	
def main():
	y0 = np.array([[1, 1]])
	t0 = 0
	A = 1
	B = 2
	n = 900
	h = 0.01
	#Рунге Кутт 1 порядка"
	name_1 = "Рунге Кутт 1 порядка"	
	y_res_1 = runge_kut_1(A, B, y0, t0, n, h)
	t_1 = np.arange(0, h * len(y_res_1), h)
	plot_runge_kut_t(y_res_1, t_1, name_1)
	#Рунге Кутт 4 порядка
	name_4 = "Рунге Кутт 4 порядка"
	y_res_4, t_4 = runge_kut_4(A, B, y0, t0, h, n)
	plot_runge_kut_t(y_res_4, t_4, name_4)
	
	#разница решений
	razn_resh_plot(A, B, y0, t0)
	
	#фазовые траектории
	faz_traect_plot(A, B, y0, t0, n, h)
	
	#print(len(y_res_4), len(y_res_1))
	#print(t_4)
	
	plt.show()
####0-0---0-0---0-0---0-0---0-0---0-0---0-0---0-0---0-0---0-0---0-0---0-0
if __name__ == "__main__":
	main()	


