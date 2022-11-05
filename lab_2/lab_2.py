import math
import matplotlib.pyplot as plt
import numpy as np

def get_jacobi_line(h, node, u):
	N = int(1 / h) + 1
	line = np.zeros(N - 2)							#without u0 and u(N - 1)

	if(node > 1):
		line[node - 2] = -1 / h / h
	
	line[node - 1] = 2 / h / h + 2 * u[node]
	
	if(node < N - 2):
		line[node] = -1 / h / h
	
	return line

def get_jacobi_matrix(h, u):
	N = int(1 / h) + 1
	j_mat = get_jacobi_line(h, 1, u)
	
	for i in range(2, N - 1):
		new_line = get_jacobi_line(h, i, u)
		j_mat = np.vstack((j_mat, new_line))
	
	return j_mat

def f_function(xi):
	if(xi > 1 or xi < 0):
		print("not needed range")
		return 0
	
	if(xi <= 0.3 and xi >= 0.2):
		return math.sin(math.pi * (xi - 0.2) / 0.1) 
	
	if(xi <= 0.7 and xi >= 0.6):
		return 0.5 * math.sin(math.pi * (xi - 0.6) / 0.1) 
	
	return 0

def F_function(x, u, h, node):
	N = int(1 / h) + 1
	#node from 0 to N - 1 	U(0) = 0 = U(N - 1)
	if(node == 0):
		return 0
	if(node == N - 1):
		return 0
		
	res = u[node] * u[node] + 2 * u[node] / h / h -  u[node - 1] / h / h - u[node + 1] / h / h - f_function(x[node])
	return res

def solve_diff_recur(h, un, x, e, u0, ite):
	N = int(1 / h) + 1
	F_un = np.zeros(N)
	un1 = np.zeros(N)
	
	for i in range(0, N):
		F_un[i] = F_function(x, un, h, i)
	
	j_mat = get_jacobi_matrix(h, un)
	inv_j_mat = np.linalg.inv(j_mat)
	
	delt = inv_j_mat.dot(F_un[1:N-1])
	
	for i in range(1, N - 1):
		un1[i] = un[i] - delt[i - 1]
	
	#print(un1, un, np.linalg.norm(un1 - u0), np.linalg.norm(un - u0))
	
	if(e < 0):
		return "error", -1
	
	if(isinstance(e, float)):						#calculate with occurasy
		ite = ite + 1
		
		if(abs(np.linalg.norm(un1 - u0) - np.linalg.norm(un - u0)) > e):
			un1, ite = solve_diff_recur(h, un1, x, e, u0, ite)
	
		return un1, ite
	
	if(isinstance(e, int)):							#calculate if you know number of iterations
		ite = un1 - un
		
		if(e - 1 > 0):
			un1, ite = solve_diff_recur(h, un1, x, e - 1, u0, ite)
	
		return un1, ite

def solve_diff_start(h, e, x):
	N = int(1 / h) + 1
	u = np.zeros(N)								#U0 - started
	
	u[N-1] = 0 
	u[0] = 0
	
	res_U = solve_diff_recur(h, u, x, e, u, 0)
	
	return res_U

def calc_converge(x, du, h, u):
	N = int(1 / h) + 1
	
	F_u = np.zeros(N)
	j_mat = get_jacobi_matrix(h, u)
	
	for i in range(0, N):
		F_u[i] = F_function(x, u, h, i)
		
	Ju = j_mat.dot(u[1:N-1])
	
	y = math.log(np.linalg.norm(Ju - F_u[1:N-1])/ np.linalg.norm(F_u), 10)
	
	return y
	

def plot_result(x, u, i, name1, name2):
	plt.figure(i)  
	plt.plot(x, u)
	plt.ylabel(name1)
	plt.xlabel(name2)
	plt.grid()

def plot_f_func(i, h):
	x = np.arange(0, 1 + h, h)
	y = np.zeros(len(x))
	
	for j in range(0, len(x)):
		y[j] =  f_function(x[j])
	
	plt.figure(i) 
	plt.plot(x, y)
	plt.ylabel("f func")
	plt.xlabel("x")
	plt.grid()

def plot_converge(x, h, N, i, param = 0):							#N - max number of iter	
	n = np.arange(1, N + 1, 1)
	#du = np.zeros((N, int(1/h + 1)))
	y = np.zeros(N)
	
	for j in range(0, N):
		u, du = solve_diff_start(h, j + 1, x)
	
		match param:
			case 0:
				y[j] = calc_converge(x, du, h, u)
				name = "converge"
		
			case _:
				y[j] = np.linalg.norm(du) 
				name = "||du||"
	
	plt.figure(i) 
	plt.plot(n, y)
	plt.ylabel(name)
	plt.xlabel("n - iterations")
	plt.grid()	
	
	

def main():
	h = 1 / 100
	
	x = np.arange(0, 1 + h, h)
	u, it = solve_diff_start(h, 0.00001, x)
	
	print("number of iterations = ", it)
	print("result vector = ", u)
	
	plot_result(x, u, 1, "U", "x")
	plot_f_func(2, h)
	
	plot_converge(x, h, 10, 3, 0)
	
	plt.show()
	
	
if __name__ == "__main__":
	main()
