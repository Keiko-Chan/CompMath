import numpy as np
import matplotlib.pyplot as plt

def get_node(i, j, h):
	N = int(1 / h + 1)
	return N * i + j

def get_x_y(node, h):
	N = int(1 / h + 1)
	i = node //  N
	j = node % N
	return i, j

def function_elements(i, j, h):
	res = np.zeros(5)
	N = int(1 / h)
	
	if(i == 0 or j == 0 or i == N or j == N):
		return res
	
	res[0] = get_node(i + 1, j, h)
	res[4] = get_node(i, j, h)
	res[2] = get_node(i - 1, j, h)
	res[3] = get_node(i, j + 1, h)
	res[1] = get_node(i, j - 1, h)
	return res

def get_matrix_line(node, h):
	N = int(1 / h + 1)
	res = np.zeros((N - 2)*(N - 2))
	
	i, j = get_x_y(node, h)
	elements = function_elements(i, j, h)
	
	for k in range(0, 5):
		koef = -1
		if(k == 4):
			koef = 4
		
		i1, j1 = get_x_y(elements[k], h)
		if(i1 != 0 and j1 != 0 and i1 != (N - 1) and j1 != (N-1)):
			a = elements[k]
			a = int((a // N - 1)*(N - 2) + a % N - 1)
			res[a] = koef / h / h
			
	return res

def get_matrix(nodes, h):
	res = get_matrix_line(nodes[0], h)
	for k in range(1, len(nodes)):
		matrix_l = get_matrix_line(nodes[k], h)
		res = np.vstack((res, matrix_l))
		
	return res
		
def draw_spy_matrix(matrix):
	fig, axs = plt.subplots(1, 1)
	axs.spy(matrix, markersize=3)
	plt.show()

def my_name_letter():
	letter = np.array([(1, 1, 0, 0, 1, 1, 0), (1, 1, 0, 1, 1, 0, 0), (1, 1, 1, 1, 0, 0, 0), (1, 1, 1, 0, 0, 0, 0), (1, 1, 1, 1, 0, 0, 0), (1, 1, 0, 1, 1, 0, 0), (1, 1, 0, 0, 1, 1, 0)])
	draw_spy_matrix(letter)	
	f = np.array([1, 1, 0, 0, 1, 1, 0, 1, 1, 0, 1, 1, 0, 0, 1, 1, 1, 1, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 1, 1, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0])
	return f

def get_node_array(h):
	N = int(1 / h) + 1
	array_size = (N - 2) * (N - 2)
	nodes = np.zeros(array_size)
	
	l = 0
	for k in range(0, N * N):
		i, j = get_x_y(k, h)
		if(i != 0 and j != 0 and i != (N - 1) and j != (N - 1)):
			nodes[l] = k
			l = l + 1
	return nodes
	
def do_matrix_form_vector(vector):
	N = int(len(vector) ** 0.5)
	matrix = np.zeros((N, N))
	
	for i in range(0, N):
		matrix[i] = vector[i * N : i * N + N]
	return matrix

def first_task():
	nodes = np.array([6, 7, 8, 11, 12, 13, 16, 17, 18])
	f = np.array([1, 1, 1, 1, 1, 1, 1, 1, 1])
	
	matrix = get_matrix( nodes, 1/4)
	print(matrix)
	draw_spy_matrix(matrix)
	
	res = np.linalg.solve(matrix, f) 
	print(res)
	
def second_task():
	h = 1/8
	nodes2 = get_node_array(h)

	matrix = get_matrix( nodes2, h)
	print(matrix)
	letter = my_name_letter()
		
	res = np.linalg.solve(matrix, letter) 
	print(res)
	
	res_matrix = do_matrix_form_vector(res)
	print(res_matrix)
	
	fig, axs = plt.subplots(1, 1)
	axs.imshow(res_matrix)
	plt.show()

def main():
	
	first_task()
	second_task()
	
	
	
if __name__ == "__main__":
	main()
		
