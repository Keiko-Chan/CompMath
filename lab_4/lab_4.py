import numpy as np
import matplotlib.pyplot as plt
from numpy.linalg import norm

def get_node(i, j, N):
    M = N - 2
    i = i - 1
    j = j - 1
    
    bord = int((M - 1) / 2)
    
    if i > bord:
        node = M * bord + (i - bord) * (bord) + j
        
    else:
        node = i * M + j
        
    return int(node)

def get_i_j(node, N):
    M = N - 2
    bord = int((M - 1) / 2)
    i0 = node //  M
	
    if i0 > bord:
        j = node % bord
        i = bord + (node - M * bord) // bord
	    
    else:
        i = i0
        j = node % M
	    
    j = node % N
    return i + 1, j + 1

def get_x_y(i, j, h):
    x0 = -1
    y0 = 1
    
    x = x0 + (i + 1) * h
    y = y0 - (j + 1) * h
    
    return x, y

def border(i, j, h):
    x, y = get_x_y(i, j, h)
    
    #on the N1
    if x >= 1 and y <= 1 and y >= 0:
        return y 
    
    #on the N2
    if x <= 0 and y <= -1 and x >= -1:
        return -x
    
    return 0
    
def function_elements(i, j, h, N):
	res = np.zeros((5, 2))
	add = 0
	
	if i + 1 >= N - 1 or (i + 1 >= (N - 1) / 2 and j >= (N - 1) / 2):
	    res[0] = (None, None)
	    add = add + border(i + 1, j, h)
	else:
	    res[0] = (i + 1, j)
	
	if i - 1 <= 0:
	    res[2] = (None, None)
	    add = add + border(i - 1, j, h)
	else:
	    res[2] = (i - 1, j)
	    
	if j - 1 <= 0:
	    res[1] = (None, None)
	    add = add + border(i, j - 1, h)
	else:
	    res[1] = (i, j - 1)
	 
	if j + 1 >= N - 1 or (j + 1 >= (N - 1) / 2 and i >= (N - 1) / 2):
	    res[3] = (None, None)
	    add = add + border(i, j + 1, h)
	else:
	    res[3] = (i, j + 1)
	    
	res[4] = (i, j)	    
	
	return res, add

def get_koef(k, i1, j1, tau, h, cur):
    a = 500
    
    x, y = get_x_y(i1, j1, h)
    
    if y <= 2/3 and y >= 1/3 and x <= -1/3 and x >= -2/3:
        a = 4000
    
    c = a / 2 / h / h * tau
    
    if cur == "current":
        match k:
            case 4:
                return 1 - 4*c
            case _:
                return c        
    
    if cur == "new":
        match k:
            case 4:
                return 1 + 4*c
            case _:
                return -c 

def get_matrix_line(i, j, h, tau, N, cur):
    res = np.zeros((N - 2)*(N - 2) - int((N - 1) / 2 * (N - 1) / 2))

    elements, add = function_elements(i, j, h, N)
	
    for k in range(0, 5):
        i1, j1 = elements[k]
		
        if not np.isnan(i1) and not np.isnan(j1):
            node = get_node(i1, j1, N)
            koef = get_koef(k, i1, j1, tau, h, cur)
            #print(node, i, j, i1, j1)
            res[node] = koef
			
    return res, add * get_koef(5, 0, 0, tau, h, cur)

def get_matrix(h, N, tau, cur):
    res, add = get_matrix_line(1, 1, h, tau, N, cur)    
    add_vec = np.array([add])
	
    for i in range(1, N - 1):
        for j in range(1, N - 1):
            if (i >= (N - 1) / 2 and j >= (N - 1) /2) or (i == 1 and j == 1):
                continue
	        
            else:
                line, add = get_matrix_line(i, j, h, tau, N, cur)

                res = np.vstack((res, line))
                add_vec = np.append(add_vec, [add])
		
    return res, add_vec
		
def draw_spy_matrix(matrix):
	fig, axs = plt.subplots(1, 1)
	axs.spy(matrix, markersize=3)
	plt.show()

def iteration(u_cur, h, N, tau, e, norm_it):
    mat_cur, add_cur = get_matrix(h, N, tau, "current")
    mat_new, add_new =  get_matrix(h, N, tau, "new")
    
    #fig, axs = plt.subplots(1, 1)
    #axs.spy(mat_cur, markersize=3)
    #plt.show()
    
    f = mat_cur.dot(u_cur) + add_cur - add_new
    u_new = np.linalg.solve(mat_new, f) 
    
    norm_cur = norm(u_cur, 1)
    norm_new = norm(u_new, 1)
    
    print(u_cur)
    
    if norm_it is None:
        norm_it = np.array([norm_new])
    
    else:
        norm_it = np.append(norm_it, norm_new)
    
    if abs(norm_new - norm_cur) < e:
        return u_new, norm_it
    
    else:
        return iteration(u_new, h, N, tau, e, norm_it)

def matrix_from_vec(answer, N, h):
    mat = np.zeros((N , N))
    
    for i in range(1, N - 1):
        for j in range(1, N - 1):
            if (i >= (N - 1) / 2 and j >= (N - 1) /2) or (i == 1 and j == 1):
                continue
	        
            else:
                node = get_node(i, j, N)
                mat[i][j] = answer[node]
    
    for i in range(1, int((N - 1) / 2)):
        x,y = get_x_y(i, j, h)
        mat[i][N - 1] = -x
    
    for j in range(1, int((N - 1) / 2)):
        x,y = get_x_y(i, j, h)
        mat[N - 1][j] = y
        
    fig, axs = plt.subplots(1, 1)
    axs.imshow(mat)
    plt.show()

def plot_norm_shod(norm):
    n = np.arange(0, len(norm), 1)
    
    plt.figure(2) 
    plt.plot(n, norm)
    plt.ylabel("norm")
    plt.xlabel("n - iterations")
    plt.grid()	
    plt.show()

def main():	
	h = 1 / 20
	N = int(2 / h + 1)
	e = 0.5
	tau = 0.0001
	
	u_cur = np.zeros((N - 2)*(N - 2) - int((N - 1) / 2 * (N - 1) / 2))
	
	answer, norm = iteration(u_cur, h, N, tau, e, None)
	matrix_from_vec(answer, N, h)
	
	plot_norm_shod(norm)
	
if __name__ == "__main__":
	main()
