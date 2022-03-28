import numpy as np
theta = np.array([[0 ,1 ,1 ,0 ,1 ,0 ,0 ,0],[1, 0, 0, 1, 0, 1, 0, 0],[1, 0, 0, 1, 0, 0, 1, 0],[0, 1, 1, 0, 0, 0, 0, 1],[1, 0, 0, 0, 0, 1, 1, 0],
                 [0, 1, 0, 0, 1, 0, 0, 1],[0, 0, 1, 0, 1, 0, 0, 1],[0, 0, 0, 1, 0, 1, 1, 0]])

def sum_ising_prob(state_vec,theta):
    return 1/2*(np.matmul(state_vec,np.matmul(theta,state_vec)))

def current_state(state_vec,theta):
    mu = abs(sum(state_vec))
    k = sum_ising_prob(state_vec,theta)
    if mu == 8:
        return 0
    elif mu == 6:
        return 1
    elif mu == 4:
        if k == 4:
            return 2
        else:
            distinc_vec = np.matmul(theta,state_vec)
            if -1 in distinc_vec:
                if np.count_nonzero(distinc_vec+3) != 6:
                    return 3
                return 4
            else:
                return 4
    elif mu == 2:
        if k == 2:
            return 5
        elif k == -2:
            return 6
        else:
            return 7
    else:
        if k == 0:
            distinc_vec = np.matmul(theta,state_vec)
            if 3 in distinc_vec:
                return 8
            else:
                return 9
        elif k == 4:
            return 10
        elif k == -4:
            distinc_vec = np.matmul(theta,state_vec)
            if 3 in distinc_vec:
                return 11
            else:
                return 12
        else:
            return 13
            
def prob_brute_force(theta,J):
    states = np.zeros(14)
    for i in range(256):
        vec = np.array([int(x) for x in list('{:08b}'.format(i))])
        aux_vec = 1 - vec;
        vec = vec*-1;
        state_vec = vec + aux_vec;
        #states[current_state(state_vec,theta)]+=1 #guardes quants cops entres a cada estat
        k = sum_ising_prob(state_vec,theta)
        states[current_state(state_vec,theta)]+= np.exp(k*J)
    states/=sum(states)
    return states

prob_brute_force(theta,1)

def gibbs_sampling():
    T = np.array([[0,1,0,0,0,0,0,0,0,0,0,0,0,0],
                  [1,0,1,1,1,0,0,0,0,0,0,0,0,0],
                  [0,1,0,0,0,1,1,0,0,0,0,0,0,0],
                  [0,1,0,0,0,1,1,1,0,0,0,0,0,0],
                  [0,1,0,0,0,0,1,0,0,0,0,0,0,0],
                  [0,0,1,1,0,0,0,0,1,1,1,1,0,0],
                  [0,0,1,1,1,0,0,0,0,1,0,1,1,0],
                  [0,0,0,1,0,0,0,0,1,0,0,1,0,1],
                  [0,0,0,0,0,1,0,1,0,0,0,0,0,0],
                  [0,0,0,0,0,1,1,0,0,0,0,0,0,0],
                  [0,0,0,0,0,1,0,0,0,0,0,0,0,0],
                  [0,0,0,0,0,1,1,0,0,0,0,0,0,0],
                  [0,0,0,0,0,0,1,0,0,0,0,0,0,0],
                  [0,0,0,0,0,0,0,1,0,0,0,0,0,0]])
    
    
    
    
    