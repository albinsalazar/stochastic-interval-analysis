import numpy as np

k_1 = 20  # Reaction rate of a first reaction.
k_2 = 1 # Reaction rate of a second reaction.
reac_rates = [k_1, k_2] # List of reaction rate constants.
U_M1 = int(input('Enter the first upper bound:'))
U_M2 = int(input('Enter the second upper bound:'))
N_B = int(input('Coordinate value of N_B:'))
N_C = int(input('Coordinate value of N_C:'))
q = np.array([N_B, N_C]) # The states values for molecules B and C, respectively.
upper_bounds = np.array([U_M1, U_M2])  # Upper bound values for the x- and y-axis, respectively.
init_states = [100, N_B, N_C] # Initial values for molecules A, B and C, respectively.
k = 100 # Number of steps.
L = False # toggle between coordinates
print(q, upper_bounds)
def gen_intervals(q, upper_bounds):
    # q := list of coordinate system, upper_bounds := list of upper bounds for each coordinate.
    '''Here, we generate a set of coordinates contained in a (2-D) Box.'''
    coordinate_list = []
    for i in range(q[0], upper_bounds[0]): # Range from value of the first coordinate to its upper bound.
        for j in range(q[1], upper_bounds[1]): # Range from value of the second coordinate to its upper bound.
            coordinate_list.append((i,j)) # Append each pair, a point in a 2-D plane.
    return coordinate_list

def U_M(q1, q2, upper_bounds, k, L):
    # q := list of coordinate system, k := number of permitted steps
    '''The goal of this function is to compute the probability to reach a first upper bound.
    Here, we have a recurrence function with two terms.'''

    # Base case: the boundary conditions for the 2-D recurrence
    # print(q1,q2)
    if L == True: # Reach the upper bound on the x-axis first.
        if q2 >= upper_bounds[1] or k <= 0: # Not upper bound of interest or no more steps.
            return 0
        elif q1 >= upper_bounds[0]: # The target goal has been met.
            return 1
        else: # A state is inside the interval of interest.
             a = (init_states[0] + init_states[1] + 2*init_states[2] - q1 - 2*q2) # Invariant for our reactions for resource competition. 
             if a<=0: # Quantity of chemical species is non-negative.
                 return 0
             else:
                 lambda_1 = reac_rates[0] # k1
                 lambda_2 = reac_rates[1]*(a-1) # k2 * invariant
                 lambda_t = lambda_1 + lambda_2
                 #print(lambda_1,lambda_2,lambda_t)
                 return (lambda_1/lambda_t) * U_M(q1+1, q2, upper_bounds, k-1, L) + (lambda_2/lambda_t) * U_M(q1,q2+1,upper_bounds, k-1, L)
    # Alternate upper bound probability
    else:
        if q1 >= upper_bounds[0] or k <= 0:
            return 0
        elif q2 >= upper_bounds[1]:
            return 1
        else:
             a = (init_states[0] + init_states[1] + 2*init_states[2] - q1 - 2*q2) 
             if a<=0:
                 return 0
             else:
                 lambda_1 = reac_rates[0] # k1
                 lambda_2 = reac_rates[1]*(a-1) # k2 * invariant
                 lambda_t = lambda_1 + lambda_2
                 #print(lambda_1,lambda_2,lambda_t)
                 return (lambda_1/lambda_t) * U_M(q1+1, q2,upper_bounds, k-1, L) + (lambda_2/lambda_t) * U_M(q1,q2+1,upper_bounds, k-1, L)

def inc_sequences(q, upper_bounds, k):
    # q := list of coordinate system, upper_bounds := list of upper bounds for each coordinate, 
    # k := number of permitted steps
    '''Starting from some initial values, we compute the probability of reaching an upper bound at along
    each (increasing) coordinate in a Box.'''

    X = [] # Vector that contains recursion solutions.
    for i in range(q[0], upper_bounds[0]): # First coordinate interval
        for j in range(q[1], upper_bounds[1]): # Second coordinate interval
            X.append(U_M(i, j, upper_bounds, k, L)) # Row-wise, we compute the probability at each coordinate in a Box.
    return np.array(X)

def tran_mat(q, upper_bounds):  
    # q := list of coordinate system
    '''Here, we want to build a transition matrix.'''
    # Zeroth, generate a list of intervals.
    interval_list = gen_intervals(q, upper_bounds)  

    # First, we compute the rate of reaction for each reaction in our network.
    def lambda1(s1, s2): # Rate of reaction for reaction one.
        return reac_rates[0]*(init_states[0] + init_states[1] + 2*init_states[2] - s1 - 2*s2)
    def lambda2(s1, s2): # Rate of reaction for reaction two.
        return reac_rates[1]*((init_states[0] + init_states[1] + 2*init_states[2] - s1 - 2*s2)-1)
    # and the total rate of reaction. 
    def lambda_t(q1, q2):
        return lambda1(q1, q2) + lambda2(q1, q2)

    # Build a dictionary associating each coordinate to an index 
    pair_index = {}
    for i in range(len(interval_list)): # The number of coordinates in an interval.
        pair_index[interval_list[i]] = i # Store in the dicionary index: coordinate pairs.

    # A square matrix (nxn dimensions).
    transition_matrix = np.zeros((len(pair_index), len(pair_index)))

    # Assign to each element in the square matrix a propensity for each respective coordinate.
    for i in range(len(interval_list) - 1): # Once we reach a goal, we stop.
        # i is the value of each index and refers to a coordinate.
        ref = (interval_list[i][0],interval_list[i][1]) # Our reference coordinate.
        if interval_list[i][0] < upper_bounds[0] and interval_list[i][1] < upper_bounds[1]:
            up = (interval_list[i][0],interval_list[i][1] + 1) # A tuple with movement "up".
            right = (interval_list[i][0]+1, interval_list[i][1]) # A tuple with movement "right".
            if up in pair_index.keys() and right in pair_index.keys(): # Ensure coordinate is in our dictionary.
                transition_matrix[pair_index[ref],pair_index[up]] = (lambda2(pair_index[ref],pair_index[up]))/lambda_t(pair_index[ref],pair_index[up]) # Increase on y (\lambda_C)
                transition_matrix[pair_index[ref],pair_index[right]] = (lambda1(pair_index[ref],pair_index[right]))/lambda_t(pair_index[ref],pair_index[right]) # Increase on x (\lambda_B)
            else:
                continue
    return transition_matrix

def cond_vector(q, upper_bounds, L):
    # upper_bounds := list of upper bounds for each coordinate.
    '''Here we create a row-wise condition vector.'''
    interval_list = gen_intervals(q, upper_bounds)
    B = [] # Vector of conditions corresponding to boundary values.
    if L == True: # Toggle on the x-axis.
        for i in interval_list: # Each coordinate in the list.
            if i[0] == upper_bounds[0] - 1: # Target goal on x-axis.
                B.append(1)
            else:
                B.append(0) # Other
    else: 
        for i in interval_list:
            if i[1] == upper_bounds[1]-1:
                B.append(1)
            else:
                B.append(0) 
    return np.array(B)

def upper_approx(q, upper_bounds, k):
    # q := list of coordinate system, upper_bounds := list of upper bounds for each coordinate, 
    # k := number of permitted steps
    A = tran_mat(q, upper_bounds) # Transition matrix.
    B = cond_vector(q, upper_bounds, L) # Cases vector.
    X_n = inc_sequences(q, upper_bounds, k) # Vector of current state probabilities.
    D = np.matmul(A, X_n) # A x X_n
    return D # A x X_n + B

print(upper_approx(q, upper_bounds, k))

