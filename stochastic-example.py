import numpy as np
import matplotlib.pyplot as plt

plt.rcParams['text.usetex'] = True

#q = [100, 5, 10] # Initial states for A, B, and C respectively.
k = [20, 1] # Kinetic rates for k1 and k2, respectively.

def propensity_function1(k_1, k_2, N_A):
      return (k_1*N_A)/((k_1*N_A) + (k_2*N_A*(N_A-1)))

def propensity_function2(k_1, k_2, N_A):
      return (k_2*N_A*(N_A-1))/((k_1*N_A) + (k_2*N_A*(N_A-1)))
      

def markov_chain(q, k):
  # Initialization of vectors containing Markov states.
  x_A = np.zeros(100)
  x_B = np.zeros(100)
  x_C = np.zeros(100)

  # Markov state initial condition.
  x_A[0] = q[0]
  x_B[0] = q[1]
  x_C[0] = q[2]

  ratio = [] # The ratio between competing molecules.
  for i in range(len(x_A)-1):
    a = (q[0] + q[1] + 2*q[2] - x_B[i] - 2*x_C[i]) # Linear invariant.
    # Reaction propensities
    lambda_1 = k[0] 
    lambda_2 = k[1]*(a-1) 
    lambda_t = lambda_1 + lambda_2
    if (0 < x_A[i]) or (x_B[i] == 0 or x_C[i] == 0): # Resource availability.
      birth_B = np.random.rand() <= (lambda_1 / lambda_t)* x_B[i] # Production of molecule B.
      birth_C = np.random.rand() <= (lambda_2 / lambda_t)* x_C[i] # Production of molecule C.
      # Store the values along each Markov chain.
      x_B[i+1] = x_B[i] + 1*birth_B 
      x_C[i+1] = x_C[i] + 1*birth_C
      x_A[i+1] = x_A[i] - 1*birth_B - 2*birth_C 
    ratio.append(x_C[i]/x_B[i]) 

  solution_A = [] # Amount of resources left.
  for i in range(len(x_A)-1):
    if x_A[i] > 0:
      solution_A.append(x_A[i])

  return [solution_A, ratio]

prop_1 = propensity_function1(20, 1, 4)
prop_2 = propensity_function2(20, 1, 4)
propensity_total = prop_1 + prop_2 
prob_trace = prop_2 / propensity_total 
print(prob_trace)

# instance_B = markov_chain([100, 20, 5], k)
# instance_C = markov_chain([100, 5, 20], k)
# plt.plot(instance_B[0], instance_B[1][:len(instance_B[0])])
# plt.plot(instance_C[0], instance_C[1][:len(instance_C[0])])
# plt.xlabel('Availability of A')
# plt.ylabel(r'$N_C / N_B$')
# plt.title('Competition for resources between molecules B and C.')
# plt.xlim(0, 100)
# plt.ylim(0, 4)
#plt.savefig('resource-comp-dynamics.png')
#plt.show()

'''Explanation for the method:
Here we examine the competition between two reactions for resources. The common reactant is molecule A, and one 
instance of A is required for the first reaction, while 2 instances are required for the second.
The propensity of each reaction gives us information about frequencies, and if done correct it will
provide information about competition between reactions.
'''