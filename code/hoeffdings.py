from random import *
from math import *

# random.uniform(a, b)
# Return a random floating point number N such that a <= N <= b for a <= b and b <= N <= a for b < a.

# The end-point value b may or may not be included in the range depending on floating-point rounding in the equation a + (b-a) * random().
def dartboard_simulation(n):
    seed()
    hits = 0
    for i in range(n):
        x = uniform(-1,1)
        y = uniform(-1,1)
        if x**2 + y**2 <= 1:
            hits += 1
    print("π ≈ " + str(hits/n * 4))
    return hits/n

def hoeffdings(p_success, delta):
    p_failure = 1.0 - p_success
    #log() function in python is natural logarithm
    n = log(2.0/p_failure)/(2.0*delta**2)
    return dartboard_simulation(ceil(n))

def hoeffdings_test(p_success, delta):
    successes = 0.0
    failures = 0.0
    for i in range(90):
        if pi/4 - delta/2 < hoeffdings(p_success, delta) < pi/4 + delta/2:
            successes += 1
        else:
            failures += 1
    print("Actual success rate is " + str(successes/(successes + failures)))


hoeffdings_test(2.0/3, .05)
        

# def bisection():
#     lower = 0
#     upper = 4
#     delta = .1 #temp
#     while(True):
#         pi_approx = estimate_pi_monte_carlo(20)
#         print(pi_approx)
#         if isLeft(lower, pi_approx, upper):
#             # print("left")
#             upper = ((lower+upper)/2)
#         else:
#             # print("right")
#             lower = ((lower+upper)/2)
           
#         if abs(upper - lower) < delta:
#             # print(str(lower) + ", " + str(upper))
#             print(pi_approx)
#             return pi_approx
   
   
# def isLeft(lower, pi_approx, upper):
#     return lower<=pi_approx<=((lower+upper)/2)

# bisection()