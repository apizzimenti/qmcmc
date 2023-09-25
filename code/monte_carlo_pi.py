import random
import math
import time
import matplotlib.pyplot as plt

def estimate_pi(num_samples):
    inside_circle = 0
    x_inside = []
    y_inside = []
    x_outside = []
    y_outside = []

    for _ in range(num_samples):
        x = random.random()
        y = random.random()
        distance = math.sqrt(x**2 + y**2)

        if distance <= 1:
            inside_circle += 1
            x_inside.append(x)
            y_inside.append(y)
        else:
            x_outside.append(x)
            y_outside.append(y)

    pi_estimate = (inside_circle / num_samples) * 4
    return pi_estimate, x_inside, y_inside, x_outside, y_outside

if __name__ == "__main__":
    #num_samples = 100000000  # Number of random samples
    num_samples = 100000

    start_time = time.time()
    pi_estimate, x_inside, y_inside, x_outside, y_outside = estimate_pi(num_samples)
    end_time = time.time()
    execution_time = end_time - start_time

    print(f"Estimated Pi: {pi_estimate}")
    print(f"Execution time: {execution_time} seconds")

    # Visualize the points
    plt.figure(figsize=(8, 8))
    plt.scatter(x_inside, y_inside, color='magenta', marker='.')
    plt.scatter(x_outside, y_outside, color='black', marker='.')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.title('Monte Carlo Estimation of Pi')
    plt.axis('equal')
    plt.show()
