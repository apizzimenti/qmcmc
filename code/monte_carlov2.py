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
    num_samples = 10000  # Number of random samples

    start_time = time.time()
    pi_estimate, x_inside, y_inside, x_outside, y_outside = estimate_pi(num_samples)
    end_time = time.time()
    execution_time = end_time - start_time

    print(f"Estimated Pi: {pi_estimate}")
    print(f"Execution time: {execution_time} seconds")

    # Create a square and a circle
    circle = plt.Circle((0, 0), 1, color='indigo', alpha=0.3)
    square = plt.Rectangle((-1, -1), 2, 2, color='b', fill=False)

    # Visualize the points, circle, and square
    fig, ax = plt.subplots(figsize=(8, 8))
    ax.add_patch(circle)
    ax.add_patch(square)
    ax.scatter(x_inside, y_inside, color='indigo', marker='.')
    ax.scatter(x_outside, y_outside, color='black', marker='.')
    ax.set_xlim([-1.5, 1.5])
    ax.set_ylim([-1.5, 1.5])
    ax.set_aspect('equal')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.title('Monte Carlo Estimation of Pi')
    plt.show()
