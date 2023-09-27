import matplotlib.pyplot as plt
import numpy as np

def calculate_pi_monte_carlo(num_samples):
    inside_circle = 0
    for _ in range(num_samples):
        x, y = np.random.random(2) * 2 - 1  # Generate random point in [-1, 1] x [-1, 1]
        if x**2 + y**2 <= 1:
            inside_circle += 1

    return (inside_circle / num_samples) * 4

def chernoff_bound(n, delta):
    return np.sqrt((1 / (2 * n)) * np.log(2 / delta))

# Number of Monte Carlo samples
num_samples = 10000

# True value of pi
true_pi = np.pi

# Calculate Monte Carlo estimate of pi
pi_estimate = calculate_pi_monte_carlo(num_samples)

# Chernoff bound calculation
n_values = np.arange(1, num_samples + 1)
chernoff_bounds = [chernoff_bound(n, 0.1) for n in n_values]

# Plot the results
plt.figure(figsize=(10, 6))
plt.axhline(y=true_pi, color='r', linestyle='--', label='True Pi Value')
plt.plot(n_values, np.full_like(n_values, true_pi), label='True Pi Value')

plt.plot(n_values, [pi_estimate] * len(n_values), label='Monte Carlo Estimate of Pi')

plt.plot(n_values, pi_estimate + np.array(chernoff_bounds), 'g--', label='Chernoff Bound')
plt.plot(n_values, pi_estimate - np.array(chernoff_bounds), 'g--')

plt.xlabel('Number of Monte Carlo Samples')
plt.ylabel('Estimated Value of Pi')
plt.title('Monte Carlo Estimation of Pi with Chernoff Bound')
plt.legend()
plt.show()
