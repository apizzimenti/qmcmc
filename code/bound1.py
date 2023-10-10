import numpy as np
import matplotlib.pyplot as plt

def estimate_pi_monte_carlo(num_samples):
    inside_circle = 0
    for _ in range(num_samples):
        x, y = np.random.rand(2)  # Generate random points (x, y)
        if x**2 + y**2 <= 1:      # Check if the point is inside the quarter circle
            inside_circle += 1

    estimated_pi = 4 * inside_circle / num_samples
    return estimated_pi

def chernoff_bound_error(epsilon, delta, num_samples):
    return np.sqrt((1 / (2 * num_samples)) * np.log(2 / delta) / (2 * epsilon))

# Parameters for Chernoff bound
epsilon = 0.1  # Desired precision
delta = 0.05   # Probability of exceeding the error

# Number of Monte Carlo samples
num_samples_list = np.arange(100, 5000, 100)

# Estimate π using Monte Carlo and calculate Chernoff bound error
estimated_pis = []
chernoff_bound_errors = []

for num_samples in num_samples_list:
    estimated_pi = estimate_pi_monte_carlo(num_samples)
    estimated_pis.append(estimated_pi)

    # Calculate Chernoff bound error
    error = chernoff_bound_error(epsilon, delta, num_samples)
    chernoff_bound_errors.append(error)

# Plot the results
plt.figure(figsize=(10, 6))
plt.plot(num_samples_list, estimated_pis, label='Estimated π')
plt.plot(num_samples_list, [np.pi] * len(num_samples_list), linestyle='--', label='True π')
plt.plot(num_samples_list, [np.pi + error for error in chernoff_bound_errors], linestyle=':', label='π + Chernoff Bound Error')
plt.plot(num_samples_list, [np.pi - error for error in chernoff_bound_errors], linestyle=':', label='π - Chernoff Bound Error')
plt.xlabel('Number of Samples')
plt.ylabel('Value')
plt.legend()
plt.title('Estimating π using Monte Carlo and Chernoff Bound')
plt.show()
