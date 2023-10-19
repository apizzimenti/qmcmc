import math

def chernoff_estimation(delta):
    count_success = 0
    total_trials = 0
    
    while True:
        total_trials += 1
        # Simulate generating a random sample (success or failure)
        # For simplicity, let's assume a success rate of 0.7
        sample_success = (total_trials % 10) < 7
        
        if sample_success:
            count_success += 1
        
        estimated_success_rate = count_success / total_trials
        upper_bound = math.exp(-2 * total_trials * ((estimated_success_rate - 2/3)**2))
        
        if upper_bound <= delta:
            return estimated_success_rate

# Testing the algorithm with delta = 0.1
delta = 0.1
estimated_rate = chernoff_estimation(delta)
print("Estimated success rate:", estimated_rate)
