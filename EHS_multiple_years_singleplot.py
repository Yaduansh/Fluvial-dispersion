import numpy as np
from scipy.special import iv
from scipy.optimize import brute
import pandas as pd
import matplotlib.pyplot as plt

# Define the EHS model pdf
def ehs_pdf(x, t, k1, k2):
    return k1 * np.exp(-(k1 * x) - (k2 * t)) * np.sqrt((k2 * t) / (k1 * x)) * iv(1, 2 * np.sqrt(k1 * x * k2 * t))

# Define the sum of squared errors between observed and predicted probabilities
def sse(params, x_bins, obs_prob):
    k1, k2 = params
    pred_prob = np.zeros_like(obs_prob)
    for i in range(1, len(x_bins)):
        x = x_bins[i]
        pred_prob[i] = ehs_pdf(x, t, k1, k2)
    pred_prob[0] = np.exp(-(k2 * t))  # Add e^(-k2*t) to the first bin
    pred_prob = pred_prob / np.sum(pred_prob)
    return np.sum((obs_prob[1:] - pred_prob[1:]) ** 2)  # Exclude the first bin when calculating SSE


# Read dispersion time from a CSV file
dispersion_df = pd.read_csv("Dispersion time.csv")
dispersion_time = dispersion_df.set_index("Year")["Time"]

# Create a list to store the results for each year
results_list = []

# Create a single figure and axis
fig, ax = plt.subplots()

years = [2007, 2008, 2009, 2010, 2011, 2013, 2014, 2015]

# Iterate over the years from 2008 to 2015
for year, year in enumerate(years):
    # Read the data for the current year
    filename = str(year) + ".csv"
    data = pd.read_csv(filename)
    x = data["Dist"]

    # Divide the travel distances into 20 bins
    x_bins = np.linspace(0, max(x), 20)
    x_centers = (x_bins[1:] + x_bins[:-1]) / 2

    # Compute observed probabilities
    obs_prob, _ = np.histogram(x, bins=x_bins)
    obs_prob = obs_prob / np.sum(obs_prob)

    # Set dispersion time for the current year
    t = dispersion_time.loc[year]

    # Perform brute force optimization
    params_range = (slice(0.001, 1, 0.001), slice(0.000001, 0.1, 0.001))
    res = brute(sse, params_range, args=(x_centers[1:], obs_prob[1:]), full_output=True, finish=None)

    # Extract estimated parameters
    if isinstance(res, tuple):
        k1_est, k2_est = res[0]
        r2 = 1 - res[1] / np.sum((obs_prob[1:] - np.mean(obs_prob[1:])) ** 2)
        print(f"Year: {year}")
        print(f"Dispersion time: {t}")
        print(f"Estimated k1: {k1_est}")
        print(f"Estimated k2: {k2_est}")
        print(f"R2 value: {r2}")
        print()

        # Append the results to the list
        results_list.append({"Year": year, "Dispersion time": t, "k1": k1_est, "k2": k2_est, "R2": r2})

        # Compute predicted probabilities using estimated parameters
        pred_prob = np.zeros_like(obs_prob)
        for i in range(len(x_centers)):
            x = x_centers[i]
            pred_prob[i] = ehs_pdf(x, t, k1_est, k2_est)
        pred_prob[0] = np.exp(-(k2_est * t))  # Add e^(-k2*t) to the first bin
        pred_prob = pred_prob / np.sum(pred_prob)

        # Plot the fitted pdf on the same figure
        ax.scatter(x_centers, obs_prob, alpha=0.5)
        ax.plot(x_centers, pred_prob, '--', label=f'{year}')

# Set plot settings and labels
ax.set_yscale('log')
ax.set_xscale('log')
ax.set_ylim([10 ** -3, 10 ** 0])
ax.set_xlim([10 ** 0, 10 ** 3])
ax.set_title("EHS")
ax.set_xlabel("Travel Distance (m)")
ax.set_ylabel("Probability")
ax.legend()
plt.tight_layout()

# Save the figure
plt.savefig("EHS_all_years.png", dpi=600)

# Show the figure
plt.show()

# Convert the results list to a pandas dataframe
results_df = pd.DataFrame(results_list)

# Save the results to a CSV file
results_df.to_csv("EHS_parameters_multiple_years.csv", index=False)
