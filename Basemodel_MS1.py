import pandas as pd
import numpy as np
from scipy import optimize
from scipy.special import iv
import random


# Define the EHS model pdf
def ehs_pdf(x, t, k1, k2):
    return k1 * np.exp(-(k1 * x) - (k2 * t)) * np.sqrt((k2 * t) / (k1 * x)) * iv(1, 2 * np.sqrt(k1 * x * k2 * t))

# Define the inverse of the cumulative distribution function (CDF) of the pdf (to draw random values for travel distances)
#inverse transform sampling method
def inverse_cdf(t, k1, k2, num_samples=1000):
    # Define the range for x values
    x_values = np.linspace(0.001, 2000, 5000)  # Adjust the range as per your requirement
    # Calculate the pdf values
    pdf_values = ehs_pdf(x_values, t, k1, k2)
    # Calculate the cumulative distribution function (CDF)
    cdf = np.cumsum(pdf_values)
    cdf = cdf/cdf[-1]  # Normalize the CDF to ensure it ranges from 0 to 1
    # Generate random samples from the uniform distribution
    u = np.random.rand(num_samples)
    # Find the corresponding x values using the inverse CDF
    x_samples = np.interp(u, cdf, x_values)
    return x_samples

# Load the probability distribution parameters CSV file
params_df = pd.read_csv('EHS_parameters_multiple_years.csv')

# List of years for which we want to calculate travel distances and positions
years = [2007, 2008, 2009, 2010, 2011, 2013, 2014, 2015]
#2012 year results in 0 value for cdf and divide by zero creates errors

# Create a DataFrame with 1000 fake tracer IDs and a position column
fake_tracers = pd.DataFrame({'ID': range(1, 1001), 'Position_2006': np.zeros(1000)})

# Number of simulations to run
num_simulations = 10000

# Loop through each simulation
for sim in range(num_simulations):
    # Shuffle the order of years randomly
    
    shuffled_years = [2007, 2008, 2009, 2010, 2011, 2013, 2014, 2015] #Using same sequence
    # Create a new DataFrame for the current simulation
    sim_results = pd.DataFrame({'ID': range(1, 1001)})
    
    # moved_both_indices = []
    # previously_moved_indices = []
    
    # Loop through each year and calculate travel distances and positions for each fake tracer ID
    for i, year in enumerate(shuffled_years):
        print(f'Simulation {sim+1} of {num_simulations}')
        
        # Extract alpha, loc, beta, and movement probability values for the current year
        k1 = params_df.loc[params_df['Year'] == year, 'k1'].values[0]
        k2 = params_df.loc[params_df['Year'] == year, 'k2'].values[0]
        t = params_df.loc[params_df['Year'] == year, 'Dispersion time'].values[0]
        movement_probability = params_df.loc[params_df['Year'] == year, 'Movement Probability'].values[0]
        
                
        # Create an array with zeros for calculated travel distance and position
        cal_distances = np.zeros(1000, dtype=float)
        cal_positions = np.zeros(1000, dtype=float)
        
        # Draw whether each tracer moved or not based on movement probability
        moved = np.random.choice([0, 1], p=[1-movement_probability, movement_probability], size=1000)
        
        
        # Map the probabilities to the desired range using the inverse transform sampling,.i.e., get distances only for the moved grains
        distances = np.where(moved == 1, inverse_cdf(t,k1,k2,num_samples=1000),0)
                
              
        # Calculate new distances and positions for each tracer
        cal_distances = distances
        if i == 0:
            cal_positions = distances
        else:
            cal_positions = sim_results[f'Position_{shuffled_years[i-1]}_{sim+1}'] + distances
        
        
        
        # Write the calculated distances and positions to the sim_results DataFrame for the current year
        sim_results[f'Distance_{year}_{sim+1}'] = cal_distances
        sim_results[f'Position_{year}_{sim+1}'] = cal_positions
    
    # Write sim_results to CSV after each simulation
    sim_results.to_csv(f'MC_EHS{sim+1}.csv', index=False)
