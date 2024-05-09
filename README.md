# Tracer-dispersion-model
EHS model fit to observed travel distances and Basemodel (MS1) code for the article "Longitudinal fluvial dispersion of coarse particles: Insights from field observations and model simulations".

The raw dataset used in this work are publicly available as supplementary information in Bradley (2017) and downloadable at https://doi.org/10.1002/2017GL075045.

The files 2007, 2008,2009,......2015 contains the tracer ID, tracer position in the current and previous year, distance travelled in the current year.

The file "dispersion time" contains the year and the corresponding dispersion times.

The filename "EHS_model_fit_multiple_years.py" uses the "dispersion time", "2007", "2008", ......"2015" to output the EHS model parameters in the file "EHS_parameters_multiple_years.csv"

The filename "EHS_parameters_multiple_years.csv" serves as input for the Monte Carlo simulations in the basemodel scenario (MS1).
