# Tracer-dispersion-model
EHS model fit to observed travel distances and Basemodel (MS1) code for the article "Longitudinal fluvial dispersion of coarse particles: Insights from field observations and model simulations".

The raw dataset used in this work are publicly available as supplementary information in Bradley (2017) and downloadable at https://doi.org/10.1002/2017GL075045.

The files 2007, 2008, 2009,......2015 contains the tracer ID, tracer position in the previous year and current year, distance travelled in the current year.

The file "dispersion time" contains the year and the corresponding dispersion times, which serves as input for the EHS model fits.

The filename "EHS_model_fit_multiple_years.py" uses the "dispersion time", and travel distances of each grain from the files "2007", "2008", ......"2015" to give the EHS model parameters in the file "EHS_parameters_multiple_years.csv"

The filename "EHS_parameters_multiple_years.csv" serves as input for the Monte Carlo simulations in the basemodel scenario (Basemodel_MS1) to give travel distance distribution of each flood event.
