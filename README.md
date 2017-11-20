Reference:
Jones D. E., Stenning D. C., Ford E. B., Wolpert R. L., Loredo T. J., Dumusque, X. (2017+). Improving Exoplanet Detection Power: Multivariate Gaussian Process Models for Stellar Activity. URL: https://arxiv.org/abs/1711.01318

This code implements the methodology of Jones et al. (2017+):
i) Customized PCA for constructing stellar activity indicators
ii) Maximum likelihood fitting of multivariate GP models for stellar activity 

Note: the code is currently being organized into a more user-friendly format. 

The main files run files are: 
- raw data and PCA procedure/make_GPCA_output_clean_final.R
- gp model fitting/fine_planets_datasets_analysis.R
- gp model fitting/null_datasets_analysis.R

Users should cite Jones et al. (2017+) in their work. 

The raw data was generated using SOAP 2.0, see:
- Dumusque, X., Boisse, I. and Santos, N. (2014). SOAP 2.0: A tool to estimate the photometric and radial velocity variations induced by stellar spots and plages. The Astrophysical Journal 796 132.
- Boisse, I., Bonfils, X. and Santos, N. (2012). SOAP: A tool for the fast computation of photometry and radial velocity induced by stellar spots. Astronomy & Astrophysics 545 A109.

