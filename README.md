# Near-real-time-prompt-assessment-for-regional-EQIL
Near-real-time prompt assessment for regional EQIL is an open-source package providing a method to perform the near-real-time prompt assessment for regional earthquake-induced landslides
It was developed in python, C++ and MATLAB language. And MATLAB code is provided by Montejo and Suarez (2013), which is compatible with MATLAB R2014a and subsequent versions

To launch Near-real-time prompt assessment for regional EQIL:
- Use python to run the script Four.py to calculate the response spectrum of the ground motion reocrds.
- Use C++ program InteStation to calculate the response spectrum of the target location.
- Use MATLAB platform to run the program Ground motion generation using CWT to generate the ground motion of the target location.
- Use C++ program NewmarkDispStation to calculate the landslide displacement of the target location. 

Preparing files:
- input.txt: This file includes the name of ground motion records that used to interpolate the ground motion at the target location.
- station.txt: This file includes the longitude and latitude of the target location.
- groundmotion: This folder includes the ground motion records that will be used for the interpolation. 

Please report any error, bug or suggestion to chengqingle@gmail.com

Reference cited:

Cheng, Q.L., Tian, Y., Lu, X.Z., Huang, Y. L., Ye, L. P., 2021. Near-real-time prompt assessment for regional earthquake-induced landslides using recorded ground motions. Computers & Geosciences. doi:10.1016/j.cageo.2021.104709.

Lu, X.Z., Cheng, Q.L., Tian, Y., Huang, Y.L., 2021. Regional ground‐motion simulation using recorded ground motions. Bulletin of the Seismological Society of America. doi:10.1785/0120200243.

Montejo, L. A., Suarez, L. E., 2013. An improved CWT-based algorithm for the generation of spectrum-compatible records. International Journal of Advanced Structural Engineering 5(1), 26. DOI: 10.1186/2008-6695-5-26.
