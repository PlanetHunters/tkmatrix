<p align="center">
  <img width="400px" src="https://github.com/martindevora/tkmatrix/blob/master/images/matrix.jpg?raw=true">
</p>

# MATRIX ToolKit
ToolKit for Multi-phAse Transits Recovery from Injected eXoplanets

MATRIX is an injection-recovery software prepared to create grids of scenarios with a set of periods, radii and epochs
of synthetic transiting exoplanet signals in a provided light curve. The typical injection-recovery execution has
consisted in the usage of 2 dimensional scenarios, where only one epoch (random or hardcoded) was used for each period
and radius. Even though this is a fair approach when the computational power is very limited, this reduces the accuracy
of the experiment as the recoverability might vary depending on the transit epoch and hence, the places where the
transits would appear in the curve. E.g: A transit would be clearly recovered falling in a quite region of the curve
but not easy to recover when it happened close to a flare or a curve gap. A multi-phase analysis can now be done with
MATRIX easily by only setting a few parameters in a configuration file and running one line of code.

Just execute the command below this text. Take into accont that the `user-properties.yaml` file needs to include 
several mandatory options. Please refer to the example files under the 
[example](https://github.com/PlanetHunters/tkmatrix/tree/master/examples/TOI-2257/matrix_properties.yaml) directory.

`python3.8 -m tkmatrix --properties user-properties.yaml`

For the complete set of properties available for use please look at the root 
[properties.yaml](https://github.com/PlanetHunters/tkmatrix/tree/master/tkmatrix/properties.yaml)


## By-products
* a_tls_report.csv: A file containing a csv formatted output given the orbital period, the radius and the epoch besides the outputs with found status, SNR and SDE of the results.
* a_tls_report.png: A file with an automatically generated plot from the csv report. You are free to build your own plot from the report if you feel like the one provided by MATRIX is not good enough for your purposes.
* Injected curves (csv files): In case you want to study the injected curves generated for the recovery, you can set a flag to the tool so that it keeps the files after it finishes. If you don't provide that flag, the files will be removed at the end of the execution.


## Citation
We are planning to write a scientific paper based on the usage of MATRIX. In the meantime, we encourage the users to cite the Software DOI in their research:
```
@MISC{2022zndo...6570831D,
       author = {{D{\'e}vora-Pajares}, Mart{\'\i}n and {Pozuelos}, Francisco J.},
        title = "{MATRIX: Multi-phAse Transits Recovery from Injected eXoplanets}",
     keywords = {exoplanets, transits, injection \& recovery, python},
 howpublished = {Zenodo},
         year = 2022,
        month = may,
          eid = {10.5281/zenodo.6570831},
          doi = {10.5281/zenodo.6570831},
      version = {0.3.17},
    publisher = {Zenodo},
       adsurl = {https://ui.adsabs.harvard.edu/abs/2022zndo...6570831D},
      adsnote = {Provided by the SAO/NASA Astrophysics Data System}
}
```

## Documentation
For more information please visit [https://tkmatrix.readthedocs.io/](https://tkmatrix.readthedocs.io/). 
