<p align="center">
  <img width="400px" src="https://github.com/martindevora/tkmatrix/blob/master/images/matrix.jpg?raw=true">
</p>

# MATRIX ToolKit
ToolKit for Multi-phAse Transits Recovery from Injected eXoplanets

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

## Main Developers
[M. Dévora Pajares](https://github.com/martindevora)

[F.J. Pozuelos](https://github.com/franpoz)


## Additional contributors
[L. Cerdeño Mota](https://github.com/LuisCerdenoMota) 

## Installation
Supported Python versions: 3.8, 3.9. Install with:

```
python3.8 -m pip install numpy==1.22.4
python3.8 -m pip install -r requirements.txt
```

You can find the requirements.txt file [here](https://github.com/PlanetHunters/tkmatrix/blob/master/requirements.txt). 

## Tests
We use [tox](https://tox.readthedocs.io) to test MATRIX under all the supported Python versions. Usage:

`tox`

## Examples
Under the [examples](https://github.com/PlanetHunters/tkmatrix/tree/master/examples) directory.

## Execution
Just execute the command below this text. Take into accont that the `user-properties.yaml` file needs to include several mandatory options. Please refer to the example file under the examples directory.

`python3.8 -m tkmatrix --properties user-properties.yaml`

## By-products
* a_tls_report.csv: A file containing a csv formatted output given the orbital period, the radius and the epoch besides the outputs with found status, SNR and SDE of the results.
* a_tls_report.png: A file with an automatically generated plot from the csv report. You are free to build your own plot from the report if you feel like the one provided by MATRIX is not good enough for your purposes.
* Injected curves (csv files): In case you want to study the injected curves generated for the recovery, you can set a flag to the tool so that it keeps the files after it finishes. If you don't provide that flag, the files will be removed at the end of the execution.
