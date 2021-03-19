<p align="center">
  <img width="400px" src="https://github.com/martindevora/tkmatrix/blob/master/images/matrix.jpg?raw=true">
</p>

# MATRIX ToolKit
ToolKit for Multi-phAse Transits Recovery from Injected eXoplanets

## Main Developers
[M. Dévora Pajares](https://github.com/martindevora)

[F.J. Pozuelos](https://github.com/franpoz)

[L. Cerdeño Mota](https://github.com/LuisCerdenoMota) 

## Additional contributors 
<i>A. Thuillier</i>

## Installation
Supported Python versions: 3.6, 3.7, 3.8. Install with:

`python3 -m pip install tkmatrix`

## Tests
We use [tox](https://tox.readthedocs.io) to test MATRIX under all the supported Python versions. Usage:

`tox`

## Examples
Under the examples directory.

## Execution
Just execute the command below this text. Take into accont that the `user-properties.yaml` file needs to include several mandatory options. Please refer to the example file under the examples directory.

`python3 -m tkmatrix --properties user-properties.yaml`

## By-products
* a_tls_report.csv: A file containing a csv formatted output given the orbital period, the radius and the epoch besides the outputs with found status, SNR and SDE of the results.
* a_tls_report.png: A file with an automatically generated plot from the csv report. You are free to build your own plot from the report if you feel like the one provided by MATRIX is not good enough for your purposes.
* Injected curves (csv files): In case you want to study the injected curves generated for the recovery, you can set a flag to the tool so that it keeps the files after it finishes. If you don't provide that flag, the files will be removed at the end of the execution.
