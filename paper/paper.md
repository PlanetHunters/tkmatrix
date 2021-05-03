---
title: '``MATRIX``: Multi-phAse Transits Recovery from Injected eXoplanets.'
tags:
  - Python
  - Astronomy
  - Exoplanets
  - Kepler
  - K2
  - TESS
authors:
  - name: Martín Dévora-Pajares(*)
    orcid: 0000-0003-2007-6144 
    affiliation: 1
  - name: Francisco J. Pozuelos(*)
    orcid: 0000-0003-1572-7707
    affiliation: "2, 3" # (Multiple affiliations must be quoted)
  - name: Antoine Thuillier
    orcid: 
    affiliation: 2 # (Multiple affiliations must be quoted)
  - name: Valérie Van Grootel
    affiliation: 2
  - name: Luis Cerdeño Mota
    orcid:
affiliations:
 - name: Dpto. Física Teórica y del Cosmos, Universidad de Granada, 18071, Granada, Spain
   index: 1
 - name: Space Sciences, Technologies and Astrophysics Research (STAR) Institute, Universitè de Liège, Allée du 6 Août 19C, B-4000 Liège, Belgium
   index: 2
 - name: Astrobiology Research Unit, Universitè de Liège, Allée du 6 Août 19C, B-4000 Liège, Belgium
   index: 3
date: 22 January 2021
bibliography: paper.bib

# Optional fields if submitting to a AAS journal too, see this blog post:
# https://blog.joss.theoj.org/2018/12/a-new-collaboration-with-aas-publishing
# aas-doi: 10.3847/xxxxx <- update this with the DOI from AAS once you know it.
# aas-journal: Astrophysical Journal <- The name of the AAS journal.
---

# Summary

The exoplanets detection is a field that has been growing exponentially in the last decades. More than
half of the known exoplanets have been discovered thanks to the usage of the transits method. It consists
on the measurement of the stellar light flux drops caused by a planet passing through the visual line
between the star and the observer. So far, several space-based missions have already been launched
as COROT, Kepler and TESS, producing huge amounts of publicly available data and growing communities
of scientists looking forward to them. Every star observation shows its own systematics, together
with the common expected ones, due to its own characteristics and its nearby field. Therefore,
it has been a common task for scientists to create exoplanet transiting models for each analyzed star
and try their search tools on them to define their detection limits. We present an extension to the 
typical model injection and recovery scenarios by adding a new dimension to the analysis, establishing
in most of the cases more confident detection thresholds.

# 1. MATRIX
The main data product that astronomers use to analyze and search for transiting exoplanets are the light
curves. These are time series for stellar flux measurements. One of the main methods used to search 
transiting exoplanets under light curves in the community is the Box Least Squares algorithm. 
This method folds the light curve with many different periods and tries to fit a squared model of a transit
to each of the folded time series. The folded signal with less residuals is chosen as the best candidate
to a transiting exoplanet and its main parameters (period, depth, epoch, duration) are returned. In the last
two years, a new method appeared named Transit Least Squares, aimed to challenge the BLS results and providing
better residuals thanks to the usage of realistic transit models for its fit stages instead of a squared one.
MATRIX will use TLS by default, though the user could switch to BLS if desired.

The user will provide a grid of periods by selecting the `MAX_PERIOD`, `MIN_PERIOD` and `STEP_PERIOD`, a grid
of planet radius by selecting `MAX_RADIUS`, `MIN_RADIUS` and `STEP_RADIUS` and the number of epochs to be explored
for each case. `MATRIX` will initially download the target star light curve (or load a `csv` file specified by the 
user), generate a model of transit injecting it into the original data and store the resultant modelled
light curve in a `csv` file for each case. That is, `MATRIX` will store a set of 
`PERIOD_GRID_SIZE` x `RADIUS_GRID_SIZE` x `EPOCHS_COUNT` files for the recovery scenario.

# 2. Scientific cases 

The ``MATRIX ToolKit`` is specially designed to help the scientists develop a robust inject and recovery 
analysis on a given target star. This is needed in almost every exploration project that aims to define
the strength of its search tools and pipelines.

## Single-phase inject and recovery
The most usual way to study the ability of a given exoplanet search tool of finding new candidates is the
launch of an inject and recovery process for a grid of periods and planet radius. For this traditional case, 
MATRIX provides an easy-to-use execution command which only needs to be fed with a YAML file including the 
scenario parameters 
(see [mono-phase.yaml](https://github.com/martindevora/matrix/blob/master/examples/mono-phase.yaml) file.). You
can appreciate that the `PHASES` property is set to `1`.

## Multi-phase inject and recovery
In many cases, a single-phase inject and recovery scenario proves to show poor results near the threshold
detection limit (when transiting exoplanet models start to be difficult to be detected) because it only assess
the detectability of a model with one sample for a given period and radius pair. But there is also one more 
reason that could reverse completely the results of a single-phase inject and recovery scenario: the epoch
of the modeled transiting exoplanets becomes crucial. In case the selected epoch makes the transit events 
appear under noisy regions or data gaps, they will become much more difficult to detect. This situation is very complicated
to correct and therefore, we have added a new dimension to the inject and recovery scenarios: a grid of epochs
for each period and radius cell. By setting the `PHASES` parameter to any value greater than `1` we will be running
an injection and recovery analysis on as many different epochs for each period and radius within the grid as phases
we determined.

HERE A PICTURE OF Single-phase VS Multi-Phase

# 3. Supported inputs for targets and tools
We have included several kind of inputs for `MATRIX`. Missions like TESS provide short and long cadence data and their study
differ in several points. Thus, the injection and recovery task might be useful for each of them separately. Sometimes the
astronomer will be provided with a file containing the target light curve and sometimes the astronomer might be
interested in specifying the target star properties manually instead of loading those from the online catalogues. Finally,
the scientist could not only want to use the `MATRIX` built-in search method but also use his cleaning and search
algorithms to execute an injection and recovery process according to his real search procedures instead of using the
ones we adopted as default. We explain some additions we included in our new tool to support all of these possibilities.

## Missions and exposure times
Support for Kepler, K2 and TESS missions online data is included. By only specifying the field `TARGET` with the appropiate
target name (e.g. KIC 1234, EPIC 1234 or TIC 1234 for Kepler, K2 or TESS targets). By default, the short cadence exposure time
will be selected (60 seconds for Kepler and K2, 120 seconds for TESS) but it can be switched to the long cadence by changing
the value of the `EXPOSURE_TIME` parameter (long cadence values are 1800 seconds for the three missions).

## Star properties
The target star properties can be loaded automatically either from the online catalogues when the `TARGET` name is specified or
providing a `STAR_INFO` root property containing child properties for each star property like mass, radius, etc.

## File input
In addition to setting the target to be searched for, a `FILE` parameter can be provisioned pointing to the absolute file path
to be loaded in `csv` format. This way, the star properties will be downloaded from the online catalogues and the light curve
used to inject the models will be loaded from disk. If no `TARGET` value is provided, the `STAR_INFO` property needs
to be included in order to be able to generate a realistic model according to the expected host star.

## Grids configuration
Both the period and radius grids can be generated linearly or logarithmically (by setting `PERIOD_GRID_GEOM` or `RADIUS_GRID_GEOM`
to `linear` or `log`) In case `linear` is selected, the grid will be generated by sampling periods or radiuses equally spaced given
`PERIOD_STEP` or `RADIUS_STEP` respectively. In case `log` is selected, `PERIOD_STEP` or `RADIUS_STEP` will define the number of
points to be hold by the respective grid.

## Rotation signal clean-up
## Grid mode
## Custom algorithms

## Missions and inputs

# 4. Performance
COMPARISON BETWEEN SINGLE vs MULTI FOR TWO TARGETS
SNR AND SDE DATA FOR SHALLOW TRANSITS



# 5. Future implementations  

- Attention to threshold regions: As we are using squared grids and the inject and recovery scenarios focus
on finding a detection threshold, some wide regions are usually showing the same detection values (found or
not found), which might represent a waste of computational power. To mitigate this, we plan to incorporate some
kind of attention mechanism into our algorithm in such a way that it could only keep testing the scenario only
near the detection limits, assuming that above them the results are true (found) and below them they are 
false (not found).

# Acknowledgements
We acknowledge support from the Spanish Ministry of Science through the project PID2019-107061GB-C64/SRA (State Research Agency/10.13039/501100011033)

# References
