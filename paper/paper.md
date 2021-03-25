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
  - name: Luis Cerdeño Mota
    orcid:
  - name: Antoine Thuillier
    orcid: 
    affiliation: 2 # (Multiple affiliations must be quoted)
  - name: Valérie Van Grootel
    affiliation: 2
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
and try their search tools on them to define their detection limits. As far as we know there is no
public software to perform such a process and thus, we developed MATRIX (Multi-phAse Transits
Recovery of Injected eXoplanets) to help the astronomers assess the recovery rates of exoplanet
transits around a given star.

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

The user will provide a grid of periods by selecting the MAX_PERIOD, MIN_PERIOD and STEP_PERIOD, a grid
of planet radius by selecting MAX_RADIUS, MIN_RADIUS and STEP_RADIUS and the number of epochs to be explored
for each case. MATRIX will initially download the target star light curve (or load a CSV file specified by the 
user), generate a model of transit injecting it into the original data and store the resultant modelled
light curve in a CSV file for each case. That is, MATRIX will store a set of 
PERIOD_GRID_SIZE x RADIUS_GRID_SIZE x EPOCHS_COUNT files for the recovery scenario.

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
appear under noisy regions, they will become much more difficult to detect. This situation is very complicated
to correct and therefore, we have added a new dimension to the inject and recovery scenarios: a grid of epochs
for each period and radius cell.

# 4. Performance
Comparison between mono-phase and multi-phase search for some target.


# 5. Future implementations  

- Attention to threshold regions: As we are using squared grids and the inject and recovery scenarios focus
on finding a detection threshold, some wide regions are usually showing the same detection values (found or
not found), which might represent a waste of computational power. To mitigate this, we plan to incorporate some
kind of attention mechanism into our algorithm in such a way that it could only keep testing the scenario only
near the detection limits, assuming that above them the results are true (found) and below them they are 
false (not found).

- Custom user search algorithms: Instead of using BLS or TLS, the user might be interested in testing his own
search pipeline. Therefore we are aimed to implement a flexible piece of code where the user could provide
the desired custom search algorithm.

# Acknowledgements


# References
