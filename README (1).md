# Dust Mapping Project

## Introduction

This project uses Gaia DR3 data to study interstellar dust extinction using the open cluster Cl Tombaugh 4.

While looking at Tombaugh 4's color-magnitude diagram, I noticed a big color shift and dimming compared to what a normal population should look like---a sign of heavy dust extinction along that sightline. That observation is what kicked off this whole project: isolate the cluster, find its red clump stars (great standard candles), and use them to measure extinction and reddening as a function of distance.

## How It Works

1. **Cluster selection**: Sky position and proper motion cuts (sigma-clipped, narrowed round by round) isolate likely Tombaugh 4 members from the surrounding field.
2. **Galactic component classification**: Stars are transformed into Galactocentric coordinates so we can check their height above the Galactic plane and classify them as thin disk, thick disk, or halo.
3. **Red clump identification**: A color-magnitude box (plus the disk classification) picks out red clump stars, which have a nearly universal absolute magnitude---making them reliable standard candles.
4. **Extinction & reddening**: Using the red clump's known intrinsic color/magnitude (Ruiz-Dern et al. 2018), we compute color excess and extinction (A_V, A_G).
5. **3D dust mapping**: Combining extinction with 3D positions gives spatial maps of dust along the line of sight.

## Important Takeaways

- The original version of this analysis went star-by-star in a messy, repetitive way (copy-pasted selection blocks for every round). It's been refactored into `cluster_selection.py`, using classes (`ClusterSelector`, `SelectionRound`) so the same selection logic gets reused instead of retyped.
- Sigma-clipping by proper motion only makes sense for populations that actually share a common motion (like a real cluster). Applying it to field-star selections (like the red clump box cut) doesn't isolate anything physical---worth knowing before slapping a pm cut on everything.
- Apparent magnitude and absolute magnitude are easy to accidentally mix up across a big notebook, especially when comparing to reference values from a paper. Always double check which one you're plotting!
- Masked Gaia columns (like `bp_rp`) can silently break plotting functions even when `np.isfinite()` says everything's fine---worth explicitly stripping masks before handing data to something like `hist2d`.

## Outputs

The code produces:
- Color-magnitude diagrams (cluster rounds, background field, red clump)
- Thin disk / thick disk / halo classification + RA/Dec maps
- Extinction vs. distance plots
- 2D and 3D dust density/extinction maps (static + interactive)

## Repository Structure

```
cluster_selection.py    # Reusable ClusterSelector / SelectionRound classes for
                         # Gaia-based kinematic selection and galactic
                         # component classification
Dust_Mapping.ipynb      # Main analysis notebook: selection, CMDs, extinction
                         # calculations, and 3D visualization
```

## Running the code!

Requires:
- numpy
- astropy
- matplotlib
- scipy
- pandas
- plotly

AND YOUR OWN GAIA DATA :) (queried via ADQL from the Gaia Archive)

## Status

Still a work in progress! Cluster selection parameters are being tuned to better match an earlier, more exploratory version of this analysis.

### Author: Eleanor Tchida
