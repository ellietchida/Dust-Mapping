# Dust Mapping Project

## Introduction

This project uses Gaia DR3 photometry to study interstellar dust extinction using the open cluster Cl Tombaugh 4 and the red clump.

While looking at Tombaugh 4's color-magnitude diagram, I noticed a big color shift and dimming compared to what a normal population should look like on the CMD, which is a sign of heavy dust extinction. That observation is what inspired this project and the steps that followed: isolate the cluster, find red clump stars (standard candles), use them to measure extinction and reddening as a function of distance.

## How It Works/Key Components

1. Cluster selection: Sky position, proper motion cuts (sigma-clipped), well-measured radial velocity/parallax masks, and optional additional masking isolate likely cluster members (in our case, Tombaugh 4).
2. Galactic component classification: Stars are transformed into Galactocentric coordinates so we can check their height above the Galactic plane and classify them as thin disk, thick disk, or halo.
3. Red clump identification: A color-magnitude box (plus the disk label) picks out red clump stars, which have a nearly universal absolute magnitude, making them good and reliable standard candles.
4. Extinction & reddening: Using the red clump's known intrinsic color/magnitude, we compute color excess and extinction.
5. 3D dust mapping: Combining extinction with 3D positions gives spatial maps of dust.

## Important Takeaways

- The original version of this analysis went star-by-star in a messy, repetitive way (copy-pasted selection blocks for every round). It's been reworked into 'cluster_selection.py', using classes ('ClusterSelector' and 'SelectionRound') so the same selection logic gets reused while keeeping rounds isolated for performance evalutation (if desired).
- Sigma-clipping by proper motion only works for populations that actually share a common motion (like a real cluster), so be careful with this selection feature and on which round.
- Effectiveness of selection rounds is a work-in-progress. Given the nature of the original design, this is a very simple project that will be expanded upon.

## Outputs

The code produces:
- Color-magnitude diagrams (cluster rounds, full background survey, red clump)
- Thin disk / thick disk / halo classification
- Many visualization plots of position, extinction, etc.
- 2D and 3D dust density/extinction maps (interactive!)

## Repository Contents

- cluster_selection.py   --> Reusable ClusterSelector / SelectionRound classes that create cluster subsets based on selection parameters. Also produces galactic component classifications, absolute magnitude calculations. 
- Dust_Mapping.ipynb   -->   Main analysis notebook, basically a playground for extinction study. Runs selection, CMD plotting, extinction calculations, and 3D visualization of dust/extinction.

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

Still a work in progress! Cluster selection parameters are being tuned, as is the case with performance evaluations. Can always improve!

### Author: Eleanor Tchida
