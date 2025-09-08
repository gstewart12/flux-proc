# Delmarva Baywatch

A home for the Palmer Lab flux data processing pipeline

> [!NOTE] Limited applicability.
>
> I designed this pipeline to be flexible, but reproducibility and other best practices fell victim to triage as I prepared to defend my dissertation. Therefore, it likely has limited applicability beyond my data. My goal is to eventually make it useful to others and other datasets. This may happen soon, later, or never.

## File directory organization:

-   archive: anything outdated or superfluous but can't rule out future value
-   data: raw data

## Before starting pipeline:

-   Check .ghg files for potential errors
-   Process raw data in EddyPro
-   Download ERA5 data (R/get-era5.R)

## Pipeline:

1.  Prepare for post-processing (R/01-merge.R)

    -   *Run once for each site/year*
    -   Combine output from different EddyPro runs
    -   Check formatting
    -   Add external data
    -   Save a single file ready for post-processing

2.  (R/02-correct.R)

    -   *Run once for each site/year*
    -   Correct variable names & convert units
    -   Flux storage correction
    -   Adjust wind direction for magnetic declination
    -   Correct for zero-offset in incoming radiation
    -   Fix known site/year-specific issues
        -   These must already be entered into data/[SITE]/flagged_periods.csv
    -   Set primary versions of replicate measurements

3.  Automatic quality control (R/03-qc-auto.R)

    -   *Run once for each site/year*
    -   Clean flux dependencies
    -   Quality checks on biomet data:
        -   Dependencies
        -   Plausible limits
        -   Spikes
        -   Multivariate comparison
    -   Quality checks on flux data:
        -   Dependencies
        -   SSITC flags 
        -   Plausible time lag
        -   Spectral correction factor
        -   Sensor signal strength
        -   Precipitation
        -   Low frequency spikes
    -   Combine results into a single aggregated QC flag for each variable

4.  Gap-fill biomet variables (R/04-biomet-gapfill.R)

    -   *Run once for each site*
    -   De-bias alternate data (ERA5 & other tower sites)

5.  Gap-fill non-CH4 fluxes (R/05-mds-gapfill.R)

    -   *Run once for each site*
    -   
