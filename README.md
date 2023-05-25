# ForecastingAndMonitoringOptimization

This is the repository to the article "The dependence of forecasts on sampling frequency as a guide to optimizing monitoring in community ecology".


DOI to be added

## Repository structure

```
ForecastingAndMonitoringOptimization
|-- README.md
|-- LICENSE.txt   
|-- ForecastingAndMonitoringOptimization.Rproj
|-- Data (to be added)
    |-- CCM_analysis_data
    |-- forecast_data
    |-- SimulationsForecastVsGrowthRate
    |-- SMap_interaction_data
    |-- Timeseries
|-- R
    |-- Functions
    |-- DataAnalysis
```

### Files and directories

- **README.md**: this file
- **LICENSE.txt**: document containing the license conditions
- **ForecastingAndMonitoringOptimization.Rproj**: R project file
- **Data**: directory containing all the data used in the study
  - **Timeseries**: directory containing the actual data
  - The other directories contain the results of the various analyses. They are **not** needed to reproduce the study, but considering that some analyses can take days or weeks, we provide their results as .RData or .csv files as well.
- **R**: directory containg all the code needed to reproduce the study
  - **Functions**: directory containing several R scripts in which several functions are implemented
  - **DataAnalysis**: directory containing the R script used in the study, i.e. to analyse the data.
