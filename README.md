# Getting Started with Cloud-Native Harmonized Landsat Sentinel-2 (HLS) Data in R

>**NOTE:** This tutorial relies on **CMR-STAC**, which will continue to evolve in order to align with the [STAC specification](https://github.com/radiantearth/stac-spec) over time. As these changes are released, the LP DAAC will update this tutorial as necessary. If you are encountering an issue with the tutorial, please review the [CHANGELOG.md](https://git.earthdata.nasa.gov/projects/LPDUR/repos/hls_tutorial_r/browse/CHANGELOG.md?at=main) and make sure that you have copied/cloned/downloaded the latest release of this tutorial. If problems persist, please use the contact information below to submit your issue to the LP DAAC User Services.  

## Objective  

NASA's Land Processes Distributed Active Archive Center (LP DAAC) archives and distributes HLS products in the LP DAAC Cumulus cloud archive as Cloud Optimized GeoTIFFs (COG). The primary objective of this tutorial is to show how to perform spatiotemporal queries using NASA’s Common Metadata Repository (CMR) SpatioTemporal Asset Catalog (STAC) to find and subset HLS data.  

## Tutorial  

- [HLS_Tutorial](https://git.earthdata.nasa.gov/projects/LPDUR/repos/hls_tutorial_r/browse/HLS_Tutorial.html?at=refs%2Fheads%2Fmain) shows how to query and subset HLS data using NASA’s CMR-STAC API. It showns how to stack the scenes, process data (quality filtering and NDVI calculation), visualize, calculate statistics for an NDVI time series, and export the outputs.  

## Use Case Description  

This tutorial demonstrates an example use case for crop monitoring over a multiple farm fields in Central Valley California in the United States. The goal of the example is to observe a time series of HLS-derived NDVI over agricultural fields in northern California without downloading the entirety of the HLS source data.  

### Products Used

1. Daily 30 meter (m) global HLS Sentinel-2 Multi-spectral Instrument Surface Reflectance - [HLSS30.002](https://doi.org/10.5067/HLS/HLSS30.002)  

    - **Science Dataset (SDS) layers:**  
      - B8A (NIR Narrow)  
      - B04 (Red)  
      - Fmask (Quality)  

2. Daily 30 meter (m) global HLS Landsat-8 OLI Surface Reflectance - [HLSL30.002](https://doi.org/10.5067/HLS/HLSL30.002)  

    - **Science Dataset (SDS) layers:**  
      - B05 (NIR)  
      - B04 (Red)  
      - Fmask (Quality)  

## Setup  

### Prerequisites  

- A [NASA Earthdata Login](https://urs.earthdata.nasa.gov/) account is required to access the data used in the `HLS_Tutorial.Rmd` Tutorial. You can create an account [here](https://urs.earthdata.nasa.gov/users/new).  
- Install R and RStudio. Details can be found [here](https://www.rstudio.com/products/rstudio/download/#download). These tutorials have been tested on Windows using R Version 4.1.0 and RStudio version 1.4.1717.  

### Get started  

1. [Clone](ssh://git@git.earthdata.nasa.gov:7999/lpdur/hls_tutorial_r.git) or [download](https://git.earthdata.nasa.gov/rest/api/latest/projects/LPDUR/repos/hls_tutorial_r/archive?format=zip) HLS_Tutorial_R repository from the LP DAAC Data User Resources Repository.  
2. Open the `HLS_Tutorial.Rproj` file to directly open the project. Next, select the R Markdown files (CMR_STAC_Tutorial.Rmd and HLS_Tutorial.Rmd) from the `Scripts` folder and open them.  

---

## Contact Information  

**Authors:** Mahsa Jami¹ and Aaron Friesz¹  
**Contact:** LPDAAC@usgs.gov  
**Voice:** +1-866-573-3222  
**Organization:** Land Processes Distributed Active Archive Center (LP DAAC)  
**Website:** <https://lpdaac.usgs.gov/>  
**Date last modified:** 09-09-2021  

¹KBR, Inc., contractor to the U.S. Geological Survey, Earth Resources Observation and Science (EROS) Center,  
 Sioux Falls, South Dakota, USA. Work performed under USGS contract G15PD00467 for LP DAAC².  
²LP DAAC Work performed under NASA contract NNG14HH33I.  
