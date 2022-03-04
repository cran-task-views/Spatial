---
name: Spatial
topic: Analysis of Spatial Data
maintainer: Roger Bivand, Jakub Nowosad
email: Roger.Bivand@nhh.no
version: 2022-01-24
source: https://github.com/cran-task-views/Spatial/
---


Base R includes many functions that can be used for reading,
visualising, and analysing spatial data. The focus in this view is on
"geographical" spatial data, where observations can be identified with
geographical locations, and where additional information about these
locations may be retrieved if the location is recorded with care.

Base R functions are complemented by contributed packages provided as 
source packages, and as ready-to-run binary packages for Windows and 
macOS (Intel 64-bit and Apple Silicon arm64 architectures). Information 
about source installs of packages using software external to R may be 
found at the end of this page. This task view covers the current status 
of contributed packages available from CRAN.

The contributed packages address two broad areas: moving spatial data
into and out of R including coordinate transformation, and analysing 
spatial data in R. Because the contributed packages constitute an evolving
ecosystem, there are several points of entry for users looking for help
and information. Two informal organisations curate websites: 
[r-spatial](https://r-spatial.org/) with a hyphen, and 
[rspatial](https://rspatial.org/) without. R-spatial is more generally 
geo-informatics based, grew from the legacy `r pkg("sp")` package and 
now clearly aligned with the modern `r pkg("sf", priority="core")` and 
`r pkg("stars", priority="core")` packages. Rspatial has grown from the 
`r pkg("raster")` package, now moving towards the modern 
`r pkg("terra", priority="core")` package. It is also worth noting the
wealth of online book projects, which may be helpful for users seeking
an introduction, including 
[Geocomputation with R](https://geocompr.robinlovelace.net/index.html).

Specific questions or issues may be raised where 
`packageDescription(<pkg>)$BugReports` returns an URL for bug 
reports or issues (where `<pkg>` is the
name of the package as a string), or directly with package
maintainers by email. Use may also be made of
the [R-SIG-Geo](https://stat.ethz.ch/mailman/listinfo/R-SIG-Geo/)
mailing-list after subscription, or of 
[stackoverflow](https://stackoverflow.com) with appropriate tags, or of
[stackexchange](https://gis.stackexchange.com/). Using the
`#rspatial` tag on [twitter](https://twitter.com) may also be worth trying,
or browsing traffic using that tag (among others).

The packages in this view can be roughly structured into the following
topics. If you think that some package is missing from the list, please
e-mail the maintainer or submit an issue or pull request in the GitHub
repository linked above.

Classes for spatial data and metadata
-------------------------------------

Many of the packages for handling and analysing spatial data use shared classes to
reduce duplication of effort. Up until 2016, the `r pkg("sp")` package provided
shared classes for spatial vector and raster data, but the representations used
preceded more modern and efficient international standards for spatial vector data.
From the release of `r pkg("sf", priority = "core")`, these modern vector
representations are to be preferred. For spatial raster data, the representations
proposed in `r pkg("stars", priority = "core")` and `r pkg("terra", priority = "core")`
suit overlapping but slightly different requirements. Conversion between objects
of classes defined by `r pkg("sf")`, `r pkg("stars")`, `r pkg("terra")` and the legacy
`r pkg("sp")` packages are available, and are described in [Conversions between different 
spatial classes in R](https://geocompr.github.io/post/2021/spatial-classes-conversion/).

Complementary initiatives are ongoing to support better handling of
geographic metadata in R.

### Spatial data - general

-   `r pkg("sf", priority = "core")` is a CRAN package for spatial vector data, 
    and is being actively developed here:
    `r github("r-spatial/sf")`, providing Simple Features
    for R, in compliance with the [OGC Simple
    Feature](http://www.opengeospatial.org/standards/sfa) standard. The
    development of the package was supported by the [R
    Consortium](https://www.r-consortium.org/). It provides simple
    features access for vector data, and as such is a modern
    implementation and standardization of parts of the legacy 
    `r pkg("sp")` package. `r pkg("sf")` is documented in an [R
    Journal](https://journal.R-project.org/archive/2018/RJ-2018-009/index.html)
    article. `r pkg("sf")` uses the PROJ, GEOS and GDAL external software
    libraries, which must be available for source installs together with
    other external software libraries that they in turn depend on.
-   `r pkg("stars", priority = "core")` is being actively developed here:
    `r github("rspatial/stars")`, and was supported by the [R
    Consortium](https://www.r-consortium.org/); it provides for
    spatiotemporal data in the form of dense arrays. It supercedes the 
    `r pkg("spacetime")` package, which 
    extended the shared classes defined in `r pkg("sp")` for
    spatio-temporal data (see [Spatio-Temporal Data in
    R](http://www.jstatsoft.org/v51/i07) ). `r pkg("stars")` uses PROJ and
    GDAL through `r pkg("sf")`.
-   The `r pkg("vapour")` package offers low-level access to GDAL functionality 
    for R packages. 
-   The `r pkg("spatstat", priority = "core")` contains classes suited to the
    analysis of point patterns, and may be coerced to and from `"sf"`, `"stars"`
    and other spatial classes.
-   The `r pkg("rcosmo")` package provides simple access to
    spherical and HEALPix data. It extends standard dataframes for
    HEALPix-type data.
-   `r pkg("inlmisc")` has followed on from Grid2Polygons
    and converts a spatial object from class SpatialGridDataFrame to
    SpatialPolygonsDataFrame among many other possibilities for legacy 
    `r pkg("sp")` classes.


### Raster data

-   `r pkg("terra", priority = "core")` is a re-implementation of
    `r pkg("raster")` functionality, linking directly to
    PROJ, GDAL and GEOS, and introducing new S4 classes for raster and
    vector data. See the [manual and
    tutorials](https://rspatial.org/terra/) to get started.
    `r pkg("terra")` is very similar to the `r pkg("raster")` package; but
    `r pkg("terra")` is simpler, better, and faster.
-   `r pkg("stars", priority = "core")` provides for spatiotemporal data in the
    form of dense arrays, with space and time being array dimensions.
    Examples include socio-economic or demographic data, environmental
    variables monitored at fixed stations, time series of satellite
    images with multiple spectral bands, spatial simulations, and
    climate model results.
-   The `r pkg("gdalcubes")` package also provides classes for data cubes,
    including proxy data cubes, it links to PROJ, GDAL and NetCDF4.

### Geographic metadata

-   `r pkg("geometa")` provides classes and methods to write
    geographic metadata following the ISO and OGC metadata standards
    (ISO 19115, 19110, 19119) and export it as XML (ISO 19139) for later
    publication into metadata catalogues. Reverserly, geometa provides a
    way to read ISO 19139 metadata into R. The package extends
    `r pkg("sf")` to provide GML (ISO 19136) representation
    of geometries. `r pkg("geometa")` is under active development on
    Github: `r github("eblondel/geometa")`.
-   `r pkg("ncdf4")` provides read and write functions for
    handling metadata (CF conventions) in the self-described NetCDF
    format.

Reading and writing spatial data
--------------------------------

**Reading and writing spatial data**

Spatial data is most often represented by one of two data models, 
vector or raster, and both models have many of their own file formats.
[GDAL (Geospatial Data Abstraction Library)](https://gdal.org/) is a (non-R)
library that provides a unified way to read and write hundreds of 
spatial data formats. Formats supported by GDAL include both OGC standard 
data formats (e.g., GeoPackage) and proprietary formats (e.g., ESRI Shapefile).
GDAL is used by a large number of GIS software and also many R packages,
such as `r pkg("sf")`, `r pkg("terra")`, and `r pkg("vapour")`.
This allows us to read and write spatial data in R from and to various 
spatial file formats. Important note: CRAN offers binary versions of packages 
`r pkg("sf")`, `r pkg("terra")`, and `r pkg("vapour")` for Windows and macOS,
that contain specific GDAL version with a subset of possible data source drivers.
If other drivers are needed, you need to either use other conversion utilities 
or install these packages from the source against a version of GDAL with the
required drivers.

In the past, `r pkg("rgdal")` and `r pkg("raster")` (through `r pkg("rgdal")`) were recommended for reading and writing of spatial data in R.
However, due to [the retirement of `r pkg("rgdal")` by the end of 2023](https://www.mail-archive.com/r-sig-geo@r-project.org/msg18468.html) new projects should not use it, and existing projects should implement migration to the packages mentioned in the previous paragraph.

### Reading and writing spatial data - data formats

Other packages provide facilities to read and write spatial data, 
dealing with open standard formats or proprietary formats.

*Open formats*

-   *Well-Known Text (WKT) / Well-Known Binary (WKB):* These standards are 
    part of the OGC Simple Feature specification. Both WKT/WKB formats are 
    supported by the `r pkg("sf")` package that implements the whole
    OGC Simple Feature specification in R. Additionally, `r pkg("wk")` and 
    `r pkg("wkutils")` may be used to parse well-known binary and
    well-known text representation of geometries to and from R-native formats.
-   *GeoJSON:* An rOpenSci [blog entry](http://ropensci.org/blog/blog/2016/11/22/geospatial-suite) 
    describes a GeoJSON-centred approach to reading GeoJSON and WKT data.
    The entry lists `r pkg("geojson")`, and `r pkg("geojsonio")`, among others.
    The GeoJSON format can also be read and write with `r pkg("sf")`,
    `r pkg("terra")`, and `r pkg("vapour")`. `r pkg("wellknown")` makes possible
    conversions from WKT to GeoJSON and GeoJSON to WKT.
-   *Geographic Markup Language (GML):* GML format can be read and written
    with `r pkg("sf")`. Additional GML native reader and writer is provided
    by `r pkg("geometa")` model with bindings to the `r pkg("sf")` classes, 
    for extension of geographic metadata with GML data and metadata 
    elements (GML 3.2.1 and 3.3) and interfacing OGC web-services 
    in `r pkg("ows4R")` package.
-   *NetCDF files:* NetCDF files can be read and write with 
    `r pkg("ncdf4")` or `r pkg("RNetCDF")`. Additionally, both `r pkg("terra")`
    and `r pkg("stars")` have capabilities for reading and writing NetCDF files.
-   *LAS / LAX:* These file formats are designed to work with lidar point
    cloud data and can be read/write with `r pkg("lidR")` or `r pkg("rLiDAR")`.

*Proprietary Data Formats*

-   *ESRI formats:* Many of spatial data saved into ESRI file formats
    can be read with GDAL, and thus also with `r pkg("sf")`, 
    `r pkg("terra")`, and `r pkg("vapour")`. Additionally, 
    `r pkg("shapefiles")` reads and writes ESRI ArcGIS/ArcView shapefiles.
    Additionally, `r pkg("maps")` (with `r pkg("mapdata")` and
    `r pkg("mapproj")`) provides a legacy tool to access to the same kinds
    of geographical databases as S.
-   *Others:* The `r pkg("gmt")` package gives a simple interface 
    between GMT map-making software and R.

### Reading and writing spatial data - GIS Software connectors

-   *PostGIS:* The `r pkg("rpostgis")` package provides additional functions
    to the `r pkg("RPostgreSQL")` package to interface R with a 
    'PostGIS'-enabled database, as well as convenient wrappers to 
    common 'PostgreSQL' queries. It is documented in an 
    [R Journal](https://journal.R-project.org/archive/2018/RJ-2018-025/index.html)
    article. `r pkg("postGIStools")` package provides functions to convert 
    geometry and 'hstore' data types from 'PostgreSQL' into standard R objects,
    as well as to simplify the import of R data frames 
    (including spatial data frames) into 'PostgreSQL'.
    `r pkg("sf")` also provides an R interface to PostGIS, 
    for both reading and writing, through GDAL.
<!--Roger, I assume we can keep this information about rgrass7 for now and update it when rgrass is on CRAN (I just discovered the name change using a GitHub search)-->
-   *GRASS GIS:* Integration with version 7.\* of the leading open source GIS,
    GRASS GIS, is provided in CRAN package `r pkg("rgrass7")`.
    For GRASS 6.\*, use `r pkg("spgrass6")`.
-   *SAGA GIS:* `r pkg("RSAGA")` and `r pkg("Rsagacmd")` offer shell-based
    wrapper for SAGA GIS commands.
-   *Quantum GIS (QGIS):* QGIS2 was supported by RQGIS 
    (`r github("r-spatial/RQGIS")`). QGIS3 (version >= 3.16) is supported by 
    `r github("paleolimbot/qgisprocess")`, which establishes an interface 
    between R and QGIS, i.e., it allows the user to access QGIS 
    functionalities from the R console. It achieves this by using 
    the qgis_process command-line utility.
-   *WhiteboxTools:* `r pkg("whitebox")` is an R frontend for 
    the WhiteboxTools software.
-   *ArcGIS:* `r pkg("RPyGeo")` is a wrapper for Python access 
    to the ArcGIS GeoProcessor. The ESRI company also offers their 
    own package (`r github("R-ArcGIS/r-bridge")`) that allows transferring data
    from ArcGIS to R.
-   Various GIS Software, including Orfeo ToolBox and SAGA GIS, can also be
    connected to R using `r pkg("link2GI")`.

### Interfaces to Spatial Web-Services

Some R packages focused on providing interfaces to web-services and web tools 
in support of spatial data management. Here follows a first tentative 
(non-exhaustive) list:

-   `r pkg("ows4R")` is a package that intends to provide an R interface 
    to OGC standard Web-Services. It is in active development at 
    `r github("eblondel/ows4R")` and currently support interfaces to the 
    Web Feature Service (WFS) for vector data access, with binding to the 
    `r pkg("sf")` package, and the Catalogue Service (CSW) for geographic
    metadata discovery and management (including transactions), with binding 
    to the `r pkg("geometa")` package.
-   `r pkg("geosapi")` is an R client for the [GeoServer](http://geoserver.org) 
    REST API, an open source implementation used widely for serving 
    spatial data.
-   `r pkg("geonapi")` provides an interface to the 
    [GeoNetwork](https://geonetwork-opensource.org/) legacy API, 
    an open source catalogue for managing geographic metadata. 
-   `r pkg("rgee")` is an [Earth Engine](https://earthengine.google.com/) 
    client library for R. All of the 'Earth Engine' API classes, modules, 
    and functions are made available. Additional functions implemented include 
    importing (exporting) of Earth Engine spatial objects, 
    extraction of time series, interactive map display, 
    assets management interface, and metadata display.

### Specific geospatial data sources of interest

-   `r pkg("rnaturalearth")` package facilitates interaction with 
    [Natural Earth](http://www.naturalearthdata.com/) map data. It includes 
    functions to download a wealth of Natural Earth vector and raster data,
    including cultural (e.g., country boundaries, airports, roads, railroads)
    and physical (e.g., coastline, lakes, glaciated areas) datasets.
-   Historical country boundaries (1886-today) can be obtained 
    from the `r pkg("cshapes")`.
-   `r pkg("marmap")` package is designed for downloading, plotting, 
    and manipulating bathymetric and topographic data in R. It allows to query 
    the ETOPO1 bathymetry and topography database hosted by the NOAA, 
    use simple latitude-longitude-depth data in ascii format, and take 
    advantage of the advanced plotting tools available in R to build 
    publication-quality bathymetric maps (see the
    [PLOS](http://www.plosone.org/article/info%3Adoi%2F10.1371%2Fjournal.pone.0073051)
    paper).
-   `r pkg("tidycensus")` provides access to US Census Bureau data in a 
    tidy format, including the option to bind the data spatially on import.
-   `r pkg("tigris")` provides access to cartographic elements provided by 
    the US Census Bureau TIGER, including cartographic boundaries, roads, 
    and water.
-   `r pkg("rgbif")` package is used to access Global Biodiversity 
    Information Facility (GBIF) occurrence data
-   `r pkg("geonames")` is an interface to the 
    [www.geonames.org](http://www.geonames.org/) service.
-   `r pkg("osmdata")` is an R package for accessing relatively 
    small datasets from OpenStreetMap (OSM), delivered via the Overpass API.
    `r pkg("osmextract")` matches, downloads, converts, and reads 
    OpenStreetMap data covering large areas, obtained from Geofabrik 
    and other providers. 
-   `r pkg("OpenStreetMap")` gives access to open street map raster images, 
    and `r pkg("osmar")` provides an infrastructure to access OpenStreetMap 
    data from different sources, to work with the data in a common R manner,
    and to convert data into available infrastructure provided by existing 
    R packages.
-   `r pkg("giscoR")` provides access to spatial elements provided by 
    GISCO - Eurostat, including boundary files of countries,  NUTS regions,
    municipalities, and other spatial objects.
-   `r pkg("chilemapas")` provides access to spatial data of political 
    and administrative divisions of Chile.
-   The former `geobr` package provided easy access to official spatial data sets of 
    Brazil for multiple geographies and years.
-   `r pkg("geouy")` loads and process geographic information for Uruguay.
-   `r pkg("rgugik")` allows to search and retrieve data from Polish Head 
    Office of Geodesy and Cartography ("GUGiK").
-   `r pkg("mapSpain")` downloads spatial boundary files of administrative 
    regions and other spatial objects of Spain.

Handling spatial data
---------------------

### Data processing - general

-   `r pkg("sf")` provides an interface to spatial geometry functions using
    the [GEOS](https://libgeos.org/) and [S2](https://s2geometry.io/) libraries.
-   `r pkg("stars")` contains tools for manipulating raster and vector
    data cubes.
-   `r pkg("terra")` package introduces many GIS methods for spatial vector
    and raster data.
-   The `r pkg("gdalUtils")` and `r pkg("gdalUtilities")` packages provide
    wrappers for the Geospatial Data Abstraction Library (GDAL) Utilities.
-   `r pkg("rmapshaper")` is a wrapper around the 'mapshaper' 'JavaScript'
    library to perform topologically-aware polygon simplification and other 
    operations such as clipping, erasing, dissolving, and converting 
    'multi-part' to 'single-part' geometries.
-   `r pkg("gdistance")`, provides functions to calculate distances and 
    routes on geographic grids.
    `r pkg("geosphere")` permits computations of distance and area to be 
    carried out on spatial data in geographical coordinates.
    `r pkg("cshapes")` package provides functions for calculating distance
    matrices (see [Mapping and Measuring Country Shapes](http://journal.R-project.org/archive/2010-1/RJournal_2010-1_Weidmann+Skrede~Gleditsch.pdf)).
-   `r pkg("magclass")` offers a data class for increased interoperability 
    working with spatial-temporal data together with corresponding functions
    and methods (conversions, basic calculations and basic data manipulation). 
-   The `r pkg("rcosmo")` package offers various tools for geometric
    transformations, computations, and statistical analysis of spherical data.
-   The `r pkg("trip")` package extends spatial classes to permit the
    accessing and manipulating of spatial data for animal tracking.
<!-- Roger, should we link here to papers such as https://besjournals.onlinelibrary.wiley.com/doi/10.1111/1365-2656.13116, https://link.springer.com/article/10.1007%2Fs40823-021-00067-y, https://link.springer.com/article/10.1007/s10109-020-00342-2, etc?-->

### Data cleaning

-   `r pkg("sf")` has a built-in functions st_is_valid to check whether 
    an sf geometry is valid and st_make_valid to fix invalid geometry (from GEOS 3.8).
-   `r pkg("lwgeom")` may also be used to facilitate handling and reporting
    of topology errors and geometry validity issues
    in sf objects.
    
### Data processing - specific

-   The `r pkg("landsat")` package with accompanying 
    [JSS paper](http://www.jstatsoft.org/v43/i04) provides tools for exploring
    and developing correction tools for remote sensing data.
-   The `r pkg("areal")` package can be used to interpolate overlapping 
    but incongruent polygons, also known as areal weighted interpolation.
-   The `r pkg("qualmap")` package can be used to digitize qualitative GIS data.

### Spatial sampling

-   `r pkg("spsurvey")` provides a range of sampling
    functions.
-   `r pkg("Spbsampling")` allows selecting probability samples well spread
    over the population of interest, in any dimension and using 
    any distance function.
-   `r pkg("spatialsample")` is a member of the tidymodel family of packages
    and contains functions and classes for spatial resampling to use with the
    `r pkg("rsample")`. 
-   `r pkg("MBHdesign")` provides spatially survey balanced designs using the
    quasi-random number method.
-   `r pkg("SpotSampling")` contains three methods for spatial and temporal
    sampling. 

Visualizing spatial data
------------------------

### Base visualization packages

-   Packages such as `r pkg("sf")`, `r pkg("stars")`, `r pkg("terra")`, 
    and `r pkg("rasterVis")` provide basic visualization methods through 
    the generic plot function.
-   `r pkg("classInt", priority = "core")` package provides functions for 
    choosing class intervals for thematic cartography.
-   `r pkg("rcosmo")` package provides several tools to interactively visualize
    HEALPix data, in particular, to plot data in arbitrary spherical windows.
-   Currently, the grDevices package (included with the R installation) 
    contains a large number of color palettes that can be accessed with 
    the hcl.colors and palette.colors functions; see also New features in this
    [blog](https://developer.r-project.org/Blog/public/2019/11/21/a-new-palette-for-r/index.html).
    Some of these color palettes can be also retrieved using separate packages,
    such as `r pkg("RColorBrewer")`, `r pkg("viridis")`, 
    or `r pkg("rcartocolor")`.

### Thematic cartography packages

-   `r pkg("tmap")` package accepts most spatial data classes and provides
    a modern basis for thematic mapping using a Grammar of Graphics syntax.
    It also allows for interactive spatial data mapping.
-   `r pkg("mapsf")` package allows various cartographic representations 
    such as proportional symbols, choropleth, or typology maps; it accepts sf
    (`r pkg("sf")`) and SpatRaster (`r pkg("terra")`) objects
-   `r pkg("ggplot2")` package has a built-in support for sf objects with the 
    geom_sf function and additional support for stars object is available
    through the geom_stars function available in the `r pkg("stars")` package.
    Its spatial visualization capabilities can be further extended with 
    `r pkg("ggspatial")`, which adds support for more spatial classes 
    (including classes from the raster package), allows adding north arrows
    and scale bars, etc.
-   The `r pkg("mapmisc")` package is a minimal, light-weight set of tools
    for producing nice-looking maps in R, with support for map projections.
-   Additional processing and mapping functions are available in 
    `r pkg("PBSmapping")` package; `r pkg("PBSmodelling")` provides 
    modelling support. In addition, `r pkg("GEOmap")` provides mapping
    facilities directed to meet the needs of geologists and uses 
    the `r pkg("geomapdata")` package.

### Packages based on web-mapping frameworks

-   `r pkg("mapview")` and `r pkg("leaflet")` packages provide methods to 
    view spatial objects interactively, usually on a web mapping base.
    Additionally, `r pkg("tmap")` has a view mode that allows for interactive 
    spatial data mapping.
-   `r pkg("mapdeck")` package provides a mechanism to plot interactive maps
    through javascript libraries 'Mapbox GL' and 'Deck.gl'.
-   `r pkg("RgoogleMaps")` package for accessing Google Maps(TM) may be
    useful if the user wishes to place a map backdrop behind other displays.
-   `r pkg("ggmap")` may be used for spatial visualization with Google Maps
    and OpenStreetMap; `r pkg("ggsn")` provides north arrows and scales 
    for such maps.
-   `r pkg("mapedit")` provides an R shiny widget based on `r pkg("leaflet")`
    for editing or creating sf geometries.

### Building Cartograms

-   `r pkg("cartogram")` package allows for constructions of a continuous 
    area cartogram by a rubber sheet distortion algorithm, non-contiguous
    area cartograms, and non-overlapping circles cartogram.
-   `r pkg("geogrid")` package turns polygons into rectangular or hexagonal
    cartograms.
-   `r pkg("micromap")` package provides linked micromaps using ggplot2.
-   `r pkg("recmap")` package provides rectangular cartograms with rectangle
    sizes reflecting for example population.
-   `r pkg("geogrid")` turns spatial polygons into regular or hexagonal grids. 
    `r pkg("statebins")` provides a simple binning approach to US states.

Analyzing spatial data
----------------------

The division of spatial statistics into three partly overlapping 
areas: point pattern analysis, geostatistics and the analysis of 
areal/lattice data, is widely accepted. However, areal data 
analysis can be split into disease mapping and spatial regression 
(also partly overlapping). In addition, ecological analyses often 
approach spatial data in particular ways, giving rise to a specific 
topical cluster of packages. All of these approaches to analysing spatial data
treat the spatial relationships between observations as a way of exploring
and making use of important sources of information about the observations over
and above what is known when assuming that they are independent of each other.

### Point pattern analysis

Point pattern analysis examines the distance relationships between observed points,
where the set of observations is expected to encompass all such entities in the 
study area.

-   `r pkg("spatstat", priority = "core")` is a family of R packages for analysing
    spatial point pattern data (and other kinds of spatial data). It has extensive
    capabilities for exploratory analysis, statistical modelling, simulation and
    statistical inference. It allows freedom in defining the region(s)
    of interest, and makes extensions to marked processes and spatial
    covariates. Its strengths are model-fitting and simulation, and it has a
    useful [homepage](http://www.spatstat.org/); it is 
    [actively developed](https://github.com/spatstat/spatstat). It is the only package
    that will enable the user to fit inhomogeneous point process models with
    interpoint interactions. 
-   The `r pkg("splancs")` package allows
    point data to be analysed within a polygonal region of interest, and
    covers many methods, including 2D kernel densities. 
-   The `r pkg("spatial")` package 
    is a recommended package shipped with base R, and contains several core functions,
    including an implementation of Khat by its author, Prof. Ripley. 
-   The `r pkg("spatgraphs")` package provides graphs, graph visualisation 
    and graph based summaries to be used with spatial point pattern analysis. 
-   The `r pkg("smacpod")` package provides various statistical
    methods for analyzing case-control point data. The methods available
    closely follow those in chapter 6 of Applied Spatial Statistics for
    Public Health Data by Waller and Gotway (2004).
-   `r pkg("ecespa")` provides wrappers, functions and data for
    spatial point pattern analysis, used in the book on Spatial Ecology of
    the ECESPA/AEET. The functions for binning points on grids in
-   `r pkg("ads")` may also be of interest. The ads package
    performs first- and second-order multi-scale analyses derived from
    Ripley's K-function. 
-   The `r pkg("dbmss")` package allows
    simple computation of a full set of spatial statistic functions of
    distance, including classical ones (Ripley's K and others) and more
    recent ones used by spatial economists (Duranton and Overman's Kd,
    Marcon and Puech's M). It relies on `r pkg("spatstat")` for core calculation.

### Geostatistics

Geostatistics uses a model fitted using the distances between observations to
interpolate values observed at point to unobserved points

-   The `r pkg("gstat", priority = "core")` package provides a
    wide range of functions for univariate and multivariate geostatistics,
    also for larger datasets.
-   `r pkg("geoR", priority = "core")` contains
    functions for model-based geostatistics. 
-   Variogram diagnostics may be
    carried out with `r pkg("vardiag")`. 
-   Automated interpolation using `r pkg("gstat")` is available 
    in `r pkg("automap")`. 
-   This family of packages is supplemented by `r pkg("intamap")` with 
    procedures for automated interpolation. 
-   A similar wide range of functions is to be found in the
    `r pkg("fields")` package, extended by `r pkg("LatticeKrig")` 
    for large spatial datasets and `r pkg("autoFRK")`.
-   The `r pkg("spatial")` package is shipped with base R, and
    contains several core geostatistical functions. 
-   The `r pkg("spBayes")` package fits Gaussian univariate and 
    multivariate models with MCMC.
-   `r pkg("ramps")` is a different Bayesian geostatistical
    modelling package. 
-   The `r pkg("geospt")` package contains
    some geostatistical and radial basis functions, including prediction and
    cross validation. Besides, it includes functions for the design of
    optimal spatial sampling networks based on geostatistical modelling. 
-   The `r pkg("rcosmo")` package offers various geostatistics
    methods for spherical data: descriptive statistics, entropy based
    methods, covariance-variogram methods, etc. Most of rcosmo features were
    developed for Cosmic Microwave Background data, but they can also be
    used for any spherical data. 
-   The `r pkg("FRK")` package is a tool for 
    spatial/spatio-temporal modelling and prediction with large
    datasets. The approach, discussed in Cressie and Johannesson (2008),
    decomposes the field, and hence the covariance function, using a fixed
    set of n basis functions, where n is typically much smaller than the
    number of data points (or polygons) m.
-   The `r pkg("RandomFields", priority = "core")` package
    provides functions for the simulation and analysis of random fields, and
    variogram model descriptions can be passed between
    `r pkg("geoR")`, `r pkg("gstat")` and this package. 
-   `r pkg("SpatialExtremes")` proposes several
    approaches for spatial extremes modelling using `r pkg("RandomFields")`. 
-   In addition, `r pkg("CompRandFld")`, `r pkg("constrainedKriging")` and
    `r pkg("geospt")` provide alternative approaches to geostatistical modelling. 
-   The `r pkg("spTimer")` package is able to fit, spatially predict and 
    temporally forecast large amounts of space-time data using \[1\] Bayesian 
    Gaussian Process (GP) Models, \[2\] Bayesian Auto-Regressive (AR) Models, 
    and \[3\] Bayesian Gaussian Predictive Processes (GPP) based AR Models. 
-   The `r pkg("rtop")` package provides functions for the
    geostatistical interpolation of data with irregular spatial support such
    as runoff related data or data from administrative units. 
-   The `r pkg("georob")` package provides functions for fitting
    linear models with spatially correlated errors by robust and Gaussian
    Restricted Maximum Likelihood and for computing robust and customary
    point and block kriging predictions, along with utility functions for
    cross-validation and for unbiased back-transformation of kriging
    predictions of log-transformed data. 
-   The `r pkg("SpatialTools")` package has an emphasis on kriging,
    and provides functions for prediction and simulation. It is extended by
    `r pkg("ExceedanceTools")`, which provides tools for
    constructing confidence regions for exceedance regions and contour
    lines. 
-   The `r pkg("gear")` package implements common
    geostatistical methods in a clean, straightforward, efficient manner,
    and is said to be a quasi reboot of `r pkg("SpatialTools")`.
-   The `r pkg("sperrorest")` package implements spatial error
    estimation and permutation-based spatial variable importance using
    different spatial cross-validation and spatial block bootstrap methods, used by
    `r pkg("mlr3spatiotempcv")`.
-   The `r pkg("sgeostat")` package is also available. Within
    the same general topical area are the
    `r pkg("deldir", priority = "core")` package for triangulation and the
    `r pkg("interp")` package for spline interpolation; the
    `r pkg("MBA")` package provides scattered data interpolation
    with multilevel B-splines. 
-   In addition, there are the
    `r pkg("spatialCovariance")` package, which supports the
    computation of spatial covariance matrices for data on rectangles, the
    `r pkg("regress")` package building in part on
    `r pkg("spatialCovariance")`, and the
    `r pkg("tgp")` package. 
-   The `r pkg("Stem")`
    package provides for the estimation of the parameters of a
    spatio-temporal model using the EM algorithm, and the estimation of the
    parameter standard errors using a spatio-temporal parametric bootstrap.
-   `r pkg("FieldSim")` is another random fields simulations
    package. 
-   The `r pkg("SSN")` is for geostatistical modeling
    for data on stream networks, including models based on in-stream
    distance. Models are created using moving average constructions. Spatial
    linear models, including covariates, can be fit with ML or REML. Mapping
    and other graphical functions are included. 
-   The `r pkg("ipdw")` provides functions to interpolate
    georeferenced point data via Inverse Path Distance Weighting. Useful
    for coastal marine applications where barriers in the landscape
    preclude interpolation with Euclidean distances.
-   `r pkg("RSurvey")` may be used as a processing program for
    spatially distributed data, and is capable of error corrections and data
    visualisation.

### Disease mapping and areal data analysis

Both point pattern analysis and geostatistics enter into disease mapping, which
is concerned with representing public health information over space and
time in a communicative and responsible way. Estimation is important to present
calculated rates that are comparable both in terms of levels and uncertainty.

-   `r pkg("DCluster", priority = "core")` is a package for the
    detection of spatial clusters of diseases. It is complemented by
    `r pkg("DClusterm")` for model-based cluster detection, and by 
    `r pkg("rflexscan")` and `r pkg("FlexScan")`, two implementations
    of flexible scan statistics.
-   `r pkg("DCluster")` extends and depends on the
    `r pkg("spdep", priority = "core")` package, which provides
    basic functions for building neighbour lists and spatial weights. 
-   `r pkg("spdep")` also provides global and local tests for spatial
    autocorrelation, including join-count tests, Moran's I, Geary's C,
    Getis-Ord G and others.
-   `r pkg("rgeoda")` is a wrapper for GeoDa and provides efficient alternatives
    for calculating global and local tests for spatial autocorrelation.
-   Some functions for fitting spatial regression models, such as SAR and CAR 
    models are in `r pkg("spatialreg", priority = "core")`, see below.
-   The `r pkg("SpatialEpi")` package provides implementations of
    cluster detection and disease mapping functions, including Bayesian
    cluster detection, and supports strata. 
-   The `r pkg("smerc")` package provides statistical methods 
    for the analysis of data areal data, with a focus on cluster detection. 
-   The `r pkg("diseasemapping")` package offers the formatting of
    population and case data, calculation of Standardized Incidence Ratios,
    and fitting the BYM model using INLA. 
-   A Markov Random Field `"mrf"` effect may be added to models in the 
    `r pkg("mgcv")` package shipped with base R, providing flexible
    modelling tools in a recommended package.
-   The `r pkg("hglm")` package also provides SAR and CAR model fitting
    approaches.
-   Regionalization of polygon objects
    is provided by `r pkg("AMOEBA")`: a function to calculate
    spatial clusters using the Getis-Ord local statistic. It searches
    for irregular clusters (ecotopes) on a map, as does `skater()` in
    `r pkg("spdep")`. 
-   The `r pkg("seg")`, `r pkg("divseg")` and
    `r pkg("OasisR")` packages provide functions for measuring
    spatial segregation; `r pkg("OasisR")` includes Monte Carlo
    simulations to test the indices. 
-   The `r pkg("lctools")` package provides
    researchers and educators with easy-to-learn user friendly tools for
    calculating key spatial statistics and to apply simple as well as
    advanced methods of spatial analysis in real data. These include: Local
    Pearson and Geographically Weighted Pearson Correlation Coefficients,
    Spatial Inequality Measures (Gini, Spatial Gini, LQ, Focal LQ), Spatial
    Autocorrelation (Global and Local Moran's I), several Geographically
    Weighted Regression techniques and other Spatial Analysis tools (other
    geographically weighted statistics). This package also contains
    functions for measuring the significance of each statistic calculated,
    mainly based on Monte Carlo simulations. The
-   `r pkg("sparr")` package provides another approach to
    relative risks. 
-   The `r pkg("CARBayes")` package implements
    Bayesian hierarchical spatial areal unit models. In such models, the
    spatial correlation is modelled by a set of random effects, which are
    assigned a conditional autoregressive (CAR) prior distribution. Examples
    of the models included are the BYM model as well as a recently developed
    localised spatial smoothing model. 
-   The `r pkg("spaMM")`
    package fits spatial GLMMs, using the Matern correlation function as the
    basic model for spatial random effects. 
-   The `r pkg("PReMiuM")` package is for profile regression, which
    is a Dirichlet process Bayesian clustering model; it provides a spatial
    CAR term that can be included in the fixed effects (which are global,
    ie. non-cluster specific, parameters) to account for any spatial
    correlation in the residuals. 
-   Spatial survival analysis is provided by the
    `r pkg("spBayesSurv")` package: Bayesian Modeling and
    Analysis of Spatially Correlated Survival Data. 
-   The `r pkg("spselect")` package provides modelling functions
    based on forward stepwise regression, incremental forward stagewise
    regression, least angle regression (LARS), and lasso models for
    selecting the spatial scale of covariates in regression models.
-   Spatial microsimulation is offered by `r pkg("rakeR")`, `r pkg("sms")`, 
    `r pkg("synthACS")`, and `r pkg("NetLogoR")` permits the building and
    running of spatially explicit agent-based models.


### Spatial regression

Many packages providing functions for fitting spatial regression models have 
already been given as they are used in disease mapping. In this subsection, 
more attention is given to the subset of methods used in spatial econometrics, and
so complements general econometric methods covered in the `r view("Econometrics")` 
Task View.

-   The choice of function for spatial regression will depend on the support
    available. If the data are characterised by point support and the
    spatial process is continuous, geostatistical methods may be used, or
    functions in the `r pkg("nlme")` package. 
-   If the support is
    areal, and the spatial process is not being treated as continuous,
    functions provided in the `r pkg("spatialreg")` package may
    be used. This package can also be seen as providing spatial econometrics
    functions. `r pkg("spdep")` provides the full range
    of local indicators of spatial association, such as local Moran's I and
    diagnostic tools for fitted linear models, including Lagrange Multiplier
    tests. Spatial regression models that can be fitted using maximum
    likelihood and Bayesian MCMC methods in
    `r pkg("spatialreg")` include spatial lag models, spatial
    error models, two parameter models, their Durbin variants and SLX models.
    For larger data sets, sparse
    matrix techniques can be used for maximum likelihood fits. In
    `r pkg("spatialreg")`, the `ME` and `SpatialFiltering`
    functions provide Moran Eigenvector model fitting, as do more modern
    functions in the `r pkg("spmoran")` package. 
-   When using the generalized method of moments (GMM), `r pkg("sphet")` can be used
    to accommodate both autocorrelation and heteroskedasticity, also with 
    instrumental variables. 
-   The `r pkg("splm")` package provides methods for fitting spatial
    panel data by maximum likelihood and GM. 
-   The `r pkg("spsur")` package provides functions to test and estimate 
    spatial seemingly unrelated regression models (spatial SUR) by maximum 
    likelihood and three-stage least squares. 
-   The two small packages `r pkg("S2sls")` and `r pkg("spanel")` provide
    alternative implementations without most of the facilities of
    `r pkg("splm")`. 
-   The former `HSAR` package provides Hierarchical
    Spatial Autoregressive Models (HSAR), based on a Bayesian Markov Chain
    Monte Carlo (MCMC) algorithm. 
-   `r pkg("spatialprobit")` makes
    possible Bayesian estimation of the spatial autoregressive probit model
    (SAR probit model). 
-   The `r pkg("ProbitSpatial")` package provides methods for
    fitting Binomial spatial probit models to larger data sets; spatial
    autoregressive (SAR) and spatial error (SEM) probit models are included.
-   The `r pkg("starma")` package provides functions to
    identify, estimate and diagnose a Space-Time AutoRegressive Moving
    Average (STARMA) model.
-   The `r pkg("spgwr")`
    package contains an implementation of geographically weighted regression
    methods for exploring possible spatial non-stationarity. The
    `r pkg("gwrr")` package fits geographically weighted
    regression (GWR) models and has tools to diagnose and remediate
    collinearity in the GWR models. It also fits geographically weighted ridge
    regression (GWRR) and geographically weighted lasso (GWL) models. The
    `r pkg("GWmodel")` package contains functions for
    computing geographically weighted (GW) models. Specifically, basic, 
    robust, local ridge, heteroskedastic, mixed, multiscale, generalised 
    and space-time GWR; GW summary statistics, GW PCA and GW discriminant
    analysis; associated tests and diagnostics; and options for a range of
    distance metrics. 

### Ecological analysis

There are many packages for analysing ecological and environmental data.
They include:

-   `r pkg("ade4")` for exploratory and Euclidean methods in
    the environmental sciences, the adehabitat family of packages for
    the analysis of habitat selection by animals
    (`r pkg("adehabitatHR")`,
    `r pkg("adehabitatHS")`,
    `r pkg("adehabitatLT")`, and
    `r pkg("adehabitatMA")`)
-   `r pkg("pastecs")` for the regulation, decomposition and
    analysis of space-time series
-   `r pkg("vegan")` for ordination methods and other useful
    functions for community and vegetation ecologists, and many other
    functions in other contributed packages. One such is
    `r pkg("tripEstimation")`, basing on the classes
    provided by `r pkg("trip")`. 
-   `r pkg("ncf")` provides a range of spatial
    nonparametric covariance functions.
-   The `r pkg("spind")` package provides functions for
    spatial methods based on generalized estimating equations (GEE) and
    wavelet-revised methods (WRM), functions for scaling by wavelet
    multiresolution regression (WMRR), conducting multi-model inference,
    and stepwise model selection.
-   The `r pkg("siplab")` package is a platform for
    experimenting with spatially explicit individual-based vegetation
    models.
-   `r pkg("ModelMap")` builds on other packages to create
    models using underlying GIS data.
-   The `r pkg("SpatialPosition")` computes spatial position
    models: Stewart potentials, Reilly catchment areas, Huff catchment
    areas.
-   The `r pkg("Watersheds")` package provides methods for
    watersheds aggregation and spatial drainage network analysis.
-   The `r pkg("ngspatial")` package provides tools for
    analyzing spatial data, especially non-Gaussian areal data. It
    supports the sparse spatial generalized linear mixed model of Hughes
    and Haran (2013) and the centered autologistic model of Caragea and
    Kaiser (2009).
-   `r pkg("landscapemetrics")` package calculates landscape
    metrics for categorical landscape patterns. It can be used as a
    drop-in replacement for
    [FRAGSTATS](https://www.umass.edu/landeco/research/fragstats/fragstats.html)
    , as it offers a reproducible workflow for landscape analysis in a
    single environment. It also provides several visualization
    functions, e.g. to show all labeled patches or the core area of all
    patches.

The `r view("Environmetrics")` Task View contains a much more
complete survey of relevant functions and packages.


Installing packages linking to PROJ, GDAL or GEOS
-------------------------------------------------

Installation of packages like `r pkg("sf")` and `r pkg("terra")` which use
external software libraries such as PROJ, GDAL or GEOS requires care. For
most users on platforms such as Windows or macOS who are not themselves 
package developers, it is always better to avoid what are known as 
source installs, because CRAN binary packages include all of the external 
software required. Because `getOption("pkgType")` on these platforms is
usually `"both"`, you may be asked to choose to install a source package
if it is more recent than the latest binary. 

Please do not be tempted to choose a source install for `r pkg("sf")` or 
`r pkg("terra")` or similar; the binary package will be generated within a
day or two. To avoid being asked, you may see from `?options` 
under options provided by the utils package that the default 
behaviour of your installation of R may be controlled
by setting options `install.packages.check.source` and
`install.packages.compile.from.source`, or by setting environment variable
`R_COMPILE_AND_INSTALL_PACKAGES`, see also this [helpful comment](https://github.com/r-spatial/sf/issues/1848#issuecomment-1005038753).

If you are a developer using Windows or macOS or installing 
from `github`, the same static-linked 
binary external software libraries, header files, etc. as those used in 
building CRAN binary packages are available from: Windows 4.0 and 4.1 
[downloaded on-the-fly](https://github.com/orgs/rwinlib/repositories), 
Windows 4.2 [forthcoming rtools42](https://www.r-project.org/nosvn/winutf8/ucrt3/) 
and macOS [both architectures](https://mac.r-project.org/tools/). These 
external software libraries have been built using the same compile and 
link settings as R itself, so avoid the risk of possible errors caused by 
mismatched binaries. 

If you are a user (or developer) on systems where `getOption("pkgType")`
is `"source"`, you will need to ensure that the external software is
available when installing source packages. Advice for some such systems may
be found [here](https://github.com/r-spatial/sf/#installing). 
[The most common reason](https://github.com/r-spatial/sf/#multiple-gdal-geos-andor-proj-versions-on-your-system) 
for failure is having multiple versions of external software installed
on your platform.


### Links
-   [R-SIG-Geo mailing list](https://stat.ethz.ch/mailman/listinfo/R-SIG-Geo/)
