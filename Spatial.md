---
name: Spatial
topic: Analysis of Spatial Data
maintainer: Roger Bivand
email: Roger.Bivand@nhh.no
version: 2021-12-29
source: https://github.com/cran-task-views/Spatial/
---


Base R includes many functions that can be used for reading,
visualising, and analysing spatial data. The focus in this view is on
"geographical" spatial data, where observations can be identified with
geographical locations, and where additional information about these
locations may be retrieved if the location is recorded with care.

Base R functions are complemented by contributed packages, some of which
are on CRAN, and others are still in development. One location is
[Github](https://github.com/) . Some key packages including
`r pkg("sf", priority = "core")` and
`r pkg("stars")` are grouped under
[r-spatial](https://github.com/r-spatial) , others including
`r pkg("raster", priority = "core")` and
`r pkg("terra")` under
[rspatial](https://github.com/rspatial) . Maintenance of the
`r pkg("sp", priority = "core")` is continuing here:
`r github("edzer/sp")`.

Another set of locations for the development and maintenance of packages
on [R-Forge](https://R-Forge.R-project.org/) , which lists "Spatial
Data and Statistics" projects in its [project
tree](https://R-Forge.R-project.org/softwaremap/trove_list.php) .
Information on R-spatial packages was until 2016 posted on the R-Forge
rspatial project [website](https://rspatial.R-Forge.R-project.org/) ,
including a visualisation gallery.

The contributed packages address two broad areas: moving spatial data
into and out of R, and analysing spatial data in R.

The [R-SIG-Geo](https://stat.ethz.ch/mailman/listinfo/R-SIG-Geo/)
mailing-list is a good place to begin for obtaining help and discussing
questions about both accessing data, and analysing it. The mailing list
is a good place to search for information about relevant courses.
Further information about courses may be found under the "Events" tab
of [this blog](http://r-spatial.org/) .

There are a number of contributed tutorials and introductions; a recent
one is [Introduction to visualising spatial data in
R](https://CRAN.R-project.org/doc/contrib/intro-spatial-rl.pdf) by Robin
Lovelace and James Cheshire.

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
suit overlapping but slightly different requirements.

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

**Reading and writing spatial data -
`r pkg("rgdal", priority = "core")`**

Maps may be vector-based or raster-based. The
`r pkg("rgdal")` package provides bindings to [GDAL
(Geospatial Data Abstraction Library)](http://www.gdal.org/) -supported
raster formats and [OGR](http://www.gdal.org/ogr/) -supported vector
formats. It contains functions to write raster and vector files in
supported formats. Formats supported by GDAL/OGR include both OGC
standard data formats (e.g. GeoJSON) and proprietary formats (e.g. ESRI
Shapefile). The package also provides [PROJ.4](https://proj4.org/)
projection support for vector objects ( [this
site](http://spatialreference.org) provides searchable online PROJ.4
representations of projections). Affine and similarity transformations
on sp objects may be made using functions in the
`r pkg("vec2dtransf")` package. The Windows and Mac OSX CRAN
binaries of `r pkg("rgdal")` include subsets of possible
data source drivers; if others are needed, use other conversion
utilities, or install from source against a version of GDAL with the
required drivers.

### Reading and writing spatial data - data formats

Other packages provide facilities to read and write spatial data,
dealing with open standard formats or proprietary formats.

*OGC Standard Data formats*

-   *Well-Known Text (WKT) / Well-Known Binary (WKB):* These standards
    are part of the OGC Simple Feature specification. Both WKT/WKB
    formats are supported by `r pkg("sf")` package that
    implements the whole OGC Simple Feature specification in R. Apart
    from the `r pkg("sf")` package, the
    `r pkg("rgeos", priority = "core")` package provides
    functions for reading and writing well-known text (WKT) geometry.
    Package `r pkg("wkb")` package provides functions for
    reading and writing well-known binary (WKB) geometry.
-   *GeoJSON:* An rOpenSci [blog
    entry](http://ropensci.org/blog/blog/2016/11/22/geospatial-suite)
    described a GeoJSON-centred approach to reading GeoJSON and WKT
    data. GeoJSON can be written and read using
    `r pkg("rgdal")`, and WKT by
    `r pkg("rgeos")`. The entry lists
    `r pkg("geojson")`, and
    `r pkg("geojsonio")`, among others.
-   *Geographic Markup Language (GML):* GML format can be read and
    writen with `r pkg("rgdal")`. Additional GML native
    reader and writer is provided by `r pkg("geometa")`
    model with bindings to the `r pkg("sf")` classes, for
    extension of geographic metadata with GML data and metadata elements
    (GML 3.2.1 and 3.3) and interfacing OGC web-services in
    `r pkg("ows4R")` package
-   *NetCDF files:* `r pkg("ncdf4")` or
    `r pkg("RNetCDF")` may be used.

*Proprietary Data Formats*

-   *ESRI formats:* `r pkg("maps")` (with
    `r pkg("mapdata")` and `r pkg("mapproj")`)
    provides access to the same kinds of geographical databases as S.
    `r pkg("maptools")` and
    `r pkg("shapefiles")` read and write ESRI ArcGIS/ArcView
    shapefiles.
-   *Others:* `r pkg("maptools")` package provides helper
    functions for writing map polygon files to be read by *WinBUGS* ,
    *Mondrian* , and the tmap command in *Stata* . The
    `r pkg("gmt")` package gives a simple interface between
    GMT map-making software and R.

### Reading and writing spatial data - GIS Software connectors

-   *PostGIS:* The `r pkg("rpostgis")` package provides
    additional functions to the `r pkg("RPostgreSQL")`
    package to interface R with a 'PostGIS'-enabled database, as well
    as convenient wrappers to common 'PostgreSQL' queries. It is
    documented in an [R
    Journal](https://journal.R-project.org/archive/2018/RJ-2018-025/index.html)
    article. `r pkg("postGIStools")` package provides
    functions to convert geometry and 'hstore' data types from
    'PostgreSQL' into standard R objects, as well as to simplify the
    import of R data frames (including spatial data frames) into
    'PostgreSQL'. `r pkg("sf")` also provides an R
    interface to Postgis, for both reading and writing, throuh GDAL.
-   *GRASS:* Integration with version 7.\* of the leading open source
    GIS, GRASS, is provided in CRAN package
    `r pkg("rgrass7")`, using `r pkg("rgdal")`
    for exchanging data. For GRASS 6.\*, use
    `r pkg("spgrass6")`.
-   *SAGA:* `r pkg("RSAGA")` is a similar shell-based
    wrapper for SAGA commands.
-   *Quantum GIS (QGIS):* QGIS2 was supported by RQGIS. QGIS3 is
    supported by `r github("r-spatial/RQGIS3")`, which
    establishes an interface between R and QGIS, i.e. it allows the user
    to access QGIS functionalities from the R console. It achieves this
    by using the QGIS Python API.
-   *ArcGIS:* `r pkg("RPyGeo")` is a wrapper for Python
    access to the ArcGIS GeoProcessor

### Interfaces to Spatial Web-Services

Some R packages focused on providing interfaces to web-services and web
tools in support of spatial data management. Here follows a first
tentative (non-exhaustive) list:

-   `r pkg("ows4R")` is a new package that intends to
    provide an R interface to OGC standard Web-Services. It is in active
    development at `r github("eblondel/ows4R")` and
    currently support interfaces to the Web Feature Service (WFS) for
    vector data access, with binding to the `r pkg("sf")`
    package, and the Catalogue Service (CSW) for geographic metadata
    discovery and management (including transactions), with binding to
    the `r pkg("geometa")` package.
-   `r pkg("geosapi")` is an R client for the
    [GeoServer](http://geoserver.org) REST API, an open source
    implementation used widely for serving spatial data.
-   `r pkg("geonapi")` provides an interface to the
    [GeoNetwork](https://geonetwork-opensource.org/) legacy API, an
    opensource catalogue for managing geographic metadata.
-   `r pkg("rgee")` ia an [Earth
    Engine](https://earthengine.google.com/) client library for R. All
    of the 'Earth Engine' API classes, modules, and functions are made
    available. Additional functions implemented include importing
    (exporting) of Earth Engine spatial objects, extraction of time
    series, interactive map display, assets management interface, and
    metadata display.

### Specific geospatial data sources of interest

-   `r pkg("rnaturalearth")` package facilitates interaction
    with [Natural Earth](http://www.naturalearthdata.com/) map data. It
    includes functions to download a wealth of Natural Earth vector and
    raster data, including cultural (e.g., country boundaries, airports,
    roads, railroads) and physical (e.g., coastline, lakes, glaciates
    areas) datasets.
-   Modern country boundaries are provided at 2 resolutions by
    `r pkg("rworldmap")` along with functions to join and
    map tabular data referenced by country names or codes. Chloropleth
    and bubble maps are supported and general functions to work on user
    supplied maps (see [A New R package for Mapping Global
    Data](http://journal.R-project.org/archive/2011-1/RJournal_2011-1_South.pdf)
    . Higher resolution country borders are available from the linked
    package `r pkg("rworldxtra")`. Historical country
    boundaries (1946-2012) can be obtained from the
    `r pkg("cshapes")`.
-   `r pkg("marmap")` package is designed for downloading,
    plotting and manipulating bathymetric and topographic data in R. It
    allows to query the ETOPO1 bathymetry and topography database hosted
    by the NOAA, use simple latitude-longitude-depth data in ascii
    format, and take advantage of the advanced plotting tools available
    in R to build publication-quality bathymetric maps (see the
    [PLOS](http://www.plosone.org/article/info%3Adoi%2F10.1371%2Fjournal.pone.0073051)
    paper).
-   `r pkg("maptools")` provides an interface to GSHHS
    shoreline databases.
-   The UScensus2000 suite of packages
    (`r pkg("UScensus2000cdp")`,
    `r pkg("UScensus2000tract")`) makes the use of data from
    the 2000 US Census more convenient.
-   `r pkg("rgbif")` package is used to access Global
    Biodiversity Information Facility (GBIF) occurence data
-   `r pkg("geonames")` is an interface to the
    [www.geonames.org](http://www.geonames.org/) service.
-   `r pkg("OpenStreetMap")` gives access to open street map
    raster images, and `r pkg("osmar")` provides
    infrastructure to access OpenStreetMap data from different sources,
    to work with the data in common R manner, and to convert data into
    available infrastructure provided by existing R packages.
-   `r pkg("tidycensus")` provides access to US Census
    Bureau data in a tidy format, including the option to bind the data
    spatially on import.
-   `r pkg("tigris")` provides access to cartographic
    elements provided by the US Census Bureau TIGER, including
    cartographic boundaries, roads, and water.
-   `r pkg("chilemapas")` provides access to spatial data of
    political and administrative divisions of Chile.
-   `r pkg("geobr")` provides easy access to official
    spatial data sets of Brazil for multiple geographies and years.
-   `r pkg("geouy")` loads and process geographic
    information for Uruguay.
-   `r pkg("rgugik")` allows to search and retrieve data
    from Polish Head Office of Geodesy and Cartography ("GUGiK").
-   `r pkg("giscoR")` provides access to spatial elements
    provided by GISCO - Eurostat, including boundary files of countries,
    NUTS regions, municipalities and other spatial objects.
-   `r pkg("mapSpain")` downloads spatial boundary files of
    administrative regions and other spatial objects of Spain.
-   `r pkg("osmextract")` matches, downloads, converts and
    reads OpenStreetMap data extracts obtained from Geofabrik and other
    providers.

Handling spatial data
---------------------

A number of packages dedicated to spatial data handling have been
written using sp classes.

### Data processing - general

-   `r pkg("rgdal")` and `r pkg("maptools")`.
    The `r pkg("rgeos")` package provides an interface to
    topology functions for `r pkg("sp")` objects using
    [GEOS](http://trac.osgeo.org/geos/) .
-   `r pkg("raster")` package introduces many GIS methods
    that now permit much to be done with spatial data without having to
    use GIS in addition to R.
-   The `r pkg("gdalUtils")` package provides wrappers for
    the Geospatial Data Abstraction Library (GDAL) Utilities.
-   `r pkg("gdistance")`, provides functions to calculate
    distances and routes on geographic grids.
    `r pkg("geosphere")` permits computations of distance
    and area to be carried out on spatial data in geographical
    coordinates. `r pkg("cshapes")` package provides
    functions for calculating distance matrices (see [Mapping and
    Measuring Country
    Shapes](http://journal.R-project.org/archive/2010-1/RJournal_2010-1_Weidmann+Skrede~Gleditsch.pdf)
    ).
-   `r pkg("spsurvey")` provides a range of sampling
    functions.
-   The `r pkg("trip")` package extends sp classes to permit
    the accessing and manipulating of spatial data for animal tracking.
-   `r pkg("magclass")` offers a data class for increased
    interoperability working with spatial-temporal data together with
    corresponding functions and methods (conversions, basic calculations
    and basic data manipulation). The class distinguishes between
    spatial, temporal and other dimensions to facilitate the development
    and interoperability of tools build for it. Additional features are
    name-based addressing of data and internal consistency checks (e.g.
    checking for the right data order in calculations).
-   `r pkg("taRifx")` is a collection of utility and
    convenience functions, and some interesting spatial functions.
-   The `r pkg("rcosmo")` package offers various tools for
    geometric transformations, computations, and statistical analysis of
    spherical data.
-   The `r pkg("areal")` package can be used to interpolate
    overlapping but incongruent polygons, also known as areal weighted
    interpolation.
-   The `r pkg("qualmap")` package can be used to digitize
    qualitative GIS data.

### Data processing - raster and imagery data

-   The `r pkg("landsat")` package with accompanying [JSS
    paper](http://www.jstatsoft.org/v43/i04) provides tools for
    exploring and developing correction tools for remote sensing data.

### Data cleaning

-   `r pkg("cleangeo")` may be used to inspect spatial
    objects, facilitate handling and reporting of topology errors and
    geometry validity issues. It may be used to reduce the likelihood of
    having issues when doing spatial data processing.
-   `r pkg("lwgeom")` may also be used to facilitate
    handling and reporting of topology errors and geometry validity
    issues.

Visualizing spatial data
------------------------

### Base visualization packages

-   Packages such as `r pkg("sp")`,
    `r pkg("sf")`, `r pkg("raster")` and
    `r pkg("rasterVis")` provide basic visualization methods
    through the generic plot function
-   `r pkg("RColorBrewer", priority = "core")` provides very
    useful colour palettes that may be modified or extended using the
    `colorRampPalette` function provided with R.
-   `r pkg("viridis")` also provides colour palettes
    designed with consideration for colorblindness and printing in
    grayscale.
-   `r pkg("classInt", priority = "core")` package provides
    functions for choosing class intervals for thematic cartography.
-   `r pkg("rcosmo")` package provides several tools to
    interactively visualize HEALPix data, in particular, to plot data in
    arbitrary spherical windows.

### Thematic cartography packages

-   `r pkg("tmap")` package provides a modern basis for
    thematic mapping optionally using a Grammar of Graphics syntax.
    Because it has a custom grid graphics platform, it obviates the need
    to fortify geometries to use with ggplot2.
-   `r pkg("quickmapr")` provides a simple method to
    visualize 'sp' and 'raster' objects, allows for basic zooming,
    panning, identifying, and labeling of spatial objects, and does not
    require that the data be in geographic coordinates.
-   `r pkg("cartography")` package allows various
    cartographic representations such as proportional symbols,
    choropleth, typology, flows or discontinuities.
-   The `r pkg("mapmisc")` package is a minimal,
    light-weight set of tools for producing nice looking maps in R, with
    support for map projections.
-   Additional processing and mapping functions are available in
    `r pkg("PBSmapping")` package;
    `r pkg("PBSmodelling")` provides modelling support. In
    addition, `r pkg("GEOmap")` provides mapping facilities
    directed to meet the needs of geologists, and uses the
    `r pkg("geomapdata")` package.

### Packages based on web-mapping frameworks

-   `r pkg("mapview")`, `r pkg("leaflet")` and
    `r pkg("leafletR")` packages provide methods to view
    spatial objects interactively, usually on a web mapping base.
-   `r pkg("RgoogleMaps")` package for accessing Google
    Maps(TM) may be useful if the user wishes to place a map backdrop
    behind other displays.
-   plotGoogleMaps package provides methods for the visualisation of
    spatial and spatio-temporal objects in Google Maps in a web browser.
-   `r pkg("ggmap")` may be used for spatial visualisation
    with Google Maps and OpenStreetMap; `r pkg("ggsn")`
    provides North arrows and scales for such maps.
-   `r pkg("mapedit")` provides an R shiny widget based on
    `r pkg("leaflet")` for editing or creating sf
    geometries.

### Building Cartograms

-   The `r pkg("micromap")` package provides linked
    micromaps using ggplot2.
-   `r pkg("recmap")` package provides rectangular
    cartograms with rectangle sizes reflecting for example population.
-   `r pkg("statebins")` provides a simpler binning approach
    to US states.
-   `r pkg("cartogram")` package allows for constructions of
    a continuous area cartogram by a rubber sheet distortion algorithm,
    non-contiguous Area Cartograms, and non-overlapping Circles
    Cartogram.
-   `r pkg("geogrid")` package turns polygons into
    rectangular or hexagonal cartograms.

Analyzing spatial data
----------------------

### Point pattern analysis

The `r pkg("spatial")` package is a recommended package
shipped with base R, and contains several core functions, including an
implementation of Khat by its author, Prof. Ripley. In addition,
`r pkg("spatstat")` allows freedom in defining the region(s)
of interest, and makes extensions to marked processes and spatial
covariates. Its strengths are model-fitting and simulation, and it has a
useful [homepage](http://www.spatstat.org/) . It is the only package
that will enable the user to fit inhomogeneous point process models with
interpoint interactions. The `r pkg("spatgraphs")` package
provides graphs, graph visualisation and graph based summaries to be
used with spatial point pattern analysis. The
`r pkg("splancs", priority = "core")` package also allows
point data to be analysed within a polygonal region of interest, and
covers many methods, including 2D kernel densities. The
`r pkg("smacpod")` package provides various statistical
methods for analyzing case-control point data. The methods available
closely follow those in chapter 6 of Applied Spatial Statistics for
Public Health Data by Waller and Gotway (2004).

`r pkg("ecespa")` provides wrappers, functions and data for
spatial point pattern analysis, used in the book on Spatial Ecology of
the ECESPA/AEET. The functions for binning points on grids in
`r pkg("ash")` may also be of interest. The ads package
perform first- and second-order multi-scale analyses derived from
Ripley's K-function. The `r pkg("dbmss")` package allows
simple computation of a full set of spatial statistic functions of
distance, including classical ones (Ripley's K and others) and more
recent ones used by spatial economists (Duranton and Overman's Kd,
Marcon and Puech's M). It relies on spatstat for core calculation.

### Geostatistics

The `r pkg("gstat", priority = "core")` package provides a
wide range of functions for univariate and multivariate geostatistics,
also for larger datasets, while
`r pkg("geoR", priority = "core")` and geoRglm contain
functions for model-based geostatistics. Variogram diagnostics may be
carried out with `r pkg("vardiag")`. Automated interpolation
using `r pkg("gstat")` is available in
`r pkg("automap")`. This family of packages is supplemented
by `r pkg("intamap")` with procedures for automated
interpolation. A similar wide range of functions is to be found in the
`r pkg("fields")` package. The
`r pkg("spatial")` package is shipped with base R, and
contains several core functions. The `r pkg("spBayes")`
package fits Gaussian univariate and multivariate models with MCMC.
`r pkg("ramps")` is a different Bayesian geostatistical
modelling package. The `r pkg("geospt")` package contains
some geostatistical and radial basis functions, including prediction and
cross validation. Besides, it includes functions for the design of
optimal spatial sampling networks based on geostatistical modelling. The
`r pkg("rcosmo")` package offers various geostatistics
methods for spherical data: descriptive statistics, entropy based
methods, covariance-variogram methods, etc. Most of rcosmo features were
developed for Cosmic Microwave Background data, but they can also be
used for any spherical data. The `r pkg("FRK")` package is a
tool for spatial/spatio-temporal modelling and prediction with large
datasets. The approach, discussed in Cressie and Johannesson (2008),
decomposes the field, and hence the covariance function, using a fixed
set of n basis functions, where n is typically much smaller than the
number of data points (or polygons) m.

The `r pkg("RandomFields", priority = "core")` package
provides functions for the simulation and analysis of random fields, and
variogram model descriptions can be passed between
`r pkg("geoR")`, `r pkg("gstat")` and this
package. `r pkg("SpatialExtremes")` proposes several
approaches for spatial extremes modelling using
`r pkg("RandomFields")`. In addition,
`r pkg("CompRandFld")`,
`r pkg("constrainedKriging")` and
`r pkg("geospt")` provide alternative approaches to
geostatistical modelling. The `r pkg("spTimer")` package is
able to fit, spatially predict and temporally forecast large amounts of
space-time data using \[1\] Bayesian Gaussian Process (GP) Models, \[2\]
Bayesian Auto-Regressive (AR) Models, and \[3\] Bayesian Gaussian
Predictive Processes (GPP) based AR Models. The
`r pkg("rtop")` package provides functions for the
geostatistical interpolation of data with irregular spatial support such
as runoff related data or data from administrative units. The
`r pkg("georob")` package provides functions for fitting
linear models with spatially correlated errors by robust and Gaussian
Restricted Maximum Likelihood and for computing robust and customary
point and block kriging predictions, along with utility functions for
cross-validation and for unbiased back-transformation of kriging
predictions of log-transformed data. The
`r pkg("SpatialTools")` package has an emphasis on kriging,
and provides functions for prediction and simulation. It is extended by
`r pkg("ExceedanceTools")`, which provides tools for
constructing confidence regions for exceedance regions and contour
lines. The `r pkg("gear")` package implements common
geostatistical methods in a clean, straightforward, efficient manner,
and is said to be a quasi reboot of `r pkg("SpatialTools")`.
The `r pkg("sperrorest")` package implements spatial error
estimation and permutation-based spatial variable importance using
different spatial cross-validation and spatial block bootstrap methods.

The `r pkg("sgeostat")` package is also available. Within
the same general topical area are the
`r pkg("deldir", priority = "core")` and
`r pkg("tripack")` packages for triangulation and the
`r pkg("akima")` package for spline interpolation; the
`r pkg("MBA")` package provides scattered data interpolation
with multilevel B-splines. In addition, there are the
`r pkg("spatialCovariance")` package, which supports the
computation of spatial covariance matrices for data on rectangles, the
`r pkg("regress")` package building in part on
`r pkg("spatialCovariance")`, and the
`r pkg("tgp")` package. The `r pkg("Stem")`
package provides for the estimation of the parameters of a
spatio-temporal model using the EM algorithm, and the estimation of the
parameter standard errors using a spatio-temporal parametric bootstrap.
`r pkg("FieldSim")` is another random fields simulations
package. The `r pkg("SSN")` is for geostatistical modeling
for data on stream networks, including models based on in-stream
distance. Models are created using moving average constructions. Spatial
linear models, including covariates, can be fit with ML or REML. Mapping
and other graphical functions are included. The
`r pkg("ipdw")` provides functions o interpolate -
georeferenced point data via Inverse Path Distance Weighting. Useful -
for coastal marine applications where barriers in the landscape -
preclude interpolation with Euclidean distances.
`r pkg("RSurvey")` may be used as a processing program for
spatially distributed data, and is capable of error corrections and data
visualisation.

### Disease mapping and areal data analysis

`r pkg("DCluster", priority = "core")` is a package for the
detection of spatial clusters of diseases. It extends and depends on the
`r pkg("spdep", priority = "core")` package, which provides
basic functions for building neighbour lists and spatial weights, tests
for spatial autocorrelation for areal data like Moran's I. Functions
for fitting spatial regression models, such as SAR and CAR models prior
to version 1.1-1 are now in
`r pkg("spatialreg", priority = "core")`. These models assume
that the spatial dependence can be described by known weights. In
`r pkg("spatialreg")`, the `ME` and `SpatialFiltering`
functions provide Moran Eigenvector model fitting, as do more modern
functions in the `r pkg("spmoran")` package. The
`r pkg("SpatialEpi")` package provides implementations of
cluster detection and disease mapping functions, including Bayesian
cluster detection, and supports strata. The `r pkg("smerc")`
package provides statistical methods for the analysis of data areal
data, with a focus on cluster detection. The
`r pkg("diseasemapping")` package offers the formatting of
population and case data, calculation of Standardized Incidence Ratios,
and fitting the BYM model using INLA. Regionalisation of polygon objects
is provided by `r pkg("AMOEBA")`: a function to calculate
spatial clusters using the Getis-Ord local statistic. It searches
irregular clusters (ecotopes) on a map, and by `skater` in
`r pkg("spdep")`. The `r pkg("seg")` and
`r pkg("OasisR")` packages provide functions for measuring
spatial segregation; `r pkg("OasisR")` includes Monte Carlo
simulations to test the indices. The `r pkg("spgwr")`
package contains an implementation of geographically weighted regression
methods for exploring possible non-stationarity. The
`r pkg("gwrr")` package fits geographically weighted
regression (GWR) models and has tools to diagnose and remediate
collinearity in the GWR models. Also fits geographically weighted ridge
regression (GWRR) and geographically weighted lasso (GWL) models. The
`r pkg("GWmodel")` package contains functions for -
computing geographically weighted (GW) models. Specifically, basic, -
robust, local ridge, heteroskedastic, mixed, multiscale, generalised -
and space-time GWR; GW summary statistics, GW PCA and GW discriminant
analysis; - associated tests and diagnostics; and options for a range of
distance metrics. The `r pkg("lctools")` package provides
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
`r pkg("sparr")` package provides another approach to
relative risks. The `r pkg("CARBayes")` package implements
Bayesian hierarchical spatial areal unit models. In such models, the
spatial correlation is modelled by a set of random effects, which are
assigned a conditional autoregressive (CAR) prior distribution. Examples
of the models included are the BYM model as well as a recently developed
localised spatial smoothing model. The `r pkg("spaMM")`
package fits spatial GLMMs, using the Matern correlation function as the
basic model for spatial random effects. The
`r pkg("PReMiuM")` package is for profile regression, which
is a Dirichlet process Bayesian clustering model; it provides a spatial
CAR term that can be included in the fixed effects (which are global,
ie. non-cluster specific, parameters) to account for any spatial
correlation in the residuals. The spacom package provides tools to
construct and exploit spatially weighted context data, and further
allows combining the resulting spatially weighted context data with
individual-level predictor and outcome variables, for the purposes of
multilevel modelling. The geospacom package generates distance matrices
from shape files and represents spatially weighted multilevel analysis
results. Spatial survival analysis is provided by the
`r pkg("spBayesSurv")` package: Bayesian Modeling and
Analysis of Spatially Correlated Survival Data. The
`r pkg("spselect")` package provides modelling functions
based on forward stepwise regression, incremental forward stagewise
regression, least angle regression (LARS), and lasso models for
selecting the spatial scale of covariates in regression models.

### Spatial regression

The choice of function for spatial regression will depend on the support
available. If the data are characterised by point support and the
spatial process is continuous, geostatistical methods may be used, or
functions in the `r pkg("nlme")` package. If the support is
areal, and the spatial process is not being treated as continuous,
functions provided in the `r pkg("spatialreg")` package may
be used. This package can also be seen as providing spatial econometrics
functions, and, as noted above, provides basic functions for building
neighbour lists and spatial weights, tests for spatial autocorrelation
for areal data like Moran's I, and functions for fitting spatial
regression models. `r pkg("spdep")` provides the full range
of local indicators of spatial association, such as local Moran's I and
diagnostic tools for fitted linear models, including Lagrange Multiplier
tests. Spatial regression models that can be fitted using maximum
likelihood and Bayesian MCMC methods in
`r pkg("spatialreg")` include spatial lag models, spatial
error models, and spatial Durbin models. For larger data sets, sparse
matrix techniques can be used for maximum likelihood fits, while spatial
two-stage least squares and generalised method of moments estimators are
an alternative. When using GMM, `r pkg("sphet")` can be used
to accommodate both autocorrelation and heteroskedasticity. The
`r pkg("splm")` package provides methods for fitting spatial
panel data by maximum likelihood and GM. 
The `r pkg("spsur")` package provides functions to test and estimate 
spatial seemingly unrelated regression models (spatial SUR) by maximum likelihood and three-stage least squares. The two small packages
`r pkg("S2sls")` and `r pkg("spanel")` provide
alternative implementations without most of the facilities of
`r pkg("splm")`. The HSAR package provides Hierarchical
Spatial Autoregressive Models (HSAR), based on a Bayesian Markov Chain
Monte Carlo (MCMC) algorithm. `r pkg("spatialprobit")` make
possible Bayesian estimation of the spatial autoregressive probit model
(SAR probit model). The ProbitSpatial package provides methods for
fitting Binomial spatial probit models to larger data sets; spatial
autoregressive (SAR) and spatial error (SEM) probit models are included.
The `r pkg("starma")` package provides functions to
identify, estimate and diagnose a Space-Time AutoRegressive Moving
Average (STARMA) model.

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
    provided by `r pkg("trip")`. `r pkg("ncf")`
    has entered CRAN recently, and provides a range of spatial
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
-   [Rcitrus](http://www.leg.ufpr.br/Rcitrus/) (off-CRAN package) is for
    the spatial analysis of plant disease incidence.
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



### Links
-   [R-SIG-Geo mailing list](https://stat.ethz.ch/mailman/listinfo/R-SIG-Geo/)
