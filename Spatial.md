---
name: Spatial
topic: Analysis of Spatial Data
maintainer: Roger Bivand, Jakub Nowosad, Krzysztof Dyba
email: Roger.Bivand@nhh.no, nowosad.jakub@gmail.com, adres7@gmail.com
version: 2025-11-11
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
is now clearly aligned with the modern `r pkg("sf", priority="core")` and 
`r pkg("stars", priority="core")` packages. Rspatial has grown from the 
`r pkg("raster")` package, now moving towards the modern 
`r pkg("terra", priority="core")` package. It is also worth noting the
wealth of online book projects, which may be helpful for users seeking
an introduction, including 
[Geocomputation with R](https://r.geocompx.org/index.html).

Specific questions or issues may be raised where 
`packageDescription(<pkg>)$BugReports` returns an URL for bug 
reports or issues (where `<pkg>` is the
name of the package as a string), or directly with package
maintainers by email. Use may also be made of
the [R-SIG-Geo](https://stat.ethz.ch/mailman/listinfo/R-SIG-Geo/)
mailing-list after subscription, or of 
[Stack Overflow](https://stackoverflow.com) with appropriate tags, or of
[Stack Exchange](https://gis.stackexchange.com/). Using the
`#rspatial` tag on [Mastodon](https://joinmastodon.org/)
may also be worth trying, or browsing traffic using that tag (among others);
BlueSky currently is emerging.

The packages in this view can be roughly structured into the following
topics. If you think that some package is missing from the list, please
e-mail the maintainer or submit an issue or pull request in the GitHub
repository linked above.

**Table of contents**

- [Classes for spatial data and metadata](#classes-for-spatial-data-and-metadata)
  - [Spatial data - general](#spatial-data---general)
  - [Raster data](#raster-data)
  - [Geographic metadata](#geographic-metadata)
- [Reading and writing spatial data](#reading-and-writing-spatial-data)
  - [Reading and writing spatial data - data formats](#reading-and-writing-spatial-data---data-formats)
  - [Reading and writing spatial data - GIS software connectors](#reading-and-writing-spatial-data---gis-software-connectors)
  - [Specific geospatial data sources of interest](#specific-geospatial-data-sources-of-interest)
  - [Interfaces to Spatial Web-Services](#interfaces-to-spatial-web-services)
- [Handling spatial data](#handling-spatial-data)
  - [Data cleaning](#data-cleaning)
  - [Data processing - general](#data-processing---general)
  - [Data processing - specific](#data-processing---specific)
  - [Remote sensing](#remote-sensing)
  - [Spatial sampling](#spatial-sampling)
- [Visualizing spatial data](#visualizing-spatial-data)
  - [Base visualization packages](#base-visualization-packages)
  - [Thematic cartography packages](#thematic-cartography-packages)
  - [Packages based on web-mapping frameworks](#packages-based-on-web-mapping-frameworks)
  - [Building cartograms](#building-cartograms)
- [Analyzing spatial data](#analyzing-spatial-data)
  - [Point pattern analysis](#point-pattern-analysis)
  - [Geostatistics](#geostatistics)
  - [Disease mapping and areal data analysis](#disease-mapping-and-areal-data-analysis)
  - [Spatial regression](#spatial-regression)
  - [Ecological analysis](#ecological-analysis)
  - [Machine learning of spatial data](#machine-learning-of-spatial-data)
- [Installing packages linking to PROJ, GDAL or GEOS](#installing-packages-linking-to-proj-gdal-or-geos)

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
spatial classes in R](https://geocompx.org/post/2021/spatial-classes-conversion/).

Complementary initiatives are ongoing to support better handling of
geographic metadata in R.

### Spatial data - general

-   `r pkg("sf", priority = "core")` is a CRAN package for spatial vector data, 
    providing Simple Features
    for R, in compliance with the [OGC Simple
    Feature](http://www.opengeospatial.org/standards/sfa) standard. The
    development of the package was supported by the [R
    Consortium](https://www.r-consortium.org/). It provides simple
    features access for vector data, and as such is a modern
    implementation and standardization of parts of the legacy 
    `r pkg("sp")` package. `r pkg("sf")` is documented in an [R
    Journal](https://journal.R-project.org/archive/2018/RJ-2018-009/index.html)
    article. `r pkg("sf")` uses the [PROJ](https://proj.org/), [GEOS](https://libgeos.org/) 
    and [GDAL](https://gdal.org/) external software
    libraries, which must be available for source installs together with
    other external software libraries that they in turn depend on.
-   `r pkg("stars", priority = "core")` is being actively developed and 
    was initially supported by the [R Consortium](https://www.r-consortium.org/); 
    it provides for
    spatiotemporal data in the form of dense arrays. It supercedes the 
    `r pkg("spacetime")` package, which 
    extended the shared classes defined in `r pkg("sp")` for
    spatio-temporal data (see [Spatio-Temporal Data in
    R](http://www.jstatsoft.org/v51/i07)). `r pkg("stars")` uses PROJ and
    GDAL through `r pkg("sf")`.
-   `r pkg("terra", priority = "core")` provides classes for spatial vector and
    raster data, linking directly to PROJ, GDAL and GEOS. 
    access to GDAL functionality for R packages. 
-   `r pkg("gdalraster")` provides API bindings to GDAL. Since
    `r pkg("gdalraster")` 2.0.0, bindings are provided for both the Raster and
    Vector [APIs](https://gdal.org/en/stable/api/index.html), the Geometry API
    ([GEOS](https://libgeos.org/) via GDAL headers) and Spatial Reference
    Systems API ([PROJ](https://proj.org/) via GDAL). Bindings are provided for
    the low-level Virtual Systems Interface
    ([VSI](https://gdal.org/en/stable/api/cpl.html#cpl-vsi-h)) which abstracts
    file system operations and binary I/O on URLs, cloud storage services,
    compressed files (.zip, .gz, .tar, .tar.gz, .7z archives) and in-memory
    files, as well as regular file systems. When built against GDAL >= 3.11.3,
    `r pkg("gdalraster")` also provides bindings to GDAL's unified command
    line interface framework, enabling
    [the use of CLI algorithms](https://usdaforestservice.github.io/gdalraster/articles/use-gdal-cli-from-r.html)
    and pipeline processing from R.
-   The `r pkg("vapour")` package also offers low-level access to GDAL
    functionality for R packages. 
-   The `r pkg("spatstat", priority = "core")` contains classes suited to the
    analysis of point patterns, and may be coerced to and from `"sf"`, `"stars"`
    and other spatial classes.


### Raster data

-   `r pkg("terra", priority = "core")` is a re-implementation of
    `r pkg("raster")` functionality and introducing new S4 classes for raster and
    vector data. See the [manual and
    tutorials](https://rspatial.org/terra/) to get started.
    `r pkg("terra")` is very similar to the `r pkg("raster")` package; but
    `r pkg("terra")` is simpler and faster.
-   `r pkg("stars", priority = "core")` provides for spatiotemporal data in the
    form of dense arrays, with space and time being array dimensions.
    Examples include socio-economic or demographic data, environmental
    variables monitored at fixed stations, time series of satellite
    images with multiple spectral bands, spatial simulations, and
    climate model results.
-   The `r pkg("gdalcubes")` package also provides classes for data cubes,
    including proxy data cubes, it links to PROJ, GDAL and NetCDF4.
-   `r pkg("gdalraster")` provides an R implementation of the [GDAL Raster Data
    Model](https://gdal.org/en/stable/user/raster_data_model.html), along with
    several utilities and algorithms for processing and analyzing raster data.

### Geographic metadata

-   `r pkg("geometa")` provides classes and methods to write
    geographic metadata following the ISO and OGC metadata standards
    (ISO 19115, 19110, 19119) and export it as XML (ISO 19139) for later
    publication into metadata catalogues. Reversely, geometa provides a
    way to read ISO 19139 metadata into R. The package extends
    `r pkg("sf")` to provide GML (ISO 19136) representation
    of geometries.
-   `r pkg("ncdf4")` provides read and write functions for
    handling metadata (CF conventions) in the self-described NetCDF
    format.
-   `r pkg("CFtime")` encapsulates the CF time coordinate and allows to deal
    with the different CF calendars

Reading and writing spatial data
--------------------------------

Spatial data is most often represented by one of two data models, 
vector or raster, and both models have many of their own file formats.
[GDAL (Geospatial Data Abstraction Library)](https://gdal.org/) is a (non-R)
library that provides a unified way to read and write hundreds of 
spatial data formats. Formats supported by GDAL include both OGC standard 
data formats (e.g., GeoPackage) and proprietary formats (e.g., ESRI Shapefile).
GDAL is used by a large number of GIS software and also many R packages,
such as `r pkg("sf")`, `r pkg("terra")`, and `r pkg("vapour")`, and the
`r pkg("gdalraster")` package explicitly implements the GDAL [Raster](https://gdal.org/en/stable/user/raster_data_model.html)
and [Vector](https://gdal.org/en/stable/user/vector_data_model.html) Data
Models. This allows us to read and write spatial data in R from and to various
spatial file formats. Important note: CRAN offers binary versions of packages
`r pkg("gdalraster")`, `r pkg("sf")`, `r pkg("terra")`, and `r pkg("vapour")`
for Windows and macOS, that contain specific GDAL version with a subset of
possible data source drivers. If other drivers are needed, you need to either
use other conversion utilities or install these packages from the source against
a version of GDAL with the required drivers.

In the past, `r rforge("rgdal")` and `r pkg("raster")` (through `r rforge("rgdal")`) 
were recommended for reading and writing of spatial data in R.
However, due to [the retirement of rgdal on 16 October 2023](https://stat.ethz.ch/pipermail/r-sig-geo/2023-October/029350.html)
new projects should not use it, and existing projects should implement migration
to the packages mentioned in the previous paragraph. In addition, `r rforge("rgeos")` 
and `r rforge("maptools")` were archived at the same time. Further details and links 
may be found in [project reports](https://r-spatial.github.io/evolution/) on the evolution
project. From October 2023, `r pkg("sp")` only uses methods 
from `r pkg("sf")` in place of those from `r rforge("rgdal")` for projection 
and access to the underlying definitions of coordinate reference systems.

### Reading and writing spatial data - data formats

Other packages provide facilities to read and write spatial data, 
dealing with open standard formats or proprietary formats.

*Open formats*

-   *Well-Known Text (WKT) / Well-Known Binary (WKB):* These standards are 
    part of the OGC Simple Feature specification. Both WKT/WKB formats are 
    supported by the `r pkg("sf")` package that implements the whole
    OGC Simple Feature specification in R. Additionally, `r pkg("wk")`
    may be used to parse well-known binary and well-known text
    representation of geometries to and from R-native formats.
-   *GeoJSON:* An rOpenSci [blog entry](http://ropensci.org/blog/blog/2016/11/22/geospatial-suite) 
    describes a GeoJSON-centred approach to reading GeoJSON and WKT data.
    The entry lists `r pkg("geojson")`, and `r pkg("geojsonio")`, among others.
    The GeoJSON format can also be read and written with `r pkg("gdalraster")`,
    `r pkg("sf")`, `r pkg("terra")`, and `r pkg("vapour")`. 
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
    cloud data and can be read/write with `r pkg("lidR")`.

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

### Reading and writing spatial data - GIS software connectors

-   *PostGIS:* The `r pkg("rpostgis")` package provides additional functions
    to the `r pkg("RPostgreSQL")` package to interface R with a 
    'PostGIS'-enabled database, as well as convenient wrappers to 
    common 'PostgreSQL' queries. It is documented in an 
    [R Journal](https://journal.R-project.org/archive/2018/RJ-2018-025/index.html)
    article.
    `r pkg("sf")` also provides an R interface to PostGIS, 
    for both reading and writing, through GDAL.
-   *GRASS GIS:* Integration with version 7.\* and 8.\* of the leading open source GIS,
    GRASS GIS, is provided in CRAN package `r pkg("rgrass")`, which uses `r pkg("terra")` for file transfer.
-   *SAGA GIS:* `r pkg("RSAGA")` and `r pkg("Rsagacmd")` offer shell-based
    wrapper for SAGA GIS commands.
-   *QGIS:* QGIS version 2 was supported by RQGIS 
    (`r github("r-spatial/RQGIS")`). Using QGIS processing algorithms is currently
    supported by `r pkg("qgisprocess")`, which uses the standalone
    'qgis_process' command-line utility from QGIS (use
    recent QGIS versions; may work since >= 3.16). Both native QGIS and third-party
    (plugin) processing providers are supported, i.e. GRASS, SAGA, GDAL, ...
-   *WhiteboxTools:* `r pkg("whitebox")` is an R frontend for 
    the WhiteboxTools software.
-   *ArcGIS:* `r pkg("RPyGeo")` is a wrapper for Python access 
    to the ArcGIS GeoProcessor. The ESRI company also offers their 
    own package (`r github("R-ArcGIS/r-bridge")`) that allows transferring data
    from ArcGIS to R.
-   *DuckDB*: `r pkg("duckspatial")` is a package that allows and simplifies the process of writing
    and reading vector data into a DuckDB database through the SPATIAL extension.
-   Various GIS Software, including Orfeo ToolBox and SAGA GIS, can also be
    connected to R using `r pkg("link2GI")`.
-   *Orfeo ToolBox* segmentation module: `r pkg("OTBsegm")` simplifies the process of applying 
    unsupervised segmentation algorithms available in *Orfeo Toolbox* via the `r pkg("link2GI")`

### Specific geospatial data sources of interest

-   `r pkg("rnaturalearth")` package facilitates interaction with 
    [Natural Earth](http://www.naturalearthdata.com/) map data. It includes 
    functions to download a wealth of Natural Earth vector and raster data,
    including cultural (e.g., country boundaries, airports, roads, railroads)
    and physical (e.g., coastline, lakes, glaciated areas) datasets.
-   `r pkg("elevatr")` provides access to elevation data from several web services. 
-   Historical country boundaries (1886-today) can be obtained 
    from the `r pkg("cshapes")`.
-   `r pkg("marmap")` package is designed for downloading, plotting, 
    and manipulating bathymetric and topographic data in R. It allows to query 
    the ETOPO1 bathymetry and topography database hosted by the NOAA, 
    use simple latitude-longitude-depth data in ASCII format, and take 
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
-   `r pkg("osmapiR")` is an R interface to the [OpenStreetMap API](https://wiki.openstreetmap.org/wiki/API_v0.6) for 
    fetching and saving all kinds of data from/to the OpenStreetMap database, 
    including map objects, GPS traces, notes, changesets, and users.
-   `r pkg("OpenStreetMap")` gives access to open street map raster images.
-   `r pkg("giscoR")` provides access to spatial elements provided by 
    GISCO - Eurostat, including boundary files of countries, NUTS regions,
    municipalities, and other spatial objects.
-   `r pkg("chilemapas")` provides access to spatial data of political 
    and administrative divisions of Chile.
-   `r pkg("geobr")` provided easy access to official spatial data sets of 
    Brazil for multiple geographies and years.
-   `r pkg("geouy")` loads and process geographic information for Uruguay.
-   `r pkg("RCzechia")` downloads spatial boundary files of administrative regions
    and other spatial objects of the Czech Republic.
-   `r pkg("rgugik")` allows to search and retrieve data from Polish Head 
    Office of Geodesy and Cartography ("GUGiK").
-   `r pkg("mapSpain")` downloads spatial boundary files of administrative 
    regions and other spatial objects of Spain.
-   `r pkg("mapme.biodiversity")` allows to download and process a number open 
    datasets related to biodiversity conservation providing efficient routines 
    and parallelization options. Datasets include among others the 
    [Global Forest Watch](https://www.globalforestwatch.org/), 
    [ESA/Copernicus Landcover](https://land.copernicus.eu/global/products/lc), 
    [WorldClim](https://www.worldclim.org/) and 
    [NASA FIRMS](https://firms.modaps.eosdis.nasa.gov/active_fire/).
-    `r pkg("terrainr")` provides an interface to the
    [United States Geological Survey's National Map services](https://apps.nationalmap.gov/services/),
    providing elevation data and orthoimagery along other basemap tiles
    for the United States.
-   `r pkg("geodata")` facilitates access to climate, elevation, soil, crop, 
    species occurrence, and administrative boundary data, and is a successor of
    the `getData()` function from the `r pkg("raster")` package.
-   `r pkg("forestdata")` allows to download forest and land cover data from
    various sources, includying forest inventory data, forest cover maps, and
    global canopy height models.

### Interfaces to Spatial Web-Services

Some R packages focused on providing interfaces to web-services and web tools 
in support of spatial data management. Here follows a first tentative 
(non-exhaustive) list:

-   `r pkg("ows4R")` is a package that intends to provide an R interface 
    to OGC standard Web-Services. It is in active development and currently 
    support interfaces to the 
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

Handling spatial data
---------------------

### Data cleaning

-   `r pkg("sf")` has a built-in functions `st_is_valid` to check whether 
    a sf geometry is valid and `st_make_valid` to fix invalid geometry (from GEOS 3.8).
-   `r pkg("lwgeom")` may also be used to facilitate handling and reporting
    of topology errors and geometry validity issues
    in sf objects.

### Data processing - general

-   `r pkg("sf")` provides an interface to spatial geometry functions using
    the [GEOS](https://libgeos.org/) and [S2](https://s2geometry.io/) libraries; S2
    is bundled in the `r pkg("s2")` package which `r pkg("sf")` methods use for
    topological predicates and operations on spherical or elliptical coordinates by default.
-   `r pkg("stars")` contains tools for manipulating raster and vector
    data cubes.
-   `r pkg("terra")` package introduces many GIS methods for spatial vector
    and raster data.
-   The `r github("cran/gdalUtils")` (see
    <https://stat.ethz.ch/pipermail/r-sig-geo/2022-April/028953.html>)
    and `r pkg("gdalUtilities")` packages provide
    wrappers for the Geospatial Data Abstraction Library (GDAL) Utilities.
-   `r pkg("gdalraster")` provides access to "traditional" GDAL utilities and
    algorithms via API bindings. When built against GDAL >= 3.11.3, bindings to
    GDAL's unified CLI framework provide access to CLI algorithms and GDAL
    facilities for pipeline processing (see
    [Using gdal CLI algorithms from R](https://usdaforestservice.github.io/gdalraster/articles/use-gdal-cli-from-r.html)).
-   The `r pkg("geos")` high-performance bindings to the GEOS library, based on
    `r pkg("libgeos")`; the latter bundles a frozen copy of GEOS, and does not
    link to system versions, which may possibly be different versions of GEOS.
-   `r pkg("rmapshaper")` is a wrapper around the 'mapshaper' 'JavaScript'
    library to perform topologically-aware polygon simplification and other 
    operations such as clipping, erasing, dissolving, and converting 
    'multi-part' to 'single-part' geometries.
-   `r pkg("gdistance")` and `r pkg("spaths")` provide functions to calculate distances and 
    routes on geographic grids.
-   `r pkg("geosphere")` permits computations of distance and area to be 
    carried out on spatial data in geographical coordinates.
-   `r pkg("dggridR", priority = "core")` provides a discrete global grid system
    via DGGRID. These grids are useful for spatial statistics because
    they tile the Earth with _equally_-sized hexagons, triangles, or diamonds.
-   `r pkg("cshapes")` package provides functions for calculating distance
    matrices (see [Mapping and Measuring Country Shapes](http://journal.R-project.org/archive/2010-1/RJournal_2010-1_Weidmann+Skrede~Gleditsch.pdf)).
-   `r pkg("magclass")` offers a data class for increased interoperability 
    working with spatial-temporal data together with corresponding functions
    and methods (conversions, basic calculations and basic data manipulation).
-   The `r pkg("trip")` package extends spatial classes to permit the
    accessing and manipulating of spatial data for animal tracking.
<!-- Roger, should we link here to papers such as https://besjournals.onlinelibrary.wiley.com/doi/10.1111/1365-2656.13116, https://link.springer.com/article/10.1007%2Fs40823-021-00067-y, https://link.springer.com/article/10.1007/s10109-020-00342-2, etc?-->

### Data processing - specific

-   The `r pkg("areal")` package can be used to interpolate overlapping 
    but incongruent polygons, also known as areal weighted interpolation.
-   The `r pkg("centerline")` package can be used for centerline (or median-axis) estimation of closed polygons,
    such as lakes, landslides, or rivers. The package also
    provides extensions to `r pkg("ggplot2")`, allowing users to place `geom_text` and `geom_label` at the center of a spatial
    polygon, which can be useful for map-making.
-   The `r pkg("qualmap")` package can be used to digitize qualitative GIS data.
-   The `r pkg("exactextractr")` for fast and accurate summary of raster values
    on polygonal areas (known as zonal statistics).

### Remote sensing

-   `r pkg("rstac")` provides functions to access, search and download
    spacetime earth observation data via [SpatioTemporal Asset Catalogs](https://stacspec.org).
    This package supports the version 1.0.0 (and older) of the 
    [STAC specification](https://github.com/radiantearth/stac-spec).
-   The `r pkg("rsi")` package provides an interface to the
    [Awesome Spectral Indices](https://github.com/awesome-spectral-indices/awesome-spectral-indices)
    project. It also provides functions to download, mask, and composite data 
    from [SpatioTemporal Asset Catalogs](https://stacspec.org), with a particular 
    focus on satellite imagery.
-   The `r pkg("RStoolbox")` is a remote sensing toolbox covering many aspects
    including data import, pre-processing, data analysis, image classification
    and graphical display.
-   The `r pkg("sits")` is an end-to-end toolkit for land use and land cover
    classification using big Earth observation data, based on machine learning
    methods applied to satellite image data cubes.
-   The `r pkg("landsat")` package with accompanying 
    [JSS paper](http://www.jstatsoft.org/v43/i04) provides tools for exploring
    and developing correction tools for remote sensing data.
-   `r github("rspatial/luna")` has tools for acquiring and processing satellite 
    remote sensing data from NASA's LANDSAT and MODIS data sources. 
    `r pkg("MODISTools")` also provides an interface to the 
    [MODIS Land Products Subsets](https://modis.ornl.gov/data/modis_webservice.html) 
    web services.
    `r pkg("modisfast")`, on its side, provides an interface to some MODIS, VIIRS 
    and GPM [OPeNDAP](https://www.opendap.org/) servers, enabling to subset the datacubes
    directly at the download phase. 
-   The `r pkg("CDSE")` provides the interface to the 'Copernicus Data Space Ecosystem'
    API (<https://dataspace.copernicus.eu/analyse/apis/sentinel-hub>), mainly for searching
    the catalogue of available data from Copernicus Sentinel missions and obtaining the
    images for the area of interest based on selected spectral bands.
-   The `r pkg("rsat")` is designed for automatically downloading, pre-processing, and
    managing time series of satellite imagery from multiple platforms (i.e., Landsat,
    Sentinel and MODIS). 
-   `r pkg("rgee")` is an [Earth Engine](https://earthengine.google.com/) 
    client library for R. All of the 'Earth Engine' API classes, modules, 
    and functions are made available. Additional functions implemented include 
    importing (exporting) of Earth Engine spatial objects, 
    extraction of time series, interactive map display, 
    assets management interface, and metadata display.
-   The `r pkg("openeo")` is an R client package that allows users to interact with
    openEO-compliant back-ends (<https://openeo.org/>) for processing Earth observation
    data (e.g., services offered by European Space Agency or VITO Remote Sensing).
    It enables to define and execute workflows remotely on cloud infrastructures
    without having to download large amounts of data locally.

### Spatial sampling

-   `r pkg("spsurvey")` provides functions to select generalized
    random-tessellation stratified (GRTS) probability samples
    and analyze survey data.
-   `r pkg("Spbsampling")` allows selecting probability samples well spread
    over the population of interest, in any dimension and using 
    any distance function.
-   `r pkg("spatialsample")` is a member of the tidymodels family of packages
    and contains functions and classes for spatial resampling to use with
    `r pkg("rsample")`. 
-   `r pkg("MBHdesign")` provides spatially survey balanced designs using the
    quasi-random number method.

Visualizing spatial data
------------------------

### Base visualization packages

-   Packages such as `r pkg("sf")`, `r pkg("stars")`, `r pkg("terra")`, 
    and `r pkg("rasterVis")` provide basic visualization methods through 
    the generic plot function.
-   `r pkg("classInt", priority = "core")` package provides functions for 
    choosing class intervals for thematic cartography.
-   Currently, the grDevices package (included with the R installation) 
    contains a large number of color palettes that can be accessed with 
    the `hcl.colors` and `palette.colors` functions; see also New features in this
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
    `geom_sf` function and additional support for stars object is available
    through the `geom_stars` function available in the `r pkg("stars")` package.
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
    and OpenStreetMap.
-   `r pkg("mapedit")` provides an R shiny widget based on `r pkg("leaflet")`
    for editing or creating sf geometries.

### Building cartograms

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
topical cluster of packages. Moreover, machine learning for spatial data 
is a growing area of interest, and there are packages that
provide tools for this.
All of these approaches to analysing spatial data
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
    Marcon and Puech's M). It relies on `r pkg("spatstat")` for core 
    calculation.
-   `r pkg("sfhotspot")` provides functions for descriptive hotspot analysis and
    mapping of spatial concentrations of points, including kernel-density 
    estimation, Getis--Ord Gi*, hotspot classification, etc. The package 
    attempts to provide sensible default values to assist with analysis
    and mapping for students and non-experts.

### Geostatistics

Geostatistics uses a model fitted using the distances between observations to
interpolate values observed at point to unobserved points

-   The `r pkg("gstat", priority = "core")` package provides a
    wide range of functions for univariate and multivariate geostatistics,
    also for larger datasets.
-   `r pkg("mbg")` offers an [interface](https://henryspatialanalysis.github.io/mbg/)
    for running model-based geostatistics with optional covariates and discrete spatial
    effects.
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
-   The `r pkg("FRK")` package is a tool for 
    spatial/spatio-temporal modelling and prediction with large
    datasets. The approach, discussed in Cressie and Johannesson (2008),
    decomposes the field, and hence the covariance function, using a fixed
    set of n basis functions, where n is typically much smaller than the
    number of data points (or polygons) m.
-   `r pkg("SpatialExtremes")` proposes several
    approaches for spatial extremes. 
-   In addition, `r pkg("constrainedKriging")` and `r pkg("geospt")` provide 
    alternative approaches to geostatistical modelling. 
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
-   The archived `Stem`
    package provided for the estimation of the parameters of a
    spatio-temporal model using the EM algorithm, and the estimation of the
    parameter standard errors using a spatio-temporal parametric bootstrap.
-   The `r pkg("SSN2")` is for geostatistical modeling
    for data on stream networks, including models based on in-stream
    distance. Models are created using moving average constructions. Spatial
    linear models, including covariates, can be fit with ML or REML. Mapping
    and other graphical functions are supported.
-   The `r pkg("ipdw")` provides functions to interpolate
    georeferenced point data via Inverse Path Distance Weighting. Useful
    for coastal marine applications where barriers in the landscape
    preclude interpolation with Euclidean distances.
-   `r pkg("sptotal")` uses Finite Population Block Kriging (FPBK) to
    provide a prediction for a quantity of interest, most commonly a
    population total or a prediction of total abundance, on a finite number
    of spatial sites.
-   `r pkg("spmodel")` fits statistical models to geostatistical
    and areal spatial data using a variety of covariance structures. 
    Additional functionality allows for prediction (Kriging), non-spatial
    random effects, anisotropy, and big data.
-   `r pkg("spStack")` uses predictive stacking to fit Bayesian spatial and 
    spatial-temporal models for point-referenced Gaussian, Poisson, binomial, 
    and binary data without using MCMC.

### Disease mapping and areal data analysis

Both point pattern analysis and geostatistics enter into disease mapping, which
is concerned with representing public health information over space and
time in a communicative and responsible way. Estimation is important to present
calculated rates that are comparable both in terms of levels and uncertainty.

-   `r pkg("DCluster", priority = "core")` is a package for the
    detection of spatial clusters of diseases. It is complemented by
    `r pkg("DClusterm")` for model-based cluster detection, and by 
    `r github("tkhrotn/rflexscan")` and `r pkg("FlexScan")`, two implementations
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
-   A Markov Random Field `"mrf"` effect may be added to models in the 
    `r pkg("mgcv")` package shipped with base R, providing flexible
    modelling tools in a recommended package.
-   The `r pkg("hglm")` package also provides SAR and CAR model fitting
    approaches.
-   Regionalization of polygon objects
    is no longer provided by archived `AMOEBA`: a function to calculate
    spatial clusters using the Getis-Ord local statistic. It searched
    for irregular clusters (ecotopes) on a map. `skater()` in
    `r pkg("spdep")` does not interpose a local statistic,
     being based on distance between features in attribute space,
     and polygon contiguity. 
-   The `r pkg("divseg")` and
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
    `r github("alexWhitworth/synthACS")` permits the building and
    running of spatially explicit agent-based models.
-   The `r pkg("geostan")` package has GLMs, SAR, proper CAR, ICAR, and eigenvector
    spatial filter (ESF) models for Bayesian disease mapping and spatial regression.
    The package uses the Stan modeling language for MCMC analysis.
    The package also contains exploratory spatial analysis tools (Moran scatter plot and
    various measures of spatial autocorrelation) and measurement error models
    designed for the use of (noisy) survey estimates as covariates. 
-   `r pkg("waywiser")` helps assess models fit to spatial data, with
    functions for calculating the spatial autocorrelation of model
    residuals, for calculating model performance statistics, for 
    assessing model performance across multiple spatial scales, and for
    calculating the "area of applicability" of a model. Functions are
    designed to be compatible with both base R and with the tidymodels
    modeling framework, and adopt `r pkg("yardstick")` classes and
    interfaces.
-   The `r pkg("gdverse")` package provides 
    an integrated and extensible toolkit for analyzing spatial stratified heterogeneity 
    and spatial associations using the geographical detector methodology. It includes 
    core functions for computing the q-statistic, performing spatial factor exploration, 
    and detecting spatial interactions between variables. The package supports various 
    geographical detector models and provides high-performance implementations for their
    efficient execution.


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
-   The two small archived packages `S2sls` and `spanel` provide alternative
    implementations without most of the facilities of `r pkg("splm")`. 
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
-   `r pkg("varycoef")` and `r pkg("spBayes")` provide implementations of
    spatially varying coefficient (SVC) models, which may be preferred to
    geographically weighted regression (GWR) models as having proper
    statistical foundations.
-   The `r pkg("gwrr")` package fits geographically weighted
    regression (GWR) models and has tools to diagnose and remediate
    collinearity in the GWR models. It also fits geographically weighted ridge
    regression (GWRR) and geographically weighted lasso (GWL) models. The
    `r pkg("GWmodel")` package contains functions for
    computing geographically weighted (GW) models. Specifically, basic, 
    robust, local ridge, heteroskedastic, mixed, multiscale, generalised 
    and space-time GWR; GW summary statistics, GW PCA and GW discriminant
    analysis; associated tests and diagnostics; and options for a range of
    distance metrics. 
-   `r pkg("waywiser")` helps assess models fit to spatial data, with
    functions for calculating the spatial autocorrelation of model
    residuals, for calculating model performance statistics, for 
    assessing model performance across multiple spatial scales, and for
    calculating the "area of applicability" of a model. Functions are
    designed to be compatible with both base R and with the tidymodels
    modeling framework, and adopt `r pkg("yardstick")` classes and
    interfaces.
-   `r pkg("spStack")` fits Bayesian spatially-temporally varying 
    coefficients (STVC) generalized linear models without MCMC using stacking of
    predictive densities.

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
<!-- -   The `r github("ArturoTorres/Watersheds")` package provides methods for -->
    <!-- watersheds aggregation and spatial drainage network analysis. -->
-   The `r pkg("ngspatial")` package provides tools for
    analyzing spatial data, especially non-Gaussian areal data. It
    supports the sparse spatial generalized linear mixed model of Hughes
    and Haran (2013) and the centered autologistic model of Caragea and
    Kaiser (2009).
-   `r pkg("landscapemetrics")` package calculates landscape
    metrics for categorical landscape patterns. It can be used as a
    drop-in replacement for
    [FRAGSTATS](https://www.umass.edu/landeco/research/fragstats/fragstats.html),
    as it offers a reproducible workflow for landscape analysis in a
    single environment. It also provides several visualization
    functions, e.g. to show all labeled patches or the core area of all
    patches.
-   `r pkg("waywiser")` helps assess models fit to spatial data, with
    functions for calculating the spatial autocorrelation of model
    residuals, for calculating model performance statistics, for 
    assessing model performance across multiple spatial scales, and for
    calculating the "area of applicability" of a model. Functions are
    designed to be compatible with both base R and with the tidymodels
    modeling framework, and adopt `r pkg("yardstick")` classes and
    interfaces.
-   `r pkg("dismo")` provides functions for species distribution modelling.

The `r view("Environmetrics")` Task View contains a much more
complete survey of relevant functions and packages.

### Machine learning of spatial data

Machine learning of spatial data requires specialized methods to account for spatial dependencies like spatial autocorrelation -- where nearby observations tend to be similar. 
Ignoring these dependencies during model training and evaluation risks information leakage.
For example, randomly splitting spatial data into training and testing subsets, without considering spatial autocorrelation, can result in test samples being spatially close to training samples. 
This violates the assumption of independence between training and test sets and can lead to inflated performance metrics and poor model generalization.

To address this, various approaches and methods to account for spatial dependencies and relationships when building models were developed.
In general, machine learning of spatial data can be performed through one of the existing machine learning frameworks in R, such as `r pkg("caret")`, `r pkg("mlr3")`, and `r pkg("tidymodels")` or through specialized spatial machine learning packages.

-    The `r pkg("caret")` package provides a consistent interface for training models but requires additional packages like `r pkg("blockCV")` or `r pkg("CAST")` to implement spatial methodologies. Functions like `CAST::knndm` and `CAST::ffs` enable spatially aware feature selection and cross-validation, while `CAST::aoa` assesses the area of applicability for spatial models. The `r pkg("mbg")` package offers convenience functions for fitting `r pkg("caret")` models with point-referenced outcomes and raster features.
-   `r pkg("mlr3")` with `r pkg("mlr3spatial")` and `r pkg("mlr3spatiotempcv")` takes an object-oriented approach with R6 classes for direct spatial object handling and cross-validation within its structured syntax.
-   `r pkg("tidymodels")` with `r pkg("spatialsample")` and `r pkg("waywiser")` introduces spatial sampling strategies and model evaluation tools following tidyverse principles, including `spatialsample::spatial_resample` and `waywiser::ww_area_of_applicability`.
-   `r pkg("RandomForestsGLS")` and `r pkg("spatialRF")` extend Random Forests to incorporate spatial dependence, offering specialized functions for spatial estimation, feature selection, and model assessment.
-   `r pkg("meteo")` implements Random Forest Spatial Interpolation by incorporating nearest observations and distances into the prediction process.
-   `r pkg("gpboost")` captures complex non-linear dependencies by combining gradient boosting with Gaussian processes.
-   `r pkg("sperrorest")` and `r pkg("blockCV")` provide frameworks for spatial resampling and validation, supporting methods like k-means clustering and block-based approaches to account for spatial dependencies in model evaluation.

Installing packages linking to PROJ, GDAL or GEOS
-------------------------------------------------

Installation of packages like `r pkg("gdalraster")`, `r pkg("sf")` and
`r pkg("terra")` which use external software libraries such as PROJ, GDAL or
GEOS requires care. For most users on platforms such as Windows or macOS who
are not themselves package developers, it is always better to avoid what are
known as source installs, because CRAN binary packages include all of the
external software required. Because `getOption("pkgType")` on these platforms
is usually `"both"`, you may be asked to choose to install a source package
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
Windows 4.2 [RTools42](https://cran.r-project.org/bin/windows/Rtools/rtools42/rtools.html), Windows 4.3 [RTools43](https://cran.r-project.org/bin/windows/Rtools/rtools43/rtools.html),
and macOS [both architectures](https://mac.r-project.org/tools/). These 
external software libraries have been built using the same compile and 
link settings as R itself, so avoid the risk of possible errors caused by 
mismatched binaries. The current versions may be updated between R releases, in order to give access to more recent versions of GDAL, GEOS or PROJ.

If you are a user (or developer) on systems where `getOption("pkgType")`
is `"source"`, you will need to ensure that the external software is
available when installing source packages. Advice for some such systems may
be found [here](https://github.com/r-spatial/sf/#installing). 
[The most common reason](https://github.com/r-spatial/sf/#multiple-gdal-geos-andor-proj-versions-on-your-system) 
for failure is having multiple versions of external software installed
on your platform.


### Links
-   [R-SIG-Geo mailing list](https://stat.ethz.ch/mailman/listinfo/R-SIG-Geo/)
-   [r-spatial](https://r-spatial.org/)
-   [rspatial](https://rspatial.org/)
-   [Geocomputation with R](https://geocompx.org/r.html)
