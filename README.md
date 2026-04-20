# ***calsimpy***: Functions to help with common CalSim data processing and analysis tasks
## *CalSim can be hard. This is my attempt to make (some) things easier*

This package includes a number of functions to work with CalSim (CalSimII, CalSim3, CalLite)
data in an organized way using Python. 
The code is primarily built around the concept of creating an object for each CalSim study
that the modeler would like to work with. 


Some of the key functions associated with these CalSim study objects are related to
reading and writing DSS files (requires a separate Python-based DSS wrapper - the `calsimpy` 
package is currently built around the `yapydss` package [https://github.com/james-m-gilbert/yapydss]()), but components to read parts of the CalSim3 modeling chain (e.g. CalSimHydro input files and CalSim3/C2VSim groundwater model files) and some plotting functions have been developed as well.
Where possible, CalSim data is handled in common Python data structures (lists, dictionaries, Pandas DataFrames, and GeoDataFrames) so that any CalSim-derived data can be written out to files that can be used later or viewed with other software (e.g., csv, GIS shapefile, etc)

Much of what is in this package has been pieced together in response to specific needs I have encountered in my work with CalSim over the last seven-ish years (and counting). 
This code has never been the _main_ project, but rather something I built to help me get to where I needed to go.
It's not a polished final product - there are many other things I'd like to add, and a lot of the code that could use
some cleanup and reorganization.
And yes - all of this needs some documentation!

**For that reason, this package should be considered very much in active development.
I'm making it available in case others may find it useful, but cannot promise it will exactly 
meet any one CalSim modeler's needs or work for every task.**

**If you find this package, test it out, and see areas for improvement or expansion, let me know 
(james.gilbert [at] ucsc.edu) - I'm always open to constructive feedback and potential collaboration!**

## Installation

A more detailed set of installation instructions are coming soon.
If you're brave and want to try installing the package in it's current state, you can `pip` to install directly from GitHub using:
`pip install git+https://github.com/james-m-gilbert/calsimpy.git`

Note - there may still be some bugs related to dependencies. I am working on a more streamlined package and hope to have updates soon.


## Concepts & Basic Workflow
Coming soon - check back for updates!


## Function Documentation

Coming soon - check back for updates!
