# ***calsimpy***: Functions to help with common CalSim data processing and analysis tasks
## *CalSim can be hard, maybe Python can make some things a little easier*

This package includes a number of functions to work with CalSim (CalSimII, CalSim3, CalLite)
data in an organized way within Python. 
The code is primarily built around the concept of creating an object for each CalSim study
that the modeler would like to work with. 
Some of the key functions associated with these CalSim study objects are related to
reading and writing DSS files (requires a separate Python-based DSS wrapper - the `calsimpy` 
package is currently built around the `yapydss` package [https://github.com/james-m-gilbert/yapydss]()), but components to read parts of the CalSim3 modeling chain (e.g. CalSimHydro input files and CalSim3/C2VSim groundwater model files) and some plotting functions have been developed as well.

Much of the  - not to mention that all of this needs some documentation!
as the need arose in my day-to-day work over the last ~6-7 years.
There are many other things I'd like to add, and probably a lot of the code that could use
some cleanup, reorganization, and certaintly documentation. 
For that reason, this package should be considered very much in active development.
I'm making it available in case others may find it useful, but cannot promise it will exactly 
meet any one CalSim modeler's needs.

## Installation

TODO


## Concepts & Basic Workflow


## Function Documentation

TODO
