MODSFabTools
============

Software tools written in MATLAB to help with testing, data analysis, and to direct fabrication efforts

###LinescanAlignment.m
This code correspondes to a GUI that allows the user to compare, manipulate, and align mechanically acquired surface data, typically from a coordinate measuring machine (CMM). All calculations are done as code within each Form Object.

The GUI looks like this:
![LinescanAlignment GUI](/src/linescanalignmentgui.png?raw=true "LinescanAlignment GUI")

###SurfaceInterpolation.m
This code reads in a series of linescans aligned with LinescanAlignment.m and outputs an interpolated colormap of the surface. Fine adjustments are made to each component linescan to make the surface contiguous and as smooth as possible.

Sample output looks like this:
![Surface Interpolation Output](/src/interpolated.jpg?raw=true "Surface Interpolation Output")

###StrokeMaker.m
Using normalized units of removal rate for a given set of grinding/polishing tools, this creates a predicted average radial removal profile and exports required dwell times for each tool and radial location, along with stroke parameters (speed, width, etc.)

There is no preview of this output.
