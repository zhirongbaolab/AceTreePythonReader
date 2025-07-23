# AceTreePythonReader
Acetree file format reader for python.  Loads Acetree .zip file of csv files into a networkX graph, preserving all fields and adding spatial fields in um when it finds xml file with pixel sizes.

Expects a triad of .zip, .xml and Auxinfo.csv files and aggregates all the data in them into one network x graph.
