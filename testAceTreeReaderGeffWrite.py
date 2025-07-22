# -*- coding: utf-8 -*-
"""
Created on Sat Jul 19 13:59:25 2025

@author: SantellA

demo of reading an acetree embryo 'file' (zip,xml,csv triad) and writing to geff
"""

#test acetreereader
from AceTreeReader import AceTreeReader
import networkx as nx
from geff import write_nx


myreader=AceTreeReader()
endtime=100
basepath='Z:/bao_data_zrc/baolab/santella_DeepLearning/StarryniteIII/pythonAcetreeSupport/20140407_JIM113_SiO-0.15_1_s3_nobifurcation_trainingversion21_edited'
graph=myreader.readFiles(basepath,endtime)
#nx.draw(graph)

#save as geff
#note some gotchas, time is in frames in acetree, with frame rate not in metadata
#often this is 1 minute (as in test data set) but it is also often 75 seconds, only can tell from microscope aquisition metadata if you have it or experimental notes

# if there is no matching xml file for the basepath,the .zip will be loaded 
#and have x,y,z,radius fields but without pixel sizes none of the _um properties will be added
write_nx(graph, basepath+'.zarr/tracking_graph',
         axis_names=['t','x_um','y_um','z_um'],
         axis_units=['minute','micrometer','micrometer','micrometer'],
         axis_types=['time','space','space','space'])
         #not sure if this is right syntax when  spec updated?
         #sphere=['radius_um']

#write pointing at the uncalibrated x,y,z
# write_nx(graph, basepath+'.zarr/tracking_graph',
#          axis_names=['t','x','y','z'],
#          axis_units=['minute','pixels','pixels','pixels'],
#          axis_types=['time','space','space','space'])
