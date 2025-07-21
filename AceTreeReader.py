# -*- coding: utf-8 -*-
"""
Created on Fri Jul 18 14:20:09 2025

@author: SantellA
"""
import networkx as nx
from zipfile import ZipFile
import pandas as pd
import xml.etree.ElementTree as ET

#acetree file reader first implementation in python to my knowlege 2025 Janelia Hackathon
#to avoid dependency on SNIII prototype creates a simple network x graph with all Acetree data as node tags
#note that it faithfully adds fields for EVERY acetree column as node attributes
#Note that this includes those specifying network topology but these are redundant
#with the networkx edges which might lead to confusion if not careful
#resulting networkx graph connectivity reflects only the predecessor ID field, replicating the 
#(undocumented) beehavior of Acetree which ignores redundant successor fields 
#when building its internal representation of the on disk data structure

#note a known bug that some AceTree generated lineages contain an edge to a 
#deleted node. This loader silently corrects any errors of this sort by 
#respecting Acetree semanitcs by refusing to make edges between invalid nodes 
# they seem to correspond to what should be roots of independent subtrees and so are treated as such, but buyer beware
# todo all xml and AuxInfo.csv tags should added as graph level attributes

class AceTreeReader:
    #pasted below is cannonical historical explanation of acetree column meanings
    #1(cell ID specific to this file, not always the same as row number), 1(validity flag: 1 valid, 0 invalid (a cell deleted in acetree)), 
    #24 (predecessor ID in previous file), 45(successor ID in next file), -1(successor 2 id in next file), 
    #178(x), 57(y), 5.2(z), 11(size), ABalpaaaaa(name), 6547(lineage marker 'weight'), 37259(reporter intensity - Boyle Units),
    #6036(summed reporter intensity), 162(voxels), (placeholder for forced Acetree name), 37259(SRI),
    #0(global correction), 0(local correction), 0(blot correction), 0("crosstalk" correction),}
    #below is an example from a real expression containing file  
    #4, 1, -1, 4, -1, 444, 343, 26.0, 25, Nuc26360, 65535, 1523178, 3308342, 2172, , 1523178, 25000, 1522646, 1522646, 0, 
    # dictionary giving the column indicies and data name in each AceTree timepoint file column
    
    # note many of these columns are historical in nature and seldom used in common practice
    #of note 'valid' which if zero indicates the detection was marked as deleted
    #forcedAceTreeName which indicates a node name that has been manually set in Acetree
    #and indictates to Acetree to skip automated naming and use this name going forwad
    #intensity* values may or may not reflect lineaging marker intensity and should not be relied on
    #expIntensity* values will either be all zero, all empty, or reflect expression within detection mask
    # some of these values are obscure in nature, and appear to be arbitrary constants but 
    #in all publications I'm aware of expression is measured as: 
    # expIntensity-expBlotCorrection 
    # this implements the subtraction of a mean background based on a halo computed by acebatch2.jar (Murray at al Nature Methods 2008  doi: 10.1038/nmeth.1228)
    # I think that the units are mean intensity in the mask times an arbitrary constant (1000?) i.e. 'Boyle' (of sainted memory) units

    aceTreeColumnDictionary={'cellIDAtTime':1,'valid':2, 'predecessorID':3, 'successorID1':4,
                            'successorID2':5, 'x':6, 'y':7, 'z':8,'diameter':9, 'name':10,
                            'meanintensity':11, 'intensity':12, 'intensitySummed':13, 'voxelCount':14,
                            'forcedAcetreeName':15, 'expIntensity':16,
                            'expIntensityGlobalCorrection':17,'expIntensityLocalCorreciton':18,
                           'explocalCorrection':19,'expBlotCorrection':20,'expCrosstalkCorrection':21}
    
    #basebath is the root path/name of an acetree file triplet
    #endtime is the (one indexed) endtime 
    #(acetree files of certain vintages have dummy empty nuclei files for later timepoints up to a 'max' and subset by time is common use )
    
    def readFiles(self,basepath,endtime):
        #read metadata
        #look in xml file for image config, location of zip file, resolution
        #should I look in auxinfo or auxinfo2 file for cannonical orientation, right now noplace to put this graph level info?
        #read zip file into network x graph
        mygraph=self.readAceTree(basepath+'.zip',endtime)
        metadata=self.readAceTreeXML(basepath+'.xml')
        #todo check for Auxinfov2.csv with different semantics? very rare use case
        embryometadata=self.readAceTreeCSV(basepath+'AuxInfo.csv')
        # to do add all metadata to graph
        return mygraph
    
    def readAceTreeXML(self,filename):
        tree = ET.parse(filename)
        
        #to do actually parse xml file, which has different tags based on version
        #just pass them back
        return 
    
    def readAceTreeCSV(self, filename):
        #to do read csv
        return None
    #to add 
    def createRowDictionary(self,row):
        rowdictionary={}
        for key in self.aceTreeColumnDictionary:
            rowdictionary[key]=row[self.aceTreeColumnDictionary[key]-1]
        return rowdictionary
    
    #function for initializing detectiondatamanger from AceTree zip file
    def readAceTree(self,file,endtime):
        acetreelineage=nx.Graph()
        addednodes={}#used to keep track of all unique in data nodes IDs and their acetree, unique per frame IDS
        #try:
        c=0 #keep track of index 
        with ZipFile(file) as myzip:   
            for t in range(endtime):
                csvfilename='nuclei/t'+format(t+1, '03')+'-nuclei'
                with myzip.open(csvfilename) as myzipcsv:
                    df = pd.read_csv(myzipcsv,header=None)
                    #create detection for each line csv file and add to detectionDataManager return the DetectionDataManager
                    for row in df.itertuples(index=False):
                            rowdict=self.createRowDictionary(list(row))
                            rowdict['t']=t+1#add t as it is not in acetree file
                            acetreelineage.add_nodes_from([(c,rowdict)])
                            thiskey=str(acetreelineage.nodes[c]['t'])+str(acetreelineage.nodes[c]['cellIDAtTime'])
                            addednodes[thiskey]=c
                            c=c+1
                
        #except Exception as e:
        #    print(f"Error reading with explicit zipfile open: {e}")
        #now that the nodes are all in graph need to add backward link for each non -1 pred, by finding 
        #unique node with t=t-1 and ID predID
        #so as not to search for them I made a list of concatenated t_ID for each cell as made it and keept a list
        for node in acetreelineage.nodes():
            pred=acetreelineage.nodes[node]['predecessorID']
            if not pred==-1:
                #key for t-1 and ID from pred
                thiskey=str(acetreelineage.nodes[node]['t']-1)+str(acetreelineage.nodes[node]['predecessorID'])
                #look up what node was added that matched this time ID signature
                if not (thiskey==None):
                    predc=addednodes[thiskey] #
                    acetreelineage.add_edge(node,predc)
                else:
                    print('should never have a pred reference a nonexistent node ID')
        return acetreelineage