# -*- coding: utf-8 -*-
"""
Created on Fri Jul 18 14:20:09 2025

@author: SantellA
"""
import networkx as nx
from zipfile import ZipFile
import pandas as pd
import csv
import xml.etree.ElementTree as ET
import warnings

#acetree file reader in python  2025 Janelia Hackathon
#to avoid dependency on SNIII prototype creates a simple network x graph with all Acetree data as node tags
#note that it faithfully adds fields for EVERY acetree column as node attributes
#Note that this includes those specifying network topology but these are redundant
#with the networkx edges which might lead to confusion if not careful
#resulting networkx graph connectivity reflects only the predecessor ID field, replicating the 
#(undocumented) beehavior of Acetree which ignores redundant successor fields 
#when building its internal representation of the on disk data structure

#note a known bug that some AceTree edited lineages contain an edge to a 
#deleted node. This loader corrects any errors of this sort with complaint
#  by refusing to make edges between invalid nodes 
# FYI these errors SEEM to correspond to what should be roots of independent subtrees and so are treated as such, but buyer beware

class AceTreeReader:
    """
      Loader that makes networkX graph from AceTree .zip file
      copies all fields unchanged and if metadata is loaded creates micron versions of all spatial fields
      see code for long commentary on obscure nature of Acetree fields
      """
    #pasted below is the cannonical historical explanation of acetree column meanings
    #1(cell ID specific to this file, often, but not always the same as row number), 1(validity flag: 1 valid, 0 invalid (a cell deleted in acetree)), 
    #24 (predecessor ID in previous file), 45(successor ID in next file), -1(successor 2 id in next file), 
    #178(x), 57(y), 5.2(z, all 3 in potentially anisotropic pixels), 11(size, (diameter)), ABalpaaaaa(name), 6547(lineage marker 'weight'), 37259(reporter intensity - Boyle Units),
    #6036(summed reporter intensity), 162(voxels, in sphere), (empty placeholder for forced Acetree name), 37259(SRI),
    #0(global correction), 0(local correction), 0(blot correction), 0("crosstalk" correction),}
    #below is an example from a real expression containing file  
    #4, 1, -1, 4, -1, 444, 343, 26.0, 25, Nuc26360, 65535, 1523178, 3308342, 2172, , 1523178, 25000, 1522646, 1522646, 0, 
    
    # note many of these columns are historical in nature and seldom used in common current practice
    #of note 'valid' which if zero indicates the detection was marked as deleted
    #forcedAceTreeName which indicates a node name that has been manually set in Acetree
    #and indictates to Acetree to skip automated naming and use this name going forward
    #intensity* values may or may not reflect lineaging marker intensity and should not be relied on
    #expIntensity* values will either be all zero, all empty, or reflect expression within detection masks
    # some of these values are obscure in nature, and appear to be arbitrary constants but 
    #in all publications I'm aware of expression is measured as: 
    # expIntensity-expBlotCorrection 
    # this implements the subtraction of a mean background based on a halo computed by acebatch2.jar (Murray at al Nature Methods 2008  doi: 10.1038/nmeth.1228)
    # I think that the units are mean intensity in the mask times an arbitrary constant (1000?) i.e. 'Boyle' (of sainted memory) units

    #similarly xml file is a mix of used and historical elements
    #the image file is most often historical exemplar single slice image, related to back compatiblity, if not found acetree searches in
    #the directory up for single timepoint per tiff files
    #these files are cropped and flipped before display automatically though newer acetree files have expliit tags
    #that independently indicate if cropping and splitting should be performed 
    #todo might be nice to implement this filename conversion
    #endtime is a max endtime that should be accurate resolution is used but endplane is not
    #blot specifies the expression correction used in displayed values, polar body info is for an obsolete acetree feature. 
    
    #auxinfocsv file (v1 or v2) are inf related to embryo orientation used in naming the 'axis' e.g. ADL for c elegans cannonical orientation is combined with centroid and angle fields to define the AP axis for aligning cannonical division angles to this embryo in naming
    #probably nothing else is used in current practice? v2 specifies 3D angle of embryo more flexibly, but is largely unused.
    # note zpixres is not the z pixel resolution, but the anisotropy factor between xy and z 
    #see https://github.com/zhirongbaolab/StarryNite/tree/master/documentation for more details about current practice
    
    # dictionary giving the column indicies and descriptive data name in each AceTree timepoint file column
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
        """
        Loads a set of AceTree files corresponding to an embryo into a networkx graph preserving all nuclei file columns as node attributes
        loads .xml (which specifies resolution and image location and auxinfo.csv (which specifies centroid and orientation for naming) attaches these as graph attributes
        Args:
            basepath (str): The path and root name of all3 files
            endtime (int): timepoint to load till
    
        Returns:
            networkx.classes.graph.Graph: a graph representing the lineage data in the acetree file with all metadata attched as graph attributes
        """
       
        #read metadata
        #look in xml file for image config, location of zip file, resolution
        #should I look in auxinfo or auxinfo2 file for cannonical orientation, right now noplace to put this graph level info?   
        metadata=self.readAceTreeXML(basepath+'.xml')
        if metadata==None:
            warnings.warn('Warning: No xml file cant find pixel resolution metadata')
        
        #else:
            #to do:
            #if have xml file lets try to dissect tiff name, its likely a slice that does not exist
            #[historical acetree back compatiblity weirdness]
            #take a stab at creating the real image name
          
        #todo check for Auxinfov2.csv with different semantics? very rare use case
        embryometadata=self.readAceTreeCSV(basepath+'AuxInfo.csv')
        if  embryometadata==None:
            warnings.warn('Warning: No auxinfo file looking for auxinfo v2')
            embryometadata=self.readAceTreeCSV(basepath+'AuxInfo_v2.csv')
        
        if  embryometadata==None:
            warnings.warn('Warning: no auxinfo (1 or 2) found, noncanonical embryo may be misnamed')
                          
        #read zip file into network x graph
        mygraph=self.readAceTree(file=basepath+'.zip',endtime=endtime,metadata=metadata)
        mygraph.graph.update(metadata)       
        mygraph.graph.update(embryometadata)
         
        return mygraph
    
    def readAceTreeXML(self,filename):
        """
        loads .xml (which specifies resolution and image location and auxinfo.csv (which specifies centroid and orientation for naming) attaches these as graph attributes
        Args:
           filename  (str): xml file to load
          
        Returns:
            dict: flattened dict of xml elements
        """
        leaf_tags={}
        tree = ET.parse(filename)
        root = tree.getroot()
        #leaf_tags = set()  # Use a set to store unique leaf tags
        for elem in root.iter():
             if not list(elem):  # Check if the element has no children
                 leaf_tags[elem.tag]=elem.attrib
        leaf_tags=dict(self.flatten(leaf_tags))
        return leaf_tags

    
    def flatten(self,obj, prefix=[]):
        if isinstance(obj, str):
            yield ('_'.join(prefix), obj)
    
        elif isinstance(obj, list):
            for o in obj:
                yield from self.flatten(o, prefix) 
        else:
            for k, v in obj.items():
                yield from self.flatten(v, prefix + [k])
                
            
    def readAceTreeCSV(self, filename):
        """
        Converts a CSV file with one data row and headers into a Python dictionary.
    
        Args:
            filepath (str): The path to the CSV file.
    
        Returns:
            dict: A dictionary representing the single data row, with headers as keys.
                  Returns an empty dictionary if the file is empty or only contains headers.
        """
        with open(filename, 'r', newline='') as csvfile:
            reader = csv.DictReader(csvfile)
            for row in reader:
                return row  # Returns the first (and only) data row as a dictionary
        return {} # Return empty dict if no data rows found


    #to add 
    def createRowDictionary(self,row):
        rowdictionary={}
        for key in self.aceTreeColumnDictionary:
            rowdictionary[key]=row[self.aceTreeColumnDictionary[key]-1]
        return rowdictionary
    
    #function for initializing detectiondatamanger from AceTree zip file
    def readAceTree(self,file=None,endtime=None,metadata=None):
        """
        parses Acetree.zip file
    
        Args:
            file (str): The path and name of the .zip file
            endtime (int): the time to load to (1 indexed)
            metadata (dict): dictionary of read xml file values used to compute spatial info in um from acetree file if not missing
    
        Returns:
            Returns:
                networkx.classes.graph.Graph: a graph representing the lineage data in the acetree file
        """
        #add nodes
        acetreelineage=nx.Graph()
        addednodes={}#used to keep track of all unique in data nodes IDs and their acetree, unique per frame IDS
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
                            rowdict['radius']=rowdict['diameter']/2.0
                            if not metadata==None:
                                #add pixel size calibrated spatial fields
                                rowdict['x_um']=rowdict['x']*float(metadata['resolution_xyRes'])
                                rowdict['y_um']=rowdict['y']*float(metadata['resolution_xyRes'])
                                rowdict['z_um']=rowdict['z']*float(metadata['resolution_zRes'])
                                rowdict['radius_um']=rowdict['radius']*float(metadata['resolution_xyRes'])
                                
                            acetreelineage.add_nodes_from([(c,rowdict)])
                            thiskey=str(acetreelineage.nodes[c]['t'])+str(acetreelineage.nodes[c]['cellIDAtTime'])
                            addednodes[thiskey]=c
                            c=c+1
                
        #now that the nodes are all in graph need to add backward link for each non -1 pred, by finding 
        #unique node with t=t-1 and ID predID
        #so as not to search for them I made a list of concatenated t_ID for each cell as made it and keept a list
        for node in acetreelineage.nodes():
            pred=acetreelineage.nodes[node]['predecessorID']
            if not pred==-1: #special no link code
           
                #key for t-1 and ID from pred
                thiskey=str(acetreelineage.nodes[node]['t']-1)+str(acetreelineage.nodes[node]['predecessorID'])
                #look up what node was added that matched this time ID signature
                if not (thiskey==None):
                    #check if weird bug where edges link to deleted nodes just drop the edge from networkX if so
                    predc=addednodes[thiskey] #
                    isvalidedge= acetreelineage.nodes[node]['valid']==1 and acetreelineage.nodes[predc]['valid']==1
                    if isvalidedge:
                        acetreelineage.add_edge(node,predc)
                    else:
                        warnings.warn('Warning: valid node pointing to invalid, this shouldnt happen but does on occasion-a longstanding mystery.')
                else:
                    raise Exception('should never have a pred reference a nonexistent node ID')
        return acetreelineage