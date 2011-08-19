#!/usr/bin/env python
import os
import Numeric
#module for triangle package
import triangulate
import triIfaceUtils
import triIfaceFileUtils

__doc__ = """ collect functionality to make a lightweight interface
for Shewchuk's Triangle package. Follow the ellipt2d example in some
respects

TODO: 

"""

class TriangleBaseMesh:
    """ This is basically a wrapper for the triangulateio interface
    that should be able to create a triangle mesh in different ways

       from .ele and .node files
       from a .poly file
       from a adhPyUtils mesh

    It should also be able to generate a adhPyUtils mesh from the
    triangle representation and allow the user to access the basic
    data arrays in triangle

    """

    def __init__(self,baseFlags="zen",nbase=0,verbose=0):
        """
        initialize the triangulateio object,
        keep track of what numbering scheme it uses,
        create a base set of flags for triangulate (e.g., z if using base 0)

        does not create a Voronoi diagram by default
        """
        self.trirep = []
        self.trirep.append(triangulate.new())
        self.trirepDefaultInit = []
        self.trirepDefaultInit.append(True)
        self.vorrep = []
        self.vorrep.append(triangulate.new()) #don't create a Voronoi diagram by default
        self.makeVoronoi = False
        self.vorrepDefaultInit = []
        self.vorrepDefaultInit.append(True)
        assert(nbase in [0,1])
        self.nbase  = nbase 
        self.baseFlags = baseFlags
        if self.nbase == 0 and self.baseFlags.find('z') == -1:
            print """WARNING TriangleMesh base numbering= %d
            reqires z in baseFlags = %s, adding """ % (self.nbase,self.baseFlags)
            self.baseFlags += "z"
        #end if
        if self.baseFlags.find('v') >= 0:
            self.makeVoronoi = True
        if verbose > 0:
            print """TriangleBaseMesh nbase=%d baseFlags= %s """ % (self.nbase,
                                                                    self.baseFlags)
        #end if
        #keep track of last set of flags used to generate rep
        self.lastFlags = None
        #to make printing out info easier
        self.infoFormat = """triangulation:
number of points        = %d
attributes per point    = %d
number of triangles     = %d
nodes per triangle      = %d
attributes per triangle = %d
number of segments      = %d
number of holes         = %d
number of regions       = %d
number of edges         = %d
        """
            
        #end init
    def resetDefaultTrirep(self,index=0,verbose=0):
        """
        reset the Triangle mesh representation if it has been set to something nontrivial
        
        """
        assert(index in range(0,len(self.trirep)))
        #delete the existing representation if it has been initialized
        if not self.trirepDefaultInit[index]:
            if verbose > -1:
                print 'TriangleIface deleting current trirep[',index,']'
            #end verbose
            del self.trirep[index]
            #get a new represenation
            self.trirep[index] = triangulate.new()
            selt.trirepDefaultInit[index] = True
        #end trirep initialized to something non trivial

    def resetDefaultVorrep(self,index=0,verbose=0):
        """
        reset the Triangle mesh Voronoi representation if it has
        been set to something nontrivial
        """
        assert(index in range(0,len(self.vorrep)))
        #delete the existing representation if it has been initialized
        if not self.vorrepDefaultInit[index]:
            if verbose > -1:
                print 'TriangleIface deleting current vorrep[',index,']'
            #end verbose
            del self.vorrep[index]
            #get a new represenation
            self.vorrep[index] = triangulate.new()
            selt.vorrepDefaultInit[index] = True
        #end trirep initialized to something non trivial

    def convertToAdhPyUtilMesh(self,verbose=0):
        """
        Generate a representation in the format expected by
        adhPyUtil.

        Need to make sure the triangle mesh has generated
           nodes
           elements (triangles)
           edges
           neighbors
           
        First set the _global arrays for
            nodes
            elements
        """
        import MeshTools
        triInfo = triangulate.getInfo(self.trirep[0])
        if verbose > 1:
            print 'generating adhPyUtilMesh:'
            print self.infoFormat % triInfo
        #end if

        #create basic adhPyUtil mesh
        meshout = MeshTools.Mesh()

        #get basic information to make sure I can proceed
        #with current mesh representation
        nNodes_global = triInfo[0]; nElems_global = triInfo[2];
        nNodes_elem   = triInfo[3]; nEdges_global = triInfo[8];
        spaceDim      = 2
        assert(nNodes_global > 0 and nElems_global > 0
               and nNodes_elem >= 3 and nEdges_global > 0)

        #subtract off base since adhPyMesh wants base 0 more or less
        nbase = self.nbase
        #should also check? base zero, 
        #get the minimum array information
        ##nodes
        nodeArray = triangulate.getPoints(self.trirep[0])
        meshout.nNodes_global     = nNodes_global
        meshout.nodeArray = Numeric.zeros((meshout.nNodes_global,3),
                                          Numeric.Float)
        for nN in range(nNodes_global):
            meshout.nodeArray[nN,:spaceDim] = nodeArray[nN,:spaceDim]-nbase
        ##elements
        elemArray  = triangulate.getTriangles(self.trirep[0])
        #ignore higher order nodes for adhPyUtil mesh
        meshout.nNodes_element    = spaceDim+1
        meshout.nElements_global  = nElems_global
        #just copy over first 3 nodes
        meshout.elementNodesArray = Numeric.zeros((nElems_global,spaceDim+1),
                                                  Numeric.Int)
        for eN in range(nElems_global):
            meshout.elementNodesArray[eN,:] = elemArray[eN,:spaceDim+1]-nbase
        #end eN

        ##adhPyMesh keeps elements per node
        nodeElementsDict={}
        for eN in range(nElems_global):
            for nN_element in range(spaceDim+1):
                nN = meshout.elementNodesArray[eN,nN_element]
                if  nodeElementsDict.has_key(nN):
                    nodeElementsDict[nN].append(eN)
                else:
                    nodeElementsDict[nN] = [eN]
                #end if
            #end for nN_element
        #end eN
        meshout.max_nElements_node = max(len(nodeElementsDict[nN]) for nN in range(meshout.nNodes_global))

        meshout.nElements_node    = Numeric.zeros((meshout.nNodes_global,),Numeric.Int)
        meshout.nodeElementsArray = Numeric.zeros((meshout.nNodes_global,
                                                   meshout.max_nElements_node),
                                                  Numeric.Int)
        for nN,elementList in nodeElementsDict.iteritems():
            meshout.nElements_node[nN] = len(elementList)
            for eN_element,eN in enumerate(elementList):
                meshout.nodeElementsArray[nN,eN_element]=eN
            #end eN_element
        #end nN

        ##now build the element <--> element boundary arrays
        meshout.nElementBoundaries_element = spaceDim+1
        #make sure Triangle keeps all edges and not just boundary ones 
        meshout.nElementBoundaries_global  = nEdges_global 

        #maps element, local edge number to global edge number
        meshout.elementBoundariesArray = Numeric.zeros((nElems_global,spaceDim+1),
                                                       Numeric.Int)

        #maps edge to global element on left and right (out of domain is 0)
        meshout.elementBoundaryElementsArray=Numeric.ones(
            (meshout.nElementBoundaries_global,2),Numeric.Int)
        meshout.elementBoundaryElementsArray*=-1
        #maps global edge to its local number on the neighboring elements
        meshout.elementBoundaryLocalElementBoundariesArray=Numeric.zeros(
            (meshout.nElementBoundaries_global,2),Numeric.Int)

        #several options for edge to "left" and "right" element neighbor,
        #could number each neighbor according to which one
        #has a given edge first in the current numbering
        #could try to make element 0 be the one that has same orientation of edge

        elementBoundaryElementsCardArray = Numeric.zeros((meshout.nElementBoundaries_global,),
                                                         Numeric.Int) #CardArray in Mesh
        #I have to generate the element-->global edge number mapping manually
        edgeArray = triangulate.getEdges(self.trirep[0])
        edgeDict  = {}
        for edgeN in range(nEdges_global):
            n0 = edgeArray[edgeN,0]-nbase; n1 = edgeArray[edgeN,1]-nbase
            edgeDict[(n0,n1)] = edgeN #global edge number
        #end for
        
        for eN in range(nElems_global):
            locNodes = elemArray[eN,:spaceDim+1]-nbase #global node numbers
            edgeOrientedSame = [False,False,False]
            #edge I is across from node I
            #0
            e0 = (locNodes[1],locNodes[2]); e0rev = (locNodes[2],locNodes[1])
            e0_global = None
            if edgeDict.has_key(e0):
                e0_global = edgeDict[e0] #same orientation
                edgeOrientedSame[0] = True
            elif edgeDict.has_key(e0rev):     #opposite orientation
                e0_global = edgeDict[e0rev]   #could make negative to denote orientation
            #end if
            assert(not e0_global == None)
            #1
            e1 = (locNodes[2],locNodes[0]); e1rev = (locNodes[0],locNodes[2])
            e1_global = None
            if edgeDict.has_key(e1):
                e1_global = edgeDict[e1] #same orientation
                edgeOrientedSame[1] = True
            elif edgeDict.has_key(e1rev):     #opposite orientation
                e1_global = edgeDict[e1rev]   #could make negative to denote orientation
            #end if
            assert(not e1_global == None)
            #2
            e2 = (locNodes[0],locNodes[1]); e2rev = (locNodes[1],locNodes[0])
            e2_global = None
            if edgeDict.has_key(e2):
                e2_global = edgeDict[e2] #same orientation
                edgeOrientedSame[2] = True
            elif edgeDict.has_key(e2rev):     #opposite orientation
                e2_global = edgeDict[e2rev]   #could make negative to denote orientation
            #end if
            assert(not e2_global == None)
            eI_global = Numeric.array([e0_global,e1_global,e2_global],
                                      Numeric.Int)
            meshout.elementBoundariesArray[eN,:] = eI_global
            for eI in eI_global:
                #edge --> element mappings
                elementBoundaryElementsCardArray[eI] += 1
            #end eI

            
            for eiloc in range(meshout.nElementBoundaries_element):
                eI     = eI_global[eiloc]
                #first visited labelling
                elneig = elementBoundaryElementsCardArray[eI]-1
                #elneig = 0 #same orientation
                #if not edgeOrientedSame[eiloc]
                #elneig = 1 #opposite
                ##endif
                meshout.elementBoundaryElementsArray[eI,elneig]=eN
                meshout.elementBoundaryLocalElementBoundariesArray[eI,elneig]=eiloc
            #end eiloc
        #end eN

        #perform some sanity checks
        for edgeN in range(nEdges_global):
            assert(elementBoundaryElementsCardArray[edgeN] in [1,2])
            eN0 = meshout.elementBoundaryElementsArray[edgeN,0]
            eN1 = meshout.elementBoundaryElementsArray[edgeN,1]
            assert(0 <= eN0 and eN0 < meshout.nElements_global)
            if elementBoundaryElementsCardArray[edgeN] == 1:
                assert(eN1 == -1)
            else:
                assert(0 <= eN1 and eN1 < meshout.nElements_global)
            #end if
        #end sanity check on edges
        #now take care of boundary edges
        #interior elements counted twice
        sumCard = Numeric.sum(elementBoundaryElementsCardArray)
        nExtBnd = 2*meshout.nElementBoundaries_global-sumCard
        nIntBnd = meshout.nElementBoundaries_global-nExtBnd
        meshout.nExteriorElementBoundaries_global=nExtBnd
        meshout.nInteriorElementBoundaries_global=nIntBnd

        #global edge numbers on exterior boundary
        meshout.exteriorElementBoundariesArray=Numeric.zeros((nExtBnd,),Numeric.Int)
        meshout.interiorElementBoundariesArray=Numeric.zeros((nIntBnd,),Numeric.Int)

        #enumerate interior and exterior
        interiorI = 0
        exteriorI = 0
        for ebN in range(meshout.nElementBoundaries_global):
            if elementBoundaryElementsCardArray[ebN] == 1:
                meshout.exteriorElementBoundariesArray[exteriorI]=ebN
                exteriorI+= 1
            else:
                meshout.interiorElementBoundariesArray[interiorI]=ebN
                interiorI+= 1
            #end if on card
        #end ebN

        #now track which nodes are on the boundary
        #maps edge --> its nodes
        #this is just the edgelist in triangle
        meshout.elementBoundaryNodesArray = Numeric.zeros((nEdges_global,spaceDim),
                                                          Numeric.Int)
        for edgeN in range(nEdges_global):
            meshout.elementBoundaryNodesArray[edgeN,:spaceDim]=edgeArray[edgeN,:spaceDim]-nbase
        #end edgeN
        #2d so edgeList is same as elementBoundaryNodesArray
        meshout.edgeNodesArray = Numeric.zeros((nEdges_global,spaceDim),
                                               Numeric.Int)
        for edgeN in range(nEdges_global):
            meshout.edgeNodesArray[edgeN,:spaceDim]=edgeArray[edgeN,:spaceDim]-nbase
        #end edgeN

        #now tell mesh that it doesn't have the list interface
        meshout.hasListInterface=False
        #compute diameters array manually
        meshout.elementDiametersArray=Numeric.zeros((meshout.nElements_global,),
                                                    Numeric.Float)
        import math
        for eN in range(meshout.nElements_global):
            elen = Numeric.zeros((meshout.nElementBoundaries_element,),
                                 Numeric.Float)
            for eloc in range(meshout.nElementBoundaries_element):
                eg = meshout.elementBoundariesArray[eN,eloc] #glocal edge number
                n0,n1 = meshout.elementBoundaryNodesArray[eg,0:2] #global node numbers number
                de = meshout.nodeArray[n0,:]-meshout.nodeArray[n1,:]
                elen[eloc] = math.sqrt(de[0]**2 + de[1]**2)
            #end eloc
            meshout.elementDiametersArray[eN]=max(elen)
        #end eN
        meshout.hasGeometricInfo=True
        return meshout
    #end convertToAdhPyUtilMesh
    def convertFromAdhPyUtilMesh(self,meshin,verbose=0):
        """
        generate a Triangle mesh representation from an adhPyUtil mesh.
        This version will copy over the nodes and elements from the
        adhPyUtil mesh.

        Deletes the existing Triangle mesh and regenerates.
        """

        #first check that the input mesh has the correct dimensions and
        #minimal information necessary
        spaceDim = 2
        assert(meshin.nElementBoundaries_element == spaceDim+1)
        assert(not meshin.nodeArray == None)
        assert(not meshin.elementNodesArray == None)

        #get a clean slate
        tri0 = triangulate.new()
        #don't set any markers for now
        #input array should be nNodes by spacedim
        nodesIn = meshin.nodeArray[:,:spaceDim]
        triangulate.setPoints(tri0,nodesIn)
        #triangle array should be nElements x 3
        triAin  = meshin.elementNodesArray[:,:spaceDim+1]
        triangulate.setTriangles(tri0,triAin)

        flagsAdd = "r" #use element and node connections
        flags = flagsAdd+self.baseFlags
        if flags.find('v') >= 0:
            self.makeVoronoi = True
            if verbose > 0:
                print 'readNodeAndEle makeVoronoi= ',self.makeVoronoi,' flags= ',flags
        #end if
        self.lastFlags = flags
        #reset data representations if necessary
        self.resetDefaultTrirep(index=0,verbose=verbose)
        if self.makeVoronoi:
            self.resetDefaultVorrep(index=0,verbose=verbose)
        #end vorrep initialized
        
        if self.makeVoronoi:
            triangulate.applyTriangulate(flags,tri0,self.trirep[0],self.vorrep[0])
            #handles no longer contain trivial representations
            self.trirepDefaultInit[0] = False
            self.vorrepDefaultInit[0] = False
            
        else:
            triangulate.applyTriangulateNoVoronoi(flags,tri0,self.trirep[0])
            #handle no longer contains trivial representations
            self.trirepDefaultInit[0] = False
        #end if

        #clean up explicitly
        del tri0
    #end convertFromAdhPyUtilMesh
    ##################################################
    #read from file routines
    ##################################################
    def readFromNodeFile(self,filebase,flagsAdd="",verbose=0):
        """
        read triangle representation from filebase.node
        files. assumes the nbase data member is set appropriately
        """
        reader = triIfaceUtils.TriangleInputFileReader()
        if verbose > 0:
            print 'readFromNodeAndEleFiles filebase= ',filebase
        #end if
        #could specify comment character too
        nodeDataInfo,nodeData = reader.readNodes(filebase,nbase=self.nbase)
        nodes,nodesA,nodesM = (nodeData['nodes'],
                               nodeData['nodeAttributes'],
                               nodeData['nodeMarkers'])

        if verbose > 3:
            print 'Nodes: nodeInfo= ',nodeDataInfo
            print """Nodes: nodes= \n%s\n nodeAttributes= \n%s\n nodeMarkers= \n%s
            """ % (nodes,nodesA,nodesM)
        #end if

        #now create an initial representation
        tri0 = triangulate.new()

        if nodesM == None:
            triangulate.setPoints(tri0,nodes)
        else:
            triangulate.setPointsAndMarkers(tri0,nodes,nodesM)
        if not nodesA == None: 
            triangulate.setPointAttributes(tri0,nodesA)
        
        #run triangulate on it using the base flags and whatever else was
        #added
        flags = flagsAdd+self.baseFlags
        if flags.find('v') >= 0:
            self.makeVoronoi = True
            if verbose > 0:
                print 'readNodeAndEle makeVoronoi= ',self.makeVoronoi,' flags= ',flags
        #end if
        self.lastFlags = flags

        #reset data representations if necessary
        self.resetDefaultTrirep(index=0,verbose=verbose)
        if self.makeVoronoi:
            self.resetDefaultVorrep(index=0,verbose=verbose)
        #end vorrep initialized


        if self.makeVoronoi:
            triangulate.applyTriangulate(flags,tri0,self.trirep[0],self.vorrep[0])
            #handles no longer contain trivial representations
            self.trirepDefaultInit[0] = False
            self.vorrepDefaultInit[0] = False
        else:
            triangulate.applyTriangulateNoVoronoi(flags,tri0,self.trirep[0])
            #handle no longer contains trivial representations
            self.trirepDefaultInit[0] = False
        #end if

        #do I need to clean up explicitly?
        del tri0
    def readFromNodeAndEleFiles(self,filebase,flagsAdd="",verbose=0):
        """
        read triangle representation from filebase.node and filebase.ele
        files. assumes the nbase data member is set appropriately
        """
        reader = triIfaceUtils.TriangleInputFileReader()
        if verbose > 0:
            print 'readFromNodeAndEleFiles filebase= ',filebase
        #end if
        nodeDataInfo,nodeData = reader.readNodes(filebase,nbase=self.nbase)
        nodes,nodesA,nodesM = (nodeData['nodes'],
                               nodeData['nodeAttributes'],
                               nodeData['nodeMarkers'])

        if verbose > 3:
            print 'Nodes: nodeInfo= ',nodeDataInfo
            print """Nodes: nodes= \n%s\n nodeAttributes= \n%s\n nodeMarkers= \n%s
            """ % (nodes,nodesA,nodesM)
        #end if

        triDataInfo,triData = reader.readTriangles(filebase,nbase=self.nbase)
        triangles,trianglesA = triData['triangles'],triData['triangleAttributes']

        if verbose > 3:
            print 'Triangles: triInfo= ',triDataInfo
            print """Triangles: elems= \n%s\n triAttributes= \n%s
            """ % (triangles,trianglesA)
        #end if

        #now create an initial representation
        tri0 = triangulate.new()
        
        if nodesM == None:
            triangulate.setPoints(tri0,nodes)
        else:
            triangulate.setPointsAndMarkers(tri0,nodes,nodesM)
        if not nodesA == None: 
            triangulate.setPointAttributes(tri0,nodesA)
        #end if
        
        triangulate.setTriangles(tri0,triangles)
        if not trianglesA == None: 
            triangulate.setTriangleAttributes(tri0,trianglesA)

        #run triangulate on it using the base flags and whatever else was
        #added
        flags = flagsAdd+self.baseFlags
        if flags.find('v') >= 0:
            self.makeVoronoi = True
            if verbose > 0:
                print 'readNodeAndEle makeVoronoi= ',self.makeVoronoi,' flags= ',flags
        #end if
        self.lastFlags = flags

        #reset data representations if necessary
        self.resetDefaultTrirep(index=0,verbose=verbose)
        if self.makeVoronoi:
            self.resetDefaultVorrep(index=0,verbose=verbose)
        #end vorrep initialized

        if self.makeVoronoi:
            triangulate.applyTriangulate(flags,tri0,self.trirep[0],self.vorrep[0])
            #handles no longer contain trivial representations
            self.trirepDefaultInit[0] = False
            self.vorrepDefaultInit[0] = False
        else:
            triangulate.applyTriangulateNoVoronoi(flags,tri0,self.trirep[0])
            #handle no longer contains trivial representations
            self.trirepDefaultInit[0] = False
        #end if

        #do I need to clean up explicitly?
        del tri0
    def readFromPolyFile(self,filebase,flagsAdd="",verbose=0):
        """
        read triangle representation from filebase.poly file
        assumes the nbase data member is set appropriately
        """
        reader = triIfaceUtils.TriangleInputFileReader()
        if verbose > 0:
            print 'readFromPolyFile filebase= ',filebase
        #end if

        polyDataInfo,polyData = reader.readPoly(filebase,nbase=self.nbase,
                                                verbose=verbose)

        if verbose > 3:    
            for type in ['node','segment','hole','region']:
                print 'Poly: ',type,'Info= \n',polyDataInfo[type],'\n'
            #end for
        #end if
            
        nodes,nodesA,nodesM = (polyData['node']['nodes'],
                               polyData['node']['nodeAttributes'],
                               polyData['node']['nodeMarkers'])
        segments,segmentsM   = (polyData['segment']['segments'],
                                polyData['segment']['segmentMarkers'])
        holes                 = polyData['hole']['holes']
        regions               = polyData['region']['regions']

        if verbose > 3:
            print """Poly file read:
            nodes   = \n%s\n nodeAttributes= \n%s\n nodeMarkers= \n%s\n 
            segments= \n%s\n segmentMarkers= \n%s\n 
            holes   = \n%s\n
            regions = \n%s\n
            """ % (nodes,nodesA,nodesM,segments,segmentsM,holes,regions)
        #end if
        tri0 = triangulate.new()
    
        if nodesM == None:
            triangulate.setPoints(tri0,nodes)
        else:
            triangulate.setPointsAndMarkers(tri0,nodes,nodesM)
        if not nodesA == None: 
            triangulate.setPointAttributes(tri0,nodesA)
        #end if
        if segmentsM == None:
            triangulate.setSegments(tri0,segments)
        else:
            triangulate.setSegmentsAndMarkers(tri0,segments,segmentsM)
        #end if
        if (not holes == None):
            triangulate.setHoles(tri0,holes)
        #end if
        if (not regions == None):
            #print 'setting trin1 regions=\n',regions2
            triangulate.setRegions(tri0,regions)
        #end if

        flags = flagsAdd+self.baseFlags
        if flags.find('p') == -1:
            print 'flags = ',flags,' must have p, appending'
            flags += "p"
        #end
        if flags.find('v') >= 0:
            self.makeVoronoi = True
            if verbose > 0:
                print 'readNodeAndEle makeVoronoi= ',self.makeVoronoi,' flags= ',flags
        #end if
        self.lastFlags = flags
        #reset data representations if necessary
        self.resetDefaultTrirep(index=0,verbose=verbose)
        if self.makeVoronoi:
            self.resetDefaultVorrep(index=0,verbose=verbose)
        #end vorrep initialized

        if self.makeVoronoi:
            triangulate.applyTriangulate(flags,tri0,self.trirep[0],self.vorrep[0])
            #handles no longer contain trivial representations
            self.trirepDefaultInit[0] = False
            self.vorrepDefaultInit[0] = False
        else:
            triangulate.applyTriangulateNoVoronoi(flags,tri0,self.trirep[0])
            #handle no longer contains trivial representations
            self.trirepDefaultInit[0] = False
        #end if

        #clean up?
        del tri0
    ##################################################
    #output routines
    ##################################################
    def writeToFile(self,filebase,verbose=0):
        """
        Just write out basic files for triangulateion
        Still need to write out Voronoi diagram 
        """
        triIfaceUtils.writeOutTriangulation(self.trirep[0],filebase,
                                            self.nbase,verbose)

    #end def

    def viewShowme(self):
        """
        just call showme for the mesh held in rep, uses meshshow-tmp file
        """
        filebase="meshshow-tmp"
        self.writeToFile(filebase)
        globshowme= triIfaceUtils.showmeCmdBase
        showmecmd = """%s  %s """ % (globshowme,filebase)
        
        failure = 0
        failure = triIfaceFileUtils.checkFileExists(globshowme)
        failure = os.system(showmecmd)

        return failure
    
    #end viewShowme
#end TriangleBaseMesh

########################################################################
#define some simple helper functions and examples on calling the code
########################################################################

def triIfaceCall3(filebase="trimesh",baseFlags="zen",
                  flagsAdd="",verbose=0):
    import TriangleIface
    nbase = 0
    if baseFlags.find('z') == -1:
        nbase=1
    mesh = TriangleIface.TriangleBaseMesh(baseFlags=baseFlags,
                                          nbase=nbase,
                                          verbose=verbose)

    if flagsAdd.find('p') >= 0:
        if verbose > 0:
            print """flagsAdd= %s trying to read %s """ % (flagsAdd,
                                                           filebase+'.poly')
        #end verbose
        mesh.readFromPolyFile(filebase,flagsAdd,verbose=verbose)

    elif flagsAdd.find('r') >= 0:
        if verbose > 0:
            print """flagsAdd= %s trying to read %s and %s """ % (flagsAdd,
                                                                  filebase+'.node',
                                                                  filebase+'.ele')
        #end verbose
        mesh.readFromNodeAndEleFiles(filebase,flagsAdd,verbose=verbose)
        
    else:
        if verbose > 0:
            print """flagsAdd= %s trying to read %s """ % (flagsAdd,
                                                           filebase+'.node')
        #end verbose
        
        mesh.readFromNodeFile(filebase,flagsAdd,verbose=verbose)
    #end if on flags

    print 'viewing mesh generated'
    mesh.viewShowme()

#end triIfaceCall3
def exAdhPyUtilMesh0(filebase="trimesh",baseFlags="zen",
                     flagsAdd="",viewMesh=1,verbose=0):
    """
    create a Triangle mesh representation
    read it in from a file and initialize
    convert to an adhPyUtilMesh
    """
    import TriangleIface
    nbase = 0
    if baseFlags.find('z') == -1:
        nbase=1
    mesh = TriangleIface.TriangleBaseMesh(baseFlags=baseFlags,
                                          nbase=nbase,
                                          verbose=verbose)

    if flagsAdd.find('p') >= 0:
        if verbose > 0:
            print """flagsAdd= %s trying to read %s """ % (flagsAdd,
                                                           filebase+'.poly')
        #end verbose
        mesh.readFromPolyFile(filebase,flagsAdd,verbose=verbose)

    elif flagsAdd.find('r') >= 0:
        if verbose > 0:
            print """flagsAdd= %s trying to read %s and %s """ % (flagsAdd,
                                                                  filebase+'.node',
                                                                  filebase+'.ele')
        #end verbose
        mesh.readFromNodeAndEleFiles(filebase,flagsAdd,verbose=verbose)
        
    else:
        if verbose > 0:
            print """flagsAdd= %s trying to read %s """ % (flagsAdd,
                                                           filebase+'.node')
        #end verbose
        
        mesh.readFromNodeFile(filebase,flagsAdd,verbose=verbose)
    #end if on flags
    if viewMesh > 0:
        print 'viewing mesh generated from .node and .ele files'
        mesh.viewShowme()
    #end if
    adhPyMesh = mesh.convertToAdhPyUtilMesh(verbose)
    adhPyMesh.writeEdgesGnuplot2("gnuMesh") #uses array interface
    if viewMesh > 0:
        adhPyMesh.viewMeshGnuplotPipe("gnuMesh")
    matmesh = "matlabMesh"
    adhPyMesh.buildMatlabMeshDataStructures(matmesh)
    
#end exAdh
def exAdhPyUtilMesh1(baseFlags="zen",viewMesh=1,verbose=0):
    """
    create an adhPyUtilMesh for a rectangle
    create a Triangle mesh representation
    convert the adhPyUtilMesh mesh to the Triangle mesh representation
    """
    import MeshTools
    import TriangleIface
    #simple domain for now
    Lx = 2.0
    Ly = 1.0
    nx = 11
    ny = 6
    adhPyMesh = MeshTools.TriangularMesh()
    adhPyMesh.constructTriangularMeshOnRectangle(Lx,Ly,nx,ny)
    adhPyMesh.writeEdgesGnuplot2("gnuMesh") #uses array interface
    if viewMesh > 0:
        adhPyMesh.viewMeshGnuplotPipe("gnuMesh")
    
    nbase = 0
    if baseFlags.find('z') == -1:
        nbase=1
    #end if
    trimesh = TriangleIface.TriangleBaseMesh(baseFlags=baseFlags,
                                             nbase=nbase,
                                             verbose=verbose)
    trimesh.convertFromAdhPyUtilMesh(adhPyMesh,verbose=verbose)

    if viewMesh > 0:
        print 'viewing mesh generated from adhPyUtil mesh'
        trimesh.viewShowme()
    #end if

#end exAdh
if __name__ == '__main__':
    #make sure python has adhPyUtils in paty
    import os,sys
    from optparse import OptionParser
    parser = OptionParser()

    #options controlling simulation behavior
    parser.add_option('-P','--adhPyDir',
                      default='/Users/mfarthin/Public/code/chrisAdhPyUtil/cek-dev/',
                      help="""where to find adhPyUtil library""")
    parser.add_option('--baseFlags',
                      default='zen',
                      help="""base flags for creating Triangle meshes""")
    parser.add_option('--filebase',
                      default='trimesh',
                      help="""base filename for reading a Triangle mesh""")
    parser.add_option('--flagsAdd',
                      default='',
                      help="""flags to add when reading/creating Triangle meshes""")

    parser.add_option('-t','--testNum',
                      default=0,
                      help="""which test to use [0]
                      0 --- triIfaceCall3 trimesh fro  .nodes and .ele file,
                            .poly file,  and .node file
                      1 --- triIfaceCall3 with spiral
                      2 --- exAdhPyUtilMesh0 with spiral
                      3 --- exAdhPyUtilMesh0 with la 
                      4 --- exAdhPyUtilMesh0 with user specified file, and flags 
                      5 --- exAdhPyUtilMesh1 convert adhPyMesh for rectangle
                      """)
    parser.add_option('-v','--verbose',
                      default=0,
                      help="""level of verbosity in simulation [0]""")

    #get options
    options, args = parser.parse_args(sys.argv[1:]) #get command line args

    #
    verbose  = int(options.verbose)
    testNum  = int(options.testNum)
    adhPyDir = str(options.adhPyDir)
    flagsAdd = str(options.flagsAdd)
    baseFlags= str(options.baseFlags)
    filebase = str(options.filebase)
    
    sys.path.insert(0,adhPyDir)

    if testNum == 0:
        #look at node and ele files
        filebase="trimesh"
        flagsAdd="r"
        triIfaceCall3(filebase=filebase,flagsAdd=flagsAdd,
                      verbose=verbose)
        #generat from poly file
        filebase="trimesh"
        flagsAdd="p"
        triIfaceCall3(filebase=filebase,flagsAdd=flagsAdd,
                      verbose=verbose)
        #generat from node file
        filebase="trimesh"
        flagsAdd=""
        triIfaceCall3(filebase=filebase,flagsAdd=flagsAdd,
                      verbose=verbose)
    elif testNum == 1:
        filebase="../examples/spiral"
        baseFlags="en" #base 1?
        flagsAdd =""   #just node
        triIfaceCall3(filebase=filebase,flagsAdd=flagsAdd,
                      baseFlags=baseFlags,verbose=verbose)
    elif testNum == 2:
        filebase="../examples/spiral"
        baseFlags="en" #base 1?
        flagsAdd =""   #just node
        exAdhPyUtilMesh0(filebase=filebase,flagsAdd=flagsAdd,
                         baseFlags=baseFlags,verbose=verbose)
    elif testNum == 3:
        filebase="../examples/la"
        baseFlags="en" #base 1?
        flagsAdd ="pa"   #from poly file
        exAdhPyUtilMesh0(filebase=filebase,flagsAdd=flagsAdd,
                         baseFlags=baseFlags,verbose=verbose)
    elif testNum == 4:
        exAdhPyUtilMesh0(filebase=filebase,flagsAdd=flagsAdd,
                         baseFlags=baseFlags,verbose=verbose)
    elif testNum == 5:
        exAdhPyUtilMesh1(baseFlags=baseFlags,verbose=verbose)
    else:
        print 'PROBLEM testNum= ',testNum,' not recognized'
    #end if
