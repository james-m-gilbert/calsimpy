# -*- coding: utf-8 -*-
"""
classes for dealing with the groudnwater model data in Calsim3
based on the 'readCS3GW.v2019xx.py files from summer of 2019

Created on Thu Oct  3 16:37:34 2019

@author: jmgilbert
"""

import cs_util as util
from collections import OrderedDict as odict

import numpy as np
import pandas as pnd
import geopandas as gpd
from fiona.crs import from_epsg
import os, sys
from shapely.geometry import Polygon, MultiPolygon, mapping, Point, LineString
from shapely.ops import transform
import pyproj
try:
    import geojson
except:
    print('geojson library not installed...')
import datetime as dt

M2_TO_FT2 = 100.*100./2.54/2.54/12./12.

# import the DSS read/write libraries
sys.path.append(r'C:\Users\jmgilbert\02_Projects\CalSim_Utilities\Python_Functions\Python_DSS')
import dss3_functions_reference as dss
import AuxFunctions as af


from collections import OrderedDict as Odict



#def progress(count, total, status=''):
#    bar_len = 60
#    filled_len = int(round(bar_len * count / float(total)))
#
#    percents = round(100.0 * count / float(total), 1)
#    bar = '=' * filled_len + '-' * (bar_len - filled_len)
#
#    sys.stdout.write('\r[%s] %s%s ...%s\r' % (bar, percents, '%', status))
#    sys.stdout.flush()  # As suggested by Rom Ruben (see: http://stackoverflow.com/questions/3173320/text-progress-bar-in-the-console/27871113#comment50529068_27871113)



class cs3gw():
    
    def __init__(self,projDir='', inFile = 'CVGroundwater.in', cs3reorg=True):
        self.ProjectDir = projDir
        if cs3reorg:
            self.GWDataDir = os.path.join(projDir,'Run', 'CVGroundwater','Data')
        else:
            self.GWDataDir = os.path.join(projDir,'Run', 'common','CVGroundwater','Data')
        self.MainFile = inFile
        self.GWMain = self.getGWMain()
        self.Nodes = {}
        self.NodesGIS = None
        self.Elements = {}
        self.ElementsGIS = None
        self.Strat = {}
        self.InitHdFile = None
        self.InitHd = {}
        self.TileDrain = {}
        self.Rivers = {}
        self.Results = {}
        self.RiverGIS_Line = None
        self.RiverGIS_Point = None
    
    
    def getGWMain(self):
        ''' read a CVGroundwater.in input file used for CS3
        
            - the approach for scraping the input file is really rough, 
             probably should improve this for better interpretability/readibilyt
        
        '''
        
        CVGWinDict  = Odict()
        inFP = os.path.join(self.GWDataDir, self.MainFile)
        
        if not os.path.exists(inFP):
            print("Groundwater model input file not found. Could not proceed")
            return
        
        with open(inFP, 'r') as of:
            ls = of.readlines()
            
        sectNum=0
        lcntr = 0
        # first section
        l = ls[lcntr]
        while l[0].upper()=='C':
            lcntr+=1
            l = ls[lcntr]
        sectNum+=1
        sectList = []
        while l[0].upper()!='C':
            lsp = l.split('/')
            print(lsp[0])
            sectList.append((lsp[1],lsp[0].strip()))
            lcntr+=1
            l = ls[lcntr]
        
        CVGWinDict[sectNum] = sectList
        
        return(CVGWinDict)
        
    def getNodes(self):
        
        node_dict = {}
        
        node_fp = os.path.join(self.GWDataDir, self.GWMain[1][1][1] )
        
       
        if not os.path.exists(node_fp):
            print("Groundwater node input file not found. Could not proceed")
            return
        
        with open(node_fp, 'r') as of:
            rl = of.readlines()
        rlnc = [i for i in rl if i[0].upper()!='C']
        nd = int(rlnc[0].strip().split()[0]) # number of nodes
        fact = float(rlnc[1].strip().split()[0]) # conversion factor for nodal coordinates
        node_dict['ND'] = nd
        node_dict['FACT'] = fact
        node_coords = {}
        for i in rlnc[2:]: # loop through rest of file - assumes the rest is lines of coordinate with ID   X   Y
            ls = i.split()
            nodeID = int(ls[0].strip())
            xcoord = float(ls[1].strip())
            ycoord = float(ls[2].strip())
            node_coords[nodeID] = [xcoord, ycoord]
        
        node_dict['NODES'] = node_coords
        
        self.Nodes = node_dict
        

    def getElements(self):
        
        elem_dict = {}
        
        elem_fp = os.path.join(self.GWDataDir, self.GWMain[1][0][1] )
        
        with open(elem_fp, 'r') as of:
            rl = of.readlines()
        rlnc = [i for i in rl if i[0].upper()!='C']
        ne = int(rlnc[0].strip().split()[0]) # number of nodes
    
        elem_dict['NE'] = ne
    
        elem_nodes = {}
        for i in rlnc[1:]: # loop through rest of file - assumes the rest is lines of coordinate with ID   X   Y
            ls = i.split()
            elemID = int(ls[0].strip())
            ide1 = int(ls[1].strip())
            ide2 = int(ls[2].strip())
            ide3 = int(ls[3].strip())
            ide4 = int(ls[4].strip())
            elem_nodes[elemID] = [ide1, ide2, ide3, ide4]
        
        elem_dict['ELEMS'] = elem_nodes   
        
        self.Elements = elem_dict
        
        
    def getDrain(self):
    
        drn_dict = {}
        
        tiledrn_fp = os.path.join(self.GWDataDir, self.GWMain[1][13][1])
        with open(tiledrn_fp, 'r') as of:
            rl = of.readlines()
        rlnc = [i for i in rl if i[0].upper()!='C']
        ntd = int(rlnc[0].strip().split()[0]) # number of tile drains
        fact = float(rlnc[1].strip().split()[0]) # conversion factor for tile drain elevations
        factcdc = float(rlnc[2].strip().split()[0]) #conversion factor for tile rain conductances - for spatial component only?
        tunit = str(rlnc[3].strip().split()[0])
        drn_dict['ND'] = ntd
        drn_dict['FACT'] = fact
        drn_dict['FACTCDC'] = factcdc
        drn_dict['TUNIT'] = tunit
        drn_bc = {}
        for i in rlnc[4:]: # loop through rest of file - assumes the rest is lines of coordinate with ID   X   Y
            if i=='\n':
                continue
            ls = i.split()
            try:
                notes = i.split('/')[1].strip()
            except:
                notes = ''
            nodeID = int(ls[0].strip())
            elevdrn = float(ls[1].strip())
            cdcdrn = float(ls[2].strip())
            istrmdrn = int(ls[3].strip())
            drn_bc[nodeID] = [elevdrn, cdcdrn, istrmdrn, notes]
        
        tdrnDF = pnd.DataFrame.from_dict(drn_bc, orient='index', columns=['DrnElev_ft','Cond','StrNode','Notes'])
        tdrnDF['NodeID'] = tdrnDF.index *-1
        drn_dict['TILEDRN'] = tdrnDF
        
        
        self.TileDrain = drn_dict
        
        
    def getStrat(self, calcBotElev=True):
        
        strat_fp = os.path.join(self.GWDataDir, self.GWMain[1][2][1])
        
        # reads in IWFM stratigraphy type file
        with open(strat_fp) as of:
            cvstrat_lines = of.readlines()
        cvstrat_lines_nocomment = [l for l in cvstrat_lines if l[0]!='C']
        
        # file starts with 2 lines that indicate 1) number of layers and 2) conversion factor for length
        # file doesn't contain number of nodes(?) - assuming this has been read in from elsewhere? - hard coding it at 1393
    
        nd = self.Nodes['ND']
        NL = int(cvstrat_lines_nocomment[0].split("/")[0].strip())
        FACT = int(cvstrat_lines_nocomment[1].split("/")[0].strip())
        
        
        # the following code assumes 3 layers and is hard-wired for this; would 
        # need to be generalized to accommodate n-number of layers
        strat_list = []
        for n in range(nd):
            ls = cvstrat_lines_nocomment[n+2].split()
            nodeID = int(ls[0])
            gw_hyd_id = "L1:GW%s" %nodeID
            strat_list.append([int(ls[0]),float(ls[1]),float(ls[2]), float(ls[3]),float(ls[4]),float(ls[5]),float(ls[6]), float(ls[7]),gw_hyd_id])
            
        strat_df = pnd.DataFrame(strat_list, columns=("NodeID","Elev","L1_bc1","L1_b1","L2_bc1","L2_b1","L3_bc1","L3_b1", "GW_Hydrograph_ID"))
        
        if calcBotElev:
            strat_df['L1bot'] = strat_df.Elev - (strat_df.L1_bc1 + strat_df.L1_b1)
            strat_df['L2bot'] = strat_df.L1bot - (strat_df.L2_bc1 + strat_df.L2_b1)
            strat_df['L3bot'] = strat_df.L2bot - (strat_df.L3_bc1 + strat_df.L3_b1)
        
        self.Strat = {'NL': NL, 'data':strat_df}
        

    def getRivers(self):
        riv_dict = odict()
        
        riv_fp = os.path.join(self.GWDataDir, self.GWMain[1][3][1])
        with open(riv_fp, 'r') as of:
            rl = of.readlines()
        rlnc = [i for i in rl if i[0].upper()!='C']
        nrh = int(rlnc[0].strip().split()[0]) # number of stream reaches
        nr = int(rlnc[1].strip().split()[0])  # number of stream nodes modeled
        nrtb = int(rlnc[2].strip().split()[0]) # number of data points in stream rating tables

        riv_dict['NRH'] = nrh
        riv_dict['NR'] = nr
        riv_dict['NRTB'] = nrtb
        riv_dict['Reach_Data'] = {}

        reach_start_row = 3
        #cntr=0
        for r in range(nrh):  # loop through all the reaches
            reach_dict = odict()
            cntr = 0
            reach_row = rlnc[reach_start_row].strip().split()
            print("Reach row: %s" %rlnc[reach_start_row])
            reach_id = int(reach_row[0])
            upstrm_node = int(reach_row[1])
            dnstrm_node = int(reach_row[2])
            outflw_node = int(reach_row[3])
            
            cntr+=1
            
            node_dict = odict()
            for rd in range(upstrm_node,dnstrm_node+1): # loop through reach nodes
                
                node_row = rlnc[reach_start_row+cntr].strip().split()
                print("Node row: %s  %s   %s" %(node_row[0], node_row[1], node_row[2]))
                strm_node = int(node_row[0])
                igw = int(node_row[1])
                irgst = int(node_row[2])
                node_dict[strm_node]= {'GW_Node':igw,'Subregion': irgst}
                cntr+=1
            
            riv_dict['Reach_Data'][reach_id] = {'UpstreamNode': upstrm_node,
                                                'DownstreamNode': dnstrm_node,
                                                'OutflowNode': outflw_node,
                                                'ReachNodes': node_dict}    
            reach_start_row += cntr   
        
        self.Rivers = riv_dict
        
    def getInitialHead(self, valPerRow=10, startRow=16, **kwargs):
      
        inithd_fp = ''
        for f in self.GWMain[1]:
            if 'INITIAL CONDITIONS DATA FILE' in f[0]:
                inithd_fp = os.path.join(self.GWDataDir, f[1])
        if inithd_fp=='':
            raise FileNotFoundError(f'Initial head file not found in CVGroundwater input - check file')
        if not os.path.exists(inithd_fp):
            raise FileNotFoundError(f'Specified initial head file {inithd_fp} not found')
        self.InitHdFile = inithd_fp
        
        if len(self.Nodes) == 0:
            # haven't called the nodes funciton yet - do that now
            self.getNodes()
        
        nnodes = self.Nodes['ND'] # set number of nodes
                      
            
    #def read_initHd(inithd_fp,nnodes, valPerRow=10, startRow=15):
        with open(inithd_fp, 'r') as of:
            all_lines = of.readlines()
        
        # count the number of comment lines at the top of the file
        initDateTime = all_lines[0].strip() #dt.datetime.strptime(all_lines[0].strip(), '%m/%d/%Y_%H:%M')
        l = all_lines[1]
        cmtLen = 0
        while l[0].upper() == 'C':
            cmtLen+=1
            l = all_lines[1+cmtLen]
        cmtLen+=1 #<-- have to adjust for the first line being the date
        layerIndicator = float(all_lines[cmtLen].strip())
        cmtLen+=1
        #nnodes = 1393
        if nnodes % valPerRow >0:
            numRows = int(nnodes//valPerRow + 1)
        else:
            numRows = int(nnodes//valPerRow)
        # this is the super naive/quick&dirty way to do this - not robust or general
        # there are 10 values per row, 1393 nodes total, so 140 rows hold the data for
        # one layer
        init_hd_dict ={}
        ncntr = 1
        row_cntr = 1
        for l in all_lines[cmtLen:cmtLen+numRows]: #[startRow:startRow+numRows]:
            ls = l.split()
            print(l)
            print(len(ls))
            if row_cntr ==numRows:
                numitems = nnodes%valPerRow
            else:
                numitems = valPerRow
            print("numitems = %s" %numitems)
            for i in range(numitems):
                init_hd_dict[ncntr] = float(ls[i])
                ncntr+=1
            row_cntr+=1
            
        #return(init_hd_dict)
        self.InitHd = init_hd_dict


        
    def node2GIS(self, **kwargs):
        '''
            convert the node information read from the input file into a 
            projected geospatial format
            
            **kwargs:   specify which other proprties to include as attributes
                        options include:  
                                        - props = ['drn'  -> include tile drains
                                                  'initHd' -> initial head
                                                  'strat' -> stratigraphy, including ground surface elvation]
            
        '''
        all_nodes_pt = []
        pts = []  # list of Shapely point objects
        
#        # this method for projecting adapted from https://gis.stackexchange.com/questions/127427/transforming-shapely-polygon-and-multipolygon-objects
#        project = lambda x, y: pyproj.transform(
#                            pyproj.Proj(init='epsg:26910'), #, preserve_flags=True),  #source
#                            pyproj.Proj(init='epsg:4326'), # preserve_flags=True),
#                            x, y)
        
        for nd in self.Nodes['NODES']:
            #print(nd)
            node_row = []

            coords = self.Nodes['NODES'][nd]
            plist=[tuple(coords)]
            node_row.append(nd)
            node_row.append(plist[0][0])
            node_row.append(plist[0][1])
            node_lab = ['NodeID','CoordY','CoordX']
            
            if 'props' in kwargs:
                if 'drn' in kwargs['props']:
                    # lookup tile drain for this node
                    if -1*nd in self.TileDrain['TILEDRN'].keys():
                        drnelev = self.TileDrain['TILEDRN'][-1*nd][0]
                        drncond = self.TileDrain['TILEDRN'][-1*nd][1]
                        drnStrNod = self.TileDrain['TILEDRN'][-1*nd][2]
                    else:
                        drnelev = np.NaN
                        drncond = np.NaN
                        drnStrNod = np.NaN
                        
                    node_row.append(drnelev)
                    node_row.append(drncond)
                    node_row.append(drnStrNod)
                    node_lab.append('DrnElv_ft')
                    node_lab.append('DrnCond')
                    node_lab.append('DrnDestSNode')
                    
                
                
                #look up assigned initial condition for node
                if 'initHd' in kwargs['props']:
                    thisInitHd = self.InitHd[nd]
                    node_row.append(thisInitHd)
                    node_lab.append('InitlHd_ft')
                    
            
                if 'strat' in kwargs['props']:
                    # look up ground surface elevation
                    thisGSElev = self.Strat['data'][self.Strat['data']['NodeID']==nd].Elev.iloc[0]
                    
                    # look up layer bottom elevations
                    l1botElev = self.Strat['data'][self.Strat['data']['NodeID']==nd].L1bot.iloc[0]
                    l2botElev = self.Strat['data'][self.Strat['data']['NodeID']==nd].L2bot.iloc[0]
                    l3botElev = self.Strat['data'][self.Strat['data']['NodeID']==nd].L3bot.iloc[0]
                    
                    node_row.append(thisGSElev)
                    node_row.append(l1botElev)
                    node_row.append(l2botElev)
                    node_row.append(l3botElev)
                    node_lab.append('SurfElev_ft')
                    node_lab.append('L1botElev_ft')
                    node_lab.append('L2botElev_ft')
                    node_lab.append('L3botElev_ft')
                    
                    
#            # calc difference between ground surface and init hd elev
#            initHdGSdiff = thisGSElev - thisInitHd
            
            #all_nodes_pt.append([nd, plist[0][0], plist[0][1], drnelev, drncond, drnStrNod, thisInitHd,thisGSElev,initHdGSdiff ])
            all_nodes_pt.append(node_row)
            
            #pts.append(transform(project,Point(plist[0])))
            pts.append(Point(plist[0]))
        
        nodeDF = gpd.GeoDataFrame(data=all_nodes_pt,columns=node_lab,
                                    geometry=pts)
        
        
        # set the projection - assuming UTM 10N (i think)
        prjStr = pyproj.Proj('+init=epsg:26910').definition_string()
        crsDict = {}
        for ps in prjStr.split():
            if ps=='+no_defs' or ps=='no_defs':
                crsDict[ps] = True
            else:
                psi = ps.split('=')
                crsDict[psi[0]] = psi[1]
                
        nodeDF.crs = crsDict
        
        self.NodesGIS = nodeDF
    
    def elem2GIS(self):
        
        nc = self.Nodes['NODES']
        all_elem_poly = []
        polys = []  # list of Shapely polygon objects
        elem_areas = []
        elem_ids = []
        
        i=0
        for e in self.Elements['ELEMS']: #loop through elem id's

            i +=1
            # provide a progress bar in the console
            util.progress(i, len(self.Elements['ELEMS']), status='Building element GIS dataset')

            elem_nodes = self.Elements['ELEMS'][e]
            plist = []
            avgTopElev=avgL1bot=avgL2bot=avgL3bot = 0.
            if elem_nodes[-1]!=0:
                num_nodes = 4
            else:
                num_nodes = 3
            for n in range(num_nodes):
                thisnode=elem_nodes[n]
                
                stratdata = self.Strat['data'].copy()
                thisstrat = stratdata[stratdata['NodeID']==thisnode]
                topelev = thisstrat['Elev'].iloc[0]
                l1b = thisstrat['L1_b1'].iloc[0] + thisstrat['L1_bc1'].iloc[0]
                l2b =thisstrat['L2_b1'].iloc[0]+ thisstrat['L2_bc1'].iloc[0]
                l3b = thisstrat['L3_b1'].iloc[0]+ thisstrat['L3_bc1'].iloc[0]
                l1elevbot = topelev-l1b
                l2elevbot = l1elevbot-l2b
                l3elevbot = l2elevbot - l3b
                avgTopElev = avgTopElev + topelev/float(num_nodes)
                avgL1bot = avgL1bot + l1elevbot/float(num_nodes)
                avgL2bot = avgL2bot + l2elevbot/float(num_nodes)
                avgL3bot = avgL3bot + l3elevbot/float(num_nodes)
                
                plist.append(tuple(nc[thisnode]))
             
            all_elem_poly.append([e, plist, avgTopElev, avgL1bot, avgL2bot, avgL3bot])
            polys.append(Polygon(plist)) #transform(project,Polygon(plist)))
            #all_elem_poly.append(plist)
            thisPoly = Polygon(plist)
            elem_areas.append(thisPoly.area)
            elem_ids.append(e)
            #polys.append(thisPoly)
            
        
        elem_polygons = MultiPolygon(polys)
        elemTopElev = [v[2] for v in all_elem_poly]
        elemL1botElev = [v[3] for v in all_elem_poly]
        elemL2botElev = [v[4] for v in all_elem_poly]
        elemL3botElev = [v[5] for v in all_elem_poly]
       
        elemGDF = gpd.GeoDataFrame(data=zip(elem_ids,elem_areas, elemTopElev, 
                                            elemL1botElev, elemL2botElev,
                                            elemL3botElev),
                                columns=['elem_id','elem_area_m2','ElevTop','ElevL1Bot','ElevL2Bot','ElevL3Bot'],geometry=polys)
        
        elemGDF['elem_area_ft'] = elemGDF['elem_area_m2']*M2_TO_FT2
        elemGDF['elem_area_ac'] = elemGDF['elem_area_ft']/42560
        
        #elemGDF.geometry = elemGDF.geometry.to_crs(epsg=26910)
        
        prjStr = pyproj.Proj('+init=epsg:26910').definition_string()
        crsDict = {}
        for ps in prjStr.split():
            if ps=='+no_defs' or ps=='no_defs':
                crsDict[ps] = True
            else:
                psi = ps.split('=')
                crsDict[psi[0]] = psi[1]
        elemGDF.crs = crsDict # {'init': 'epsg:26910'}
        
        self.ElementsGIS = elemGDF
        
    def riv2GIS(self):
        
        nc = self.Nodes['NODES']
        nrh = self.Rivers['NRH']
        nr = self.Rivers['NR']
        
        all_riv_pts = []
        lines = []  # list of Shapely line objects
        strnodeID = []
        gwnodeID = []
        reachID = []
        reachID_lines=[]
        subrgn = []
        
        reachUp = []
        reachDwn = []
        reachOut = []
        
        i=0
        for r in self.Rivers['Reach_Data']: #loop through river reaches

            i +=1
            # provide a progress bar in the console
            util.progress(i, nrh, status='Building rivers GIS dataset')

            reach_nodes = self.Rivers['Reach_Data'][r]['ReachNodes']
            
            llist = []
            plist=[]

            for n in reach_nodes:  # loop through gw nodes in reach
                
                gw_node = reach_nodes[n]['GW_Node']
                
                all_riv_pts.append(Point(tuple(nc[gw_node])))
                
                llist.append(tuple(nc[gw_node]))
                strnodeID.append(n)
                gwnodeID.append(gw_node)
                subrgn.append(reach_nodes[n]['Subregion'])
                reachID.append(r)

             
            #all_riv_line.append([r, llist, avgTopElev, avgL1bot, avgL2bot, avgL3bot])
            lines.append(LineString(llist)) #transform(project,Polygon(plist)))
            reachID_lines.append(r)
            reachUp.append(self.Rivers['Reach_Data'][r]['UpstreamNode'])
            reachDwn.append(self.Rivers['Reach_Data'][r]['DownstreamNode'])
            reachOut.append(self.Rivers['Reach_Data'][r]['OutflowNode'])
       
        riv_Pt_GDF = gpd.GeoDataFrame(data=zip(strnodeID,gwnodeID, subrgn, reachID),
                                columns=['str_id','gw_id','Subregion','ReachID'],geometry=all_riv_pts)
        
        riv_Line_GDF = gpd.GeoDataFrame(data=zip(reachID_lines, reachUp, reachDwn,
                                                 reachOut),
                                        columns=['reach_id','UpNode','DwnNode','OutNode'],
                                        geometry=lines)
        
        prjStr = pyproj.Proj('+init=epsg:26910').definition_string()
        crsDict = {}
        for ps in prjStr.split():
            if ps=='+no_defs':
                crsDict[ps] = True
            else:
                psi = ps.split('=')
                crsDict[psi[0]] = psi[1]
        riv_Pt_GDF.crs = crsDict # {'init': 'epsg:26910'}
        riv_Line_GDF.crs = crsDict
        
        self.RiverGIS_Point = riv_Pt_GDF
        self.RiverGIS_Line = riv_Line_GDF

        
class gwResults(cs3gw):

    def __init__(self, cs3gw, ModelRunID=''):
        self.RunID = ModelRunID
        self.Model = cs3gw
        self.Pumping = {}
        self.DeepPerc = {}
        self.Head = {}
        self.StorageChange = None
        self.StreamGW = None
        self.LateralFlow = None
        
        
    def getGWPumping(self, gw_startDate = None, gw_endDate = None, ctime24=True):
        
        # set up path to results file that has element-by-element results
        resultsFP = os.path.join(self.Model.ProjectDir, 'CONV','DSS',self.Model.GWMain[1][17][1])
        
        if not os.path.exists(resultsFP):
            print("Specified results file %s does not exist...nothing can be done now..." %resultsFP)
            return
        
        # create the DSS file object
        gwresDSS = util.dssFile(fp=resultsFP)

        gwresDSS.openDSS()
            
        
        numElems = self.Model.Elements['NE']
        NL = self.Model.Strat['NL']
        lyrs = ['L%d' %(x+1) for x in range(NL)] 
        lyrs.append('LT')
        
        Apt = 'IWFM'
        Cpt = 'PUMPING'
        Dpt = ''
        Ept = '1MON'
        Fpt = 'PUMPING_AT_ALL_ELEMENTS'
        

#        df = pnd.DataFrame()
#        df['Datetime'] = pnd.to_datetime(dateList)
#        df.set_index('Datetime', inplace=True)
        i=0
        # get gw pumping first
        for ne in range(1,numElems+1):

            for ly in lyrs:
                # provide a progress bar in the console
                util.progress(i, (NL+1)*numElems, status='Retrieving %s GW Pumping Time Series' %((NL+1)*numElems))                
                
                Bpt = ly+':E'+str(ne)
                #print('retrieving time series data for: %s' %Bpt)
#                var = Bpt
                cpath = '/'+'/'.join([Apt,Bpt, Cpt,Dpt, Ept, Fpt]) + '/'
                
                
                i +=1

            
                thisvar = util.dssVar(gwresDSS, cpath=cpath)
            
                thisvar.getRTS(gw_startDate,ctime24=True, endDateTime=gw_endDate)

                thispump = thisvar.RTS
                
                self.Pumping[Bpt] = [thispump, thisvar.Units, thisvar.RecordType]

            
        gwresDSS.closeDSS()
        
    def getGWDeepPerc(self, gw_startDate = None, gw_endDate = None, ctime24=True):
        
        # set up path to results file that has element-by-element results
        resultsFP = os.path.join(self.Model.ProjectDir, 'CONV','DSS',self.Model.GWMain[1][18][1])
        
        if not os.path.exists(resultsFP):
            print("Specified results file %s does not exist...nothing can be done now..." %resultsFP)
            return
        
        # create the DSS file object
        gwresDSS = util.dssFile(fp=resultsFP)

        gwresDSS.openDSS()
            
        
        numElems = self.Model.Elements['NE']
        NL = self.Model.Strat['NL']
        lyrs = ['L%d' %(x+1) for x in range(NL)] 
        lyrs.append('LT')
        
        Apt = 'IWFM'
        Cpt = 'DEEP_PERC'
        Dpt = ''
        Ept = '1MON'
        Fpt = 'DEEP_PERC_AT_ALL_ELEMENTS'
        

#        df = pnd.DataFrame()
#        df['Datetime'] = pnd.to_datetime(dateList)
#        df.set_index('Datetime', inplace=True)
        i=0
        # get gw pumping first
        for ne in range(1,numElems+1):
            # provide a progress bar in the console
            util.progress(i,numElems, status='Retrieving %s GW Deep Perc Time Series' %(numElems))

            Bpt = 'E'+str(ne)
            #print('retrieving time series data for: %s' %Bpt)
#                var = Bpt
            cpath = '/'+'/'.join([Apt,Bpt, Cpt,Dpt, Ept, Fpt]) + '/'
            
            
            i +=1

        
            thisvar = util.dssVar(gwresDSS, cpath=cpath)
        
            thisvar.getRTS(gw_startDate,ctime24=True, endDateTime=gw_endDate)

            thisDP = thisvar.RTS
            
            self.DeepPerc[Bpt] = [thisDP, thisvar.Units, thisvar.RecordType]

            
        gwresDSS.closeDSS()
        
    def getGWHead(self, gw_startDate = None, gw_endDate = None, ctime24=True):
        
        # set up path to results file that has element-by-element results
        resultsFP = os.path.join(self.Model.ProjectDir, 'CONV','DSS',self.Model.GWMain[1][16][1])
        
        if not os.path.exists(resultsFP):
            print("Specified results file %s does not exist...nothing can be done now..." %resultsFP)
            return
        
        # create the DSS file object
        gwresDSS = util.dssFile(fp=resultsFP)

        gwresDSS.openDSS()
            
        
        numNodes = self.Model.Nodes['ND']
        NL = self.Model.Strat['NL']
        lyrs = ['L%d' %(x+1) for x in range(NL)] 

        
        Apt = 'IWFM'
        Cpt = 'HEAD'
        Dpt = ''
        Ept = '1MON'
        Fpt = 'GW_HEAD_AT_ALL_NODES'
        

#        df = pnd.DataFrame()
#        df['Datetime'] = pnd.to_datetime(dateList)
#        df.set_index('Datetime', inplace=True)
        i=0
        # get gw pumping first
        for nd in range(1,numNodes+1):
            # provide a progress bar in the console
            util.progress(i, (NL+1)*numNodes, status='Retrieving %s GW Head Time Series' %((NL+1)*numNodes))
            for ly in lyrs:
                Bpt = ly+':GW'+str(nd)
                #print('retrieving time series data for: %s' %Bpt)
#                var = Bpt
                cpath = '/'+'/'.join([Apt,Bpt, Cpt,Dpt, Ept, Fpt]) + '/'
                
                
                i +=1

            
                thisvar = util.dssVar(gwresDSS, cpath=cpath)
            
                thisvar.getRTS(gw_startDate,ctime24=True, endDateTime=gw_endDate)

                thishead = thisvar.RTS
                
                self.Head[Bpt] = [thishead, thisvar.Units, thisvar.RecordType]

            
        gwresDSS.closeDSS()