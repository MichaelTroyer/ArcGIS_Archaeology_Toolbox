# -*- coding: utf-8 -*-

"""
Arch general Functions

Michael D. Troyer

mtroyer@blm.gov
719-269-8587
"""

import arcpy
import copy
import datetime
import getpass
import re
import os
import sys
import traceback

import pandas as pd

from arcpy import mapping
from arcpy import env

from collections import defaultdict

env.addOutputsToMap = False
arcpy.env.overwriteOutput = True


##---Classes---------------------------------------------------------------------------------------

class pyt_log(object):
    """A custom logging class that can simultaneously write to the console - AddMessage,
       write to an optional logfile, and/or a production report..."""
       
    def __init__(self, report_path, log_path, log_active=True):
        self.report_path = report_path
        self.log_path = log_path
        self.log_active = log_active
    
    def _write_arg(self, arg, path, starting_level=0):
        """Accepts a [path] txt from open(path) and unpacks the data [arg]"""
        level = starting_level
        txtfile = open(path, 'a')
        if level == 0:
            txtfile.write("_"*80)
        if type(arg) == dict:
            txtfile.write("\n"+(level*"\t")+(str(arg))+"\n")
            txtfile.write((level*"\t")+str(type(arg))+"\n")
            for key, value in arg.items():
                txtfile.write((level*"\t\t")+(str(key))+": "+(str(value))+"\n")
                if hasattr(value, '__iter__'):
                    txtfile.write((level*"\t")+"Values:"+"\n")
                    txtfile.close()
                    for val in value:
                        self._write_arg(val, path, starting_level=level+1)
                txtfile.close()
        else:
            txtfile.write("\n"+(level*"\t")+(str(arg))+"\n")
            txtfile.write((level*"\t")+str(type(arg))+"\n")
            if hasattr(arg, '__iter__'): #does not include strings
                txtfile.write((level*"\t")+"Iterables:"+"\n")
                txtfile.close()
                for a in arg:
                    self._write_arg(a, path, starting_level=level+1)

    def _writer(self, msg, path, *args):
        """A writer to write the msg, and unpacked variable"""
        if os.path.exists(path):
            write_type = 'a'
        else:
            write_type = 'w'
        with open(path, write_type) as txtfile:
            txtfile.write("\n"+msg+"\n")
            txtfile.close()
            if args:
                for arg in args:
                    self._write_arg(arg, path)

    def console(self, msg):
        """Print to console only"""
        arcpy.AddMessage(msg)

    def report(self, msg):
        """Write to report only"""
        self._writer(msg, path=self.report_path)

    def logfile(self, msg, *args):
        """Write to logfile only"""
        if self.log_active:
            path = self.log_path
            self._writer(msg, path, *args)
            
    def log_report(self, msg, *args):
        """Write to logfile and report only"""
        self.report(msg)
        self.logfile(msg, *args)
        
    def log_all(self, msg, *args):
        """Write to all"""
        self.console(msg)
        self.report(msg)
        self.logfile(msg, *args)


##---Functions-------------------------------------------------------------------------------------

def blast_my_cache():
    """Delete in memory tables and feature classes
       reset to original worksapce when done"""

    # get the original workspace location
    orig_workspace = arcpy.env.workspace
    
    # Set the workspace to in_memory
    arcpy.env.workspace = "in_memory"
    # Delete all in memory feature classes
    fcs = arcpy.ListFeatureClasses()
    if len(fcs) > 0:
        for fc in fcs:
            arcpy.Delete_management(fc)
    # Delete all in memory tables
    tbls = arcpy.ListTables()
    if len(tbls) > 0:
        for tbl in tbls:
            arcpy.Delete_management(tbl)

    # Reset the workspace
    arcpy.env.workspace = orig_workspace


def buildWhereClauseFromList(table, field, valueList):
    """Takes a list of values and constructs a SQL WHERE
    clause to select those values within a given field and table."""

    # Add DBMS-specific field delimiters
    fieldDelimited = arcpy.AddFieldDelimiters(arcpy.Describe(table).path, field)

    # Determine field type
    fieldType = arcpy.ListFields(table, field)[0].type

    # Add single-quotes for string field values
    if str(fieldType) == 'String':
        valueList = ["'%s'" % value for value in valueList]

    # Format WHERE clause in the form of an IN statement
    whereClause = "%s IN(%s)" % (fieldDelimited, ', '.join(map(str, valueList)))
    return whereClause


# Constants
degToRad = math.pi/180.0
radToDeg = 180.0/math.pi


# Functions
def p1p2Dist(pnt1, pnt2):
    """
    Returns the distance between two different points
    """
    dY = pnt2[1] - pnt1[1]
    dX = pnt2[0] - pnt1[0]
    dist = math.hypot(dX,dY)
    return dist
  
  
def p1p2Angle(pnt1, pnt2):
    """
    Returns the angle between two different points relative to the x axis
    """
    dY = pnt2[1] - pnt1[1]
    dX = pnt2[0] - pnt1[0]
    theAngle = math.atan2(dY,dX)*radToDeg
    return theAngle


def extentPnts(pnts):
    """Returns the min, max X, Y of input points"""
    xList = []; yList = []
    for pnt in pnts:
        xList.append(pnt[0]); yList.append(pnt[1])
    L = min(xList); R = max(xList); B = min(yList); T = max(yList)
    LL = [L,B]; UL = [L,T]; UR = [R,T]; LR = [R,B]
    return [LL, UL, UR, LR]


def extentCenter(pnts):
    """Returns the average and median X, Y coordinates of a series
       of input points"""
    xList = []; yList = []
    for pnt in pnts:
        xList.append(pnt[0]); yList.append(pnt[1])
    L = min(xList); R = max(xList); B = min(yList); T = max(yList)
    Xcent = (R - L)/2.0 + L; Ycent = (T - B)/2.0 + B
    return [Xcent, Ycent]


def polyAngles(pnts):
    """
    Requires:  a list of points forming a polygon
    Returns:   a list containing from pnt, to pnt, anAngle, aDistance
    """
    pnts2 = pnts[:]
    if pnts2[0] != pnts2[-1]:
        N = len(pnts2)
        # Add the first point to the list
        pnts2.append(pnts2[0])
    else:
        N = len(pnts2) - 1
    angleList = []  
    for i in range(1, len(pnts2)):
        pnt1 = pnts2[i-1];  pnt2 = pnts2[i]
        if pnt1 != pnt2:
            theAngle = p1p2Angle(pnt1, pnt2)
            theDist = p1p2Dist(pnt1, pnt2)
            if i < N:
                angleList.append([i-1, i, theAngle, theDist])
            else:
                angleList.append([i-1, 0, theAngle, theDist])
    return angleList

  
def transRotatePnts(pnts, pntCent, angle):
    """
    Requires a list of points, an extent point and an angle in degrees.
    Translates and rotates points about the origin using the negative angle
    formed between the points.
    """
   
    X0 = pntCent[0]; Y0 = pntCent[1]
    # Check for a duplicate closure point
    if pnts[0] != pnts[-1]:
        N = len(pnts)
    else:
        N = len(pnts)-1
    #translate and rotate shapes and determine area
    angle = angle*(-1.0) #reverse the rotation
    cosXY = math.cos(degToRad * angle) 
    sinXY = math.sin(degToRad * angle)
    rotPnts = []
    for j in range(0, N):
        X1 = pnts[j][0] - X0; Y1 = pnts[j][1] - Y0
        X = (X1 * cosXY) - (Y1 *sinXY)
        Y = (X1 * sinXY) + (Y1 *cosXY)
        pnt = [X, Y]
        rotPnts.append(pnt)
    #Return the rotated points and the centre
    return rotPnts
  
  
def minRect(pnts):
    """
    Determines the minimum area rectangle for a shape represented
    by a list of points
    """
    areaDict = {}
    angleList = polyAngles(pnts)  #determine the angles
    pntCent = extentCenter(pnts)  #determine center of the extent
    xCent = pntCent[0]; yCent = pntCent[1]
    for angle in angleList:
        rotPnts = transRotatePnts(pnts, pntCent, angle[2])  #slice the angle
        Xs = []; Ys = []
        for pnt in rotPnts:
            Xs.append(pnt[0]); Ys.append(pnt[1])
        #Determine the area of the rotated hull
        Xmin = min(Xs); Xmax = max(Xs); Ymin = min(Ys); Ymax = max(Ys)
        area = (max(Xs) - min(Xs))*(max(Ys) - min(Ys))
        areaDict.update({area:[pntCent, Xmin, Xmax, Ymin, Ymax,angle[2]]})
    #Get the minimum rectangle centred about the origin
    #Rotate the rectangle back
    minArea = min(areaDict.keys())
    a = areaDict.get(minArea)
    Xmin = a[1];  Xmax = a[2]
    Ymin = a[3];  Ymax = a[4]
    angle = a[5] * (-1.0)
    rectPnts = [[Xmin,Ymin], [Xmin,Ymax], [Xmax,Ymax], [Xmax,Ymin]]
    originPnt = [0.0,0.0]
    rotPnts = transRotatePnts(rectPnts, originPnt, angle)
    outList = []
    xyPnts = []
    for pnt in rotPnts:
        XY = [pnt[0] + xCent, pnt[1] + yCent]
        xyPnts.append(XY)
    dx = Xmax - Xmin       
    dy = Ymax - Ymin
    outList = [xyPnts, angle, dx, dy]
    #return the points
    return outList


##---Variables-------------------------------------------------------------------------------------

start_time = datetime.datetime.now()

user = getpass.getuser()

##---Settings--------------------------------------------------------------------------------------

arcpy.env.addOutputsToMap = False

arcpy.env.overwriteOutput = True


class Toolbox(object):   
    def __init__(self):
        self.label = "Arch_General_Functions_v2_2"
        self.alias = "Arch_General_Functions"
        
        # List of tool classes associated with this toolbox
        self.tools = [Arch_General_Functions_v2_2]


class Arch_General_Functions_v2_2(object):
    def __init__(self):
        self.label = "Arch_General_Functions_v2_2"
        self.description = ""
        self.canRunInBackground = True
        
    def getParameterInfo(self):
        
        #Input Target Shapefile
        param0=arcpy.Parameter(
            displayName="Input Feature Class",
            name="Input_Shape",
            datatype="Feature Class",
            parameterType="Required",
            direction="Input")
        
        param1=arcpy.Parameter(
            displayName="Selection based on case value",
            name="Select_Boolean",
            datatype="Boolean",
            parameterType="Optional",
            direction="Input",
            enabled = "False")
        
        param2=arcpy.Parameter(
            displayName="Select Feature Case Field",
            name="Select_Field",
            datatype="String",
            parameterType="Optional",
            direction="Input",
            enabled = "False")
        
        param3=arcpy.Parameter(
            displayName="Select Feature Case Value",
            name="Select_Value",
            datatype="String",
            parameterType="Optional",
            direction="Input",
            enabled = "False")
        
        #Output Location and Name
        param4=arcpy.Parameter(
            displayName="Output Workspace and File Naming Convention",
            name="Out_Name",
            datatype="File",
            parameterType="Required",
            direction="Output")
        
###
# Map Options
###

        #Create Map Document
        param5 = arcpy.Parameter(
            displayName="Create a Map Document",
            name="Create_Map",
            datatype="Boolean",
            parameterType="Optional",
            category = "Map Options",
            direction="Input")
        
        #Select Template
        param6 = arcpy.Parameter(
            displayName="Select Map Template",
            name="Map_Template",
            datatype="string",
            parameterType="Optional",
            category = "Map Options",            
            direction="Input")
        
        #Cultural Resources Report Number
        param7=arcpy.Parameter(
            displayName="Culural Resource Report Number",
            name="Input_CRNum",
            datatype="String",
            parameterType="required", 
            direction="Input")
        
        #Cultural Resources Report Name
        param8=arcpy.Parameter(
            displayName="Cultural Resources Report Name",
            name="Input_CRName",
            datatype="String",
            parameterType="Optional",
            category = "Map Options", 
            direction="Input")
        
        #Author
        param9=arcpy.Parameter(
            displayName="Author",
            name="Input_Author",
            datatype="String",
            parameterType="Optional",
            category = "Map Options",
            direction="Input")
        
###
# SHPO data prep
###

        #Populate SHPO Layers
        param10 = arcpy.Parameter(
            displayName="Populate SHPO Layers",
            name="SHPO_Option",
            datatype="boolean",
            parameterType="Optional",
            category = "SHPO Options",
            direction="Input")
        
        #Select SHPO Layers
        param11 = arcpy.Parameter(
            displayName="Select SHPO Layer to Populate",
            name="SHPO_Layers",
            datatype="string",
            parameterType="Optional",
            enabled = "False",
            category = "SHPO Options",
            direction="Input")
        
        #Input Site/Survey ID
        param12 = arcpy.Parameter(
            displayName="Input SHPO Site or Survey ID",
            name="Input_ID",
            datatype="string",
            parameterType="Optional",
            enabled = "False",
            category = "SHPO Options",
            direction="Input")
        
###
# Data options
###

        #Output table of spatial information in PLSS
        param13 = arcpy.Parameter(
            displayName="Write PLSS Spatial Location Data to Table",
            name="Output_PLSS",
            datatype="Boolean",
            parameterType="Optional",
            category = "Data Options",
            direction="Input")
        
        #Clip Elevation
        param14 = arcpy.Parameter(
            displayName="Clip Degital Elevation Model (DEM)",
            name="Output_DEM",
            datatype="Boolean",
            parameterType="Optional",
            category = "Data Options",
            direction="Input")
        
        #Clip and Dissolve Sites
        param15 = arcpy.Parameter(
            displayName="Clip Sites",
            name="Output_Sites",
            datatype="Boolean",
            parameterType="Optional",
            category = "Data Options",
            direction="Input")
        
        #Clip and Dissolve Surveys
        param16 = arcpy.Parameter(
            displayName="Clip Surveys",
            name="Output_Surveys",
            datatype="Boolean",
            parameterType="Optional",
            category = "Data Options",
            direction="Input")

        #Create Lit Search Tables
        param17 = arcpy.Parameter(
            displayName="Copy Lit Search Tables",
            name="Lit_Tables",
            datatype="Boolean",
            parameterType="Optional",
            category = "Data Options",
            direction="Input")
        
        params = [param0, param1, param2, param3, param4, param5, param6,
                  param7, param8, param9, param10, param11, param12, param13,
                  param14, param15, param16, param17]
        
        return params


    def isLicensed(self):
        return True


    def updateParameters(self, params):
        #Param0 - Shape In and handle subselections
        params[0].filter.list = ["Polygon"]
        
        if params[0].value:
            params[1].enabled = "True"
        else:
            params[1].enabled = "False"
            
        if params[1].value == 1:
            params[1].enabled = "True"
        else:
            params[2].value = ""
            params[2].enabled = "False"
        
        if params[1].value == 1:
            desc = arcpy.Describe(params[0].value)
            fields = desc.fields
            featurefieldList = [field.name for field in fields
                                if field.type in ["String", "Integer", "SmallInteger"]]
                                
            params[2].enabled = "True"
            params[2].filter.type = "ValueList"
            params[2].filter.list = featurefieldList
            
        if params[2].value:
            params[3].enabled = "True"
        else:
            params[3].value = ""
            params[3].enabled = "False"
            
        if params[2].value:
            field_select = params[2].value
            arcpy.Frequency_analysis(params[0].value, "in_memory\\field_freq", field_select)
            
            
            for field in fields:
                if field.name == field_select:
                    type = field.type
                    if type in ("Integer", "SmallInteger"):
                        where = '"{}" IS NOT NULL'.format(field_select)
                    elif type == "String":
                        where = '"{}" <> \'\' and "{}" IS NOT NULL'.format(
                            field_select, field_select)

            featurevalueList = [row[0] for row in arcpy.da.SearchCursor(
                "in_memory\\field_freq", [field_select], where)]        
                    
            featurevalueList.sort()
            
            params[3].enabled = "True"
            params[3].filter.type = "ValueList"
            params[3].filter.list = featurevalueList

        #param5 - Create Map Document
        if not params[5].value == 1:
            params[6].enabled = "False"
            params[8].enabled = "False"
            params[9].enabled = "False"
        else:
            params[6].enabled = "True"
            params[8].enabled = "True"
            params[9].enabled = "True"
            
        templateList = []
        
        #Hard Code MAP TEMPLATES
        path = r'T:\CO\GIS\gistools\tools\Cultural\Templates\Map_Templates'
        
        #Hard Code MAP TEMPLATES
        dirList = os.listdir(path)
        
        for fname in dirList:
            if "mxd" in fname:
                templateList.append(fname)
                
        params[6].filter.type = "ValueList"
        params[6].filter.list = templateList
        
        #param6 - map Template
        if not params[6].altered:
            params[6].value = templateList[0]
            
        #param7 - CR Report Number
        if not params[7].altered:
            params[7].value = "CR-RG-18-xxx X"
            
        #param8 - CR Report Name
        if not params[8].altered:
            params[8].value = ""

        #param9 - Author
        if not params[9].altered:
            params[9].value = ""

        #param10 - SHPO Layers and values
        if not params[10].value == 1:
            params[11].enabled = "False"
            params[12].enabled = "False"
        else:
            params[11].enabled = "True"
            params[12].enabled = "True"
            
        #param11 - Select Layer
        populateList = ["Site", "Survey"]
        params[11].filter.type = "ValueList"
        params[11].filter.list = populateList
        if not params[11].altered:
            params[11].value = "Survey"
            
        return


    def updateMessages(self, params):    
        #param10 - Populate SHPO Layers
        if params[10].value == 1 and params[11].value == "Survey":
            surveyString = str(params[12].ValueAsText)
            surveyCount = surveyString.count(".")
            if not surveyCount == 2:
                params[12].setErrorMessage(
                    "Survey ID must be in format xx.LM.xxXXX (e.g. MC.LM.NR999)")
                
        if params[10].value == 1 and params[11].value == "Site":
            siteString = str(params[12].ValueAsText)
            siteCount = siteString.count(".")
            if not siteCount >= 1:
                params[12].setErrorMessage(
                    "Site ID must be in format 5xx.XXXX (e.g. 5FN.9999)")                             
        return


    def execute(self, params, messages):
        blast_my_cache()

        #Hard Codes
        Sites               = r'T:\ReferenceState\CO\CorporateData\data_sharing\Sites.lyr'
        Surveys             = r'T:\ReferenceState\CO\CorporateData\data_sharing\Survey.lyr'
        DEM                 = r'T:\ReferenceState\CO\CorporateData\topography\dem\Elevation 10 Meter Zunits Feet.lyr'
        mxd                 = arcpy.mapping.MapDocument(
                              r'T:\CO\GIS\gistools\tools\Cultural\Templates\Map_Templates\{}'\
                              .format(params[6].valueAsText))
        sections            = r'T:\ReferenceState\CO\CorporateData\cadastral\Sections.lyr'
        GCDB                = r'T:\ReferenceState\CO\CorporateData\cadastral\Survey Grid.lyr'
        counties            = r'T:\ReferenceState\CO\CorporateData\admin_boundaries\County Boundaries.lyr'
        quad                = r'T:\ReferenceState\CO\CorporateData\cadastral'\
                              r'\24k USGS Quad Index.lyr'
        shpoSurveyTarget    = r'T:\CO\GIS\gistools\tools\Cultural\Templates\SHPO_Templates\survey_ply_tmp.shp'
        shpoSiteTarget      = r'T:\CO\GIS\gistools\tools\Cultural\Templates\SHPO_Templates\site_ply_tmp.shp'

        try:
            date_time_stamp = re.sub('[^0-9]', '', str(datetime.datetime.now())[5:16])

            #Identify Workspace
            baseName = params[4].valueAsText
            envPath = os.path.dirname(baseName) + '\\' + '{}_env_data'.format(
                      os.path.basename(baseName))
            if not os.path.exists(envPath) and params[14].value == 1:
                os.mkdir(envPath)
                
            # Create the logger
            report_path = baseName + r'_Report.txt'
            logfile_path = r'T:\CO\GIS\gistools\tools\Cultural\z_logs\logfile.txt'
            logger = pyt_log(report_path, logfile_path)

            # Start logging
            logger.log_all("Arch General Functions "+str(datetime.datetime.now()))
            logger.log_report("_"*120+"\n")
            logger.log_all("Running environment: Python - {}\n".format(sys.version))
            logger.log_all("User: "+user+"\n")
            logger.log_all("Output Location:\n\t" + params[4].valueAsText+'\n')

            for i, param in enumerate(params):
                logger.logfile("param {} - {}: {}".format(i, param.displayName, param.valueAsText))
            
            #Hold Param0 and check for multi-part and select case values where appropriate
            # If case selct
            if params[1].value == 1:
                desc = arcpy.Describe(params[0].value)
                fields = desc.fields
                
                for field in fields:
                    if field.name == params[2].valueAsText:
                        ftype = field.type
                        field = params[2].value
                        field_select = '"'+params[2].valueAsText+'"'
                        title = params[3].valueAsText
                        if ftype in ("Integer", "SmallInteger"):
                            where = field_select + ' = ' + title
                        elif ftype == "String":
                            where = field_select + ' = ' + "'" + title + "'"
                        break
                    
                arcpy.MakeFeatureLayer_management(params[0].value, "in_memory\\selected", where)
                
            else:  # take everything
                arcpy.MakeFeatureLayer_management(params[0].value, "in_memory\\selected")
            
            polyParts = int(arcpy.GetCount_management("in_memory\\selected").getOutput(0))
            logger.logfile('poly parts: {}'.format(polyParts))
            
            if polyParts >1:
                arcpy.Dissolve_management("in_memory\\selected", 'in_memory\\Poly')    
            
            else:
                arcpy.MakeFeatureLayer_management("in_memory\\selected", 'in_memory\\Poly')
                
            poly = 'in_memory\\Poly'

            arcpy.MakeFeatureLayer_management(Sites, "in_memory\\siteLayer")
            arcpy.MakeFeatureLayer_management(Surveys, "in_memory\\surveyLayer")
            
            # Peel back input polygon 10 meters to prevent extranneous boundary overlap - i.e. PLSS
            # If poly(s) is too small and gets erased, keep original(s)
            arcpy.MakeFeatureLayer_management(poly, "in_memory\\polycopy")
            arcpy.PolygonToLine_management("in_memory\\polycopy", "in_memory\\polylines")
            arcpy.Buffer_analysis("in_memory\\polylines", "in_memory\\polybuffer", 10)
            
            arcpy.Erase_analysis(
                "in_memory\\polycopy", "in_memory\\polybuffer", "in_memory\\PLSSpoly")
            
            inResult  = int(arcpy.GetCount_management("in_memory\\polycopy").getOutput(0))
            outResult = int(arcpy.GetCount_management("in_memory\\PLSSpoly").getOutput(0))
            logger.logfile('Trim - in: {}'.format(inResult))
            logger.logfile('Trim - out: {}'.format(outResult))
            
            if not inResult == outResult:
                arcpy.Delete_management("in_memory\\PLSSpoly")
                arcpy.MakeFeatureLayer_management(poly, "in_memory\\PLSSpoly")
                
            arcpy.Delete_management("in_memory\\polycopy")
            arcpy.Delete_management("in_memory\\polylines")
            arcpy.Delete_management("in_memory\\polybuffer")
          
            #Clip elevation
            if params[14].value == 1:
                param14Name = envPath+"\\dem"
                arcpy.Clip_management(
                    in_raster=DEM, out_raster=param14Name, in_template_dataset=poly,
                    clipping_geometry="ClippingGeometry")
                
            #Map Processes 
            if params[5].value == 1:
                logger.log_all("Calculating PLSS Info\n")
                
                df1 = arcpy.mapping.ListDataFrames(mxd)[0]
                df2 = arcpy.mapping.ListDataFrames(mxd)[1]

                now = datetime.datetime.now()
                
                #Intersect survey poly, GCDB Survey Grid, counties, and quad
                arcpy.Intersect_analysis(
                    ["in_memory\\PLSSpoly", GCDB, counties, quad], "in_memory\\survey", "NO_FID")
                
                arcpy.Frequency_analysis("in_memory\\survey", "in_memory\\County", "COUNTY")
                
                countyList = [str(row[0]).title() + " County"
                              for row in arcpy.da.SearchCursor("in_memory\\County", ["COUNTY"])]
                            
                arcpy.Frequency_analysis("in_memory\\survey", "in_memory\\Quad", "QUAD_NAME")
                
                quadList =   [str(row[0]).title() + " 7.5'" 
                              for row in arcpy.da.SearchCursor("in_memory\\Quad", ["QUAD_NAME"])]

                logger.log_all('Counties:\n{}\n'.format(countyList))
                logger.log_all('Quads:\n{}\n'.format(quadList))
                
                #Select sections that intersect a feature class
                arcpy.MakeFeatureLayer_management(sections, "in_memory\\sections_join")
                arcpy.SelectLayerByLocation_management(
                    "in_memory\\sections_join", "INTERSECT", "in_memory\\PLSSpoly")                
                
                #Get unique PLSSID's
                arcpy.Frequency_analysis(
                    "in_memory\\sections_join", "in_memory\\PLSSID_freq", ["PLSSID"])              
                
                PLSSIDlist = \
                    [row[0] for row in arcpy.da.SearchCursor("in_memory\\PLSSID_freq", ["PLSSID"])]
                meridianlist_orig = [PLSSID[2:4]for PLSSID in PLSSIDlist] 
              
                #Message the console
                infolist=[]
                meridianlist_unique = list(set(meridianlist_orig))
                
                for meridian in meridianlist_unique:
                    if meridian == "06":
                        meridian_name = "6th Principal Meridian"
                    elif meridian == "31":
                        meridian_name = "Ute Principal Meridian"
                    elif meridian == "23":
                        meridian_name = "New Mexico Principal Meridian"
                        
                    infolist.append(meridian_name)
                    
                    for PLSSID in PLSSIDlist:
                        if PLSSID[2:4] == meridian:
                            twnsp = str(int(PLSSID[5:7]))
                            twnspdir = PLSSID[8]
                            plss_range = str(int(PLSSID[9:12]))
                            rangedir = PLSSID[13]                         
                            plss_where = 'PLSSID = {}'.format(PLSSID)

                            trDesc = "T. {} {}., R. {} {}.".format(
                                twnsp, twnspdir, plss_range, rangedir)
                         
                            infolist.append(trDesc)
                            
                legalDesc = '\n'.join(infolist)
                logger.log_all(legalDesc+'\n')

                #Pull and populate elevation in feet
                logger.log_all("Creating Map Document")
                arcpy.FeatureToPoint_management(poly, "in_memory\\centroid", "CENTROID")
                arcpy.sa.ExtractValuesToPoints(
                    "in_memory\\centroid", DEM, "in_memory\\centValue", "NONE", "VALUE_ONLY")
                
                ElePoint = str([row[0] for row
                                in arcpy.da.SearchCursor("in_memory\\centValue", "RASTERVALU")])
                
                Elevation = ElePoint.split(".")
                
                Elev1 = Elevation[0]
                ElevPrint = Elev1[1:]

                #Update report elements with dict - k: element_ID, v: text                
                map_elements = {"ProjectID" : params[7].valueAsText,
                                "Title"     : params[8].valueAsText,
                                "Author"    : params[9].valueAsText,
                                "Date"      : str(now.month)+"\\"+str(now.day)+"\\"+str(now.year),
                                "Location"  : legalDesc,
                                "County"    : "\n".join(countyList),
                                "Quad"      : "\n".join(quadList),
                                "Elevation" : ElevPrint+" feet"}
                                                                
                for item in arcpy.mapping.ListLayoutElements(mxd):
                    # Not all will be in dict..
                    try:
                        ePX = item.elementPositionX
                        item.text = map_elements[item.name]
                        item.elementPositionX = ePX
                    except: pass
                    
                #Identify and select intersected counties
                arcpy.MakeFeatureLayer_management(counties, "in_memory\\NewCounty")       
                arcpy.SelectLayerByLocation_management(
                    "in_memory\\NewCounty", "INTERSECT", "in_memory\\PLSSpoly") 
                    
                arcpy.MakeFeatureLayer_management("in_memory\\NewCounty", "Inset_Cty")
                
                countyLayer = baseName+"_InsetCounty.lyr"
                arcpy.SaveToLayerFile_management("Inset_Cty", countyLayer)
                addLayer = arcpy.mapping.Layer(countyLayer)
                arcpy.mapping.AddLayer(df2, addLayer, "TOP")
                arcpy.Delete_management(countyLayer)
                
                #Add Layers
##                polyLyr = arcpy.mapping.Layer(params[0].valueAsText)
##                polyLyr = arcpy.mapping.Layer(poly)
##                arcpy.mapping.AddLayer(df1, polyLyr, "TOP")
                
                #Set visible layers
                for item in arcpy.mapping.ListLayers(mxd, "", df1):
                    if item.supports("VISIBLE") == "True":
                        item.visible = "False"

                arcpy.RefreshActiveView
                arcpy.RefreshTOC           
                
                #Save as new mxd
                saveName = ((params[4].valueAsText)+".mxd")
                mxd.saveACopy(saveName)
                
                #set scale and pan to extent
                mxd = arcpy.mapping.MapDocument(saveName)
                df1 = arcpy.mapping.ListDataFrames(mxd)[0]
                surDesc = arcpy.Describe(poly)
                newExtent = surDesc.extent
                df1.extent = newExtent
                df1.scale = 24000
                mxd.save()
            
            #SHPO Processes
            if params[10].value == 1:
                shpoName = str(params[12].ValueAsText)
                
                #Survey
                if params[11].ValueAsText == "Survey":
                    logger.log_all("Populating SHPO Survey Layer")
                    shpoName2 = shpoName.replace(".", "_")
                    output = os.path.dirname(params[4].valueAsText)+"\\"+shpoName2+"_SHPO.shp"
                    
                    #Clear template for appending
                    arcpy.DeleteRows_management(shpoSurveyTarget)
                    
                    #Append feature
                    arcpy.Append_management(poly, shpoSurveyTarget, "NO_TEST")
                    arcpy.MakeFeatureLayer_management(shpoSurveyTarget, 'in_memory\\outPoly')
                    
                    #Calculate geometry (area, perimeter, acres)
                    arcpy.CalculateField_management(
                        'in_memory\\outPoly', "AREA", "!shape.area@SQUAREMETERS!", "PYTHON_9.3")
                    
                    arcpy.CalculateField_management(
                        'in_memory\\outPoly', "PERIMETER", "!shape.length@METERS!", "PYTHON_9.3")
                    
                    arcpy.CalculateField_management(
                        'in_memory\\outPoly', "ACRES", "!shape.area@ACRES!", "PYTHON_9.3")
                    
                    #Calculate location (x, y)
                    arcpy.CalculateField_management(
                        'in_memory\\outPoly', "X","!SHAPE.CENTROID.X!", "PYTHON_9.3")
                    
                    arcpy.CalculateField_management(
                        'in_memory\\outPoly', "Y","!SHAPE.CENTROID.Y!", "PYTHON_9.3")
                    
                    #Add other info (DOC_, CONF, ZONE, AGENCY_, SOURCE, DATE)
                    fields = ('DOC_', 'CONF','ZONE', 'AGENCY_', 'DATE')   
                    
                    with arcpy.da.UpdateCursor("in_memory\\outpoly", fields) as cursor:
                        for row in cursor:
                            row[0] = params[12].valueAsText 
                            row[1] = "HC"
                            row[2] = "13"
                            row[3] = params[7].valueAsText
                            row[4] = datetime.datetime.now()
                            
                            cursor.updateRow(row)   
                            
                    #output shape and save
                    arcpy.CopyFeatures_management("in_memory\\outpoly", output)
                    
                    #Clear template for future use
                    arcpy.DeleteRows_management(shpoSurveyTarget)
                    
                #Site
                if params[11].ValueAsText == "Site":
                    logger.log_all("Populating SHPO Site Layer")
                    shpoName2 = shpoName.replace(".", "_")
                    output = os.path.dirname(params[4].valueAsText)+"\\"+shpoName2+"_SHPO.shp"
                    
                    #Clear template for appending
                    arcpy.DeleteRows_management(shpoSiteTarget)
                    
                    #Append feature
                    arcpy.Append_management(poly, shpoSiteTarget, "NO_TEST")
                    arcpy.MakeFeatureLayer_management(shpoSiteTarget, 'in_memory\\outPoly')
                    
                    #Calculate geometry (area, perimeter, acres)
                    arcpy.CalculateField_management(
                        'in_memory\\outPoly', "AREA", "!shape.area@SQUAREMETERS!", "PYTHON_9.3")
                    
                    arcpy.CalculateField_management(
                        'in_memory\\outPoly', "PERIMETER", "!shape.length@METERS!", "PYTHON_9.3")
                    
                    arcpy.CalculateField_management(
                        'in_memory\\outPoly', "ACRES", "!shape.area@ACRES!", "PYTHON_9.3")
                    
                    #Calculate location (x, y)
                    arcpy.CalculateField_management(
                        'in_memory\\outPoly', "X","!SHAPE.CENTROID.X!", "PYTHON_9.3")
                    
                    arcpy.CalculateField_management(
                        'in_memory\\outPoly', "Y","!SHAPE.CENTROID.Y!", "PYTHON_9.3")
                    
                    #Add other info (SITE_, BND_CMPLT, DATE, LINEAR, ZONE, CONF)  
                    fields = ('SITE_','BND_CMPLT','DATE', 'LINEAR', 'ZONE', 'CONF')
                    
                    with arcpy.da.UpdateCursor("in_memory\\outpoly", fields)as cursor:
                        for row in cursor:
                            row[0] = params[12].valueAsText
                            row[1] = "Y"
                            row[2] = datetime.datetime.now()
                            row[3] = '0'
                            row[4] = "13"
                            row[5] = "HC"
                            
                            cursor.updateRow(row)    
                            
                    #output shape and save
                    arcpy.CopyFeatures_management("in_memory\\outpoly", output)
                    
                    #Clear template for future use
##                    arcpy.DeleteRows_management(shpoSiteTarget)
            
            #Output PLSS
            if params[13].value == 1:
                logger.log_all("Writing PLSS Location Data")
                
                #Intersect survey poly, GCDB Survey Grid, counties, and quad
                
                arcpy.Intersect_analysis(
                    ["in_memory\\PLSSpoly", GCDB, counties, quad], "in_memory\\survey", "NO_FID")
                
                freqFields = ["PLSSID","FRSTDIVNO","QQSEC"]
                arcpy.Frequency_analysis("in_memory\\survey", "in_memory\\Loc", freqFields)
                
                fieldNames = ["PM", "TWN", "RNG", "SEC", "QQ1", "QQ2"]
                
                plss = defaultdict(list)
                
                with arcpy.da.SearchCursor("in_memory\\Loc", freqFields) as cursor:
                    for row in cursor:
                        pm  = str(row[0])[2:4]
                        twn = str(row[0])[5:7]+str(row[0])[8]
                        rng = str(row[0])[10:12]+str(row[0])[13]
                        sec = str(row[1])
                        qq1 = str(row[2])[0:2]
                        qq2 = str(row[2])[2:4]
                        
                        plss[pm +'-'+ twn +'-'+ rng +'-'+ sec].append([qq1, qq2])
                            
                print_list = []

                for section in plss.keys():
                    splits = section.split('-')
                    
                    if len(plss[section]) == 16:
                        splits.extend(['ENTIRE', 'SECTION'])
                        print_list.append(splits)
                        continue
                        
                    quarters = defaultdict(list)
                    
                    for quarter in plss[section]:
                        quarters[quarter[1]].append(quarter)
                        
                    for quarter, q_list in quarters.items():
                        if len(q_list) == 4:
                            splits_ = copy.copy(splits)
                            splits_.extend(['ENTIRE', quarter])
                            print_list.append(splits_)
                        else: 
                            for q in q_list:
                                splits_ = copy.copy(splits)
                                splits_.extend(q)
                                print_list.append(splits_)

                out_csv = baseName +'_PLSS.csv'                     

                df = pd.DataFrame(print_list, columns=fieldNames)

                # sort into proper (non-lexigraphic) order
                df.PM  = df.PM.astype(int)
                df.SEC = df.SEC.astype(int)
                df['TWN_v'] = df.TWN.str[:-1]
                df['TWN_d']= df.TWN.str[-1]
                df['RNG_v'] = df.RNG.str[:-1]
                df['RNG_d']= df.RNG.str[-1]

                for col in df.columns:
                    try:
                        df[col] = df[col].astype(int)
                    except: pass

                df.sort(['PM', 'TWN_v', 'TWN_d', 'RNG_v', 'RNG_d', 'SEC', 'QQ2', 'QQ1'],
                        inplace=True)
                df.drop(['TWN_v', 'TWN_d', 'RNG_v', 'RNG_d'], axis=1, inplace=True)
                df.to_csv(out_csv, index=False)
            
            #Output DEM
            if params[14].value == 1:
                logger.log_all("Clipping Digital Elevation Model")
                # Not really - already did it 
                
            #Output Sites
            if params[15].value == 1:
                logger.log_all("Clipping Sites")
                param20Name = baseName+"_Sites.shp"
                
                arcpy.SelectLayerByLocation_management(
                    "in_memory\\siteLayer", "INTERSECT", poly, "", "NEW_SELECTION")
                
                siteResult=int(arcpy.GetCount_management("in_memory\\siteLayer").getOutput(0)) 
                
                if siteResult == 0:
                    logger.log_all(
                        "####"+'\n'+'\n'+"There are no sites within this polygon"+'\n'+'\n'+"####")          
                else:
                    arcpy.CopyFeatures_management("in_memory\\siteLayer", param20Name) 
                    
            #Output Surveys
            if params[16].value == 1:
                logger.log_all("Clipping Surveys")
                param21Name = baseName+"_Surveys.shp"
                
                arcpy.SelectLayerByLocation_management(
                    "in_memory\\surveyLayer", "INTERSECT", poly, "", "NEW_SELECTION")
                
                surveyResult=int(arcpy.GetCount_management("in_memory\\surveyLayer").getOutput(0))           
                
                if surveyResult == 0:
                    logger.log_all(
                      "####"+'\n'+'\n'+"There are no surveys within this polygon"+'\n'+'\n'+"####")
                else:
                    arcpy.Clip_analysis("in_memory\\surveyLayer", poly, param21Name)
                    try:
                        arcpy.AddField_management(param21Name, "ACRES", "DOUBLE", 15, 2)
                    except:
                        try:
                            arcpy.CalculateField_management( param21Name, "ACRES", "!shape.area@ACRES!", "PYTHON_9.3")
                        except:
                            pass
            #Lit
            if params[17].value == 1:
                # sites
                arcpy.SelectLayerByLocation_management(
                    in_layer="in_memory\\siteLayer",
                    overlap_type='WITHIN_A_DISTANCE',
                    select_features=poly,
                    search_distance='1 Mile',
                    selection_type='NEW_SELECTION')
                
                if int(arcpy.GetCount_management("in_memory\\siteLayer").getOutput(0)) <> 0:
                    arcpy.TableToTable_conversion("in_memory\\siteLayer",
                        os.path.dirname(params[4].valueAsText),
                        os.path.basename(params[4].valueAsText+"_Sites_lit.csv"))
                
                # surveys
                arcpy.SelectLayerByLocation_management(
                    in_layer="in_memory\\surveyLayer",
                    overlap_type='WITHIN_A_DISTANCE',
                    select_features=poly,
                    search_distance='1 Mile',
                    selection_type='NEW_SELECTION')
                
                if int(arcpy.GetCount_management("in_memory\\surveyLayer").getOutput(0)) <> 0:
                    arcpy.TableToTable_conversion("in_memory\\surveyLayer",
                        os.path.dirname(params[4].valueAsText),
                        os.path.basename(params[4].valueAsText+"_Surveys_lit.csv"))
            
            # Get the extent and long-axis measurements      
            try:
                # polygon to points
                pnts = [(pnt.X, pnt.Y) for row in arcpy.da.SearchCursor(
                     poly, ["OID@", "SHAPE@"]) for part in row[1] for pnt in part if part]
                # get extent
                extent = extentPnts(pnts)
                # minRect returns [xyPnts, angle, dx, dy]
                min_boundary = minRect(pnts)

                logger.log_all('\nMaximum Extent Points:')
                for ext in extent: logger.log_all(str(ext))
                logger.log_all('\nMinimum Bounding Rectangle Points:')
                for mbrp in min_boundary[0]: logger.log_all(str(mbrp))
                logger.log_all('\nBounding Rectangle Angle: \n{:.5} degrees'.format(min_boundary[1]))
                logger.log_all('\nAxis 1 Extent:\n{:.5}\n'.format(min_boundary[2]))
                logger.log_all('\nAxis 2 Extent:\n{:.5}\n'.format(min_boundary[3]))
            except:
                logger.log_all('Error calculating bounds')
        except:
            logger.log_all(str(traceback.format_exc()))
           
        finally:
            try:
                blast_my_cache()
                
                # sweep up the stray .xmls
                for f in os.listdir(os.path.dirname(baseName)):
                    n, e = os.path.splitext(f)
                    if e == '.xml':
                        os.remove(os.path.join(os.path.dirname(baseName), f))
            except:
                logger.log_all('Error in clean up..')
                logger.log_all(str(traceback.format_exc()))
            
        return                      
