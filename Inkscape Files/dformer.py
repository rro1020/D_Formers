#!/usr/bin/env python 
'''
This extension module contains all helper functions for all 
dformer connectors.  This project was assigned to us for 
Texas A&M University's CSCE 482 Course.  

Copyright (C) 2010 Alvin Penner
Copyright (C) 2006 Georg Wiora
Copyright (C) 2006 Nathan Hurst
Copyright (C) 2005 Aaron Spike, aaron@ekips.org

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version. 

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

''' 

# standard library
import locale
# local library
import inkex
import simplestyle
import simpletransform
import cubicsuperpath
import bezmisc
import addnodes
import math

from simpletransform import fuseTransform 
from simplepath import parsePath, formatPath
from simplestyle import * 


inkex.localize()
locale.setlocale(locale.LC_ALL, '')

# third party
try:
    import numpy
except:
    inkex.errormsg(_("Failed to import the numpy modules. These modules are required by this extension. Please install them and try again.  On a Debian-like system this can be done with the command, sudo apt-get install python-numpy."))
    exit()

mat_area   = numpy.matrix([[  0,  2,  1, -3],[ -2,  0,  1,  1],[ -1, -1,  0,  2],[  3, -1, -2,  0]])
mat_cofm_0 = numpy.matrix([[  0, 35, 10,-45],[-35,  0, 12, 23],[-10,-12,  0, 22],[ 45,-23,-22,  0]])
mat_cofm_1 = numpy.matrix([[  0, 15,  3,-18],[-15,  0,  9,  6],[ -3, -9,  0, 12],[ 18, -6,-12,  0]])
mat_cofm_2 = numpy.matrix([[  0, 12,  6,-18],[-12,  0,  9,  3],[ -6, -9,  0, 15],[ 18, -3,-15,  0]])
mat_cofm_3 = numpy.matrix([[  0, 22, 23,-45],[-22,  0, 12, 10],[-23,-12,  0, 35],[ 45,-10,-35,  0]])

def numsegs(csp):
    return sum([len(p)-1 for p in csp])
def interpcoord(v1,v2,p):
    return v1+((v2-v1)*p)
def interppoints(p1,p2,p):
    return [interpcoord(p1[0],p2[0],p),interpcoord(p1[1],p2[1],p)]
def pointdistance((x1,y1),(x2,y2)):
    return math.sqrt(((x2 - x1) ** 2) + ((y2 - y1) ** 2))
def bezlenapprx(sp1, sp2):
    return pointdistance(sp1[1], sp1[2]) + pointdistance(sp1[2], sp2[0]) + pointdistance(sp2[0], sp2[1])
def tpoint((x1,y1), (x2,y2), t = 0.5):
    return [x1+t*(x2-x1),y1+t*(y2-y1)]
def cspbezsplit(sp1, sp2, t = 0.5):
    m1=tpoint(sp1[1],sp1[2],t)
    m2=tpoint(sp1[2],sp2[0],t)
    m3=tpoint(sp2[0],sp2[1],t)
    m4=tpoint(m1,m2,t)
    m5=tpoint(m2,m3,t)
    m=tpoint(m4,m5,t)
    return [[sp1[0][:],sp1[1][:],m1], [m4,m,m5], [m3,sp2[1][:],sp2[2][:]]]
def cspbezsplitatlength(sp1, sp2, l = 0.5, tolerance = 0.001):
    bez = (sp1[1][:],sp1[2][:],sp2[0][:],sp2[1][:])
    t = bezmisc.beziertatlength(bez, l, tolerance)
    return cspbezsplit(sp1, sp2, t)
def cspseglength(sp1,sp2, tolerance = 0.001):
    bez = (sp1[1][:],sp1[2][:],sp2[0][:],sp2[1][:])
    return bezmisc.bezierlength(bez, tolerance)    
def csplength(csp):
    total = 0
    lengths = []
    for sp in csp:
        lengths.append([])
        for i in xrange(1,len(sp)):
            l = cspseglength(sp[i-1],sp[i])
            lengths[-1].append(l)
            total += l
    return lengths, total

#Defined Functions 
def verifyDocument(document):
    document.set('width', '18in')
    document.set('height', '32in')
    
    margin = inkex.etree.Element(inkex.addNS('rect', 'svg'))
    margin.set('y', '14.940762')
    margin.set('x', '-301.96371')
    margin.set('width', '1022.4807')
    margin.set('height', '1840.9845')
    margin.set('stroke', 'red')
    margin.set('stroke-width', '1')
    margin.set('fill', 'none')
    
    document.append(margin)

def rotateSlope((x, y), degrees):
    degrees = degrees % 360
    radians = math.radians(degrees)
    new_x = x * math.cos(radians) - y * math.sin(radians)
    new_y = x * math.sin(radians) + y * math.cos(radians)
    return new_x, new_y 
    
def computeSlope((x1, y1), (x2, y2)):
    return x2 - x1, y2 - y1
    
def negateSlope((dx, dy)): 
    return -dx, -dy
	
def perpendicularSlope((dx, dy)):
    return -dy, dx
    
def computeLineIntercept((x, y), slope): 
    return y - slope * x 

    
def findPointInList((x, y), pointList): 
    for p in pointList: 
        if p[0] == x and p[1] == y:
            index = pointList.index(p)
    return index 

def crossProductDirection((x1, y1), (x2, y2), (x3, y3)):
    vector1 = (x1 - x2), (y1 - y2), 0
    vector2 = (x3 - x2), (y3 - y2), 0 
    direction = True
    
    product = vector1[1] * vector2[2] - vector1[2] * vector2[1], vector1[2] * vector2[0] - vector1[0] * vector2[2], vector1[0] * vector2[1] - vector1[1] * vector2[0]
    if (product[0] > 0 or product[1] > 0 or product[2] > 0): 
        direction = False 
    return direction 
    
def getDirection(index, pointList):
    #need to consider cases with a small amount of points 
    if (index + 1 == len(pointList) - 1): 
        firstPoint = pointList[index]
        secondPoint = pointList[index + 1]
        thirdPoint = pointList[0] 
        index = -2 
    elif (index == len(pointList) - 1): 
        firstPoint = pointList[index]
        secondPoint = pointList[0]
        thirdPoint = pointList[1] 
        index = -1 
    else: 
        firstPoint = pointList[index]
        secondPoint = pointList[index + 1] 
        thirdPoint = pointList[index + 2]
    
    midPoint = getCenterPoint(pointList) 
    inkex.errormsg("1stPoint = " + str(firstPoint) + "\t2ndPoint = " + str(secondPoint)) 
    inkex.errormsg("3rdPoint = " + str(thirdPoint) + "\tmidPoint = " + str(midPoint)) 
    
    isInward = True 
    if (firstPoint[1] > midPoint[1]): 
        #starting point is above midpoint
        inkex.errormsg("firstPoint is above midPoint")
        if (firstPoint[0] > secondPoint[0] and secondPoint[0] > thirdPoint[0]): 
            inkex.errormsg("if")
            isInward = True 
        elif (firstPoint[0] > secondPoint[0] and secondPoint[0] < thirdPoint[0]):
            #redo with different points 
            inkex.errormsg("elif")
            isInward = getDirection(index + 3, pointList)
        else: 
            inkex.errormsg("else")
            isInward = False 
    else: 
        #starting point is at or below midpoint 
        inkex.errormsg("firstPoint is below midPoint")
        if (firstPoint[0] < secondPoint[0] and secondPoint[0] < thirdPoint[0]): 
            inkex.errormsg("if")
            isInward = True 
        elif (firstPoint[0] < secondPoint[0] and secondPoint[0] > thirdPoint[0]): 
            #redo with different points 
            inkex.errormsg("elif")
            isInward = getDirection(index + 3, pointList)
        else: 
            inkex.errormsg("else")
            isInward = False 
    
    return isInward 
            
            
def getDirectionOfPoint(p, pointList):
    #index = findPointInList((x, y), pointList)
    index = pointList.index(p)
    pAfter = pointList[index]
    pBefore = pointList[index]
    
    if (index == 0): 
        #inkex.errormsg("index is 0")
        pBefore = pointList[len(pointList) - 1]
        pAfter = pointList[index + 1]
    elif (index == len(pointList) - 1):
        #inkex.errormsg("index is at the end")
        pBefore = pointList[index - 1] 
        pAfter = pointList[0] 
    else: 
        pBefore = pointList[index - 1]
        pAfter = pointList[index + 1]
    
    inkex.errormsg("Index = " + str(index)) 
    inkex.errormsg("PBefore = " + str(pBefore)) 
    inkex.errormsg("P = " + str(p)) 
    inkex.errormsg("PAfter = " + str(pAfter)) 
    
    direction = crossProductDirection(pBefore, p, pAfter)
    
    return direction 
    
def computeMidpointIntersection((x1, y1), (x2, y2), (dx1, dy1), (dx2, dy2)):
    # slope1 = 0 
    # slope2 = 0 
    # if (dy1 != 0): 
        # slope1 = dx1/dy1 
    # if (dy2 != 0): 
        # slope2 = dx2/dy2
    if (dx1 == 0): 
        b2 = computeLineIntercept((x2, y2), (dy2/dx2))
        x = x1 
        y = (dy2/dx2) * x + b2
    elif (dx2 == 0):
        b1 = computeLineIntercept((x1, y1), (dy1/dx1))
        x = x2 
        y = (dy1/dx1) * x + b1
    elif (dx1 == 0 and dx2 == 0): 
        y = 0 
    else: 
        b1 = computeLineIntercept((x1, y1), (dy1/dx1))
        b2 = computeLineIntercept((x2, y2), (dy2/dx2))
    
        x = (b1 - b2) / ((dy2/dx2) - (dy1/dx1)) 
        y = (dy1/dx1) * x + b1
    
    #inkex.errormsg("slope1 = " + str((dy1/dx1)) + "\tslope2 = " + str((dy2/dx2)))
    #inkex.errormsg("b1 = " + str(b1) + "\tb2 = " + str(b2))
    
    return (x, y)
    
def computePointAlongLine((dx, dy), (x, y), distance):
    #compute unit vector
    #inkex.errormsg("dx = " + str(dx) + "  dy = " + str(dy))
    norm = (dx ** 2 + dy ** 2) ** 0.5
    
    #if slope is undefined, new point is the x_coord + distance  
    new_x = 0 
    new_y = 0 
    if (norm == 0): 
        #Subtraction only works for vertical slopes on the right side
        #isInside = getDirection((dx, dy))
        new_x = x - distance 
        new_y = y
    else: 
        u_x = dx / norm
        u_y = dy / norm
    
        new_x = x + distance * u_x
        new_y = y + distance * u_y
    
    return new_x, new_y

def computePointAlongLine2((dx, dy), (x, y), distance, dirOfPoint):
    #compute unit vector
    #inkex.errormsg("dx = " + str(dx) + "  dy = " + str(dy))
    norm = (dx ** 2 + dy ** 2) ** 0.5
    
    #if slope is undefined, new point is the x_coord + distance  
    new_x = 0 
    new_y = 0 
    if (norm == 0): 
        #Subtraction only works for vertical slopes on the right side
        inkex.errormsg("Direction = " + str(dirOfPoint))
        if (dirOfPoint):
            new_x = x + distance 
            new_y = y
        else: 
            new_x = x - distance 
            new_y = y 
    else: 
        u_x = dx / norm
        u_y = dy / norm
    
        new_x = x + distance * u_x
        new_y = y + distance * u_y
    
    return new_x, new_y
    
def replaceSegmentWith(path, start, end, subpath):
    new_path = ""
    
    split_path = path.replace("L", "C").split("C")
    
    curves = split_path[1:]
    new_path += split_path[0] + curves[0]
    i = 1
    while i < len(curves):
        points = curves[i - 1].split()
        inkex.errormsg(str(start) + " " + str(end) + " " + curves[i] + '\n')
        if start == (float(points[-2]), float(points[-1])):
            points = curves[i].split()
            while end != (points[-2], points[-1]):
                i += 1
                points = curves[i].split()
            new_path += subpath
        else:
            new_path += "C" + curves[i]
        
        i += 1
    
    return new_path

def findPointPairs(points):
    pairs = []
    
    for i, elem in enumerate(points):
        if i % 2 != 0:
            pairs += [(points[i - 1], points[i])]
    
    return pairs

def slidePointPairs(points): 
    firstPoint = points[0] 
    for i, j in enumerate(points):
        if (i+1 > len(points) - 1): 
            points[i] = firstPoint
        else: 
            points[i] = points[i+1]
    return findPointPairs(points)


#Based off the Bounding Box Method 
def getCenterPoint(points): 
    px_min = points[0] 
    px_max = points[1] 
    py_min = points[0] 
    py_max = points[1]
    
    for p in points: 
        if p[0] < px_min[0]:
            px_min = p 
        elif p[0] > px_max[0]: 
            px_max = p 

        if p[1] < py_min[1]: 
            py_min = p 
        elif p[1] > py_max[1]:
            py_max = p 
        
    x_slope = computeSlope(px_min, px_max) 
    y_slope = computeSlope(py_min, py_max) 
    
    return computeMidpointIntersection(px_min, py_min, x_slope, y_slope)
    
    
def getMidPoint(p1, p2): 
    x = (p1[0] + p2[0])/2 
    y = (p1[1] + p2[1])/2

    return (x, y)
    
def checkRelativeDifference(val1, val2, tolerance=0.9999):
    return (min(val1,val2)/max(val1,val2)) >= tolerance

def rescalePath(path, n, oriCenterPoint, offset, diameter, document): 
    #oriPointList = generatePoints(path, n)
    #oriCenterPoint = getCenterPoint(oriPointList)
    
    #originalOrigin = originParse(path)
    
    perimeter = read_stored_info("pathlength", path)
    path.set('transform', 'scale(' + str(1 + offset/100) + ' ' + str(1 + offset/100) +')')
    fuseTransform(path) 
    
    newPath = cubicsuperpath.parsePath(path.get('d'))
    #newOrigin = originParse(path)
    
    #mat = simpletransform.composeParents(node, [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]])
    #node.set('d', cubicsuperpath.formatPath(parsedPath))
    #simpletransform.applyTransformToPath(mat, parsedPath)
    
    slengths, stotal = csplength(newPath)
    #newPoints = generatePoints(newPath, n)
    
    
    #save the path length and segment lengths in the document 
    path.set('pathlength', str(stotal))
    tmp = ''
    for slen in slengths[0][::-1]:
        tmp += str(slen) + ' '
    tmp = tmp[:-1] #remove the space at the end
    path.set('segmentlengths', tmp)
    
    newPointList = generatePoints(path, n)
    newCenterPoint = getCenterPoint(newPointList)
    drawCircle(newCenterPoint, diameter/2, document)
    
    for p in newPointList: 
        drawCircle(p, diameter/2, document)
    
    transformOrigin = [(oriCenterPoint[0] - newCenterPoint[0]), (oriCenterPoint[1] - newCenterPoint[1])]
    inkex.errormsg("transformOrigin = " + str(transformOrigin))
    path.set('transform', 'translate(' + str(transformOrigin[0]) + ' ' + str(transformOrigin[1]) +')')
    fuseTransform(path)
    
def read_stored_info(type, obj):
    if type == 'pathlength':
        return float(obj.get(type))
    elif type == 'segmentlengths':
        tmp = obj.get(type).split(' ')
        tmp = map(float, tmp)
        #raise Exeception(str(tmp))
        return tmp
    return None

def generatePoints(obj, offset, n):
    pointList = []  
    targetDist = offset * read_stored_info("pathlength", obj)/ n 
    segmentLengths = [0.0] + read_stored_info("segmentlengths", obj)[::-1]
    
    mat = simpletransform.composeParents(obj, [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]])
    p = cubicsuperpath.parsePath(obj.get('d'))
    simpletransform.applyTransformToPath(mat, p)
    
    segments = p[0]
    new_path = segments
    
    i = 0
    count = 0
    total = 0
    
    while i < len(segments): 
        if segmentLengths[i] == targetDist:
            pointList += [segments[i][1]]
            targetDist = read_stored_info("pathlength", obj)/ n

        elif targetDist > segmentLengths[i] or segmentLengths[i] == 0:
            targetDist -= segmentLengths[i]

        else:
            t = targetDist / float(segmentLengths[i]) 
            
            t1 = segments[i-1][1]
            t2 = segments[i][1]
            
            pathList = addnodes.cspbezsplitatlength(segments[i-1], segments[i], t,tolerance=0.00001)
            pointList += [pathList[1][1]]
            prev = segments[:i - 1]
            next = segments[i + 1:]
            segments = prev + pathList + next
            
            len1 = targetDist
            len2 = segmentLengths[i] - targetDist
            
            if abs((len1 + len2) - segmentLengths[i]) >= 0.0001:
                raise Exception("There is an issue with the bezier split function. (" + str(len1 + len2) + " vs. " + str(segmentLengths[i]) + ")")
                
            
            #inkex.errormsg("Split Point = " + str(pointList[-1]))
            #inkex.errormsg("len1 = " + str(len1) + "\nlen2 = " + str(len2) + " \nsegmentLengths[i] = " + str(segmentLengths[i]))
            
            prev = segmentLengths[:i]
            next = segmentLengths[i+1:]
            #inkex.errormsg("Pathlist = " + str(pathList))
            segmentLengths = prev + [len1, len2] + next
        
            #obj.set('d', cubicsuperpath.formatPath([segments]))
            #segments = cubicsuperpath.parsePath(obj.get('d'))[0]
            targetDist = read_stored_info("pathlength", obj)/n
            
            
            count += 1
        if(len(pointList) == n):
            break;
        i += 1
    
    pointList += [segments[0][1]]
    
    if abs(sum(segmentLengths) - read_stored_info("pathlength", obj)) >= 0.001:
        raise Exception("Internal Error: The total length of the new path does not equal the original path.")
    obj.set('d', cubicsuperpath.formatPath([segments[:-1]])+'Z')
    obj.set('segmentlengths', str(segmentLengths))
      
    return pointList 

def drawColoredCircle(tarPoint, radius, doc, color):
#draws a circle at a point
    circ = inkex.etree.Element(inkex.addNS('circle', 'svg'))
    circ.set('cx', str(tarPoint[0]))
    circ.set('cy', str(tarPoint[1]))
    circ.set('r', str(radius))
    circ.set('stroke', 'black')
    circ.set('stroke-width', '1')
    circ.set('fill', color)
    
    doc.append(circ)
    
def drawCircle(tarPoint, radius, doc):
#draws a circle at a point
    circ = inkex.etree.Element(inkex.addNS('circle', 'svg'))
    circ.set('cx', str(tarPoint[0]))
    circ.set('cy', str(tarPoint[1]))
    circ.set('r', str(radius))
    circ.set('stroke', 'black')
    circ.set('stroke-width', '1')
    circ.set('fill', 'white')
    
    doc.append(circ)

def printValue(value, self):
#creates text object based on given value (debugging purposes)   
    what = value

    # Get access to main SVG document element and get its dimensions.
    svg = self.document.getroot()
    # or alternatively
    # svg = self.document.xpath('//svg:svg',namespaces=inkex.NSS)[0]

    # Again, there are two ways to get the attibutes:
    width  = self.unittouu(svg.get('width'))
    height = self.unittouu(svg.attrib['height'])

    # Create a new layer.
    layer = inkex.etree.SubElement(svg, 'g')
    layer.set(inkex.addNS('label', 'inkscape'), 'Hello %s Layer' % (what))
    layer.set(inkex.addNS('groupmode', 'inkscape'), 'layer')
    # Create text element
    text = inkex.etree.Element(inkex.addNS('text','svg'))
    text.text = 'Value: %s' % (what)

    # Set text position to center of document.
    text.set('x', str(width / 2))
    text.set('y', str(height / 2))
    # Center text horizontally with CSS style.
    style = {'text-align' : 'center', 'text-anchor': 'middle'}
    text.set('style', formatStyle(style))

    # Connect elements together.
    layer.append(text)
    
def originParse(pathTarget):
# find the starting point of a path
    dString = pathTarget.get('d')
    
    dToken = 1
    xString = ""
    yString = ""
    symbols = "0123456789.-"
    
    while dString[dToken] not in symbols:
        dToken += 1
    
    while dString[dToken] in symbols:
        xString += dString[dToken]
        dToken += 1
    
    while dString[dToken] not in symbols:
        dToken += 1
    
    while dString[dToken] in symbols:
        yString += dString[dToken]
        dToken += 1
    
    xf = float(xString)
    yf = float(yString)

    return [xf,yf]    
   
class Length(inkex.Effect):
    def __init__(self):
        inkex.Effect.__init__(self)
        self.OptionParser.add_option("--type",
                        action="store", type="string", 
                        dest="type", default="length",
                        help="Type of measurement") 
                        
        self.OptionParser.add_option("-u", "--unit",
                        action="store", type="string", 
                        dest="unit", default="mm",
                        help="The unit of the measurement")
                        
        self.OptionParser.add_option("-p", "--points",
                        action="store", type="int", 
                        dest="points", default=20,
                        help="How Many points of connection are present")
                        
        self.OptionParser.add_option("-o", "--offset",
                        action="store", type="float", 
                        dest="offset", default=5,
                        help="Offset")
                        
        self.OptionParser.add_option("-x", "--paraStitch",
                        action="store", type="float", 
                        dest="paraStitch", default=5,
                        help="Scale Factor (Drawing:Real Length)")

        self.OptionParser.add_option("-y", "--paraTooth",
                        action="store", type="float", 
                        dest="paraTooth", default=5,
                        help="Scale Factor (Drawing:Real Length)")
                        
        self.OptionParser.add_option("-z", "--paraLeaf",
                        action="store", type="float", 
                        dest="paraLeaf", default=5,
                        help="Scale Factor (Drawing:Real Length)")
                        
        self.OptionParser.add_option("-r", "--radioScale",
                        action="store", type="string", 
                        dest="radioScale", default="B2S",
                        help="Radio Button Scaling")
                        
        self.OptionParser.add_option("-b", "--radioBuild",
                        action="store", type="string", 
                        dest="radioBuild", default="none",
                        help="Radio Button Builb Type")
                        
        self.OptionParser.add_option("--tab",
                        action="store", type="string", 
                        dest="tab", default="desc",
                        help="The selected UI-tab when OK was pressed")
                        
        self.OptionParser.add_option("--measurehelp",
                        action="store", type="string", 
                        dest="measurehelp", default="",
                        help="dummy")

    def effect(self):
        # get number of digits
        #prec = int(self.options.precision)
        scale = self.unittouu('1px')    # convert to document units
        factor = 1.0
        doc = self.document.getroot()
        if doc.get('viewBox'):
            [viewx, viewy, vieww, viewh] = doc.get('viewBox').split(' ')
            factor = self.unittouu(doc.get('width'))/float(vieww)
            if self.unittouu(doc.get('height'))/float(viewh) < factor:
                factor = self.unittouu(doc.get('height'))/float(viewh)
            factor /= self.unittouu('1px')
        
        # loop over all selected paths
        obj_lengths = []  
        obj_ids = []
        obj_nodes = []
        for id, node in self.selected.iteritems():
            if node.tag == inkex.addNS('path','svg'):
                mat = simpletransform.composeParents(node, [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]])
                p = cubicsuperpath.parsePath(node.get('d'))
                node.set('d', cubicsuperpath.formatPath(p))
                simpletransform.applyTransformToPath(mat, p)
                factor *= scale/self.unittouu('1'+self.options.unit)
                if self.options.type == "length":
                    slengths, stotal = csplength(p)
                    
                    #save the path length and segment lengths in the document 
                    node.set('pathlength', str(stotal))
                    tmp = ''
                    for slen in slengths[0][::-1]:
                        tmp += str(slen) + ' '
                    tmp = tmp[:-1] #remove the space at the end
                    node.set('segmentlengths', tmp)
                    
                    # self.group = inkex.etree.SubElement(node.getparent(),inkex.addNS('text','svg'))
                    obj_lengths += [stotal]
                    obj_ids += [id]
                    obj_nodes += [node] 
        
        id_min = 0
        id_max = 1
        minOrigin = []
        maxOrigin = []
        if obj_lengths[id_min] > obj_lengths[id_max]:
            id_min = 1
            id_max = 0
            
        minOrigin = originParse(obj_nodes[id_min])
        maxOrigin = originParse(obj_nodes[id_max])
        
        if self.options.radioScale == "B2S":
            ratio = obj_lengths[id_min] / obj_lengths[id_max]
            obj_ori = []
            ori_trans = []      
            #obj_nodes[id_max].set('transform', 'scale(' + str(ratio * 2) + ' ' + str(ratio * 2) +')')
            obj_nodes[id_max].set('transform', 'scale(' + str(ratio) + ' ' + str(ratio) +')')      
            #obj_nodes[id_min].set('transform', 'scale(' + str(2) + ' ' + str(2) +')')      
            fuseTransform(obj_nodes[id_max])
            #fuseTransform(obj_nodes[id_min])
            
            obj_ori = originParse(obj_nodes[id_max])
            ori_trans = [(maxOrigin[0] - obj_ori[0]), (maxOrigin[1] - obj_ori[1])]
            obj_nodes[id_max].set('transform', 'translate(' + str(ori_trans[0]) + ' ' + str(ori_trans[1]) +')')
            # obj_ori = originParse(obj_nodes[id_min])
            # ori_trans = [(minOrigin[0] - obj_ori[0]), (minOrigin[1] - obj_ori[1])]
            # obj_nodes[id_min].set('transform', 'translate(' + str(ori_trans[0]) + ' ' + str(ori_trans[1]) +')')
            fuseTransform(obj_nodes[id_max])
            #fuseTransform(obj_nodes[id_min])
            
        elif self.options.radioScale == "S2B":
            ratio = obj_lengths[id_max] / obj_lengths[id_min]
            obj_ori = []
            ori_trans = []
            obj_nodes[id_min].set('transform', 'scale(' + str(ratio) + ' ' + str(ratio) +')')
            fuseTransform(obj_nodes[id_min])
            
            obj_ori = originParse(obj_nodes[id_min])
            ori_trans = [(minOrigin[0] - obj_ori[0]), (minOrigin[1] - obj_ori[1])]
            obj_nodes[id_min].set('transform', 'translate(' + str(ori_trans[0]) + ' ' + str(ori_trans[1]) +')')
            fuseTransform(obj_nodes[id_min])
        
        #verifyDocument(doc)
        obj_lengths = []  
        obj_ids = []
        obj_nodes = []
        for id, node in self.selected.iteritems():
            if node.tag == inkex.addNS('path','svg'):
                mat = simpletransform.composeParents(node, [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]])
                p = cubicsuperpath.parsePath(node.get('d'))
                simpletransform.applyTransformToPath(mat, p)
                node.set('d', cubicsuperpath.formatPath(p))
                factor *= scale/self.unittouu('1'+self.options.unit)
                if self.options.type == "length":
                    slengths, stotal = csplength(p)
                    
                    #save the path length and segment lengths in the document 
                    node.set('pathlength', str(stotal))
                    tmp = ''
                    for slen in slengths[0][::-1]:
                        tmp += str(slen) + ' '
                    tmp = tmp[:-1] #remove the space at the end
                    node.set('segmentlengths', tmp)
                    
                    # self.group = inkex.etree.SubElement(node.getparent(),inkex.addNS('text','svg'))
                    obj_lengths += [stotal]
                    obj_ids += [id]
                    obj_nodes += [node] 
                # Format the length as string
                # lenstr = locale.format("%(len)25."+str(prec)+"f",{'len':round(stotal*factor*self.options.scale,prec)}).strip()
        
        points = []
        
        if self.options.tab == "\"stitch\"":
            addStitching(obj_nodes[id_min], self.options.points, self.options.offset, self.options.paraStitch, doc)
            addStitching(obj_nodes[id_max], self.options.points, self.options.offset, self.options.paraStitch, doc)
        elif self.options.tab == "\"tooth\"":
            addNotches(obj_nodes[id_min], self.options.points, self.options.offset, False, self.options.paraTooth,doc)
            addNotches(obj_nodes[id_max], self.options.points, self.options.offset, True, self.options.paraTooth,doc)
        elif self.options.tab == "\"leaf\"": 
            addLeaves(obj_nodes[id_min], self.options.points, self.options.offset, self.options.paraLeaf, doc)
            addLeaves(obj_nodes[id_max], self.options.points, self.options.offset, self.options.paraLeaf, doc)
        


    
if __name__ == '__main__':
    e = Length()
    e.affect()


# vim: expandtab shiftwidth=4 tabstop=8 softtabstop=4 fileencoding=utf-8 textwidth=99


# Site for Extension GUI: http://wiki.inkscape.org/wiki/index.php/Extensions