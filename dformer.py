#!/usr/bin/env python 
'''
This extension module can measure arbitrary path and object length
It adds text to the selected path containing the length in a given unit.
Area and Center of Mass calculated using Green's Theorem:
http://mathworld.wolfram.com/GreensTheorem.html

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

TODO:
 * should use the standard attributes for text
 * Implement option to keep text orientation upright 
    1. Find text direction i.e. path tangent,
    2. check direction >90 or <-90 Degrees
    3. rotate by 180 degrees around text center
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
    
def perpendicularSlope((dx, dy)):
    return -dy, dx
    
def computePointAlongLine((dx, dy), (x, y), distance):
    #compute unit vector
    norm = (dx ** 2 + dy ** 2) ** 0.5
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

def check_relative_difference(val1, val2, tolerance=0.9999):
    return (min(val1,val2)/max(val1,val2)) >= tolerance
 
def addStitching(path, n, offset, diameter, document): 
    points = generatePoints(path, n) 
    pointPairs = findPointPairs(points)
    
    # for p1, p2, p3, p4 in points: 
        # #x, y = bezierslopeatt(((bx0,by0),(bx1,by1),(bx2,by2),(bx3,by3)), p) 
        # x1, y1 = bezierslopeatt((p1, p2, p3, p4), p1) 
        # x2, y2 = bezierslopeatt((p1, p2, p3, p4), p2)
        # x3, y3 = bezierslopeatt((p1, p2, p3, p4), p3)
        # x4, y4 = bezierslopeatt((p1, p2, p3, p4), p4)
        
        # newPoint1 = computePointAlongLine((x1, y1), p1, offset)
        # newPoint2 = computePointAlongLine((x2, y2), p2, offset)
        # newPoint3 = computePointAlongLine((x3, y3), p3, offset)
        # newPoint4 = computePointAlongLine((x4, y4), p4, offset)
        # drawCircle(newPoint1, diameter/2, document)
        # drawCircle(newPoint2, diameter/2, document)
        # drawCircle(newPoint3, diameter/2, document)
        # drawCircle(newPoint4, diameter/2, document)
    
    for p1, p2 in pointPairs:
        pSlope = perpendicularSlope(computeSlope(p1, p2)) 
        newPoint1 = computePointAlongLine(pSlope, p1, offset)
        newPoint2 = computePointAlongLine(pSlope, p2, offset)
        drawCircle(newPoint1, diameter/2, document)
        drawCircle(newPoint2, diameter/2, document)
    
    # Attempt to rescale path by offset 
    # originalOrigin = originParse(path)
    # path.set('transform', 'scale(' + str(offset) + ' ' + str(offset) +')')
    # newOrigin = originParse(path)
    # transformOrigin = [(originalOrigin[0] - newOrigin[0]), (originalOrigin[1] - newOrigin[1])]
    # path.set('transform', 'translate(' + str(transformOrigin[0]) + ' ' + str(transformOrigin[1]) +')')
    #useTransform(path)
    
def addNotches(path, n, offset, angle, document):
    points = generatePoints(path, 2 * n)
    pt_pairs = findPointPairs(points)
    #inkex.errormsg("pt_pairs = " + str(pt_pairs))
    
    p = cubicsuperpath.parsePath(path.get('d'))[0]
    i = 1
    
    for a, b in pt_pairs:
        slope = computeSlope(a, b)
        p_slope = perpendicularSlope(slope)
        
        ac_slope = rotateSlope(p_slope, 360 - angle)
        bd_slope = rotateSlope(p_slope, angle)
        
        tmp = offset * math.tan(math.radians(angle))
        dist = math.sqrt(tmp ** 2 + offset ** 2)
        
        c = computePointAlongLine(ac_slope, a, -dist)
        d = computePointAlongLine(bd_slope, b, -dist)
        #inkex.errormsg("a = " + str(cubicsuperpath.parsePath('m 1,2C 3,4 5,6 7,8L 9,10L11,12L13,14C15,16 17,18 19,20')[0]) + "\n")
        #raise Exception('IGNORE ME!!!')
        #inkex.errormsg("before i = " + str(i) + "\n")
        start_pt = p[i - 1][1]
        
        while not(check_relative_difference(start_pt[0], a[0]) and check_relative_difference(start_pt[1], a[1])):
            #inkex.errormsg("rejected: a = " + str(list(a)) + " " + str(start_pt) + "\n")
            #inkex.errormsg("rejected: b = " + str(list(b)) + " " + str(end_pt) + "\n")
            i += 1
            start_pt = p[i - 1][1]
        begin_idx = i
        end_pt = p[i][1]
        while not(check_relative_difference(end_pt[0], b[0]) and check_relative_difference(end_pt[1], b[1])):
            i += 1
            end_pt = p[i][1]
        end_idx = i
        
        before = p[:begin_idx]
        after = p[end_idx+1:]
        
        before[-1][1] = list(a)
        before[-1][2] = list(a)
        line_to_c = [list(c)] * 3
        line_to_d = [list(d)] * 3
        line_to_b = [list(b)] * 3
        
        p = before + [line_to_c] + [line_to_d] + [line_to_b] + after
        
        #inkex.errormsg("final: a = " + str(list(a)) + " " + str(start_pt) + "\n")
        #inkex.errormsg("final: b = " + str(list(b)) + " " + str(end_pt) + "\n")
        
       
        
        # line_ac = inkex.etree.Element(inkex.addNS('line', 'svg'))
        # line_ac.set('x1', str(a[0]))
        # line_ac.set('y1', str(a[1]))
        # line_ac.set('x2', str(c[0]))
        # line_ac.set('y2', str(c[1]))
        # line_ac.set('stroke', 'red')
        
        # line_cd = inkex.etree.Element(inkex.addNS('line', 'svg'))
        # line_cd.set('x1', str(c[0]))
        # line_cd.set('y1', str(c[1]))
        # line_cd.set('x2', str(d[0]))
        # line_cd.set('y2', str(d[1]))
        # line_cd.set('stroke', 'red')
        
        # line_db = inkex.etree.Element(inkex.addNS('line', 'svg'))
        # line_db.set('x1', str(d[0]))
        # line_db.set('y1', str(d[1]))
        # line_db.set('x2', str(b[0]))
        # line_db.set('y2', str(b[1]))
        # line_db.set('stroke', 'red')
        
        # layer = inkex.etree.SubElement(document, 'g')
        # layer.set(inkex.addNS('label','inkscape'), 'New Layer')
        
        # layer.append(line_ac)
        # layer.append(line_cd)
        # layer.append(line_db)
        
        # document.append(layer)
        
        #replaceSegmentWith(path.get('d'), a, b, '')
    path.set('d', cubicsuperpath.formatPath([p]) + 'Z')
    
def read_stored_info(type, obj):
    if type == 'pathlength':
        return float(obj.get(type))
    elif type == 'segmentlengths':
        tmp = obj.get(type).split(' ')
        tmp = map(float, tmp)
        #raise Exeception(str(tmp))
        return tmp
    return None

# generate points but points are not correct!     
def generatePoints(obj, n):
    pointList = []  
    targetDist = read_stored_info("pathlength", obj)/ n 
    segmentLengths = [0.0] + read_stored_info("segmentlengths", obj)[::-1]
    
    mat = simpletransform.composeParents(obj, [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]])
    p = cubicsuperpath.parsePath(obj.get('d'))
    simpletransform.applyTransformToPath(mat, p)
    
    segments = p[0]
    new_path = segments
    
    i = 0
    count = 0
    total = 0
    #inkex.errormsg("Len(segments) = " + str(len(segments)))
    
    while i < len(segments) - 1: 
        #inkex.errormsg("Target Distance = " + str(targetDist) + "  " + str(i))
        #inkex.errormsg(str(len(segments)) + " " + str(len(segmentLengths)))
        #inkex.errormsg("segmentLengths = " + str(segmentLengths))
        #inkex.errormsg("Curve = " + str((segments[i - 1][1], segments[i - 1][2], segments[i][0], segments[i][1])))
        #inkex.errormsg("record length = " + str(segmentLengths[i]))
        #inkex.errormsg("measured length = " + str(bezmisc.bezierlength((segments[i - 1][1], segments[i - 1][2], segments[i][0], segments[i][1]),tolerance=0.00001)))
        if segmentLengths[i] == targetDist:
            pointList += [segments[i][1]]
            #inkex.errormsg("0.001 confirmed\n added point " + str(pointList[-1]))
            targetDist = read_stored_info("pathlength", obj)/ n
            #raise Exception(str(pointList))

        elif targetDist > segmentLengths[i] or segmentLengths[i] == 0:
            #raise Exception(str(targetDist))
            targetDist -= segmentLengths[i]

        else:
            t = targetDist / float(segmentLengths[i]) 
            
            #inkex.errormsg("t = " + str(t))
            #inkex.errormsg("targetPt = " + str(bezmisc.bezierpointatt((segments[i - 1][1], segments[i - 1][2], segments[i][0], segments[i][1]),t)))
             
            t1 = segments[i-1][1]
            t2 = segments[i][1]
            
            #sub1, sub2 = bezmisc.beziersplitatt((segments[i-1][-2], segments[i-1][-1], segments[i][0], segments[i][1]), t)
            pathList = addnodes.cspbezsplitatlength(segments[i-1], segments[i], t,tolerance=0.00001)
            pointList += [pathList[1][1]]
            #inkex.errormsg("added point " + str(pointList[-1]))
            prev = segments[:i - 1]
            next = segments[i + 1:]
            segments = prev + pathList + next
            #inkex.errormsg("Segments[\n" + str(segments) + "\n]")
            
            #len1 = bezmisc.bezierlength((pathList[0][1][:], pathList[0][2][:], pathList[1][0][:], pathList[1][1][:]),tolerance=0.00001)
            #len2 = bezmisc.bezierlength((pathList[1][1][:], pathList[1][2][:], pathList[2][0][:], pathList[2][1][:]),tolerance=0.00001)
            
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
            
            
            #inkex.errormsg("prev = " + str(prev))
            #inkex.errormsg("next = " + str(next))
            #inkex.errormsg("SegmentLengths[\n" + str(segmentLengths) + "\n]")
            #inkex.errormsg('\n\n')
            count += 1
        i += 1
    
    pointList += [segments[0][1]]
    
    if abs(sum(segmentLengths) - read_stored_info("pathlength", obj)) >= 0.001:
        #inkex.errormsg(str(sum(segmentLengths) - read_stored_info("pathlength", obj)))
        raise Exception("Internal Error: The total length of the new path does not equal the original path.")
    #raise Exception(str(segments)+ "\n\n" + str(new_path))
    #raise Exception(pointList)
    #inkex.errormsg(str(pointList))
    #raise Exception(cubicsuperpath.formatPath([segments[:-1]]))
    obj.set('d', cubicsuperpath.formatPath([segments[:-1]]))
    #inkex.errormsg("segmentLengths = " + str(len(segmentLengths)))
    #inkex.errormsg(str(len(segments)))
    #inkex.errormsg(str(count))
    
    #if len(pointList) != n:
    #    raise Exception("Internal Error: The algorithm did not find the required number of points (" + str(len(pointList)) + " out of " + str(n) + ").")
    
    return pointList 
    
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
                # Format the length as string
                # lenstr = locale.format("%(len)25."+str(prec)+"f",{'len':round(stotal*factor*self.options.scale,prec)}).strip()
        
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
            obj_nodes[id_max].set('transform', 'scale(' + str(ratio) + ' ' + str(ratio) +')')
            fuseTransform(obj_nodes[id_max])
            
            obj_ori = originParse(obj_nodes[id_max])
            ori_trans = [(maxOrigin[0] - obj_ori[0]), (maxOrigin[1] - obj_ori[1])]
            obj_nodes[id_max].set('transform', 'translate(' + str(ori_trans[0]) + ' ' + str(ori_trans[1]) +')')
            fuseTransform(obj_nodes[id_max])
            
        elif self.options.radioScale == "S2B":
            ratio = obj_lengths[id_max] / obj_lengths[id_min]
            obj_ori = []
            ori_trans = []
            obj_nodes[id_min].set('transform', 'scale(' + str(ratio) + ' ' + str(ratio) +')')
            fuseTransform(obj_nodes[id_min])
            
            obj_ori = originParse(obj_nodes[id_min])
            #drawCircle(obj_ori, self.options.paraStitch, doc)
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
        #points = generatePoints(obj_nodes[id_min], self.options.points)
        #points = generatePoints(obj_nodes[id_max], self.options.points)
        #addNotches(obj_nodes[id_min], self.options.points, self.options.offset, self.options.paraTooth,doc)
        
        #printValue(self.options.tab, self)
        
        if self.options.tab == "\"stitch\"":
            addStitching(obj_nodes[id_min], self.options.points, self.options.offset, self.options.paraStitch, doc)
            addStitching(obj_nodes[id_max], self.options.points, self.options.offset, self.options.paraStitch, doc)
            #printValue(self.options.tab, self)
        elif self.options.tab == "\"tooth\"":
            addNotches(obj_nodes[id_min], self.options.points, self.options.offset, self.options.paraTooth,doc)
            #printValue(self.options.tab, self)
            addNotches(obj_nodes[id_max], self.options.points, self.options.offset, self.options.paraTooth,doc)
            #printValue(self.options.tab, self)
        
        #inkex.errormsg(str(points))
        #inkex.errormsg(str(len(points)))
        
        #buildType Radio Button Here
        
    #q = inkex.getElementById(doc, obj_nodes[id_min]); 
        

    # def addCross(self, node, x, y, scale):
        # l = 3*scale         # 3 pixels in document units
        # node.set('d', 'm %s,%s %s,0 %s,0 m %s,%s 0,%s 0,%s' % (str(x-l), str(y), str(l), str(l), str(-l), str(-l), str(l), str(l)))
        # node.set('style', 'stroke:#000000;fill:none;stroke-width:%s' % str(0.5*scale))

    
if __name__ == '__main__':
    e = Length()
    e.affect()


# vim: expandtab shiftwidth=4 tabstop=8 softtabstop=4 fileencoding=utf-8 textwidth=99


# Site for Extension GUI: http://wiki.inkscape.org/wiki/index.php/Extensions