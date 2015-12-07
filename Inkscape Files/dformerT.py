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
import dformer

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
    
def addNotches(path, n, offset, slideRatio, angle, document):  
    points = dformer.generatePoints(path, slideRatio, 2 * n)
    #if slideRatio > 0: 
    #    points = dformer.generatePoints(path, 2 * n)
    #else: 
    #    points = dformer.generatePointsWithOffset(path, slideRatio, 2 * n)    
    
    pt_pairs = dformer.findPointPairs(points)
    #inkex.errormsg("pt_pairs = " + str(len(pt_pairs)))
    
    p = cubicsuperpath.parsePath(path.get('d'))[0]
    i = 1
    
    for a, b in pt_pairs:
        slope = dformer.computeSlope(a, b)
        p_slope = dformer.perpendicularSlope(slope)
        
        ac_slope = dformer.rotateSlope(p_slope, 360 - angle)
        bd_slope = dformer.rotateSlope(p_slope, angle)
        
        tmp = offset * math.tan(math.radians(angle))
        dist = math.sqrt(tmp ** 2 + offset ** 2)
        
        c = dformer.computePointAlongLine(ac_slope, a, -dist)
        d = dformer.computePointAlongLine(bd_slope, b, -dist)
        start_pt = p[i - 1][1]
        
        while not(dformer.checkRelativeDifference(start_pt[0], a[0]) and dformer.checkRelativeDifference(start_pt[1], a[1])):
            i += 1
            start_pt = p[i - 1][1]
        begin_idx = i
        end_pt = p[i][1]
        while not(dformer.checkRelativeDifference(end_pt[0], b[0]) and dformer.checkRelativeDifference(end_pt[1], b[1])):
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
        
    path.set('d', cubicsuperpath.formatPath([p]) + 'Z')
   
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
                              
        self.OptionParser.add_option("-q", "--pointsT",
                        action="store", type="int", 
                        dest="pointsT", default=20,
                        help="How Many points of connection are present(teeth)")

        self.OptionParser.add_option("-n", "--offsetT",
                        action="store", type="float", 
                        dest="offsetT", default=5,
                        help="Offset")

        self.OptionParser.add_option("-y", "--paraTooth",
                        action="store", type="float", 
                        dest="paraTooth", default=5,
                        help="Scale Factor (Drawing:Real Length)")
                        
        self.OptionParser.add_option("-l", "--slideT",
                        action="store", type="float", 
                        dest="slideT", default="0",
                        help="Adjust Slide")

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
        if obj_lengths[id_min] > obj_lengths[id_max]:
            id_min = 1
            id_max = 0
            
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
        
        addNotches(obj_nodes[id_min], self.options.pointsT, self.options.offsetT, self.options.slideT, self.options.paraTooth, doc)
        addNotches(obj_nodes[id_max], self.options.pointsT, self.options.offsetT, self.options.slideT, self.options.paraTooth, doc)

        
if __name__ == '__main__':
    e = Length()
    e.affect()


# vim: expandtab shiftwidth=4 tabstop=8 softtabstop=4 fileencoding=utf-8 textwidth=99


# Site for Extension GUI: http://wiki.inkscape.org/wiki/index.php/Extensions