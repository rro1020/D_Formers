var balfax = 0;
var balfbx = 0;
var balfcx = 0;
var balfay = 0;
var balfby = 0;
var balfcy = 0;

function balf(t)
{
    var retval = Math.pow((balfax * (t * t) + balfbx * t + balfcx), 2) + Math.pow((balfay * (t * t) + balfby * t + balfcy), 2);
    return Math.sqrt(retval);
}

function bezierparameterize(start, ctrl_1, ctrl_2, end)
{
    var x0 = start[0];
    var y0 = start[1];
    
    var cx = 3 * (ctrl_1[0] - x0);
    var bx = 3 * (ctrl_2[0] - ctrl_1[0]) - cx;
    var ax = end[0] - x0 - cx - bx;
    
    var cy = 3 * (ctrl_1[1] - y0);
    var by = 3 * (ctrl_2[1] - ctrl_1[1]) - cy;
    var ay = end[1] - y0 - cy - by;
    
    return [ax, ay, bx, by, cx, cy, x0, y0];
}

function simpson(a, b)
{
    var n_limit = 4096;
    var tolerance = 0.00001;
    var n = 2;
    var multiplier = (b - a)/6.0;
    var endsum = balf(a) + balf(b);
    var interval = (b - a)/2.0;
    var asum = 0.0;
    var bsum = balf(a + interval);
    var est1 = multiplier * (endsum + (2.0 * asum) + (4.0 * bsum));
    var est0 = 2.0 * est1;
    
    while(n < n_limit && Math.abs(est1 - est0) > tolerance)
    {
        n = n * 2;
        multiplier /= 2;
        interval /= 2.0;
        asum += bsum;
        bsum = 0.0;
        est0 = est1;
        for(var i = 1; i < n; i += 2)
        {
            bsum += balf(a + (i * interval));
            est1 = multiplier * (endsum + (2.0 * asum) + (4.0 * bsum))
        }
    }
    
    return est1;
}

function segmentLength(a, b)
{
    var tmp = bezierparameterize(a[0], a[2], b[1], b[0]);
    balfax = 3 * tmp[0];
    balfbx = 2 * tmp[2];
    balfcx = tmp[4];
    balfay = 3 * tmp[1];
    balfby = 2 * tmp[3];
    balfcy = tmp[5];
    
    return simpson(0.0,1.0);
}

function bezieratlength(a, b, l)
{
    var tmp = bezierparameterize(a[0], a[2], b[1], b[0]);
    balfax = 3 * tmp[0];
    balfbx = 2 * tmp[2];
    balfcx = tmp[4];
    balfay = 3 * tmp[1];
    balfby = 2 * tmp[3];
    balfcy = tmp[5];
    
    var t = 1.0;
    var tdiv = t;
    var curlen = simpson(0.0, t);
    var targetlen = l * curlen;
    var diff = curlen - targetlen;
    var tolerance = 0.00001;
    
    while(Math.abs(diff) > tolerance)
    {
        tdiv /= 2.0;
        if(diff < 0)
        {
            t += tdiv;
        }
        else
        {
            t -= tdiv;
        }
        curlen = simpson(0.0, t);
        diff = curlen - targetlen;
    }
    
    return t;
}


/* function segmentLength(a, b)
{
    var tmpPath = app.activeDocument.pathItems.add();
    tmpPath.setEntirePath([a[0],b[0]]);
    
    tmpPath.pathPoints[0].pointType = a[3];
    tmpPath.pathPoints[0].leftDirection = a[1];
    tmpPath.pathPoints[0].rightDirection = a[2];
    
    tmpPath.pathPoints[1].pointType = b[3];
    tmpPath.pathPoints[1].leftDirection = b[1];
    tmpPath.pathPoints[1].rightDirection = b[2];
    
    var length = tmpPath.length;
    tmpPath.remove();
    return length;
} */

function segmentLengths(path)
{
    var lengths = [];
    for(var j = 0; j < path.length; j++)
    {
        var next = (j + 1) % path.length;
        lengths.push(segmentLength(path[j], path[next]));
    }
    return lengths;
}

function parse_pnt(pnt)
{
    with(pnt)
    {
        if(pointType == PointType.CORNER) return [anchor, anchor, anchor, PointType.SMOOTH];
        return [anchor, leftDirection, rightDirection, pointType];
    }
}

function parse_path(path)
{
    pnts = [];
    for(var i = 0; i < path.pathPoints.length; i++)
    {
        pnts.push(parse_pnt(path.pathPoints[i]));
    }
    
    return pnts;
}

function unparse_path(path, p_pnts)
{
    anchors = [];
    for(var i = 0; i < p_pnts.length; i++)
    {
        anchors.push(p_pnts[i][0]);
    }
    //alert("anchors = " + anchors);
    path.setEntirePath(anchors);
    
    for(var j = 0; j < path.pathPoints.length; j++)
    {
        path.pathPoints[j].leftDirection = p_pnts[j][1];
        path.pathPoints[j].rightDirection = p_pnts[j][2];
        path.pathPoints[j].pointType = p_pnts[j][3];
    }
}

function split_curve(a, b, t)
{
    var start_pt = a[0];
    var ctrl_1 = a[2];
    var ctrl_2 = b[1];
    var end_pt = b[0];
    
    var c12 = [((ctrl_1[0] - start_pt[0]) * t + start_pt[0]), ((ctrl_1[1] - start_pt[1]) * t + start_pt[1])];
    var c23 = [((ctrl_2[0] - ctrl_1[0]) * t + ctrl_1[0]), ((ctrl_2[1] - ctrl_1[1]) * t + ctrl_1[1])];
    var c34 = [((end_pt[0] - ctrl_2[0]) * t + ctrl_2[0]), ((end_pt[1] - ctrl_2[1]) * t + ctrl_2[1])];
    var c123 = [((c23[0] - c12[0]) * t + c12[0]), ((c23[1] - c12[1]) * t + c12[1])];
    var c234 = [((c34[0] - c23[0]) * t + c23[0]), ((c34[1] - c23[1]) * t + c23[1])];
    var c1234 = [((c234[0] - c123[0]) * t + c123[0]), ((c234[1] - c123[1]) * t + c123[1])];
    
    return [[start_pt,a[1],c12, a[3]], [c1234, c123, c234, PointType.SMOOTH], [end_pt,c34,b[2], b[3]]]
}

function generate_points(p, n, slide)
{
    var seg_lengths = segmentLengths(p);
    
    var path_length = 0;
    for(var i in seg_lengths) {path_length += seg_lengths[i];}
    
    var target_dist = (path_length/n) * slide;
    //alert(target_dist);
    var point_list = [];
    var test = 0.0;
    //alert("p length = " + p.length);
    for(var i = 0; i < p.length + 1; i++)
    {
        i %= p.length;
        var next = (i + 1) % p.length;
        var curr_seg_len = segmentLength(p[i], p[next]);
        //alert("current = " + curr_seg_len);
        if(target_dist == 0.0)
        {
            point_list.push(p[i]);
            target_dist = path_length / n;
            i -= 1;
            continue;
        }
        
        else if(target_dist >= curr_seg_len)
        {
            target_dist -= curr_seg_len;
        }
        
        else
        {
            //alert("before split = " + p.length);
            var t =  bezieratlength(p[i], p[next], target_dist/curr_seg_len);
            
            tmp = split_curve(p[i], p[next], t);
            //alert("original = " + seg_lengths[i]);    
            p[i][2] = tmp[0][2];
            p[next][1] = tmp[2][1];
            p.splice(i + 1,0, tmp[1]);
            
            point_list.push(tmp[1]);

            //seg_lengths[i] = segmentLength(tmp[0], tmp[1]);
            //seg_lengths.splice(i + 1, 0, segmentLength(tmp[1], tmp[2]));
            
            //alert("new = " + (segmentLength(tmp[0], tmp[1]) + segmentLength(tmp[1], tmp[2])));
            
            target_dist = path_length / n;
            //alert("after split = " + p.length);
            //i += 1;
        }
        test += curr_seg_len;
        
        if(point_list.length == n)
            break;
    }
    return [p, point_list]
}

function computeSlope(pt1, pt2)
{
    var x1 = pt1[0][0];
    var y1 = pt1[0][1];
    var x2 = pt2[0][0];
    var y2 = pt2[0][1];
    
    return [(x2 - x1), (y2 - y1)];
}

function perpendicularSlope(slope)
{
    return [(-slope[1]), (slope[0])];
}

function rotateSlope(pts, degrees)
{
    var degrees = degrees % 360;
    var radians = degrees * (Math.PI / 180);
    
    var new_x = pts[0] * Math.cos(radians) - pts[1] * Math.sin(radians);
    var new_y = pts[0] * Math.sin(radians) + pts[1] * Math.cos(radians);

    return [new_x, new_y];
}
    
function findPointPairs(pts)
{
    var pairs = [];
    
    for(var i = 0; i < pts.length; i++)
    {
        if(i % 2 != 0)
        {
            pairs.push([pts[i - 1], pts[i]]);
        }
    }
    
    return pairs;
}
    
function computePointAlongLine(slope, pt, distance)
{
    var norm = Math.sqrt(((slope[0] * slope[0]) + (slope[1] * slope[1])));
    var u_x = slope[0] / norm;
    var u_y = slope[1] / norm;
    
    var new_x = pt[0][0] + distance * u_x;
    var new_y = pt[0][1] + distance * u_y;

    return [new_x, new_y];
}

function check_relative_difference(val1, val2, tolerance)
{
    return ((Math.min(val1, val2) / Math.max(val1, val2)) >= tolerance);
}
    
function add_leaves(p, n, o, slide, invert)
{
    tmp = generate_points(p, n * 2, slide);
    
    p = tmp[0];
    pts = tmp[1];
    //alert("found pts = " + pts.length);
    pair_pts = findPointPairs(pts);
    
    //alert("pair lengths = " + pair_pts.length);
    
    var j = 0;
    
    for(var i = 0; i < pair_pts.length; i++)
    {

        var a = pair_pts[i][0];
        var b = pair_pts[i][1];
        
        var slope = computeSlope(a, b);
        var p_slope = perpendicularSlope(slope);

        var ac_slope = rotateSlope(p_slope, 0);
        var bd_slope = rotateSlope(p_slope, 0);
        
        var tmp = o * Math.tan(0 * (Math.PI / 180));
        var dist = Math.sqrt(tmp * tmp + o * o);
        
        //alert("dist = " + dist);
        
        var polarity = 1;
        
        if(invert)
        {
            polarity = -1;
        }
        
        dist = dist * polarity;
        
        var c = computePointAlongLine(ac_slope, a, dist);
        
        var d = computePointAlongLine(bd_slope, b, dist);
        var midpt = [(c[0] + d[0])/2, (c[1] + d[1])/2];
        
        var a1 = 45;
        var a2 = 315;
        
        if(invert)
        {
            a1 = 315;
            a2 = 45;
        }
        
        var mc_slope = rotateSlope(ac_slope, a1);
        var md_slope = rotateSlope(ac_slope, a2);
        var m_c = computePointAlongLine(mc_slope, a, dist/2);
        var m_d = computePointAlongLine(md_slope, b, dist/2);
        
        var start_pt = p[j];
        //alert("j = " + j + " " + start_pt);
        while(start_pt[0] != a[0])
            {
                j += 1;
                j %= p.length;
                //alert("-j = " + j);
                start_pt = p[j];
            }
            
        var begin_idx = j;
        
        var end_pt = p[j];
        
        while(end_pt[0] != b[0])
            {
                j += 1;
                j %= p.length;
                end_pt = p[j];
            }
            
        var end_idx = j;
        //alert("len = " + p.length);
        var e = end_idx - begin_idx - 1;
        
        var before;
        var after;
        
        if(e < 0)
        {
            before = p.splice(1, begin_idx + 1);
            after = p.splice(0, 1);
        }
        else
        {
            before = p.splice(0, begin_idx + 1);
            after = p.splice(e, p.length);
        }
        
        
        before[before.length - 1][2] = m_c;
        
        var pt_m = [midpt, c, d, PointType.SMOOTH];
        
        after[0][1] = m_d;
        
        //alert("c = " + pt_c);
        
        before.push.apply(before, [pt_m]);
        before.push.apply(before, after);
        
        p = before
        
        //alert("len1 = " + p.length);
    }
    
    return p;
}

var args = [];
var skip = false;

function setUpUI()
{
    var ui = new Window("dialog", "Add Leaves to Path");
    
    skip = true;
    
    var g1 = ui.add("group", undefined, "");
    g1.alignment = "column"
    var st1 = g1.add("statictext", undefined, "Number of leaves: ");
    var quantity = g1.add("edittext", undefined, "0");
    quantity.characters = 5;
    
    var g2 = ui.add("group", undefined, "");
    var st2 = g2.add("statictext", undefined, "Length of each leaf: ");
    var notch_length = g2.add("edittext", undefined, "0.00");
    notch_length.characters = 5;
    
    var g4 = ui.add("group", undefined, "");
    var st4 = g4.add("statictext", undefined, "Slide along path: ");
    var slide = g4.add("slider", undefined, "");
    var txt_slide = g4.add("edittext", undefined, "0.00");
    txt_slide.characters = 5;
    
    var g5 = ui.add("group", undefined, "");
    var invert_check = g5.add("checkbox", undefined, "Invert");
    invert_check.characters = 5;
    
    var g6 = ui.add("group", undefined, "");
    var closeBtn = g6.add("button", undefined, "Apply");
    var skipBtn = g6.add("button", undefined, "Skip");
    
    slide.onChanging = function(){
        txt_slide.text = "" + (slide.value / 100.0).toFixed(2);
    }
    txt_slide.onChanging = function(){
        var v = parseFloat(txt_slide.text);
        slide.value = v * 100;
    }
    closeBtn.onClick = function(){
        args = [];
        skip = false;
        args.push(parseInt(quantity.text));
        args.push(parseFloat(notch_length.text));
        args.push(parseFloat((slide.value/100).toFixed(2)));
        args.push(invert_check.value);
        ui.close();
    }
    skipBtn.onClick = function()
    {
        skip = true;
        ui.close();
    }
    
    ui.show();
}

function isInteger(x) { return Math.floor(x) === x; }
function isFloat(x) { return !!(x % 1); }
function isPositive(x) { return x >= 0.0; }


if(app.documents.length > 0)
{
    var allPaths = app.activeDocument.pathItems;
    
    if (allPaths.length > 0)
    {
        var selected_paths = []
        for(var i = 0; i < allPaths.length; i++)
        {
            if(allPaths[i].selected)
            {
                selected_paths.push(allPaths[i]);
            }
        }
        
        for(var j = 0; j < selected_paths.length; j++)
        {
            
            
            var invert = false
            
            if(selected_paths[j].polarity == PolarityValues.NEGATIVE)
            {
                selected_paths.reverse = true;
            }
            p = parse_path(selected_paths[j]);
            var valid = false;
            var num_nodes;
            var node_length;
            var slide;
            var input;
            
            var first_time = true;
            skip = false;
            var error_msg = "";
            while(!valid && !skip)
            {
                if(!first_time)
                {
                    alert(error_msg);
                }
                error_msg = "";
                setUpUI();
                valid = true;
                num_nodes = args[0];
                if(!isInteger(num_nodes) || !isPositive(num_nodes) || num_nodes == 0)
                {
                    valid = false;
                    error_msg += "Number of Leaves: Please enter a positive integer.\n";
                }
                node_length = args[1];
                if((!isInteger(node_length) && !isFloat(node_length)) || !isPositive(node_length))
                {
                    valid = false;
                    error_msg += "Length of each leaf: Please enter a non-negative number.\n";
                }
                slide = args[2];
                if((!isInteger(slide) && !isFloat(slide)) || !isPositive(slide))
                {
                    valid = false;
                    error_msg += "Slide: Please enter a non-negative number.\n";
                }
                input = args[3];
                first_time = false;
            
            }

            if(skip)
            {
                continue;
            }
            
            if(input)
            {
                invert = !invert;
            }
            
            p = add_leaves(p, num_nodes, node_length,slide, invert);
            
            unparse_path(selected_paths[j], p);
        }
        if(selected_paths.length == 0)
        {
            alert("Please select one or more path objects.")
        }
    }
}
else
{
    alert("Please create a new document.");
}