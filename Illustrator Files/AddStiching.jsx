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
    
    for(var j = 0; j < pnts.length; j++)
    {
        if(pnts[j][3] == PointType.CORNER)
        {
            var prev = (j - 1 + pnts.length) % pnts.length;
            var next = (j + 1 + pnts.length) % pnts.length;
            
            var ctrl_x = (pnts[j][0][0] - pnts[prev][0][0]) * 0.75 + pnts[prev][0][0];
            var ctrl_y = (pnts[j][0][1] - pnts[prev][0][1]) * 0.75 + pnts[prev][0][1]; 
            
            pnts[j][1] = [ctrl_x, ctrl_y];
            
            var ctrl_x_r = (pnts[next][0][0] - pnts[j][0][0]) * 0.25 + pnts[j][0][0];
            var ctrl_y_r = (pnts[next][0][1] - pnts[j][0][1]) * 0.25 + pnts[j][0][1]; 
            
            pnts[j][2] = [ctrl_x_r, ctrl_y_r];
        }
        if(pnts[j][1][0] == pnts[j][0][0] && pnts[j][1][1] == pnts[j][0][1])
        {
            var prev = (j - 1 + pnts.length) % pnts.length;
            var ctrl_x = (pnts[j][0][0] - pnts[prev][0][0]) * 0.75 + pnts[prev][0][0];
            var ctrl_y = (pnts[j][0][1] - pnts[prev][0][1]) * 0.75 + pnts[prev][0][1]; 
            
            pnts[j][1] = [ctrl_x, ctrl_y];
        }
        if(pnts[j][0][0] == pnts[j][2][0] && pnts[j][0][1] == pnts[j][2][1])
        {
            var next = (j + 1 + pnts.length) % pnts.length;
            var ctrl_x_r = (pnts[next][0][0] - pnts[j][0][0]) * 0.25 + pnts[j][0][0];
            var ctrl_y_r = (pnts[next][0][1] - pnts[j][0][1]) * 0.25 + pnts[j][0][1]; 
            
            pnts[j][2] = [ctrl_x_r, ctrl_y_r];
        }
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

function derivative_bcurve(start, ctrl1, ctrl2, end, t)
{
    var d_x = ctr1[0] + 2 * ctrl2[0] * t + 3 * end[0] * t * t;
    var d_y = ctr1[1] + 2 * ctrl2[1] * t + 3 * end[1] * t * t;
    
    return [d_x, d_y];
}

function add_stiches(p, n, o, diameter, slide, invert)
{
    tmp = generate_points(p, n, slide);
    
    //p = tmp[0];
    pts = tmp[1];
    //alert("found pts = " + pts.length);
    
    for(var i = 0; i < pts.length; i++)
    {
        var top = pts[i][0][1] + (diameter/2);
        var left = pts[i][0][0] - (diameter/2);

        app.activeDocument.pathItems.ellipse(top, left, diameter, diameter);
    }
}

function offset_entire_path(p, o, invert)
{
    if(!invert)
        {
            o = -1 * o;
        }

    for(var j = 0; j < p.length; j++)
    {
        var dx_r = p[j][2][0] - p[j][0][0];
        var dy_r = p[j][2][1] - p[j][0][1];

        var dx_l = p[j][1][0] - p[j][0][0];
        var dy_l = p[j][1][1] - p[j][0][1];
        
        var op_dx_r = -1 * dx_r;
        var op_dy_r = -1 * dy_r;
        
        var unit_vec_a = [(dx_l/(Math.sqrt(dx_l * dx_l + dy_l * dy_l))), (dy_l/(Math.sqrt(dx_l * dx_l + dy_l * dy_l)))];
        var unit_vec_b = [(op_dx_r/(Math.sqrt(op_dx_r * op_dx_r + op_dy_r * op_dy_r))), (op_dy_r/(Math.sqrt(op_dx_r * op_dx_r + op_dy_r * op_dy_r)))];
        
        var tangent = [(unit_vec_a[0] + unit_vec_b[0]), (unit_vec_a[1] + unit_vec_b[1])];
        
        var p_slope = perpendicularSlope(tangent);
        var new_pt = computePointAlongLine(p_slope, p[j], o);
        
        var prev = ((j - 1) + p.length) % p.length;
        var next = (j + 1) % p.length;
        
        /* var vec_p_c = [(p[j][1][0] - p[prev][2][0]), (p[j][1][1] - p[prev][2][1])];
        var norm1 = Math.sqrt(((vec_p_c[0] * vec_p_c[0]) + (vec_p_c[1] * vec_p_c[1])));
        var norm2 = Math.sqrt(((dx_l * dx_l) + (dy_l * dy_l)));
        var dot_prod = vec_p_c[0] * dx_l + vec_p_c[1] * dy_l;
        var angle = Math.acos(dot_prod / (norm1 * norm2)) * 0.5;
        var left_ctrl_slope = rotateSlope(vec_p_c, 360 - (angle * (180 / Math.PI)));
        var new_left_ctrl = computePointAlongLine(left_ctrl_slope, [p[j][1]],o);
        
        var vec1 = [(new_left_ctrl[0] - new_pt[0]), (new_left_ctrl[1] - new_pt[1])];
        var norm1 = Math.sqrt(((vec1[0] * vec1[0]) + (vec1[1] * vec1[1])));
        new_left_ctrl = computePointAlongLine([dx_l, dy_l], [new_pt], norm1);
        
        var vec_n_c = [(p[j][2][0] - p[next][1][0]), (p[j][2][1] - p[next][1][1])];
        var norm1 = Math.sqrt(((vec_n_c[0] * vec_n_c[0]) + (vec_n_c[1] * vec_n_c[1])));
        var norm2 = Math.sqrt(((dx_r * dx_r) + (dy_r * dy_r)));
        var dot_prod = vec_n_c[0] * dx_r + vec_n_c[1] * dy_r;
        var angle = Math.acos(dot_prod / (norm1 * norm2));
        var right_ctrl_slope = rotateSlope(vec_n_c, (angle * (180 / Math.PI)));
        var new_right_ctrl = computePointAlongLine(right_ctrl_slope, [p[j][2]],o);
        
        var vec1 = [(new_right_ctrl[0] - new_pt[0]), (new_right_ctrl[1] - new_pt[1])];
        var norm1 = Math.sqrt(((vec1[0] * vec1[0]) + (vec1[1] * vec1[1])));
        new_right_ctrl = computePointAlongLine([dx_r, dy_r], [new_pt], norm1); */
        
        p[j][0] = new_pt;
        p[j][1] = [new_pt[0] + dx_l, new_pt[1] + dy_l];
        p[j][2] = [new_pt[0] + dx_r, new_pt[1] + dy_r];
    }
    
    return p;
}

var args = [];
var skip = false;

function setUpUI()
{
    var ui = new Window("dialog", "Add Stitching to Path");
    
    skip = true;
    
    var g1 = ui.add("group", undefined, "");
    g1.alignment = "column"
    var st1 = g1.add("statictext", undefined, "Number of Holes: ");
    var quantity = g1.add("edittext", undefined, "0");
    quantity.characters = 5;
    
    var g2 = ui.add("group", undefined, "");
    var st2 = g2.add("statictext", undefined, "Diameter of Holes: ");
    var hole_diameter = g2.add("edittext", undefined, "0.00");
    hole_diameter.characters = 5;
    
    var g3 = ui.add("group", undefined, "");
    var st3 = g3.add("statictext", undefined, "Distance from Path: ");
    var dist = g3.add("edittext", undefined, "0.00");
    dist.characters = 5;
    
    var g4 = ui.add("group", undefined, "");
    var st4 = g4.add("statictext", undefined, "Slide along Path: ");
    var slide = g4.add("slider", undefined, "0.00");
    var txt_slide = g4.add("edittext", undefined, "0.00");
    txt_slide.characters = 5;
    
    var g5 = ui.add("group", undefined, "");
    var invert_check = g5.add("checkbox", undefined, "Invert? ");
    invert_check.characters = 5;
    
    var g6 = ui.add("group", undefined, "");
    var closeBtn = g6.add("button", undefined, "Apply");
    var skipBtn = g6.add("button", undefined, "Skip");
    
    slide.onChanging = function(){
        txt_slide.text = "" + (slide.value / 100.0).toFixed(2);
        if(slide.value == 0 || slide.value == 100)
        {
            txt_slide.text += ".00"
        }
    }
    txt_slide.onChanging = function(){
        var v = parseFloat(txt_slide.text);
        slide.value = v * 100;
    }
    closeBtn.onClick = function(){
        args = [];skip = false;
        args.push(parseInt(quantity.text));
        args.push(parseFloat(hole_diameter.text));
        args.push(parseFloat(dist.text));
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
                invert = true;
            }
            
            var p = parse_path(selected_paths[j]);
            var valid = false;
            var num_nodes;
            var node_length;
            var diameter;
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
                    error_msg += "Number of Holes: Please enter a positive integer.\n";
                }
                node_length = args[2];
                if((!isInteger(node_length) && !isFloat(node_length)) || !isPositive(node_length))
                {
                    valid = false;
                    error_msg += "Length from path: Please enter a non-negative number.\n";
                }
                diameter = args[1];
                if((!isInteger(diameter) && !isFloat(diameter)))
                {
                    valid = false;
                    error_msg += "Diameter of Holes: Please enter a non-negative number.\n";
                }
                slide = args[3];
                if((!isInteger(slide) && !isFloat(slide)) || !isPositive(slide))
                {
                    valid = false;
                    error_msg += "Slide: Please enter a non-negative number.\n";
                }
                input = args[4];
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
            
            add_stiches(p, num_nodes, node_length, diameter, slide);
            
            var new_p = parse_path(selected_paths[j]);
            new_p = offset_entire_path(p, node_length, invert);
            
            unparse_path(selected_paths[j], new_p);
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