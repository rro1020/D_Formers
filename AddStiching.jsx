function segmentLength(a, b)
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
        return [anchor, leftDirection, rightDirection, pointType];
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
    
    var target_dist = slide % (path_length / n);
    
    var point_list = [];
    var test = 0.0;
    //alert("p length = " + p.length);
    for(var i = 0; i < p.length; i++)
    {
        var next = (i + 1) % p.length;
        if(target_dist == 0.0)
        {
            point_list.push(p[i]);
            target_dist = path_length / n;
            i -= 1;
            continue;
        }
        
        else if(target_dist >= seg_lengths[i])
        {
            target_dist -= seg_lengths[i];
        }
        
        else
        {
            //alert("before split = " + p.length);
            t = target_dist / seg_lengths[i];
            
            tmp = split_curve(p[i], p[next], t);
                
            p[i][2] = tmp[0][2];
            p[next][1] = tmp[2][1];
            p.splice(i + 1,0, tmp[1]);
            
            point_list.push(tmp[1]);

            seg_lengths[i] = segmentLength(tmp[0], tmp[1]);
            seg_lengths.splice(i + 1, 0, segmentLength(tmp[1], tmp[2]));

            target_dist = path_length / n;
            //alert("after split = " + p.length);
            //i += 1;
        }
        test += seg_lengths[i];
        
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
        
        var tangent = [(unit_vec_a[0] + unit_vec_b[0])/2, (unit_vec_a[1] + unit_vec_b[1])/2];
        
        var p_slope = perpendicularSlope(tangent);
        
        var new_pt = computePointAlongLine(p_slope, p[j], o);
        
        p[j][0] = new_pt;
        p[j][1] = [new_pt[0] + dx_l, new_pt[1] + dy_l];
        p[j][2] = [new_pt[0] + dx_r, new_pt[1] + dy_r];
    }
    
    return p;
}

var args = [];

function setUpUI()
{
    var ui = new Window("dialog", "Add Notches to Path");
    
    var g1 = ui.add("group", undefined, "");
    var st1 = g1.add("statictext", undefined, "Number of Notches: ");
    var quantity = g1.add("edittext", undefined, "");
    
    //add more arguements
    
    var closeBtn = ui.add("button", undefined, "Apply");
    
    args = [];
    
    closeBtn.onClick = function(){
        ui.close();
    }
    ui.show();
}

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
            p = parse_path(selected_paths[j]);
            
            var invert = false
            
            if(selected_paths[j].polarity == PolarityValues.NEGATIVE)
            {
                invert = true;
            }
            
            //setUpUI();
            
            var num_nodes = 20;
            var node_length = 20;
            var diameter = 10;
            var slide = 0;
            var input = false;
            
            //alert(selected_paths[j].name);
            
            if(input)
            {
                invert = !invert;
            }
            
            p = add_stiches(p, num_nodes, node_length, diameter, slide, invert);
            
            unparse_path(selected_paths[j], p);
        }
    }
}