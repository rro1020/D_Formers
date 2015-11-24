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

function generate_points(p, n)
{
    var seg_lengths = segmentLengths(p);
    
    var path_length = 0;
    for(var i in seg_lengths) {path_length += seg_lengths[i];}
    
    var target_dist = 0.0;
    
    var point_list = [];

    for(var i = 0; i < p.length; i++)
    {
        alert(target_dist);
        var next = (i + 1) % p.length;
        if(target_dist == 0.0)
        {
            point_list.push(p[i]);
            target_dist = path_length / n;
        }
        
        else if(target_dist >= seg_lengths[i])
        {
            target_dist -= seg_lengths[i];
        }
        
        else
        {
            t = target_dist / seg_lengths[i];
            
            tmp = split_curve(p[i], p[next], t);
                
            p[i] = tmp[0];
            p[next] = tmp[2];
            p.splice(next,0, tmp[1]);
            
            point_list.push(tmp[1]);

            seg_lengths[i] = segmentLength(tmp[0], tmp[1]);
            seg_lengths.splice(next, 0, segmentLength(tmp[1], tmp[2]));

            target_dist = path_length / n;
        }
    }
    
    return [p, point_list]
}

if(app.documents.length > 0)
{
    var allPaths = app.activeDocument.pathItems;
    
    if (allPaths.length > 0)
    {
        var selectedPaths = []
        for(var i = 0; i < allPaths.length; i++)
        {
            if(allPaths[i].selected)
            {
                alert(allPaths[i].length)
                p = parse_path(allPaths[i]);
                
                tmp = generate_points(p, 8);
                p = tmp[0];
                alert(tmp[1])
                
                unparse_path(allPaths[i], p);
            }
        }
    }
}