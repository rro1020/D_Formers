
function comparePathsByLength(p1, p2)
{
    if(p1.length < p2.length)
        return -1;
    if(p1.length > p2.length)
        return 1;
    return 0;
}

var smallToBig = true;

function setUpUI()
{
    var ui = new Window("dialog", "Normalize Path Perimeters");
    
    var scaleOpt = ui.add('group', undefined, "Normalize all paths to:");
    scaleOpt.orientation = 'column'
    var btnSL = scaleOpt.add('radiobutton', undefined, 'Largest path length');
    var btnLS = scaleOpt.add('radiobutton', undefined, 'Smallest path length');
    btnSL.value = true;
    
    var closeBtn = ui.add("button", undefined, "Apply");
    
    closeBtn.onClick = function(){
        if(btnSL.value == true)
            smallToBig = true;
        if(btnLS.value == true)
            smallToBig = false;
        ui.close();
    }
    ui.show();
}

if(app.documents.length > 0)
{
    var allPaths = app.activeDocument.pathItems;
    
    if (allPaths.length > 0)
    {
        setUpUI();
        
        var selectedPaths = []
        for(var i = 0; i < allPaths.length; i++)
        {
            if(allPaths[i].selected)
            {
                selectedPaths.push(allPaths[i]);
            }
        }
        
        selectedPaths.sort(comparePathsByLength)
        
        if(!smallToBig)
        {
            selectedPaths.reverse();
        }
        
        targetLength = selectedPaths[selectedPaths.length - 1].length;
        
        for(var i = 0; i < selectedPaths.length; i++)
        {
            scale = (targetLength / selectedPaths[i].length) * 100;
            selectedPaths[i].resize(scale, scale);
        }
    }
}