function app=FigKeyPressApp(figHandle,eventdata,app) %#ok<*INUSD>
    if ~isempty(app.selectTile)
        dropImage = false;
        key = eventdata.Key;
        shifts = app.currShift; % could also have stored this in handles
        switch(key)
            case('w') %want to move top image up
                shifts.upDownOffset = shifts.upDownOffset - 1;
            case('s') %want to move top image down
                shifts.upDownOffset = shifts.upDownOffset + 1;
            case('a') %want to move top image left
                shifts.leftRightOffset = shifts.leftRightOffset - 1;
            case('d') %want to move top image right
                shifts.leftRightOffset = shifts.leftRightOffset + 1;
             case('i') %want to move top image up
                shifts.upDownOffset = shifts.upDownOffset - 10;
            case('k') %want to move top image down
                shifts.upDownOffset = shifts.upDownOffset + 10;
            case('j') %want to move top image left
                shifts.leftRightOffset = shifts.leftRightOffset - 10;
            case('l') %want to move top image right
                shifts.leftRightOffset = shifts.leftRightOffset + 10;
            case('e')
                shifts.angle = shifts.angle + shifts.angleStep;
            case('r')
                shifts.angle = shifts.angle - shifts.angleStep;
            case('x')
                dropImage = true;
            case(28) %For right and left arrows,
            case(29) %will add this soon
        otherwise
        end
        app.currShift = shifts; % could also stick into handles;
        % app.lastKey = key;
        if ~dropImage
            app = UpdateOverlay(app);
        else
            app = UpdateTilePos(app);
        end
    end
end