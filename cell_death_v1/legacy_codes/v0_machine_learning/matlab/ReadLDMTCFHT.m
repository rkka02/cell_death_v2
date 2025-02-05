function data = ReadLDMTCFHT(tcfFilePath, timeFrameIndex)
    dataPath = sprintf('/Data/3D/%06d', timeFrameIndex);

    dataSizeX = h5readatt(tcfFilePath,dataPath, 'DataSizeX');
    dataSizeY = h5readatt(tcfFilePath,dataPath, 'DataSizeY');
    dataSizeZ = h5readatt(tcfFilePath,dataPath, 'DataSizeZ');

    tileCount = h5readatt(tcfFilePath,dataPath, 'NumberOfTiles');

    minValue = h5readatt(tcfFilePath, dataPath, 'RIMin');

    data = zeros(dataSizeY, dataSizeX, dataSizeZ, 'single');
    
    scalarType = h5readatt(tcfFilePath, dataPath, 'ScalarType');

    for tileIndex = 1 : tileCount
        tileString = sprintf('%02d',tileIndex-1);
        tilePath = strcat(dataPath, '/TILE_', tileString);
        
        samplingStep = h5readatt(tcfFilePath, tilePath, 'SamplingStep');
        
        if (samplingStep ~= 1)
            continue;
        end

        if (scalarType == 0)
            tileData = single(permute(h5read(tcfFilePath, tilePath), [2,1,3]))/10000;
        else
            tileData = single(permute(h5read(tcfFilePath, tilePath), [2,1,3]))/1000 + minValue;
        end
        
        offsetX = h5readatt(tcfFilePath, tilePath, 'DataIndexOffsetPointX') + 1;
        offsetY = h5readatt(tcfFilePath, tilePath, 'DataIndexOffsetPointY') + 1;
        offsetZ = h5readatt(tcfFilePath, tilePath, 'DataIndexOffsetPointZ') + 1;
        
        lastX = h5readatt(tcfFilePath, tilePath, 'DataIndexLastPointX') + 1;
        lastY = h5readatt(tcfFilePath, tilePath, 'DataIndexLastPointY') + 1;
        lastZ = h5readatt(tcfFilePath, tilePath, 'DataIndexLastPointZ') + 1;

        mappingRangeX = offsetX:lastX;
        mappingRangeY = offsetY:lastY;
        mappingRangeZ = offsetZ:lastZ;
        
        dataRangeX = mappingRangeX - offsetX + 1;
        dataRangeY = mappingRangeY - offsetY + 1;
        dataRangeZ = mappingRangeZ - offsetZ + 1;
        
        data(mappingRangeY, mappingRangeX, mappingRangeZ) = tileData(dataRangeY, dataRangeX, dataRangeZ);
    end
end

