function [MusclesOI,CoordinatesOI] = getMuscleStructure( Model_OS,CoordinatesOI_Names )

% MusclesOI:
    % Structure that holds information of the muscles that actuate the
    % degrees of freedom (or coordinates) Of Interest
    % ---> Name,Index,CoordinatesActuated(Name,Index),PropertiesOfMuscle(related to the Hill model)
% MusclesOI:
    % Name, Index of the coordinates Of Interest

% NOTE: "Index" refers to the index in the OpenSim model 



% Get all Coordinate Names in the model
Nr_Coordinates = Model_OS.getCoordinateSet.getSize;
for i = 1:Nr_Coordinates
    CoordinateNames(i).name = Model_OS.getCoordinateSet.get(i-1).getName.toCharArray';
end

% Assign index to our coordinates of interest
for i = 1:length(CoordinatesOI_Names)
    for j = 1:length(CoordinateNames)
        if CoordinateNames(j).name == string(CoordinatesOI_Names(i))
            CoordinatesOI(i).index = j;
            break
        end
    end
end

% We get the information on the muscles of Interest
state = Model_OS.initSystem;
Model_OS.equilibrateMuscles(state); 

Nr_Muscles = Model_OS.getMuscles.getSize;
incM = 0;
MuscleNr = 1;
for m = 1:Nr_Muscles
    muscle = Model_OS.getMuscles.get(m-1);
    Coordinate_Nr = 1;
    for i = 1:length(CoordinatesOI)
        coordinateIndex = CoordinatesOI(i).index;
        coordinate = Model_OS.getCoordinateSet.get(coordinateIndex-1);
        if abs(muscle.computeMomentArm(state, coordinate)) > 0.0001
            MusclesOI(MuscleNr).name = muscle.getName;
            MusclesOI(MuscleNr).index = m;
            MusclesOI(MuscleNr).coordinate(Coordinate_Nr).coordinateName = coordinate.getName.toCharArray';
            MusclesOI(MuscleNr).coordinate(Coordinate_Nr).coordinateIndex = coordinateIndex;
            incM = 1;
            Coordinate_Nr = Coordinate_Nr + 1;
        end
    end
    if incM == 1
        MuscleNr = MuscleNr + 1;
    end
    inc = 0;
end

% Put the Muscle Properties into structure
% These are the properties related to the hill model (isometric force,
% tendon slack length, etc...)
MusclesOI = getMuscleProperties(Model_OS,MusclesOI);



end

