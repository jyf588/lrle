function MusclesOI = getMuscleProperties(osimModel,MusclesOI)
nMuscles = osimModel.getMuscles.getSize;
state = osimModel.initSystem;
Muscle_nr = 1;
for i = 1:length(MusclesOI)
    muscle = osimModel.getMuscles.get(MusclesOI(i).index-1);
    disabled = muscle.isDisabled(state);
    if disabled
    else
        
        MusclesOI(Muscle_nr).maxIsoForce  = muscle.getMaxIsometricForce;
        MusclesOI(Muscle_nr).optFL  = muscle.getOptimalFiberLength;
        MusclesOI(Muscle_nr).tendonSL  = muscle.getTendonSlackLength;
        MusclesOI(Muscle_nr).alpha  = muscle.getPennationAngleAtOptimalFiberLength;
        MusclesOI(Muscle_nr).maxConVel  = muscle.getMaxContractionVelocity;
        Muscle_nr = Muscle_nr + 1;
    end  
end
end
