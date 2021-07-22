//---------------------------------Spheral++----------------------------------//
// SpheralCTimers -- Declare the Timers for use in profiling SpheralC
//----------------------------------------------------------------------------//
#include "Utilities/Timer.hh"

Timer TIME_SpheralC                                ("Root Timer for SpheralC");
Timer TIME_SpheralC_initialize                     ("SpheralPsuedoScript::initialize", TIME_SpheralC);
Timer TIME_SpheralC_updateState                    ("SpheralPsuedoScript::updateState", TIME_SpheralC);
Timer TIME_SpheralC_initializeBoundariesAndPhysics ("SpheralPsuedoScript::initializeBoundariesAndPhysics", TIME_SpheralC);
Timer TIME_SpheralC_initializeStep                 ("SpheralPsuedoScript::initializeStep", TIME_SpheralC);
Timer TIME_SpheralC_evaluateDerivatives            ("SpheralPsuedoScript::evaluateDerivatives", TIME_SpheralC);
Timer TIME_SpheralC_iterateHfield                  ("SpheralPsuedoScript::iterateHfield", TIME_SpheralC);
Timer TIME_SpheralC_computeFragmentID              ("SpheralPsuedoScript::computeFragmentID", TIME_SpheralC);
Timer TIME_SpheralC_sampleLatticeMesh              ("SpheralPsuedoScript::sampleLatticeMesh", TIME_SpheralC);
Timer TIME_SpheralC_polyhedralMesh                 ("SpheralPsuedoScript::polyhedralMesh", TIME_SpheralC);
