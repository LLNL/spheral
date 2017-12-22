//---------------------------------Spheral++----------------------------------//
// SpheralTimers -- Declare the Timers for use in profiling Spheral
//
// Created by J. Michael Owen, Sat Dec  1 15:03:19 PST 2001
//----------------------------------------------------------------------------//

#include "Timer.hh"

// Must initialize the static list defined in Timer.hh
#ifdef TIMER
list<Timer*> Timer::TimerList(0); 
#endif

//------------------------------------------------------------------------------
// Root Timers
//------------------------------------------------------------------------------
Timer TimeSpheral               ("Root Timer                ");
// Timer TimeNestedGridNeighbor    ("Root NestedGridNeighbor   ");
// Timer TimePhysics               ("All physics derivatives   ");
// Timer TimeSphNodeList           ("Root Sph NodeList         ");
// Timer TimeMashNodeList          ("Root MASH NodeList        ");
// Timer TimeDataBase              ("Root DataBase             ");
// Timer TimeDistributedBound      ("Root DistributedBoundary  ");
// Timer TimeNestedDistributedBound("Root NestedDistribBound   ");

//------------------------------------------------------------------------------
// clipFacetedVolume timers
//------------------------------------------------------------------------------
Timer TIME_clipFacetedVolumeByPlanes3d("clipFacetedVolumeByPlanes3d");
Timer TIME_convertfrom("Convert from Spheral::Polyhedron", TIME_clipFacetedVolumeByPlanes3d);
Timer TIME_clipverts("Clip vertices", TIME_clipFacetedVolumeByPlanes3d);
Timer TIME_clipedges("Clip edges", TIME_clipFacetedVolumeByPlanes3d);
Timer TIME_clipfaces("Clip faces", TIME_clipFacetedVolumeByPlanes3d);
Timer   TIME_clipEdgesInFace("Clip edges in face", TIME_clipfaces);
Timer   TIME_eraseface("Erase face", TIME_clipfaces);
Timer   TIME_erasedups("Erase node duplicates", TIME_clipfaces);
Timer   TIME_sortline("Sort new face nodes in line", TIME_clipfaces);
Timer   TIME_newedges("Construct new edges in face", TIME_clipfaces);
Timer   TIME_newloops("Construct new face loops", TIME_clipfaces);
Timer   TIME_deactivate("Mark clipped vertices inactive", TIME_clipfaces);
Timer TIME_cap("Cap new faces", TIME_clipFacetedVolumeByPlanes3d);
Timer TIME_convertto("Convert back to Spheral::Polyhedron", TIME_clipFacetedVolumeByPlanes3d);

// //------------------------------------------------------------------------------
// // Second order predictor corrector integrator
// //------------------------------------------------------------------------------
// Timer TimePC2               ("Root Pred-Correct 2nd or", TimeSpheral);
// Timer TimePC2Boundaries1    ("Set boundaries 1        ", TimePC2);
// Timer TimePC2Dt             ("Set dt                  ", TimePC2);
// Timer TimePC2StepInitialize1("Pre-step initialize 1   ", TimePC2);
// Timer TimePC2ReadDB         ("Read DataBase FieldLists", TimePC2);
// Timer TimePC2CopyFields     ("Copy initial Derivs     ", TimePC2);
// Timer TimePC2PredStep       ("Predictor step          ", TimePC2);
// Timer TimePC2Weights1       ("Update Weights 1        ", TimePC2);
// Timer TimePC2StepInitialize2("Pre-step initialize 2   ", TimePC2);
// Timer TimePC2Boundaries2    ("Set boundaries 2        ", TimePC2);
// Timer TimePC2Zero1          ("Zero derivs  1          ", TimePC2);
// Timer TimePC2CopyDvDx1      ("Copy DvDx to temp       ", TimePC2);
// Timer TimePC2EvalDerivs1    ("Eval derivs mid step    ", TimePC2);
// Timer TimePC2CorStep        ("Corrector step          ", TimePC2);
// Timer TimePC2CopyDvDx2      ("Copy DvDx to DB         ", TimePC2);
// Timer TimePC2Weights2       ("Update Weights 2        ", TimePC2);
// Timer TimePC2Boundaries3    ("Set boundaries 3        ", TimePC2);
// Timer TimePC2StepInitialize3("Pre-step initialize 3   ", TimePC2);
// Timer TimePC2Zero2          ("Zero derivs  2          ", TimePC2);
// Timer TimePC2EvalDerivs2    ("Eval derivs end step    ", TimePC2);
// Timer TimePC2CopyDvDx3      ("Copy DvDx to DB again   ", TimePC2);
// Timer TimePC2SumDensity     ("Sum end mass density    ", TimePC2);
// Timer TimePC2Weights3       ("Update Weights 3        ", TimePC2);

// //------------------------------------------------------------------------------
// // RK2 Synchronous integrator
// //------------------------------------------------------------------------------
// Timer TimeSyncRK2               ("Root Sync RK2Integrator ", TimeSpheral);
// Timer TimeSyncRK2Boundaries1    ("Set boundaries 1        ", TimeSyncRK2);
// Timer TimeSyncRK2Dt             ("Set dt                  ", TimeSyncRK2);
// Timer TimeSyncRK2StepInitialize1("Pre-step initialize 1   ", TimeSyncRK2);
// Timer TimeSyncRK2ReadDB         ("Read DataBase FieldLists", TimeSyncRK2);
// Timer TimeSyncRK2CopyFields     ("Copy initial Derivs     ", TimeSyncRK2);
// Timer TimeSyncRK2MidStep        ("Advance to mid-step     ", TimeSyncRK2);
// Timer TimeSyncRK2Weights1       ("Update Weights 1        ", TimeSyncRK2);
// Timer TimeSyncRK2StepInitialize2("Pre-step initialize 2   ", TimeSyncRK2);
// Timer TimeSyncRK2Boundaries2    ("Set boundaries 2        ", TimeSyncRK2);
// Timer TimeSyncRK2Zero1          ("Zero derivs  1          ", TimeSyncRK2);
// Timer TimeSyncRK2CopyDvDx1      ("Copy DvDx to temp       ", TimeSyncRK2);
// Timer TimeSyncRK2EvalDerivs1    ("Eval derivs mid step    ", TimeSyncRK2);
// Timer TimeSyncRK2EndStep        ("Advance to end of step  ", TimeSyncRK2);
// Timer TimeSyncRK2CopyDvDx2      ("Copy DvDx to DB         ", TimeSyncRK2);
// Timer TimeSyncRK2Weights2       ("Update Weights 2        ", TimeSyncRK2);
// Timer TimeSyncRK2Boundaries3    ("Set boundaries 3        ", TimeSyncRK2);
// Timer TimeSyncRK2StepInitialize3("Pre-step initialize 3   ", TimeSyncRK2);
// Timer TimeSyncRK2Zero2          ("Zero derivs  2          ", TimeSyncRK2);
// Timer TimeSyncRK2EvalDerivs2    ("Eval derivs end step    ", TimeSyncRK2);
// Timer TimeSyncRK2CopyDvDx3      ("Copy DvDx to DB again   ", TimeSyncRK2);
// Timer TimeSyncRK2SumDensity     ("Sum end mass density    ", TimeSyncRK2);
// Timer TimeSyncRK2Weights3       ("Update Weights 3        ", TimeSyncRK2);

// //------------------------------------------------------------------------------
// // Cheap RK2 Synchronous integrator
// //------------------------------------------------------------------------------
// Timer TimeCheapRK2               ("Root Cheap RK2Integrator", TimeSpheral);
// Timer TimeCheapRK2Boundaries1    ("Set boundaries 1        ", TimeCheapRK2);
// Timer TimeCheapRK2Dt             ("Set dt                  ", TimeCheapRK2);
// Timer TimeCheapRK2StepInitialize1("Pre-step initialize 1   ", TimeCheapRK2);
// Timer TimeCheapRK2ReadDB         ("Read DataBase FieldLists", TimeCheapRK2);
// Timer TimeCheapRK2CopyFields     ("Copy initial Derivs     ", TimeCheapRK2);
// Timer TimeCheapRK2MidStep        ("Advance to mid-step     ", TimeCheapRK2);
// Timer TimeCheapRK2Weights1       ("Update Weights 1        ", TimeCheapRK2);
// Timer TimeCheapRK2Boundaries2    ("Set boundaries 2        ", TimeCheapRK2);
// Timer TimeCheapRK2StepInitialize2("Pre-step initialize 2   ", TimeCheapRK2);
// Timer TimeCheapRK2Zero1          ("Zero derivs  1          ", TimeCheapRK2);
// Timer TimeCheapRK2CopyDvDx1      ("Copy DvDx to temp       ", TimeCheapRK2);
// Timer TimeCheapRK2EvalDerivs1    ("Eval derivs mid step    ", TimeCheapRK2);
// Timer TimeCheapRK2EndStep        ("Advance to end of step  ", TimeCheapRK2);
// Timer TimeCheapRK2CopyDvDx2      ("Copy DvDx to DB         ", TimeCheapRK2);
// Timer TimeCheapRK2SumDensity     ("Sum end mass density    ", TimeCheapRK2);
// Timer TimeCheapRK2Weights2       ("Update Weights 2        ", TimeCheapRK2);
// Timer TimeCheapRK2Boundaries3    ("Set boundaries 3        ", TimeCheapRK2);

// //------------------------------------------------------------------------------
// // NestedGridNeighbor
// //------------------------------------------------------------------------------
// Timer TimeNestedMaster          ("Set the master list (node)", TimeNestedGridNeighbor);
// Timer TimeNestedRefine          ("Set the refine list (node)", TimeNestedGridNeighbor);
// Timer TimeNestedMasterPlane     ("Set the master list(plane)", TimeNestedGridNeighbor);
// Timer TimeNestedUpdateNodes     ("Update the node info      ", TimeNestedGridNeighbor);
// Timer TimeNestedUpdateNodes1    ("Update the given node info", TimeNestedGridNeighbor);
// Timer TimeNestedPrecull         ("Precull nodes             ", TimeNestedRefine);

// //------------------------------------------------------------------------------
// // Hydro
// //------------------------------------------------------------------------------
// Timer TimeHydro                 ("Root Hydro                ", TimePhysics);
// Timer TimeHydroFlagNodes        ("Initialize flag nodes     ", TimeHydro);
// Timer TimeHydroReadState        ("Read the initial state    ", TimeHydro);
// Timer TimeHydroSetMaster        ("Set master nodes          ", TimeHydro);
// Timer TimeHydroSetRefine        ("Set the refine nodes      ", TimeHydro);
// Timer TimeHydroEvalDerivs       ("Set derivs for node       ", TimeHydro);

// //------------------------------------------------------------------------------
// // Sph NodeList
// //------------------------------------------------------------------------------
// Timer TimeSphSum                ("Sph Mass Summation        ", TimeSphNodeList);

// Timer TimeSphDerivs             ("Base Sph Derivs           ", TimeSphNodeList);
// Timer TimeSphNodeIState         ("Get state for node i      ", TimeSphDerivs);
// Timer TimeSphNodeJState         ("Get state for node j      ", TimeSphDerivs);
// Timer TimeSphRij                ("Calc rij, etaij, ...      ", TimeSphDerivs);
// Timer TimeSphKernel             ("Evaluate Wij              ", TimeSphDerivs);
// Timer TimeSphQ                  ("Evaluate Qij              ", TimeSphDerivs);
// Timer TimeSphIncDerivs          ("Increment the derivs      ", TimeSphDerivs);
// Timer TimeSphHDeriv             ("Evaluate the H deriv      ", TimeSphDerivs);
// Timer TimeSphHSmooth            ("Smoothing contrib to H    ", TimeSphDerivs);

// //------------------------------------------------------------------------------
// // MASH NodeList
// //------------------------------------------------------------------------------
// Timer TimeMashDerivs            ("Base Mash Derivs          ", TimeMashNodeList);
// Timer TimeMashNodeIState        ("Get state for node i      ", TimeMashDerivs);
// Timer TimeMashNodeJState        ("Get state for node j      ", TimeMashDerivs);
// Timer TimeMashRij               ("Calc rij, etaij, ...      ", TimeMashDerivs);
// Timer TimeMashKernel            ("Evaluate Wij              ", TimeMashDerivs);
// Timer TimeMashQ                 ("Evaluate Qij              ", TimeMashDerivs);
// Timer TimeMashIncDerivs         ("Increment the derivs      ", TimeMashDerivs);
// Timer TimeMashHDeriv            ("Evaluate the H deriv      ", TimeMashDerivs);
// Timer TimeMashHSmooth           ("Smoothing contrib to H    ", TimeMashDerivs);

// //------------------------------------------------------------------------------
// // DataBase
// //------------------------------------------------------------------------------
// Timer TimeDataBaseRho           ("Sum fluid density         ", TimeDataBase);
// Timer TimeDataBaseGetState      ("Read state FieldLists     ", TimeDataBaseRho);
// Timer TimeDataBaseInitFlags     ("Init node done flags      ", TimeDataBaseRho);
// Timer TimeDataBaseSetMaster     ("Set master nodes          ", TimeDataBaseRho);
// Timer TimeDataBaseSetRefine     ("Set refine nodes          ", TimeDataBaseRho);
// Timer TimeDataBaseSumNode       ("Sum rho for single node   ", TimeDataBaseRho);
// Timer TimeDataBasePostConditions("Verify post conditions    ", TimeDataBaseRho);

// //------------------------------------------------------------------------------
// // DistributedBoundary
// //------------------------------------------------------------------------------
// Timer TimeDBExchange            ("Generic Field exchange    ", TimeDistributedBound);

// //------------------------------------------------------------------------------
// // NestedGridDistributedBoundary
// //------------------------------------------------------------------------------
// Timer TimeNDSetGhost            ("Set ghost nodes           ", TimeNestedDistributedBound);
// Timer TimeNDReset               ("Clear existing info       ", TimeNDSetGhost);
// Timer TimeNDFlatten             ("Flatten occupied gridcells", TimeNDSetGhost);
// Timer TimeNDReduceGridCells     ("Reduce occupied gridcells ", TimeNDSetGhost);
// Timer TimeNDBuildCommMaps       ("Build comm maps on master ", TimeNDSetGhost);
// Timer TimeNDDistributeCommMaps  ("Distribute comm maps      ", TimeNDSetGhost);
// Timer TimeNDBuildGhost          ("Build ghost nodes         ", TimeNDSetGhost);
// Timer TimeNDExchangeMinimal     ("Exchange r, H             ", TimeNDSetGhost);

// //------------------------------------------------------------------------------
// // Golden rule grad div vector field.
// //------------------------------------------------------------------------------
// Timer TimeGoldenSecond          ("Golden rule grad div      ");
// Timer TimeGoldenSecondReserve   ("Reserve return Fields     ", TimeGoldenSecond);
// Timer TimeGoldenSecondMaster    ("Select master nodes       ", TimeGoldenSecond);
// Timer TimeGoldenSecondRefine    ("Select refine nodes       ", TimeGoldenSecond);
// Timer TimeGoldenSecondNodeI     ("State for node I          ", TimeGoldenSecond);
// Timer TimeGoldenSecondNodeJ     ("State for node J          ", TimeGoldenSecond);
// Timer TimeGoldenSecondRefineLoop("Main loop over refine nei ", TimeGoldenSecond);
// Timer TimeGoldenSecondW         ("Calculate W, gradW, grad2W", TimeGoldenSecond);
// Timer TimeGoldenSecondAccumulate("Accumulate components     ", TimeGoldenSecond);
// Timer TimeGoldenSecondFinalize  ("Determine final answer    ", TimeGoldenSecond);
// Timer TimeGoldenSecondCheck     ("Completeness check        ", TimeGoldenSecond);
