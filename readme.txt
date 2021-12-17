To launch a problematic instance with alternance between pricing and cutting use the following command : "./SCIP_Coloring_BP ../test/le450_15a.col ca"
                                                                                           or this one : "./SCIP_Coloring_BP ../test/4-FullIns_3.col ca"

The pricing class used is : "SCIP_pricer_alternate"

The separation class used is : "ChvatalG_M_A_sep"

note that separation and pricing are both solving a MILP.
