#ifndef _SCIP_Alternate_Pricer
#define _SCIP_Alternate_Pricer

#include "objscip/objscip.h"
#include "scip/pub_var.h"

#include <vector>
#include <list>
#include <chrono>

#include "Graph.h"
#include "SCIP_master.h"
#include "Find_Stable_set_with_Cplex.h"
#include "ChvatalG_M_A_sep.h"
//#include "constraints.h"


using namespace std;
using namespace scip;


/** pricer class */
class PricerAlternate : public ObjPricer{
public:

   C_Graph *G;

   C_master_coloring *M;

   Cplex_pricing_algo A_cplex;

   int nbNewConstaints=0;
   int begin;
   ChvatalGMAsep* separator;
   int nbPricingRound;
   bool costIsNegative;

   /** Constructs the pricer object with the data needed */
   PricerAlternate(
      SCIP*                               scip,        /**< SCIP pointer */
      const char*                         p_name,      /**< name of pricer */
      C_Graph*                      G_,
      C_master_coloring *           M_
   
      );


   /** Destructs the pricer object. */
   virtual ~PricerAlternate();

   /** initialization method of variable pricer (called after problem was transformed) */
   virtual SCIP_DECL_PRICERINIT(scip_init);

   /** reduced cost pricing method of variable pricer for feasible LPs */
   virtual SCIP_DECL_PRICERREDCOST(scip_redcost);

   /** perform pricing */
   void coloring_pricing(
      SCIP*              scip               /**< SCIP data structure */
      );

   bool isSolInteger(SCIP* scip);
	SCIP_Real sumValue(SCIP* scip);


   
};

#endif
