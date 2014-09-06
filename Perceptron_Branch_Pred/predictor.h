#ifndef _PREDICTOR_H_
#define _PREDICTOR_H_

#include "utils.h"
#include "tracer.h"


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

class PREDICTOR{

  // The state is defined for Gshare, change for your design

 private:
  UINT64  ghr;           // global history register
  UINT64  *pht;          // pattern history table
  UINT64  historyLength; // history length
  UINT64  numPhtEntries; // entries in pht 
  UINT64  numPercEntries; // entries in perceptron
  int     **perc_table;  // perceptron table
  int     global_total;
 public:

  // The interface to the four functions below CAN NOT be changed

  PREDICTOR(void);
  bool    GetPrediction(UINT32 PC);  
  void    UpdatePredictor(UINT32 PC, bool resolveDir, bool predDir, UINT32 branchTarget);
  void    TrackOtherInst(UINT32 PC, OpType opType, UINT32 branchTarget);

  // Contestants can define their own functions below

};



/***********************************************************/
#endif

