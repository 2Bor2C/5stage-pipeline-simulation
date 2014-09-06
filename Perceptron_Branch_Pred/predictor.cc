#include "predictor.h"
#include <climits>

#define PHT_CTR_MAX  3
#define PHT_CTR_INIT 2
#define HIST_LEN   127
//#define original  
/////////////// STORAGE BUDGET JUSTIFICATION ////////////////
// Total storage budget: 32KB + 17 bits
// Total PHT counters: 2^17 
// Total PHT size = 2^17 * 2 bits/counter = 2^18 bits = 32KB
// GHR size: 17 bits
// Total Size = PHT size + GHR size
/////////////////////////////////////////////////////////////



/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

#ifdef original 
PREDICTOR::PREDICTOR(void){

  historyLength    = HIST_LEN;
  ghr              = 0;
  numPhtEntries    = (1 << HIST_LEN);

  pht = new UINT64[numPhtEntries];

  for(UINT32 ii=0; ii< numPhtEntries; ii++){
    pht[ii]=PHT_CTR_INIT; 
  }
  
}
#else
PREDICTOR::PREDICTOR(void){

  historyLength    = HIST_LEN;
  ghr              = 0;
  numPercEntries   = 256;

  perc_table = new int*[numPercEntries];
  
  for(UINT32 ii=0; ii< numPercEntries; ii++)
  {
	 perc_table[ii] =  new int[HIST_LEN+1];
	 
  }
 
  for(UINT32 j=0; j < numPercEntries; j++)
  {
 	for(UINT32 ii=0; ii < HIST_LEN + 1; ii++)
 	{
 	   perc_table[j][ii] = 0;	
	} 
  }  
}
#endif
/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
#ifdef original
bool   PREDICTOR::GetPrediction(UINT32 PC){

  UINT32 phtIndex   = (PC^ghr) % (numPhtEntries);
  UINT32 phtCounter = pht[phtIndex];
  
  if(phtCounter > PHT_CTR_MAX/2){
    return TAKEN; 
  }else{
    return NOT_TAKEN; 
  }
  
}
#else
bool   PREDICTOR::GetPrediction(UINT32 PC){

  UINT32 percIndex   = (PC^ghr) % (numPercEntries);
  //UINT32 phtCounter = pht[phtIndex];
  int local_total = 0;
  int local_ghr = ghr;
  for(UINT32 ii=0; ii < HIST_LEN + 1; ii++)
  {
	if(ii != HIST_LEN)
	{
		if((local_ghr >> ii) & 1)
		{
			local_total += perc_table[percIndex][ii]; 
		}
		else
		{
			local_total -= perc_table[percIndex][ii];
		}
	}
	else
	{
		local_total += perc_table[percIndex][ii];   ///bias weight
	}
  }
  global_total = local_total;
  if(local_total >= 0)
	return 1;
  else
	return 0;
}

#endif

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
#ifdef original
void  PREDICTOR::UpdatePredictor(UINT32 PC, bool resolveDir, bool predDir, UINT32 branchTarget){

  UINT32 phtIndex   = (PC^ghr) % (numPhtEntries);
  UINT32 phtCounter = pht[phtIndex];

  // update the PHT

  if(resolveDir == TAKEN){
    pht[phtIndex] = SatIncrement(phtCounter, PHT_CTR_MAX);
  }else{
    pht[phtIndex] = SatDecrement(phtCounter);
  }

  // update the GHR
  ghr = (ghr << 1);

  if(resolveDir == TAKEN){
    ghr++; 
  }

}
#else
void  PREDICTOR::UpdatePredictor(UINT32 PC, bool resolveDir, bool predDir, UINT32 branchTarget){

  UINT32 percIndex   = (PC^ghr) % (numPercEntries);
  //UINT32 phtCounter = pht[phtIndex];
  float threshold_f = (1.93 * HIST_LEN) + 14;
  int threshold_i = ((int)(threshold_f - (((int)threshold_f) % 1))) + 1;
  int local_ghr = ghr;
  // update the PHT
  // Modulus of total required 
  int local_total = global_total;
  if(local_total < 0)
	{
		local_total = (0 - local_total);
	}
  if((local_total < threshold_i) || (predDir != resolveDir)){
  
		for( int ii = 0; ii < HIST_LEN + 1; ii++)
		{
			// Check if w are not saturated +127 to -128 or INT_MAX INT_MIN ?? (@piazza)
			
				if(((local_ghr >> ii) & 1) == resolveDir)
				{
					if(perc_table[percIndex][ii] < 127)
						perc_table[percIndex][ii] += 1; 
				}
				else
				{
					if(perc_table[percIndex][ii] > -128)
						perc_table[percIndex][ii] -= 1;
				}
			
		}
	}
  // update the GHR
  ghr = (ghr << 1);

  if(resolveDir == TAKEN){
    ghr++; 
  }

}
#endif

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

void    PREDICTOR::TrackOtherInst(UINT32 PC, OpType opType, UINT32 branchTarget){

  // This function is called for instructions which are not
  // conditional branches, just in case someone decides to design
  // a predictor that uses information from such instructions.
  // We expect most contestants to leave this function untouched.

  return;
}

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
