#include "common.h"
#include "sim.h"
#include "trace.h"
#include <stdlib.h>
/*******************************************************************/
/* Simulator frame */
/*******************************************************************/

bool run_a_cycle();
void init_structures();

//// MY GLOBAL Variables ////
//Different type of Bubbles //
int latency_stall = FALSE;
int control_hazard = FALSE;
int data_hazard = FALSE;
//Others //
int latency = 0;  // Counter for execution time for ops
int num_reg_inv = 0; // Change if tracking of particular register required.

//// Source Regs having data dependancy ////
int source[2];
//// End my  Global variables ////

/* uop_pool related variables */

uint32_t free_op_num;
uint32_t active_op_num;
Op *op_pool;
Op *op_pool_free_head = NULL;

/* simulator related local functions */

bool icache_access(uint32_t addr);
bool dcache_access(uint32_t addr);
void init_latches(void);
void init_registers(void);

#include "knob.h"
#include "all_knobs.h"

// knob variables
KnobsContainer *g_knobsContainer; /* < knob container > */
all_knobs_c    *g_knobs; /* < all knob variables > */

gzFile g_stream;

void init_knobs(int argc, char** argv)
{
  // Create the knob managing class
  g_knobsContainer = new KnobsContainer();

  // Get a reference to the actual knobs for this component instance
  g_knobs = g_knobsContainer->getAllKnobs();

  // apply the supplied command line switches
  char* pInvalidArgument = NULL;
  g_knobsContainer->applyComandLineArguments(argc, argv, &pInvalidArgument);

  g_knobs->display();
}

void read_trace_file(void)
{
  g_stream = gzopen((KNOB(KNOB_TRACE_FILE)->getValue()).c_str(), "r");
}

// simulator main function is called from outside of this file

void simulator_main(int argc, char** argv)
{
  init_knobs(argc, argv);

  // trace driven simulation
  read_trace_file();
  init_structures();
  run_a_cycle();

}
int op_latency[NUM_OP_TYPE];

void init_op_latency(void)
{
  op_latency[OP_INV]   = 1;
  op_latency[OP_NOP]   = 1;
  op_latency[OP_CF]    = 1;
  op_latency[OP_CMOV]  = 1;
  op_latency[OP_LDA]   = 1;
  op_latency[OP_LD]    = 1;
  op_latency[OP_ST]    = 1;
  op_latency[OP_IADD]  = 1;
  op_latency[OP_IMUL]  = 2;
  op_latency[OP_IDIV]  = 4;
  op_latency[OP_ICMP]  = 2;
  op_latency[OP_LOGIC] = 1;
  op_latency[OP_SHIFT] = 2;
  op_latency[OP_BYTE]  = 1;
  op_latency[OP_MM]    = 2;
  op_latency[OP_FMEM]  = 2;
  op_latency[OP_FCF]   = 1;
  op_latency[OP_FCVT]  = 4;
  op_latency[OP_FADD]  = 2;
  op_latency[OP_FMUL]  = 4;
  op_latency[OP_FDIV]  = 16;
  op_latency[OP_FCMP]  = 2;
  op_latency[OP_FBIT]  = 2;
  op_latency[OP_FCMP]  = 2;
  op_latency[OP_END]   = 1;   // My addition for end of simulation
}

void init_op(Op *op)
{
  op->num_src               = 0;
  op->src[0]                = -1;
  op->src[1]                = -1;
  op->dst                   = -1;
  op->opcode                = 0;
  op->is_fp                 = false;
  op->cf_type               = NOT_CF;
  op->mem_type              = NOT_MEM;
  op->write_flag             = 0;
  op->inst_size             = 0;
  op->ld_vaddr              = 0;
  op->st_vaddr              = 0;
  op->instruction_addr      = 0;
  op->branch_target         = 0;
  op->actually_taken        = 0;
  op->mem_read_size         = 0;
  op->mem_write_size        = 0;
  op->valid                 = FALSE;
  /* you might add more features here */
}


void sim_end_op(Op *op)			//// SD's addition for end of simulation
{
  op->num_src               = 0;
  op->src[0]                = -1;
  op->src[1]                = -1;
  op->dst                   = -1;
  op->opcode                = OP_END;
  op->is_fp                 = false;
  op->cf_type               = NOT_CF;
  op->mem_type              = NOT_MEM;
  op->write_flag             = 0;
  op->inst_size             = 0;
  op->ld_vaddr              = 0;
  op->st_vaddr              = 0;
  op->instruction_addr      = 0;
  op->branch_target         = 0;
  op->actually_taken        = 0;
  op->mem_read_size         = 0;
  op->mem_write_size        = 0;
  op->valid                 = FALSE;
  /* you might add more features here */
}


void init_op_pool(void)
{
  /* initialize op pool */
  op_pool = new Op [1024];
  free_op_num = 1024;
  active_op_num = 0;
  uint32_t op_pool_entries = 0;
  int ii;
  for (ii = 0; ii < 1023; ii++) {

    op_pool[ii].op_pool_next = &op_pool[ii+1];
    op_pool[ii].op_pool_id   = op_pool_entries++;
    init_op(&op_pool[ii]);
  }
  op_pool[ii].op_pool_next = op_pool_free_head;
  op_pool[ii].op_pool_id   = op_pool_entries++;
  init_op(&op_pool[ii]);
  op_pool_free_head = &op_pool[0];
}


Op *get_free_op(void)
{
  /* return a free op from op pool */

  if (op_pool_free_head == NULL || (free_op_num == 1)) {
    std::cout <<"ERROR! OP_POOL SIZE is too small!! " << endl;
    std::cout <<"please check free_op function " << endl;
    assert(1);
    exit(1);
  }

  free_op_num--;
  assert(free_op_num);

  Op *new_op = op_pool_free_head;
  op_pool_free_head = new_op->op_pool_next;
  assert(!new_op->valid);
  init_op(new_op);
  active_op_num++;
  return new_op;
}

void free_op(Op *op)
{
  free_op_num++;
  active_op_num--;
  op->valid = FALSE;
  op->op_pool_next = op_pool_free_head;
  op_pool_free_head = op;
}



/*******************************************************************/
/*  Data structure */
/*******************************************************************/

typedef struct pipeline_latch_struct {
  Op *op; /* you must update this data structure. */
  bool op_valid;
   /* you might add more data structures. But you should complete the above data elements */
}pipeline_latch;


typedef struct Reg_element_struct{
  bool valid;
  int same_reg_inv_cnt; // This is required because ... multiple statements in the pipeline may have the same destination. So, when the 1st of these group instructions reaches WB stage it sets the reg valid false, whereas the last instruction in the group if its in execute stage it still requires it to be invalid. Now the data hazard is missed.
  // data is not needed
  /* you might add more data structures. But you should complete the above data elements */
}REG_element;

REG_element register_file[NUM_REG];


/*******************************************************************/
/* These are the functions you'll have to write.  */
/*******************************************************************/

void FE_stage();
void ID_stage();
void EX_stage();
void MEM_stage();
void WB_stage();

/*******************************************************************/
/*  These are the variables you'll have to write.  */
/*******************************************************************/

bool sim_end_condition = FALSE;     /* please complete the condition. */
UINT64 retired_instruction = 0;    /* number of retired instruction. (only correct instructions) */
UINT64 cycle_count = 0;            /* total number of cycles */
UINT64 data_hazard_count = 0;
UINT64 control_hazard_count = 0;
UINT64 icache_miss_count = 0;      /* total number of icache misses. for Lab #2 and Lab #3 */
UINT64 dcache_miss_count = 0;      /* total number of dcache  misses. for Lab #2 and Lab #3 */
UINT64 l2_cache_miss_count = 0;    /* total number of L2 cache  misses. for Lab #2 and Lab #3 */

pipeline_latch *MEM_latch;
pipeline_latch *EX_latch;
pipeline_latch *ID_latch;
pipeline_latch *FE_latch;
UINT64 ld_st_buffer[LD_ST_BUFFER_SIZE];
UINT64 next_pc;

/*******************************************************************/
/*  Print messages  */
/*******************************************************************/

void print_stats() {
  std::ofstream out((KNOB(KNOB_OUTPUT_FILE)->getValue()).c_str());
  /* Do not modify this function. This messages will be used for grading */
  out << "Total instruction: " << retired_instruction << endl;
  out << "Total cycles: " << cycle_count << endl;
  float ipc = (cycle_count ? ((float)retired_instruction/(float)cycle_count): 0 );
  out << "IPC: " << ipc << endl;
  out << "Total I-cache miss: " << icache_miss_count << endl;
  out << "Total D-cache miss: " << dcache_miss_count << endl;
  out << "Total L2-cache miss: " << l2_cache_miss_count << endl;
  out << "Total data hazard: " << data_hazard_count << endl;
  out << "Total control hazard : " << control_hazard_count << endl;
  out.close();
}

/*******************************************************************/
/*  Support Functions  */
/*******************************************************************/

bool get_op(Op *op)
{
  static UINT64 unique_count = 0;
  Trace_op trace_op;
  bool success = FALSE;
  // read trace
  // fill out op info
  // return FALSE if the end of trace
  success = (gzread(g_stream, &trace_op, sizeof(Trace_op)) >0 );
  if (KNOB(KNOB_PRINT_INST)->getValue()) dprint_trace(&trace_op);

  /* copy trace structure to op */
  if (success) {
    copy_trace_op(&trace_op, op);

    op->inst_id  = unique_count++;
    op->valid    = TRUE;
  }
  return success;
}
/* return op execution cycle latency */

int get_op_latency (Op *op)
{
  assert (op->opcode < NUM_OP_TYPE);
  return op_latency[op->opcode];
}

/* Print out all the register values */
void dump_reg() {
  for (int ii = 0; ii < NUM_REG; ii++) {
    std::cout << cycle_count << ":register[" << ii  << "]: V:" << register_file[ii].valid << endl;
  }
}

void print_pipeline() {
  std::cout << "--------------------------------------------" << endl;
  std::cout <<"cycle count : " << dec << cycle_count << " retired_instruction : " << retired_instruction << endl;
  std::cout << (int)cycle_count << " FE: " ;
  if (FE_latch->op) {
    Op *op = FE_latch->op;
    cout << (int)op->inst_id ;
  }
  else {
    cout <<"####";
  }
  std::cout << " ID: " ;
  if (ID_latch->op_valid) {
    Op *op = ID_latch->op;
    cout << (int)op->inst_id ;
  }
  else {
    cout <<"####";
  }
  std::cout << " EX: " ;
  if (EX_latch->op_valid) {
    Op *op = EX_latch->op;
    cout << (int)op->inst_id ;
  }
  else {
    cout <<"####";
  }


  std::cout << " MEM: " ;
  if (MEM_latch->op_valid) {
    Op *op = MEM_latch->op;
    cout << (int)op->inst_id ;
  }
  else {
    cout <<"####";
  }
  cout << endl;
  //  dump_reg();
  std::cout << "--------------------------------------------" << endl;
}

void print_heartbeat()
{
  static uint64_t last_cycle ;
  static uint64_t last_inst_count;
  float temp_ipc = float(retired_instruction - last_inst_count) /(float)(cycle_count-last_cycle) ;
  float ipc = float(retired_instruction) /(float)(cycle_count) ;
  /* Do not modify this function. This messages will be used for grading */
  cout <<"**Heartbeat** cycle_count: " << cycle_count << " inst:" << retired_instruction << " IPC: " << temp_ipc << " Overall IPC: " << ipc << endl;
  last_cycle = cycle_count;
  last_inst_count = retired_instruction;
}
/*******************************************************************/
/*                                                                 */
/*******************************************************************/

bool run_a_cycle(){

  for (;;) {
    if (((KNOB(KNOB_MAX_SIM_COUNT)->getValue() && (cycle_count >= KNOB(KNOB_MAX_SIM_COUNT)->getValue())) ||
      (KNOB(KNOB_MAX_INST_COUNT)->getValue() && (retired_instruction >= KNOB(KNOB_MAX_INST_COUNT)->getValue())) ||  (sim_end_condition))) {
        // please complete sim_end_condition
        // finish the simulation
        print_heartbeat();
        print_stats();
        return TRUE;
    }
    cycle_count++;
    if (!(cycle_count%5000)) {
      print_heartbeat();
    }
    WB_stage();
    MEM_stage();
    EX_stage();
    ID_stage();
    FE_stage();
    if (KNOB(KNOB_PRINT_PIPE_FREQ)->getValue() && !(cycle_count%KNOB(KNOB_PRINT_PIPE_FREQ)->getValue())) print_pipeline();
  }
  return TRUE;
}


/*******************************************************************/
/* Complete the following fuctions.  */
/* You can add new data structures and also new elements to Op, Pipeline_latch data structure */
/*******************************************************************/

void init_structures(void)
{
	init_op_pool();
	init_op_latency();
	/* please initialize other data stucturs */
	/* you must complete the function */
	init_latches();
	init_registers();

}

void WB_stage()
{
//	cout << "WB stage" << endl;
	/* You MUST call free_op function here after *  an op is retired */
	/* you must complete the function */
	if(MEM_latch->op_valid == FALSE)
	{
		return;
	}
	
	else
	{
		Op *op = MEM_latch->op;
		
	//// END OF SIMULATION CHECK FOR OPCODE: OP_END ////
		if(op->opcode == OP_END)
		{
//			cout << "SIM END CONDITION" <<endl;
			sim_end_condition = 1;	
			return;
		}
		if((op->cf_type)>=CF_BR)   // if true its a branch. hold fetch, but continue passing ops
		{
			control_hazard = FALSE;
		}

		

	////Data hazard busy register to be free here ////
		if((int)(op->dst) != -1)
		{	
			register_file[(int)(op->dst)].same_reg_inv_cnt--;
			if(register_file[(int)(op->dst)].same_reg_inv_cnt == 0)
			register_file[(int)(op->dst)].valid = TRUE;
			
			for(int i=0; i < 2; i++)
			{
				
				if(source[i] == (int)(op->dst))
				{
				num_reg_inv--;
		//cout << "Number of Dependancies:" << num_reg_inv << endl;
		//cout << "Reg:" << (int)(op->dst) << " was one of the dependancies, now freed" << endl; 
				source[i] = 0;
				if(num_reg_inv == 0)
				{
					data_hazard = FALSE;
//			cout << "DATA HAZARD REMOVED" << endl;
				}

				}
			}
			//cout << "Register :" << (int) op->dst << "freed right now" << endl;
			
		}
	retired_instruction++;
	free_op(op);
	}
}

void MEM_stage()
{
//	cout << "MEM Stage " << endl;
	if(EX_latch->op_valid == FALSE)
	{
		MEM_latch->op_valid = FALSE;
		MEM_latch->op = NULL;
		return;
	}
	else
	{
	Op *op = EX_latch->op;
	MEM_latch->op = op;
	MEM_latch->op_valid = TRUE;
	}
  /* you must complete the function */
}


void EX_stage()
{
//	cout << "EX_stage" <<endl;
	/* you must complete the function */

	if(ID_latch->op_valid == FALSE)
	{
	  EX_latch->op_valid = FALSE;
	  EX_latch->op = NULL;
            return;
	}
    if(latency != 0 )
	{
		latency--;
//		cout<< "latency is "<<latency<<endl;;
		EX_latch->op_valid = FALSE;
		EX_latch->op = NULL;
		latency_stall = TRUE;
		if (latency != 0)
		return;
	}
//	cout << "latency = " << latency << "latency stall = " << latency_stall <<endl;
	if(latency == 0 && latency_stall == 1)
	{

	Op *op = ID_latch->op;
	EX_latch->op = op;
	EX_latch->op_valid = TRUE;
	latency_stall = FALSE;
//	cout << "NO more latency stall for 1 cycle. .. let the next instr flow in !" <<endl;
	}

	else
	{
	Op *op = ID_latch->op;
	latency = get_op_latency(op) - 1;
//	cout<< "issue latency is "<<latency<<endl;;
	if (latency != 0)
	{
	latency_stall = TRUE;
	EX_latch->op_valid = FALSE;
	EX_latch->op = NULL;
	}
	else
	{
	latency_stall = FALSE;
	EX_latch->op = op;
	EX_latch->op_valid = TRUE;
	}
	}
}

void ID_stage()
{
//	cout << " INVALID REG CHECK " << endl;
//	for(int i = 0; i < 32 ; i++)
//	{
//		if(register_file[i].valid == FALSE)
//		cout << "Register :" << i << "  is still invalid" << endl;
//	}
//	cout<< " ______________________ " <<endl;
//	cout << "ID Stage "   << endl;		
//	cout<< " ______________________ " <<endl;
//	cout << "Number of Dependancies:" << num_reg_inv << endl;
//	cout << "Source 1: "   << source[0] << endl;
//	cout << "source 2: "   << source[1] << endl;
   	if(latency_stall == TRUE)
	return;
	else if(FE_latch->op_valid == TRUE && data_hazard == FALSE)
	{
//		cout<<" Instruction Decode stage entered, Because data hazard was false " << endl;
		Op *op = FE_latch->op;
//// Data Hazard DBZ///
		for(int i = 0; i < op->num_src; i++)
		{
			if(register_file[(int)op->src[i]].valid == FALSE)
			{	
//				cout << "Register :" <<(int) op->src[i] << "found invalid" << endl;
				source[i] = (int) op->src[i];
//				cout << "Number of Dependancies:" << num_reg_inv << endl;
				num_reg_inv++;
				data_hazard = TRUE;
//				cout << "DATA Hazard Found" << endl;
//				cout << "Number of Dependancies:" << num_reg_inv << endl;
			}
		}
			
		if(data_hazard == TRUE)
		{
			ID_latch->op_valid = FALSE;
			ID_latch->op = NULL;
			data_hazard_count++;
			return;
		}	
			
			
    		if((int) op->dst != -1)
		{
//			cout << " HI, I am invalidating Destination Register: " << (int)op->dst << endl; 
			register_file[(int)op->dst].valid = FALSE;
			register_file[(int)op->dst].same_reg_inv_cnt++;
		}
////Cheking for Branch OP here. Following FAQ from website ////
		if((op->cf_type)>=CF_BR)   // if true its a branch. hold fetch, but continue passing ops
		{	
			control_hazard = TRUE;  // This will be reset in the memory function, when the branch is resolved.
			control_hazard_count++;
//			cout << "control hazard detected" << endl;
		}

		
			ID_latch->op_valid = TRUE;
			ID_latch->op = op;
		
	}
	else
	{
		if(data_hazard == TRUE)
		data_hazard_count++;
		ID_latch->op_valid = FALSE;
		ID_latch-> op = NULL;
        	return;
	}
  /* you must complete tihe function */
}


void FE_stage()
{
  /* only part of FE_stage function is implemented */
  /* please complete the rest of FE_stage function */
    
//cout << "DATA Hazard: " << data_hazard;
 //   cout << "FE Stage" <<endl;
	if(latency_stall == TRUE || data_hazard == TRUE)
		return;
	if(control_hazard == TRUE)
	{
		FE_latch->op_valid = FALSE;
		FE_latch->op = NULL;
	}
	else
	{
		Op *op = get_free_op();
		if( get_op(op) == TRUE)
       		 {
            		FE_latch->op = op;    ////Pass on the OP to the decode stage via latch
            		FE_latch->op_valid = TRUE;
        	}
		else {
			sim_end_op(op);
			FE_latch->op = op;
			FE_latch->op_valid = TRUE;
		}
	}
  //   next_pc = pc + op->inst_size;  // you need this code for building a branch predictor

}

void init_registers()
{
        for(int i = 0; i < NUM_REG; i++)
        {
            register_file[i].valid = TRUE;
	    register_file[i].same_reg_inv_cnt = 0;

        }
}

void  init_latches()
{
  MEM_latch = new pipeline_latch();
  EX_latch = new pipeline_latch();
  ID_latch = new pipeline_latch();
  FE_latch = new pipeline_latch();

  MEM_latch->op = NULL;
  EX_latch->op = NULL;
  ID_latch->op = NULL;
  FE_latch->op = NULL;

  /* you must set valid value correctly  */
  MEM_latch->op_valid = false;
  EX_latch->op_valid = false;
  ID_latch->op_valid = false;
  FE_latch->op_valid = false;

}

bool icache_access(uint32_t addr) {

  /* For Lab #1, you assume that all I-cache hit */
  bool hit = FALSE;
  if (KNOB(KNOB_PERFECT_ICACHE)->getValue()) hit = TRUE;
  return hit;
}



bool dcache_access(uint32_t addr) {
  /* For Lab #1, you assume that all D-cache hit */
  bool hit = FALSE;
  if (KNOB(KNOB_PERFECT_DCACHE)->getValue()) hit = TRUE;
  return hit;
}
