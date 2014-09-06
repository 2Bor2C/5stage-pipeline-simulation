#include "common.h"
#include "sim.h"
#include "trace.h" 
#include "cache.h"  /**** NEW-LAB2*/ 
#include "memory.h" // NEW-LAB2 
#include <stdlib.h>
#include <ctype.h> /* Library for useful character operations */
/*******************************************************************/
/* Simulator frame */ 
/*******************************************************************/

bool run_a_cycle(memory_c *m ); // please modify run_a_cycle function argument  /** NEW-LAB2 */ 
void init_structures(memory_c *m); // please modify init_structures function argument  /** NEW-LAB2 */ 

list<Op*> temp_ops;
int temp_ops_list_empty = 0;
int dcache_latency = 0;
int dcache_latency_stall = FALSE;
int mshr_full_stall =  FALSE;
/* uop_pool related variables */ 
//New variables added by apkarande
UINT64 reg_writing_ops[NUM_REG];
UINT64 last_inst_id = 0;
int simulation_end_flag;
uint32_t free_op_num;
uint32_t active_op_num; 
Op *op_pool; 
Op *op_pool_free_head = NULL; 

void init_register_file();
/* simulator related local functions */ 

bool icache_access(ADDRINT addr); /** please change uint32_t to ADDRINT NEW-LAB2 */  
bool dcache_access(ADDRINT addr); /** please change uint32_t to ADDRINT NEW-LAB2 */   
void  init_latches(void);

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

  /** NEW-LAB2 */ /* just note: passing main memory pointers is hack to mix up c++ objects and c-style code */  /* Not recommended at all */ 
  memory_c *main_memory = new memory_c();  // /** NEW-LAB2 */ 


  init_structures(main_memory);  // please modify run_a_cycle function argument  /** NEW-LAB2 */ 
  run_a_cycle(main_memory); // please modify run_a_cycle function argument  /** NEW-LAB2 */ 


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
  bool stage_stall;
  list<Op*> queue_ops;
   /* you might add more data structures. But you should complete the above data elements */ 
}pipeline_latch; 


typedef struct Reg_element_struct{
  bool valid;
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
void MEM_stage(memory_c *main_memory); // please modify MEM_stage function argument  /** NEW-LAB2 */ 
void WB_stage(memory_c *main_memory); 

/*******************************************************************/
/*  These are the variables you'll have to write.  */ 
/*******************************************************************/

bool sim_end_condition = FALSE;     /* please complete the condition. */ 
UINT64 retired_instruction = 0;    /* number of retired instruction. (only correct instructions) */ 
UINT64 cycle_count = 0;            /* total number of cycles */ 
UINT64 data_hazard_count = 0;  
UINT64 control_hazard_count = 0; 
UINT64 icache_miss_count = 0;      /* total number of icache misses. for Lab #2 and Lab #3 */ 
UINT64 dcache_hit_count = 0;      /* total number of dcache  misses. for Lab #2 and Lab #3 */ 
UINT64 dcache_miss_count = 0;      /* total number of dcache  misses. for Lab #2 and Lab #3 */ 
UINT64 l2_cache_miss_count = 0;    /* total number of L2 cache  misses. for Lab #2 and Lab #3 */  
UINT64 dram_row_buffer_hit_count = 0; /* total number of dram row buffer hit. for Lab #2 and Lab #3 */   // NEW-LAB2
UINT64 dram_row_buffer_miss_count = 0; /* total number of dram row buffer hit. for Lab #2 and Lab #3 */   // NEW-LAB2
UINT64 store_load_forwarding_count = 0;  /* total number of store load forwarding for Lab #2 and Lab #3 */  // NEW-LAB2


pipeline_latch *MEM_latch;  
pipeline_latch *EX_latch;
pipeline_latch *ID_latch;
pipeline_latch *FE_latch;
UINT64 ld_st_buffer[LD_ST_BUFFER_SIZE];   /* this structure is deprecated. do not use */ 
UINT64 next_pc; 

Cache *data_cache;  // NEW-LAB2 

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
  out << "Total D-cache hit: " << dcache_hit_count << endl; 
  out << "Total D-cache miss: " << dcache_miss_count << endl; 
  out << "Total L2-cache miss: " << l2_cache_miss_count << endl; 
  out << "Total data hazard: " << data_hazard_count << endl;
  out << "Total control hazard : " << control_hazard_count << endl; 
  out << "Total DRAM ROW BUFFER Hit: " << dram_row_buffer_hit_count << endl;  // NEW-LAB2
  out << "Total DRAM ROW BUFFER Miss: "<< dram_row_buffer_miss_count << endl;  // NEW-LAB2 
  out <<" Total Store-load forwarding: " << store_load_forwarding_count << endl;  // NEW-LAB2 

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
  bool new_success = success && (trace_op.opcode < NUM_OP_TYPE) && (trace_op.opcode > 0);
  /* copy trace structure to op */ 
  if (new_success) { 
    copy_trace_op(&trace_op, op); 

    op->inst_id  = unique_count++;
    op->valid    = TRUE; 
    last_inst_id = op->inst_id + 1; 
 }
  else
  {
	last_inst_id = unique_count - 1;
  }
  return new_success; 
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
  if (FE_latch->op_valid) {
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
    //Op *op = MEM_latch->op; 
    //cout << (int)op->inst_id ;
    list<Op*>::iterator cii;
  for(cii = MEM_latch->queue_ops.begin(); cii != MEM_latch->queue_ops.end(); cii++)
  {
 	Op *op = (*cii);
	cout << (int)op->inst_id;
  }
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

bool run_a_cycle(memory_c *main_memory){   // please modify run_a_cycle function argument  /** NEW-LAB2 */ 

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

    main_memory->run_a_cycle();          // *NEW-LAB2 
    
    WB_stage(main_memory); 
    MEM_stage(main_memory);  // please modify MEM_stage function argument  /** NEW-LAB2 */ 
    mem_latch_queue_ops();
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

void init_structures(memory_c *main_memory) // please modify init_structures function argument  /** NEW-LAB2 */ 
{
  init_op_pool(); 
  init_op_latency();
  /* please initialize other data stucturs */ 
  /* you must complete the function */
  init_latches();
  	init_register_file();
	//Memory Initialisation
	main_memory->init_mem();
	
	// Cache initialisation
	data_cache = new Cache;
	cache_size = KNOB(KNOB_DCACHE_SIZE)->getValue();	
	block_size = KNOB(KNOB_BLOCK_SIZE)->getValue();
	assoc = KNOB(KNOB_DCACHE_WAY)->getValue();
	cache_init(data_cache,cache_size,block_size,assoc,"Data Cache");
}

void WB_stage(memory_c *main_memory)
{
  /* You MUST call free_op function here after an op is retired */ 
  /* you must complete the function */
  	if(MEM_latch->op_valid==true)
	{
		list<Op*>::iterator cii;
	for(cii = MEM_latch->queue_ops.begin(); cii != MEM_latch->queue_ops.end(); cii++)
	{
		Op *op = (*cii);
		if(op->cf_type>=CF_BR)
		{
			FE_latch->stage_stall = false;
		}
		if((op->dst!=-1) && (reg_writing_ops[(op)->dst] == op->inst_id))
		{
			register_file[op->dst].valid=true;
		}		
		
		if(op->inst_id == last_inst_id)
		{
			simulation_end_flag = 1;
		}
		free_op(op);
	}
		retired_instruction = retired_instruction + MEM_latch->queue_ops.size();
		MEM_latch->queue_ops.clear();
		MEM_latch->op_valid=false;
	}

		if(simulation_end_flag == 1 && main_memory->m_mshr.empty() && temp_ops_list_empty == 1) 
		{
			sim_end_condition = true;
		}
}

void MEM_stage(memory_c *main_memory)  // please modify MEM_stage function argument  /** NEW-LAB2 */ 
{
  /* you must complete the function */
	static int dcache_latency = 0;
	static int op_inst_id_store = -1;
  	if(EX_latch->op_valid==true)
	{
		/* Do everything*/
				
		if(MEM_latch->stage_stall==false) 		//check if stage is stalled
		{
			
			Op* op = EX_latch->op;
			if(op->mem_type > 0 )
			{
				if(dcache_latency == 0)
				{
					dcache_latency = KNOB(KNOB_DCACHE_LATENCY)->getValue();
				//Memory Instruction
				}
				dcache_latency--;
				if(dcache_latency == 0)
				{
					ADDRINT addr;
					if(op->mem_type == MEM_ST)
					{
						addr = op->st_vaddr;
					}					
					else
					{
					addr = op->ld_vaddr;
					}
					if(dcache_access(addr))
					{
						MEM_latch->queue_ops.push_back(op);
						MEM_latch->op_valid = TRUE;
						dcache_hit_count++;
						EX_latch->op_valid=false;
						EX_latch->stage_stall=false;
						
					}
					else
					{
						if(op_inst_id_store != op->inst_id)
						{
							dcache_miss_count++;
						}
						op_inst_id_store = op->inst_id;
						if(main_memory->store_load_forwarding(op))
						{
							MEM_latch->queue_ops.push_back(op);
							MEM_latch->op_valid = TRUE;
							EX_latch->op_valid=false;
							EX_latch->stage_stall=false;
							store_load_forwarding_count++;
						}
						else if(main_memory->store_store_forwarding(op))
						{
							MEM_latch->queue_ops.push_back(op);
							MEM_latch->op_valid = TRUE;
							EX_latch->op_valid=false;
							EX_latch->stage_stall=false;
							
						}
						else if(main_memory->check_piggyback(op))
						{
							EX_latch->op_valid=false;
							EX_latch->stage_stall=false;
						}
						else if(main_memory->insert_mshr(op))
						{
							EX_latch->op_valid=false;
							EX_latch->stage_stall=false;
						}
						else
						{
							EX_latch->stage_stall = true;
							return;
						}
					}
				}
				else
				{
					EX_latch->stage_stall = true;
				}
				
				
		}
		else	
		{
		EX_latch->op_valid=false;
		EX_latch->stage_stall=false;
		MEM_latch->op_valid=true;
		MEM_latch->queue_ops.push_back(EX_latch->op);
		}
	}
	else
	{
		EX_latch->stage_stall = true;
	}
}
}

void EX_stage() 
{
	static int ex_latency = 0;
  	if(ID_latch->op_valid==true)
	{
	

		if(EX_latch->stage_stall==false) 
		{
			
			if(ex_latency==0)
			{
				ex_latency=get_op_latency(ID_latch->op);	//check execution latency of op
			}
			ex_latency--;
			if(ex_latency==0)
			{
				ID_latch->op_valid=false;		//set op as invalid for previous stage since op is already transferred
				ID_latch->stage_stall=false;	//if stage stalled remove stall
				EX_latch->op=ID_latch->op;		//transfer pointer to next stage
				EX_latch->op_valid=true;		// set op as valid for the stage

			}
			else
			{
				ID_latch->stage_stall=true;
			}
		}
		else
		{
			ID_latch->stage_stall=true;
		}
	}
}

void ID_stage()
{
  /* you must complete the function */
#define NUM_SRC_1			((FE_latch->op)->num_src==1)
#define NUM_SRC_2			((FE_latch->op)->num_src==2)
#define SRC(x)				(FE_latch->op)->src[x]
#define SRC_INVALID(x)		!(register_file[SRC(x)].valid)




	if(FE_latch->op_valid == true)		//check if inst lies in the FE latch
	{
		if((NUM_SRC_1 && (SRC_INVALID(0))) || (NUM_SRC_2 && (SRC_INVALID(0) || SRC_INVALID(1))))	//check for source-dest clash
		{
			if(ID_latch->stage_stall == false)
			{
				if((FE_latch->op)->cf_type>=CF_BR)
				{
					control_hazard_count++;
				}
				data_hazard_count++;
			}
			FE_latch->stage_stall = true;
		}
		else
		{
			FE_latch->stage_stall = false;

			if(ID_latch->stage_stall == false)		//transfer op only if ID stage is not stalled
			{ 
				FE_latch->op_valid = false;
				ID_latch->op = FE_latch->op;
				ID_latch->op_valid = true;
				if((FE_latch->op)->dst != -1)			//if dest is valid then mark the specific register as busy/invalid and indicate which op writes to the reg
				{
					register_file[(FE_latch->op)->dst].valid = false;
					reg_writing_ops[(FE_latch->op)->dst] = FE_latch->op->inst_id; // in WB stage, clear valid flag of register only if this is the op
				}
				//checking branch in the end is important
				if((FE_latch->op)->cf_type>=CF_BR)
				{
					control_hazard_count++;
					FE_latch->stage_stall = true;
				}
			}
			else
			{
				FE_latch->stage_stall = true;
			}
		}	
	}
	
}


void FE_stage()
{
  /* only part of FE_stage function is implemented */ 
  /* please complete the rest of FE_stage function */ 

	if(FE_latch->stage_stall==false)	//check if no pipeline stalled
	{
		Op* new_op = get_free_op();	//get a placeholder op out of the pool
		if(get_op(new_op))		//copy op from trace into the placeholder op
		{
			FE_latch->op=new_op;
			FE_latch->op_valid=true;
		}
		else
		{
			free_op(new_op);
		}
	}
 
  //   next_pc = pc + op->inst_size;  // you need this code for building a branch predictor 

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
  
  MEM_latch->stage_stall = false;
  EX_latch->stage_stall = false;
  ID_latch->stage_stall = false;
  FE_latch->stage_stall = false;
}

// My functions //

void mem_latch_queue_ops()
{
	        while(!temp_ops.empty())
        	{
                      Op* op = temp_ops.front();
                      temp_ops.pop_front();
                      MEM_latch->queue_ops.push_back(op);
                      MEM_latch->op_valid = TRUE;
		      
        	}
	temp_ops_list_empty = 1;
	return;
}
bool icache_access(ADDRINT addr) {   /** please change uint32_t to ADDRINT NEW-LAB2 */ 

  /* For Lab #1, you assume that all I-cache hit */     
  bool hit = FALSE; 
  if (KNOB(KNOB_PERFECT_ICACHE)->getValue()) hit = TRUE; 
  return hit; 
}


bool dcache_access(ADDRINT addr) { /** please change uint32_t to ADDRINT NEW-LAB2 */ 
  /* For Lab #1, you assume that all D-cache hit */     
  /* For Lab #2, you need to connect cache here */   // NEW-LAB2 
  bool hit = FALSE;
  if (KNOB(KNOB_PERFECT_DCACHE)->getValue())
  {
	 hit = TRUE;
  return hit;
  }
  else
  {
	
	return cache_access(data_cache, addr);
	
  }
}



// NEW-LAB2 
void dcache_insert(ADDRINT addr)  // NEW-LAB2 
{                                 // NEW-LAB2 
  /* dcache insert function */   // NEW-LAB2 
  cache_insert(data_cache, addr) ;   // NEW-LAB2 
 
}                                       // NEW-LAB2 

void broadcast_rdy_op(Op* op)             // NEW-LAB2 
{                                          // NEW-LAB2 
  /* you must complete the function */     // NEW-LAB2 
  // mem ops are done.  move the op into WB stage   // NEW-LAB2 
  temp_ops.push_back(op);
	temp_ops_list_empty = 0;  
}      // NEW-LAB2 



/* utility functions that you might want to implement */     // NEW-LAB2 
int64_t get_dram_row_id(ADDRINT addr)    // NEW-LAB2 
{  // NEW-LAB2 
 // NEW-LAB2 
/* utility functions that you might want to implement */     // NEW-LAB2 
/* if you want to use it, you should find the right math! */     // NEW-LAB2 
/* pleaes carefull with that DRAM_PAGE_SIZE UNIT !!! */     // NEW-LAB2 
  // addr >> 6;   // NEW-LAB2 
  return ((addr)/((KNOB(KNOB_DRAM_PAGE_SIZE)->getValue())*1024));;   // NEW-LAB2 
}  // NEW-LAB2 

int get_dram_bank_id(ADDRINT addr)  // NEW-LAB2 
{  // NEW-LAB2 
 // NEW-LAB2 
/* utility functions that you might want to implement */     // NEW-LAB2 
/* if you want to use it, you should find the right math! */     // NEW-LAB2 

  // (addr >> 6);   // NEW-LAB2 
  return ((addr/((KNOB(KNOB_DRAM_PAGE_SIZE)->getValue())*1024))%(KNOB(KNOB_DRAM_BANK_NUM)->getValue())) ;   // NEW-LAB2 
}  // NEW-LAB2 

void init_register_file()
{
	for(int i=0;i<NUM_REG;i++)
		register_file[i].valid=true;		//all registers are initially available
}



