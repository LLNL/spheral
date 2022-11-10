#ifndef EXECSTRATEGY_HH
#define EXECSTRATEGY_HH
class ExecutionStrategy {
  using T = unsigned int;
  using ExecCtx = RAJA::ExecPlace;

public:
  ExecutionStrategy(T n_p, T d_sz, ExecCtx ctx = ExecCtx::HOST) {
    n_pairs = n_p;
    data_sz = d_sz;
    exec_context = ctx;

    //// Host attempts to balance based on number of omp threads being used.
    //// and adjusts block size accordingly.
    //if (ctx == ExecCtx::HOST){
    //  blocks_per_pool = default_blocks_per_pool;
    //  //n_blocks = omp_get_max_threads();

    //  block_sz = RAJA_DIVIDE_CEILING_INT(n_pairs, n_blocks);
    //  n_data_pools = RAJA_DIVIDE_CEILING_INT(n_blocks, blocks_per_pool);

    //  if (n_blocks == 1) sequential = true;
    //}

    //// Device favours defined block sizes and multiple blocks per data
    //// pool. 
    //if (ctx == ExecCtx::DEVICE){
    //  block_sz = default_blocks_sz;
    //  n_blocks = RAJA_DIVIDE_CEILING_INT(n_pairs, block_sz);
    //  n_data_pools = default_n_data_pools;
    //  blocks_per_pool = n_blocks / n_data_pools;
    //}
  
  }

  void print(){
    std::cout << "**********************************\n";
    std::cout << "Execution Strategy ";
    if (exec_context == ExecCtx::HOST) std::cout << "-- HOST --\n";
    if (exec_context == ExecCtx::DEVICE) std::cout << "-- DEVICE --\n";
    std::cout << "**********************************\n";
    std::cout << " - n_pairs  : " << n_pairs << "\n";
    std::cout << " - data_sz  : " << data_sz << "\n";
    std::cout << std::endl;
    std::cout << " - block_sz : " << block_sz << "\n";
    std::cout << " - n_blocks : " << n_blocks << "\n";
    std::cout << std::endl;
    std::cout << " - blk/pool : " << blocks_per_pool << "\n";
    std::cout << " - n_pools  : " << n_data_pools << "\n";
    std::cout<< "**********************************\n";
  }

  T n_pairs, data_sz;
  
  T block_sz, n_blocks;
  T blocks_per_pool, n_data_pools;

  bool sequential = false;

  ExecCtx exec_context = ExecCtx::HOST;

private:
  static constexpr T default_blocks_per_pool = 1;
  static constexpr T default_blocks_sz = 512;
  static constexpr T default_n_data_pools = 16;

};


#endif //  EXECSTRATEGY_HH
