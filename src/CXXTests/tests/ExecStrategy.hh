#ifndef EXECSTRATEGY_HH
#define EXECSTRATEGY_HH

class ExecutionStrategy {
  using T = unsigned int;
  using Platform = RAJA::Platform;

public:
  ExecutionStrategy(T n_p, T d_sz, Platform p = Platform::host) :
    n_pairs(n_p),
    data_sz(d_sz),
    platform(p) {

    if (platform == Platform::cuda) {
      block_sz = 1024;
      n_blocks = RAJA_DIVIDE_CEILING_INT(n_pairs, block_sz);

      n_data_pools = 16;
      pool_block_sz = RAJA_DIVIDE_CEILING_INT(n_pairs, n_data_pools);
    
    } else { // Platform::host
      block_sz = n_pairs / omp_get_max_threads();
      n_blocks = RAJA_DIVIDE_CEILING_INT(n_pairs, block_sz);

      n_data_pools = n_blocks;
      pool_block_sz = block_sz;
    }
  
  }

  void print(){
    std::cout << "**********************************\n";
    std::cout << "Execution Strategy ";
    if (platform == Platform::host) std::cout << "-- HOST --\n";
    if (platform == Platform::cuda) std::cout << "-- DEVICE --\n";
    std::cout << "**********************************\n";
    std::cout << " - n_pairs       : " << n_pairs << "\n";
    std::cout << " - data_sz       : " << data_sz << "\n";
    std::cout << std::endl;
    std::cout << " - block_sz      : " << block_sz << "\n";
    std::cout << " - n_blocks      : " << n_blocks << "\n";
    std::cout << std::endl;
    std::cout << " - pool_block_sz : " << pool_block_sz << "\n";
    std::cout << " - n_pools       : " << n_data_pools << "\n";
    std::cout<< "**********************************\n";
  }

  T n_pairs, data_sz;
  
  T block_sz, n_blocks;
  T pool_block_sz, n_data_pools;

  Platform platform = Platform::host;
  Platform host_platform = Platform::host;
};

#endif //  EXECSTRATEGY_HH
