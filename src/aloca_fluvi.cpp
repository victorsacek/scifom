#include <stdio.h>
#include <stdlib.h>

double **Aloc_matrix_real (long p, long n);
long **Aloc_matrix_long (long p, long n);
double *Aloc_vector_real (long n);
long *Aloc_vector_long (long n);

extern long nodes_max_aloca;

extern double *Df;
extern long *ordem_fluvi;
extern long *direc_fluvi;
extern double *dist_fluvi;
extern double *Qr;
extern double *Qr_prov;
extern double *Qf;
extern long *basins;
extern long **global_basin;
extern long **global_index;
extern double **h_global_basin;
extern long *n_global;
extern long *global_flag;
extern double *h_min;
extern long *pos_min;

extern long *stack_fluvi;

extern double *h_topo_prov;

extern long *lagos;

void aloca_fluvi()
{
     Df = Aloc_vector_real (nodes_max_aloca);
     
     ordem_fluvi = Aloc_vector_long (nodes_max_aloca);
     direc_fluvi = Aloc_vector_long (nodes_max_aloca);
     dist_fluvi = Aloc_vector_real (nodes_max_aloca);
     
     Qr = Aloc_vector_real (nodes_max_aloca);
     Qr_prov = Aloc_vector_real (nodes_max_aloca);	
     Qf = Aloc_vector_real (nodes_max_aloca); 
     
     basins = Aloc_vector_long (nodes_max_aloca);
     global_basin = Aloc_matrix_long (nodes_max_aloca,200);
	 global_index = Aloc_matrix_long (nodes_max_aloca,200);
	 n_global = Aloc_vector_long (nodes_max_aloca);	
	 global_flag = Aloc_vector_long (nodes_max_aloca);	
	 h_global_basin = Aloc_matrix_real (nodes_max_aloca,200);
	
	stack_fluvi = Aloc_vector_long(nodes_max_aloca);
	
	lagos = Aloc_vector_long(nodes_max_aloca);
	
	h_topo_prov = Aloc_vector_real(nodes_max_aloca);
	
     h_min = Aloc_vector_real (nodes_max_aloca);
     pos_min = Aloc_vector_long (nodes_max_aloca);
     
}
