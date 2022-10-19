package comparison.kernel.graph;

import comparison.kernel.Kernel;
import graph.LGraph;

/**
 * Base  class for kernels based on product graphs.
 * 
 * @author kriege
 *
 * @param <V> vertex label type
 * @param <E> edge label type
 * @param <P> product graph type
 */
//TODO support composite kernels (most product graph kernels are composite kernels!)
// Note: This class should also be used by the subgraph matching 
// kernel family
public abstract class ProductGraphKernel<V, E, P extends ProductGraph> 
	implements Kernel<LGraph<V, E>>
{
	public long computations = 0;
	public long productGraphComputationTime = 0;
	public long kernelComputationTime = 0;
	
	ProductGraphBuilder<V, E, P> pgb;
	
	public ProductGraphKernel(ProductGraphBuilder<V, E, P> pgb) {
		this.pgb = pgb;
	}
	
	public Kernel<? super V> getVertexKernel() {
		return pgb.getVertexKernel();
	}
	
	public Kernel<? super E> getEdgeKernel() {
		return pgb.getEdgeKernel();
	}

	public double compute(LGraph<V, E> g1, LGraph<V, E> g2) {
		computations++;
		
		long startTime = System.nanoTime();
		P pg = pgb.build(g1, g2);
		productGraphComputationTime += System.nanoTime() - startTime;

		startTime = System.nanoTime();
		double r = computeKernel(pg);
		kernelComputationTime += System.nanoTime() - startTime;
		
		return r;
	}
	
	/**
	 * Computes the kernel value based on the given weighted product 
	 * graph.
	 * 
	 * Note: Implementing classes should use both, vertex and edge
	 * weights, to compute a weight matrix, where each element
	 * encodes the product of an edge kernel result and two vertex
	 * kernel results.
	 * 
	 * @param pg weighted product graph encoding the vertex and edge 
	 * kernel results for computation
	 * @return kernel value
	 */
	abstract double computeKernel(P pg);
	
	@Override
	public String toString() {
		return "Total number of computations: "+computations+"\n"+
			"Total computation time: "+(productGraphComputationTime+kernelComputationTime)/1000/1000+" [ms]\n"+
			"Total product graph time: "+productGraphComputationTime/1000/1000+" [ms]\n"+
			"Total kernel computation time: "+kernelComputationTime/1000/1000+" [ms]\n";
	}

}
