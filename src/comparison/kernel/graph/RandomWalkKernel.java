package comparison.kernel.graph;

import comparison.kernel.Kernel;

/**
 * Base class for random walk kernels inspired by the approach proposed by 
 * Gaertner et al. (2003), which originally requires vertex and edge labels 
 * to be identical.
 * This version is extended to arbitrary vertex and edge kernels by using
 * a weighted direct product graph.
 * 
 * @author kriege
 *
 * @param <V> vertex label type
 * @param <E> edge label type
 */
public abstract class RandomWalkKernel<V,E> extends ProductGraphKernel<V, E, ProductGraph> {
	
	/**
	 * @param edgeKernel edge comparison function used to compute weight matrix
	 * @param vertexKernel vertex comparison function used to compute weight matrix
	 */
	public RandomWalkKernel(Kernel<? super E> edgeKernel, Kernel<? super V> vertexKernel) {
//		super(new WeightedDirectProductGraphBuilderDense<V, E>(edgeKernel, vertexKernel));
		super(new WeightedDirectProductGraphBuilderSparse<V, E>(edgeKernel, vertexKernel));
	}

	/**
	 * {@inheritDoc}
	 */
	abstract double computeKernel(ProductGraph pg);
	

}
