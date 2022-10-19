package comparison.kernel.graph;

import comparison.kernel.Kernel;
import graph.LGraph;

public abstract class ProductGraphBuilder<V,E,O extends ProductGraph> {
	
	Kernel<? super E> edgeKernel;
	Kernel<? super V> vertexKernel;
	
	/**
	 * Creates a new instance that can be used to subsequently build multiple
	 * product graphs. 
	 * @param edgeKernel kernel used to compare edges
	 * @param vertexKernel kernel used to compare vertices
	 */
	public ProductGraphBuilder(Kernel<? super E> edgeKernel, Kernel<? super V> vertexKernel) {
		this.edgeKernel = edgeKernel;
		this.vertexKernel = vertexKernel;
	}
	
	public void setVertexKernel(Kernel<? super V> vertexKernel) {
		this.vertexKernel = vertexKernel;		
	}
	
	public Kernel<? super V> getVertexKernel() {
		return vertexKernel;
	}
	
	public void setEdgeKernel(Kernel<? super E> edgeKernel) {
		this.edgeKernel = edgeKernel;		
	}

	public Kernel<? super E> getEdgeKernel() {
		return edgeKernel;
	}	
	
	public abstract O build(LGraph<V, E> lg1, LGraph<V, E> lg2);
	

	// TODO still required!?
//	protected double[] initVertexWeights(int size, boolean fillZero) {
//		if (vertexWeights != null && vertexWeights.length >= size) {
//			if (fillZero) {
//				Arrays.fill(vertexWeights, 0d);
//			}
//			return vertexWeights;
//		} else {
//			return new double[size];
//		}
//	}
//
//	protected double[][] initEdgeWeights(int size, boolean fillZero) {
//		if (edgeWeights != null && edgeWeights.length >= size) {
//			if (fillZero) {
//				for (int i=0; i<size; i++) {
//					Arrays.fill(edgeWeights, 0d);
//				}
//			}
//			return edgeWeights;
//		} else {
//			return new double[size][size];
//		}
//	}
}
