package comparison.kernel.graph;


import java.util.ArrayList;
import java.util.List;
import java.util.Map.Entry;

import algorithm.shortestpath.betweenness.BrandesEdgeLengthHopper;
import comparison.kernel.ExplicitMappingKernel;
import comparison.kernel.Kernel;
import datastructure.FeatureVector;
import datastructure.Matrix;
import datastructure.MatrixTool;
import datastructure.SparseFeatureVector;
import datastructure.Triple;
import graph.Graph.Vertex;
import graph.LGraph;
import graph.properties.EdgeArray;
import graph.properties.VertexArray;

/**
 * Computes the GraphHopper kernel for attributed graphs. Our implementation
 * exploits ideas from centrality computation and differs from the authors
 * implementation, but gives the same results.
 * 
 * Feragen, A.; Kasenburg, N.; Petersen, J.; Bruijne, M. D. & Borgwardt
 * Scalable kernels for graphs with continuous attributes 
 * NIPS 2013
 * 
 * @author Nils Kriege
 *
 * @param <V>
 * @param <E>
 */
public class GraphHopper<V,E> implements Kernel<LGraph<V, E>>, ExplicitMappingKernel<LGraph<V, E>, Triple<Integer,Integer,?>>{
	
	private Kernel<? super V> vertexKernel;
	private boolean useEdgeLength;
	
	/**
	 * 
	 * @param vertexKernel
	 * @param useEdgeLength if true edge labels are used to compute path length; 
	 * supports integer and double edge length; doubles are discretized to avoid 
	 * numerical problems, see {@link BrandesEdgeLengthHopper#discretizeDouble(graph.Graph, EdgeArray)}}
	 */
	public GraphHopper(Kernel<? super V> vertexKernel, boolean useEdgeLength) {
		this.vertexKernel = vertexKernel;
		this.useEdgeLength = useEdgeLength;
	}

	@Override
	public double compute(LGraph<V, E> g1, LGraph<V, E> g2) {
		BrandesEdgeLengthHopper belh = new BrandesEdgeLengthHopper();
		
		EdgeArray<Integer> ea1=getEdgeLength(g1);
		EdgeArray<Integer> ea2=getEdgeLength(g2);

		VertexArray<Matrix<Long>> M1 = belh.compute(g1.getGraph(), ea1);
		VertexArray<Matrix<Long>> M2 = belh.compute(g2.getGraph(), ea2);
		
		return compute(g1, g2, M1, M2);
	}
	
	@SuppressWarnings("unchecked")
	public EdgeArray<Integer> getEdgeLength(LGraph<V, ?> g) {
		if (!useEdgeLength) return null;
		if (g.getGraph().getEdgeCount()==0) return null;
		if (g.getEdgeLabel().get(g.getGraph().getEdge(0)) instanceof Integer) {
			return (EdgeArray<Integer>) g.getEdgeLabel();
		}
		if (g.getEdgeLabel().get(g.getGraph().getEdge(0)) instanceof Double) {
			// discretize
			return BrandesEdgeLengthHopper.discretizeDouble(g.getGraph(), 
					(EdgeArray<Double>) g.getEdgeLabel());
		}
		throw new IllegalArgumentException("Unsupported edge length type! Use Integer or Double!");
	}
	
	public double compute(LGraph<V, E> g1, LGraph<V, E> g2,
			VertexArray<Matrix<Long>> M1, VertexArray<Matrix<Long>> M2) {
		
		double val = 0; 
		
		for (Vertex v : g1.getGraph().vertices()) {
			Matrix<Long> M1v = M1.get(v);
			for (Vertex w : g2.getGraph().vertices()) {
				Matrix<Long> M2w = M2.get(w);
				int rowDim = Math.min(M1v.getRowDimension(), M2w.getRowDimension());
				int colDim = Math.min(M1v.getColumnDimension(), M2w.getColumnDimension());
				long hits = 0;
				for (int i=0; i<rowDim; i++) {
					for (int j=0; j<colDim; j++) {
						hits += M1v.get(i, j)*M2w.get(i, j);
					}
				}
				//vertexkernel
				val += vertexKernel.compute(g1.getVertexLabel().get(v), g2.getVertexLabel().get(w))*hits;
			}
		}
		
		return val;
	}


	@Override
	public String getID() {
		return "GH_"+vertexKernel.getID();
	}
	
	
	public double[][] computeAll(List<? extends LGraph<V, E>> graphs) {
		
		BrandesEdgeLengthHopper belh = new BrandesEdgeLengthHopper();
		
		ArrayList<VertexArray<Matrix<Long>>> allM = new ArrayList<VertexArray<Matrix<Long>>>(graphs.size());
		for (LGraph<V, ?> lg : graphs) {
			EdgeArray<Integer> ea=getEdgeLength(lg);
			VertexArray<Matrix<Long>> M = belh.compute(lg.getGraph(), ea);
			for (Vertex v : lg.getGraph().vertices()) {
				M.set(v, MatrixTool.truncate(M.get(v)));
			}
			allM.add(M);
		}
		
		int n = graphs.size();
		double[][] gram = new double[n][n];
		
		for (int i=0; i<n; i++) {
			LGraph<V, E> g1 = graphs.get(i);
			for (int j=i; j<n; j++) {
				LGraph<V, E> g2 = graphs.get(j);
				gram[i][j] = gram[j][i] = compute(g1, g2, allM.get(i), allM.get(j));
			}
		}
		
		return gram;
	}

	@Override
	public FeatureVector<Triple<Integer,Integer,?>> getFeatureVector(LGraph<V, E> t) throws IllegalStateException {
		if (!(vertexKernel instanceof ExplicitMappingKernel<?, ?>)) {
			throw new IllegalStateException("Explicit mapping not supported by vertex kernel!");
		}
		
		BrandesEdgeLengthHopper belh = new BrandesEdgeLengthHopper();
		
		EdgeArray<Integer> ea=getEdgeLength(t);
		VertexArray<Matrix<Long>> M = belh.compute(t.getGraph(), ea);
		
		// build the feature map
		SparseFeatureVector<Triple<Integer,Integer,?>> r = new SparseFeatureVector<Triple<Integer,Integer,?>>();
		@SuppressWarnings("unchecked")
		ExplicitMappingKernel<? super V,E> vk = (ExplicitMappingKernel<? super V, E>) vertexKernel;
		VertexArray<V> va = t.getVertexLabel();
		for (Vertex v : t.getGraph().vertices()) {
			for (Entry<?, Double> fu : vk.getFeatureVector(va.get(v)).nonZeroEntries()) {
				Matrix<Long> Mv = M.get(v);
				for (int i=0; i<Mv.getRowDimension(); i++) {
					for (int j=0; j<Mv.getRowDimension(); j++) {
						r.increase(new Triple<>(i, j, fu.getKey()), fu.getValue()*Mv.get(i, j));
					}
				}
			}
		}

		return r;
	}
	
}
