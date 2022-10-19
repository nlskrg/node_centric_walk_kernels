package comparison.kernel.graph;

import java.util.HashMap;
import java.util.Map.Entry;

import comparison.kernel.ExplicitMappingKernel;
import comparison.kernel.Kernel;
import comparison.kernel.basic.DiracKernel;
import comparison.kernel.graph.ProductGraph.PEdge;
import comparison.kernel.graph.ProductGraph.PVertex;
import datastructure.FeatureVector;
import datastructure.SparseFeatureVector;
import graph.Graph;
import graph.Graph.Edge;
import graph.Graph.Vertex;
import graph.LGraph;
import graph.properties.EdgeArray;
import graph.properties.VertexArray;

/**
 * Implementation of a fixed-length random walk kernel especially suitable
 * for sparse product graphs. A walk (v0, e0, v1, ..., en-1, vn) contributes
 * weight vk(v0) * ek(e0) * vk(v1) * ... * ek(en-1) + vk(vn)
 * to the kernel value, where vertices and edges of the product graph 
 * correspond to matched vertices and edges of the two factor graphs.
 * 
 * This class supports explicit and implicit feature mapping for kernel computation.
 * 
 * Note: Explicit mapping is only supported for Dirac vertex and edge kernels;
 * labels must be uniquely identified by their toString() method. 
 * 
 * This class implements the two methods of computation presented in
 * Kriege et al. 2014, ICDM. 
 * 
 * @author kriege
 *
 * @param <V> vertex label type
 * @param <E> edge label type
 */
public class FixedLengthRandomWalkKernel<V,E> 
	extends RandomWalkKernel<V, E> 
	implements ExplicitMappingKernel<LGraph<V, E>, String> 
{
	
	public static final char DELIMITER = '|';
	
	int k;
	
	/**
	 * @param k maximum size of the random walk
	 * @param edgeKernel
	 * @param vertexKernel
	 */
	public FixedLengthRandomWalkKernel(int k, Kernel<? super E> edgeKernel, Kernel<? super V> vertexKernel) {
		super(edgeKernel, vertexKernel);
		this.k = k;
	}

	@Override
	double computeKernel(ProductGraph pg) {

		// no need to copy values in this case
		if (k==0) {
			double r = 0;
			for (PVertex v : pg.vertices()) {
				r += v.getWeight();
			}
			return r;
		}
		
		double[] weights = new double[pg.getVertexCount()];
		double[] newWeights = new double[pg.getVertexCount()];
		for (PVertex v : pg.vertices()) {
			weights[v.getIndex()] = v.getWeight(); 
		}
		
		for (int i=1; i<=k; i++) {
			// compute new weights
			for (PVertex v : pg.vertices()) {
				newWeights[v.getIndex()] = computeNewWeight(v, weights);
			}
			// swap weight arrays
			double[] tmp = weights;
			weights = newWeights;
			newWeights = tmp;
		}
		
		double r = 0;
		for (PVertex v : pg.vertices()) {
			r += weights[v.getIndex()];
		}

		return r;
	}
	
	private double computeNewWeight(PVertex v, double[] weights) {
		double r = 0;
		for (PEdge e : v.edges()) {
			PVertex w = e.getOppositeVertex(v);
			// TODO edge vertex product can be computed once and then stored
			// use a directed graph as product graph
			r += weights[w.getIndex()] * e.getWeight() * v.getWeight(); 
		}
		
		return r;
	}
	
	public String getID() {
		if (pgb.getVertexKernel() instanceof DiracKernel && pgb.getEdgeKernel() instanceof DiracKernel) {
			return "FLRW_"+k;
		} else {
			return "FLRW_"+k+"_"+pgb.getVertexKernel().getID()+"_"+pgb.getEdgeKernel().getID();
		}
	}

	//--------------------------------------------
	// EXPLICIT COMPUTATION
	//--------------------------------------------
	
	@Override
	public FeatureVector<String> getFeatureVector(LGraph<V, E> lg) throws IllegalStateException {
		if (!(pgb.getVertexKernel() instanceof DiracKernel) || !(pgb.getEdgeKernel() instanceof DiracKernel)) {
			throw new IllegalStateException("Explicit mapping requires both, vertex and edge kernel, to be dirac kernels!");
		}
		
//		WalkIntegerConversion wic = new WalkIntegerConversion();
		return computeLabeledWalkCount(lg/*, wic*/);
	}

	private SparseFeatureVector<String> computeLabeledWalkCount(LGraph<V,E> lg/*, WalkIntegerConversion wic*/) {
		Graph g = lg.getGraph();
		VertexArray<V> va = lg.getVertexLabel();
		EdgeArray<E> ea = lg.getEdgeLabel();
		VertexArray<SparseFeatureVector<String>> walks = new VertexArray<SparseFeatureVector<String>>(g);
		VertexArray<SparseFeatureVector<String>> newWalks = new VertexArray<SparseFeatureVector<String>>(g);
		
		// initialize with walks of size 0
		for (Vertex v : g.vertices()) {
			SparseFeatureVector<String> f = new SparseFeatureVector<String>();
			walks.set(v, f);
			String walk = String.valueOf(va.get(v));
			f.increaseByOne(/*wic.getInteger(*/walk/*)*/);
		}
		
		for (int i=1; i<=k; i++) {
			// compute new walks
			newWalks = new VertexArray<SparseFeatureVector<String>>(g);
			for (Vertex v : g.vertices()) {
				// compute new walks
				SparseFeatureVector<String> f = new SparseFeatureVector<String>();
				newWalks.set(v, f);
				for (Edge e : v.edges()) {
					Vertex w = e.getOppositeVertex(v);
					for (Entry<String,Double> wie : walks.get(w).nonZeroEntries()) {
						String wString = wie.getKey();/*wic.getWalk(wi);*/
						String newWString = String.valueOf(va.get(v))
								+DELIMITER+String.valueOf(ea.get(e))
								+DELIMITER+wString;
						f.increase(/*wic.getInteger(*/newWString/*)*/, wie.getValue());
					}
				}
			}
			walks = newWalks;
		}
		
		SparseFeatureVector<String> r = new SparseFeatureVector<String>();
		for (Vertex v : g.vertices()) {
			r.add(walks.get(v));
		}

		return r;
	}
	
	
	public static class WalkIntegerConversion {
		HashMap<String, Integer> w2i = new HashMap<String, Integer>();
		HashMap<Integer, String> i2w = new HashMap<Integer, String>();

		public String getWalk(Integer i) {
			return i2w.get(i);
		}
		
		public Integer getInteger(String walk) {
			Integer i = w2i.get(walk);
			if (i == null) {
				i = w2i.size();
				w2i.put(walk, i);
				i2w.put(i, walk);
			}
			return i;
		}
		
		public void clear() {
			w2i.clear();
			i2w.clear();
		}
	}

	
}
