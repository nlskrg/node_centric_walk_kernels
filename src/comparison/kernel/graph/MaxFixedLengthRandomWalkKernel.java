package comparison.kernel.graph;

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
 * FLRW until a fixed length; this kernel is equivalent to the sum of
 * {@link FixedLengthRandomWalkKernel} for 0..k, using a more efficient
 * method of computation.
 * 
 * 
 * @author kriege
 *
 * @param <V> vertex label type
 * @param <E> edge label type
 */
public class MaxFixedLengthRandomWalkKernel<V,E> 
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
	public MaxFixedLengthRandomWalkKernel(int k, Kernel<? super E> edgeKernel, Kernel<? super V> vertexKernel) {
		super(edgeKernel, vertexKernel);
		this.k = k;
	}

	@Override
	double computeKernel(ProductGraph pg) {

		double r = 0;
		
		// no need to copy values in this case
		if (k==0) {
			for (PVertex v : pg.vertices()) {
				r += v.getWeight();
			}
			return r;
		}
		
		double[] weights = new double[pg.getVertexCount()];
		double[] newWeights = new double[pg.getVertexCount()];
		for (PVertex v : pg.vertices()) {
			weights[v.getIndex()] = v.getWeight(); 
			r += v.getWeight();
		}
		
		for (int i=1; i<=k; i++) {
			// compute new weights
			for (PVertex v : pg.vertices()) {
				double w =  computeNewWeight(v, weights);
				newWeights[v.getIndex()] = w;
				r += w;
			}
			// swap weight arrays
			double[] tmp = weights;
			weights = newWeights;
			newWeights = tmp;
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
			return "MFLRW_"+k; 
		} else {
			return "MFLRW_"+k+"_"+pgb.getVertexKernel().getID()+"_"+pgb.getEdgeKernel().getID();
		}
	}
	
	//--------------------------------------------
	// EXPLICIT COMPUTATION
	//--------------------------------------------
	
	@Override
	public FeatureVector<String> getFeatureVector(LGraph<V, E> in) throws IllegalStateException {
		if (!(pgb.getVertexKernel() instanceof DiracKernel) || !(pgb.getEdgeKernel() instanceof DiracKernel)) {
			throw new IllegalStateException("Explicit mapping requires both, vertex and edge kernel, to be dirac kernels!");
		}
		
		return computeLabeledWalkCount(in);
	}
	
	private SparseFeatureVector<String> computeLabeledWalkCount(LGraph<V,E> lg/*, WalkIntegerConversion wic*/) {
		Graph g = lg.getGraph();
		VertexArray<V> va = lg.getVertexLabel();
		EdgeArray<E> ea = lg.getEdgeLabel();
		VertexArray<SparseFeatureVector<String>> walks = new VertexArray<SparseFeatureVector<String>>(g);
		VertexArray<SparseFeatureVector<String>> newWalks = new VertexArray<SparseFeatureVector<String>>(g);

		SparseFeatureVector<String> r = new SparseFeatureVector<String>();
		
		// initialize with walks of size 0
		for (Vertex v : g.vertices()) {
			SparseFeatureVector<String> f = new SparseFeatureVector<String>();
			walks.set(v, f);
			String walk = String.valueOf(va.get(v));
			f.increaseByOne(walk);
			r.add(f);;
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
						String wString = wie.getKey();
						String newWString = String.valueOf(va.get(v))
								+DELIMITER+String.valueOf(ea.get(e))
								+DELIMITER+wString;
						f.increase(newWString, wie.getValue());
					}
				}
				r.add(f);
			}
			walks = newWalks;
		}

		return r;
	}


}
