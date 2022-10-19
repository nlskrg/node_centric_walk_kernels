package comparison.kernel.graph;

import java.util.ArrayList;
import java.util.List;
import java.util.Map.Entry;

import comparison.kernel.Kernel;
import comparison.kernel.basic.DiracKernel;
import comparison.kernel.graph.ProductGraph.PEdge;
import comparison.kernel.graph.ProductGraph.PVertex;
import datastructure.SparseFeatureVector;
import graph.Graph;
import graph.LGraph;
import graph.Graph.Edge;
import graph.Graph.Vertex;
import graph.properties.EdgeArray;
import graph.properties.VertexArray;


/**
 * The node centric walk kernel.
 * @author kriege
 *
 * @param <V>
 * @param <E>
 */
public class NodeCentricWalkKernel<V,E> implements Kernel<LGraph<V, E>> {
	ProductGraphBuilder<V, E, ProductGraph> pgb;
	int k;
	double alpha;
	double beta;
	boolean wl;
	Standardizer std;
	
	public static final char DELIMITER = '|';
	
	/**
	 * Create the node centric walk kernel.
	 * 
	 * @param k number of iterations/walk length
	 * @param alpha neighborhood strictness
	 * @param beta walk count influence
	 * @param edgeKernel
	 * @param vertexKernel
	 * @param wl enable WL expressive power
	 */
	public NodeCentricWalkKernel(int k, double alpha, double beta, Kernel<? super E> edgeKernel, Kernel<? super V> vertexKernel, boolean wl) {
		pgb = new WeightedDirectProductGraphBuilderSparse<V, E>(edgeKernel, vertexKernel);
		this.k = k;
		this.wl = wl;
		if (alpha==Double.POSITIVE_INFINITY && beta == 0d) {
			this.std = new Standardizer.NodeKernelAlphaInftyBetaZero();
		} else if (alpha==Double.POSITIVE_INFINITY && beta == 1d) {
			this.std = new Standardizer.NodeKernelAlphaInftyBetaOne();
		} else if (alpha==Double.POSITIVE_INFINITY) {
			this.std = new Standardizer.NodeKernelAlphaInfty(beta);
		} else {
			this.std= new Standardizer.NodeKernel(alpha,beta);
		}
	}

	private double computeNewWeight(PVertex v, double[] weights) {
		double r = 0;
		for (PEdge e : v.edges()) {
			PVertex w = e.getOppositeVertex(v);
			r += weights[w.getIndex()] * e.getWeight() * v.getWeight(); 
		}
		
		return r;
	}
	
	private double[] computeNewWeights(ProductGraph pg, double[] weights) {
		
		double[] newWeights = new double[pg.getVertexCount()];
		// compute new weights
		for (PVertex v : pg.vertices()) {
			double w =  computeNewWeight(v, weights);
			newWeights[v.getIndex()] = w;
		}
		
		return newWeights;
	}
	
	private double[] initialWeights(ProductGraph pg) {
		double[] weights = new double[pg.getVertexCount()];
		for (PVertex v : pg.vertices()) {
			weights[v.getIndex()] = v.getWeight(); 
		}
		return weights;
	}
	
	private int[] selfSimIndex(ProductGraph pg, int n) {
		int[] selfSimIndex = new int[n];
		for (PVertex pv : pg.vertices()) {
			if (pv.getFirst() == pv.getSecond()) {
				selfSimIndex[pv.getFirst().getIndex()] = pv.getIndex();
			}
		}
		return selfSimIndex;
	}
	
	@Override
	public double[][] computeAll(List<? extends LGraph<V, E>> set) {

		int n = set.size();

		// compute all self-similarities
		ArrayList<double[][]> selfSims = new ArrayList<double[][]>(n);
		for (int i=0; i<n; i++) {
			selfSims.add(computeSelfSimilarities(set.get(i)));
		}
		
		// compute kernel values
		double[][] r = new double[n][n];
		for (int i=0; i<n; i++) {
			LGraph<V, E> e1 = set.get(i);
			for (int j=i; j<n; j++) {
				LGraph<V, E> e2 = set.get(j);
				r[i][j] = r[j][i] = this.compute(e1, e2, selfSims.get(i), selfSims.get(j));
			}
		}
		return r;
	}
	
	public double[][] computeSelfSimilarities(LGraph<V, E> lg) {
		int n = lg.getGraph().getVertexCount();
		ProductGraph pg = pgb.build(lg, lg);
		int[] selfSimIdx = selfSimIndex(pg, n);
		double[][] rW = new double[k+1][pg.getVertexCount()];
		double[] W = initialWeights(pg);
		double[] Wp = W.clone();
		rW[0] = W;
		for (int l=1; l<=k; l++) {
			// new weights
			W = computeNewWeights(pg, W);
			for (int z=0; z<pg.getVertexCount(); z++) { Wp[z] += W[z]; } // add up
			rW[l] = Wp.clone();
			if (wl) W = modifySelf(pg, selfSimIdx, Wp, W, std);
		}
		// collect self-similarities
		double[][] selfSims = new double[k+1][n];
		for (int l=0; l<=k; l++) {
			for (int i=0; i<n; i++) {
				selfSims[l][i] = rW[l][selfSimIdx[i]];
			}
		}
		
		return selfSims;
	}
	
	private double compute(LGraph<V, E> lg1, LGraph<V, E> lg2, double[][] selfSims1, double[][] selfSims2) {
		ProductGraph pg = pgb.build(lg1, lg2);
		
		double[][] rW = new double[k+1][pg.getVertexCount()];
		
		// initial weights
		double[] W = initialWeights(pg);
		double[] Wp = W.clone();
		rW[0] = gnwk(pg, selfSims1[0], selfSims2[0], Wp, W, std);
		
		for (int l=1; l<=k; l++) {
			// new weights
			W = computeNewWeights(pg, W);
			for (int z=0; z<pg.getVertexCount(); z++) { Wp[z] += W[z]; } // add up
			rW[l] = gnwk(pg, selfSims1[l], selfSims2[l], Wp, W, std);
			if (wl) W = rW[l];
		}
		
		double r = 0;
		for (PVertex v : pg.vertices()) {
			int idx = v.getIndex();

			for (int l=0; l<=k; l++) {
				r += rW[l][idx];
			}
			
		}
		
		return r;
	}
	
	public double compute(LGraph<V, E> lg1, LGraph<V, E> lg2) {
		
		double[][] selfSims1 = computeSelfSimilarities(lg1);
		double[][] selfSims2 = computeSelfSimilarities(lg2);
		
		return compute(lg1, lg2, selfSims1, selfSims2);
	}

	
	private double[] modifySelf(ProductGraph pg, int[] selfSim, double[] rNp, double[] rN, Standardizer std) {
		double[] R = new double[pg.getVertexCount()];
		for (PVertex v : pg.vertices()) {
			int idx = v.getIndex();
			int idx1 = selfSim[v.getFirst().getIndex()];
			int idx2 = selfSim[v.getSecond().getIndex()];
			
			double g = std.standardize(rNp[idx], rNp[idx1], rNp[idx2], rN[idx]);
			R[idx] = g;
		}
		return R;
	}
	
	private double[] gnwk(ProductGraph pg, double[] selfSims1, double[] selfSims2, double[] rNp, double[] rN, Standardizer std) {
		double[] R = new double[pg.getVertexCount()];
		for (PVertex v : pg.vertices()) {
			int idx = v.getIndex();
//			System.out.print(rNp[idx]+" "+selfSims1[v.getFirst().getIndex()]+" "+selfSims2[v.getSecond().getIndex()]+" "+rN[idx]);
			R[idx] = std.standardize(rNp[idx], selfSims1[v.getFirst().getIndex()], selfSims2[v.getSecond().getIndex()], rN[idx]);
//			System.out.println("->"+R[idx]);
		}
		return R;
	}
	
	@Override
	public String getID() {
		if (pgb.getVertexKernel() instanceof DiracKernel && pgb.getEdgeKernel() instanceof DiracKernel) {
			return "NCW"+(wl? "WL":"")+"_"+k+"_"+std.getID(); 
		} else {
			return "NCW"+(wl? "WL":"")+"_"+k+"_"+std.getID()+"_"+pgb.getVertexKernel().getID()+"_"+pgb.getEdgeKernel().getID();
		}
	}
	
	public static interface Standardizer {
		public double standardize(double abp, double aap, double bbp, double ab);
		public String getID();
		
		/**
		 * Gaussian and Polynomial kernel
		 */
		public class NodeKernel implements Standardizer {
			double alpha, beta;
			public NodeKernel(double alpha, double beta) {this.alpha = alpha; this.beta = beta;}
			public double standardize(double abp, double aap, double bbp, double ab) {
				// compute squared distance using kernel trick
				double d = aap - 2*abp + bbp;
				// compute Gaussian RBF kernel
				double g = Math.exp(-alpha*d);
				return g*Math.pow(ab, beta);
			}
			public String getID() { return "A"+alpha+"B"+beta; }
		}
		
		public class NodeKernelAlphaInftyBetaZero implements Standardizer {
			public NodeKernelAlphaInftyBetaZero() {}
			public double standardize(double abp, double aap, double bbp, double ab) {
				return (aap == abp && aap == bbp) ? 1d : 0d;
			}
			public String getID() { return "AinfB0"; }
		}
		
		public class NodeKernelAlphaInftyBetaOne implements Standardizer {
			public NodeKernelAlphaInftyBetaOne() {}
			public double standardize(double abp, double aap, double bbp, double ab) {
				return (aap == abp && aap == bbp) ? ab : 0d;
			}
			public String getID() { return "AinfB1"; }
		}
		public class NodeKernelAlphaInfty implements Standardizer {
			double beta;
			public NodeKernelAlphaInfty(double beta) {this.beta = beta;}
			public double standardize(double abp, double aap, double bbp, double ab) {
				return (aap == abp && aap == bbp) ? Math.pow(ab, beta) : 0d;
			}
			public String getID() { return "AinfB"+beta; }
		}

	}
	
	
	//--------------------------------------------
	// Computation using walk counts of input graphs 
	//--------------------------------------------
	public double[][] computeAllOnInputGraphs(List<? extends LGraph<V, E>> set) {
//		if (!(pgb.getVertexKernel() instanceof DiracKernel) || !(pgb.getEdgeKernel() instanceof DiracKernel) || (wl) ) {
//			throw new IllegalStateException("Explicit mapping requires both, vertex and edge kernel, to be dirac kernels and does not support wl expressive power!");
//		}
		
		ArrayList<ArrayList<VertexArray<SparseFeatureVector<String>>>> allLabeledWalkCounts = new ArrayList<ArrayList<VertexArray<SparseFeatureVector<String>>>>();

		long startTime = System.nanoTime();
		for (LGraph<V, E> lg : set) {
			allLabeledWalkCounts.add(computeLabeledWalkCounts(lg));
		}
		long duration = System.nanoTime() - startTime;
		System.out.println("Computed vector space embedding "+duration/1000/1000+" [ms]");
		
		int n = set.size();
		double[][] r = new double[n][n];
		for (int i=0; i<n; i++) {
			ArrayList<VertexArray<SparseFeatureVector<String>>> e1 = allLabeledWalkCounts.get(i);
			for (int j=i; j<n; j++) {
				ArrayList<VertexArray<SparseFeatureVector<String>>> e2 = allLabeledWalkCounts.get(j);
				r[i][j] = r[j][i] = computeFromWalkCounts(e1, e2);
			}
		}
		
		return r;
	}
	
	private double computeFromWalkCounts(ArrayList<VertexArray<SparseFeatureVector<String>>> A, ArrayList<VertexArray<SparseFeatureVector<String>>> B) {
		
		Graph G = A.get(0).getGraph();
		Graph H = B.get(0).getGraph();
		double r=0;

		for (Vertex a : G.vertices()) {
			for (Vertex b : H.vertices()) {
				double abp = 0;
				double aap = 0;
				double bbp = 0;
				for (int l=0; l<=k; l++) {
					VertexArray<SparseFeatureVector<String>> VA = A.get(l);
					VertexArray<SparseFeatureVector<String>> VB = B.get(l);
					aap += VA.get(a).dotProduct(VA.get(a));
					bbp += VB.get(b).dotProduct(VB.get(b));
					double ab = VA.get(a).dotProduct(VB.get(b));
					abp += ab;
					r+=std.standardize(abp, aap, bbp, ab);
				}
			}
		}
		
		return r;
	}


	private ArrayList<VertexArray<SparseFeatureVector<String>>> computeLabeledWalkCounts(LGraph<V,E> lg) {
		Graph g = lg.getGraph();
		VertexArray<V> va = lg.getVertexLabel();
		EdgeArray<E> ea = lg.getEdgeLabel();
		ArrayList<VertexArray<SparseFeatureVector<String>>> walkCounts = new ArrayList<VertexArray<SparseFeatureVector<String>>>();
		VertexArray<SparseFeatureVector<String>> walks = new VertexArray<SparseFeatureVector<String>>(g);
		VertexArray<SparseFeatureVector<String>> newWalks = new VertexArray<SparseFeatureVector<String>>(g);
		
		// initialize with walks of size 0
		for (Vertex v : g.vertices()) {
			SparseFeatureVector<String> f = new SparseFeatureVector<String>();
			walks.set(v, f);
			String walk = String.valueOf(va.get(v));
			f.increaseByOne(walk);
		}
		walkCounts.add(walks);
		
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
			}
			walkCounts.add(newWalks);
			
			// prepare for next iteration
			walks = newWalks;
		}
		
		return walkCounts;
	}
	
}
