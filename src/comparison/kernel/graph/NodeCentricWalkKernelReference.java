package comparison.kernel.graph;

import java.util.ArrayList;
import comparison.kernel.Kernel;
import comparison.kernel.basic.DiracKernel;
import comparison.kernel.graph.ProductGraph.PEdge;
import comparison.kernel.graph.ProductGraph.PVertex;
import datastructure.Pair;
import graph.LGraph;
import graph.LGraphTools;
import jeigen.DenseMatrix;
import jeigen.SparseMatrixLil;


/**
 * This implementation closely follows the pseudo code of the publication,
 * but does not implement any optimizations.
 * 
 * @author kriege
 *
 * @param <V>
 * @param <E>
 */
public class NodeCentricWalkKernelReference<V,E> implements Kernel<LGraph<V, E>> {
	
	ProductGraphBuilder<V, E, ProductGraph> pgb;
	int k;
	double alpha;
	double beta;
	boolean wl;
	
	
	/**
	 * Create the node centric walk kernel.
	 * 
	 * @param k number of iterations/walk length
	 * @param alpha neighborhood strictness
	 * @param beta walk counts
	 * @param edgeKernel
	 * @param vertexKernel
	 * @param wl enable WL expressive power (for suitable alpha/beta)
	 */
	public NodeCentricWalkKernelReference(int k, double alpha, double beta, Kernel<? super E> edgeKernel, Kernel<? super V> vertexKernel, boolean wl) {
//		if (!(vertexKernel instanceof OneKernel) || !(edgeKernel instanceof OneKernel)) throw new IllegalArgumentException("Supports only OneKernel");
		pgb = new WeightedDirectProductGraphBuilderSparse<V, E>(edgeKernel, vertexKernel);
		
		this.k = k;
		this.alpha = alpha;
		this.beta = beta;
		this.wl = wl;
	}
	
	@Override
	public double compute(LGraph<V, E> g1, LGraph<V, E> g2) {
		int n1 = g1.getGraph().getVertexCount();
		int n2 = g2.getGraph().getVertexCount();
		
		LGraph<V,E> G = LGraphTools.copyLGraph(g1);
		LGraphTools.copyLGraph(g2, G, false);
		
		Pair<ArrayList<double[][]>,ArrayList<double[][]>> KP = computeGeneralizedLWalkKernel(G);
		ArrayList<double[][]> lK = KP.getFirst();
		ArrayList<double[][]> lKhat = KP.getSecond();
		
		double r=0;
		for (int i=0; i<lK.size(); i++) {
			double[][] K = lK.get(i);
			double[][] Khat = lKhat.get(i);
//			System.out.println(n1+" "+n2+" "+K.length);
			for (int u=0; u<n1; u++) {
				for (int v=n1; v<n1+n2; v++) {
					r+=Khat[u][v]*Math.pow(K[u][v], beta);
				}
			}
		}
		
		return r;
	}

	@Override
	public String getID() {
		if (pgb.getVertexKernel() instanceof DiracKernel && pgb.getEdgeKernel() instanceof DiracKernel) {
			return "NCWRef"+(wl? "WL":"")+"_"+k; 
		} else {
			return "NCWRef"+(wl? "WL":"")+"_"+k+"_"+pgb.getVertexKernel().getID()+"_"+pgb.getEdgeKernel().getID();
		}
	}
	
	private Pair<ArrayList<double[][]>,ArrayList<double[][]>> computeGeneralizedLWalkKernel(LGraph<V, E> lg) {
		ArrayList<double[][]> K = new ArrayList<double[][]>();
		ArrayList<double[][]> Khat = new ArrayList<double[][]>();
		int n = lg.getGraph().getVertexCount();

		ProductGraph PG = pgb.build(lg, lg);
		// create adjacency matrix with original indices
		SparseMatrixLil AX = new SparseMatrixLil(n*n, n*n);
		for (PEdge e : PG.edges()) {
			int i = e.getFirstVertex().getFirst().getIndex()*n+e.getFirstVertex().getSecond().getIndex();
			int j = e.getSecondVertex().getFirst().getIndex()*n+e.getSecondVertex().getSecond().getIndex();
			AX.append(i, j, 1d);
			AX.append(j, i, 1d);
		}
		
//		SparseMatrixLil A = LinAlgTools.createAdjacencyMatrix(G);
//		SparseMatrixLil AX = LinAlgTools.kroneckerProduct(A, A);
		
		DenseMatrix w = new DenseMatrix(AX.rows, 1);
		double[][] r = new double[n][n];
		for (PVertex v : PG.vertices()) {
			w.set(v.getFirst().getIndex()*n+v.getSecond().getIndex(), 0, 1d);
			r[v.getFirst().getIndex()][v.getSecond().getIndex()] = 1d;
		}
		DenseMatrix wp = w;
		
		// K^(1)
		K.add(r);
		Khat.add(r);
		
		for (int i=1; i <= k; i++) {
			w = AX.mmul(w);
			wp=wp.add(w);
			r = new double[n][n]; 
			double[][] rHat = new double[n][n];
			for (int s=0; s<n; s++) {
				for (int t=0; t<n; t++) {
					r[s][t] = w.get(s*n+t, 0);
//					System.out.print(wp.get(s*n+s, 0)+"\t"+wp.get(t*n+t, 0)+"\t"+wp.get(s*n+t, 0)+"\t"+w.get(s*n+t, 0));
					rHat[s][t] = comp(alpha, wp.get(s*n+s, 0), wp.get(t*n+t, 0), wp.get(s*n+t, 0));
					if (wl) w.set(s*n+t, 0, rHat[s][t]);
//					System.out.println(" -> "+r[s][t]);
				}
			}
			K.add(r);
			Khat.add(rHat);
		}
		
		return new Pair<ArrayList<double[][]>, ArrayList<double[][]>>(K, Khat);
	}
	
	private static double comp(double alpha,double aap, double bbp, double abp) {
		// compute squared distance using kernel trick
		double d = aap + bbp - 2*abp;
		// compute Gaussian RBF kernel
		return Math.exp(-alpha*d);
	}
        
}
