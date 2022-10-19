package algorithm.shortestpath.betweenness;

import java.util.ArrayList;
import java.util.List;
import java.util.Stack;

import org.teneighty.heap.BinaryHeap;
import org.teneighty.heap.Heap;
import org.teneighty.heap.Heap.Entry;

import algorithm.shortestpath.SPTools;
import datastructure.FullMatrix;
import datastructure.Matrix;
import graph.Graph;
import graph.Graph.Edge;
import graph.Graph.Vertex;
import graph.properties.EdgeArray;
import graph.properties.VertexArray;

/**
 * Computes matrices for GraphHopper kernel inspired by Brandes
 * betweenness algorithm.
 * 
 * @author Nils Kriege
 *
 */
public class BrandesEdgeLengthHopper {

	public static int MAX_INT = 1000;
	

	/**
	 * 
	 * @param g
	 * @param lambda null for uniform edge length
	 * @return
	 */
	public VertexArray<Matrix<Long>> compute(Graph g, EdgeArray<Integer> lambda) {
		
		// we do not determine the diameter here, but shrink matrices later for 
		// the graph hopper kernel
		int diam = SPTools.getDiameter(g); // TODO this is not correct for weighted graphs
		
		// initialize output
		VertexArray<Matrix<Long>> cB = new VertexArray<Matrix<Long>>(g);
		for (Vertex v : g.vertices()) {
			Matrix<Long> M = new FullMatrix<Long>(diam, diam);
			M.fill(0L);
			cB.set(v, M);
		}
		if (diam==0) return cB;

		// initialize data for loop; cleaned in each iteration
		Heap<Integer,Vertex> Q = new BinaryHeap<Integer, Vertex>();
		Stack<Vertex> S = new Stack<Vertex>();

		VertexArray<Integer> dist = new VertexArray<Integer>(g);
		VertexArray<List<Vertex>> Pred = new VertexArray<List<Vertex>>(g);
		for (Vertex v : g.vertices()) {
			Pred.set(v, new ArrayList<Vertex>());
		}
		VertexArray<Long> sigma = new VertexArray<Long>(g);
		VertexArray<Double> delta = new VertexArray<Double>(g);
		VertexArray<Entry<Integer, Vertex>> queueEntries = new VertexArray<Entry<Integer,Vertex>>(g);
		VertexArray<Matrix<Long>> Hops = new VertexArray<Matrix<Long>>(g);
		VertexArray<Matrix<Long>> backHops = new VertexArray<Matrix<Long>>(g);
		for (Vertex v : g.vertices()) {
			Hops.set(v, new FullMatrix<Long>(diam, diam));
			backHops.set(v, new FullMatrix<Long>(diam, diam));
		}
		
		for (Vertex s : g.vertices()) {
			// single-source shortest-paths problem
			
			// initialize data structures
			for (Vertex v : g.vertices()) {
				Pred.get(v).clear();
				dist.set(v, Integer.MAX_VALUE);
				sigma.set(v, 0L);
				Hops.get(v).fill(0L);
				backHops.get(v).fill(0L);
			}
			dist.set(s, 0);
			sigma.set(s, 1L);
			Hops.get(s).set(0, 0, 1L);
			Entry<Integer, Vertex> qE = Q.insert(0, s);
			queueEntries.set(s, qE);
			
			while (!Q.isEmpty()) {
				Vertex v = Q.extractMinimum().getValue();
				queueEntries.set(v, null);
				S.push(v);
				for (Edge e : v.edges()) {
					Vertex w = e.getOppositeVertex(v);
					// path discovery --- w found for the first time?
					int newdist = dist.get(v) + ((lambda == null) ? 1 : lambda.get(e));
					if (dist.get(w) > newdist) {
						dist.set(w, newdist);
						if (queueEntries.get(w) != null) {
							Q.decreaseKey(queueEntries.get(w), newdist);
						} else {
							qE = Q.insert(newdist, w);
							queueEntries.set(w, qE);
						}
						sigma.set(w, 0L);
						Hops.get(w).fill(0L);
						Pred.get(w).clear();
					}
					// path counting
//					System.out.println(dist.get(w) +" "+newdist+" >? "+(dist.get(w) > newdist)+ " ==? "+(dist.get(w) == newdist));
					if (dist.get(w) == newdist) {
						sigma.set(w, sigma.get(w)+sigma.get(v));
						addShifted(Hops.get(w), Hops.get(v));
						Pred.get(w).add(v);
					}
					
				}
			}
			
			// accumulation --- back-propagation of dependencies
			for (Vertex v : g.vertices()) {
				delta.set(v, 1d);
				backHops.get(v).set(0, 0, 1L);
			}
			while (!S.isEmpty()) {
				Vertex w = S.pop();
				for (Vertex v : Pred.get(w)) {
					delta.set(v, delta.get(v)+delta.get(w));
					addShifted(backHops.get(v), backHops.get(w));
				}
				for (int i=0; i<diam; i++) {
					for (int j=i; j<diam; j++) {
						long value = Hops.get(w).get(i, 0);
						value *= backHops.get(w).get(j-i, 0);
						cB.get(w).set(i, j, cB.get(w).get(i, j)+value);
					}
				}
			}
		}
		
		return cB;
	}
	
	private void addShifted(Matrix<Long> A, Matrix<Long> B) {
		for (int i=1; i<A.getRowDimension(); i++) {
			for (int j=0; j<A.getColumnDimension(); j++) {
				A.set(i, j, A.get(i, j)+B.get(i-1, j));
			}
		}
	}
	
	
	/**
	 * Works with double weights which are discretized in order to make
	 * equal-length path identifiable.
	 * @param g
	 * @param lambda
	 * @return
	 */
	public static EdgeArray<Integer> discretizeDouble(Graph g, EdgeArray<Double> lambda) {
		
		double maxD = Double.NEGATIVE_INFINITY;
		for (Edge e : g.edges()) {
			maxD = Math.max(maxD, lambda.get(e));
		}
		
		EdgeArray<Integer> lambda2 = new EdgeArray<Integer>(g);
		for (Edge e : g.edges()) {
			int value = (int)Math.floor(lambda.get(e)/maxD*MAX_INT);
			lambda2.set(e, value);
		}
		
		return lambda2;
	}

	
}
