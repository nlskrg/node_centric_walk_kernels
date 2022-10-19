package comparison.kernel.graph;

import comparison.kernel.Kernel;
import comparison.kernel.graph.ProductGraph.PVertex;
import graph.Graph;
import graph.Graph.Edge;
import graph.Graph.Vertex;
import graph.LGraph;
import graph.properties.EdgeArray;
import graph.properties.VertexArray;

/**
 * @inheritDocs
 * The direct product graphs connects vertex pairs with common adjacency only. 
 */
public class WeightedDirectProductGraphBuilderSparse<V,E> extends ProductGraphBuilder<V, E, ProductGraph> {
	
	/**
	 * @inheritDocs
	 */
	public WeightedDirectProductGraphBuilderSparse(Kernel<? super E> edgeKernel, Kernel<? super V> vertexKernel) {
		super(edgeKernel, vertexKernel);
	}

	/**
	 * @inheritDocs
	 */
	public ProductGraph build(LGraph<V,E> lg1, LGraph<V,E> lg2) {
		Graph g1 = lg1.getGraph();
		Graph g2 = lg2.getGraph();
		VertexArray<V> va1 = lg1.getVertexLabel();
		VertexArray<V> va2 = lg2.getVertexLabel();
		EdgeArray<E> ea1 = lg1.getEdgeLabel();
		EdgeArray<E> ea2 = lg2.getEdgeLabel();
		int n1 = g1.getVertexCount();
		int n2 = g2.getVertexCount();
		
		ProductGraph pg = new ProductGraph(g1, g2);
		// compute compatible vertices
		PVertex[][] m = new PVertex[n1][n2];
		for (Vertex v1 : g1.vertices()) {
			V l1 = va1.get(v1);
			for (Vertex v2 : g2.vertices()) {
				V l2 = va2.get(v2);
				double value = vertexKernel.compute(l1, l2);
				if (value != 0) {
					PVertex pv = pg.createVertex(v1, v2, value);
					m[v1.getIndex()][v2.getIndex()] = pv;
				}
			}
		}

		for (PVertex pv : pg.vertices()) {
			Vertex v1 = pv.getFirst();
			Vertex v2 = pv.getSecond();
			for (Edge e1 : v1.edges()) {
				Vertex w1 = e1.getOppositeVertex(v1);
				// handle each edge only from one vertex
				if (v1.getIndex() > w1.getIndex()) continue;
				for (Edge e2 : v2.edges()) {
					Vertex w2 = e2.getOppositeVertex(v2);
					PVertex pw = m[w1.getIndex()][w2.getIndex()];
					if (pw != null) {
						double value = edgeKernel.compute(ea1.get(e1), ea2.get(e2));
						if (value != 0) {
							pg.createEdge(pv, pw, value);
						}
					}
				}
			}
		}
		
		return pg;
	}

}
