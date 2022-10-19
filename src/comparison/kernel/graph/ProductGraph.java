package comparison.kernel.graph;

import java.util.ArrayList;
import java.util.Arrays;

import datastructure.Pair;
import graph.AdjListGraph.AdjListEdge;
import graph.Graph;
import graph.ImplementationHelper;


/**
 * A light-weight implementation of the graph interface used for 
 * efficient product graph creation, which is essential for kernel 
 * computation.
 * 
 * Note: Does not support attached properties, all relevant properties
 * are contained as attributes!
 * 
 * @author kriege
 *
 */
public class ProductGraph implements Graph {
	
	ArrayList<PVertex> vertices;
	ArrayList<PEdge> edges;
	Graph g1;
	Graph g2;
	
	double[] vertexWeights;
	double[][] edgeWeights;
	
	
	public ProductGraph(Graph g1, Graph g2) {
		this.g1 = g1;
		this.g2 = g2;
		vertices = new ArrayList<PVertex>();
		edges = new ArrayList<PEdge>();
	}
	
	public Iterable<PEdge> edges() {
		return edges;
	}

	public Iterable<PVertex> vertices() {
		return vertices;
	}

	public PVertex getVertex(int index) {
		return vertices.get(index);
	}
	
	public PEdge getEdge(int index) {
		return edges.get(index);
	}

	/**
	 * Note: This method has runtime O( min{deg(u),deg(v)} )
	 */
	public AdjListEdge getEdge(Vertex u, Vertex v) {
		return (AdjListEdge)ImplementationHelper.getEdge(u, v);
	}
	
	/**
	 * @see #getEdge(graph.Graph.Vertex, graph.Graph.Vertex)
	 */
	public boolean hasEdge(Vertex u, Vertex v) {
		return getEdge(u, v) != null;		
	}
	
	public int getVertexCount() {
		return vertices.size();
	}
	
	public int getEdgeCount() {
		return edges.size();
	}
	
	public PVertex createVertex(Vertex v1, Vertex v2, double weight) {
		return new PVertex(v1, v2, weight);
	}
	
	public PEdge createEdge(Vertex u, Vertex v, double weight) {
		return new PEdge((PVertex)u, (PVertex)v, weight);
	}
	
	public String toString() {
		return ImplementationHelper.toString(this);
	}
	
	public Graph getG1() {
		return g1;
	}
	
	public Graph getG2() {
		return g2;
	}
	
	/**
	 * Note: This method is only provided for convenience and is not
	 * light-weight, i.e. the first call might require significant 
	 * computation time. 
	 */
	public double[] getVertexWeights() {
		if (vertexWeights != null) {
			return vertexWeights;
		} else {
			return getVertexWeights(new double[getVertexCount()]);
		}
	}
	
	/**
	 * Fills the given array with the vertex weights. Note that for the
	 * provided data structure <code>data.length >= n</code> must hold,
	 * where <code>n=getVertexCount()</code>. Only the first n elements 
	 * are modified by this method.
	 * @see #getVertexWeights()
	 * @param data the provided data structure
	 * @return data array with vertex weights filled in
	 */
	public double[] getVertexWeights(double[] data) {
		for (PVertex v : vertices) {
			data[v.getIndex()] = v.getWeight();
		}
		vertexWeights = data;
		return vertexWeights;
	}
	
	/**
	 * Note: This method is only provided for convenience and is not
	 * light-weight, i.e. the first call might require significant 
	 * computation time. 
	 */
	public double[][] getEdgeWeights() {
		if (edgeWeights != null) {
			return edgeWeights;
		} else {
			int n = getVertexCount();
			double[][] data = new double[n][n];
			return getEdgeWeights(data, false);
		}
	}
	
	/**
	 * Fills the given array with the edge weights. Note that for the
	 * provided data structure <code>data.length >= n</code> and
	 * <code>data[i].length >= n</code>, for all 0<=i<n, must hold,
	 * where <code>n=getVertexCount()</code>. Only the first n elements 
	 * are modified by this method.
	 * @see #getVertexWeights()
	 * @param data the provided data structure
	 * @param fill fill elemnts with zero first
	 * @return data array with vertex weights filled in
	 */
	public double[][] getEdgeWeights(double[][] data, boolean fill) {
		int n = getVertexCount();
		if (fill) {
			for (int i=0; i<n; i++) {
				Arrays.fill(data[i], 0, n, 0d);
			}
		}
		for (PEdge e : edges()) {
			int iU = e.getFirstVertex().getIndex();
			int iV = e.getSecondVertex().getIndex();
			data[iV][iU] = data[iU][iV] = e.getWeight();
		}
		edgeWeights = data;
		return edgeWeights;
	}


	
	public class PVertex extends Pair<Vertex,Vertex> implements Vertex {
		
		int index;
		double weight;
		ArrayList<PEdge> edges;
		
		public PVertex(Vertex v1, Vertex v2, double weight) {
			super(v1, v2);
			this.weight = weight;
			this.index = vertices.size();
			vertices.add(this);
			edges = new ArrayList<PEdge>();
		}
		
		public ArrayList<PEdge> edges() {
			return edges;
		}
		
		public Iterable<PVertex> neighbors() {
			return ImplementationHelper.createNeighborIterator(this);
		}
		
		public int getDegree() {
			return edges.size();
		}
		
		public int getIndex() {
			return index;
		}
		
		protected void addEdge(PEdge e) {
			edges.add(e);
		}
		
		public double getWeight() {
			return weight;
		}
		
		public void setWeight(double weight) {
			this.weight = weight;
		}
		
	}
	
	public class PEdge implements Edge {
		int index;
		double weight;
		private PVertex u, v;
		
		public PEdge(PVertex u, PVertex v, double weight) {
			this.u = u;
			this.v = v;
			this.weight = weight;
			this.index = edges.size();
			edges.add(this);
			u.addEdge(this);
			v.addEdge(this);
		}
		
		public PVertex getFirstVertex() {
			return u;
		}
		
		public PVertex getSecondVertex() {
			return v;
		}
		
		public PVertex getOppositeVertex(Vertex w) {
			return (w == u) ? v : u;
		}
		
		public int getIndex() {
			return index;
		}
		
		public double getWeight() {
			return weight;
		}
		
		public void setWeight(double weight) {
			this.weight = weight;
		}

	}

	public String[] getProperties() {
		throw new UnsupportedOperationException("This class does not support properties.");
	}

	public Object getProperty(String key) {
		//throw new UnsupportedOperationException("This class does not support properties.");
		return null;
	}

	public void setProperty(String key, Object value) {
		throw new UnsupportedOperationException("This class does not support properties.");
	}

}
