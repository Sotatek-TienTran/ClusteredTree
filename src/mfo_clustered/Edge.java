package mfo_clustered;

public class Edge {
	private int startVertex;
	private int endVertex;

	public Edge(int s_Vertex, int e_Vertex) {
		startVertex = s_Vertex;
		endVertex = e_Vertex;
	}

	// getter
	public int getStartVertex() {
		return startVertex;
	}

	// setter
	public void setStartVertex(int value) {
		this.startVertex = value;
	}

	// getter
	public int getEndVertex() {
		return endVertex;
	}

	// setter
	public void setEndVertex(int value) {
		this.endVertex = value;
	}
}
