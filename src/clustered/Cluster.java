package clustered;

import java.util.ArrayList;

import clustered.Edge;

public class Cluster {
	private ArrayList<Integer> cluster = new ArrayList<Integer>();
	private ArrayList<Edge> listEdge = new ArrayList<Edge>();

	public ArrayList<Integer> getCluster() {
		return cluster;
	}

	public void setCluster(ArrayList<Integer> cluster) {
		this.cluster = cluster;
	}

	public ArrayList<Edge> getListEdge() {
		return listEdge;
	}

	public void setListEdge(ArrayList<Edge> listEdge) {
		this.listEdge = listEdge;
	}

	public Cluster(ArrayList<Integer> cluster, ArrayList<Edge> newListEdge) {
		this.cluster = cluster;
		this.listEdge = newListEdge;
	}

	public Cluster() {

	}
}
