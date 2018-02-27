package graph;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;

import clustered.Cluster;
import clustered.Edge;
import iofile.ReadWriteFile;


public class DijkstraHeap {
	private double[] d = new double[ReadWriteFile.numberOfCity];
	private int[] p = new int[ReadWriteFile.numberOfCity];
	private int[] K = new int[ReadWriteFile.numberOfCity];

	public DijkstraHeap() {
		Arrays.fill(d, 0);
		Arrays.fill(K, 0);
		Arrays.fill(p, 0);
	}

	public int parent(int i) {
		return (i + 1) / 2 - 1;
	}

	public int lChild(int i) {
		int x = i + 1;
		int y = x * 2;
		return y - 1;
	}

	public int rChild(int i) {
		int x = i + 1;
		int y = x * 2 + 1;
		return y - 1;
	}

	private void Build_Min_Heap(ArrayList<Double> Q, Map<Integer, Integer> pos) {
		int n = Q.size();
		for (int i = n / 2 - 1; i > -1; i--) {
			Min_Heapify(Q, i, n, pos);
		}
	}

	private void Min_Heapify(ArrayList<Double> A, int i, int n, Map<Integer, Integer> pos) {
		int l = lChild(i);
		int r = rChild(i);
		int minimum;
		if (l < n && A.get(l) < A.get(i)) {
			minimum = l;
		} else {
			minimum = i;
		}

		if (r < n && A.get(r) < A.get(minimum)) {
			minimum = r;
		}

		if (minimum != i) {

			double tmp = A.get(minimum);
			A.set(minimum, A.get(i));
			A.set(i, tmp);
			int tmp1 = pos.get(i);
			pos.put(i, pos.get(minimum));
			pos.put(minimum, tmp1);
			Min_Heapify(A, minimum, n, pos);
		}

	}

	private int Extract_Min(ArrayList<Double> Q, Map<Integer, Integer> pos) {
		int n = Q.size();

		int k = pos.get(0);
		Q.set(0, Q.get(n - 1));
		int tmp = pos.get(n - 1);
		pos.put(n - 1, pos.get(0));
		pos.put(0, tmp);
		Q.remove(n - 1);
		Min_Heapify(Q, 0, n - 1, pos);
		return k;
	}

	private void Decrease_Key(ArrayList<Double> Q, int m, double d, Map<Integer, Integer> pos) {
		Q.set(m, d);
		while (m > 0) {
			if (Q.get(parent(m)) > Q.get(m)) {
				int c = parent(m);
				double tmp = Q.get(c);
				Q.set(c, Q.get(m));
				Q.set(m, tmp);
				int tmp1 = pos.get(m);
				pos.put(m, pos.get(c));
				pos.put(c, tmp1);
				m = c;

			} else {
				break;
			}
		}

	}

	public ArrayList<Integer> DsKe(int v, Cluster clt1, Cluster clt2) {
		ArrayList<Integer> L1 = new ArrayList<>();
		for (Edge e : clt1.getListEdge()) {
			if (e.getStartVertex() == v && !L1.contains(e.getEndVertex())) {
				L1.add(e.getEndVertex());
			} else if (e.getEndVertex() == v && !L1.contains(e.getStartVertex())) {
				L1.add(e.getStartVertex());
			}
		}
		for (Edge e : clt2.getListEdge()) {
			if (e.getStartVertex() == v && !L1.contains(e.getEndVertex())) {
				L1.add(e.getEndVertex());
			} else if (e.getEndVertex() == v && !L1.contains(e.getStartVertex())) {
				L1.add(e.getStartVertex());
			}
		}
		return L1;
	}

	public Cluster Al_Dijkstra_Heap(Cluster clt1, Cluster clt2, int start) {

		for (int u : clt1.getCluster()) {
			d[u] = Double.MAX_VALUE;
			p[u] = 0;
			K[u] = 0;
		}
		d[start] = 0;
		ArrayList<Double> Q = new ArrayList<>();
		for (int i = 0; i < clt1.getCluster().size(); i++) {
			Q.add(0.0);
		}

		int j = 0;
		Map<Integer, Integer> pos = new HashMap<Integer, Integer>();
		for (int i : clt1.getCluster()) {
			Q.set(j, d[i]);
			pos.put(j, i);
			j++;
		}
		Build_Min_Heap(Q, pos);
		while (Q.size() != 0) {
			int u = Extract_Min(Q, pos);
			K[u] = 1;
			ArrayList<Integer> Keof_u = DsKe(u, clt1, clt2);
			for (int v : Keof_u) {
				if (K[v] == 0 && d[v] > d[u] + ReadWriteFile.weightMatrix[u][v]) {
					d[v] = d[u] + ReadWriteFile.weightMatrix[u][v];
					p[v] = u;
					int m = 0;
					for (int i : pos.keySet()) {
						if (pos.get(i) == v) {
							m = i;
							break;
						}
					}
					Decrease_Key(Q, m, d[v], pos);
				}
			}

		}

		ArrayList<Edge> NewListEdge = new ArrayList<Edge>();
		ArrayList<Integer> T = new ArrayList<Integer>();
		for (int i = 0; i < clt1.getCluster().size(); i++) {
			T.add(clt1.getCluster().get(i));
		}
		// (ArrayList) clt1.getCluster().clone();
		T.remove(T.indexOf(start));
		for (int i : T) {
			NewListEdge.add(new Edge(p[i], i));
		}
		return new Cluster(clt1.getCluster(), NewListEdge);
	}

	public double[][] convertClusterToArray(Cluster clt, double[][] edgeMatrix) {
		for (int i = 0; i < clt.getListEdge().size(); i++) {
			edgeMatrix[clt.getListEdge().get(i).getStartVertex()][clt.getListEdge().get(i).getEndVertex()] = 1;
			edgeMatrix[clt.getListEdge().get(i).getEndVertex()][clt.getListEdge().get(i).getStartVertex()] = 1;
		}
		return edgeMatrix;
	}

	public ArrayList<Cluster> convertArrayToCluster(double[][] edgeMatrix, int[][] vertexInCluster) {
		ArrayList<Cluster> clusters = new ArrayList<Cluster>();
		for (int i = 0; i < vertexInCluster.length; i++) {

			ArrayList<Edge> edge = new ArrayList<Edge>();
			for (int j = 0; j < vertexInCluster[i].length; j++) {
				for (int k = j + 1; k < vertexInCluster[i].length; k++) {
					if (edgeMatrix[vertexInCluster[i][j]][vertexInCluster[i][k]] > 0) {
						Edge e = new Edge(vertexInCluster[i][j], vertexInCluster[i][k]);
						edge.add(e);
					}
				}
			}
			ArrayList<Integer> temp = new ArrayList<Integer>();
			for(int j = 0; j< vertexInCluster[i].length; j++){
				temp.add(vertexInCluster[i][j]);
			}
			Cluster c = new Cluster(temp, edge);
			clusters.add(c);
		}
		return clusters;
	}
}
