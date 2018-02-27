package evaluate;

import java.util.ArrayList;
import java.util.LinkedList;
import java.util.Queue;

import chromo_activities.Mutation;
import chromosome.InitializeChromosome;
import graph.GraphMethods;

public class Evaluate {
	GraphMethods graph = new GraphMethods();
	Mutation mutation = new Mutation();
	InitializeChromosome initChromo = new InitializeChromosome();

	public double evaluation(double[][] tree, double[][] weightMatrix, int num_vertex, int startVertex) {
		double[] distances = new double[num_vertex];// distance between root and
													// the others
		double sum = 0;
		distances[startVertex] = 0;
		boolean[] mark = new boolean[num_vertex];
		Queue<Integer> queue = new LinkedList<>();
		for (int i = 0; i < num_vertex; i++) {
			mark[i] = true;
		}
		queue.add(startVertex);
		while (!queue.isEmpty()) {
			int u = queue.poll();
			mark[u] = false;
			for (int i = 0; i < num_vertex; i++) {
				if (tree[u][i] > 0 && mark[i]) {
					queue.add(i);
					mark[i] = false;
					distances[i] = distances[u] + weightMatrix[u][i];
					sum += distances[i];
				}
			}
		}
		return sum;
	}

	public double evaluation2(double[][] edge_Matrix, double[][] weigh_Matrix, int num_Vertex, int start_Vertex) {
		double path_Length = 0;
		int[] pre = new int[num_Vertex];

		ArrayList<Integer> path_1;
		// Tạo ma trận trọng của cây khung từ ma trận canh và ma trận trọng số
		// để tìm đường đi ngắn nhất bằng dijstra
		double[][] temp_Matrix = new double[num_Vertex][num_Vertex];
		for (int i = 0; i < num_Vertex; i++) {
			for (int j = 0; j < num_Vertex; j++) {
				if ((edge_Matrix[i][j] > 0) && (edge_Matrix[i][j] < Double.MAX_VALUE)) // neu
																						// co
																						// canh
				{
					temp_Matrix[i][j] = weigh_Matrix[i][j];
				} else {
					temp_Matrix[i][j] = 0;
				}
			}
		}
		pre = graph.dijkstra(temp_Matrix, num_Vertex, start_Vertex);
		for (int i = 1; i < num_Vertex; i++) {
			path_1 = graph.printPath(start_Vertex, i, pre);
			for (int l = path_1.size() - 1; l > 0; l--) {
				path_Length = path_Length + temp_Matrix[path_1.get(l)][path_1.get(l - 1)];
			}
		}
		return path_Length;
	}
	
	
}
