package chromo_activities;

import java.util.ArrayList;
import java.util.Random;

import chromosome.InitializeChromosome;
import graph.GraphMethods;
import iofile.ReadWriteFile;

public class Mutation {
	GraphMethods graph = new GraphMethods();

	public double[][] mutationTree(double[][] par, int num_Vertex, double mutation_Rate, int num_Cluster,
			int[][] vertex_In_Cluster, Random rnd) {
		double[][] child = new double[num_Vertex][num_Vertex];
		for (int i = 0; i < num_Vertex; i++) {
			for (int j = 0; j < num_Vertex; j++) {
				child[i][j] = par[i][j];
			}
		}
		for (int i = 0; i < num_Vertex; i++) {
			double r = 0 + (1 - 0) * rnd.nextDouble();
			if (r < mutation_Rate) {
				child = edgeClusteredTreeMutation(par, num_Vertex, num_Cluster, vertex_In_Cluster, rnd);
			}
		}
		return child;
	}

	public double[][] edgeClusteredTreeMutation(double[][] par, int num_Vertex, int num_Cluster,
			int[][] vertex_In_Cluster, Random rnd) {
		InitializeChromosome init_Chrome = new InitializeChromosome();

		double[][] child = new double[num_Vertex][num_Vertex];
		for (int i = 0; i < num_Vertex; i++) {
			for (int j = 0; j < num_Vertex; j++) {
				child[i][j] = par[i][j];
			}
		}

		int idx_Cluster = rnd.nextInt(num_Cluster);
		while (vertex_In_Cluster[idx_Cluster].length < 3) {
			idx_Cluster = rnd.nextInt(num_Cluster);

		}
		// 01. Chuyển ma trận của cluster thành ma trận cây để áp dụng đột biến
		double[][] weight_Cluster_Matrix;
		double[][] spanning_Tree_of_Cluster;
		int num_Vertex_in_Cluster = vertex_In_Cluster[idx_Cluster].length;
		weight_Cluster_Matrix = init_Chrome.createWeightMatrixForCluster(num_Vertex, num_Vertex_in_Cluster, par,
				vertex_In_Cluster[idx_Cluster]);
		spanning_Tree_of_Cluster = edgeMutation2(weight_Cluster_Matrix, num_Vertex_in_Cluster, rnd);

		// Chuyen ra cay khung cua do thi G
		int[] cluster = new int[num_Vertex_in_Cluster];
		for (int i = 0; i < num_Vertex_in_Cluster; i++) {
			cluster[i] = vertex_In_Cluster[idx_Cluster][i];
		}
		for (int k = 0; k < num_Vertex_in_Cluster; k++) {
			for (int j = 0; j < num_Vertex_in_Cluster; j++) {
				// child[cluster[k], cluster[j]] = 0;
				child[cluster[k]][cluster[j]] = spanning_Tree_of_Cluster[k][j];
			}
		}
		return child;
	}

	public double[][] edgeMutation(double[][] par, int num_Vertex, Random rnd) {
		double[][] child = new double[num_Vertex][num_Vertex];
		for (int i = 0; i < num_Vertex; i++) {
			for (int j = 0; j < num_Vertex; j++) {
				child[i][j] = par[i][j];
			}
		}
		int start_Vertex, end_Vertex; // 2 dinh cua canh ngau nhien them vao
		start_Vertex = rnd.nextInt(num_Vertex);
		end_Vertex = rnd.nextInt(num_Vertex);
		while ((start_Vertex == end_Vertex) || (child[start_Vertex][end_Vertex] > 0)) {
			end_Vertex = rnd.nextInt(num_Vertex);
			start_Vertex = rnd.nextInt(num_Vertex);
		}

		// Tìm duong ti tu dinh start_Vertex -> end_Vertex
		Boolean[] visited = new Boolean[num_Vertex];
		int[] pre = new int[num_Vertex];
		for (int i = 0; i < num_Vertex; i++) {
			visited[i] = false;
			pre[i] = -1;
		}
		// Tìm đường đi nối từ đỉnh start_Vertex --> end_Vertex
		graph.findCyclic(start_Vertex, end_Vertex, child, num_Vertex, visited, pre);
		ArrayList<Integer> path = graph.printPath(start_Vertex, end_Vertex, pre);
		// Xóa cạnh bất kỳ trên đường đi.
		int del_idx_1 = rnd.nextInt(path.size() - 1);
		// int del_idx_2 = rnd.Next(path.Count);
		int del_idx_2 = del_idx_1 + 1;
		child[path.get(del_idx_1)][path.get(del_idx_2)] = 0f;
		child[path.get(del_idx_2)][path.get(del_idx_1)] = 0f;

		child[path.get(del_idx_1)][path.get(del_idx_2)] = 0f;
		child[path.get(del_idx_2)][path.get(del_idx_1)] = 0f;
		// Dat canh moi
		child[start_Vertex][end_Vertex] = 1f;
		child[end_Vertex][start_Vertex] = 1f;
		return child;
	}

	public double[][] edgeMutation2(double[][] par, int num_Vertex, Random rnd) {
		double[][] child = new double[num_Vertex][num_Vertex];
		for (int i = 0; i < num_Vertex; i++) {
			for (int j = 0; j < num_Vertex; j++) {
				child[i][j] = par[i][j];
			}
		}
		int start_Vertex, end_Vertex; // 2 dinh cua canh ngau nhien them vao
		start_Vertex = rnd.nextInt(num_Vertex);
		end_Vertex = rnd.nextInt(num_Vertex);
		while ((start_Vertex == end_Vertex) || (child[start_Vertex][end_Vertex] > 0)) {
			end_Vertex = rnd.nextInt(num_Vertex);
			start_Vertex = rnd.nextInt(num_Vertex);
		}

		// if (rnd.nextDouble() < mutation_Rate) {
		// Tìm duong ti tu dinh start_Vertex -> end_Vertex
		Boolean[] visited = new Boolean[num_Vertex];
		int[] pre = new int[num_Vertex];
		for (int i = 0; i < num_Vertex; i++) {
			visited[i] = false;
			pre[i] = -1;
		}
		// Tìm đường đi nối từ đỉnh start_Vertex --> end_Vertex
		graph.findCyclic(start_Vertex, end_Vertex, child, num_Vertex, visited, pre);
		ArrayList<Integer> path = graph.printPath(start_Vertex, end_Vertex, pre);
		// Xóa cạnh dài nhất trên đường đi.
		double max = 0;
		int del_idx_1 = 0, del_idx_2 = 1;
		for (int i = 0; i < path.size() - 1; i++) {
			if (ReadWriteFile.weightMatrix[path.get(i)][path.get(i + 1)] > max) {
				del_idx_1 = i;
				del_idx_2 = i + 1;
				max = ReadWriteFile.weightMatrix[path.get(i)][path.get(i + 1)];
			}
		}
		child[path.get(del_idx_1)][path.get(del_idx_2)] = 0f;
		child[path.get(del_idx_2)][path.get(del_idx_1)] = 0f;
		// Dat canh moi
		child[start_Vertex][end_Vertex] = 1f;
		child[end_Vertex][start_Vertex] = 1f;
		// }
		return child;
	}

}
