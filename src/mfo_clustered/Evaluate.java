package mfo_clustered;

import java.util.LinkedList;
import java.util.Queue;

public class Evaluate {
	GraphMethods graph_Method_Class = new GraphMethods();
	Mutation mutation_Class = new Mutation();
	InitializeChromosome init_Chromo = new InitializeChromosome();

	public double evaluation(double [][] tree, double[][] weightMatrix, int num_vertex, int startVertex) {
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

	/*********************************************************************************************************************************************
	 * Decoding cho bai toan cay khung Giam tu n dinh ve thanh m dinh 01. Xóa
	 * các đỉnh và cạnh liên thuộc với nó > m 02. Tìm các thành phần liên thông
	 * 03. Nối các thành phần liên thông thứ i -> i + 1 + Chọn ngẫu nhiên ở mỗi
	 * thành phần liên thông 1 đỉnh + Chọn đỉnh đầu tiên
	 ********************************************************************************************************************************************/
	public double[][] decodingMFOVertexInSubGraph(double[][] ind_Matrix, int max_Genes, int num_Gen_of_Task_j) {
		double[][] tmp_Matrix = new double[max_Genes][max_Genes];
		// 01. Tìm các đỉnh lớn hơn num_Gen_of_Task_j và xóa nó cùng cạnh liên
		// thuộc
		for (int i = 0; i < num_Gen_of_Task_j; i++) {
			for (int j = 0; j < num_Gen_of_Task_j; j++) {
				if ((ind_Matrix[i][j] > 0) && (ind_Matrix[i][j] < Double.MAX_VALUE)) {
					tmp_Matrix[i][j] = 1.0f;
				} else {
					tmp_Matrix[i][j] = 0.0f;
				}
			}
		}
		// 02. Tìm các thành phần liên thông va lay moi tp lien thong 1 dinh
		int[] tp_LT = graph_Method_Class.get_Vertex_In_Each_SubGraph(tmp_Matrix, num_Gen_of_Task_j);

		// 03. Nối các thành phần liên thông thứ i -> i + 1
		for (int i = 0; i < tp_LT.length - 1; i++) {
			tmp_Matrix[tp_LT[i]][tp_LT[i + 1]] = 1.0f;
			tmp_Matrix[tp_LT[i + 1]][tp_LT[i]] = 1.0f;
		}

		// Tao ra ma tran ket qua
		double[][] final_Matrix = new double[num_Gen_of_Task_j][num_Gen_of_Task_j];
		for (int i = 0; i < num_Gen_of_Task_j; i++) {
			for (int j = 0; j < num_Gen_of_Task_j; j++) {
				final_Matrix[i][j] = (int) tmp_Matrix[i][j];
			}
		}
		return final_Matrix;
	}

}
