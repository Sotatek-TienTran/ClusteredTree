package mfo_clustered;

import java.util.ArrayList;

public class Evaluate {
	GraphMethods graph_Method_Class = new GraphMethods();
	Mutation mutation_Class = new Mutation();
	InitializeChromosome init_Chromo = new InitializeChromosome();

	public double clusteredTreeEvaluate(double[][] edge_Matrix, double[][] weigh_Matrix, int num_Vertex,
			int start_Vertex) {
		double path_Length = 0;
		int[] pre = new int[num_Vertex];

		ArrayList<Integer> path_1 = new ArrayList<Integer>();
		// Tạo ma trận trọng của cây khung từ ma trận canh và ma trận trọng số
		// để tìm đường đi ngắn nhất bằng dijstra
		double[][] temp_Matrix = new double[num_Vertex][num_Vertex];
		for (int i = 0; i < num_Vertex; i++) {
			for (int j = 0; j < num_Vertex; j++) {
				if ((edge_Matrix[i][j] > 0) && (edge_Matrix[i][j] < Double.MAX_VALUE))// neu
																						// co
																						// canh
				{
					temp_Matrix[i][j] = weigh_Matrix[i][j];
				} else {
					temp_Matrix[i][j] = 0;
				}
			}
		}

		pre = graph_Method_Class.dijkstra(temp_Matrix, num_Vertex, start_Vertex);
		for (int i = 1; i < num_Vertex; i++) {
			path_1 = graph_Method_Class.printPath(start_Vertex, i, pre);
			for (int l = path_1.size() - 1; l > 0; l--) {
				path_Length = path_Length + temp_Matrix[path_1.get(l)][path_1.get(l - 1)];
			}
		}
		return path_Length;
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
